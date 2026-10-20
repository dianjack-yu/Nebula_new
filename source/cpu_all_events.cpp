#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"
#include "common/cli_params.h"
#include "common/work_pool.h"
#include "common/time_log.h"
#include "io/output_stream.h"
#include "io/load_tri_file.h"
#include "io/load_pri_file.h"

#include "drivers/cpu/all_events_cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>

// Main typedefs
using geometry_t = nbl::geometry::octree<false>;
using material_t = material<scatter_physics<false>>;

using driver = nbl::drivers::all_events_cpu_driver<
	scatter_physics<false>,
	intersect_t,
	geometry_t
>;

int main(int argc, char **argv)
{
	// Print version information
	std::clog << "This is Nebula All Events 32-bit version "
		<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n\n"
		"Physics models:\n";
	scatter_physics<false>::print_info(std::clog);
	intersect_t::print_info(std::clog);
	std::clog << "\n" << std::string(80, '-') << "\n\n";

	// Settings
	cli_params p("[options] <geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]");
	p.add_option("energy-threshold", "Lowest energy to simulate", 0);
	p.add_option("seed", "Random seed", 0x78c7e39b);
	p.add_option("detect-filename", "Filename for detected electrons", "detected.bin");
	p.add_option("events-filename", "Filename for event records", "stdout");
	// Dian Yu: Variable phonon domain, skip all or specify one for each material
	p.add_multi_option("phonon-domain", "Phonon scattering energy domains per material (eV)", std::string(""));
	//Dian Yu: for variable starting material, vacuum by default
	p.add_option("initial-material", "Material index where primaries start (-123=vacuum, -126=detector, >=0=material)", -123);

	p.parse(argc, argv);

	const real energy_threshold = p.get_flag<real>("energy-threshold");
	const unsigned int seed = p.get_flag<unsigned int>("seed");
	const std::string detect_filename = p.get_flag<std::string>("detect-filename");
	const std::string events_filename = p.get_flag<std::string>("events-filename");
	// Dian Yu: fetch initial material
	const int initial_material = p.get_flag<int>("initial-material");
	// Dian Yu: Parse phonon domains, fetch as a vector
	std::vector<real> phonon_domains = p.get_flag_vector<real>("phonon-domain");


	// Setup time logging
	time_log timer;

	// Interpret command-line options
	std::vector<std::string> pos_flags = p.get_positional();
	if (pos_flags.size() < 3)
	{
		p.print_usage(std::clog);
		return 1;
	}

	std::mt19937 random_generator(seed);

	// Load geometry
	std::clog << "Loading geometry..." << std::endl;
	timer.start();
	std::vector<triangle> triangles = nbl::load_tri_file(pos_flags[0]);
	timer.stop("Loading triangles");

	if (triangles.empty())
	{
		std::clog << "Error: could not load triangles!\n";
		p.print_usage(std::clog);
		return 1;
	}
	// Sanity check with number of materials
	{
		int max_material = -1;
		for (triangle const & tri : triangles)
		{
			if (tri.material_in > max_material)
				max_material = tri.material_in;
			if (tri.material_out > max_material)
				max_material = tri.material_out;
		}

		if (max_material > int(pos_flags.size())-3)
		{
			std::clog << "Error: not enough materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
			return 1;
		}
		if (max_material < int(pos_flags.size())-3)
		{
			std::clog << "Warning: too many materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
		}
	}

	timer.start();
	geometry_t geometry = geometry_t::create(triangles);
	timer.stop("Building acceleration structure");


	// Load materials
	std::clog << "Loading materials..." << std::endl;
	timer.start();
	nbl::cpu_material_manager<material_t> materials;
	for (size_t parameter_idx = 2; parameter_idx < pos_flags.size(); ++parameter_idx)
	{
		nbl::hdf5_file material(pos_flags[parameter_idx]);
		
		// Dian Yu: Get phonon domain threshold for this material (default 200 eV if not specified)
		size_t mat_idx = parameter_idx - 2;
		real phonon_domain = (mat_idx < phonon_domains.size()) //  && phonon_domains[mat_idx] > 0
							? phonon_domains[mat_idx] 
							: 200.0;
		
		materials.add(material, phonon_domain);  // Dian Yu: Pass phonon domain

		std::clog << "  Material " << mat_idx << ":\n"
			"    Name: " << material.get_property_string("name") << "\n"
			"    cstool version: " << material.get_property_string("cstool_version") << "\n"
			"    Phonon domain: " << phonon_domain << " eV\n";  // Dian Yu: Print phonon domain
	}
	timer.stop("Loading materials");


	// Load primaries
	std::clog << "Loading primary electrons..." << std::endl;
	timer.start();
	std::vector<particle> primaries; std::vector<int2> pixels;
	std::tie(primaries, pixels) = nbl::load_pri_file(pos_flags[1], geometry.AABB_min(), geometry.AABB_max(), materials.get_max_energy());
	timer.stop("Loading primary electrons");

	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n";
		p.print_usage(std::clog);
		return 1;
	}

	// The driver only accepts uint32 tags. So we make a map: simulation tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1

	// This manages the work to be done (thread-safe).
	work_pool pool(primaries.data(), gpu_tags.data(), primaries.size());


	intersect_t inter;

	// Print debug data
	std::clog << "\n"
		<< "Loaded " << triangles.size() << " triangles.\n"
		<< "  min = {" << geometry.AABB_min().x << ", " << geometry.AABB_min().y << ", " << geometry.AABB_min().z << "}\n"
		<< "  max = {" << geometry.AABB_max().x << ", " << geometry.AABB_max().y << ", " << geometry.AABB_max().z << "}\n"
		<< "Loaded " << primaries.size() << " primaries.\n"
		<< "Loaded " << materials.size() << " materials.\n\n" << std::flush;

	// Prepare output file
	output_stream detect_out(detect_filename);
	output_stream events_out(events_filename);

	// Simulation loop
	auto sim_loop = [&pool, &detect_out, &events_out, &pixels,
		&geometry, &inter, &materials, energy_threshold, initial_material](uint32_t seed) // Dian Yu: pass initial material
		{
		driver d(
			inter, materials, geometry,
			energy_threshold, materials.get_max_energy(), seed, initial_material); //Dian Yu: pass initial material
		output_buffer detect_buff(detect_out,
			1024*(8*sizeof(float) + 3*sizeof(int))); //Dian Yu: include origin_z and lineage tag
		output_buffer events_buff(events_out,
			1024*1024*(6*sizeof(float) + 4*sizeof(int))); //Dian Yu: include origin_z, lineage tag and event ID

		for (;;)
		{
			// Push new particles
			auto work_data = pool.get_work(1);

			if (std::get<2>(work_data) == 0)
				break;

			d.push(
				std::get<0>(work_data),  // particle*
				std::get<1>(work_data),  // tag*
				std::get<2>(work_data)); // number

			// Simulate a little, output event data
			d.simulate_to_end([&events_buff, &pixels](particle const & before, particle const & after, uint32_t tag)
			                  {
				events_buff.add(std::array<float, 6>{
					after.pos.x, after.pos.y, after.pos.z,
					before.kin_energy, before.kin_energy - after.kin_energy, after.origin_z});
				events_buff.add(std::array<int, 4>{ // Jonathan Choi & Dian Yu edit
					pixels[tag].x, pixels[tag].y, after.tag_lineage,
				    after.event_id_int}); }); // Jonathan Choi & Dian Yu edit

			// Flush output data
			d.flush_detected([&detect_buff, &pixels](particle const &p, uint32_t t)
							 {
				detect_buff.add(std::array<float, 8>{
					p.pos.x, p.pos.y, p.pos.z,
					p.dir.x, p.dir.y, p.dir.z, p.kin_energy,
				    p.origin_z}); // Jonathan Choi & Dian Yu edit
				detect_buff.add(std::array<int, 3>{ 
					pixels[t].x, pixels[t].y, p.tag_lineage}); }); // Jonathan Choi & Dian Yu edit, no p.event_id_int
			
		}

		detect_buff.flush();
		events_buff.flush();
    };

	// Simulation
	//const auto n_threads = std::thread::hardware_concurrency();
	const auto n_threads = 1;
	std::clog << "Creating " << n_threads << " CPU drivers" << std::endl;
	std::vector<std::thread> threads;

	timer.start();
	for (unsigned int i = 0; i < n_threads; ++i)
		threads.push_back(std::thread(sim_loop, random_generator()));

	// Progress indicator
	for (;;)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		auto primaries_to_go = pool.get_primaries_to_go();
		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size())) << "%";
		if (primaries_to_go == 0)
			break;
	}

	for (auto& t : threads)
		t.join();

	timer.stop("Simulation");

	std::clog << "\n\n";
	timer.print(std::clog);
	return 0;
}