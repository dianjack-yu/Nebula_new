namespace nbl { namespace drivers {

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::cpu_driver(
	intersect_t intersect,
	material_manager_t const & materials,
	geometry_manager_t const & geometry,
	real min_energy, real max_energy,
	seed_t seed,
	int initial_material  //Dian Yu: initial material for pushed particles
) :
	_min_energy(min_energy), _max_energy(max_energy),
	_particles(particle_manager_t::create()),
	_materials(materials),
	_geometry(geometry),
	_intersect(intersect),
	_initial_material(initial_material), //Dian Yu: set initial material

	rand_state(seed)
{
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::~cpu_driver()
{
	particle_manager_t::destroy(_particles);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::push(
	particle* particles,
	primary_tag_t* tags,
	particle_index_t N
) -> particle_index_t
{
	return _particles.push(particles, tags, N, _initial_material); // Dian Yu: variable initial material
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_running_count() const
-> particle_index_t
{
	return _particles.get_running_count();
}
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_detected_count() const
-> particle_index_t
{
	return _particles.get_detected_count();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
template<typename detect_function>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::flush_detected(
	detect_function func)
{
	_particles.flush_detected(func);
	_particles.flush_terminated();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::init(particle_index_t particle_idx)
{
	// Get data from memory
	auto this_particle = _particles[particle_idx];

	// If not in domain, terminate
	if (!_geometry.in_domain(this_particle.pos))
	{
		_particles.terminate(particle_idx);
		return;
	}

	// Next scattering event
	scatter_event next_scatter{
		0,
		_geometry.get_max_extent()
	};

	// If not in a vacuum, get next scatter event
	auto this_material_idx = _particles.get_material_index(particle_idx);
	if (_materials.is_physical(this_material_idx))
	{
		const auto this_material = _materials[this_material_idx];

		// Terminate if we are below the energy threshold
		// (which is with respect to the vacuum energy)
		if (this_particle.kin_energy < this_material.barrier + _min_energy ||
			this_particle.kin_energy > this_material.barrier + _max_energy)
		{
			_particles.terminate(particle_idx);
			return;
		}

		// Sample next scattering event
		next_scatter = this_material.sample_path(this_particle, rand_state);
		// Dian Yu: handle case of no scattering events
		// If no scattering event is available (type == 0 and distance is invalid),
		// the particle should propagate without scattering until it hits a boundary
		if (next_scatter.type == 0 && next_scatter.distance < 0)
		{
			// Set distance to maximum extent so particle can propagate to geometry boundary
			next_scatter.distance = _geometry.get_max_extent();
		}

	}

	// Move particle to next event, unless there is a triangle in the way
	normalise(this_particle.dir);
	intersect_event next_intersect = _geometry.propagate(
		this_particle.pos, this_particle.dir, next_scatter.distance,
		_particles.get_last_triangle(particle_idx),
		_particles.get_material_index(particle_idx)
	);

	if (next_intersect.isect_triangle == nullptr)
	{
		// No triangle intersection: move to scattering position.
		// Scatter there later (after sorting)
		_particles.set_scatter_event(particle_idx, next_scatter);
	}
	else
	{
		// Triangle intersection: move to triangle position.
		_particles.set_intersect_event(particle_idx, next_intersect);
	}
}


template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::intersect(particle_index_t particle_idx)
{
	// ignore all particles except those with an intersect event.
	if (!_particles.next_intersect(particle_idx))
		return;

	_intersect.execute(_materials, _particles, particle_idx, rand_state);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::scatter(particle_index_t particle_idx)
{
	// ignore all particles except those with an inelastic event.
	if (!_particles.next_scatter(particle_idx))
		return;

	// forget last intersected triangle. This event might cause us to scatter
	// back into that triangle and we don't want to ignore that triangle if so.
	_particles.forget_last_triangle(particle_idx);

	_materials[_particles.get_material_index(particle_idx)].execute(
		_particles.get_next_scatter(particle_idx), _particles, particle_idx, rand_state);
}

}} // namespace nbl::drivers
