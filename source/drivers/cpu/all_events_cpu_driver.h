#ifndef __ALL_EVENTS_CPU_DRIVER_H_
#define __ALL_EVENTS_CPU_DRIVER_H_

#include "cpu_driver.h"

namespace nbl { namespace drivers {

/**
 * \brief CPU driver for tracking energy deposits during scattering.
 * \see cpu_driver
 */

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class all_events_cpu_driver
	: public cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>
{
public:
	using base_t = cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>;

	using typename base_t::particle_index_t;
	using typename base_t::material_t;
	using typename base_t::material_manager_t;
	using typename base_t::seed_t;

	/// Constructor
	all_events_cpu_driver(
		intersect_t intersect,
		material_manager_t const & materials,
		geometry_manager_t const & geometry,
		real min_energy, real max_energy,
		seed_t seed = util::random_generator<false>::default_seed,
		int initial_material = -123) //Dian Yu: default initial material
		: base_t(intersect, materials, geometry, min_energy, max_energy, seed, initial_material)
	{}

	/**
	 * \brief Perform a single iteration of the simulation for all particles.
	 *
	 * \param function Callback function to be called if a particle loses energy
	 *                 in a scattering event. Should have signature
	 *                 `void(particle const & before, particle const & after, uint32_t tag)`.
	 */
	template<typename deposit_function>
	void do_iteration(deposit_function function)
	{
		const auto particle_count = this->_particles.get_total_count();
		for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
		{
			if (!this->_particles.active(particle_idx))
				continue;

			this->init(particle_idx);

			// Dian Yu: Tracking both scattering and intersection
			bool track_scatter = this->_particles.next_scatter(particle_idx);
			bool track_intersect = this->_particles.next_intersect(particle_idx);
			bool track = track_scatter || track_intersect;
			const particle before = this->_particles[particle_idx];

			const auto count_before = this->_particles.get_total_count();

            this->intersect(particle_idx);
            this->scatter(particle_idx);

            const auto count_after = this->_particles.get_total_count();

            if (count_after > count_before)
            {
	            this->emit_secondary_creation_events(
		            count_before,
		            count_after,
		            function
	            );
            }

			const particle after = this->_particles[particle_idx];

			if (track) // Dian Yu: removed && before.kin_energy != after.kin_energy
				function(before, after,
					this->_particles.get_primary_tag(particle_idx));
		}
	}

	/// Keep simulating until there are no particles left.
	template<typename deposit_function>
	void simulate_to_end(deposit_function&& func)
	{
		for (particle_index_t particle_idx = 0;
			particle_idx < this->_particles.get_total_count(); ++particle_idx)
		{
			// Dian Yu: Record initial state (event ID 0)
			{
				const particle initial = this->_particles[particle_idx];
				particle initial_event = initial;
				initial_event.event_id_int = 0;  // Starting event
				func(initial_event, initial_event, this->_particles.get_primary_tag(particle_idx));
			}

			while (this->_particles.active(particle_idx))
			{
				this->init(particle_idx);

				// Dian Yu: Tracking both scattering and intersection
			    bool track_scatter = this->_particles.next_scatter(particle_idx);
			    bool track_intersect = this->_particles.next_intersect(particle_idx);
			    bool track = track_scatter || track_intersect;
				const particle before = this->_particles[particle_idx];

				const auto count_before = this->_particles.get_total_count();

				this->intersect(particle_idx);
				this->scatter(particle_idx);

				const auto count_after = this->_particles.get_total_count();

				if (count_after > count_before)
				{
					this->emit_secondary_creation_events(
						count_before,
						count_after,
						func
					);
				}

				const particle after = this->_particles[particle_idx];
			
				if (track) // Dian Yu: removed && before.kin_energy != after.kin_energy
					func(before, after,
						this->_particles.get_primary_tag(particle_idx));


			}
		}
	}

private:
	template<typename event_function>
	void emit_secondary_creation_events(
		particle_index_t first_new_idx,
		particle_index_t last_new_idx,
		event_function& func)
	{
		for (particle_index_t i = first_new_idx; i < last_new_idx; ++i)
		{
			const particle& p = this->_particles[i];

			// Creation is a zero-step event: before == after
			func(
				p,
				p,
				this->_particles.get_primary_tag(i)
			);
		}
	}

};
}} // namespace nbl::drivers

#endif // __ALL_EVENTS_CPU_DRIVER_H_


