#ifndef __PARTICLE_H_
#define __PARTICLE_H_

/// Represents a particle (currently, always an electron) in the simulation.
struct particle
{
	vec3 pos;		 ///< Position
	vec3 dir;		 ///< Direction (unnormalised).
	real kin_energy; ///< Kinetic energy
    float origin_z;      // Dian Yu electron escape depth stats, records z-position where this particle's current state originated
	// Jonathan Choi ____ new fields _____
	int tag_x = 0;		 // beam pixel X index
	int tag_y = 0;		 // beam pixel Y index
	int tag_lineage = 0; // electron lineage tag (0, 1, 11, 21, 12, 22, 13, 23)
	int event_id_int = 0; // Dian Yu scatter event ID
};

#endif // __PARTICLE_H_
