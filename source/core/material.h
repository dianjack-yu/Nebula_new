#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include "scatter_list.h"
#include "../io/hdf5_file.h"

template<typename...>
struct material;

/**
 * \brief Material class.
 *
 * Template parameter is a ::scatter_list, i.e. a list of physical scattering
 * mechanisms. These scattering mechansims contain the physics of what actually
 * goes on inside the material.
 */
template<typename... scatter_types>
struct material<scatter_list<scatter_types...>>
	: public scatter_list<scatter_types...>
{
public:
	using base_t = scatter_list<scatter_types...>;
	
	/**
	 * \brief Dian Yu: Read from a HDF5 material file with default phonon domain.
	 *
	 * \param mat HDF5 file with the data
	 */
	static CPU material create(nbl::hdf5_file const & mat);
	
	/**
	 * \brief Dian Yu: Read from a HDF5 material file with custom phonon domain.
	 *
	 * \param mat HDF5 file with the data
	 * \param phonon_domain Energy threshold for phonon scattering (eV)
	 */
	static CPU material create(nbl::hdf5_file const & mat, real phonon_domain);
	
	/**
	 * \brief Clone from a different material.
	 */
	template<typename... source_scatter_types>
	static CPU material create(material<scatter_list<source_scatter_types...>> const & source);
	
	/**
	 * \brief Deallocate all data
	 */
	static CPU void destroy(material & mat);
	
	/**
	 * \brief Barrier energy, that is, the work function plus the Fermi energy.
	 */
	real barrier = 0;
};

#include "material.inl"

#endif // __MATERIAL_H_