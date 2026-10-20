#include <initializer_list>

namespace detail
{
	template<typename source_type, typename... target_scatter_types, size_t... is>
	scatter_list<target_scatter_types...> create_scatter_list(
		source_type const & source,
		nbl::index_sequence<is...>)
	{
		return scatter_list<target_scatter_types...>{
			target_scatter_types::create(nbl::tuple::get<is>(source))...
		};
	}
	
	template<typename... scatter_types, size_t... is>
	void destroy_scatter_list(
		scatter_list<scatter_types...> & list,
		nbl::index_sequence<is...>)
	{
		(void)std::initializer_list<int>{
			(scatter_types::destroy(nbl::tuple::get<is>(list)), 0)... };
	}
	
	// Dian Yu: Helper structure to detect if a scatter type has a two-parameter create method
	template<typename scatter_type>
	struct has_phonon_create
	{
	private:
		template<typename T>
		static auto test(int) -> decltype(T::create(std::declval<nbl::hdf5_file const &>(), std::declval<real>()), std::true_type{});
		
		template<typename>
		static std::false_type test(...);
		
	public:
		static constexpr bool value = decltype(test<scatter_type>(0))::value;
	};
	
	// Dian Yu: Create scatter with phonon domain if supported
	template<typename scatter_type>
	CPU typename std::enable_if<has_phonon_create<scatter_type>::value, scatter_type>::type
	create_single_scatter(nbl::hdf5_file const & mat, real phonon_domain)
	{
		return scatter_type::create(mat, phonon_domain);
	}
	
	// Dian Yu: Create scatter without phonon domain if not supported
	template<typename scatter_type>
	CPU typename std::enable_if<!has_phonon_create<scatter_type>::value, scatter_type>::type
	create_single_scatter(nbl::hdf5_file const & mat, real phonon_domain)
	{
		(void)phonon_domain;
		return scatter_type::create(mat);
	}
	
	// Dian Yu: Create scatter list with phonon domain support
	template<typename... target_scatter_types, size_t... is>
	scatter_list<target_scatter_types...> create_scatter_list_with_phonon(
		nbl::hdf5_file const & mat,
		real phonon_domain,
		nbl::index_sequence<is...>)
	{
		return scatter_list<target_scatter_types...>(
			create_single_scatter<target_scatter_types>(mat, phonon_domain)...
		);
	}
}

// Dian Yu: Original create method - maintains backward compatibility
template<typename... scatter_types>
CPU auto material<scatter_list<scatter_types...>>::create(nbl::hdf5_file const & mat)
	-> material
{
	material<base_t> target;
	target.barrier = static_cast<real>(mat.get_property_quantity("barrier") / nbl::units::eV);
	target.base_t::operator=(base_t(scatter_types::create(mat)...));
	return target;
}

// Dian Yu: New create method with phonon domain parameter
template<typename... scatter_types>
CPU auto material<scatter_list<scatter_types...>>::create(nbl::hdf5_file const & mat, real phonon_domain)
	-> material
{
	material<base_t> target;
	target.barrier = static_cast<real>(mat.get_property_quantity("barrier") / nbl::units::eV);
	target.base_t::operator=(
		detail::create_scatter_list_with_phonon<scatter_types...>(
			mat, phonon_domain, nbl::make_index_sequence<sizeof...(scatter_types)>{})
	);
	return target;
}

// Clone from another material
template<typename... scatter_types>
template<typename... source_scatter_types>
CPU auto material<scatter_list<scatter_types...>>::create(
	material<scatter_list<source_scatter_types...>> const & source)
	-> material
{
	material<base_t> target;
	target.barrier = source.barrier;
	target.base_t::operator=(
		detail::create_scatter_list<scatter_list<source_scatter_types...>, scatter_types...>(
			source, nbl::make_index_sequence<sizeof...(scatter_types)>{})
	);
	return target;
}

// Destroy material
template<typename... scatter_types>
CPU void material<scatter_list<scatter_types...>>::destroy(material & mat)
{
	detail::destroy_scatter_list(mat,
		nbl::make_index_sequence<sizeof...(scatter_types)>{});
}