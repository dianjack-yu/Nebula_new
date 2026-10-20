#ifndef __CLI_PARAMS_H_
#define __CLI_PARAMS_H_
/**
 * \brief Very simple command-line parameter parser.
 *
 * Types of arguments supported are:
 *  * `--key=value`
 *  * `--key value` (single value)
 *  * `--key val1 val2 val3` (multiple values, use get_flag_vector)
 *  * Anything that doesn't match this pattern are "positional arguments".
 */

#include <string>
#include <vector>
#include <map>
#include <sstream>

class cli_params
{
public:
	/**
	 * \brief Constructor.
	 *
	 * \param param_string String to be printed for the documentation listing.
	 */
	cli_params(std::string const & param_string = "[options]")
		: param_string(param_string)
	{}
	
	/**
	 * \brief Add key-value option pair.
	 */
	template<typename T>
	void add_option(
		std::string const & key,
		std::string const & description,
		T const & default_value)
	{
		std::stringstream def;
		def << default_value;
		options[key] = { description, def.str(), false };
	}
	
	/**
	 * \brief Dian Yu: Add key-value option pair that accepts multiple values.
	 * 
	 * When parsing --key val1 val2 val3, all values until the next flag
	 * will be collected. Use get_flag_vector<T>() to retrieve.
	 */
	template<typename T>
	void add_multi_option(
		std::string const & key,
		std::string const & description,
		T const & default_value)
	{
		std::stringstream def;
		def << default_value;
		options[key] = { description, def.str(), true };
	}
	
	/**
	 * \brief Parse command-line options.
	 *
	 * Parameters are the same as what you get in the C++ `main()` function.
	 *
	 * This function parses the arguments and stores them in private variables.
	 *
	 * \param argc Number of command-line arguments.
	 * \param argv The command-line arguments themselves.
	 */
	void parse(int argc, char**argv)
	{
		program_name = argv[0];
		for (int i = 1; i < argc; )  // REMOVE ++i from here!
		{
			std::string item = argv[i];
			++i;  // Increment at START of loop instead

			if (item.size() >= 2 && item[0] == '-' && item[1] == '-')
			{
				std::string key;
				std::string value;
				auto split = item.find("=", 2);
				
				if (split == std::string::npos)
				{
					// Format: --key value(s)
					key = item.substr(2);

					auto map_it = options.find(key);
					if (map_it == options.end())
						throw std::runtime_error("Unexpected key '" + key + "'");
					
					// Check if this is a multi-value option
					if (map_it->second.is_multi)
					{
						// Collect all values until next flag
						std::stringstream value_stream;
						bool first = true;
						while (i < argc)  // Changed from i + 1
						{
							std::string next_item = argv[i];

							// Stop if we hit another flag
							if (next_item.size() >= 2 && next_item[0] == '-' && next_item[1] == '-')
								break;
							
							if (!first)
								value_stream << " ";
							value_stream << next_item;
							first = false;
							++i;  // Consume this value
						}
						value = value_stream.str();

						// If no values collected, use default
						if (value.empty())
							value = map_it->second.value;
					}
					else
					{
						// Single value only
						if (i < argc)  // Changed from i + 1
						{
							std::string next_item = argv[i];

							// Check if next item is a flag
							if (next_item.size() >= 2 && next_item[0] == '-' && next_item[1] == '-')
							{
								// Next is a flag, use default value
								value = map_it->second.value;
							}
							else if (next_item == "=")
							{
								// Handle: --key = value
								++i;  // Skip "="
								if (i < argc)
								{
									value = argv[i];
									++i;  // Consume value
								}
								else
								{
									value = map_it->second.value;
								}
							}
							else
							{
								// Normal case: --key value
								value = argv[i];
								++i;  // Consume value
							}
						}
						else
						{
							// No more arguments, use default
							value = map_it->second.value;
						}
					}

					map_it->second.value = value;
				}
				else
				{
					// Format: --key=value
					key = item.substr(2, split - 2);
					value = item.substr(split + 1);

						auto map_it = options.find(key);
					if (map_it == options.end())
						throw std::runtime_error("Unexpected key '" + key + "'");
					map_it->second.value = value;
				}
			}
			else
			{
				positional.push_back(item);
			}
		}
	}
	
	/**
	 * \brief Print documentation string, updated by Dian Yu
	 */
	void print_usage(std::ostream & out)
	{
		out << "Usage: " << program_name << " " << param_string << "\n"
			<< "Options:\n";
		for (auto const & opt : options)
		{
			out << "\t--" << opt.first << ": " << opt.second.description
				<< " [" << opt.second.value << "]";
			if (opt.second.is_multi)
				out << " (multi-value)";
			out << "\n";
		}
	}
	
	/**
	 * \brief Get the program name, i.e. argv[0].
	 */
	std::string const & get_program_name() const
	{
		return program_name;
	}
	
	/**
	 * \brief Get a single-value flag.
	 *
	 * \param key Flag name
	 */
	template<typename T>
	T get_flag(std::string const & key) const
	{
		return options.at(key).get_as<T>();
	}
	
	/**
	 * \brief Dian Yu: Get a multi-value flag as a vector.
	 *
	 * For options added with add_multi_option(), this parses the space-separated
	 * values into a vector of the requested type.
	 *
	 * \param key Flag name
	 */
	template<typename T>
	std::vector<T> get_flag_vector(std::string const & key) const
	{
		return options.at(key).get_as_vector<T>();
	}
	
	/**
	 * \brief Get the positional arguments (in the correct order).
	 *
	 * The "positional" arguments are arguments that don't follow the --key=value format.
	 */
	std::vector<std::string> const & get_positional() const
	{
		return positional;
	}

private:
	struct option
	{
		std::string description;
		std::string value;
		bool is_multi;  // Dian Yu: Whether this option accepts multiple values
		
		template<typename T>
		T get_as() const
		{
			T ret;
			std::stringstream(value) >> ret;
			return ret;
		}
		
		template<typename T>
		std::vector<T> get_as_vector() const
		{
			std::vector<T> ret;
			std::istringstream iss(value);
			T val;
			while (iss >> val)
				ret.push_back(val);
			return ret;
		}
	};
	
	std::string program_name;
	std::string param_string;
	std::map<std::string, option> options;
	std::vector<std::string> positional;
};

#endif // __CLI_PARAMS_H_
