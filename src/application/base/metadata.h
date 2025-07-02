#ifndef CUBIQUITY_METADATA_H
#define CUBIQUITY_METADATA_H

#include "types.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>

// Forward declaration of toml::ordered_value
namespace toml
{
	struct ordered_type_config;

	template<typename TypeConfig>
	class basic_value;

	using ordered_value = basic_value<ordered_type_config>;
}

// Metadata is kept as TOML objects and converted to real types on demand. 
// When working with metadata, one option would be to parse the TOML and
// populate our own data structues, use/modify that data, and later write it
// back out as TOML if required.
// 
// However, we would also ideally keep track of which values were defaulted
// (rather than being specified) so that we could avoid writing them out,
// otherwise a terse input file might turn into a verbose output file. We'd also
// lose additional information (such as comments or other pieces of data) if
// the user had hand-crafted a TOML file.
//
// Therefore we instead keep the data as TOML types and convert on demand. I'm
// not really sure how valuable this but it does address the concerns above.
class Metadata
{
public:
	Metadata();
	~Metadata();

	// Non-copyable (deleted) as we hold a raw pointer
	Metadata(const Metadata&) = delete;
	Metadata& operator=(const Metadata&) = delete;

	// But moveable (defaulted)
	Metadata(Metadata&& m) = default;
	Metadata& operator=(Metadata&& m) = default;

	void load(const std::filesystem::path& path);
	void save(const std::filesystem::path& path) const;

	// Bounds
	std::optional<std::array<int, 3>> lower_bound() const;
	void set_lower_bound(const std::array<int, 3>& lower_bound);

	std::optional<std::array<int, 3>> upper_bound() const;
	void set_upper_bound(const std::array<int, 3>& upper_bound);

	// Materials
	int material_count() const;
	void clear_all_materials();

	void clear_material(int index);
	void set_material(int index, std::string name, vec3 base_color);

	std::string material_name(int index) const;
	void set_material_name(int index, std::string name);

	vec3 material_base_color(int index) const;
	void set_material_base_color(int index, vec3 color);

	// Explicitly adding enties for these makes the files more readable
	// (rather than relying on missing entries falling back on default values).
	void set_material_to_default(int index);
	void set_material_to_empty_space(int index);

private:

	void ensure_bounds_exists();
	void ensure_material_exists(int index);

	// It should not be necessary to call this externally
	// as the materials array is resized on demand.
	void resize_materials(int size);

	// Using a pointer to a forward-declared type so that we don't have to
	// include toml.hpp from this header and bloat compilation times. It is a 
	// raw pointer as usage is very simple, and as there are some constraints
	// when using an std::unique_ptr with only a forward declaration available.
	toml::ordered_value* m_toml_data;
};

#endif // CUBIQUITY_METADATA_H
