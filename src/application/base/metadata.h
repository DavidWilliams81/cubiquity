#ifndef CUBIQUITY_METADATA_H
#define CUBIQUITY_METADATA_H

#include "types.h"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

class Material
{
	friend class Metadata; // For loading/saving materials

public:

	static const Material Default;
	static const Material EmptySpace;

	Material(std::optional<std::string> name = std::nullopt,
		     std::optional<vec3> base_color  = std::nullopt)
		:m_name(name), m_base_color(base_color) {}

	std::string name() const { return m_name.value_or(""); }
	void set_name(std::optional<std::string> name) { m_name = name; }

	vec3 base_color() const { return m_base_color.value_or(DefaultColor); }
	void set_base_color(std::optional<vec3> base_color) { m_base_color = base_color; }

private:
	inline static const vec3 DefaultColor = { 0.8f, 0.8f, 0.8f };

	std::optional<std::string> m_name;
	std::optional<vec3> m_base_color;
};

class Metadata
{
public:
	Metadata() {};
	Metadata(const std::filesystem::path& path) { load(path); };

	void load(const std::filesystem::path& path);
	void save(const std::filesystem::path& path) const;

public:
	std::optional<ivec3>  dimensions;
	std::vector<Material> materials;
};

#endif // CUBIQUITY_METADATA_H
