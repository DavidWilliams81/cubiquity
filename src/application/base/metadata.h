#ifndef CUBIQUITY_METADATA_H
#define CUBIQUITY_METADATA_H

#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "json.hpp"

#include <cstdint>
#include <filesystem>

typedef std::array<float, 3> Col; // To be renamed when we can avoid conflicts.;

struct Material
{
	std::string name;
	Col diffuse = { 1.0f, 1.0f, 1.0f };
};

struct Metadata
{
	static const Material EmptySpace;
	static const Material Default;
	static const Material Warning;

	std::vector<Material> materials;
};

Metadata loadMetadataForVolume(const std::filesystem::path& volumePath);
void saveMetadataForVolume(const Metadata& metadata, const std::filesystem::path& volumePath);

#endif // CUBIQUITY_METADATA_H