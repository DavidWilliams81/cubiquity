#include "metadata.h"

#include "logging.h"

#include <fstream>

using namespace nlohmann; // For JSON

const Material Metadata::EmptySpace = { "_EmptySpace", { 0.0f, 0.0f, 0.0f} }; // FIXME - Should we actually add a diffuse color for this? Or just an opacity of zero?
const Material Metadata::Default    = { "_Default",    { 0.8f, 0.8f, 0.8f} }; // Grey
const Material Metadata::Warning    = { "_Warning",    { 1.0f, 0.0f, 1.0f} }; // Bright purple

std::filesystem::path getMetadataPath(std::filesystem::path volumePath)
{
	return volumePath.replace_extension(".json");
}

Metadata loadMetadataForVolume(const std::filesystem::path& volumePath)
{
	Metadata metadata;

	std::ifstream metadataFile(getMetadataPath(volumePath));
	if (!metadataFile)
	{
		log_error("Error: Failed to open metadata '{}'", getMetadataPath(volumePath));
		return metadata;
	}

	json metadataAsJSON;
	try {
		metadataAsJSON = json::parse(metadataFile);
	} catch (json::parse_error& ex) {
		log_error("Error: JSON parse error at byte {}", ex.byte);
		return metadata;
	}

	metadata.materials.resize(metadataAsJSON["materials"].size());
	for (int i = 0; i < metadataAsJSON["materials"].size(); i++)
	{
		metadata.materials[i].name =
			metadataAsJSON["materials"][i]["name"].template get<std::string>();
		metadata.materials[i].diffuse =
			metadataAsJSON["materials"][i]["diffuse"].template get<std::array<float, 3>>();
	}

	// FIXME - We could potentially check that a material has
	// been found for each material which exists in the volume.

	return metadata;
}

void saveMetadataForVolume(const Metadata& metadata, const std::filesystem::path& volumePath)
{
	ordered_json metadataAsJSON;
	for (int i = 0; i < metadata.materials.size(); i++)
	{
		metadataAsJSON["materials"][i]["name"] = metadata.materials[i].name;
		metadataAsJSON["materials"][i]["diffuse"] = metadata.materials[i].diffuse;
	}

	std::ofstream metdataFile(getMetadataPath(volumePath));
	metdataFile << std::setw(4) << metadataAsJSON << std::endl; // setw triggers pretty-printing.
}
