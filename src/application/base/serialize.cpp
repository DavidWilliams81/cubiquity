#include "serialize.h"

#include "base/bounds.h"
#include "base/logging.h"

#include "cubiquity.h"

bool checkInputFileIsValid(const std::filesystem::path& inputFile)
{
	if (!std::filesystem::exists(inputFile)) {
		log_error("Input '{}' does not exist!", inputFile);
		return false;
	}
	if (!std::filesystem::is_regular_file(inputFile)) {
		log_error("Input '{}' is not a regular file!", inputFile);
		return false;
	}
	return true;
}

bool checkOutputDirIsValid(const std::filesystem::path& outputDir)
{
	if (!std::filesystem::exists(outputDir)) {
		log_error("Output '{}' does not exist!", outputDir);
		return false;
	}
	if (!std::filesystem::is_directory(outputDir)) {
		log_error("Output '{}' is not a directory!", outputDir);
		return false;
	}
	return true;
}

std::filesystem::path getMetadataPath(std::filesystem::path volumePath)
{
	return volumePath.replace_extension(".txt");
}

std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata> loadVolume(const std::filesystem::path& volume_path)
{
	std::unique_ptr<Cubiquity::Volume> volume =
		std::make_unique<Cubiquity::Volume>(volume_path.string());

	Metadata metadata;
	metadata.load(getMetadataPath(volume_path));

	// A .dag file does not need dimensions in the metadata and will ignore any
	// that are present. However, they usually *are* present to facilitate
	// sharing metadata between .dag and .bin files. It seems prudent to at
	// least warn the user if the values do not appear to be correct (though
	// perhaps there is a valid reason, such as a border of empty voxels?).
	if (metadata.dimensions.has_value()) {
		auto [lower, upper] = find_bounds(*volume);
		ivec3 real_dimensions = (upper - lower) + ivec3(1);
		log_warning_if(metadata.dimensions != real_dimensions,
			"Dimensions in metadata {} do not match computed values {}",
			metadata.dimensions.value(), real_dimensions);
	}

	// TODO - We could potentially check that a material has been found for each
	// material which exists in the volume, and perhaps that there are no unused
	// materials? Can also check the dimensions of the volume are as expected.

	return { std::move(volume), metadata };
}

void saveVolume(const std::filesystem::path& volume_path,
	            Cubiquity::Volume& volume, Metadata& metadata)
{
	// A .dag file does not actually need dimensions to be known as they can
	// easily be computed. However, they are needed for a .bin file. We choose
	// to require them for a .dag file anyway because it simplifies the process
	// of exporting a .dag to a .bin (otherwise the dimensions need to be added,
	// overwritting the TOML file with the same name, and we should prompt the
	// user as to whether this is what they really wanted).
	if (metadata.dimensions.has_value() == false) {
		throw std::runtime_error("Metadata is missing dimensions.");
	}

	volume.save(volume_path.string());
	metadata.save(getMetadataPath(volume_path));
}
