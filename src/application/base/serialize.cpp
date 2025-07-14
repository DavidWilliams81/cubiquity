#include "serialize.h"

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

	// TODO - We could potentially check that a material has
	// been found for each material which exists in the volume. 

	return { std::move(volume), std::move(metadata) };
}

void saveVolume(const std::filesystem::path& volume_path,
	            Cubiquity::Volume& volume, Metadata& metadata)
{
	volume.save(volume_path.string());
	metadata.save(getMetadataPath(volume_path));
}
