#include "serialize.h"

#include "base/bounds.h"
#include "base/logging.h"

#include "cubiquity.h"

std::ifstream make_safe_ifstream(const std::string& path,
								 std::ios::openmode mode,
								 std::ios_base::iostate except_flags)
{
	std::ifstream file;
	file.exceptions(except_flags);
	file.open(path, mode);
	return file;
}

std::ofstream make_safe_ofstream(const std::string& path,
								 std::ios::openmode mode,
								 std::ios_base::iostate except_flags)
{
	std::ofstream file;
	file.exceptions(except_flags);
	file.open(path, mode);
	return file;
}


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

std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata> loadVolume(const std::filesystem::path& volume_path)
{
	std::unique_ptr<Cubiquity::Volume> volume =
		std::make_unique<Cubiquity::Volume>(volume_path.string());

	std::filesystem::path metadata_path = volume_path;
	metadata_path.replace_extension(".toml");
	Metadata metadata(metadata_path);

	// Metadata for a DAG should not contain dimensions.
	if (metadata.dimensions.has_value()) {
		log_warning("Metadata for a DAG should not contain dimensions. "
			        "They will be ignored");
		metadata.dimensions.reset();
	}

	// TODO - We could potentially check that a material has been found for each
	// material which exists in the volume, and perhaps that there are no unused
	// materials? Can also check the dimensions of the volume are as expected.

	return { std::move(volume), metadata };
}

void saveVolume(const std::filesystem::path& volume_path,
	            Cubiquity::Volume& volume, Metadata& metadata)
{
	// Metadata for a DAG should not contain dimensions. Throw an exception as
	// this function is only called by our own code, which should be fixed.
	if (metadata.dimensions.has_value()) {
		throw std::runtime_error(
			"Metadata for a DAG should not contain dimensions.");
	}

	std::filesystem::path metadata_path = volume_path;
	metadata_path.replace_extension(".toml");

	volume.save(volume_path.string());
	std::ofstream file(metadata_path);
	metadata.save(file);
}
