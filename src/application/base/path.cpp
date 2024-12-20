#include "paths.h"

#include "base/logging.h"

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
