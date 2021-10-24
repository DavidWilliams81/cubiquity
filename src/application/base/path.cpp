#include "paths.h"

#include "base/logging.h"

bool checkInputFileIsValid(const std::filesystem::path& inputFile)
{
	if (!std::filesystem::exists(inputFile)) {
		log(Error, "Input \'", inputFile, "\' does not exist!");
		return false;
	}
	if (!std::filesystem::is_regular_file(inputFile)) {
		log(Error, "Input \'", inputFile, "\' is not a regular file!");
		return false;
	}
	return true;
}

bool checkOutputDirIsValid(const std::filesystem::path& outputDir)
{
	if (!std::filesystem::exists(outputDir)) {
		log(Error, "Output \'", outputDir, "\' does not exist!");
		return false;
	}
	if (!std::filesystem::is_directory(outputDir)) {
		log(Error, "Output \'", outputDir, "\' is not a directory!");
		return false;
	}
	return true;
}
