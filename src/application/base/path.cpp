#include "paths.h"

#include "base/logging.h"

bool checkInputFileIsValid(const std::filesystem::path& inputFile)
{
	if (!std::filesystem::exists(inputFile)) {
		log(Error, "Input \'%s\' does not exist!", inputFile.c_str());
		return false;
	}
	if (!std::filesystem::is_regular_file(inputFile)) {
		log(Error, "Input \'%s\' is not a regular file!", inputFile.c_str());
		return false;
	}
	return true;
}

bool checkOutputDirIsValid(const std::filesystem::path& outputDir)
{
	if (!std::filesystem::exists(outputDir)) {
		log(Error, "Output \'%s\' does not exist!", outputDir.c_str());
		return false;
	}
	if (!std::filesystem::is_directory(outputDir)) {
		log(Error, "Output \'%s\' is not a directory!", outputDir.c_str());
		return false;
	}
	return true;
}
