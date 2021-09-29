#ifndef CUBIQUITY_PATHS_H
#define CUBIQUITY_PATHS_H

#include <filesystem>

bool checkInputFileIsValid(const std::filesystem::path& inputFile);
bool checkOutputDirIsValid(const std::filesystem::path& outputDir);

#endif // CUBIQUITY_PATHS_H


