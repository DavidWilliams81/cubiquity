#ifndef CUBIQUITY_PATHS_H
#define CUBIQUITY_PATHS_H

#include "storage.h"

#include "metadata.h"

#include <filesystem>
#include <utility>

bool checkInputFileIsValid(const std::filesystem::path& inputFile);
bool checkOutputDirIsValid(const std::filesystem::path& outputDir);

// FIXME - Ideally we wouldn't need this declaration (could be hidden in the
// .cpp file), as only loadVolume and saveVolume should call it.
std::filesystem::path getMetadataPath(std::filesystem::path volumePath);

// FIXME - I'd rather return the volume by value but
// I need to make a working move constructor first.
std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata>
    loadVolume(const std::filesystem::path& vol_path);

void saveVolume(const std::filesystem::path& volume_path, 
                Cubiquity::Volume& volume, Metadata& metadata);

#endif // CUBIQUITY_PATHS_H


