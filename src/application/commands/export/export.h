#ifndef CUBIQUITY_VOXELIZE_EXPORT_H
#define CUBIQUITY_VOXELIZE_EXPORT_H

#include "base/materials.h"

#include "utility.h"

#include "flags.h"

bool exportVolume(const flags::args& args);
void saveVolumeAsImages(Cubiquity::Volume& volume, const MaterialSet& materials, const std::string& filename);

#endif // CUBIQUITY_VOXELIZE_EXPORT_H