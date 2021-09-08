#ifndef CUBIQUITY_VOXELIZE_EXPORT_H
#define CUBIQUITY_VOXELIZE_EXPORT_H

#include "utility.h"

void saveVolumeAsImages(Cubiquity::Volume& volume, const std::string& filename, Cubiquity::ProgressMonitor* progMon = nullptr);

#endif // CUBIQUITY_VOXELIZE_EXPORT_H