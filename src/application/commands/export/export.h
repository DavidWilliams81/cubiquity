#ifndef CUBIQUITY_EXPORT_H
#define CUBIQUITY_EXPORT_H

#include "base/metadata.h"

#include "utility.h"

#include "flags.h"

bool exportVolume(const flags::args& args);
void saveVolumeAsImages(Cubiquity::Volume& volume, const Metadata& metadata, const std::string& filename);

#endif // CUBIQUITY_EXPORT_H
