#ifndef CUBIQUITY_APP_GENERATE_H
#define CUBIQUITY_APP_GENERATE_H

#include "base/types.h"

#include <filesystem>

bool generateVolume(const std::filesystem::path& output_path, uint size_exp);

#endif // CUBIQUITY_APP_GENERATE_H