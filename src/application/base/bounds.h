#ifndef CUBIQUITY_APP_UTILS_H
#define CUBIQUITY_APP_UTILS_H

#include "cubiquity.h"
#include "storage.h"

#include "metadata.h"

#include <filesystem>
#include <utility>

std::pair<ivec3, ivec3> find_bounds(Cubiquity::Volume& volume);
ivec3 find_dimensions(Cubiquity::Volume& volume);

#endif // CUBIQUITY_APP_UTILS_H


