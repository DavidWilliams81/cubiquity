#ifndef CUBIQUITY_APP_GENERATE_H
#define CUBIQUITY_APP_GENERATE_H

#include "base/types.h"

#include <filesystem>

enum class Algorithm {
    fractal_noise,
    menger_sponge,
    worley_noise,
};

bool generateVolume(Algorithm algorithm,
                    const std::filesystem::path& output_path,
                    int size);

#endif // CUBIQUITY_APP_GENERATE_H