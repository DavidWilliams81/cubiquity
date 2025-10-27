#ifndef CUBIQUITY_SIMPLEX_NOISE_H_
#define CUBIQUITY_SIMPLEX_NOISE_H_

#include "base/types.h"
#include "storage.h"

#include <cstdint>

float fractal_noise(int x, int y, int z, int octaves);

u8 worley_noise(int x, int y, int z, int cell_size, int border = 0);

#endif // CUBIQUITY_SIMPLEX_NOISE_H_