/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/

#include "generate.h"

#include "base/logging.h"
#include "base/metadata.h"
#include "base/noise.h"
#include "base/serialize.h"
#include "base/types.h"

#include "storage.h"
#include "utility.h"

// Work around missing std::execution support.
#ifdef CUBIQUITY_USE_POOLSTL
	#define POOLSTL_STD_SUPPLEMENT
	#define POOLSTL_STD_SUPPLEMENT_FORCE
	#include "poolstl.hpp"
#else
	#include <execution>
#endif // CUBIQUITY_USE_POOLSTL

#include <filesystem>
#include <mutex>
#include <numeric>

using Cubiquity::Volume;
using Cubiquity::uint;

using namespace std;

u8 mengerSponge(int x, int y, int z)
{
	int levels = 8;
	bool mat_id = 1;

	for (int i = 0; i < levels; i++)
	{
		int count = 0;
		if ((x - 1) % 3 == 0) count++;
		if ((y - 1) % 3 == 0) count++;
		if ((z - 1) % 3 == 0) count++;

		if (count >= 2)
		{
			mat_id = 0;
			break;
		}

		x /= 3;
		y /= 3;
		z /= 3;
	}

	return mat_id;
}

float f(float n, float h, float s, float v) {
	float k = fmod(n + (h / 60.0f), 6.0f);
	return v - v * s * std::max(0.0f, std::min({ k, 4.0f - k, 1.0f }));
}

// All parameters are in the range[0.0, 1.0].
std::tuple<float, float, float> hsv_to_rgb( float h, float s, float v ) {
	h = h * 360.0; // To degrees
	float r = f(5.0, h, s, v);
	float g = f(3.0, h, s, v);
	float b = f(1.0, h, s, v);
	return std::make_tuple(r, g, b);
}

bool generateVolume(Algorithm algorithm,
	const std::filesystem::path& output_path,
	int size)
{
	// Fractal noise
	int octaves = log2(size);

	// Worley noise
	int cell_size = 100;
	int border = 2;

	std::function<u8(int, int, int)> func;
	switch (algorithm)
	{
	case Algorithm::fractal_noise:
		func = [=](int x, int y, int z) {
			if (fractal_noise(x, y, z, octaves) > 0.0f) {
				// Tempory hack to combine fractal and worley noise,
				// until we add CSG-like combining of volumes.
				return worley_noise(x, y, z, cell_size, border);
			} else {
				return u8(0);
			}
		};
		break;
	case Algorithm::menger_sponge:
		func = mengerSponge;
		break;
	case Algorithm::worley_noise:
		func = [=](int x, int y, int z) {
			return worley_noise(x, y, z, cell_size, border);
		};
		break;
	}

	Cubiquity::Timer timer;
	Volume volume;
	Metadata metadata;
	metadata.materials.push_back(Material::EmptySpace);
	metadata.materials.push_back(Material("CellBorder", vec3(0.1f, 0.1f, 0.1f)));

	int colour_count = 254;
	float h = 0.0;
	float s = 1.0;
	float v = 0.8;

	for (int i = 0; i < colour_count; i++) {
		Material material;
		auto [r, g, b] = hsv_to_rgb(h, s, v);
		material.set_base_color(vec3(r, g, b));
		metadata.materials.push_back(material);
		h += 1.0f / colour_count;
	}

	while (metadata.materials.size() < 256) {
		metadata.materials.push_back(Material("LightGrey", vec3(0.8f, 0.8f, 0.8f)));
	}

	Cubiquity::MaterialId matId = 1;

	std::vector<int> z_vals(size);
	std::iota(z_vals.begin(), z_vals.end(), 0);

	//for (int z = 0; z < size; z++)
	std::mutex m;
	std::for_each(std::execution::par, z_vals.begin(), z_vals.end(), [&](int z) {

		auto slice = std::make_unique<u8[]>(size * size);
		// FIXME - We need to think how do drive a proper progress bar from a parallel for loop.
		log_info("{} of {}", z + 1, size);
		for (int y = 0; y < size; y++)
		{
			for (int x = 0; x < size; x++)
			{
				slice.get()[y * size + x] = func(x, y, z);
			}
		}

		std::scoped_lock lock{m};

		for (int y = 0; y < size; y++)
		{
			for (int x = 0; x < size; x++)
			{
				volume.setVoxel(x, y, z, slice.get()[y * size + x]);
			}
		}
	});
	log_info("Generated in {} seconds", timer.elapsedTimeInSeconds());

	saveVolume(output_path, volume, metadata);

	return true;
}
