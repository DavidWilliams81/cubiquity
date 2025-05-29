#include "fractal_noise.h"

extern "C"
{
	#include "simplexnoise1234.h"
}

#include "storage.h"

#include <cfloat>
#include <climits>

using Cubiquity::MaterialId;

// Based on 'stb_perlin_fbm_noise3' (public domain).
float fractal_noise(int x, int y, int z, int octaves)
{
	const float lacunarity = 0.5;
	const float gain = 2.0;

	float frequency = 1.0f;
	float amplitude = 1.0f;

	float sum = 0.0f;
	// High frequency components first
	for (int i = 0; i < octaves; i++)
	{
		sum += snoise3(x*frequency, y*frequency, z*frequency)*amplitude;
		frequency *= lacunarity;
		amplitude *= gain;
	}

	// Not normalised, but we (usually) only care about the sign
	return sum;
}

ivec3 centre(ivec3 cell, int cell_size)
{
	// Pick a random position within the cell to become its centre
	std::minstd_rand rng(std::hash<ivec3>{}(cell));
	ivec3 local_pos = ivec3(rng(), rng(), rng()) % cell_size;
	ivec3 global_pos = (cell * cell_size) + local_pos;
	return global_pos;
}

u8 material(ivec3 cell)
{
	// Pick a random material for the cell (material 1 reserved for border)
	return std::hash<ivec3>{}(cell) % 254 + 2; // Range 2 - 255
}

// Worley noise is also known as Voronoi Noise.
// See https://www.ronja-tutorials.com/post/028-voronoi-noise/
u8 worley_noise(int x, int y, int z, int cell_size, int border)
{
	// Current position and cell
	const ivec3 sample_pos({ x, y, z });
	const ivec3 sample_cell = sample_pos / cell_size;

	// Best values found so far
	int closest_cell_dist_sq = INT_MAX;
	MaterialId closest_cell_material = 0;
	ivec3 closest_cell_centre;

	// Iterate over surrounding cells
	for (int nz = -1; nz <= 1; nz++) {
		for (int ny = -1; ny <= 1; ny++) {
			for (int nx = -1; nx <= 1; nx++) {

				// Determine the cell's centre
				ivec3 cell = sample_cell + ivec3({ nx, ny, nz });
				ivec3 cell_centre = centre(cell, cell_size);

				// Find distance from sample point to centre of cell
				ivec3 sample_to_centre = cell_centre - sample_pos;
				int cell_dist_sq = dot(sample_to_centre, sample_to_centre);

				// Use this cell if it is the closest.
				if (cell_dist_sq < closest_cell_dist_sq) {
					closest_cell_dist_sq = cell_dist_sq;
					closest_cell_material = material(cell);
					closest_cell_centre = cell_centre;
				}
			}
		}
	}

	float min_dist_to_edge = FLT_MAX;

	// Iterate over surrounding cells again
	for (int nz = -1; nz <= 1; nz++) {
		for (int ny = -1; ny <= 1; ny++) {
			for (int nx = -1; nx <= 1; nx++) {

				// Cell centre is the same as in previous pass above
				ivec3 cell = sample_cell + ivec3({ nx, ny, nz });
				ivec3 cell_centre = centre(cell, cell_size);

				// Skip comparing cell with itself
				if (cell_centre == closest_cell_centre) continue;

				// Distance calculation most easily done with floats
				vec3 midpoint =
					(vec3(cell_centre) + vec3(closest_cell_centre)) * 0.5f;
				vec3 closest_to_current =
					normalize(vec3(cell_centre) - vec3(closest_cell_centre));
				vec3 sample_pos_to_midpoint = midpoint - vec3(sample_pos);
				float dist_to_edge =
					dot(closest_to_current, sample_pos_to_midpoint);

				// Use edge from this cell if it is the closest.
				min_dist_to_edge = std::min(dist_to_edge, min_dist_to_edge);
			}
		}
	}

	// Apply border if wanted
	u8 border_material = 1;
	float half_border = static_cast<float>(border) * 0.5f;
	return min_dist_to_edge < half_border ?
		border_material : closest_cell_material;
}
