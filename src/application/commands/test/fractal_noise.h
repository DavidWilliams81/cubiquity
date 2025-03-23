#ifndef CUBIQUITY_SIMPLEX_NOISE_H_
#define CUBIQUITY_SIMPLEX_NOISE_H_

#include "geometry.h"
#include "storage.h"

#include <cstdint>

// Fractal noise geometry with materials assigned based on Voronoi cells.
// Useful for producing detailed test volumes of potentially unlimited size.
// Note that the highest frequency component is always present. Adding more
// octaves increases the size of the largest features but the high-frequency
// detail always goes down to sub-voxel level. 
class FractalNoise
{
public:
	FractalNoise(int octaves, int offsetX = 0, int offsetY = 0, int offsetZ = 0);
	Cubiquity::MaterialId  operator()(int x, int y, int z);

private:

	float fractalNoise(int x, int y, int z);
	Cubiquity::MaterialId voronoiCell(int x, int y, int z);

	Cubiquity::vec3i chooseCentre(Cubiquity::vec3i cell);
	Cubiquity::MaterialId chooseMaterial(Cubiquity::vec3i cell);

	int mOctaves;
	int mOffsetX = 0;
	int mOffsetY = 0;
	int mOffsetZ = 0;

	const int mCellSize = 32;
};

#endif // CUBIQUITY_SIMPLEX_NOISE_H_