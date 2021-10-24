#include "fractal_noise.h"

#include "geometry.h"
#include "storage.h"

using namespace Cubiquity;

extern "C"
{
	#include "simplexnoise1234.h"
}

FractalNoise::FractalNoise(int octaves, int offsetX, int offsetY, int offsetZ)
	: mOctaves(octaves), mOffsetX(offsetX), mOffsetY(offsetY), mOffsetZ(offsetZ) {}

MaterialId  FractalNoise::operator()(int x, int y, int z)
{
	MaterialId material = 0;

	if (fractalNoise(x, y, z) > 0.0f)
	{
		material = voronoiCell(x, y, z);
	}

	return material;
}

// Derived from stb_perlin.h, but modified to use Simplex noise instead. Note that the stb_perlin docs state:
//
//		Typical values to start playing with:
//			octaves    =   6     -- number of "octaves" of noise3() to sum
//			lacunarity = ~ 2.0   -- spacing between successive octaves (use exactly 2.0 for wrapping output)
//			gain       =   0.5   -- relative weighting applied to each successive octave
//
// However, the suggested values for gain and lacunarity appear to be the wrong way around.
// The code generates the base layer first (frequency and amplitude are 1.0) and each successive
// layer need to have lower frequency and bigger amplitude.
float FractalNoise::fractalNoise(int x, int y, int z)
{
	float lacunarity = 0.5;
	float gain = 2.0;

	int i;
	float frequency = 2.0f; // Possibly higher than needed, but catches the high-freq comonents.
	float amplitude = 1.0f;
	float sum = 0.0f;

	x += mOffsetX;
	y += mOffsetY;
	z += mOffsetZ;

	for (i = 0; i < mOctaves; i++)
	{
		sum += snoise3(x*frequency, y*frequency, z*frequency)*amplitude;
		frequency *= lacunarity;
		amplitude *= gain;
	}

	// The return value is not normalised, but this doesn't matter
	// as long as we only test whether it is above or below zero.
	return sum;
}

MaterialId FractalNoise::voronoiCell(int x, int y, int z)
{
	Vector3i pos(x, y, z);
	const Vector3i cell = pos / mCellSize;

	MaterialId closestMaterial = 0;
	int closestCellDistanceSquared = 1000000000;

	for (int nz = -1; nz <= 1; nz++)
	{
		for (int ny = -1; ny <= 1; ny++)
		{
			for (int nx = -1; nx <= 1; nx++)
			{
				Vector3i neighbourCellOffset(nx, ny, nz);
				Vector3i neighbourCell = cell + neighbourCellOffset;
				Vector3i neighbourCentre = chooseCentre(neighbourCell);
				Vector3i toNeighbour = neighbourCentre - pos;
				int neighbourDistSquared = dot(toNeighbour, toNeighbour);

				if (neighbourDistSquared < closestCellDistanceSquared)
				{
					closestCellDistanceSquared = neighbourDistSquared;
					closestMaterial = chooseMaterial(neighbourCell);
				}
			}
		}
	}

	return closestMaterial;
}

Vector3i FractalNoise::chooseCentre(Vector3i cell)
{
	// When hashing for centre we use a differnt seed than when hashing for material.
	uint32_t cellHash = Internals::murmurHash3(cell.data, sizeof(cell.data), 17);

	int32_t x = cellHash % mCellSize;
	cellHash = Internals::mix(cellHash);
	int32_t y = cellHash % mCellSize;
	cellHash = Internals::mix(cellHash);
	int32_t z = cellHash % mCellSize;

	return cell * mCellSize + Vector3i(x, y, z);
}

MaterialId FractalNoise::chooseMaterial(Vector3i cell)
{
	// When hashing for material we use a differnt seed than when hashing for centre.
	uint32_t cellHash = Internals::murmurHash3(cell.data, sizeof(cell.data), 65);

	const uint32_t limit = 50; // Increase this to get more cells which are clamped to max material.
	cellHash %= limit; // Constrain to range 0 - (limit-1)
	cellHash += 1;     // Constrain to range 1 - limit
	cellHash = std::min(cellHash, uint32_t(7)); // Cells between 7 and limit get set to 7.

	// Most cells get mat id of 7, with some getting a smaller id.
	MaterialId material = cellHash;
	return material;
}
