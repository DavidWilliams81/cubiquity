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
#include "geometry.h"
#include "storage.h"
#include "voxelization.h"

#include <math.h>

#include "utility.h"
#include "geometry.h"
#include "storage.h"

#include <iostream>
#include <limits>
#include <list>
#include <queue>
#include <thread>
#include <mutex>
#include <unordered_set>

using namespace Cubiquity;
using namespace Internals;

int voxelizeTerrain(const std::string& heightMap, const std::string& colourMap);

int main(int /*argc*/, char** /*argv*/)
{
	//voxelizeTerrain("", "");
	//return 0;

	// Create a volume
	int targetSize = 256;
	Volume volume;

	// Load a mesh to voxelise
	const char* path = "../data";
	const char* filename = "shapes.obj";

	auto splitTriangles = loadObjFile(path, filename);

	std::cout << std::endl;
	std::cout << "Loaded " << filename << ":" << std::endl;
	std::cout << std::endl;

	// Flip mesh for Quake test.
	/*for (auto& vertex : mesh.mVertices)
	{
		std::swap(vertex.data[1], vertex.data[2]);

		vertex.data[2] = -vertex.data[2];
	}*/
	//mesh.computeBounds(); // As we have modified the vertices since loading.

	//Box3f bounds = computeBounds(triangles);

	Box3f bounds = computeBounds(splitTriangles);

	// Scale the mesh such that it fills volume.
	Vector3f minBound = bounds.lower();
	Vector3f maxBound = bounds.upper();

	// Move to the origin
	// FIXME - It seems that removing this translation changes the number of inside voxels by a suprisingly large amount
	// (approx 10% for the e1m1 map at 256 resolution). But why? I think it probably looks worse without the translation
	// but I'm not sure. Why does it have an effect? Because the translation lines up the value with integer boundaries
	// perhaps? Or floating point errors when further from the origin? Needs more investigation.
	//translate(splitTriangles, minBound * -1.0f);

	// The same scale factor is applied to all axes to avoid distorting the model.
	Vector3f dims = maxBound - minBound;
	float scaleFactor = (targetSize - 5) / dims.maxComponentValue(); // Couple of voxels of border.

	scale(splitTriangles, scaleFactor);
	//translate(splitTriangles, Vector3f(2.0f));

	// Perform the voxelization
    Timer timer;
	
	//insideOutsideVoxelize(volume, splitTriangles);
	TerminalProgressBar progressBar;
	voxelize(volume, splitTriangles, true, 0, &progressBar);
	std::cout << "Voxelised in " << timer.elapsedTimeInSeconds() << " seconds" << std::endl;
	std::cout << "Node count before merging = " << volume.countNodes() << std::endl;

	// Save the result
	std::cout << "Saving volume...";
	volume.save("../data/voxelized.vol");
	std::cout << "done." << std::endl;

	// Write out images for visualization.
	//saveVolumeAsImages(volume, "..", &progressBar);

	auto result = estimateBounds(volume);
	Box3i estimatedBounds = result.second;
	std::cout << estimatedBounds.lower() << " " << estimatedBounds.upper() << std::endl;
	Histogram histogram = computeHistogram(volume, estimatedBounds);
	printHistogram(histogram);
}

int voxelizeTerrain(const std::string& heightMap, const std::string& colourMap)
{
	Timer timer;
	Volume volume;

	uint32 width = 256;
	uint32 height = 256;

	bool doneTidy = false;

	for (int x = 0; x < width; x++)
	{
		std::cout << "x = " << x << std::endl;
		for (int y = 0; y < height; y++)
		{
			//std::cout << "y = " << y << std::endl;
			for (int z = 0; z < 256; z++)
			{
				//std::cout << "z = " << z << std::endl;
				volume.setVoxel(x, y, z, 0x0FFF);
				if (z == 100 && !doneTidy)
				{
					//volume.tidyHashMap();
					doneTidy = true;
				}
			}
		}
		//volume.tidyEdits();
	}

	volume.save("../data/terrain.vol");
	std::cout << "Finished in " << timer.elapsedTimeInSeconds() << " seconds" << std::endl;

	return 0;
}
