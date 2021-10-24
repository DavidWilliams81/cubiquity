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

#include "mesh.h"
#include "voxelize.h"

#include "base/logging.h"

#include "geometry.h"
#include "storage.h"
#include "voxelization.h"
#include "utility.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <list>
#include <cmath>
#include <queue>
#include <thread>
#include <mutex>
#include <unordered_set>

using namespace Cubiquity;
using namespace Internals;

using namespace std;

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

bool voxelizeMesh(const flags::args& args)
{
	const auto inputPath{ args.positional().at(1) };
	const auto outputPath = args.get<std::filesystem::path>("output", "output.vol");
	const auto scaleFactor = args.get<float>("scale", 32.0);

	if (!std::filesystem::exists(inputPath))
	{
		cerr << "Path \'" << inputPath << "\' does not exist!" << endl;
		return false;
	}
	// Load a mesh to voxelise
	auto objects = loadObjFile(inputPath);

	// Perform the voxelization
    Timer timer;
	Volume volume;
	MaterialSet materials;
	uint objectIndex = 0;
	for(auto& object : objects)
	{
		objectIndex++;
		std::cout << "Processing '" << object.name << "' ("
			<< objectIndex << " of " << objects.size() << ")..." << std::endl;

		if (object.triangles.empty()) { continue; } // Often happens for default object.

		Surface surface;
		surface.name = object.name;

		for (auto& tri : object.triangles)
		{
			Vector3f vertices[3];
			vertices[0] = Vector3f(tri.vertices[0][0], tri.vertices[0][1], tri.vertices[0][2]);
			vertices[1] = Vector3f(tri.vertices[1][0], tri.vertices[1][1], tri.vertices[1][2]);
			vertices[2] = Vector3f(tri.vertices[2][0], tri.vertices[2][1], tri.vertices[2][2]);
			Triangle triangle(vertices[0], vertices[1], vertices[2]);

			MaterialId matId = materials.findOrInsert(tri.colour);
			if (!matId)
			{
				log(Error, "Failed to find or insert material");
				matId = Cubiquity::Internals::MaxMaterial;
			}

			surface.addTriangle(triangle, matId);
		}	
		scale(surface.triangles, scaleFactor);
		surface.build();

		std::cout << "Voxelising " << surface.name << "..." << std::endl;
		Volume temp;
		MaterialId* internalMaterialOverride = nullptr;
		MaterialId* externalMaterialOverride = nullptr;
		voxelize(temp, surface, true, internalMaterialOverride, externalMaterialOverride);
		volume.addVolume(temp);
		volume.bake();

		std::cout << std::endl;
	}
	std::cout << "Voxelised in " << timer.elapsedTimeInSeconds() << " seconds" << std::endl;
	std::cout << "Node count before merging = " << volume.countNodes() << std::endl;

	// Save the result
	std::cout << "Saving volume as \'" << outputPath << "\'...";
	volume.save(outputPath.string());
	materials.save(getMaterialsPath(outputPath));
	std::cout << "done." << std::endl;	

	auto result = estimateBounds(volume);
	Box3i estimatedBounds = result.second;
	std::cout << estimatedBounds.lower() << " " << estimatedBounds.upper() << std::endl;
	Histogram histogram = computeHistogram(volume, estimatedBounds);
	printHistogram(histogram);

	return true;
}

