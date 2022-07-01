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

#include "stb_image.h"

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

int imageWidth, imageHeight;
float* imageData;
float heightmapScale = 0.0f;

int32 myHeightFunc(int32 x, int32 y)
{
	if (x >= 0 && x < imageWidth && y >= 0 && y < imageHeight)
	{
		return static_cast<int32>(lroundf(imageData[x + y * imageWidth] * heightmapScale));
	}
	else
	{
		return 0;
	}
}

int voxelizeTerrain(const std::filesystem::path& inputPath, Volume& volume, MaterialSet& materials, float scaleFactor)
{
	heightmapScale = scaleFactor;

	// Prevent stb_image from applying a gamma mapping when loading an LDR image 
	// (which doesn't really make sense when the data represents a heightmap).
	stbi_ldr_to_hdr_gamma(1.0f);

	// Use the float interface for consistancy between 8-bit, 16bit and float images.
	// stb_image always gives data in the range 0.0-1.0 when using this interface.
	int channels;
	const int desiredChannels = 1;
	imageData = stbi_loadf(inputPath.string().c_str(), &imageWidth, &imageHeight, &channels, desiredChannels);
	if (channels > desiredChannels)
	{
		log(Warning, "Input heightmap contains multiple channels (convert it to a greyscale image without alpha)");
	}

	MaterialId surface = materials.findOrInsert({ 0.0, 0.5, 0.0 });
	MaterialId underground = materials.findOrInsert({ 0.75, 0.75, 0.75 });

	int32 minBounds[] = { 0, 0, 0 };
	int32 maxBounds[] = { imageWidth - 1, imageHeight - 1, static_cast<int32>(lroundf(heightmapScale)) };
	voxelizeHeightmap(volume, &myHeightFunc, minBounds, maxBounds, underground, surface);

	return 0;
}

bool voxelizeMesh(const std::filesystem::path& inputPath, Volume& volume, MaterialSet& materials, float scaleFactor)
{
	// Load a mesh to voxelise
	auto objects = loadObjFile(inputPath);
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
			vertices[0] = Vector3f({ tri.vertices[0][0], tri.vertices[0][1], tri.vertices[0][2] });
			vertices[1] = Vector3f({ tri.vertices[1][0], tri.vertices[1][1], tri.vertices[1][2] });
			vertices[2] = Vector3f({ tri.vertices[2][0], tri.vertices[2][1], tri.vertices[2][2] });
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

	return true;
}

bool voxelise(const flags::args& args)
{
	const std::filesystem::path inputPath(args.positional().at(1));
	const auto outputPath = args.get<std::filesystem::path>("output", "output.vol");
	const auto scaleFactor = args.get<float>("scale", 32.0);

	if (!std::filesystem::exists(inputPath))
	{
		cerr << "Path \'" << inputPath << "\' does not exist!" << endl;
		return false;
	}

	// Perform the voxelization
	Timer timer;
	Volume volume;
	MaterialSet materials;

	std::string extension = inputPath.extension().string();

	if(extension == ".obj")
	{ 
		voxelizeMesh(inputPath, volume, materials, scaleFactor);
	}
	else if ((extension == ".png") || (extension == ".jpg"))
	{
		voxelizeTerrain(inputPath, volume, materials, scaleFactor);
	}
	else
	{
		log(Error, "Unrecognised extension \'", extension, "\'");
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

	// This chunk of code can be useful for debuging and validation, but is slow.
	// Therefore only run it for mesh voxelisation as these are usually smaller.
	if (extension == ".obj")
	{
		Histogram histogram = computeHistogram(volume, estimatedBounds);
		printHistogram(histogram);
	}

	return true;
}
