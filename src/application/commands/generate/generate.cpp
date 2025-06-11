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
#include "base/types.h"

#include "storage.h"

#include <filesystem>

using Cubiquity::Volume;
using Cubiquity::uint;

using namespace std;

bool mengerSponge(int x, int y, int z)
{
	int levels = 8;
	bool occupied = true;

	for (int i = 0; i < levels; i++)
	{
		int count = 0;
		if ((x - 1) % 3 == 0) count++;
		if ((y - 1) % 3 == 0) count++;
		if ((z - 1) % 3 == 0) count++;

		if (count >= 2)
		{
			occupied = false;
			break;
		}

		x /= 3;
		y /= 3;
		z /= 3;
	}

	return occupied;
}

bool generateVolume(const std::filesystem::path& output_path, uint size_exp)
{
	Volume volume;
	Cubiquity::MaterialId matId = 1;
	Metadata metadata;
	metadata.materials.push_back(Metadata::EmptySpace);
	metadata.materials.push_back(Metadata::Default);

	u32 size = 1;
	for (uint i = 0; i < size_exp; i++)
	{
		size *= 3;
	}

	for (int z = 0; z < size; z++)
	{
		log_info("{} of {}", z+1, size);
		for (int y = 0; y < size; y++)
		{
			for (int x = 0; x < size; x++)
			{
				bool occupied = mengerSponge(x, y, z);
				volume.setVoxel(x, y, z, occupied ? matId : 0);
			}
		}
	}	

	// Save the result
	log_info("Saving volume as '{}'", output_path);
	volume.save(output_path.string());
	saveMetadataForVolume(metadata, output_path);

	return true;
}

