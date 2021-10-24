#include "export.h"

#include "base/paths.h"
#include "base/logging.h"
#include "base/materials.h"
#include "base/progress.h"

#include "storage.h"

#include "stb_image_write.h"

#include <algorithm>
#include <fstream>

using namespace Cubiquity;

void saveVolumeAsImages(Volume& volume, const MaterialSet& materials, const std::string& dirName)
{
	Box3i bounds = estimateBounds(volume).second;

	uint32_t width = bounds.sideLength(0);
	uint32_t height = bounds.sideLength(1);

	for (int z = bounds.lower().z(); z <= bounds.upper().z(); z += 1)
	{
		// Note that the filenames start at zero (they are never negative). Using +/- symbols in the filenames is problematic,
		// at least because when sorting by name the OS lists '+' before'-', and also larger-magnitiude negative number after
		// smaller-magnitude negative numbers. This makes it more difficult to scroll through the slices.
		char filepath[256];
		std::snprintf(filepath, sizeof(filepath), "%s/%06d.png", dirName.c_str(), z - bounds.lower().z());

		//Image image(width, height);
		std::vector<uint8> imageData;
		for (int y = bounds.lower().y(); y <= bounds.upper().y(); y++)
		{
			for (int x = bounds.lower().x(); x <= bounds.upper().x(); x++)
			{
				MaterialId matId = volume.voxel(x, y, z);

				uint8 r = std::clamp(std::lround(materials[matId][0] * 255.0f), 0L, 255L);
				uint8 g = std::clamp(std::lround(materials[matId][1] * 255.0f), 0L, 255L);
				uint8 b = std::clamp(std::lround(materials[matId][2] * 255.0f), 0L, 255L);

				imageData.push_back(r);
				imageData.push_back(g);
				imageData.push_back(b);
				imageData.push_back(matId > 0 ? 255 : 0);
			}
		}

		//image.save(filepath);
		int result = stbi_write_png(filepath, width, height, 4, imageData.data(), width * 4);
		if (result == 0)
		{
			log(Error, "Failed to write PNG image.");
		}

		// A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar. 
		cubiquityProgressHandler("Saving volume as images", bounds.lower().z(), z, bounds.upper().z());
	}
}

bool exportVolume(const flags::args& args)
{
	std::filesystem::path inputPath(args.positional().at(1));
	if (!checkInputFileIsValid(inputPath)) return false;

	const auto outputPath = args.get<std::string>("output_path", ".");
	if (!checkOutputDirIsValid(outputPath)) return false;

	Volume volume(inputPath.string());
	MaterialSet materials;
	materials.load(getMaterialsPath(inputPath));

	saveVolumeAsImages(volume, materials, ".");
	return true;
}
