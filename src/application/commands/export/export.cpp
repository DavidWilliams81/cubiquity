#include "export.h"

#include "base/paths.h"
#include "base/logging.h"
#include "base/progress.h"

#include "stb_image_write.h"

using namespace Cubiquity;

Colour rgbFromMaterialId(MaterialId matId)
{
	uint8 red = (matId >> 8) & 0xF;
	uint8 green = (matId >> 4) & 0xF;
	uint8 blue = (matId) & 0xF;

	return Colour(red * 16.0f, green * 16.0f, blue * 16.0f);
}

void saveVolumeAsImages(Volume& volume, const std::string& dirName)
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
				Colour colour = rgbFromMaterialId(matId);
				colour.alpha = matId > 0 ? 255 : 0;

				imageData.push_back(colour.red);
				imageData.push_back(colour.green);
				imageData.push_back(colour.blue);
				imageData.push_back(colour.alpha);
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
	saveVolumeAsImages(volume, ".");
	return true;
}