#include "export.h"

#include "base/paths.h"
#include "base/logging.h"
#include "base/progress.h"

#include "cubiquity.h"
#include "storage.h"

#include "volume_vox_writer.h"

#include "stb_image_write.h"

#include <algorithm>
#include <fstream>

// Hack for testing example code from main project
/*#define main run_vox_writer_example
#include "vox_writer/example.cpp"
#undef main*/

using namespace Cubiquity;

void saveVolumeAsImages(Volume& volume, const Metadata& metadata, const std::string& dirName)
{
	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	// Expand the bounds in case we have a scene with a solid exterior, as in
	// this case it is useful to see some of the exterior in the image slices.
	const int border = 5;
	lower_x -= border; lower_y -= border; lower_z -= border;
	upper_x += border; upper_y += border; upper_z += border;

	for (int z = lower_z; z <= upper_z; z += 1)
	{
		// Note that the filenames start at zero (they are never negative). Using +/- symbols in the filenames is problematic,
		// at least because when sorting by name the OS lists '+' before'-', and also larger-magnitiude negative number after
		// smaller-magnitude negative numbers. This makes it more difficult to scroll through the slices.
		char filepath[256];
		std::snprintf(filepath, sizeof(filepath), "%s/%06d.png", dirName.c_str(), z - lower_z);

		//Image image(width, height);
		std::vector<uint8> imageData;
		for (int y = lower_y; y <= upper_y; y++)
		{
			for (int x = lower_x; x <= upper_x; x++)
			{
				MaterialId matId = volume.voxel(x, y, z);

				Col diffuse = metadata.materials.at(matId).diffuse;

				float gamma = 1.0f / 2.2f;
				diffuse[0] = pow(diffuse[0], gamma);
				diffuse[1] = pow(diffuse[1], gamma);
				diffuse[2] = pow(diffuse[2], gamma);

				uint8 r = std::clamp(std::lround(diffuse[0] * 255.0f), 0L, 255L);
				uint8 g = std::clamp(std::lround(diffuse[1] * 255.0f), 0L, 255L);
				uint8 b = std::clamp(std::lround(diffuse[2] * 255.0f), 0L, 255L);

				imageData.push_back(r);
				imageData.push_back(g);
				imageData.push_back(b);
				imageData.push_back(matId > 0 ? 255 : 0);
			}
		}

		int width  = (upper_x - lower_x) + 1;
		int height = (upper_y - lower_y) + 1;
		int result = stbi_write_png(filepath, width, height, 4, imageData.data(), width * 4);
		if (result == 0)
		{
			log(Error, "Failed to write PNG image.");
		}

		// A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar. 
		cubiquityProgressHandler("Saving volume as images", lower_z, z, upper_z);
	}
}

void saveVolumeAsVox(Volume& volume, const Metadata& metadata, const std::string& fileName)
{	
	// Hack for testing example code from main project
	/*run_vox_writer_example();
	return;*/

	Timer timer;
	try {		
		volume_vox_writer writer(volume, metadata);
		writer.write("output.vox", false);
		std::cout << "Exported .vox in " << timer.elapsedTimeInSeconds() << " seconds." << std::endl;
	} catch (std::exception& e) {;
		std::cerr << "Failed to write .vox file (" << e.what() << ")." << std::endl;
	}
}

bool exportVolume(const flags::args& args)
{
	std::filesystem::path inputPath(args.positional().at(1));
	if (!checkInputFileIsValid(inputPath)) return false;

	Volume volume(inputPath.string());
	Metadata metadata = loadMetadataForVolume(inputPath);

	const auto outputPath = args.get<std::filesystem::path >("output_path", ".");
	if(outputPath.extension() == ".vox") {
		saveVolumeAsVox(volume, metadata, "output.vox");
	} else {
		if (!checkOutputDirIsValid(outputPath)) return false;
		saveVolumeAsImages(volume, metadata, ".");
	}

	return true;
}
