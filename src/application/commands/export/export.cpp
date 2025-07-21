#include "export.h"

#include "base/bounds.h"
#include "base/logging.h"
#include "base/progress.h"
#include "base/serialize.h"

#include "cubiquity.h"
#include "storage.h"
#include "utility.h"

#include "volume_vox_writer.h"

#include "stb_image_write.h"

#include <algorithm>
#include <fstream>

// Hack for testing example code from main project
/*#define main run_vox_writer_example
#include "vox_writer/example.cpp"
#undef main*/

using Cubiquity::Volume;

void export_as_bin(Volume& volume, const Metadata& metadata,
                   const std::filesystem::path& output_path)
{
	std::ofstream file(output_path, std::ios::out | std::ios::binary);

	auto [lower_bound, upper_bound] = find_bounds(volume);

	for (int z = lower_bound.z; z <= upper_bound.z; z++)
	{
		for (int y = lower_bound.y; y <= upper_bound.y; y++)
		{
			for (int x = lower_bound.x; x <= upper_bound.x; x++)
			{
				Cubiquity::MaterialId matId = volume.voxel(x, y, z);
				file.write(reinterpret_cast<const char*>(&matId), sizeof(matId));
			}
		}

		// A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar. 
		cubiquityProgressHandler("Saving volume as raw 3D array",
			lower_bound.z, z, upper_bound.z);
	}
}

void export_as_images(Volume& volume, const Metadata& metadata,
	                  const std::filesystem::path& output_path)
{
	// Including a border gives useful context and is aesthetically pleasing.
	// It makes it clear that the whole object has been exported (not clipped),
	// and if the exterior is solid then some of that is captured too.
	// We only offer this border for image export (not other export types)
	// because this is intended for debug and visualisation, whereas other types
	// might be reimported and tampering with the bounds might complicate this.
	const int border = 5;
	auto [lower_bound, upper_bound] = find_bounds(volume);
	lower_bound -= ivec3(border);
	upper_bound += ivec3(border);

	for (int z = lower_bound.z; z <= upper_bound.z; z += 1)
	{
		// Note that the filenames start at zero (they are never negative). Using +/- symbols in the filenames is problematic,
		// at least because when sorting by name the OS lists '+' before'-', and also larger-magnitiude negative number after
		// smaller-magnitude negative numbers. This makes it more difficult to scroll through the slices.
		char filepath[256];
		std::snprintf(filepath, sizeof(filepath), "%s/%06d.png", output_path.string().c_str(), z - lower_bound.z);

		//Image image(width, height);
		std::vector<u8> imageData;
		for (int y = lower_bound.y; y <= upper_bound.y; y++)
		{
			for (int x = lower_bound.x; x <= upper_bound.x; x++)
			{
				Cubiquity::MaterialId matId = volume.voxel(x, y, z);

				vec3 base_color = metadata.materials[matId].base_color();

				float gamma = 1.0f / 2.2f;
				base_color[0] = pow(base_color[0], gamma);
				base_color[1] = pow(base_color[1], gamma);
				base_color[2] = pow(base_color[2], gamma);

				u8 r = std::clamp(std::lround(base_color[0] * 255.0f), 0L, 255L);
				u8 g = std::clamp(std::lround(base_color[1] * 255.0f), 0L, 255L);
				u8 b = std::clamp(std::lround(base_color[2] * 255.0f), 0L, 255L);

				imageData.push_back(r);
				imageData.push_back(g);
				imageData.push_back(b);
				imageData.push_back(matId > 0 ? 255 : 0);
			}
		}

		int width  = (upper_bound.x - lower_bound.x) + 1;
		int height = (upper_bound.y - lower_bound.y) + 1;
		int result = stbi_write_png(filepath, width, height, 4, imageData.data(), width * 4);
		if (result == 0)
		{
			log_error("Failed to write PNG image");
		}

		// A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar. 
		cubiquityProgressHandler("Saving volume as images",
		                         lower_bound.z, z, upper_bound.z);
	}
}

void export_as_vox(Volume& volume, const Metadata& metadata,
	               const std::filesystem::path& output_path)
{	
	// Hack for testing example code from main project
	/*run_vox_writer_example();
	return;*/

	Cubiquity::Timer timer;
	try {		
		volume_vox_writer writer(volume, metadata);
		writer.write(output_path.string(), false);
		log_info("Exported .vox in {} seconds", timer.elapsedTimeInSeconds());
	} catch (std::exception& e) {;
		log_error("Failed to write .vox file ({}).", e.what());
	}
}

bool export_as(ExportFormat           format,
         const std::filesystem::path& input_path,
               std::filesystem::path  output_path)
{
	auto [volume, metadata] = loadVolume(input_path.string());

	switch (format)
	{
	case ExportFormat::bin:
		if (output_path.empty()) {
			output_path = input_path.filename().replace_extension(".bin");
		}
		export_as_bin(*volume, metadata, output_path);
		break;
	case ExportFormat::pngs:
		if (output_path.empty()) {
			output_path = ".";
		}
		if (!checkOutputDirIsValid(output_path)) return false;
		export_as_images(*volume, metadata, output_path);
		break;
	case ExportFormat::vox:
		if (output_path.empty()) {
			output_path = input_path.filename().replace_extension(".vox");
		}
		export_as_vox(*volume, metadata, output_path);
		break;
	default:
		log_error("Unknown export format");
	}

	return true;
}
