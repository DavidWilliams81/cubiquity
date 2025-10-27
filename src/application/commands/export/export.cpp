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
#include <iostream>
#include <fstream>
#include <memory>
#include <filesystem>

// Hack for testing example code from main project
/*#define main run_vox_writer_example
#include "vox_writer/example.cpp"
#undef main*/

using Cubiquity::Volume;

const std::string bin_comment =
R"(# Overview
# --------
# This text file stores volume metadata such as dimensions and materials, using
# the TOML markup format (see https://toml.io for specs and parsers). The actual
# voxels are stored in a separate binary file with a '.bin' extension:
#
# The format of a .bin file is very simple. The voxels are stored as a raw 3D
# array of 8-bit unsigned integers. The array dimensions are given below. The
# array is laid out with the 'x' component changing most quickly and the 'z'
# component changing most slowly. Each voxel represents an index (0-255) into
# the material array (also given below) describing the voxel.
)";

void write_volume_as_bin(Volume* volume, // FIXME - Should ideally be ref not ptr
	                     const std::filesystem::path& output_path)
{
	OutputHandle file(output_path.string());

	auto [lower_bound, upper_bound] = find_bounds(*volume);

	for (int z = lower_bound.z; z <= upper_bound.z; z++)
	{
		for (int y = lower_bound.y; y <= upper_bound.y; y++)
		{
			for (int x = lower_bound.x; x <= upper_bound.x; x++)
			{
				Cubiquity::MaterialId matId = volume->voxel(x, y, z);
				file->write(reinterpret_cast<const char*>(&matId), sizeof(matId));
			}
		}

		// A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar.
		cubiquityProgressHandler("Saving volume as raw 3D array",
			lower_bound.z, z, upper_bound.z);
	}
}

void write_metadata_as_txt(Volume* volume,  // FIXME - Should ideally be ref not ptr
                           const Metadata& metadata,
	                       const std::filesystem::path& metadata_path)
{
	OutputHandle metadata_file(metadata_path.string());
	Metadata bin_metadata = metadata;
	bin_metadata.header = bin_comment;
	bin_metadata.dimensions = find_dimensions(*volume);
	bin_metadata.save(metadata_file.get());
}

void export_as_bin(Volume& volume, const Metadata& metadata,
	               const std::filesystem::path& output_path,
	               const std::filesystem::path& metadata_path)
{
	// We write the files from background threads in case we are outputting to
	// named pipes. In this case the writing may block until there is an active
	// reader, and we don't control in which order an external reader will read
	// the files. If we get it wrong we'll have a deadlock, so it is safest to
	// write both files in parallel.
	std::thread t1(write_volume_as_bin, &volume, output_path);
	std::thread t2(write_metadata_as_txt, &volume, metadata, metadata_path);

	t1.join();
	t2.join();
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
               std::filesystem::path  output_path,
               std::filesystem::path  output_metadata_path)
{
	auto [volume, metadata] = loadVolume(input_path.string());

	switch (format)
	{
	case ExportFormat::bin:
		if(output_path == "-") {
			if (output_metadata_path.empty()) {
				throw std::runtime_error(
					"Output metadata path must be specified when writing to stdout");
			}
		} else {
			if (output_path.empty()) {
				output_path = input_path.filename(); // Output to current dir
				output_path.replace_extension(".bin");
			}
			if (output_metadata_path.empty()) {
				// We use the .txt extension for metadata associated with a .bin file partly
				// because it avoids a name conflict when importing/exporting a .bin to/from
				// a .dag, and partly to make it obvious that the metadata for a .bin is
				// human-readable (I'm not sure everyone knows what a .toml extension means).
				output_metadata_path = output_path;
				output_metadata_path.replace_extension(".txt");
			}
		}
		export_as_bin(*volume, metadata, output_path, output_metadata_path);
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
