#ifndef CUBIQUITY_EXPORT_H
#define CUBIQUITY_EXPORT_H

#include <filesystem>

enum class ExportFormat {
    bin,  // Raw 3D array
    pngs, // PNG slices
	vox,  // MagicaVoxel
};

bool export_as(ExportFormat           format,
         const std::filesystem::path& input_path,
               std::filesystem::path  output_path);

#endif // CUBIQUITY_EXPORT_H
