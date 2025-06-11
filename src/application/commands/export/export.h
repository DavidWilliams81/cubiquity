#ifndef CUBIQUITY_EXPORT_H
#define CUBIQUITY_EXPORT_H

#include <filesystem>

enum class ExportFormat {
	vox, // MagicaVoxel
	pngs // PNG slices
};

bool exportVolume(ExportFormat           format,
	        const std::filesystem::path& input_path,
	              std::filesystem::path  output_path);

#endif // CUBIQUITY_EXPORT_H
