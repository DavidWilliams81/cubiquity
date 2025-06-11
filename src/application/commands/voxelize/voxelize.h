#ifndef CUBIQUITY_VOXELIZE_H
#define CUBIQUITY_VOXELIZE_H

#include <filesystem>
#include <optional>

static const int VoxelizeDefaultSize = 500;

// Note: Scale and size are mutually exclusive. It is an error to provide both.
bool voxelize(const std::filesystem::path& in_path,
	          const std::filesystem::path& out_path,
	          std::optional<float> scale, std::optional<int> size);

#endif // CUBIQUITY_VOXELIZE_H
