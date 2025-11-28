#ifndef CUBIQUITY_VIEW_H
#define CUBIQUITY_VIEW_H

// Sanity check that build system is setting things up correctly.
#ifndef CUBIQUITY_APP_ENABLE_VIEW
#error "CUBIQUITY_APP_ENABLE_VIEW should be set when building view command"
#endif

#include <filesystem>
#include <string>

enum class ViewMode {
	cpu_pathtracing,
	gpu_pathtracing,
	instancing
};

bool viewVolume(ViewMode mode, const std::filesystem::path& input_path,
	            int width, int height, int duration);

#endif // CUBIQUITY_VIEW_H

