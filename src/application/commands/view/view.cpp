#include "view.h"
#include "gpu_pathtracing_viewer.h"
#include "instancing_demo.h"
#include "pathtracing_demo.h"

#include "base/logging.h"
#include "base/serialize.h"

bool viewVolume(ViewMode mode, const std::filesystem::path& input_path,
	            int width, int height, int duration)
{
	if (!checkInputFileIsValid(input_path)) return false;

	if (mode == ViewMode::instancing)
	{
		InstancingDemo app(input_path.string());
		app.show(width, height, duration);
	}
	else if (mode == ViewMode::cpu_pathtracing)
	{
		PathtracingDemo pathtracingDemo(input_path.string());
		pathtracingDemo.show(width, height, duration);
	}
	else if (mode == ViewMode::gpu_pathtracing)
	{
		GPUPathtracingViewer gpuPathtracingViewer(input_path.string());
		gpuPathtracingViewer.show(width, height, duration);
	}

	return true;
}
