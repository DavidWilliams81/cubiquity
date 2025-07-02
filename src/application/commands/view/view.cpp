#include "view.h"
#include "gpu_pathtracing_viewer.h"
#include "instancing_demo.h"
#include "pathtracing_demo.h"

#include "base/logging.h"
#include "base/serialize.h"

bool viewVolume(ViewMode mode, const std::filesystem::path& input_path)
{
	if (!checkInputFileIsValid(input_path)) return false;

	if (mode == ViewMode::instancing)
	{
		InstancingDemo app(input_path.string());
		app.show(1600, 1200);
	}
	else if (mode == ViewMode::cpu_pathtracing)
	{
		PathtracingDemo pathtracingDemo(input_path.string());
		pathtracingDemo.show(800, 600);
	}
	else if (mode == ViewMode::gpu_pathtracing)
	{
		GPUPathtracingViewer gpuPathtracingViewer(input_path.string());
		gpuPathtracingViewer.show(1200, 900);
	}

	return true;
}
