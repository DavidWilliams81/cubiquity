#include "view.h"
#include "gpu_pathtracing_viewer.h"
#include "instancing_demo.h"
#include "pathtracing_demo.h"

#include "base/logging.h"
#include "base/paths.h"

bool viewVolume(const flags::args& args)
{
	std::filesystem::path inputPath(args.positional().at(1));
	if (!checkInputFileIsValid(inputPath)) return false;

	const auto mode = args.get<std::string>("mode", "gpu-pathtracing");
	if (mode == "instancing")
	{
		InstancingDemo app(inputPath.string());
		app.show(1600, 1200);
	}
	else if (mode == "gpu-pathtracing")
	{
		GPUPathtracingViewer gpuPathtracingViewer(inputPath.string());
		gpuPathtracingViewer.show(1200, 900);
	}
	else if (mode == "pathtracing")
	{
		PathtracingDemo pathtracingDemo(inputPath.string());
		pathtracingDemo.show(800, 600);
	}
	else
	{
		log_error("Unrecognised value '{}' for --mode", mode);
	}
	return true;
}
