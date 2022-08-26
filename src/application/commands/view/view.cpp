#include "view.h"
#include "instancing_demo.h"
#include "pathtracing_demo.h"

#include "base/paths.h"

bool viewVolume(const flags::args& args)
{
	std::filesystem::path inputPath(args.positional().at(1));
	if (!checkInputFileIsValid(inputPath)) return false;

	const auto mode = args.get<std::string>("mode", "instancing");
	if (mode == "instancing")
	{
		InstancingDemo app(inputPath.string());
		app.show(1600, 1200);
	}
	else if (mode == "pathtracing")
	{
		PathtracingDemo pathtracingDemo(inputPath.string());
		pathtracingDemo.show(1600, 1200);
	}
	else
	{
		std::cerr << "Unrecognised value \'" << mode << "'for --mode" << std::endl;
	}
	return true;
}