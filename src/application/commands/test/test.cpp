#include "test.h"

#include "test_base.h"
#include "test_linear_algebra.h"
#include "test_rendering.h"
#include "test_visibility.h"
#include "test_volume.h"
#include "test_voxelisation.h"

#include "base/logging.h"

bool test(const flags::args& args)
{
	const std::string feature{ args.positional().at(0) };
	//testBase();

	//testLinearAlgebra();

	/*if (!testRasterization())
	{
		log_error("TEST FAILED!");
	}*/

	if (!testRaytracingBehaviour())
	{
		log_error("TEST FAILED!");
	}

	if (!testRaytracingPerformance())
	{
		log_error("TEST FAILED!");
	}

	if (!testVisibility())
	{
		log_error("TEST FAILED!");
	}

	//testVoxelization();

	testVolume();

	//testMerge();

	return true;
}
