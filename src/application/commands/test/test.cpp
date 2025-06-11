#include "test.h"

#include "test_base.h"
#include "test_linear_algebra.h"
#include "test_rendering.h"
#include "test_volume.h"

#include "base/logging.h"

bool test(const Test& tests)
{
	if ((tests == Test::all) || (tests == Test::base)) {
		testBase();
	}

	//testLinearAlgebra();

	//testVoxelization();


	if ((tests == Test::all) || (tests == Test::volume)) {
		testVolume();
	}

	if ((tests == Test::all) || (tests == Test::raytracing)) {
		if (!testRaytracingBehaviour())
		{
			log_error("TEST FAILED!");
		}

		if (!testRaytracingPerformance())
		{
			log_error("TEST FAILED!");
		}
	}

	return true;
}
