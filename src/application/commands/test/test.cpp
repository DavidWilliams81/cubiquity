#include "test.h"

#include "test_base.h"
#include "test_linear_algebra.h"
#include "test_rendering.h"
#include "test_visibility.h"
#include "test_volume.h"
#include "test_voxelisation.h"

#include <iostream>

bool test(const flags::args& args)
{
	const std::string feature{ args.positional().at(0) };
	//testBase();

	//testLinearAlgebra();

	/*if (!testRasterization())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}*/

	if (!testRaytracingBehaviour())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}

	if (!testRaytracingPerformance())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}

	/*if (!testVisibility())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}

	testVoxelization();

	testVolume();*/

	//testMerge();

	return true;
}