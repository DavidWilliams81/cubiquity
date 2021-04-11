#include "test_base.h"
#include "test_linear_algebra.h"
#include "test_rendering.h"
#include "test_utility.h"
#include "test_visibility.h"
#include "test_volume.h"
#include "test_voxelisation.h"

#include <iostream>

int main(int argc, char *argv[])
{
	//testBase();

	//testUtility();

	//testLinearAlgebra();

	/*if (!testRasterization())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}*/

	/*if (!testRaytracingBehaviour())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}

	if (!testRaytracingPerformance())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}*/

	if (!testVisibility())
	{
		std::cout << "TEST FAILED!" << std::endl;
	}

	//testVoxelization();

	testVolume();

	//testMerge();
}
