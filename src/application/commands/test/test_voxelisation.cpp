#include "test_voxelisation.h"

#include "framework.h"

#include "geometry.h"
#include "voxelization.h"

#include <math.h>
#include <random>
#include <sstream>

using namespace Cubiquity;
using namespace std;

bool testWindingNumberBehaviour(const TriangleList& triangles, const Box3f& bounds)
{
	/*std::cout << "Testing winding number behaviour" << std::endl;
	std::cout << "--------------------------------" << std::endl;

	uint32 internalSamplesBruteForce = 0;
	uint32 internalSamplesHierarchical = 0;

	double bruteForceVsHierarchicalError = 0.0;

	Patch patchNodeNormal(triangles);

	Box3fSampler Box3fSampler(bounds);
	const uint32 samples = 10000;
	for (int i = 0; i < samples; i++)
	{
		auto v = Box3fSampler.next();

		double bruteForce = computeWindingNumber(v, triangles);
		double hierarchicalNormal = computeWindingNumber(v, patchNodeNormal);

		if (bruteForce > 0.5) { internalSamplesBruteForce++; }
		if (hierarchicalNormal > 0.5) { internalSamplesHierarchical++; }

		bruteForceVsHierarchicalError += std::abs(bruteForce - hierarchicalNormal);
	}

	bruteForceVsHierarchicalError /= samples;

	assert(internalSamplesBruteForce         == 4157);
	assert(internalSamplesHierarchicalNormal == 4157);

	std::cout << "Internal samples count (brute force) =          "
		      << internalSamplesBruteForce << std::endl;
	std::cout << "Internal samples count (hierarchical) = "
		      << internalSamplesHierarchical << std::endl;
	std::cout << std::endl;
	std::cout << "Brute-force vs. hierarchical average error = "
			  << bruteForceVsHierarchicalError << std::endl;

			  */
	return true;
}

bool testWindingNumberPerformance(const TriangleList& triangles, const Box3f& bounds)
{
	/*
	std::cout << "Testing winding number performance" << std::endl;
	std::cout << "----------------------------------" << std::endl;

	Patch patchNode(triangles);

	Box3fSampler Box3fSampler(bounds);

	uint32_t insideCount = 0;
	uint32_t outsideCount = 0;

	Timer timer;

	const uint32 samples = 100000;
	for (int i = 0; i < samples; i++)
	{
		auto v = Box3fSampler.next();
		computeWindingNumber(v, patchNode) > 0.5f ? insideCount++ : outsideCount++;
	}

	std::cout << "Interior voxels = " << insideCount  << " / " << samples << std::endl;
	std::cout << "Exterior voxels = " << outsideCount << " / " << samples << std::endl;
	std::cout << "Time elapsed    = " << timer.elapsedTimeInSeconds() << " seconds" << std::endl;

	check(insideCount,  41679);
	check(outsideCount, 58321);

	*/
	return true;
}

/*bool testVoxelize(Geometry& geometry, bool preserveSurfaceMaterials = false,
	uint16 internalMaterialOveride = 0, bool useBruteForce = false)
{
	std::cout << "Testing voxelisation" << std::endl;
	std::cout << "--------------------" << std::endl;

	Volume volume;
	TerminalProgressBar progressBar;

	for (auto& object : geometry)
	{
		std::cout << "Voxelising " << object.name << "..." << std::endl;
		voxelize(volume, object, preserveSurfaceMaterials, internalMaterialOveride, &progressBar, useBruteForce);
		std::cout << std::endl;
	}

	auto result = estimateBounds(volume);
	Box3i bounds = result.second;
	Histogram histogram = computeHistogram(volume, bounds);

	// Print details
	std::cout << "Lower bound = (" << bounds.lower().x() << ", " << bounds.lower().y() << ", " << bounds.lower().z() << ")" << std::endl;
	std::cout << "Upper bound = (" << bounds.upper().x() << ", " << bounds.upper().y() << ", " << bounds.upper().z() << ")" << std::endl;
	printHistogram(histogram);

	// Validate details
	check(bounds.lower(),  Vector3i(-99, -154, -99));
	check(bounds.upper(),  Vector3i(99, 161, 108));
	check(histogram[0],    7640947);
	if (!preserveSurfaceMaterials && !internalMaterialOveride)
	{
		check(histogram[15], 506923);
		check(histogram[240], 1639206);
		check(histogram[3840], 2753352);
		check(histogram[4080], 539444);
	}
	else if (preserveSurfaceMaterials && internalMaterialOveride)
	{
		check(histogram[15], 28300);
		check(histogram[240], 160822);
		check(histogram[3840], 182581);
		check(histogram[4080], 52437);
		check(histogram[internalMaterialOveride], 5014785);
	}
	else
	{
		std::cout << "No reference results for configuration" << std::endl;
	}

	// Save in case we need to inspect the results
	std::stringstream ss;
	ss << "../test-voxelize-" << preserveSurfaceMaterials <<
		"-" << internalMaterialOveride << "-" << useBruteForce << ".vol";
	volume.save(ss.str());
	saveVolumeAsImages(volume, "..", &progressBar);

	return true;
}*/

bool testVoxelization()
{
	/*const char* filename = "shapes.obj";
	auto geometry = loadObjFile("../../Data", filename);

	// The 'shapes' model is relatively small so we need to scale it up to be representitive
	// of a real voxelisation. This particularly matters when snapping the vertex positions 
	// to a sub-voxel grid (for improved robustness) as this operation assumes that the mesh
	// is large relative to the voxel size.
	scale(geometry, 100.0f);
	TriangleList triangles = mergedTriangles(geometry);
	Box3f bounds = computeBounds(triangles);

	//testWindingNumberBehaviour(triangles, bounds);
	std::cout << std::endl;
	//testWindingNumberPerformance(triangles, bounds);
	std::cout << std::endl;

	//testVoxelize(geometry);
	testVoxelize(geometry, false, 0);
	std::cout << std::endl;
	testVoxelize(geometry, true, 0x0777);
	std::cout << std::endl;
	testVoxelize(geometry, false, 0, true); // Brute force, slow.
	//testVoxelizeWithPreserveSurface(geometry);
	*/
	return true;
}