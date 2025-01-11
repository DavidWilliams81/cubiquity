#include "test_voxelisation.h"

#include "framework.h"

#include "geometry.h"
#include "voxelization.h"

#include <math.h>
#include <random>

using namespace Cubiquity;
using namespace std;

bool testWindingNumberBehaviour(const TriangleList& triangles, const Box3f& bounds)
{
	/*log_info("Testing winding number behaviour");
	log_info("--------------------------------");

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

	log_info("Internal samples count (brute force)  = {}", internalSamplesBruteForce);
	log_info("Internal samples count (hierarchical) = {}", internalSamplesHierarchical);
	log_info("");
	log_info("Brute-force vs. hierarchical average error = {}", bruteForceVsHierarchicalError);

			  */
	return true;
}

bool testWindingNumberPerformance(const TriangleList& triangles, const Box3f& bounds)
{
	/*
	log_info("Testing winding number performance");
	log_info("----------------------------------");

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

	log_info("Interior voxels = {} / {}", insideCount, samples);
	log_info("Exterior voxels = {} / {}", outsideCount, samples);
	log_info("Time elapsed    = {} seconds", timer.elapsedTimeInSeconds());

	check(insideCount,  41679);
	check(outsideCount, 58321);

	*/
	return true;
}

/*bool testVoxelize(Geometry& geometry, bool preserveSurfaceMaterials = false,
	uint16 internalMaterialOveride = 0, bool useBruteForce = false)
{
	log_info("Testing voxelisation");
	log_info("--------------------");

	Volume volume;
	TerminalProgressBar progressBar;

	for (auto& object : geometry)
	{
		log_info("Voxelising {} ...", object.name);
		voxelize(volume, object, preserveSurfaceMaterials, internalMaterialOveride, &progressBar, useBruteForce);
		log_info("");
	}

	auto result = estimateBounds(volume);
	Box3i bounds = result.second;
	Histogram histogram = computeHistogram(volume, bounds);

	// Print details
	log_info("Lower bound = ({},{},{})", bounds.lower().x(), bounds.lower().y(), bounds.lower().z());
	log_info("Upper bound = ({},{},{})", bounds.upper().x(), bounds.upper().y(), bounds.upper().z());
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
		log_info("No reference results for configuration");
	}

	// Save in case we need to inspect the results
	std::stringstream ss;
	ss << "../test-voxelize-" << preserveSurfaceMaterials <<
		"-" << internalMaterialOveride << "-" << useBruteForce << ".dag";
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
	log_info("");
	//testWindingNumberPerformance(triangles, bounds);
	log_info("");

	//testVoxelize(geometry);
	testVoxelize(geometry, false, 0);
	log_info("");
	testVoxelize(geometry, true, 0x0777);
	log_info("");
	testVoxelize(geometry, false, 0, true); // Brute force, slow.
	//testVoxelizeWithPreserveSurface(geometry);
	*/
	return true;
}
