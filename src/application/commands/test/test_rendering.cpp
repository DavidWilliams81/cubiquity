#include "test_rendering.h"

#include "framework.h"

#include "geometry.h"
#include "rendering.h"
#include "utility.h"

#include <functional>
#include <iostream>
#include <random>

using namespace Cubiquity;
using namespace Cubiquity::Internals;
using namespace std;

bool testRasterization()
{
	std::cout << "Running rasterization test..." << std::endl;

	const uint32_t maskSize = 1024;

	VisibilityMask* mask = new VisibilityMask(maskSize, maskSize);
	mask->clear();

	std::minstd_rand simple_rand;

	Timer timer;

	// Note, this is a slightly distorted cube because it was originally just an arbitrary hexagon (plus a couple of points 
	// inside) in 2d, which I then tried to map to the corners of a cube as I moved from using drawConvexPolygon to drawQuads().
	Vector2f corners2DFloat[8] =
	{
		{ -1.0f, -0.4f },
		{ -0.1f, -1.0f },
		{ -0.9f,  0.4f },
		{ -0.1f, -0.1f },
		{  0.2f,  0.2f },
		{  0.9f, -0.4f },
		{  0.1f,  1.0f },
		{  1.0f,  0.4f }
	};

	int iterations = 1000;
	for (int iter = 0; iter < iterations; iter++)
	{
		mask->clear();
		simple_rand.seed(42);

		for (int ct = 0; ct < 5000; ct++)
		{
			Vector2i centre({ simple_rand() % maskSize, simple_rand() % maskSize });

			// Note that the base polygon has size of approx two (-1.0 to + 1.0)
			// FIXME - It would be nice to test with zero-size (or very timy) polygons, but
			// these currently don't match between the reference renderer and the fancy bitwise one.
			// I'm not sure what they should do, but they should probably at least be consistent. 
			const float scaleFactor = static_cast<float>(simple_rand() % 81) / 10.0f + 1.0f; // From 1.0 to 9.0

			PolygonVertexArray corners2D;
			for (int i = 0; i < 8; i++)
			{
				Vector2f vertexAsFloat = corners2DFloat[i] * scaleFactor;
				Vector2i vertexAsInt({ static_cast<int>(vertexAsFloat.x() + 0.5f), static_cast<int>(vertexAsFloat.y() + 0.5f) }); // Add half and cast

				corners2D[i] = centre + vertexAsInt;
			}

			// Normally we would determine the set of front faces from the camera position.
			// We don't have camera information for this test, so just draw all faces.
			FrontFaces frontFaces = { true, true, true, true, true, true };

			mask->drawNode(corners2D, frontFaces, true);
		}
	}

	float elapsedTime = timer.elapsedTimeInMilliSeconds();

	cout << "\tTime elapsed = " << elapsedTime << "ms." << endl;

	saveVisibilityMaskAsImage(*mask, "TestRasterizationMask.ppm");

	// Tile size affects memory layout and hence hash.
	uint32_t expectedHash = 0;
	if (VisibilityMask::TileSize == 4)
	{
		expectedHash = 1819352790;
	}
	else if (VisibilityMask::TileSize == 8)
	{
		expectedHash = 3463318368;
	}

	check(mask->hash(), expectedHash);

	bool result = mask->hash() == expectedHash;

	delete mask;

	return result;
}

class IntersectionFinder
{
public:
	IntersectionFinder(const Ray3f& ray) : mRay(ray)
	{
		mIntersection.material = 0;
		mIntersection.distance = 10000000000.0f;
	}

	bool operator()(NodeDAG& nodes, uint32 nodeIndex, Box3i bounds)
	{
		Box3f dilatedBounds = static_cast<Box3f>(bounds);
		dilatedBounds.dilate(0.5f);
		RayBoxIntersection intersection = intersect(mRay, dilatedBounds);

		if (!intersection) { return false; } // Stop traversal if the ray missed the node.

		if ((isMaterialNode(nodeIndex)) && (nodeIndex > 0)) // Non-empty leaf node
		{
			// Discard nodes behind the start point, but include the node we start in.
			if (intersection.exit > 0.0)
			{
				if (intersection.entry < mIntersection.distance) // Is it closer than any other interection we foud?
				{
					mIntersection.distance = intersection.entry;
					mIntersection.material = static_cast<MaterialId>(nodeIndex);
				}
			}
		}

		// Should early out here if we miss the node?
		return true;
	}

public:
	Ray3f mRay;
	RayVolumeIntersection mIntersection;
};

RayVolumeIntersection traceRayRef(Volume& volume, Ray3f ray)
{
	IntersectionFinder intersectionFinder(ray);
	traverseNodesRecursive(volume, intersectionFinder);

	return intersectionFinder.mIntersection;
}

bool testRaytracingBehaviour()
{
	/*Volume volume;
	volume.setVoxel(0, 0, 0, 1);

	Vector3f origin(-1, -1, -1);
	Vector3f dir = normalize(Vector3f(1, 1, 1));
	Ray3f ray(origin, dir);

	RayVolumeIntersection intersection = intersectVolume(volume, ray);
	std::cout << intersection.distance << std::endl;

	return true;*/

	Volume volume;
	volume.load("../data/axis.vol");

	Box3f bounds = static_cast<Box3f>(estimateBounds(volume).second);
	Box3fSampler sampler(bounds);

	uint hitCount = 0;
	uint hitCountRef = 0;
	float maxError = 0.0f;

	Timer timer;

	const uint rayCount = 1000;
	for (uint i = 0; i < rayCount; i++)
	{
		Vector3d origin = static_cast<Vector3d>(sampler.next());
		Vector3d target = static_cast<Vector3d>(sampler.next());
		Vector3d dir = target - origin;
		dir = normalize(dir);
		Ray3d ray(origin, dir);

		RayVolumeIntersection intersection = intersectVolume(volume, ray);
		if (intersection) { hitCount++; }

		RayVolumeIntersection intersectionRef = traceRayRef(volume, static_cast<Ray3f>(ray));
		if (intersectionRef) { hitCountRef++; }

		if (intersection && intersectionRef)
		{
			float error = std::abs(intersection.distance - intersectionRef.distance);
			maxError = std::max(error, maxError);
		}

		check(intersection.material, intersectionRef.material);
	}

	std::cout << timer.elapsedTimeInSeconds() << " seconds\n";

	std::cout << "Hit count = " << hitCount << " out of " << rayCount << std::endl;
	check(hitCount, hitCountRef);

	std::cout << "Max error = " << maxError << std::endl;
	check((maxError < 0.001f), true);

	return (hitCount == hitCountRef) && ((maxError < 0.001f));
}

bool testRaytracingPerformance()
{
	Volume volume;
	volume.load("../data/axis.vol");

	Box3f bounds = static_cast<Box3f>(estimateBounds(volume).second);
	Box3fSampler sampler(bounds);

	uint hitCount = 0;

	Timer timer;

	const uint rayCount = 1000000;
	for (uint i = 0; i < rayCount; i++)
	{
		Vector3d origin = static_cast<Vector3d>(sampler.next());
		Vector3d target = static_cast<Vector3d>(sampler.next());
		Vector3d dir = target - origin;
		dir = normalize(dir);
		Ray3d ray(origin, dir);

		RayVolumeIntersection intersection = intersectVolume(volume, ray);
		if (intersection) { hitCount++; }
	}

	std::cout << "Traced " << rayCount << " rays in " << timer.elapsedTimeInSeconds() << " seconds\n";
	std::cout << "Hit count = " << hitCount << " out of " << rayCount << std::endl;
	check(hitCount, 124000);

	return true;
}