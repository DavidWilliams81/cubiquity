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

        std:: minstd_rand simple_rand;
        simple_rand.seed(42);

	Timer timer;

        //auto randCorner = std::bind(std::uniform_float_distribution<float>(-1.0,1.0), mt19937(12345));

        Vector2f corners2DFloat[8] =
        {
            Vector2f(0.2f, 0.2f),
            Vector2f(0.1f, 1.0f),
            Vector2f(1.0f, 0.4f),
            Vector2f(0.9f, -0.4f),
            Vector2f(-0.1f, -1.0f),
            Vector2f(-1.0f, -0.4f),
            Vector2f(-0.9f, 0.4f),
            Vector2f(-0.1f, -0.1f)
        };

        // Segments here are in-order
        // std::vector<int> segments = {13,6,5,5,4,4,3,3,2,2,1,1,6};

        // Segments here are out of order. The idea that in order segments are not so useful because
        // if a test passes one segment then there is a good chance it will pass it's neighbour too.
        // So better to try a completely different segment to maximise chances of early out.
        // Could maybe improve this further?
        const PolygonEdgeArray segments = PolygonEdgeArray({Edge(6,5),Edge(4,3),Edge(2,1),Edge(5,4),Edge(3,2),Edge(1,6) });

		BoundingIndices boundingIndices(5, 4, 2, 1);

        // In order of increasing Y
        const int li[6] = { 4,4,5,6,1,1 }; // Longer than required to allow reading past 'end'.
        const int ri[6] = { 4,4,3,2,1,1 }; // Longer than required to allow reading past 'end'.

	for (int ct = 0; ct < 1000; ct++)
	{
                Vector2i centre(simple_rand() % maskSize, simple_rand() % maskSize);

				const float polygonSize = simple_rand() % 14 + 2; // From 2 to 15

				PolygonVertexArray corners2D;
                for (int i = 0; i < 8; i++)
                {
                        Vector2f vertexAsFloat = corners2DFloat[i] * polygonSize;
                        Vector2i vertexAsInt(vertexAsFloat.x() + 0.5f, vertexAsFloat.y() + 0.5f); // Add half and cast

                        corners2D[i] = centre + vertexAsInt;
                }

                // Bounds are trivial to compute for our hard-coded polygon
                const Vector2i min_bounds(corners2D[5].x(), corners2D[4].y());
                const Vector2i max_bounds(corners2D[2].x(), corners2D[1].y());

				Bounds bounds;
				bounds.lower = min_bounds;
				bounds.upper = max_bounds;

                // Rasterize
                mask->drawConvexPolygon(corners2D, segments, bounds);
	}

	float elapsedTime = timer.elapsedTimeInMilliSeconds();

	cout << "\tTime elapsed = " << elapsedTime << "ms." << endl;

	saveVisibilityMaskAsImage(*mask, "TestRasterizationMask.ppm");

	uint32_t expectedHash = 3347657948;

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

	RayVolumeIntersection intersection = ray_parameter(volume, ray);
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
	for(uint i = 0; i < rayCount; i++)
	{
		Vector3d origin = static_cast<Vector3d>(sampler.next());
		Vector3d target = static_cast<Vector3d>(sampler.next());
		Vector3d dir = target - origin;
		dir = normalize(dir);
		Ray3d ray(origin, dir);

		RayVolumeIntersection intersection = ray_parameter(volume, ray);
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

		RayVolumeIntersection intersection = ray_parameter(volume, ray);
		if (intersection) { hitCount++; }
	}

	std::cout << "Traced " << rayCount << " rays in " << timer.elapsedTimeInSeconds() << " seconds\n";
	std::cout << "Hit count = " << hitCount << " out of " << rayCount << std::endl;
	check(hitCount, 123989);

	return true;
}
