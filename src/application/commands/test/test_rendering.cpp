#include "test_rendering.h"

#include "framework.h"

#include "cubiquity.h"
#include "geometry.h"
#include "extraction.h"
#include "raytracing.h"
#include "utility.h"

#include <cfloat>
#include <functional>
#include <random>

using namespace Cubiquity;
using namespace Cubiquity::Internals;
using namespace std;

class IntersectionFinder
{
public:
	IntersectionFinder(const Ray3f& ray) : mRay(ray)
	{
		mIntersection.hit = false;
		mIntersection.material = 0;
		mIntersection.distance = DBL_MAX; // Note: Might change type to float in the future?
	}

	bool operator()(NodeDAG& nodes, uint32 nodeIndex, const Box3i& bounds)
	{
		Box3d dilatedBounds = static_cast<Box3d>(bounds);
		dilatedBounds.dilate(0.5);
		RayBoxIntersection intersection = intersect(mRay, dilatedBounds);

		if (!intersection) { return false; } // Stop traversal if the ray missed the node.

		if ((isMaterialNode(nodeIndex)) && (nodeIndex > 0)) // Non-empty leaf node
		{
			// Discard nodes behind the start point, but include the node we start in.
			if (intersection.exit > 0.0)
			{
				if (intersection.entry < mIntersection.distance) // Is it closer than any other interection we foud?
				{
					mIntersection.hit = true;
					mIntersection.distance = intersection.entry;
					mIntersection.material = static_cast<MaterialId>(nodeIndex);
				}
			}
		}

		// Should early out here if we miss the node?
		return true;
	}

public:
	Ray3d mRay;
	RayVolumeIntersection mIntersection;
};

RayVolumeIntersection traceRayRef(Volume& volume, Ray3f ray)
{
	IntersectionFinder intersectionFinder(ray);
	visitVolumeNodes(volume, intersectionFinder);

	return intersectionFinder.mIntersection;
}

bool testRaytracingBehaviour()
{
	Volume volume;
	volume.load("../data/tests/axis.dag");

	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	Box3f bounds(
		{ static_cast<float>(lower_x), static_cast<float>(lower_y), static_cast<float>(lower_z) },
		{ static_cast<float>(upper_x), static_cast<float>(upper_y), static_cast<float>(upper_z) });
	Box3fSampler sampler(bounds);

	uint hitCount = 0;
	uint hitCountRef = 0;
	float maxError = 0.0f;

	Timer timer;

	SubDAGArray subDAGs = findSubDAGs(
		Internals::getNodes(volume).nodes(), getRootNodeIndex(volume));

	const uint rayCount = 1000;
	for (uint i = 0; i < rayCount; i++)
	{
		vec3f origin = sampler.next();
		vec3f target = sampler.next();
		vec3f dir = target - origin;
		dir = normalize(dir);
		Ray3f ray(origin, dir);

		RayVolumeIntersection intersection = intersectVolume(volume, subDAGs, ray, true);
		if (intersection.hit) { hitCount++; }

		RayVolumeIntersection intersectionRef = traceRayRef(volume, static_cast<Ray3f>(ray));
		if (intersectionRef.hit) { hitCountRef++; }

		if (intersection.hit && intersectionRef.hit)
		{
			float error = std::abs(intersection.distance - intersectionRef.distance);
			maxError = std::max(error, maxError);

			check(intersection.material, intersectionRef.material);
		}
	}

	log_info("{} seconds", timer.elapsedTimeInSeconds());

	log_info("Hit count = {} out of {}", hitCount, rayCount);
	check(hitCount, hitCountRef);

	log_info("Max error = {}", maxError);
	check((maxError < 0.001f), true);

	return (hitCount == hitCountRef) && ((maxError < 0.001f));
}

bool testRaytracingPerformance()
{
	Volume volume;
	volume.load("../data/tests/axis.dag");

	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	Box3f bounds(
		{ static_cast<float>(lower_x), static_cast<float>(lower_y), static_cast<float>(lower_z) },
		{ static_cast<float>(upper_x), static_cast<float>(upper_y), static_cast<float>(upper_z) });
	Box3fSampler sampler(bounds);

	uint hitCount = 0;

	Timer timer;

	SubDAGArray subDAGs = findSubDAGs(
		Internals::getNodes(volume).nodes(), getRootNodeIndex(volume));

	const uint rayCount = 1000000;
	for (uint i = 0; i < rayCount; i++)
	{
		vec3f origin = sampler.next();
		vec3f target = sampler.next();
		vec3f dir = target - origin;
		dir = normalize(dir);
		Ray3f ray(origin, dir);

		RayVolumeIntersection intersection = intersectVolume(volume, subDAGs, ray, false);
		if (intersection.hit) { hitCount++; }
	}

	log_info("Traced {} rays in {} seconds", rayCount, timer.elapsedTimeInSeconds());
	log_info("Hit count = {} out of {}", hitCount, rayCount);
	check(hitCount, 124084);

	return true;
}
