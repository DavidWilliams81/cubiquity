#include "test_rendering.h"

#include "framework.h"
#include "base/random3d.h"
#include "base/ray.h"

#include "cubiquity.h"
#include "extraction.h"
#include "raytracing.h"
#include "utility.h"

#include <cfloat>
#include <functional>
#include <random>

using Cubiquity::Volume;
using Cubiquity::uint;

using namespace std;

bool testRaytracingBehaviour()
{
	Volume volume;
	volume.load("../data/tests/axis.dag");

	uint8_t outside_material;
	int32_t lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	uniform_vec3f_distribution rand_vec3f(vec3(lower_x, lower_y, lower_z),
		vec3(upper_x, upper_y, upper_z));

	uint hitCount = 0;
	uint hitCountRef = 0;
	float maxError = 0.0f;

	Cubiquity::Timer timer;

	Cubiquity::SubDAGArray subDAGs = Cubiquity::findSubDAGs(
		Cubiquity::Internals::getNodes(volume).nodes(), Cubiquity::getRootNodeIndex(volume));

	const uint rayCount = 1000;
	for (uint i = 0; i < rayCount; i++)
	{
		vec3 origin = rand_vec3f();
		vec3 target = rand_vec3f();
		vec3 dir = target - origin;
		dir = normalize(dir);
		Ray3f ray(origin, dir);

		Cubiquity::RayVolumeIntersection intersection = intersectVolume(volume, subDAGs,
			ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
			ray.mDir.x, ray.mDir.y, ray.mDir.z, true);
		if (intersection.hit) { hitCount++; }

		Cubiquity::RayVolumeIntersection intersectionRef = traceRayRef(volume,
			ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
			ray.mDir.x, ray.mDir.y, ray.mDir.z);
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

	uint8_t outside_material;
	int32_t lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);

	uniform_vec3f_distribution random_vec3(vec3(lower_x, lower_y, lower_z),
		vec3(upper_x, upper_y, upper_z));

	uint hitCount = 0;

	Cubiquity::Timer timer;

	Cubiquity::SubDAGArray subDAGs = Cubiquity::findSubDAGs(
		Cubiquity::Internals::getNodes(volume).nodes(), Cubiquity::getRootNodeIndex(volume));

	const uint rayCount = 1000000;
	for (uint i = 0; i < rayCount; i++)
	{
		vec3 origin = random_vec3();
		vec3 target = random_vec3();
		vec3 dir = target - origin;
		dir = normalize(dir);
		Ray3f ray(origin, dir);

		Cubiquity::RayVolumeIntersection intersection = intersectVolume(volume, subDAGs,
			ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
			ray.mDir.x, ray.mDir.y, ray.mDir.z, false);
		if (intersection.hit) { hitCount++; }
	}

	log_info("Traced {} rays in {} seconds", rayCount, timer.elapsedTimeInSeconds());
	log_info("Hit count = {} out of {}", hitCount, rayCount);
	check(hitCount, 124084);

	return true;
}
