#ifndef CUBIQUITY_APP_RANDOM3D_H
#define CUBIQUITY_APP_RANDOM3D_H

#include "base/types.h"

class uniform_vec3i_distribution
{
public:
	uniform_vec3i_distribution(const ivec3& lower, const ivec3& upper)
		: eng(0)
		, dist_x(lower.x, upper.x)
		, dist_y(lower.y, upper.y)
		, dist_z(lower.z, upper.z) {}

	ivec3 operator()()
	{
		return ivec3({ dist_x(eng), dist_y(eng), dist_z(eng) });
	}

private:
	std::minstd_rand eng;
	std::uniform_int_distribution<int32_t> dist_x;
	std::uniform_int_distribution<int32_t> dist_y;
	std::uniform_int_distribution<int32_t> dist_z;
};

class uniform_vec3f_distribution
{
public:
	uniform_vec3f_distribution(const vec3& lower, const vec3& upper)
		: eng(0)
		, dist_x(lower.x, upper.x)
		, dist_y(lower.y, upper.y)
		, dist_z(lower.z, upper.z) {
	}

	vec3 operator()()
	{
		return vec3({ dist_x(eng), dist_y(eng), dist_z(eng) });
	}

private:
	std::minstd_rand eng;
	std::uniform_real_distribution<float> dist_x;
	std::uniform_real_distribution<float> dist_y;
	std::uniform_real_distribution<float> dist_z;
};

#endif // CUBIQUITY_APP_RANDOM3D_H