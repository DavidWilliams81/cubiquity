/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2023 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#ifndef CUBIQUITY_RAYTRACING_H
#define CUBIQUITY_RAYTRACING_H

#include "geometry.h"
#include "storage.h"

#include "geometry.h"

#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <algorithm>
#include <cstdint>
#include <numeric>

namespace Cubiquity
{
	using namespace Internals;

	typedef Vector3f vec3;
	typedef Vector4f vec4;
	typedef Vector3i ivec3;
	typedef Vector4i ivec4;
	typedef Vector3u uvec3;
	typedef Vector4u uvec4;
	typedef Vector3b bvec3;
	typedef Vector4b bvec4;

	// This structure contains an explicit 'hit' member rather than relying on a sentinal value
	// such as a negative distance or a material of zero. I have found this simplifies the code
	// and we might want to skip material identification anyway (if traversal is stopped early
	// due to LOD). There is some redundancy storing both the distance and the position, but we
	// get the former basically for free as we traverse and the latter is very often useful.
	struct RayVolumeIntersection
	{
		bool hit;
		double distance;
		uint material;
		vec3 position;
		vec3 normal;
	};

	struct SubDAG
	{
		ivec3 lowerBound;
		int32 nodeHeight;

		uint padding0;
		uint nodeIndex;
		uint padding1, padding2;
	};

	typedef std::array<SubDAG, 8> SubDAGArray;

	SubDAGArray findSubDAGs(const Internals::NodeStore& nodes, uint32 rootNodeIndex);

	const float MAX_FOOTPRINT_DISABLED = -1.0f;
	RayVolumeIntersection intersectVolume(const Volume& volume, const SubDAGArray& subDAGs, Ray3f ray, bool computeSurfaceProperties, float maxFootprint = MAX_FOOTPRINT_DISABLED);
}

#endif // CUBIQUITY_RAYTRACING_H
