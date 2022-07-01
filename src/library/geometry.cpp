/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#include "geometry.h"

#include <algorithm>
#include <limits>

namespace Cubiquity
{
	template <typename T> int sign(T val) {
		return static_cast<int>(T(0) < val) - static_cast<int>(val < T(0));
	}

	float clamp(float n, float lower, float upper) {
		return std::max(lower, std::min(n, upper));
	}

	// Explanation: http://www.iquilezles.org/www/articles/triangledistance/triangledistance.htm
	float distance(const Vector3f& point, const Triangle& triangle)
	{
		auto dotWithSelf = [](const Vector3f& vector) { return dot(vector, vector); };

		// Precompute vectors
		Vector3f edges[3], verticesToPoint[3];
		for (uint i = 0; i < 3; i++)
		{
			edges[i] = triangle.vertices[(i+1)%3] - triangle.vertices[i]; 
			verticesToPoint[i] = point - triangle.vertices[i];
		}
		Vector3f normal = cross(edges[0], edges[2]);

		// Determine whether the point projects inside of the triangle
		bool inside =
			dot(cross(edges[0], normal), verticesToPoint[0]) > 0.0f &&
			dot(cross(edges[1], normal), verticesToPoint[1]) > 0.0f &&
			dot(cross(edges[2], normal), verticesToPoint[2]) > 0.0f;

		float distanceSquared = std::numeric_limits<float>::max();
		if (inside)
		{
			// If it does project inside then the required distance is
			// computed between the query point and the projected point.
			distanceSquared = dot(normal, verticesToPoint[0]) *
				dot(normal, verticesToPoint[0]) / dotWithSelf(normal);
		}
		else
		{
			// But if it projects outside then we have to project onto the
			// three edges and find the one nearest to out query point.
			for (uint i = 0; i < 3; i++)
			{
				distanceSquared =
					std::min(dotWithSelf(edges[i] * clamp(dot(edges[i], verticesToPoint[i]) /
					dotWithSelf(edges[i]), 0.0f, 1.0f) - verticesToPoint[i]), distanceSquared);
			}
		}

		return sqrt(distanceSquared);
	}

	void Triangle::flip()
	{
		std::swap(vertices[1], vertices[2]);
	}

	void Triangle::translate(const Cubiquity::Vector3f& dir)
	{
		for (auto& vertex : vertices)
		{
			vertex += dir;
		}
	}

	void Triangle::scale(float factor)
	{
		for (auto& vertex : vertices)
		{
			vertex *= factor;
		}
	}

	float Triangle::sideLength(uint32 index) const
	{
		assert(index < 3);
		return length(vertices[index] - vertices[(index + 1) % 3]);
	}

	Vector3f Triangle::computeNormal() const
	{
		// FIXME - Does this normal need to be flipped?
		return normalize(cross(vertices[2] - vertices[0], vertices[1] - vertices[0]));
	}

	float Triangle::area() const
	{
		// Heron's formula - May be inaccurate for long, thin
		// triangles: https://math.stackexchange.com/a/1483219
		float a = length(vertices[0] - vertices[1]);
		float b = length(vertices[1] - vertices[2]);
		float c = length(vertices[2] - vertices[0]);
		float s = 0.5f * (a + b + c);
		float area = sqrt(s*(s-a)*(s-b)*(s-c));
		return area;
	}

	Box3f computeBounds(const TriangleList& triangles)
	{
		Cubiquity::Box3f bounds;

		for (const Triangle& triangle : triangles)
		{
			bounds.accumulate(triangle.vertices[0]);
			bounds.accumulate(triangle.vertices[1]);
			bounds.accumulate(triangle.vertices[2]);
		}

		return bounds;
	}

	void translate(TriangleList& triangles, const Cubiquity::Vector3f& dir)
	{
		for (auto& triangle : triangles)
		{
			triangle.translate(dir);
		}
	}

	void scale(TriangleList& triangles, float factor)
	{
		for (auto& triangle : triangles)
		{
			triangle.scale(factor);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	//   _____      _                          _   _                                              //
	//  |_   _|    | |                        | | (_)                                             //
	//    | | _ __ | |_ ___ _ __ ___  ___  ___| |_ _  ___  _ __  ___                              //
	//    | || '_ \| __/ _ \ '__/ __|/ _ \/ __| __| |/ _ \| '_ \/ __|                             //
	//   _| || | | | ||  __/ |  \__ \  __/ (__| |_| | (_) | | | \__ \                             //
	//   \___/_| |_|\__\___|_|  |___/\___|\___|\__|_|\___/|_| |_|___/                             //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////

	// See https://tavianator.com/fast-branchless-raybounding-box-intersections/
	RayBoxIntersection intersect(const Ray3f& ray, const Box3f& box)
	{
		// Inverse direction could be precomputed and stored in the ray
		// if we find we often intersect the same ray with multiple boxes.
		const Vector3f invDir = Vector3f({ 1.0f, 1.0f, 1.0f }) / ray.mDir;

		const Vector3f lower = (box.lower() - ray.mOrigin) * invDir;
		const Vector3f upper = (box.upper() - ray.mOrigin) * invDir;

		const Vector3f minCorner = min(lower, upper);
		const Vector3f maxCorner = max(lower, upper);

		RayBoxIntersection intersection;
		intersection.entry = *(std::max_element(minCorner.begin(), minCorner.end()));
		intersection.exit  = *(std::min_element(maxCorner.begin(), maxCorner.end()));
		return intersection;
	}
}
