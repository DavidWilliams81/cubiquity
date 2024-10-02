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

	// Implementation of geometric solution of ray triangle intersection. Could also
	// consider the Moller-Trumbore algorithm in the future though it looks more complex.
	// See https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-
	//     rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
	bool intersect(const Ray3f& ray, const Triangle& triangle, float& t)
	{
		const auto& verts = triangle.vertices;
		// Compute the normal (we could pass this in if intersecting multiple
		// rays with the same triangle). Does not need to be normalised.
		const Vector3f normal = cross(verts[1] - verts[0], verts[2] - verts[0]);

		// No intersection if ray and triangle are parallel.
		const float epsilon = 1e-6f;
		const float normalDotRayDir = dot(normal, ray.mDir);
		if (fabs(normalDotRayDir) < epsilon) { return false; } // Miss

		// Calculate d value in equation of a plane and then t (distance along ray). This
		// calculation of t seems more complex than the version used in ray-disc intersection
		// and I'm not quite sure why: https://www.scratchapixel.com/lessons/3d-basic-rendering/
		//     minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection.html
		float d = -dot(normal, verts[0]);
		t = -(dot(normal, ray.mOrigin) + d) / normalDotRayDir;

		// No intersection if triangle is behind the ray.
		if (t < 0) { return false; } // Miss

		// Compute the intersection point.
		Vector3f p = ray.mOrigin + (ray.mDir * t);

		// Check if the point is inside the triangle;
		for (int i = 0; i < 3; i++)
		{
			Vector3f edge = verts[(i+1)%3] - verts[i];
			Vector3f vp = p - verts[i];
			Vector3f c = cross(edge, vp); // Perpendicular to triangle
			if (dot(normal, c) < 0) {
				return false; // Point is outside triangle (miss)
			}
		}

		return true; // Hit
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

	float Triangle::sideLength(int index) const
	{
		return length(vertices[(index + 1) % 3] - vertices[index]);
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

	Cubiquity::Vector3f Triangle::centre() const
	{
		return (vertices[0] + vertices[1] + vertices[2]) / 3.0f;
	}

	// Note: Could templatize on container type.
	Box3f computeBounds(const std::array<Cubiquity::Vector3f, 3>& points)
	{
		Cubiquity::Box3f bounds;
		for (const Cubiquity::Vector3f& point : points)	{
			bounds.accumulate(point);
		}
		return bounds;
	}

	Box3f computeBounds(ConstTriangleSpan triangles)
	{
		Cubiquity::Box3f bounds;
		for (const Triangle& triangle : triangles) {
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
}
