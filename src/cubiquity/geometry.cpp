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

#include <limits>

namespace Cubiquity
{
	Matrix4x4f frustum_matrix(float x0, float x1, float y0, float y1, float n, float f)
	{
		const float s = -1.0f, o = n;

		Matrix4x4f result;
		result[0] = Vector4f(2*n/(x1-x0),0,0,0);
		result[1] = Vector4f(0,2*n/(y1-y0),0,0);
		result[2] = Vector4f(-s*(x0+x1)/(x1-x0),-s*(y0+y1)/(y1-y0),s*(f+o)/(f-n),s);
		result[3] = Vector4f(0,0,-(n+o)*f/(f-n),0);
		return result;
	}

	Matrix4x4f perspective_matrix(float fovy, float aspect, float n, float f)
	{
		float y = n*std::tan(fovy / 2), x = y*aspect;

		return frustum_matrix(-x, x, -y, y, n, f);
	}

	Matrix4x4f lookAtRH(const Vector3f& eye, const Vector3f& center, const Vector3f& up)
	{
		const Vector3f f(normalize(center - eye));
		const Vector3f s(normalize(cross(f, up)));
		const Vector3f u(cross(s, f));

		Matrix4x4f result;
		result[0][0] = s.x();
		result[1][0] = s.y();
		result[2][0] = s.z();
		result[0][1] = u.x();
		result[1][1] = u.y();
		result[2][1] = u.z();
		result[0][2] = -f.x();
		result[1][2] = -f.y();
		result[2][2] = -f.z();
		result[3][0] = -dot(s, eye);
		result[3][1] = -dot(u, eye);
		result[3][2] = dot(f, eye);
		return result;
	}

	Matrix4x4f translation_matrix(const Vector3f& pos)
	{
		Matrix4x4f result;
		result[3] = Vector4f(pos.x(), pos.y(), pos.z(), 1.0f);
		return result;
	}

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

	Box3f computeBounds(const Object& object)
	{
		Cubiquity::Box3f bounds;
		for (auto& subObject : object.subObjects)
		{
			Box3f triangleBounds = computeBounds(subObject.second);
			bounds.accumulate(triangleBounds);
		}
		return bounds;
	}

	Box3f computeBounds(const Geometry& geometry)
	{
		Cubiquity::Box3f bounds;
		for (auto& object : geometry)
		{
			for (auto& subObject : object.subObjects)
			{
				Box3f triangleBounds = computeBounds(subObject.second);
				bounds.accumulate(triangleBounds);
			}
		}
		return bounds;
	}

	void translate(Geometry& geometry, const Cubiquity::Vector3f& dir)
	{
		for (auto& object : geometry)
		{
			for (auto& subObject : object.subObjects)
			{
				TriangleList& triList = subObject.second;
				translate(triList, dir);
			}
		}
	}

	void scale(Geometry& geometry, float factor)
	{
		for (auto& object : geometry)
		{
			for (auto& subObject : object.subObjects)
			{
				scale(subObject.second, factor);
			}
		}
	}

	SubObject mergeSubObjects(const SubObjectList& subObjects, uint16 resultingMaterial)
	{
		TriangleList triangles;
		for (const SubObject& subObject : subObjects)
		{
			triangles.insert(triangles.end(), subObject.second.begin(), subObject.second.end());
		}

		SubObject subObject;
		subObject.first = resultingMaterial;
		subObject.second = triangles;
		return subObject;
	}

	Object mergeObjects(const Geometry& geometry, std::string resultingObjectName, uint16 resultingMaterial)
	{
		TriangleList triangles;
		for (const Object& object : geometry)
		{
			for (const SubObject& subObject : object.subObjects)
			{
				triangles.insert(triangles.end(), subObject.second.begin(), subObject.second.end());
			}
		}

		SubObject subObject;
		subObject.first = resultingMaterial;
		subObject.second = triangles;

		Object object;
		object.name = resultingObjectName;
		object.subObjects.push_back(subObject);
		return object;
	}

	TriangleList mergedTriangles(const Object& object)
	{
		TriangleList result;
		for (const auto& subObject : object.subObjects)
		{
			result.insert(result.end(), subObject.second.begin(), subObject.second.end());
		}
		return result;
	}

	TriangleList mergedTriangles(const Geometry& geometry)
	{
		TriangleList result;
		for (const auto& object : geometry)
		{
			TriangleList mergedObjectTris = mergedTriangles(object);
			result.insert(result.end(), mergedObjectTris.begin(), mergedObjectTris.end());
		}
		return result;
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
		const Vector3f invDir = Vector3f(1.0f) / ray.mDir;

		const Vector3f lower = (box.lower() - ray.mOrigin) * invDir;
		const Vector3f upper = (box.upper() - ray.mOrigin) * invDir;

		RayBoxIntersection intersection;
		intersection.entry = min(lower, upper).maxComponentValue();
		intersection.exit  = max(lower, upper).minComponentValue();
		return intersection;
	}
}
