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
#include "voxelization.h"

#include "geometry.h"

#include "storage.h"

#include <algorithm>
#include <cassert>
#include <cmath>

// g++ currently doesn't support std::execution.
#if defined (_MSC_VER) && _MSC_VER >= 1900
#include <execution>
#endif


#include <mutex>
#include <sstream>
#include <thread>

#ifdef CUBIQUITY_USE_AVX
	#include <immintrin.h>
#endif

namespace Cubiquity
{
	using namespace Internals;

	////////////////////////////////////////////////////////////////////////////////////////////////
	//   _    _ _           _ _               _   _                 _                             //
	//  | |  | (_)         | (_)             | \ | |               | |                            //
	//  | |  | |_ _ __   __| |_ _ __   __ _  |  \| |_   _ _ __ ___ | |__   ___ _ __ ___           //
	//  | |/\| | | '_ \ / _` | | '_ \ / _` | | . ` | | | | '_ ` _ \| '_ \ / _ \ '__/ __|          //
	//  \  /\  / | | | | (_| | | | | | (_| | | |\  | |_| | | | | | | |_) |  __/ |  \__ \          //
	//   \/  \/|_|_| |_|\__,_|_|_| |_|\__, | \_| \_/\__,_|_| |_| |_|_.__/ \___|_|  |___/          //
	//                                 __/ |                                                      //
	//                                |___/                                                       //
	////////////////////////////////////////////////////////////////////////////////////////////////
	typedef std::pair<Vector3f, Vector3f> TriangleEdge;

	// Theorectcally the threshold should be 0.5f. However, triangles often pass *exactly*
	// through the centre of a voxel due to e.g. being modelled on a grid in a DCC package,
	// and/or being scaled by interger amounts. Rounding errors mean that such triangles can
	// end up with a noisy surface, and the epsilon below helps to avoid this.
	const float gWindingNumberThreshold = 0.5 + 0.001f;

	TriangleEdge flip(const TriangleEdge& edge)
	{
		return std::make_pair(edge.second, edge.first);
	}

#ifdef CUBIQUITY_USE_AVX

	inline void cross(__m256 ax, __m256 ay, __m256 az, __m256 bx, __m256 by, __m256 bz, __m256& resx, __m256& resy, __m256& resz)
	{
		resx = _mm256_sub_ps(_mm256_mul_ps(ay, bz), _mm256_mul_ps(az, by));
		resy = _mm256_sub_ps(_mm256_mul_ps(az, bx), _mm256_mul_ps(ax, bz));
		resz = _mm256_sub_ps(_mm256_mul_ps(ax, by), _mm256_mul_ps(ay, bx));
	}

	inline __m256 dot(__m256 ax, __m256 ay, __m256 az, __m256 bx, __m256 by, __m256 bz)
	{
		__m256 axbx = _mm256_mul_ps(ax, bx);
		__m256 ayby = _mm256_mul_ps(ay, by);
		__m256 azbz = _mm256_mul_ps(az, bz);

		__m256 result = _mm256_add_ps(axbx, ayby);
		result = _mm256_add_ps(result, azbz);
		return result;
	}

	inline __m256 length(__m256 x, __m256 y, __m256 z)
	{
		return _mm256_sqrt_ps(dot(x, y, z, x, y, z));
	}

	float computeWindingNumber(const Vector3f& queryPoint, const TriangleList& triangles)
	{
		if (triangles.empty())
		{
			// FIXME - what is allowing this to happen?
			return 0.0f;
		}

		const uint32 avxWidth = 8; // Eight floats to each 256-bit AVX operation.
		const uint32 elementsPerTriangle = 9;

		// The triangle data comes in with the floats in an interleaved (packed) format. The order is
		// [ax, ay, az, bx, by, bz, cx ,cy, cz, ax, ay, az, ...] which is not optimal for SIMD processing.
		// Because AVX does not support 'gather' operations (requires AVX2) we restructure the data into
		// a planar format where the order is [ax, ax, ..., ay, ay, ..., ..., cz, cz, ...]. We also round
		// the number of triangles up to the nearest multiple of the AVX width (in terms of floats) to
		// simplify later logic. All this repacking will probably disappear in the future if/when I get a
		// CPU supporting the newer AVX2 instruction set.
		const uint32_t triCount = triangles.size();
		const uint32 triCountRoundedUp = (triCount + avxWidth - 1) & -avxWidth;

		// Reordered version of the data. FIXME - Use aligned_alloc when my compiler supports it.
		float* const planar = new float[triCountRoundedUp * elementsPerTriangle];
		/*if (!isAligned(planar, 32))
		{
		throw std::runtime_error("Incorrect memory alignment");
		}*/

		float* dst = planar;
		float* const dstEnd = dst + triCountRoundedUp * elementsPerTriangle;

		// For each element type (ax, ay, ... cy, cy)
		for (uint32 elementIndex = 0; elementIndex < elementsPerTriangle; elementIndex++)
		{
			const float* src = reinterpret_cast<const float*>(&(triangles[0])) + elementIndex;
			const float* const srcEnd = (src + triCount * elementsPerTriangle);

			// Only the end of *this particular set of elements* (not the whole planer buffer).
			float* const dstEnd = dst + triCountRoundedUp;

			// Gather the source elements into the destination buffer
			while (src < srcEnd)
			{
				*dst = *src;
				dst++;
				src += elementsPerTriangle;
			}

			// Fill any remaining elements with zeros;
			while (dst < dstEnd)
			{
				*dst = 0.0f;
				dst++;
			}
		}

		float windingNumber = 0.0f;

		__m256 qx = _mm256_set1_ps(queryPoint.x());
		__m256 qy = _mm256_set1_ps(queryPoint.y());
		__m256 qz = _mm256_set1_ps(queryPoint.z());

		for (uint32_t index = 0; index < triCountRoundedUp; index += 8)
		{
			// The code below is an (untested) implementation using  AVX2
			// gather instructions. It is saved here for future reference.
			//__m256i vindex = _mm256_set_epi32(0, 9, 18, 27, 36, 45, 54, 63); // Move out of loop
			//float* firstTriangle = (float*)(&(paddedTriangles[index]));

			//__m256 qax = _mm256_i32gather_ps(firstTriangle + 0, vindex, 4);
			//__m256 qay = _mm256_i32gather_ps(firstTriangle + 1, vindex, 4);
			//__m256 qaz = _mm256_i32gather_ps(firstTriangle + 2, vindex, 4);

			//__m256 qbx = _mm256_i32gather_ps(firstTriangle + 3, vindex, 4);
			//__m256 qby = _mm256_i32gather_ps(firstTriangle + 4, vindex, 4);
			//__m256 qbz = _mm256_i32gather_ps(firstTriangle + 5, vindex, 4);

			//__m256 qcx = _mm256_i32gather_ps(firstTriangle + 6, vindex, 4);
			//__m256 qcy = _mm256_i32gather_ps(firstTriangle + 7, vindex, 4);
			//__m256 qcz = _mm256_i32gather_ps(firstTriangle + 8, vindex, 4);

			__m256 qax = _mm256_loadu_ps(planar + (triCountRoundedUp * 0) + index);
			__m256 qay = _mm256_loadu_ps(planar + (triCountRoundedUp * 1) + index);
			__m256 qaz = _mm256_loadu_ps(planar + (triCountRoundedUp * 2) + index);

			__m256 qbx = _mm256_loadu_ps(planar + (triCountRoundedUp * 3) + index);
			__m256 qby = _mm256_loadu_ps(planar + (triCountRoundedUp * 4) + index);
			__m256 qbz = _mm256_loadu_ps(planar + (triCountRoundedUp * 5) + index);

			__m256 qcx = _mm256_loadu_ps(planar + (triCountRoundedUp * 6) + index);
			__m256 qcy = _mm256_loadu_ps(planar + (triCountRoundedUp * 7) + index);
			__m256 qcz = _mm256_loadu_ps(planar + (triCountRoundedUp * 8) + index);

			qax = _mm256_sub_ps(qax, qx);
			qay = _mm256_sub_ps(qay, qy);
			qaz = _mm256_sub_ps(qaz, qz);

			qbx = _mm256_sub_ps(qbx, qx);
			qby = _mm256_sub_ps(qby, qy);
			qbz = _mm256_sub_ps(qbz, qz);

			qcx = _mm256_sub_ps(qcx, qx);
			qcy = _mm256_sub_ps(qcy, qy);
			qcz = _mm256_sub_ps(qcz, qz);

			// In AVX we could actually compute the(approx) reciprocal length via the rsqrt()
			// intrinsic and do a multiply below instead of a divide. This might be faster, but
			// the results then deviate from our reference implementation.
			__m256 alength = length(qax, qay, qaz);
			__m256 blength = length(qbx, qby, qbz);
			__m256 clength = length(qcx, qcy, qcz);

			// Note that the reference C++ version of the code checks that the lengths are not zero
			// before normalizing. I've skipped that here as it's a bit more tricky when using SIMD.
			// Should watch out for Nans/Infs though?
			qax = _mm256_div_ps(qax, alength);
			qay = _mm256_div_ps(qay, alength);
			qaz = _mm256_div_ps(qaz, alength);

			qbx = _mm256_div_ps(qbx, blength);
			qby = _mm256_div_ps(qby, blength);
			qbz = _mm256_div_ps(qbz, blength);

			qcx = _mm256_div_ps(qcx, clength);
			qcy = _mm256_div_ps(qcy, clength);
			qcz = _mm256_div_ps(qcz, clength);

			// Subtracting qa from qb and qc is not strictly required,
			// but gives a more stable result if the triangle is far away.
			__m256 qb_qax = _mm256_sub_ps(qbx, qax);
			__m256 qb_qay = _mm256_sub_ps(qby, qay);
			__m256 qb_qaz = _mm256_sub_ps(qbz, qaz);

			__m256 qc_qax = _mm256_sub_ps(qcx, qax);
			__m256 qc_qay = _mm256_sub_ps(qcy, qay);
			__m256 qc_qaz = _mm256_sub_ps(qcz, qaz);

			// Numerator
			__m256 crossx, crossy, crossz;
			cross(qb_qax, qb_qay, qb_qaz, qc_qax, qc_qay, qc_qaz, crossx, crossy, crossz);
			__m256 numerator = dot(qax, qay, qaz, crossx, crossy, crossz);


			// Denominator
			__m256 dqaqb = dot(qax, qay, qaz, qbx, qby, qbz);
			__m256 dqaqc = dot(qax, qay, qaz, qcx, qcy, qcz);
			__m256 dqbqc = dot(qbx, qby, qbz, qcx, qcy, qcz);

			__m256 denominator = _mm256_set1_ps(1.0f);
			denominator = _mm256_add_ps(denominator, dqaqb);
			denominator = _mm256_add_ps(denominator, dqaqc);
			denominator = _mm256_add_ps(denominator, dqbqc);

			// Unfortunately AVX does not provide us with a way to compute the atan2 of an __m256
			// register. More accurately, the 'a_mm256_atan2_ps' intrinsic is part of the Intel
			// SVML library which is not freely available. I have found other SIMD atan2()
			// implementation but I think they implement and single atan2() via SIMD, rather than
			// operating on a whole register at once. Agner Fog's vector library might be an
			// exception here, but the code is complex and not public domain.
			//
			// Therefor I implement the atan2() part of the process by pulling the values out of the
			// __m256 registers and just using the normal atan2() function.
			alignas(32) float numerators[8];
			alignas(32) float denominators[8];

			_mm256_store_ps(numerators, numerator);
			_mm256_store_ps(denominators, denominator);

			for (int i = 0; i < 8; i++)
			{
				// Numerator of zero means we are on the surface (treat as no solid angle).
				if (numerators[i] != 0)
				{
					assert(std::isnormal(numerators[i]));
					assert(std::isnormal(denominators[i]));
					windingNumber += 2.0f * ::atan2(numerators[i], denominators[i]);
				}
			}
		}

		delete[] planar;

		// Normalise to [-1.0f, 1.0f] range.
		windingNumber /= 4.0 * 3.14159265358979323846;

		// Make sure we got a 'valid' result.
		assert(std::isnormal(windingNumber) || (windingNumber == 0));

		return windingNumber;
	}

#else // CUBIQUITY_USE_AVX

	float computeWindingNumber(const Vector3f& queryPoint, const TriangleList& triangles)
	{
		float windingNumber = 0.0f;

		for (const Triangle& triangle : triangles)
		{
			const auto& a = triangle.vertices[0];
			const auto& b = triangle.vertices[1];
			const auto& c = triangle.vertices[2];

			Vector3f qa = a - queryPoint;
			Vector3f qb = b - queryPoint;
			Vector3f qc = c - queryPoint;

			const float alength = length(qa);
			const float blength = length(qb);
			const float clength = length(qc);

			// Avoid potential NaNs/Infs in the result. If any of these values are
			// zero then the triangle does not contribute to the winding number.
			if (alength != 0 && blength != 0 && clength != 0)
			{
				// Normalize the vectors
				qa /= alength;
				qb /= blength;
				qc /= clength;

				// Subtracting qa from qb and qc is not strictly required,
				// but gives a more stable result if the triangle is far away.
				const float numerator = dot(qa, cross(qb - qa, qc - qa));
				const float denominator = 1.0f + dot(qa, qb) + dot(qa, qc) + dot(qb, qc);

				// Numerator of zero means we are on the surface (treat as no solid angle).
				if (numerator != 0)
				{
					assert(std::isnormal(numerator));
					assert(std::isnormal(denominator));
					windingNumber += 2.0f * ::atan2(numerator, denominator);
				}
			}
		}

		// Normalise to [-1.0f, 1.0f] range.
		windingNumber /= 4.0 * 3.14159265358979323846;

		// Make sure we got a 'valid' result.
		assert(std::isnormal(windingNumber) || (windingNumber == 0));

		return windingNumber;
	}

#endif // CUBIQUITY_USE_AVX

	float computeWindingNumber(const Vector3f& queryPoint, const ClosedTriangleTree& meshNode)
	{
		// If one child exists they must both exist
		assert((meshNode.mChildren[0] == nullptr) == (meshNode.mChildren[1] == nullptr));

		float windingNumber = 0.0f;
		bool hasChildren = meshNode.mChildren[0] && meshNode.mChildren[1];

		if (!hasChildren)
		{
			assert(meshNode.mTriangles.size() > 0);
			windingNumber += computeWindingNumber(queryPoint, meshNode.mTriangles);
			return windingNumber;

		}
		else if (!meshNode.mBounds.contains(queryPoint))
		{
			//assert(meshNode.mClosingTriangles.size() > 0);
			windingNumber += computeWindingNumber(queryPoint, meshNode.mClosingTriangles);
			return windingNumber;
		}
		else
		{
			//assert(meshNode.mTriangles.size() == 0 && meshNode.mClosingTriangles.size() == 0);
			windingNumber += computeWindingNumber(queryPoint, *(meshNode.mChildren[0]));
			windingNumber += computeWindingNumber(queryPoint, *(meshNode.mChildren[1]));

			return windingNumber;
		}

		return windingNumber;
	}

	std::vector<TriangleEdge> computeExteriorEdges(const TriangleList& triangles)
	{
		std::unordered_map<TriangleEdge, int32_t, Internals::MurmurHash3<TriangleEdge> > edgeCounts;
		for (uint ct = 0; ct < triangles.size(); ct++)
		{
			Vector3f i0 = triangles[ct].vertices[0];
			Vector3f i1 = triangles[ct].vertices[1];
			Vector3f i2 = triangles[ct].vertices[2];

			// FIXME - This assumes i0, i1 and i2 are all different
			// values (which they should be for a decent mesh).
			if (i0 < i1)
			{
				edgeCounts[TriangleEdge(i0, i1)]++;
			}
			else
			{
				edgeCounts[TriangleEdge(i1, i0)]--;
			}

			if (i1 < i2)
			{
				edgeCounts[TriangleEdge(i1, i2)]++;
			}
			else
			{
				edgeCounts[TriangleEdge(i2, i1)]--;
			}

			if (i2 < i0)
			{
				edgeCounts[TriangleEdge(i2, i0)]++;
			}
			else
			{
				edgeCounts[TriangleEdge(i0, i2)]--;
			}
		}

		std::vector<TriangleEdge> exteriorEdges;

		for (const auto& edgeCount : edgeCounts)
		{
			if (edgeCount.second > 0)
			{
				exteriorEdges.push_back(edgeCount.first);
			}

			if (edgeCount.second < 0)
			{
				exteriorEdges.push_back(flip(edgeCount.first));
			}
		}

		return exteriorEdges;
	}

	TriangleList computeClosingTriangles(const TriangleList& triangles, std::vector<TriangleEdge> exteriorEdges)
	{
		TriangleList closingTriangles;
		for (const auto& edge : exteriorEdges)
		{
			// Arbitrary vertex, but must be inside the bounding box.
			// Using one of the original vertices is easiest.
			Vector3f arbitraryVertex = triangles[0].vertices[0];

			Triangle closingTriangle(edge.first, edge.second, arbitraryVertex);
			closingTriangles.push_back(closingTriangle);
		}

		return closingTriangles;
	}

	void splitTriangles(const TriangleList& triangles, TriangleList& triangles0, TriangleList& triangles1)
	{
		Cubiquity::Box3f input;

		Vector3f mean(0.0f);

		for (uint ct = 0; ct < triangles.size(); ct++)
		{
			Vector3f i0 = triangles[ct].vertices[0];
			Vector3f i1 = triangles[ct].vertices[1];
			Vector3f i2 = triangles[ct].vertices[2];

			Vector3f triangleCentre = (i0 + i1 + i2) / 3.0f;

			mean += triangleCentre;

			input.accumulate(triangleCentre);
		}

		mean /= triangles.size();

		float width = input.upper().x() - input.lower().x();
		float height = input.upper().y() - input.lower().y();
		float depth = input.upper().z() - input.lower().z();

		// I'd rather use the mean for the threshold, but it gives
		// different results and I'm not yet sure it is better.
		Vector3f threshold = input.centre();

		for (const auto triangle : triangles)
		{
			Vector3f triangleCentre = (triangle.vertices[0] + 
				triangle.vertices[1] + triangle.vertices[2]) / 3.0f;

			if ((width > height) && (width > depth))
			{
				if (triangleCentre.x() <= threshold.x())
				{
					triangles0.push_back(triangle);
				}
				else
				{
					triangles1.push_back(triangle);
				}
			}
			else if (height > depth)
			{
				if (triangleCentre.y() <= threshold.y())
				{
					triangles0.push_back(triangle);
				}
				else
				{
					triangles1.push_back(triangle);
				}
			}
			else
			{
				if (triangleCentre.z() <= threshold.z())
				{
					triangles0.push_back(triangle);
				}
				else
				{
					triangles1.push_back(triangle);
				}
			}
		}
	}

	//! \param triangles A list of triangles to be converted into a ClosedTriangleTree.
	//! \param tidyInput Specifies whether the triangles should be preprocessed by snapping them to
	//!                  a sub-voxel grid and removing degnerates. This process improves the
	//!                  robustness of the algoithm but may introduce small deviations compared to
	//!                  the reference result. It is recommended to set it to the default value of
	//!                  'true' as it is mostly provided for internal use.
	ClosedTriangleTree::ClosedTriangleTree(const TriangleList &triangles, bool tidyInput)
		: mTriangles(triangles)
	{
		assert(mTriangles.size() > 0);

		if(tidyInput)
		{
			// Moves each vertex of each triangle so that it lies exactly on a grid. This helps to
			// reduce errors if two supposely-identical positions are actually slightly different.
			// This is important because if two edges are supposed to exactly line up, but they
			// don't due to floating point errors, then those edges may be incorrectly identified
			// as exterior edges and hence the generated closing triangles will be sub-optimal.
			const int gridUnitsPerVoxel = 16;
			for (Triangle& triangle : mTriangles)
			{
				for (Vector3f& vertex : triangle.vertices)
				{
					vertex *= gridUnitsPerVoxel;

					vertex[0] = std::round(vertex[0]);
					vertex[1] = std::round(vertex[1]);
					vertex[2] = std::round(vertex[2]);

					vertex /= gridUnitsPerVoxel;
				}
			}

			// Remove degenerate (zero-area) triangles. We've already snapped vertex positions
			// to a fine grid, so these exact floating-point comparisons should be well behaved.
			auto newEnd = std::remove_if(mTriangles.begin(), mTriangles.end(),
			[](auto const x)
			{
				return (x.vertices[0] == x.vertices[1]) ||
					   (x.vertices[1] == x.vertices[2]) ||
					   (x.vertices[2] == x.vertices[0]);
			});

			ptrdiff_t removedCount = std::distance(newEnd, mTriangles.end());
			mTriangles.erase(newEnd, mTriangles.end());

			// This is just a debug message really, as there may have been tiny (but valid)
			// triangles which only became degenerate after snapping-to-grid. We need a better
			// logging mechanism really.
			if (removedCount > 0)
			{
				std::cout << "WARNING - Removed " << removedCount << " degenerate triangles" << std::endl;
			}
		}

		// The paper doesn't mention it, but presumably we now need to expand the bounds to match
		// the triangles it contains? Triangles were added if their centre was inside the Box, but
		// their vertices may not be.
		for (const auto& triangle : mTriangles)
		{
			mBounds.accumulate(triangle.vertices[0]);
			mBounds.accumulate(triangle.vertices[1]);
			mBounds.accumulate(triangle.vertices[2]);
		}

		// The paper suggests a value of 100 as a threshold for when to stop splitting, but
		// basic testing has indiciated that a smaller value gives better runtime performance.
		const int triangleThreshold = 20;

		// From this point on we are implementing Algorithm 1 from the paper.
		auto exteriorEdges = computeExteriorEdges(mTriangles);
		mClosingTriangles = computeClosingTriangles(mTriangles, exteriorEdges);
		assert(exteriorEdges.size() == mClosingTriangles.size());

		if (mTriangles.size() < triangleThreshold || mClosingTriangles.size() >= mTriangles.size())
		{
			return;
		}

		std::vector<Triangle> triangles0, triangles1;
		splitTriangles(mTriangles, triangles0, triangles1);
		assert(triangles0.size() + triangles1.size() == mTriangles.size());

		// The original triangles aren't needed now that we have the split lists.
		mTriangles.clear();

		// Note that we do not need to call the 'tidyInput' process when recursively
		// building the tree, because it should already have been done for the root node.
		mChildren[0] = std::unique_ptr<ClosedTriangleTree>(new ClosedTriangleTree(triangles0, false));
		mChildren[1] = std::unique_ptr<ClosedTriangleTree>(new ClosedTriangleTree(triangles1, false));
	}

	Surface::Surface()
		:mainMaterial(0) {}

	void Surface::addTriangle(const Triangle& tri, MaterialId matId)
	{
		triangles.push_back(tri);
		materials.push_back(matId);
	}

	void Surface::build()
	{
		assert(triangles.size() == materials.size());

		// The mean winding number of the samples can be used to determine how well-formed the surface is.
		//
		// For a well-formed mesh we expect samples to have a winding number of +/-1.0 inside the surface
		// and 0.0 outside the surface. If we assume that there are approximately the same number of inside
		// vs. outside samples then we would expect to see a mean winding number of +/-0.5 for a well formed 
		// mesh, but it would be closer zero for a inconsistently wound mesh (as there would be a mix of
		// positive and negative winding numbers).
		//
		// In practice there may not be an equal number of inside vs. outside samples. If the mesh is long
		// and thin then it may actually miss most of the samples, which will push the mean winding number
		// towards zero and incorrectly suggest that it is inconsistently-wound. On the other hand, if it
		// is a large cubic mesh it may catch most of the samples and push the winding number towards one.
		// Note that this is less of a problem, because it should still detect an inconsistenly-wound
		// surfaces as they will still have a mix of positive and negative values.
		//
		// To resolve the above we only include samples in the mean if their absolute value is above a
		// threshold. A well formed mesh will then have amean winding number of 1.0 (as all the included
		// samples will have a value of 1.0). If the mean winding number has a smaller magnitude it is likely
		// to suggest that the surface is either inconsistenly-wound or is not properly closed (both of these
		// make it difficult to fill the interiour). A mean winding number greater than 1.0 has been seen to
		// occur when an otherwise-closed surface has extra details or protrusions seperately modelled on it.
		// Anecdotally these do not see to cause such a problem for voxelisation (possibly just because they
		// are typically small compred to the surface itself?)
		const int sampleCount = 100;
		int validSamples = 0;
		meanWindingNumber = 0.0f;

		bounds = computeBounds(triangles);
		Cubiquity::Box3fSampler sampler(bounds);
		for (uint32_t i = 0; i < sampleCount; i++)
		{
			Vector3f point = sampler.next();
			float windingNumber = computeWindingNumber(point, triangles);

			const float inclusionThreshold = 0.5f; // No real basis for this value. It can be tweaked as required.
			if (std::abs(windingNumber) > inclusionThreshold)
			{
				meanWindingNumber += windingNumber;
				validSamples++;
			}
		}
		meanWindingNumber /= validSamples;

		// Warn the user if the surface does not appear to be well formed.
		const float windingNumberTolerance = 0.001;
		if(std::abs(meanWindingNumber) < (1.0f - windingNumberTolerance))
		{
			std::cerr << "Warning: Surface " << name << " has a mean winding number of "
				<< meanWindingNumber << ", which is closer to zero than expected. Perhaps it "
				<< "contains inconsistently-wound faces or is not properly closed?" << std::endl;
		}

		if (std::abs(meanWindingNumber) > (1.0f + windingNumberTolerance))
		{
			std::cerr << "Warning: Surface " << name << " has a mean winding number of "
				<< meanWindingNumber << ", which is further from zero than expected. Perhaps "
				<< "it contains doubled-up geometry or separate surface details?" << std::endl;
		}

		// Find the main material of the surface as the material which covers the greatest area. 
		// Attempts to use winding numbers have less obvious behaviour for open/inverted surfaces.
		std::vector<float> areas(MaterialCount);
		std::fill(areas.begin(), areas.end(), 0.0f);
		for (uint i = 0; i < triangles.size(); i++)
		{
			areas[materials[i]] += triangles[i].area();
		}
		mainMaterial = std::distance(areas.begin(), std::max_element(areas.begin(), areas.end()));

		// Build the closed triangle tree representaton from the list of triangles.
		closedTriangleTree = std::make_unique<ClosedTriangleTree>(triangles);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	//   ___________   _____                   _____                               _              //
	//  |____ |  _  \ /  ___|                 /  __ \                             (_)             //
	//      / / | | | \ `--.  ___ __ _ _ __   | /  \/ ___  _ ____   _____ _ __ ___ _  ___  _ __   //
	//      \ \ | | |  `--. \/ __/ _` | '_ \  | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \  //
	//  .___/ / |/ /  /\__/ / (_| (_| | | | | | \__/\ (_) | | | \ V /  __/ |  \__ \ | (_) | | | | //
	//  \____/|___/   \____/ \___\__,_|_| |_|  \____/\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_| //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                            //
	// This function scan-converts a 3D triangle into a corresponding set of '3D fragments'       //
	// (voxel positions). It does this by iterating over all voxels in the triangle's bounding    //
	// box and outputting those within a certain distance threshold of the triangle. Unless       //
	// specified the threshold is chosen as described in the paper "An Accurate Method To         //
	// Voxelize Polygonal Meshes" by Huang et al.                                                 //
	//                                                                                            //
	// The result should in theory be minimal, but I haven't given much thought to the corners or //
	// edges and so there may be a few extra fragments output where different triangles meet. The //
	// above paper has further ideas on this in Section 5. For our application we don't care much //
	// about extra fragments (and would rather have too many than too few) because we run them    //
	// throuh the inside-outside test anyway, so it's just a little extra work but still a        //
	// correct result.                                                                            //
	//                                                                                            //
	// The 'thickness' of a surface is defined by its separability which (in perhaps a very loose //
	// sense) is conceptually the opposite of its connectivity. It is apparantly hard to define   //
	// the connectivity of a surface and so separability is used instead. A 6-separating surface  //
	// is one which cannot be penetrated by a 6-connected line, and it is thicker than a 26-      //
	// separating surface (which cannot be penetrated by a 26-connected line but *can* be         //
	// penetrated by a 6-connected line). See the sections 2 and 3 of the paper for further       //
	// details.                                                                                   //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////
	float computeMinThickness(const Triangle& triangle)
	{
		// Note: There seems to be some errors/inconsistencies in the paper. Perhaps it is a pre-print? I'm looking at this copy:
		//
		//		https://pdfs.semanticscholar.org/2827/c00e687c870f84e5006b851aa038ca9f971a.pdf
		//
		// Looking at the 2D case (section 4) they state:
		//
		// "For generating a 4-separating line we use t = ..., denoted by t8. For generating 8-separating lines we define t4 = ..."
		//
		// So t*8* is used for *4*-seperating lines, and vice-versa?
		//
		// But in the 3D case (section 5) they state:
		//
		// "For generating a 26-separating surface we use t = ..., denoted by t26. For generating 6-separating line we define t6 = ..."
		//
		// So now t*26* is used for *26*-seperating surfaces, t*6* is used for *6*-seperating surfaces (note, they actually said 'line' here...)?
		//
		// The 2D and 3d cases seem to be opposite? It looks like there is some bad copy/paste going on.
		//
		// Unless I find a better copy I think I'll just do what I think is correct, ensuring that a seperabilty of 6 gives the thicker surface.

		const uint seperation = 26; // Amount of seperation, must be '6' (thicker) or '26' (thinner).
		const float L = 1.0f; // From figure 10 in the paper, our voxels have unity side length.

		const Vector3f a = triangle.vertices[0];
		const Vector3f b = triangle.vertices[1];
		const Vector3f c = triangle.vertices[2];

		// Compute threshold based on Section 4.2 of "An Accurate Method To Voxelize Polygonal Meshes" by Huang et al.
		Vector3f planeNormal = cross(b - a, a - c);

		// For thick surfaces
		if (seperation == 6)
		{
			// Compute a 'normal' vector to the appropriate corner.
			// This roughly corresponds to 'K' in Figure 10 of the paper.
			Vector3f cornerNormal;
			cornerNormal[0] = planeNormal[0] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal[1] = planeNormal[1] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal[2] = planeNormal[2] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal = normalize(cornerNormal);

			float cosAlpha = dot(cornerNormal, planeNormal) / (length(cornerNormal) * length(planeNormal));

			return (L / 2) * sqrt(3.0f) * cosAlpha;
		}

		// For thin surfaces.
		if (seperation == 26)
		{
			// Compute the normal to the face which is best
			// aligned (most co-planar) with the input triangle;
			Vector3f faceNormal(0.0f);
			auto maxIndex = abs(planeNormal).maxComponentIndex();
			faceNormal[maxIndex] = planeNormal[maxIndex] >= 0.0f ? 1.0f : -1.0f;

			const float cosBeta = dot(faceNormal, planeNormal) / (length(faceNormal) * length(planeNormal));

			// The paper says "For generating 6-separating line we define t6 = (L/2)cos[Beta]". I think
			// they  mean 'surface' instead of 'line', to match the previous sentence in the paper?
			return (L / 2) * cosBeta;
		}

		assert(false && "Seperation must be '6' or '26'");
		return 0;
	}

	void drawFragment(int32 x, int32 y, int32 z, MaterialId matId, MaterialId externalMaterial, Volume& volume, uint pass)
	{
		const MaterialId checkerboard = ((x & 0x1) ^ (y & 0x1) ^ (z & 0x1)) + 1; // '1' or '2'.

		switch (pass)
		{
		case 0:
			// Write a checkerboard of values which are different from the exteriour material.
			volume.setVoxel(x, y, z, externalMaterial + checkerboard);
			break;
		case 1:
			// Draw surface only for interiour voxels.
			if (volume.voxel(x, y, z) != externalMaterial)
			{
				volume.setVoxel(x, y, z, matId);
			}
			break;
		default:
			assert(false && "Invalid pass");
		}
	}

	void scanConvert3D(const Triangle& triangle, MaterialId matId, MaterialId externalMaterial, Volume& volume, Thickness thickness, float multiplier, uint pass)
	{
		// Now that we have the threshold we iterate over all voxels in the
		// triangle's bounding box and output those which are close enough to it.

		const Vector3f a = triangle.vertices[0];
		const Vector3f b = triangle.vertices[1];
		const Vector3f c = triangle.vertices[2];

		// Compute threshold based on Section 4.2 of "An Accurate Method To Voxelize Polygonal Meshes" by Huang et al.
		Vector3f planeNormal = cross(b - a, a - c);

		// Not sure if/how we should handle tiny triangles?
		//assert(length(planeNormal) > 0.01f);
		if (length(planeNormal) <= 0.01f)
		{
			return;
		}

		float distThreshold = computeMinThickness(triangle);

		distThreshold *= multiplier;

		distThreshold += 0.01f;


		// Uncomment the line below to manually specify the threshold.
		// distThreshold = sqrt(3.0f) / 2.0f; // Half the diagonal of a voxel

		// FIXME - Handle rounding
		Vector3i boundsMin = static_cast<Vector3i>(min(a, min(b, c))); // Round down
		Vector3i boundsMax = static_cast<Vector3i>(max(a, max(b, c)) + Vector3f(1.0f)); // Round up

		// Expand
		boundsMin -= Vector3i(distThreshold + 0.5f);
		boundsMax += Vector3i(distThreshold + 0.5f);

		for (int32_t z = boundsMin.z(); z <= boundsMax.z(); z++)
		{
			for (int32_t y = boundsMin.y(); y <= boundsMax.y(); y++)
			{
				for (int32_t x = boundsMin.x(); x <= boundsMax.x(); x++)
				{
					// Using the distance here is slightly incuurate. We should actually test if
					// the point is inside the (very thin) trianglular prism formed by extending
					// slightly along the triangle normal in both directions. The distance function
					// is similar but also has rounded edges, meaning it generates fragments which
					// are slightly outside the triangle. This doesn't matter for our case as later
					// we run all fragments through our winding number test, and it may even help
					// fix and cracks between triangles.
					float dist = distance(Vector3f(x, y, z), Triangle(a, b, c));

					const float epsilon = 0.001f;
					if (dist <= (distThreshold + epsilon))
					{
						drawFragment(x, y, z, matId, externalMaterial, volume, pass);
					}
				}
			}
		}
	}

	void scanConvert3DRecursive(const Triangle& triangle, MaterialId matId, MaterialId externalMaterial, Volume& volume, Thickness thickness, float multiplier, uint pass)
	{
		const float splitThreshold = 20.0f;

		uint32_t longestSideIndex = 0;
		if (triangle.sideLength(1) > triangle.sideLength(0)) { longestSideIndex = 1; }
		if (triangle.sideLength(2) > triangle.sideLength(1)) { longestSideIndex = 2; }

		float alignmentThreshold = 0.9f;
		Vector3f absNormal = abs(triangle.computeNormal());
		bool isAxisAligned = absNormal.x() > alignmentThreshold || absNormal.y() > alignmentThreshold || absNormal.z() > alignmentThreshold;

		// If a triangle is large and not aligned to an axis then its bounding box will overlap a
		// lot of voxels. This makes it slow to rasterize, therefore we split such triangles in half.
		if ((!isAxisAligned) && (triangle.sideLength(longestSideIndex) > splitThreshold))
		{
			Vector3f longestSideMidpoint = (triangle.vertices[longestSideIndex] + triangle.vertices[(longestSideIndex + 1) % 3]) / 2.0f;

			Triangle output0 = triangle;
			Triangle output1 = triangle;

			output0.vertices[longestSideIndex] = longestSideMidpoint;
			output1.vertices[(longestSideIndex + 1) % 3] = longestSideMidpoint;

			scanConvert3DRecursive(output0, matId, externalMaterial, volume, thickness, multiplier, pass);
			scanConvert3DRecursive(output1, matId, externalMaterial, volume, thickness, multiplier, pass);
		}
		else
		{
			scanConvert3D(triangle, matId, externalMaterial, volume, thickness, multiplier, pass);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	//   _____ _   _                     _          __  __                                        //
	//  |  _  | | | |                   | |        / _|/ _|                                       //
	//  | | | | |_| |__   ___ _ __   ___| |_ _   _| |_| |_                                        //
	//  | | | | __| '_ \ / _ \ '__| / __| __| | | |  _|  _|                                       //
	//  \ \_/ / |_| | | |  __/ |    \__ \ |_| |_| | | | |                                         //
	//   \___/ \__|_| |_|\___|_|    |___/\__|\__,_|_| |_|                                         //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////

	void voxelizeShell(Volume& volume, const TriangleList& triList, const MaterialList& materialList, MaterialId externalMaterial, uint pass, ProgressMonitor* progMon)
	{
		uint32_t triCount = triList.size();

		uint32_t trianglesDrawn = 0;
		if (progMon) { progMon->startTask("Performing 3D scan conversion"); }

		assert(triList.size() == materialList.size());

		//for (const auto triangle : triList)
		for (uint i = 0; i < triList.size(); i++)
		{
			const Triangle& triangle = triList[i];
			const MaterialId& matId = materialList[i];
			if ((trianglesDrawn % 1000 == 0))
			{
				if (progMon) { progMon->setProgress(0, trianglesDrawn, triCount); };
			}

			/*if (index % 300 == 0)
			{
			std::cout << "Writing triangle... " << index / 3 << std::endl;
			}*/

			//uint32_t i0 = mesh.triangles()[index].indices[0];
			//uint32_t i1 = mesh.triangles()[index].indices[1];
			//uint32_t i2 = mesh.triangles()[index].indices[2];

			Triangle scaled(triangle);
			const float scaleFactor = 1.01f;
			Vector3f centre = (scaled.vertices[0] + scaled.vertices[1] + scaled.vertices[2]) / 3.0f;

			scaled.translate(centre * -1.0f);
			scaled.scale(scaleFactor);
			scaled.translate(centre);

			// We want to draw the thinnest triangles we can (while still seperating nodes on
			// one side from nodes on the other) in order to minimise the number of nodes we
			// have to classify. Intuitively this should mean 26-seperating triangles. But this
			// is not enough because some (roughly half) of the generated 'fragments' actually 
			// fall outside of the object. Once we remove these during the subsequent inside/outside
			// test we end up with holes in our surface. I believe this is a problem both for
			// seperating the inside vs. outside, and also if we later colour the surface (as we
			// need to make sure underlying voxels don't show through).
			//
			// Therefore we need to draw a thicker surface. It can be as thick as we like
			// (and perhaps we should let the user configure it) but must be at least enough
			// to avoid the aforementioned holes. I don't have a mathematical derivation for
			// this but have simply found that scaling the thickness by the multiplier below
			// gives visually satisfactory results. Any less and we start to see holes.
			//
			// We could also use 6-seperating triangles in this case a multiplier is still
			// needed, but the slightly smaller value of 1.9 seemed to work), but overall I
			// prefered the qualty and consistancy of the 26-seperating triangles.
			float multiplier = pass == 0 ? 1.0f : 2.1f;
			scanConvert3DRecursive(scaled, matId, externalMaterial, volume, Thickness::Separate26, multiplier, pass);
			trianglesDrawn++;
		}

		if (progMon) { progMon->finishTask(); }
	}

	bool isInside(const Vector3f& queryPoint, const Surface& surface, bool useHierarchicalEval)
	{
		float sum = 0.0f;
		if (useHierarchicalEval)
		{
			sum = computeWindingNumber(queryPoint, *(surface.closedTriangleTree));
		}
		else
		{
			sum = computeWindingNumber(queryPoint, surface.triangles);
		}
		return std::abs(sum) > gWindingNumberThreshold;
	}

	void doBruteForceVoxelisation(Volume& volume, Surface& surface, MaterialId internalMaterial,
		MaterialId externalMaterial, bool useHierarchicalMeshEval, ProgressMonitor* progMon)
	{
		Box3i voxelisationBounds = static_cast<Box3i>(surface.bounds);
		voxelisationBounds.dilate(2); // Make sure we really cover everything

		Vector3i minBound = voxelisationBounds.lower();
		Vector3i maxBound = voxelisationBounds.upper();

		// Iterate over each voxel within the bounds and classify as inside or outside
		if (progMon) { progMon->startTask("Classifying voxels (Brute Force!)"); }
		for (int32 volZ = minBound.z(); volZ <= maxBound.z(); volZ++)
		{
			if (progMon){ progMon->setProgress(minBound.z(), volZ, maxBound.z()); }
			for (int32 volY = minBound.y(); volY <= maxBound.y(); volY++)
			{
				for (int32 volX = minBound.x(); volX <= maxBound.x(); volX++)
				{
					if (isInside(Vector3f(volX, volY, volZ), surface, useHierarchicalMeshEval))
					{
						volume.setVoxel(volX, volY, volZ, internalMaterial);
					}
					else
					{
						volume.setVoxel(volX, volY, volZ, externalMaterial);
					}
				}
			}			
		}
		if (progMon) { progMon->finishTask(); }
	}

	struct NodeToTest
	{
		uint32_t index;
		uint32_t childId;
		Vector3f centre;
		bool result;
	};

	class NodeFinder
	{
	public:
		bool operator()(NodeDAG& nodes, uint32 nodeIndex, Box3i nodeBounds)
		{
			// Signal to stop traversing parts of the tree which do not overlap our voxelised object.
			if (!overlaps(nodeBounds, mBounds)) { return false; }

			// If this is a non-leaf (internal) node then check if any of its children are leaves.
			// If so, add them to the list of nodes on which the inside/outside test will be performed.
			if (!isMaterialNode(nodeIndex))
			{
				for (unsigned int childId = 0; childId < 8; childId++)
				{
					uint32 childNodeIndex = nodes[nodeIndex][childId];
					if (isMaterialNode(childNodeIndex))
					{
						Box3i childNodeBounds = childBounds(nodeBounds, childId);
						if (overlaps(childNodeBounds, mBounds))
						{
							NodeToTest toTest;
							toTest.index = nodeIndex;
							toTest.childId = childId;
							toTest.centre = static_cast<Vector3f>(childNodeBounds.lower() + childNodeBounds.upper()) * 0.5f;
							mNodes.push_back(toTest);
						}
					}
				}
			}

			return true;
		}

		std::vector<NodeToTest> nodes()
		{
			return mNodes;
		}

	private:
		std::vector<NodeToTest> mNodes;

	public:
		Box3i mBounds;

	};

	std::vector<NodeToTest> findNodes(Volume& volume, Box3i bounds)
	{
		NodeFinder nodeFinder;
		nodeFinder.mBounds = bounds;
		traverseNodesRecursive(volume, nodeFinder);
		return nodeFinder.nodes();
	}

	void classifyNodes(std::vector<NodeToTest>& nodesToTest, NodeDAG& /*nodeData*/, Surface& surface, ProgressMonitor* progMon, bool useHierarchicalMeshEval)
	{
		int i = 0;
		std::stringstream ss;
		ss << "Classifying " << nodesToTest.size() << " nodes";
		if (progMon) { progMon->startTask(ss.str()); }

		std::mutex m;

#if defined (_MSC_VER) && _MSC_VER >= 1900
		std::for_each(std::execution::par, nodesToTest.begin(), nodesToTest.end(), [&](auto&& nodeToTest)
#else
		std::for_each(nodesToTest.begin(), nodesToTest.end(), [&](auto&& nodeToTest)
#endif
		{
			nodeToTest.result = isInside(nodeToTest.centre, surface, useHierarchicalMeshEval);
			std::lock_guard<std::mutex> guard(m);

			if (progMon && i % 1000 == 0) { progMon->setProgress(0, i, nodesToTest.size()); }
			i++;
		});


		if (progMon) { progMon->finishTask(); }
	}

	uint32 prune(NodeDAG& nodes, uint32 nodeIndex)
	{
		Node node = nodes[nodeIndex];
		for (int i = 0; i < 8; i++)
		{
			uint32 nodeChildIndex = node[i];
			if (!isMaterialNode(nodeChildIndex))
			{
				node[i] = prune(nodes, nodeChildIndex);
			}
		}

		return nodes.isPrunable(node) ? node[0] : nodes.insert(node);
	}

	void prune(Volume& volume)
	{
		uint32 rootNodeIndex = getRootNodeIndex(volume);
		NodeDAG& nodes = getNodes(volume);
		uint32 newRootNodeIndex = prune(nodes, rootNodeIndex);
		volume.setRootNodeIndex(newRootNodeIndex);
	}

	void doHierarchicalVoxelisation(Volume& volume, Surface& surface, MaterialId internalMaterial,
		MaterialId externalMaterial, bool useHierarchicalMeshEval, ProgressMonitor* progMon)
	{
		Box3i voxelisationBounds = computeBounds(volume, externalMaterial);

		// Find all the nodes - both the single-voxel nodes which are part of the
		// voxelised surface and (hopefully larger) nodes which are either side of it.
		auto nodesToTest = findNodes(volume, voxelisationBounds);

		// Classify all nodes according to which side of the surface they are on.
		NodeDAG& nodeData = Internals::getNodes(volume);
		classifyNodes(nodesToTest, nodeData, surface, progMon, useHierarchicalMeshEval);

		// Apply the results to the volume
		for (auto& toTest : nodesToTest)
		{
			MaterialId resultingMaterial = toTest.result ? internalMaterial : externalMaterial;
			nodeData.nodes().setNodeChild(toTest.index, toTest.childId, resultingMaterial);
		}

		// The voxelisation process can cause the volume to become unpruned, which we consider to be an invalid state.
		// This might be because the shell voxelisaion is too thick, though I think we have avoided that. But even so,
		// the mesh might represent a small box touching eight voxels, all of which are inside. These would be
		// individually classified and would all get the same value, so they should be pruned and replaced by the parent.
		prune(volume);
	}

	void voxelize(Volume& volume, Surface& surface, bool useSurfaceMaterials,
		MaterialId* internalMaterialOverride, MaterialId* externalMaterialOverride,
		ProgressMonitor* progMon)
	{
		const bool useHierarchicalMeshEval = true;
		// Perform the classification per-octree-node instead of per-voxel.This is
		// much faster and for a well formed mesh the result should be the same, but it 
		// may have problems on scenes with open meshes, ground planes, windows, etc.
		const bool useHierarchicalNodeEval = true;

		MaterialId internalMaterial = internalMaterialOverride ? *internalMaterialOverride : surface.mainMaterial;
		MaterialId externalMaterial = externalMaterialOverride ? *externalMaterialOverride : 0;

		// A negative mean winding number indcates that the mesh is inside-out.
		if (surface.meanWindingNumber < 0.0f) { std::swap(internalMaterial, externalMaterial); }

		volume.fill(externalMaterial);

		voxelizeShell(volume, surface.triangles, surface.materials, externalMaterial, 0, progMon);

		if(std::abs(surface.meanWindingNumber) > 0.9f)
		{
			if (useHierarchicalNodeEval)
			{
				doHierarchicalVoxelisation(volume, surface, internalMaterial, externalMaterial, useHierarchicalMeshEval, progMon);
			}
			else
			{
				doBruteForceVoxelisation(volume, surface, internalMaterial, externalMaterial, useHierarchicalMeshEval, progMon);
			}
		}
		else
		{
			std::cout << "Warning: Skipping fill for badly-formed surface \'" << surface.name << "\'" << std::endl;
		}

		if (useSurfaceMaterials)
		{
			// Write the surface voxels with their correct materials as specified in the mesh. 
			voxelizeShell(volume, surface.triangles, surface.materials, externalMaterial, 1, progMon);
		}
	}
}