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

#include <cassert>
#include <cmath>
#include <thread>
#include <unordered_map>

#ifdef CUBIQUITY_USE_AVX
	#include <immintrin.h>
#endif

namespace Cubiquity
{
	using namespace Internals;

	static const int gCheckerBoardBitShift = 15;
	static const MaterialId gCheckerBoardBit = 1 << gCheckerBoardBitShift;

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
	void scanConvert3D(const Triangle& triangle, MaterialId matId, Volume& volume, Thickness thickness, float multiplier)
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

		const float L = 1.0f; // From figure 10 in the paper, our voxels have unity side length.
		float distThreshold = 0.0f;

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

		// For thick surfaces
		if (thickness == Thickness::Separate6)
		{
			// Compute a 'normal' vector to the appropriate corner.
			// This roughly corresponds to 'K' in Figure 10 of the paper.
			Vector3f cornerNormal;
			cornerNormal[0] = planeNormal[0] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal[1] = planeNormal[1] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal[2] = planeNormal[2] >= 0.0f ? 1.0f : -1.0f;
			cornerNormal = normalize(cornerNormal);

			float cosAlpha = dot(cornerNormal, planeNormal) / (length(cornerNormal) * length(planeNormal));

			distThreshold = (L / 2) * sqrt(3.0f) * cosAlpha;
		}

		// For thin surfaces.
		if (thickness == Thickness::Separate26)
		{
			// Compute the normal to the face which is best
			// aligned (most co-planar) with the input triangle;
			Vector3f faceNormal(0.0f);
			auto maxIndex = abs(planeNormal).maxComponentIndex();
			faceNormal[maxIndex] = planeNormal[maxIndex] >= 0.0f ? 1.0f : -1.0f;

			const float cosBeta = dot(faceNormal, planeNormal) / (length(faceNormal) * length(planeNormal));

			// The paper says "For generating 6-separating line we define t6 = (L/2)cos[Beta]". I think
			// they  mean 'surface' instead of 'line', to match the previous sentence in the paper?
			distThreshold = (L / 2) * cosBeta;
		}

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
						// Logical OR the mateiral id with a checkerboard pattern, in order to  
						// prevent neighbouring voxels from being identical and therefore pruned.
						MaterialId checkerboard = (x & 0x1) ^ (y & 0x1) ^ (z & 0x1);
						checkerboard <<= gCheckerBoardBitShift;
						volume.setVoxel(x, y, z, matId | checkerboard);
					}
				}
			}
		}
	}

	void scanConvert3DRecursive(const Triangle& triangle, MaterialId matId, Volume& volume, Thickness thickness, float multiplier)
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

			scanConvert3DRecursive(output0, matId, volume, thickness, multiplier);
			scanConvert3DRecursive(output1, matId, volume, thickness, multiplier);
		}
		else
		{
			scanConvert3D(triangle, matId, volume, thickness, multiplier);
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

	// Finds the main material for an object by looking for the subobject with the biggest
	// winding number with respect to the object centre, and returning the material from that.
	MaterialId findMainMaterial(const Object& object)
	{
		Box3f bounds = computeBounds(object);
		Vector3f queryPoint = bounds.centre();

		MaterialId materialId = 0;
		float biggestWindingNumber = std::numeric_limits<float>::lowest();

		for (const auto& subObject : object.subObjects)
		{
			float windingNumber = computeWindingNumber(queryPoint, subObject.second);
			if (windingNumber > biggestWindingNumber)
			{
				biggestWindingNumber = windingNumber;
				materialId = subObject.first;
			}
		}

		return materialId;
	}

	float threshold = 0.0f;

	MaterialId determineMaterial(const Vector3f& queryPoint, ClosedTriangleTreeList& triangleMap)
	{
		float sum = 0.0f;
		float biggestWindingNumber = std::numeric_limits<float>::lowest();
		MaterialId biggestMaterial = 0x0;

		for (auto& triangleContainer : triangleMap)
		{
			float windingNumber = computeWindingNumber(queryPoint, triangleContainer.second);

			sum += windingNumber;

			// It is quite common for a voxel to be inside two closed objects, meaning it could get the material.
			// from either. The epsilon below causes it to choose the first unless the second is signifcantly larger.
			float epsilon = 0.5f;
			if (windingNumber > biggestWindingNumber + epsilon)
			{
				biggestWindingNumber = windingNumber;
				biggestMaterial = triangleContainer.first;
			}
		}

		return sum > threshold ? biggestMaterial : 0;
	}

	void fillInsideMeshBruteForce(Volume& volume, ClosedTriangleTreeList& closedTriangleTreeList, Box3f bounds, bool preserveSurfaceMaterials, ProgressMonitor* progMon)
	{
		// In the current implementation the bouds we receive are actually for the mesh, but some 
		// dilation occurs due to the threshold used in the 3D scan conversion. Therefore we need
		// to dilate the bounds slightly to ensure we capture all voxels set during this process.
		bounds.dilate(3.0f);

		Vector3f minBound = bounds.lower();
		Vector3f maxBound = bounds.upper();

		//ProgressBar progressBar(maxBound.z() - minBound.z(), 70);

		if (progMon) { progMon->startTask("Classifying nodes (Brute Force!)"); }
		for (int32 volZ = minBound.z(); volZ <= maxBound.z(); volZ++)
		{
			if (progMon){ progMon->setProgress(minBound.z(), volZ, maxBound.z()); }
			for (int32 volY = minBound.y(); volY <= maxBound.y(); volY++)
			{
				for (int32 volX = minBound.x(); volX <= maxBound.x(); volX++)
				{
					// Read the current material, while also clearing the checkerboard pattern.
					MaterialId matId = volume.voxel(volX, volY, volZ) & (~gCheckerBoardBit);

					// Determine the suggested new material
					MaterialId newMatId = determineMaterial(
						Vector3f(volX, volY, volZ), closedTriangleTreeList);

					// This condition preserves the materials generated by scan conversion. When
					// preserving, we can only write to a voxel if it is currently empty or if it
					// is being deleted. We cannot replace one valid material with another.
					//
					// This flag might cause differences at material boundaries, e.g. if two 
					// differnt polygons affect a voxel then with this condition the result will be
					// whatever was drawn last, rather than which has the biggest winding number.
					if ((!preserveSurfaceMaterials) || (matId == 0) || (newMatId == 0))
					{
						matId = newMatId;
					}

					// Write back result (clearing the checkerboard if nothing else).
					volume.setVoxel(volX, volY, volZ, matId);
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
		MaterialId result; // FIXME - Should this actually be a bool?
	};

	class NodeFinder
	{
	public:
		void operator()(NodeDAG& nodes, uint32 nodeIndex, Box3i bounds)
		{
			if (!isMaterialNode(nodeIndex))
			{
				// Octree should not be merged here, so all ref counts should be '1' (except for the roor, which should always be zero).

				//FIXME - Used to use ref count, wht to do here?
				//const uint32 expectedRefCount = nodeIndex == RootNodeIndex ? 0 : 1;
				//(void)(expectedRefCount); // Suppress 'unused' warnings when not using asserts.
				//assert(nodes[nodeIndex].mRefCount == expectedRefCount && "Octree cannot be merged for this operation.");

				uint32 nodeSideLength = bounds.upper().x() - bounds.lower().x() + 1;
				uint32 childSideLength = nodeSideLength / 2;
				for (unsigned int childId = 0; childId < 8; childId++)
				{
					uint32 childIndex = nodes[nodeIndex][childId];
					if (isMaterialNode(childIndex))
					{
						uint32_t childX = (childId >> 0) & 0x01;
						uint32_t childY = (childId >> 1) & 0x01;
						uint32_t childZ = (childId >> 2) & 0x01;

						Vector3i childLowerCorner;

						// childX/Y/Z are all zero or one.
						childLowerCorner[0] = bounds.lower()[0] + (childSideLength * childX);
						childLowerCorner[1] = bounds.lower()[1] + (childSideLength * childY);
						childLowerCorner[2] = bounds.lower()[2] + (childSideLength * childZ);

						Vector3i childUpperCorner = childLowerCorner + Vector3i(childSideLength - 1);

						NodeToTest toTest;
						toTest.index = nodeIndex;
						toTest.childId = childId;
						toTest.centre = static_cast<Vector3f>(childUpperCorner + childLowerCorner) * 0.5f;
						mNodes.push_back(toTest);
					}
				}
			}
		}

		std::vector<NodeToTest> nodes()
		{
			return mNodes;
		}

	private:
		std::vector<NodeToTest> mNodes;

	};

	std::vector<NodeToTest> findNodes(Volume& volume)
	{
		NodeFinder nodeFinder;
		traverseNodes(volume, nodeFinder);
		return nodeFinder.nodes();
	}


	void classifyNodesWorker(std::vector<NodeToTest>& nodesToTest, std::list<std::pair<MaterialId, ClosedTriangleTree> >& nodes, int threadId, int threadCount, ProgressMonitor* progMon)
	{
		if (progMon) { progMon->startTask("Classifying nodes"); }
		for (uint32_t index = threadId; index < nodesToTest.size(); index += threadCount)
		{
			// FIXME - index % 1000 might never be hit as we increament by threadCount.
			if (threadId == 0 && (index % 1000 == 0))
			{
				if (progMon) { progMon->setProgress(threadId, index, nodesToTest.size()); }
			}

			auto& toTest = nodesToTest[index];

			toTest.result = determineMaterial(toTest.centre, nodes);
		}
		if (threadId == 0 && progMon) { progMon->finishTask(); }
	}

	void classifyNodes(std::vector<NodeToTest>& nodesToTest, NodeDAG& /*nodeData*/, ClosedTriangleTreeList& closedTriangleTreeList, ProgressMonitor* progMon)
	{
		const int threadCount = 4;
		std::vector<std::thread> threads;
		for (int i = 0; i < threadCount; i++)
		{
			threads.push_back(std::thread(classifyNodesWorker, std::ref(nodesToTest), std::ref(closedTriangleTreeList), i, threadCount, progMon));
		}

		for (auto& thread : threads)
		{
			thread.join();
		}
	}

	void voxelizeShell(Volume& volume, const Geometry& splitTriangles, bool preserveSurfaceMaterials, ProgressMonitor* progMon)
	{
		uint32_t triCount = 0;
		for (auto& object : splitTriangles)
		{
			for (auto& subObject : object.subObjects)
			{
				triCount += subObject.second.size();
			}
		}

		std::vector<NodeToTest> boundaryNodesToTest;

		uint32_t trianglesDrawn = 0;
		if (progMon) { progMon->startTask("Performing 3D scan conversion"); }
		for (auto& object : splitTriangles)
		{
			for (auto& subObject : object.subObjects)
			{
				for (const auto triangle : subObject.second)
				{
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

					MaterialId materialId = subObject.first;

					// When *not* preserving surface materials we draw the thinnest triangles we
					// can (while still seperating nodes on one side from nodes on the other)
					// which is are 26-seperating triangles. But this is not enough if we are
					// trying to preserve the surface material, because some (roughly half) of the
					// generated 'fragments' actually fall outside of the object. Once we remove
					// these during the subsequent inside/outside test we end up with holes in our
					// surface.
					//
					// Therefore we need to draw a thicker surface. It can be as thick as we like
					// (and perhaps we should let the user configure it) but must be at least enough
					// to avoid the aforementioned holes. I don't have a mathematical derivation for
					// this but have simply found that scaling the thickness by the multiplier below
					// gives visually satisfactory results. Any less and we start to see holes.
					//
					// We could also use 6-seperating triangles in this case  a multiplier is still
					// needed, but the slightly smaller value of 1.9 seemed to work), but overall I
					// prefered the qualty and consistancy of the 26-seperating triangles.
					float multiplier = preserveSurfaceMaterials ? 2.1f : 1.0f;
					scanConvert3DRecursive(scaled, materialId, volume, Thickness::Separate26, multiplier);
					trianglesDrawn++;
				}
			}
		}
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

	void fillInsideMeshNodeBased(Volume& volume, ClosedTriangleTreeList& closedTriangleTreeList, bool preserveSurfaceMaterials, ProgressMonitor* progMon)
	{
		// The volume starts off full of zeros before the boudaries are written. A boundary should never be written.
		// as zeros. Therefore any zero nodes at this point do not correspond to boundaries, and should be classified.
		auto nodesToTest = findNodes(volume);

		NodeDAG& nodeData = Internals::getNodes(volume);

		std::cout << "Classifying " << nodesToTest.size() << " nodes..." << std::endl;
		classifyNodes(nodesToTest, nodeData, closedTriangleTreeList, progMon);
		// Apply the results to the volume
		for (auto& toTest : nodesToTest)
		{
			Node myNode = nodeData[toTest.index];

			// Read the current material, while also clearing the checkerboard pattern.
			MaterialId matId = myNode[toTest.childId] & (~gCheckerBoardBit);

			// Determine the suggested new material
			MaterialId newMatId = toTest.result;

			// This condition preserves the materials generated by scan conversion. When
			// preserving, we can only write to a voxel if it is currently empty or if it
			// is being deleted. We cannot replace one valid material with another.
			//
			// This flag might cause differences at material boundaries, e.g. if two 
			// differnt polygons affect a voxel then with this condition the result will be
			// whatever was drawn last, rather than which has the biggest winding number.
			if ((!preserveSurfaceMaterials) || (matId == 0) || (newMatId == 0))
			{
				matId = newMatId;
			}

			// Write back result (clearing the checkerboard if nothing else).
			nodeData.nodes().setNodeChild(toTest.index, toTest.childId, matId);
		}
	}

	float findThreshold(const Geometry& geometry)
	{
		// Determine whether the threshold should be positive or negative (+0.5 or -0.5).
		//
		// A given mesh may represent a solid object in the world or it may represent a space to be
		// carved out of the world. Assuming the triangle windings are correctly and consistently
		// defined, a solid object should have winding number of +1.0 on the inside with the
		// surrounding winding number of 0.0 representing empty space, while a hollow object would have
		// a winding number of -1.0 in the inside with the surrounding winding number of 0.0
		// representing solid material. In other words, a winding number of 0.0 can represent either
		// empty of solid space depending on the nature of the scene.
		//
		// A typical scene will contain multiple solid object which may be intersecting or placed
		// inside of eachother (for example to give the centre of an object a different material to the
		// rest of it). These solid objects may be placed inside an object representing empty space,
		// such the surrounding room. It is unlikely that multiple objects representing empty space
		// will intersect of be placed inside eachother, because this does not make sense from a
		// physical point of view. You can't carve empty space out of other emty space. Or from another
		// point of view, Cubiquity only has one representation of empty space (a material id of zero)
		// but many representation of solid space (any other material id). 
		//
		// We therefore attempt to determine the threshold automatically by sampling numerous random points, some of which will
		// hopefully be inside the mesh and some of which will be outside. Ideally we simply check the sign of all non-
		// zero values, find they are consistant and that the magnitude is an integer, and we have our answer. In practice, many
		// meshes are not well formed due to holes and inconsistant windings, and so we get a range of generalised winding
		// numbers.
		//
		// Therefore we need to do a bit of statistical analysis on the results. At a basic level we could just look at the
		// sign of the average or count the number of positive vs. negative results. More advanced approaches could explore
		// the distribution further but we'll leave that until we encounter a need for it.
		std::vector<float> samples;

		Cubiquity::Box3f inputBounds = computeBounds(geometry);

		Cubiquity::Box3fSampler sampler(inputBounds);

		// Generate the samples.
		const int sampleCount = 100;
		for (uint32_t i = 0; i < sampleCount; i++)
		{
			auto v = sampler.next();
			float windingNumber = 0.0f;
			for (const auto& object : geometry)
			{
				for (const auto& subobject : object.subObjects)
				{
					windingNumber += computeWindingNumber(v, subobject.second);
				}
			}
			samples.push_back(windingNumber);
		}

		// Compute some basic statistics about them.
		uint32_t positiveCount = 0;
		uint32_t negativeCount = 0;
		for (auto sample : samples)
		{
			const float inclusionThreshold = 0.5f; // No real basis for this value. It can be tweaked as required.

			// I believe that samples with a samll magnitude are more likely to be errornous (unless they are zero, but
			// those don't tell us anything about the winding anyway?). Let's focus on analysing the bigger samples.
			if (std::abs(sample) > inclusionThreshold)
			{
				sample > 0.0f ? positiveCount++ : negativeCount++;
			}
		}

		std::cout << "Positive count = " << positiveCount << std::endl;
		std::cout << "Negative count = " << negativeCount << std::endl;

		// For a good quality mesh we expect that all non-zero samples have the same sign. Warn the user if this is not the
		// case. Note that this won't detect all poor-quality meshes, it depends on the sampling and inclusion threshold.
		if (positiveCount > 0 && negativeCount > 0)
		{
			std::cout << "WARNING - Mixed face orientations in mesh?" << std::endl;
		}

		// If the mesh is watertight then then all generalised winding number should be integers.
		for (auto sample : samples)
		{
			if (fabsf(roundf(sample) - sample) > 0.001f)
			{
				std::cout << "WARNING - Non-integer sample found. Is the mesh open or with missing faces?" << std::endl;
				break;
			}
		}

		// Theorectcally the thresholds should be +/-0.5f. However, triangles often pass *exactly*
		// through the centre of a voxel due to e.g. being modelled on a grid in a DCC package,
		// and/or being scaled by interger amounts. Rounding errors mean that such triangles can
		// end up with a noisy surface, and the epsilon below helps to avoid this.
		float epsilon = 0.001f;
		float threshold = 0.5f + epsilon;
		return positiveCount >= negativeCount ? threshold : -threshold;

		//if (negativeCount > positiveCount)
		//{
		//	std::cout << "WARNING - Flipping triangles" << std::endl;
		//	for (auto& triangle : triangles)
		//	{
		//		triangle.flip();
		//	}
		//}
	}

	void voxelize(Volume& volume, Geometry& splitTriangles, bool preserveSurfaceMaterials,
		uint16 internalMaterialOveride,	ProgressMonitor* progMon, bool useBruteForce)
	{
		threshold = findThreshold(splitTriangles);

		// Voxelise the shell even in brute force mode to support surface material preservation.
		voxelizeShell(volume, splitTriangles, preserveSurfaceMaterials, progMon);

		// We may now merge objects and/or sub-objects before doing the inside vs. outside testing.
		// If the user has requested that a single material be used to fill all objects then we no
		// longer care about the individual objects and they can be merged into one. It's not yet
		// clear if this makes the process fastr - on one hand the original objects are quite
		// possibly closed or require less closing triangles (which makes for a nice grouping), but
		// on the other hand there may be a lot (and they have no hierarchy) so they may lose out
		// to a single ClosedTriangleTree culling more aggressively.
		Geometry mergedGeometry;
		if (internalMaterialOveride)
		{
			Object mergedObject = mergeObjects(splitTriangles, "Merged Object", internalMaterialOveride);
			mergedGeometry.push_back(mergedObject);
		}
		else
		{

			// A single object in a scene is often closed, but if it consists of multiple subobjects
			// then there is a good chance that each of those is more open. Triangle lists which are 
			// 'more closed' need less closing triagles and so are faster to process, and there will 
			// also be less triangle lists if subobjects have been merged. However, it may result in
			// a less optimal choice of material allocation because we have to choose a single
			// material for the whole interior of the object. It is not yet clear what the best
			// approach is, so for now it can be toggled via the constant below.
			const bool doMergeSubObjects = true;
			if (doMergeSubObjects)
			{
				for (auto& inputObject : splitTriangles)
				{
					MaterialId mainMaterialId = findMainMaterial(inputObject);
					SubObject mergedSubObject = mergeSubObjects(inputObject.subObjects, mainMaterialId);

					Object outputObject;
					outputObject.subObjects.push_back(mergedSubObject);

					mergedGeometry.push_back(outputObject);
				}
			}
			else
			{
				mergedGeometry = splitTriangles;
			}
		}

		// We now build a list of ClosedTriangleTrees with associated materials. There will be one
		// entry for each object or sub-object in the scene (which might mean one entry in total if
		// objects/sub-objects were merged above). These are used for the point-in-mesh tests.
		ClosedTriangleTreeList closedTriangleTreeList;
		for (auto& object : mergedGeometry)
		{
			for (auto& subObject : object.subObjects)
			{
				closedTriangleTreeList.push_back(std::make_pair(subObject.first, ClosedTriangleTree(subObject.second)));
			}
		}

		// Brute-force mode is only really useful for testing and debugging purposes. It applies
		// the inside/outside test directly on each voxel, rather than on the nodes resulting
		// from dividing up the volume during scan conversion.
        if (useBruteForce)
        {
			Cubiquity::Box3f bounds = computeBounds(splitTriangles);
			fillInsideMeshBruteForce(volume, closedTriangleTreeList, bounds, preserveSurfaceMaterials, progMon);
        }
        else
        {
            fillInsideMeshNodeBased(volume, closedTriangleTreeList, preserveSurfaceMaterials, progMon);

            // The voxelisation process can cause the volume to become unpruned, which we consider to be an invalid state. I
            // suspect it happens mostly when the voxelised shell thickness is to great, though perhaps it also happens with
            // difficult geometry? Taking the former case, we might submit too many nodes for classification and a group of
            // eight children might all fall the same side of a boundary and so get the same value. These are then unpruned.
			prune(volume);
			//volume.validate();

            // Sanity check that there are no unpruned nodes. Should be moved into utility/test function?

			//FIXME - Used to use ref count, wht to do here?
            /*NodeDAG& nodeData = Internals::getNodes(volume);
            for (uint32 index = RootNodeIndex; index < nodeData.size(); index++)
            {
                const Node& node = nodeData[index];
                if (node.refCount() > 1)
                {
                    std::cerr << "Unexpected ref count!" << std::endl;
                }

                //if (node.refCount() > 0)
                {
                    if (nodeData.isMaterialNode(node.child(0)))
                    {
                        bool allChildrenMatch = true;
                        for (uint32 childId = 1; childId < 8; childId++)
                        {
                            if (node.child(0) != node.child(childId))
                            {
                                allChildrenMatch = false;
                            }
                        }

                        if (allChildrenMatch)
                        {
                            std::cerr << "All children match! " << index << ", " << node.child(0) << std::endl;
                        }
                    }
                }
            }*/
        }
	}
}