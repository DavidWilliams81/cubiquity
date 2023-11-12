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
#include "utility.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>
#include <string>

#include <sstream>

using namespace std;

namespace Cubiquity
{
	using namespace Internals;

	class NodeCounter
	{
	public:
		void operator()(NodeDAG& /*nodes*/, uint32 nodeIndex, Box3i /*bounds*/)
		{
			mUniqueNodes.insert(nodeIndex);
		}
		uint64_t count() { return mUniqueNodes.size(); }
	private:
		std::set<uint32> mUniqueNodes;
	};

	// FIXME - I think it probably makes sense to remove the usage of external material here, as it simplifies some code. 
	// From the perspective of positioning the camera we can probably assume the externaml material is just zero.
	// However, the voxelisation code is currently dpendant on the concept of extenal materials. This can probably be reviewed.
	class BoundsCalculator
	{
	public:
		BoundsCalculator(MaterialId externalMaterial)
		{
			mExternalMaterial = externalMaterial;
			mBounds.invalidate();
		}

		bool operator()(NodeDAG& nodes, uint32 nodeIndex, const Box3i& bounds)
		{
			if (isMaterialNode(nodeIndex))
			{
				MaterialId matId = static_cast<MaterialId>(nodeIndex);
				if (matId != mExternalMaterial)
				{
					mBounds.lower() = Cubiquity::min(mBounds.lower(), bounds.lower());
					mBounds.upper() = Cubiquity::max(mBounds.upper(), bounds.upper());
				}
			}

			return true; // Signal to process children
		}

		Box3i bounds()
		{
			return mBounds;
		}

	private:
		Box3i mBounds;
		MaterialId mExternalMaterial;
	};

	Box3i computeBounds(Cubiquity::Volume& volume, MaterialId externalMaterial)
	{
		BoundsCalculator boundsCalculator(externalMaterial);
		visitVolumeNodes(volume, boundsCalculator);
		return boundsCalculator.bounds();

	}

	std::pair<uint16_t, Box3i> estimateBounds(Volume& volume)
	{
		// Take a guess at what the outside material is by looking at the most extreme voxels
		// Note: This line could be shorter if we had a version of setVoxel which took a Vector3i, should add that.
		uint16_t mostNegativeVoxel = volume.voxel(std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest());
		uint16_t mostPositiveVoxel = volume.voxel(std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max());

		if (mostNegativeVoxel != mostPositiveVoxel)
		{
			log(WARN, "Unable to accurately determine outside voxel");
		}
		uint16_t outsideMaterialId = mostPositiveVoxel;

		log(DBG, "Outside material = ", outsideMaterialId);

		Box3i bounds = computeBounds(volume, outsideMaterialId);

		// If the bounds are the whole volume then we probably choose the wrong material to compute bounds for.
		// FIXME - Also check for Invalid here (max < min). A volume could be *all* empty or *not at all* empty.
		if (bounds == Box3i::invalid())
		{
			log(WARN, "Bounds are invalid, volume filled with outside material (i.e. it is empty).");
		}

		if (bounds == Box3i::max())
		{
			log(WARN, "Bounds are maxed out, did something go wrong?");
		}

		return std::make_pair(outsideMaterialId, bounds);
	}

	class HistogramCalculator
	{
	public:
		HistogramCalculator()
		{
		}

		// Add with overflow handling (see https://stackoverflow.com/a/33948556)
		bool add(uint64 a, uint64 b, uint64& result)
		{
			result = a + b;
			return ((a + b) < a); // True in case of overflow
		}

		// Multiply with overflow handling (see https://stackoverflow.com/a/1815371)
		bool multiply(uint64 a, uint64 b, uint64& result)
		{
			result = a * b;
			return(a != 0) && (result / a != b); // True in case of overflow
		}

		bool operator()(NodeDAG& nodes, uint32 nodeIndex, const Box3i& bounds)
		{
			if (isMaterialNode(nodeIndex))
			{
				MaterialId matId = static_cast<MaterialId>(nodeIndex);
				HistogramEntry& histEntry = mHistogram[matId];

				uint64 width  = static_cast<uint64>(std::max(bounds.width(),  INT64_C(0)));
				uint64 height = static_cast<uint64>(std::max(bounds.height(), INT64_C(0)));
				uint64 depth  = static_cast<uint64>(std::max(bounds.depth(),  INT64_C(0)));

				uint64 voxelsInNode = 0;
				histEntry.overflow |= multiply(width, height, voxelsInNode);
				histEntry.overflow |= multiply(voxelsInNode, depth, voxelsInNode);
				histEntry.overflow |= add(histEntry.count, voxelsInNode, histEntry.count);
			}
			return true; // Signal to process children
		}

		Histogram histogram()
		{
			return mHistogram;
		}

	private:
		Histogram mHistogram;
	};

	Histogram computeHistogram(Cubiquity::Volume& volume)
	{
		HistogramCalculator histogramCalculator;
		visitVolumeNodes(volume, histogramCalculator);
		return histogramCalculator.histogram();

	}

	void printHistogram(const Histogram& histogram)
	{
		for (const auto& entry : histogram)
		{
			if (entry.second.overflow)
			{
				log(INF, "Material ", static_cast<uint16_t>(entry.first), ": Too many to count! (64-bit overflow)");
			}
			else
			{
				// Static cast to avoid 8-bit matId being treated as ASCII char.
				log(INF, "Material ", static_cast<uint16_t>(entry.first), ": ", entry.second.count, " voxels");
			}
		}
	}
}