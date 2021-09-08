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

#include <iostream>
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

	class BoundsCalculator
	{
	public:
		BoundsCalculator(MaterialId externalMaterial)
		{
			mExternalMaterial = externalMaterial;
			mBounds.invalidate();
		}

		bool operator()(NodeDAG& nodes, uint32 nodeIndex, Box3i bounds)
		{
			if (isMaterialNode(nodeIndex))
			{
				MaterialId matId = static_cast<MaterialId>(nodeIndex);
				if (matId != mExternalMaterial)
				{
					mBounds.accumulate(bounds.lower());
					mBounds.accumulate(bounds.upper());
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
		traverseNodesRecursive(volume, boundsCalculator);
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
			std::cout << "Warning: Unable to accurately determine outside voxel" << std::endl;
		}
		uint16_t outsideMaterialId = mostPositiveVoxel;

		std::cout << "Outside material = " << outsideMaterialId << std::endl;

		Box3i bounds = computeBounds(volume, outsideMaterialId);

		// If the bounds are the whole volume then we probably choose the wrong material to compute bounds for.
		// FIXME - Also check for Invalid here (max < min). A volume could be *all* empty or *not at all* empty.
		if (bounds == Box3i::invalid())
		{
			std::cout << "Warning: Bounds are invalid, volume filled with outside material (i.e. it is empty)." << std::endl;
		}

		if (bounds == Box3i::max())
		{
			std::cout << "Warning: Bounds are maxed out, did something go wrong?" << std::endl;
		}

		return std::make_pair(outsideMaterialId, bounds);
	}

	// Not for general boxes - only works for node bounds which are cubic,
	// powers-of-two sizes, aligned to power-of-two boundaries, etc.
	Box3i childBounds(Box3i nodeBounds, uint childId)
	{
		uint childX = (childId >> 0) & 0x01;
		uint childY = (childId >> 1) & 0x01;
		uint childZ = (childId >> 2) & 0x01;
		Vector3i childOffset(childX, childY, childZ); // childOffset holds zeros or ones.

		// Careful ordering of operations to avoid signed integer overflow. Note that child
		// node dimensions might max-out the signed integer type but should not exceed it.
		const Vector3i childNodeDimsInCells = ((nodeBounds.upper() - Vector3i(1)) / 2) - (nodeBounds.lower() / 2);
		Vector3i childLowerBound = nodeBounds.lower() + (childNodeDimsInCells * childOffset) + childOffset;
		Vector3i childUpperBound = childLowerBound + childNodeDimsInCells;
		return Box3i(childLowerBound, childUpperBound);
	}

	Histogram computeHistogram(Volume& volume, const Box3i& bounds)
	{
		Histogram histogram;
		//std::fill(histogram.begin(), histogram.end(), 0);

		for (int32 z = bounds.lower().z(); z <= bounds.upper().z(); z++)
		{
			for (int32 y = bounds.lower().y(); y <= bounds.upper().y(); y++)
			{
				for (int32 x = bounds.lower().x(); x <= bounds.upper().x(); x++)
				{
					histogram[volume.voxel(x, y, z)]++;
				}
			}
		}

		return histogram;
	}

	void printHistogram(const Histogram& histogram)
	{
		for (const auto& entry : histogram)
		{
			std::cout << "Material " << entry.first << ": " << entry.second << " voxels" << std::endl;
		}
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
	}

	// Should match shader code
	float mix(float x, float y, float a)
	{
		return x * (1.0f - a) + y * a;
	}

	// See http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl
	Vector3f hsv2rgb(Vector3f c)
	{
		Vector4f K = Vector4f(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
		Vector3f p = abs(fract(c.xxx() + K.xyz()) * 6.0f - K.www());
		return mix(K.xxx(), clamp(p - K.xxx(), 0.0f, 1.0f), c.y()) * c.z();
	}

	Vector3f colourFromMaterialId(MaterialId matId)
	{
		// Use the golden ratio to pick well-spaced colour in HSV space.
		// http://gofiguremath.org/natures-favorite-math/the-golden-ratio/the-golden-angle/
		float golden_ratio = 1.618;
		float angle = golden_ratio / 1 + golden_ratio;

		double intpart;
		float hue = std::modf(matId * angle, &intpart);

		// Saturation and value are fixed.
		float sat = 1.0;
		float val = 1.0;

		// Convert to RGB.
		Vector3f hsv = Vector3f(hue, sat, val);
		return hsv2rgb(hsv);
	}

	GaloisLFSR::GaloisLFSR(uint32_t mask, uint32_t startState)
		: mMask(mask), mState(startState) {}

	void GaloisLFSR::next()
	{
		uint32_t lsb = mState & 1;
		mState >>= 1;
		if (lsb)
		{
			mState ^= mMask;
		}
	}

	uint32_t GaloisLFSR::state()
	{
		return mState;
	}

	// Note that the LSFR will never hit zero. The period of a
	// maximal length LSFR is therefore 2^n-1 instead of 2^n.
	uint64_t GaloisLFSR::computePeriod()
	{
		uint64_t period = 0;
		uint32_t startState = mState;

		do
		{
			period++;
			next();
		} while (state() != startState);

		return period;
	}

	ShuffledSequence::ShuffledSequence(uint32 sequenceLength)
		: mSequenceLength(sequenceLength)
		, mGaloisLFSR(maximulLengthMask(logBase2(roundUpToPowerOf2(mSequenceLength + 1))))
	{
	}

	void ShuffledSequence::next()
	{
		do
		{
			mGaloisLFSR.next();
		} while (mGaloisLFSR.state() > mSequenceLength);
	}

	uint32_t ShuffledSequence::state()
	{
		// Subtract one to include zero in outputs.
		return mGaloisLFSR.state() - 1;
	}

	uint32_t ShuffledSequence::maximulLengthMask(int sizeInBits)
	{
		if (sizeInBits < 4 || sizeInBits > 32)
		{
			throw std::out_of_range("Invalid mask selection");
		}
		return MaximulLengthMasks[sizeInBits];
	}

	// These values are taken from here:
	//
	//     https://users.ece.cmu.edu/~koopman/lfsr/index.html
	//
	// For each bitsize the page provides a list of mask values which
	// correspond to maximul-length LSFRs. I've taken the first value
	// from each list and placed it in the array below.
	// 
	// The position in the array corresponds to the size (in bits) of
	// the LFSR. There is no such thing as a zero-bit version so that
	// entry is invalid, and the page above does not contain mask for
	// 1-3 bits. Therefore only elements 4-32 in the array are valid.
	const uint32_t ShuffledSequence::MaximulLengthMasks[33] =
	{
		0x00000000u, // Element zero is invalid

		0x00000000u,0x00000000u,0x00000000u,0x00000009u,
		0x00000012u,0x00000021u,0x00000041u,0x0000008Eu,
		0x00000108u,0x00000204u,0x00000402u,0x00000829u,
		0x0000100Du,0x00002015u,0x00004001u,0x00008016u,
		0x00010004u,0x00020013u,0x00040013u,0x00080004u,
		0x00100002u,0x00200001u,0x00400010u,0x0080000Du,
		0x01000004u,0x02000023u,0x04000013u,0x08000004u,
		0x10000002u,0x20000029u,0x40000004u,0x80000057u,
	};
}