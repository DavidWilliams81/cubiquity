#include "morton_position_enumerator.h"

#include <cassert>

namespace Cubiquity
{
	// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
	uint32_t Compact1By2(uint32_t x)
	{
		x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		x = (x ^ (x >> 2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x >> 4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x >> 8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
		return x;
	}

	uint32_t DecodeMorton3X(uint32_t code)
	{
		return Compact1By2(code >> 0);
	}

	uint32_t DecodeMorton3Y(uint32_t code)
	{
		return Compact1By2(code >> 1);
	}

	uint32_t DecodeMorton3Z(uint32_t code)
	{
		return Compact1By2(code >> 2);
	}

	MortonPositionEnumerator::MortonPositionEnumerator(const Box3i& box)
		: mOffset(box.lower())
		, mWidth(box.sideLength(0)), mHeight(box.sideLength(1)), mDepth(box.sideLength(2))
		, mXPos(0), mYPos(0), mZPos(0)
		, mMortonCode(0)
	{
		assert(mWidth == mHeight);
		assert(mWidth == mDepth);
	}

	bool MortonPositionEnumerator::next()
	{
		if (x() == (mWidth - 1) && y() == (mHeight - 1) && z() == (mDepth - 1))
		{
			return false;
		}

		mMortonCode++;

		mXPos = DecodeMorton3X(mMortonCode);
		mYPos = DecodeMorton3Y(mMortonCode);
		mZPos = DecodeMorton3Z(mMortonCode);

		return true;
	}

	uint32_t MortonPositionEnumerator::x()
	{
		return mOffset.x() + mXPos;
	}

	uint32_t MortonPositionEnumerator::y()
	{
		return mOffset.y() + mYPos;
	}

	uint32_t MortonPositionEnumerator::z()
	{
		return mOffset.z() + mZPos;
	}
}