#ifndef CUBIQUITY_POSITION_ENUMERATOR_H
#define CUBIQUITY_POSITION_ENUMERATOR_H

#include "base/types.h"
#include "utility.h"

#include <cstdint>

enum EnumerationMode { Linear, Morton, Random };

// Warning - This class enumerates over the whole 32-bit integer space, and
// therfore it is slow (taking several seconds to complete). It behaves in this
// way even if the specified domain is small. Therefore the iteration time may
// becomes dominant when operating over a small domain and it should not be used
// for benchmarking in this scenario. Benchmarking on a large domain should be ok.
template <int mode>
class PositionEnumerator
{
public:

	PositionEnumerator(const ivec3& lower, const ivec3& upper)
		:mOrigin(lower),
			mPosition(mStartPos),
			mWidth(upper.x - lower.x + 1),
			mHeight(upper.y - lower.y + 1),
			mDepth(upper.z - lower.z + 1)
	{
		if (mode == Morton)
		{
			// Not sure is this constraint is necessary but
			// I haven't thought about it enough nor tested.
			assert(Cubiquity::Internals::isPowerOf2(mWidth));
			assert(Cubiquity::Internals::isPowerOf2(mHeight));
			assert(Cubiquity::Internals::isPowerOf2(mDepth));
		}

		updatePositionsFromState();
	}

	bool next()
	{
		// Iterate over positions (either in pseudorandom order or linearly)
		// until one is found which is valid for the current bounds.
		u32 offset = Random ? 7753    : 1;
		u32 scale  = Random ? 2976221 : 1;
		do 			{
			mPosition = mPosition * scale + offset;
		} while (mPosition >= (mWidth * mHeight * mDepth));

		updatePositionsFromState();

		// Used to indicate that there is nothing further to traverse.
		return mPosition != mStartPos;
	}

	u32 x() { return mOrigin.x + mXPos; }
	u32 y() { return  mOrigin.y + mYPos; }
	u32 z() { return  mOrigin.z + mZPos; }

private:

	// See https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
	u32 Compact1By2(u32 x)
	{
		x &= 0x09249249;
		x = (x ^ (x >> 2)) & 0x030c30c3;
		x = (x ^ (x >> 4)) & 0x0300f00f;
		x = (x ^ (x >> 8)) & 0xff0000ff;
		x = (x ^ (x >> 16)) & 0x000003ff;
		return x;
	}

	void updatePositionsFromState()
	{
		if (mode == Morton)
		{
			mXPos = Compact1By2(mPosition >> 0);
			mYPos = Compact1By2(mPosition >> 1);
			mZPos = Compact1By2(mPosition >> 2);
		}
		else
		{
			u32 position = mPosition;
			mXPos = position % mWidth;
			position /= mWidth;
			mYPos = position % mHeight;
			position /= mHeight;
			mZPos = position;
		}
	}

	static const u32 mStartPos = 0;

	// Domain
	ivec3 mOrigin;
	u32 mWidth;
	u32 mHeight;
	u32 mDepth;

	// Encoded position
	u32 mPosition;

	// Decoded components
	u32 mXPos;
	u32 mYPos;
	u32 mZPos;
};

typedef PositionEnumerator<Linear> LinearPositionEnumerator;
typedef PositionEnumerator<Morton> MortonPositionEnumerator;
typedef PositionEnumerator<Random> RandomPositionEnumerator;

#endif // CUBIQUITY_RANDOM_POSITION_ENUMERATOR_H