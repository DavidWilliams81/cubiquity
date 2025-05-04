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
		uint32_t offset = Random ? 7753    : 1;
		uint32_t scale  = Random ? 2976221 : 1;
		do 			{
			mPosition = mPosition * scale + offset;
		} while (mPosition >= (mWidth * mHeight * mDepth));

		updatePositionsFromState();

		// Used to indicate that there is nothing further to traverse.
		return mPosition != mStartPos;
	}

	uint32_t x() { return mOrigin.x + mXPos; }
	uint32_t y() { return  mOrigin.y + mYPos; }
	uint32_t z() { return  mOrigin.z + mZPos; }

private:

	// See https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
	uint32_t Compact1By2(uint32_t x)
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
			uint32_t position = mPosition;
			mXPos = position % mWidth;
			position /= mWidth;
			mYPos = position % mHeight;
			position /= mHeight;
			mZPos = position;
		}
	}

	static const uint32_t mStartPos = 0;

	// Domain
	ivec3 mOrigin;
	uint32_t mWidth;
	uint32_t mHeight;
	uint32_t mDepth;

	// Encoded position
	uint32_t mPosition;

	// Decoded components
	uint32_t mXPos;
	uint32_t mYPos;
	uint32_t mZPos;
};

typedef PositionEnumerator<Linear> LinearPositionEnumerator;
typedef PositionEnumerator<Morton> MortonPositionEnumerator;
typedef PositionEnumerator<Random> RandomPositionEnumerator;

#endif // CUBIQUITY_RANDOM_POSITION_ENUMERATOR_H