#include "random_position_enumerator.h"

#include <cassert>

namespace Cubiquity
{
	using namespace Internals;

	RandomPositionEnumerator::RandomPositionEnumerator(const Box3i& box)
		:mOffset(box.lower()),
		mXBits(logBase2(box.sideLength(0))), mYBits(logBase2(box.sideLength(1))), mZBits(logBase2(box.sideLength(2)))
		, mTotalBits(mXBits + mYBits + mZBits)
		, mGaloisLFSR(mMasks[mTotalBits])
	{
		// We don't have masks for really small sizes.
		assert(mTotalBits >= 4);

		// FIXME - We should also check side lengths are power of two?

		mInitialState = mGaloisLFSR.state();

		updatePositionsFromState();
	}

	bool RandomPositionEnumerator::next()
	{
		// If the state is zero then we cannot move on to another state.
		if (mGaloisLFSR.state() == 0)
		{
			return false;
		}

		mGaloisLFSR.next();

		// If we have returned to our initial state the we have competed a full
		// cycle, visiting all states except zero. Zero is special, as an LFSR
		// will get stuck in that state if it ever enters it. We use this
		// property and set zero as the final state we visit.
		if (mGaloisLFSR.state() == mInitialState)
		{
			mGaloisLFSR = GaloisLFSR(0, 0);
		}

		updatePositionsFromState();

		return true;
	}

	uint32_t RandomPositionEnumerator::x()
	{
		return mOffset.x() + mXPos;
	}

	uint32_t RandomPositionEnumerator::y()
	{
		return  mOffset.y() + mYPos;
	}

	uint32_t RandomPositionEnumerator::z()
	{
		return  mOffset.z() + mZPos;
	}

	void RandomPositionEnumerator::updatePositionsFromState()
	{
		uint32_t state = mGaloisLFSR.state();

		assert(state >> mTotalBits == 0);

		uint32_t xMask = ~(~0 << mXBits);
		uint32_t yMask = ~(~0 << mYBits);
		uint32_t zMask = ~(~0 << mZBits);

		mXPos = state & xMask;
		state = state >> mXBits;
		mYPos = state & yMask;
		state = state >> mYBits;
		mZPos = state & zMask;
	}

	bool RandomPositionEnumerator::test()
	{
		bool result = true;

		for (int maskIndex = 4; maskIndex <= 32; maskIndex++)
		{
			uint32_t mask = mMasks[maskIndex];
			GaloisLFSR lfsr(mask);

			uint64_t period = lfsr.computePeriod();
			uint64_t expectedPeriod = (0x1ull << maskIndex) - 1;

			if (period != expectedPeriod)
			{
				result = false;
			}

			// Uncomment this line to view progress
			//std::cout << maskIndex << " : Mask = " << mask << ", Period = " << period << " (expected " << expectedPeriod << ")" << std::endl;
		}

		return result;
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
	uint32_t RandomPositionEnumerator::mMasks[33] =
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
