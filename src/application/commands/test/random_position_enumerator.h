#ifndef CUBIQUITY_RANDOM_POSITION_ENUMERATOR_H
#define CUBIQUITY_RANDOM_POSITION_ENUMERATOR_H

#include "geometry.h"
#include "utility.h"

#include <cstdint>

namespace Cubiquity
{
	class RandomPositionEnumerator
	{
	public:
		RandomPositionEnumerator(const Box3i& box);

		bool next();

		uint32_t x();
		uint32_t y();
		uint32_t z();

		static bool test();

	private:

		void updatePositionsFromState();

		uint32_t mXBits;
		uint32_t mYBits;
		uint32_t mZBits;
		uint32_t mTotalBits;

		GaloisLFSR mGaloisLFSR;

		uint32_t mInitialState;
		uint32_t mXPos;
		uint32_t mYPos;
		uint32_t mZPos;

		static uint32_t mMasks[33];

		Vector3i mOffset;
	};
}

#endif // CUBIQUITY_RANDOM_POSITION_ENUMERATOR_H