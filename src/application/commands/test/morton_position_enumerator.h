#ifndef CUBIQUITY_MORTON_POSITION_ENUMERATOR_H
#define CUBIQUITY_MORTON_POSITION_ENUMERATOR_H

#include "geometry.h"

#include <cstdint>

namespace Cubiquity
{
	class MortonPositionEnumerator
	{
	public:
		MortonPositionEnumerator(const Box3i& box);

		bool next();

		uint32_t x();
		uint32_t y();
		uint32_t z();

	private:

		uint32_t mMortonCode;

		uint32_t mWidth;
		uint32_t mHeight;
		uint32_t mDepth;

		uint32_t mXPos;
		uint32_t mYPos;
		uint32_t mZPos;

		Vector3i mOffset;
	};
}

#endif // CUBIQUITY_MORTON_POSITION_ENUMERATOR_H