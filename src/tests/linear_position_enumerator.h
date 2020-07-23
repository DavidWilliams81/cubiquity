#ifndef CUBIQUITY_LINEAR_POSITION_ENUMERATOR_H
#define CUBIQUITY_LINEAR_POSITION_ENUMERATOR_H

#include "geometry.h"

#include <cstdint>

namespace Cubiquity
{
	class LinearPositionEnumerator
	{
	public:
		LinearPositionEnumerator(const Box3i& box);

		bool next();

		uint32_t x();
		uint32_t y();
		uint32_t z();

	private:

		Box3i mBounds;

		uint32_t mXPos;
		uint32_t mYPos;
		uint32_t mZPos;
	};
}

#endif // CUBIQUITY_LINEAR_POSITION_ENUMERATOR_H