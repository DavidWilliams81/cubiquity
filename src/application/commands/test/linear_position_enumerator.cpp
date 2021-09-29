#include "linear_position_enumerator.h"

#include <cassert>

namespace Cubiquity
{
	LinearPositionEnumerator::LinearPositionEnumerator(const Box3i& box)
		:mBounds(box), mXPos(0), mYPos(0), mZPos(0)
	{
	}

	bool LinearPositionEnumerator::next()
	{
		if (mXPos < mBounds.upper().x())
		{
			mXPos++;
		}
		else
		{
			if (mYPos < mBounds.upper().y())
			{
				mXPos = mBounds.lower().x();
				mYPos++;
			}
			else
			{
				if (mZPos < mBounds.upper().z())
				{
					mXPos = mBounds.lower().x();
					mYPos = mBounds.lower().y();
					mZPos++;
				}
				else
				{
					return false;
				}
			}
		}

		return true;
	}

	uint32_t LinearPositionEnumerator::x()
	{
		return mXPos;
	}

	uint32_t LinearPositionEnumerator::y()
	{
		return mYPos;
	}

	uint32_t LinearPositionEnumerator::z()
	{
		return mZPos;
	}
}
