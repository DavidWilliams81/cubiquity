#include "cubiquity.h"

#include "utility.h"

using namespace Cubiquity;

void cubiquity_estimate_bounds(Volume* volume, uint8_t* outside_material,
							   int32_t* lower_x, int32_t* lower_y, int32_t* lower_z,
							   int32_t* upper_x, int32_t* upper_y, int32_t* upper_z)
{
	auto result = estimateBounds(*volume);
	*outside_material = result.first;
	*lower_x = result.second.lower().x();
	*lower_y = result.second.lower().y();
	*lower_z = result.second.lower().z();
	*upper_x = result.second.upper().x();
	*upper_y = result.second.upper().y();
	*upper_z = result.second.upper().z();
}

// I'm not sure this shouold really be part of the public API. We only really use it for validating that the
// voxelisation hasn't changed, and we shouldn't need to do that once the voxelisation work is finished.
// Maybe it should be considered debug code which eventually shouldn't be part of Cubiquity at all?
void cubiquity_compute_histogram(Cubiquity::Volume* volume, int64_t histogram[256])
{
	Histogram h = computeHistogram(*volume);
	for (int i = 0; i < 256; i++)
	{
		HistogramEntry he = h[i];

		// If the uint64_t overflowed, or if it would cause an int64_t to overflow.
		if (he.overflow || he.count > static_cast<uint64_t>(std::numeric_limits<int>::max()))
		{
			histogram[i] = -1; // Indicate overflow
		}
		else
		{
			histogram[i] = static_cast<int64_t>(he.count);
		}
	}
}
