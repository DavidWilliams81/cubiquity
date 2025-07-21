#include "bounds.h"

std::pair<ivec3, ivec3> find_bounds(Cubiquity::Volume& volume)
{
	u8 outside_material = 0;
	i32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&volume, &outside_material,
		&lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);
	return { {lower_x, lower_y, lower_z}, {upper_x, upper_y, upper_z} };
}