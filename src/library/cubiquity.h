/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#ifndef CUBIQUITY_H
#define CUBIQUITY_H

#include <cstdint>

namespace Cubiquity {
	class Volume;
}

void cubiquity_estimate_bounds(Cubiquity::Volume* volume, uint8_t* outside_material,
							   int32_t* lower_x, int32_t* lower_y, int32_t* lower_z,
							   int32_t* upper_x, int32_t* upper_y, int32_t* upper_z);

void cubiquity_compute_histogram(Cubiquity::Volume* volume, int64_t histogram[256]);

#endif // CUBIQUITY_H