/*******************************************************************************
* Cubiquity - A micro-voxel engine for games and other applications.           *
*                                                                              *
* Written in 2024 by David Williams                                            *
*                                                                              *
* To the extent possible under law, the author(s) have dedicated all copyright *
* and related and neighboring rights to this software to the public domain     *
* worldwide. This software is distributed without any warranty.                *
*                                                                              *
* You should have received a copy of the CC0 Public Domain Dedication along    *
* with this software.                                                          *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.               *
*******************************************************************************/
#ifndef CUBIQUITY_APP_TYPEDEFS_H
#define CUBIQUITY_APP_TYPEDEFS_H

// Cubiquity has some basic maths types but they are for internal use only.
// Therefore this demo application pulls in a third party maths library.
#include "external/linalg.h"

// Integer typedefs
typedef int8_t  i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef unsigned int uint;

// GLSL-like vector and matrix typedefs.
using  vec3 = linalg::vec<float,  3>;
using dvec3 = linalg::vec<double, 3>;
using ivec3 = linalg::vec<int,    3>;
using uvec3 = linalg::vec<uint,   3>;
using bvec3 = linalg::vec<bool,   3>;

using  mat4 = linalg::mat<float,  4, 4>;
using dmat4 = linalg::mat<double, 4, 4>;

#endif // CUBIQUITY_APP_TYPEDEFS_H
