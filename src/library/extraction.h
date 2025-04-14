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
#ifndef CUBIQUITY_RENDERING_H
#define CUBIQUITY_RENDERING_H

#include "geometry.h"
#include "storage.h"

#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <algorithm>
#include <numeric>

namespace Cubiquity
{
	using namespace Internals;

	enum class NormalEstimation { None, FromChildren, FromNeighbours };

	// Was in Glyph.h

	struct Glyph
	{
		vec3f position;
		float size;
		vec3f normal;
		float material;
	};

	static_assert(sizeof(Glyph) == sizeof(float) * 8); // Check tightly packed

	uint32_t extractGlyphs(Volume& volume, bool subdivideMaterialNodes, Glyph* glyphs, uint32_t maxGlyphCount);

}

#endif //CUBIQUITY_RENDERING_H
