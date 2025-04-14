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
#include "extraction.h"

#include "utility.h"
#include "storage.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstring>
#include <stack>
#include <vector>

namespace Cubiquity
{
	using namespace Internals;

	static const uint32_t EmptyNodeIndex = 0;

	bool is_surrounded(const Volume& volume, const vec3i& pos)
	{
		return
			volume.voxel(pos + vec3i(-1,  0,  0)) != EmptyNodeIndex &&
			volume.voxel(pos + vec3i (1,  0,  0)) != EmptyNodeIndex &&
			volume.voxel(pos + vec3i( 0, -1,  0)) != EmptyNodeIndex &&
			volume.voxel(pos + vec3i( 0,  1,  0)) != EmptyNodeIndex &&
			volume.voxel(pos + vec3i( 0,  0, -1)) != EmptyNodeIndex &&
			volume.voxel(pos + vec3i( 0,  0,  1)) != EmptyNodeIndex;
	}

	vec3f get_normal_from_neighbours(const Volume& volume, const vec3i& pos)
	{
		vec3f normal = { 0.0f, 0.0f, 0.0f };
		for (int z = -1; z <= 1; z++)
		{
			for (int y = -1; y <= 1; y++)
			{
				for (int x = -1; x <= 1; x++)
				{
					if (volume.voxel(pos + vec3i(x, y, z)) == EmptyNodeIndex)
					{
						vec3f contribution = normalize(vec3f(x, y, z));
						normal += contribution;
					}
				}
			}
		}
		return normalize(normal);
	}

	class GlyphExtractor
	{
	public:
		GlyphExtractor(const Volume& volume, bool subdivideMaterialNodes, Glyph* glyphs, uint32_t maxGlyphCount)
			: volume(volume)
			, glyphs(glyphs), glyphCount(0), maxGlyphCount(maxGlyphCount)
			, mSubdivideMaterialNodes(subdivideMaterialNodes)
		{
		}

		bool operator()(NodeDAG& nodes, uint32 nodeIndex, const Box3i& bounds)
		{
			if ((isMaterialNode(nodeIndex)) && (nodeIndex > 0)) // Non-empty leaf node
			{
				if (mSubdivideMaterialNodes)
				{
					for (int z = bounds.lower().z; z <= bounds.upper().z; z++)
					{
						for (int y = bounds.lower().y; y <= bounds.upper().y; y++)
						{
							for (int x = bounds.lower().x; x <= bounds.upper().x; x++)
							{
								if ((x == bounds.lower().x) || (x == bounds.upper().x) ||
									(y == bounds.lower().y) || (y == bounds.upper().y) ||
									(z == bounds.lower().z) || (z == bounds.upper().z))
								{
									if (is_surrounded(volume, vec3i(x, y, z))) continue;
									Glyph glyph;
									glyph.position = vec3f(x, y, z);
									glyph.size = 1.0f;

									glyph.normal = get_normal_from_neighbours(volume, vec3i(x, y, z));
									glyph.material = nodeIndex; // Material

									// FIXME - If this condition fails we should exit loop early?
									if (glyphCount < maxGlyphCount)
									{
										glyphs[glyphCount] = glyph;
										glyphCount++;
									}
								}
							}
						}
					}
				}
				else
				{
					Box3f floatBounds = static_cast<Box3f>(bounds);
					Glyph glyph;
					glyph.position = floatBounds.centre();
					glyph.size = (floatBounds.upper().x - floatBounds.lower().x) + 1.0f; // All dimensions should be the same.;

					glyph.normal = vec3f(0.0f, 0.0f, 0.0f);
					glyph.material = nodeIndex; // Material

					// FIXME - If this condition fails we should also prevent traversal of sibling nodes?
					if (glyphCount < maxGlyphCount)
					{
						glyphs[glyphCount] = glyph;
						glyphCount++;
					}

					return false;
				}
			}

			return true;
		}

	public:
		const Volume& volume;
		Glyph* glyphs;
		uint32_t glyphCount = 0;
		uint32_t maxGlyphCount;
		bool mSubdivideMaterialNodes;
	};

	uint32_t extractGlyphs(Volume& volume, bool subdivideMaterialNodes, Glyph* glyphs, uint32_t maxGlyphCount)
	{
		GlyphExtractor glyphExtractor(volume, subdivideMaterialNodes, glyphs, maxGlyphCount);
		visitVolumeNodes(volume, glyphExtractor);
		return glyphExtractor.glyphCount;
	}
}
