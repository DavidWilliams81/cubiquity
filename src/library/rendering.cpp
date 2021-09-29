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
#include "rendering.h"

#include "utility.h"
#include "storage.h"

#include <algorithm>
#include <bitset>
#include <stack>
#include <vector>

namespace Cubiquity
{
	using namespace Internals;

	static const uint32_t EmptyNodeIndex = 0;

	// Was in VisibilityMask.cpp

	VisibilityMask::VisibilityMask(uint32_t width, uint32_t height)
	{
		mWidth = width;
		mHeight = height;

		mWidthInTiles = (mWidth / 8) + 1; // FIXME - handle rounding
		mHeightInTiles = (mHeight / 8) + 1;

		mTiles = new uint64_t[mWidthInTiles * mHeightInTiles];

		clear();
	}

	VisibilityMask::~VisibilityMask()
	{
		delete[] mTiles;
		mTiles = nullptr;
	}

	void VisibilityMask::clear()
	{
		for (uint32_t y = 0; y < mHeightInTiles; y++)
		{
			for (uint32_t x = 0; x < mWidthInTiles; x++)
			{
				mTiles[x + y * mWidthInTiles] = 0;
			}
		}
	}

	void VisibilityMask::setOpaque()
	{
		for (uint32_t y = 0; y < mHeightInTiles; y++)
		{
			for (uint32_t x = 0; x < mWidthInTiles; x++)
			{
				mTiles[x + y * mWidthInTiles] = 0xffffffffffffffff;
			}
		}
	}

	uint32_t VisibilityMask::hash()
	{
		uint32_t result = Internals::murmurHash3(mTiles, mWidthInTiles * mHeightInTiles, 42);
		return result;
	}

	uint32_t VisibilityMask::getFaceSize()
	{
		assert(mWidth == mHeight);
		return mWidth;
	}

	uint64_t VisibilityMask::hashCubeCorners(const PolygonVertexArray& vertices)
	{
		uint64_t result = 0;

		for (int i = 0; i < 8; i++)
		{
			result |= vertices[i].x() & 0xf;
			result = result << 4;
			result |= vertices[i].y() & 0xf;
			result = result << 4;
		}

		return result;
	}

	// Test whether all edges in edge array are diagonal.
	bool VisibilityMask::allEdgesAreDiagonal(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges)
	{
		for (const Edge& edge : edges)
		{
			// If the two end points of an edge have matching x or y coordinates the the edge is not diagonal.
			if (vertices[edge.p0].x() == vertices[edge.p1].x() || vertices[edge.p0].y() == vertices[edge.p1].y())
			{
				return false;
			}
		}

		return true;
	}

	bool VisibilityMask::pointInRect(const Vector2i& c, const Vector2i& clippedLowerLeft, const Vector2i& clippedUpperRight)
	{
		return c.x() >= clippedLowerLeft.x() && c.y() >= clippedLowerLeft.y() && c.x() <= clippedUpperRight.x() && c.y() <= clippedUpperRight.y();
	}

	bool VisibilityMask::pointInPolygon(const Vector2i& pointToTest, const PolygonVertexArray& vertices, const PolygonEdgeArray& edges)
	{
		const Vector2i& c = pointToTest; // Rename for conciseness

										 // Check the point aginst each segment. There are some decisions to make here:
										 //      (1) Should we check against every segment, or early out when one fails? Early out seems
										 //          faster, but if we want to SIMD this code then the conditionals might cause issues.
										 //      (2) What order should we check segments in? If a pixel passes a particular segment then
										 //          it is likely that it will pass the neighbouring segment. Should jump to the opposite segment instead?
										 //      (3) On the other hand, if it *fails* a particular segment then the next pixel will probaly fail it too.
										 //          Should we always start by testing the last failed segment? Or sort the segments into the order in
										 //          which we encounter them when traversing left-to-right and top-to-bottom?
										 //      (4) I think 2 and 3 are not mutually exclusive. Start with the last segment, then test the opposite one.
										 //      (5) Are we better off trianglulating the polygon, so that there are
										 //          only 3 segments to test against? But the triangle bounding boxes would then overlap?
										 //for (uint32_t edgeIndex = 0; edgeIndex < diagonalEdgeCount; edgeIndex++)
		for (const auto edge : edges)
		{
			// This assert should not be commented out! It currently is because the PolygonEdgeArray
			// is currently initialised with two dummy elements at the end. This should be fixed.
			assert(edge.p0 != edge.p1); // Should never have added this edge.

										// Question - What should we do if a == b? This is not a valid edge, and we can't use
										// it to eliminate pixels. Ideally all pixels with pass the indice test when a == b,
										// but if this doesn't happen we might need to test for and skip that edge.
			const Vector2i& a = vertices[edge.p0];
			const Vector2i& b = vertices[edge.p1];

			bool inside = (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x()) >= 0;

			if (!inside) return false;
		}

		return true;
	}

	bool VisibilityMask::isConvex(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges)
	{
		for (auto edge : edges)
		{
			if (edge.p0 == edge.p1)
			{
				continue;
			}

			const Vector2i& a = vertices[edge.p0];
			const Vector2i& b = vertices[edge.p1];

			if (a == b)
			{
				continue;
			}

			for (auto c : vertices)
			{
				bool inside = (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x()) >= 0;

				if (!inside) return false;
			}
		}

		return true;
	}

	void VisibilityMask::drawPixel(uint32_t x, uint32_t y)
	{
		assert(mTiles);
		assert(x < mWidth && y < mHeight);

		uint32_t tileX = x / 8;
		uint32_t tileY = y / 8;

		uint32_t offsetX = x % 8;
		uint32_t offsetY = y % 8;

		uint32_t offset = offsetX + offsetY * 8;

		uint64_t mask = uint64_t(0x01) << offset;

		uint64_t& tile = mTiles[tileY * mWidthInTiles + tileX];

		tile = tile | mask;
	}

	bool VisibilityMask::testPixel(uint32_t x, uint32_t y)
	{
		assert(mTiles);
		assert(x < mWidth && y < mHeight);

		uint32_t tileX = x / 8;
		uint32_t tileY = y / 8;

		uint32_t offsetX = x % 8;
		uint32_t offsetY = y % 8;

		uint32_t offset = offsetX + offsetY * 8;

		uint64_t mask = uint64_t(0x01) << offset;

		uint64_t& tile = mTiles[tileY * mWidthInTiles + tileX];

		uint64_t result = tile & mask;

		return result != 0;
	}

	void VisibilityMask::drawConvexPolygon(const PolygonVertexArray& vertices)
	{
		PolygonEdgeArray segments;
		segmentsFromCorners(vertices, segments);

		Bounds bounds = computeBounds(vertices);

		drawConvexPolygon(vertices, segments, bounds);
	}

	void VisibilityMask::drawConvexPolygon(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds)
	{
		//const Vector2i lowerLeft = bounds.lower;
		//const Vector2i upperRight = bounds.upper;

		//Vector2i lowerLeftTile = lowerLeft / 8;
		//Vector2i upperRightTile = upperRight / 8;

		//const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		//const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		//int polyWidthInTiles = upperRightTile.x() - lowerLeftTile.x() + 1;
		//int polyHeightInTiles = upperRightTile.y() - lowerLeftTile.y() + 1;

		int polyWidth = bounds.upper.x() - bounds.lower.x() + 1;
		int polyHeight = bounds.upper.y() - bounds.lower.y() + 1;

		if (polyWidth > 12 && polyHeight > 12 && polyWidth <= 16 && polyHeight <= 16)
			//if (polyWidthInTiles <= 2 && polyHeightInTiles <= 2)
		{
			drawConvexPolygonSmall(vertices, segments, bounds);
		}
		else
		{
			drawConvexPolygonLarge(vertices, segments, bounds);
		}

		//drawConvexPolygonReference(vertices, segments, boundingIndices);
		//drawConvexPolygonTiny(vertices, segments, boundingIndices);
	}

	void VisibilityMask::drawConvexPolygonTiny(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i lowerLeftTile = clippedLowerLeft / 8;
		Vector2i upperRightTile = clippedUpperRight / 8;

		for (int tileY = lowerLeftTile.y(); tileY <= upperRightTile.y(); tileY++)
		{
			for (int tileX = lowerLeftTile.x(); tileX <= upperRightTile.x(); tileX++)
			{
				uint64_t& tile = mTiles[tileY * mWidthInTiles + tileX];
				if (tile == 0xffffffffffffffff)
				{
					continue;
				}

				const Vector2i clippedLowerLeftTileSpace = clippedLowerLeft - Vector2i(tileX * 8, tileY * 8);
				const Vector2i clippedUpperRightTileSpace = clippedUpperRight - Vector2i(tileX * 8, tileY * 8);

				const int firstOffsetX = std::max(0, clippedLowerLeftTileSpace.x());
				const int firstOffsetY = std::max(0, clippedLowerLeftTileSpace.y());

				const int lastOffsetX = std::min(7, clippedUpperRightTileSpace.x());
				const int lastOffsetY = std::min(7, clippedUpperRightTileSpace.y());

				uint64_t bitToDraw = 0x0000000000000001;
				bitToDraw <<= firstOffsetY * 8 + firstOffsetX;

				for (int offsetY = firstOffsetY; offsetY <= lastOffsetY; offsetY++)
				{
					for (int offsetX = firstOffsetX; offsetX <= lastOffsetX; offsetX++)
					{
						if ((bitToDraw & tile) == 0)
						{
							if (pointInPolygon(Vector2i(tileX * 8 + offsetX, tileY * 8 + offsetY), vertices, edges))
							{
								tile |= bitToDraw;
							}
						}

						bitToDraw <<= 1;
					}
					bitToDraw <<= 7 - lastOffsetX;
					bitToDraw <<= firstOffsetX;
				}
			}
		}
	}

	void VisibilityMask::drawConvexPolygonSmall(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		PolygonVertexArray localVertices;
		for (int i = 0; i < 8; i++)
		{
			localVertices[i] = vertices[i] - lowerLeft;
		}
		uint64_t hash = hashCubeCorners(localVertices);

		Entry<std::array<std::array<bool, 16>, 16>>& cachedPolygonEntry = mSmallPolygonCache.find(hash);

		if (cachedPolygonEntry.key != hash)
		{
			// FIXME - Should only iterate up to width and height here.
			for (int y = 0; y < 16; y++)
			{
				for (int x = 0; x < 16; x++)
				{
					cachedPolygonEntry.value[y][x] = pointInPolygon(Vector2i(x, y), localVertices, edges);
				}
			}
			cachedPolygonEntry.key = hash;
		}

		Vector2i lowerLeftTile = clippedLowerLeft / 8;
		Vector2i upperRightTile = clippedUpperRight / 8;

		for (int tileY = lowerLeftTile.y(); tileY <= upperRightTile.y(); tileY++)
		{
			for (int tileX = lowerLeftTile.x(); tileX <= upperRightTile.x(); tileX++)
			{
				uint64_t& tile = mTiles[tileY * mWidthInTiles + tileX];
				if (tile == 0xffffffffffffffff)
				{
					continue;
				}

				const Vector2i clippedLowerLeftTileSpace = clippedLowerLeft - Vector2i(tileX * 8, tileY * 8);
				const Vector2i clippedUpperRightTileSpace = clippedUpperRight - Vector2i(tileX * 8, tileY * 8);

				const int firstOffsetX = std::max(0, clippedLowerLeftTileSpace.x());
				const int firstOffsetY = std::max(0, clippedLowerLeftTileSpace.y());

				const int lastOffsetX = std::min(7, clippedUpperRightTileSpace.x());
				const int lastOffsetY = std::min(7, clippedUpperRightTileSpace.y());

				uint64_t bitToDraw = 0x0000000000000001;
				bitToDraw <<= firstOffsetY * 8 + firstOffsetX;

				for (int offsetY = firstOffsetY; offsetY <= lastOffsetY; offsetY++)
				{
					for (int offsetX = firstOffsetX; offsetX <= lastOffsetX; offsetX++)
					{
						Vector2i c(tileX * 8 + offsetX, tileY * 8 + offsetY);

						const Vector2i localPos = c - lowerLeft;

						if ((bitToDraw & tile) == 0)
						{
							if (cachedPolygonEntry.value[localPos.y()][localPos.x()])
							{
								tile |= bitToDraw;
							}

							/*if (pointInPolygon(c, vertices, edges))
							{
							tile |= bitToDraw;
							}*/
						}

						bitToDraw <<= 1;
					}
					bitToDraw <<= 7 - lastOffsetX;
					bitToDraw <<= firstOffsetX;
				}
			}
		}
	}

	void VisibilityMask::drawConvexPolygonLarge(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		drawConvexPolygonReference(vertices, edges, bounds);
	}

	bool VisibilityMask::testConvexPolygon(const PolygonVertexArray& vertices)
	{
		// Check whether any of the eight cube vertices are visible.  If so we can skip testing
		// the fragments. Should verify how common this is and whether it really helps speed.
		for (const auto& vertex : vertices)
		{
			// Can we eliminate this bounds test?
			if (vertex.x() >= 0 && vertex.x() < static_cast<int32>(mWidth) &&
				vertex.y() >= 0 && vertex.y() < static_cast<int32>(mHeight))
			{
				if (testPixel(vertex.x(), vertex.y()) == false)
				{
					return true;
				}
			}
		}

		Bounds bounds = computeBounds(vertices);
		return testConvexPolygon(vertices, bounds);
	}

	bool VisibilityMask::testConvexPolygon(const PolygonVertexArray& vertices, const Bounds& bounds)
	{
		PolygonEdgeArray segments;
		segmentsFromCorners(vertices, segments);

		return testConvexPolygon(vertices, segments, bounds);
	}

	bool VisibilityMask::testConvexPolygon(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds)
	{
		bool convex = isConvex(vertices, segments);
		assert(convex);

		if (!convex)
		{
			return false;
		}

		/*const Vector2i lowerLeft(vertices[boundingIndices.minX].x, vertices[boundingIndices.minY].y);
		const Vector2i upperRight(vertices[boundingIndices.maxX].x, vertices[boundingIndices.maxY].y);

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i lowerLeftTile = clippedLowerLeft / 8;
		Vector2i upperRightTile = clippedUpperRight / 8;*/

		//int polyWidth = bounds.upper.x - bounds.lower.x + 1;
		//int polyHeight = bounds.upper.y - bounds.lower.y + 1;

		//if (polyWidth <= 4 && polyHeight <= 4)
		{
			return testConvexPolygonLarge(vertices, segments, bounds);
		}
		/*else //if (polyWidth <= 128 && polyHeight <= 128)
		{
		//return testConvexPolygonSmall(vertices, segments, boundingIndices);
		}*/
		/*else
		{
		return testConvexPolygonSmall(vertices, segments, bounds);
		}*/

		//return testConvexPolygonReference(vertices, segments, boundingIndices);
	}

	bool VisibilityMask::testConvexPolygonTiny(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i lowerLeftTile = clippedLowerLeft / 8;
		Vector2i upperRightTile = clippedUpperRight / 8;

		for (int tileY = lowerLeftTile.y(); tileY <= upperRightTile.y(); tileY++)
		{
			for (int tileX = lowerLeftTile.x(); tileX <= upperRightTile.x(); tileX++)
			{
				const uint64_t& tile = mTiles[tileY * mWidthInTiles + tileX];
				if (tile == 0xffffffffffffffff)
				{
					continue;
				}

				uint64_t bitToTest = 0x0000000000000001;
				for (int offsetY = 0; offsetY < 8; offsetY++)
				{
					for (int offsetX = 0; offsetX < 8; offsetX++)
					{
						Vector2i c(tileX * 8 + offsetX, tileY * 8 + offsetY);

						if ((bitToTest & tile) == 0)
						{
							if (pointInPolygon(c, vertices, edges))
							{
								return true;
							}
						}

						bitToTest <<= 1;
					}
				}
			}
		}
		return false;
	}

	bool VisibilityMask::testConvexPolygonSmall(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i lowerLeftTile = clippedLowerLeft / 8;
		Vector2i upperRightTile = clippedUpperRight / 8;

		for (int tileY = lowerLeftTile.y(); tileY <= upperRightTile.y(); tileY++)
		{
			for (int tileX = lowerLeftTile.x(); tileX <= upperRightTile.x(); tileX++)
			{
				uint64_t tile = mTiles[tileY * mWidthInTiles + tileX];
				if (tile == 0xffffffffffffffff)
				{
					continue;
				}

				const Vector2i clippedLowerLeftTileSpace = clippedLowerLeft - Vector2i(tileX * 8, tileY * 8);
				const Vector2i clippedUpperRightTileSpace = clippedUpperRight - Vector2i(tileX * 8, tileY * 8);

				const int firstOffsetX = std::max(0, clippedLowerLeftTileSpace.x());
				const int firstOffsetY = std::max(0, clippedLowerLeftTileSpace.y());

				const int lastOffsetX = std::min(7, clippedUpperRightTileSpace.x());
				const int lastOffsetY = std::min(7, clippedUpperRightTileSpace.y());

				uint64_t bitToDraw = 0x0000000000000001;
				bitToDraw <<= firstOffsetY * 8 + firstOffsetX;

				for (int offsetY = firstOffsetY; offsetY <= lastOffsetY; offsetY++)
				{
					for (int offsetX = firstOffsetX; offsetX <= lastOffsetX; offsetX++)
					{
						if ((bitToDraw & tile) == 0)
						{
							if (pointInPolygon(Vector2i(tileX * 8 + offsetX, tileY * 8 + offsetY), vertices, edges))
							{
								return true;
							}
						}

						bitToDraw <<= 1;
					}
					bitToDraw <<= 7 - (lastOffsetX - firstOffsetX);
					//bitToDraw <<= firstOffsetX;
				}

				mTiles[tileY * mWidthInTiles + tileX] = tile;
			}
		}
		return false;
	}

	bool VisibilityMask::testConvexPolygonLarge(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		//const Vector2i lowerLeft = bounds.lower;
		//const Vector2i upperRight = bounds.upper;

		Bounds clippedBounds;
		clippedBounds.lower = max(bounds.lower, Vector2i(0, 0));
		clippedBounds.upper = min(bounds.upper, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i lowerLeftTile = clippedBounds.lower / 8;
		Vector2i upperRightTile = clippedBounds.upper / 8;

		Vector2i tilePos = lowerLeftTile * 8;

		for (int tileY = lowerLeftTile.y(); tileY <= upperRightTile.y(); tileY++, tilePos.data[1] += 8)
		{
			tilePos.data[0] = lowerLeftTile.x() * 8;
			for (int tileX = lowerLeftTile.x(); tileX <= upperRightTile.x(); tileX++, tilePos.data[0] += 8)
			{
				assert(tilePos.x() == tileX * 8);
				assert(tilePos.y() == tileY * 8);

				/*if(tilePos.x != tileX * 8)
				{
				std::cout << "Error" << std::endl;
				}
				else
				{
				std::cout << "Fine" << std::endl;
				}*/

				// Could take a refernce here, but it seems slightly faster to do a copy?
				// Maybe the tile data gets pulled into a register or something?
				const uint64_t tile = mTiles[tileY * mWidthInTiles + tileX];
				if (tile == 0xffffffffffffffff)
				{
					continue;
				}

				// Should perhaps apply these optimisations in the case of large polygons.
				if (false)
				{
					Vector2i p00(tileX * 8 + 0, tileY * 8 + 0);
					Vector2i p08(tileX * 8 + 0, tileY * 8 + 8);
					Vector2i p80(tileX * 8 + 8, tileY * 8 + 0);
					Vector2i p88(tileX * 8 + 0, tileY * 8 + 8);

					bool b00 = pointInPolygon(p00, vertices, edges) && pointInRect(p00, clippedBounds.lower, clippedBounds.upper);
					bool b08 = pointInPolygon(p08, vertices, edges) && pointInRect(p08, clippedBounds.lower, clippedBounds.upper);
					bool b80 = pointInPolygon(p80, vertices, edges) && pointInRect(p80, clippedBounds.lower, clippedBounds.upper);
					bool b88 = pointInPolygon(p88, vertices, edges) && pointInRect(p88, clippedBounds.lower, clippedBounds.upper);

					// Note this test doesn't seem to help performance?
					if (b00 && b08 && b80 && b88) // All corners are inside
					{
						// This tile should be filled. We already checked it has holes in it, so it is visible.
						//std::cout << "returning" << std::endl;
						return true;
					}

					if (!b00 && !b08 && !b80 && !b88)
					{
						// We could potentially early out here. If the corners of the tile are all outside the polygon then
						// the rest of the tile is outside the polygon too, *except* in the case where a corner of the polygon
						// is inside the tile. For a large polygon this will be relatively rare (there are many tiles, but
						// only a few corners of the polygon), but I still need to think how to detect these efficiently.
					}
				}

				const Vector2i clippedLowerLeftTileSpace = clippedBounds.lower - Vector2i(tileX * 8, tileY * 8);
				const Vector2i clippedUpperRightTileSpace = clippedBounds.upper - Vector2i(tileX * 8, tileY * 8);

				const int firstOffsetX = std::max(0, clippedLowerLeftTileSpace.x());
				const int firstOffsetY = std::max(0, clippedLowerLeftTileSpace.y());

				const int lastOffsetX = std::min(7, clippedUpperRightTileSpace.x());
				const int lastOffsetY = std::min(7, clippedUpperRightTileSpace.y());

				const Vector2i firstPixel(tileX * 8 + firstOffsetX, tileY * 8 + firstOffsetY);
				const Vector2i lastPixel(tileX * 8 + lastOffsetX, tileY * 8 + lastOffsetY);

				uint64_t bitToTest = 0x0000000000000001;
				bitToTest <<= firstOffsetY * 8 + firstOffsetX;

				Vector2i pixel;
				for (pixel.data[1] = firstPixel.y(); pixel.data[1] <= lastPixel.y(); pixel.data[1]++)
				{
					for (pixel.data[0] = firstPixel.x(); pixel.data[0] <= lastPixel.x(); pixel.data[0]++)
					{
						if ((bitToTest & tile) == 0)
						{
							if (pointInPolygon(pixel, vertices, edges))
							{
								return true;
							}
						}

						bitToTest <<= 1;
					}
					bitToTest <<= 7 - (lastOffsetX - firstOffsetX);
				}
			}
		}
		return false;
	}

	void VisibilityMask::drawConvexPolygonReference(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i c;
		for (c.data[1] = clippedLowerLeft.y(); c.data[1] <= clippedUpperRight.y(); c.data[1]++)
		{
			for (c.data[0] = clippedLowerLeft.x(); c.data[0] <= clippedUpperRight.x(); c.data[0]++)
			{
				// Note: We have calls to both testPixel() and drawPixel(), which means the logic to find the tile within the image and
				// the pixel within the tile actually gets executed twice. I tried pulling this outside but it actually made performance
				// worse. I don't know why, note that drawPixel() is only sometimes called but that doesn't really explain it?
				if (testPixel(c.x(), c.y()) == false) // Useful to avoid point-in-polygon tests for pixels which are already drawn.
				{
					if (pointInPolygon(c, vertices, edges))
					{
						drawPixel(c.x(), c.y());
					}
				}
			}
		}
	}

	bool VisibilityMask::testConvexPolygonReference(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges, const Bounds& bounds)
	{
		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		// Our edge list should only contain diagonal edges. Testing against
		// horizontal or vertical edges is a waste of time, because those are
		// handled by restricting drawing to the convex polygon's bounding box.
		assert(allEdgesAreDiagonal(vertices, edges));

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i(0, 0));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i(mWidth - 1, mHeight - 1));

		Vector2i c;
		for (c.data[1] = clippedLowerLeft.y(); c.data[1] <= clippedUpperRight.y(); c.data[1]++)
		{
			for (c.data[0] = clippedLowerLeft.x(); c.data[0] <= clippedUpperRight.x(); c.data[0]++)
			{
				if (pointInPolygon(c, vertices, edges))
				{
					// If a pixel has not been set then nothing has been drawn to it yet. If any pixel
					// in our cube meets this criteria then the cube is visible and should be drawn.
					if (testPixel(c.x(), c.y()) == false)
					{
						return true;
					}
				}
			}
		}
		return false;
	}

	Bounds computeBounds(const PolygonVertexArray& vertices)
	{
		// See http://www.randygaul.net/2015/01/08/computing-aabb-trick/
		Bounds bounds;
		bounds.lower = vertices[0];
		bounds.upper = vertices[0];

		for (uint32_t ct = 1; ct < 8; ct++)
		{
			// FIXME - For some reason this version is slower (on Windows)
			// than the expanded version below. Should investigate why.
			//bounds.lower = (min)(bounds.lower, vertices[ct]);
			//bounds.upper = (max)(bounds.upper, vertices[ct]);

			bounds.lower.data[0] = (std::min)(bounds.lower.x(), vertices[ct].x());
			bounds.lower.data[1] = (std::min)(bounds.lower.y(), vertices[ct].y());
			bounds.upper.data[0] = (std::max)(bounds.upper.x(), vertices[ct].x());
			bounds.upper.data[1] = (std::max)(bounds.upper.y(), vertices[ct].y());
		}

		return bounds;
	}

	// Calculates the area of the parallelogram of the three points.
	// This is actually the same as the area of the triangle defined by the three points, multiplied by 2.
	// The sign of the return value can be used to determine whether 'c' is to the left or the right of the line 'ab'.
	//
	// See here: http://www.sunshine2k.de/coding/java/PointOnLine/PointOnLine.html#step2
	// See also: https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
	int32_t perpDotProduct(const Vector2i& a, const Vector2i& b, const Vector2i& c)
	{
		// Make sure the result doesn't overflow. See the section 'Integer Overflows':
		// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
		assert(a.x() >= -16384 && a.x() <= 16383);
		assert(b.x() >= -16384 && b.x() <= 16383);
		assert(c.x() >= -16384 && c.x() <= 16383);

		// Should SIMD if we pack a and b into a four-element vector?
		return (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x());
	}

	float perpDotProduct(const Vector4f& a, const Vector4f& b, const Vector4f& c)
	{
		// Make sure the result doesn't overflow. See the section 'Integer Overflows':
		// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
		assert(a.x() >= -16384 && a.x() <= 16383);
		assert(b.x() >= -16384 && b.x() <= 16383);
		assert(c.x() >= -16384 && c.x() <= 16383);

		// Should SIMD if we pack a and b into a four-element vector?
		return (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x());
	}

	int32_t distSquared(const Vector2i& a, const Vector2i& b)
	{
		int32_t xDiff = a.x() - b.x();
		int32_t yDiff = a.y() - b.y();

		return xDiff * xDiff + yDiff * yDiff;
	}

	/*
	* This function determines the set segments (pairs of vertices) which define the footprint of the cube.
	* It does this using a relatively brute-force approach of checking all vertices against each candidate
	* segment, to determine which sements have all vetices on one side or the other. Although this is
	* brute-force, there are only a limited number of vertices (8) and candidate segments (12) to check.
	* Furthermore, it should be quite straight-forward to paralleze this with SIMD by checking multiple
	* vertices against each candidate segment at a time.
	*
	* The algorithm replaces a previous implementation using 'gift-wrapping' (AKA Javis March). This
	* is more general but suffers with various degenerate cases which require special handling
	* (http://wcipeg.com/wiki/Convex_hull). Specifically we had to elimiate duplicate input vertices which
	* caused confusion on several occasions when I later wondered why they were missing. The approach also
	* seems more difficult to parallize.
	*
	* It is worth noting that the algorithm below can occasionally fail. Given a cube with eight vertices
	* one might assume that the footprint would always be convex (and usually a hexagon), but this is not
	* always true. Specifically, when the footprint get very small (just a few pixels) the rounding errors
	* which are introduced by mapping our floating-point 3D positions to integer 2D positions can cause the
	* footprint to become concave.
	*
	* The upshot of this is that some segments do not get output, so the fragments do not get discarded and
	* extra pixels end up being drawn/tested. Of course, the drawing/testing is stil constrained to the
	* bounding box of the footprint but is none-the-less incorrect. Note that the gift- wrapping approach
	* can also fail in the presence of concave footprints but in practice worked slightly better (it found
	* the convex hull, whereas the current approach can draw outside of that).
	*
	* I have not fully evaluated the consequences, but I prefer the simplicity of the current approach and
	* the fact that it is conceptually similar to the barycentric approach which is used to draw the footprints.
	* Looking towards the future I expect to call this function quite rarely as child nodes will often use
	* the same set of segments as their parents (which will be large and hence error-free), and so it might
	* not be a problem in practice.
	*
	* It might be worth having a way to detect when a set of segments do not form a closed shape (the failure
	* mode described above) so that some kind of warning can be issued, but I haven't looked in to this. If we
	* detected a problem we could even choose to skip drawing/testing that footprint.
	*
	*/
	void segmentsFromCorners(const PolygonVertexArray& verticesIn, PolygonEdgeArray& edgesOut)
	{
		// I think this is ok, because we should only fille each cached PolygonEdgeArray once.
		assert(edgesOut.size() == 0);

		// Cube edges are stored (and hence tested) in an order designed to jump around as much as possible. The idea
		// is that if a point is outside a polygon then we want to detect this with as few edge tests as possible,
		// rather than finding that all the tests pass except the last one. Therefore I think that if a given edge
		// passes, it makes sense to test the 'most different' edge next, rather than testing its neighbour.
		const Edge cubeEdges[] = { { 0,1 },{ 5,7 },{ 6,2 },{ 2,3 },{ 4,6 },{ 1,5 },{ 1,3 },{ 6,7 },{ 0,4 },{ 0,2 },{ 4,5 },{ 7,3 } };

		for (const Edge& edge : cubeEdges)
		{
			int positiveCount = 0;
			int negativeCount = 0;

			const Vector2i& a = verticesIn[edge.p0];
			const Vector2i& b = verticesIn[edge.p1];

			// If the two vertices are (almost) the same then skip this edge, as we'll probably have
			// precision issues. Don't know what the threshold should be but we can change it if we need to.
			if (a == b)
			{
				//std::cout << "Warning: Skipping edge with coincident points" << std::endl;
				continue;
			}

			if (a.x() == b.x() || a.y() == b.y())
			{
				continue;
			}

			// Check each point against the edge to see which side it is. This can be improved - we could skip points
			// 'a' and 'b' as e know they are on the edge, and we could early out as we have one point each side of
			// the line. But I also think we could SIMD this loop, so I won't complicate it with logic for now.
			for (int cornerIndex = 0; cornerIndex < 8; cornerIndex++)
			{
				const Vector2i& c = verticesIn[cornerIndex];

				int result = (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x());

				if (result > 0) positiveCount++;
				if (result < 0) negativeCount++;
			}

			assert(!(positiveCount == 0 && negativeCount == 0));

			if (negativeCount == 0)
			{
				if (edgesOut.size() >= edgesOut.capacity())
				{
					log(ERR, "Too many edges!");
					return; // Return the ones we have found so far.
				}

				Edge edgeToAdd(edge.p0, edge.p1);
				edgesOut.push_back(edgeToAdd);
			}
			if (positiveCount == 0)
			{
				if (edgesOut.size() >= edgesOut.capacity())
				{
					log(ERR, "Too many edges!");
					return; // Return the ones we have found so far.
				}

				Edge edgeToAdd(edge.p1, edge.p0);
				edgesOut.push_back(edgeToAdd);
			}
		}
	}

	// Compute the 2D convex hull of a set of 2D input points.
	//
	// Note that the number of output indices is never greater than the number of input vertices.
	// That is, every vertex is used at most once. This might seem obvious, but an exception could
	// be a set of colinear points (a,b,c,d) where the computed hull was (a,b,c,d,c,b,a). Fortunatly
	// we handle such colinear cases and so don't have this issue.
	//
	// This implmementation uses the 'gift-wrapping' (AKA Javis March) which is described here:
	//     https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
	//     http://wcipeg.com/wiki/Convex_hull
	// Note that the second link gives some ideas on handling degenerate cases.
	//
	// The slightly odd syntax for the array parameter means the array is passed by reference. This
	// means the size is known (rather than just being a pointer) so we can use range-based for loops.
	void computeConvexHull(const PolygonVertexArray& verticesIn, PolygonEdgeArray& edgesOut)
	{
		// Begin by clearing the index array which we will later fill.
		int32_t indicesOut[8];
		uint32_t indexCount = 0; // Number of indices we have added to our array, also serves as the position the next index should be added.
		const int32_t INVALID_INDEX = -1; // Invalid indices are marked with '-1', so we clear the array to this value.
		std::fill_n(indicesOut, 8, INVALID_INDEX);

		// The gift-wrapping algorithm does not properly handle duplicate vertices, so we eliminate them here.
		// When a duplicate is found it is deleted and the next vertex takes it's place. The number of non-
		// duplicated vertices is given by uniqueVertexCount, and entries past this point in the array
		// should be ignored. Based on http://stackoverflow.com/a/1533394
		bool isDuplicate[8];
		int uniqueVertexCount = 0;
		for (uint32_t i = 1; i < 8; i++)
		{
			// Check whether 'i' has any duplicates.
			isDuplicate[i] = false;
			for (uint32_t j = i + 1; j < 8; j++)
			{
				if (verticesIn[i] == verticesIn[j])
					isDuplicate[i] = true;
				break;
			}

			if (isDuplicate[i] == false)
			{
				uniqueVertexCount++;
			}
		}

		// If we only have one unique vertex then all other vertices are duplicates of it.
		// Return no edges in this case, we want to render just a single pixel.
		if (uniqueVertexCount == 1)
		{
			return;
		}

		// We start our convex hull with the left-most point (smallest x coordinate). If we have multiple points
		// with the same smallest x coordinate then we use the 'top' one (with the smallest y coordinate). It is
		// important to choose the top left-most point and not an arbitrary left-most point as otherwise the
		// algorithm may 'jump-over' it when completing the hull and never actually terminate.
		uint32_t topLeftPointIndex = 0;
		for (uint32_t index = 0; index < 8; index++)
		{
			// Our start point should not be a duplicate as those shold be ignored.
			if (isDuplicate[index])
			{
				continue;
			}

			if (verticesIn[index].x() < verticesIn[topLeftPointIndex].x())
			{
				topLeftPointIndex = index;
			}
			else if (verticesIn[index].x() == verticesIn[topLeftPointIndex].x())
			{
				if (verticesIn[index].y() < verticesIn[topLeftPointIndex].y())
				{
					topLeftPointIndex = index;
				}
			}
		}

		// Start at the top left point, which also becomes the last (and only) point on the hull.
		const uint32_t startPointIndex = topLeftPointIndex;
		uint32_t lastKnownPointOnHullIndex = startPointIndex;

		// Add points to the convex hull until we get back to our start point.
		do
		{
			// We can't have more points in our convex hull than we have input points, so if we hit this assert
			// then there is a problem with the algorithm which has caused the loop to execute too many times.
			//assert(indexCount < uniqueVertexCount);

			// Add the last known point to our index list
			indicesOut[indexCount] = lastKnownPointOnHullIndex;
			indexCount++;

			// From the last known point on the hull we iterate over all points to find the one which is furthest outside.
			int32_t bestEndPointIndexSoFar = INVALID_INDEX; // No best point to start with
			for (uint32_t index = 0; index < 8; index++)
			{
				if (isDuplicate[index])
				{
					continue;
				}

				// Don't bother checking the point against itself
				if (index == lastKnownPointOnHullIndex)
				{
					continue;
				}

				// If we don't yet have a best point (first iteration) then this must be it.
				if (bestEndPointIndexSoFar == INVALID_INDEX)
				{
					bestEndPointIndexSoFar = index;
					continue;
				}

				// Having references to the points makes the code more readable.
				const auto& bestEndPointSoFar = verticesIn[bestEndPointIndexSoFar];
				const auto& lastKnownPointOnHull = verticesIn[lastKnownPointOnHullIndex];

				// Determine how the point to test is positioned relative to the vector from
				// the last known point on  the hull to the best (outermost) point so far.
				int32_t orientation = perpDotProduct(lastKnownPointOnHull, bestEndPointSoFar, verticesIn[index]);
				if (orientation < 0)
				{
					// If the current point is to the outside of the best point so far then it becomes the new best point.
					bestEndPointIndexSoFar = index;
				}
				else if (orientation == 0)
				{
					// If the current point is colinear with the best point then the situation is a little more complex, and
					// we decide which point to use based on it's distance. We can choose to always use the closest of any
					// given pair or the most distant. If we always choose the closet then we maximize the number of segments
					// in our hull (probably bad for our purposes) and have to be very careful not to go back on ourselves
					// (the closest point could be behind us if there are many coliear point in a row). If we choose the most
					// distant point then we minimize the number of segments in our convex hull by jumping over intermediate
					// points in a colinear group. However, this mans we can also jump over our start point (if it is inside
					// a colinear group) so the algorithm might not terminate. In our case we handle this by ensuring that
					// the starting point is not inside a colinear group by choosing the start point with the smallest 'y' if
					// there are many point with the same smallest 'x'. Hence we choose the most distant point below.
					int32_t distToCurrent = distSquared(lastKnownPointOnHull, verticesIn[index]);
					int32_t distToBestSoFar = distSquared(lastKnownPointOnHull, bestEndPointSoFar);

					if (distToCurrent > distToBestSoFar)
					{
						bestEndPointIndexSoFar = index;
					}
				}
			}

			// After looking at all the points we use the best one which we found
			lastKnownPointOnHullIndex = bestEndPointIndexSoFar;

		} while (lastKnownPointOnHullIndex != startPointIndex); // If we get back to the start point then our convex hull is complete.

		for (uint i = 0; i < indexCount - 1; i++)
		{
			Edge edgeToAdd(indicesOut[i], indicesOut[i + 1]);
			edgesOut.push_back(edgeToAdd);
		}

		if (indexCount)
		{
			Edge edgeToAdd(indicesOut[indexCount - 1], indicesOut[0]);
			edgesOut.push_back(edgeToAdd);
		}
	}

	bool isInside(const Vector3f& pos, float x, float y, float z, float size)
	{
		float halfSize = size * 0.5f;

		// less/greater or less-than/greater-than?
		return (pos.x() > (x - halfSize))
			&& (pos.x() < (x + halfSize))
			&& (pos.y() >(y - halfSize))
			&& (pos.y() < (y + halfSize))
			&& (pos.z() >(z - halfSize))
			&& (pos.z() < (z + halfSize));
	}

	bool isInsideInt(const Vector3f& pos, int32_t x, int32_t y, int32_t z, int64_t size)
	{
		Vector3i intPos(pos.x() + 0.5f, pos.y() + 0.5f, pos.z() + 0.5f);

		// less/greater or less-than/greater-than?
		return (intPos.x() >= (x))
			&& (intPos.x() <  (x + size))
			&& (intPos.y() >= (y))
			&& (intPos.y() <  (y + size))
			&& (intPos.z() >= (z))
			&& (intPos.z() <  (z + size));
	}

	// Was in Renderer.cpp

	// These are the polygon edges used in the omnidirectional mode. There are up to six edges
	// in a cube footprint, but in omnidirectional mode at least four are horizontal or vertical.
	// We don't need to store these because they correspond to the bounding box.
	const PolygonEdgeArray seg[9]
	{
		PolygonEdgeArray({ Edge(1,5),Edge(6,2) }),
		PolygonEdgeArray({ Edge(3,7),Edge(6,2) }),
		PolygonEdgeArray({ Edge(3,7),Edge(4,0) }),
		PolygonEdgeArray({ Edge(1,5),Edge(7,3) }),
		PolygonEdgeArray(),
		PolygonEdgeArray({ Edge(2,6),Edge(4,0) }),
		PolygonEdgeArray({ Edge(0,4),Edge(7,3) }),
		PolygonEdgeArray({ Edge(0,4),Edge(5,1) }),
		PolygonEdgeArray({ Edge(5,1),Edge(2,6) }),
	};

	// Data is uninitialed buy initial size is zero.
	PolygonEdgeArray segmentCache[27];

	const BoundingIndices boundingIndices[9]
	{
		BoundingIndices(0,0,7,7),
		BoundingIndices(0,0,3,7),
		BoundingIndices(4,0,3,7),
		BoundingIndices(0,0,7,3),
		BoundingIndices(0,0,3,3),
		BoundingIndices(4,0,3,3),
		BoundingIndices(0,4,7,3),
		BoundingIndices(0,4,3,3),
		BoundingIndices(4,4,3,3),
	};

	int getSectionIndex(const Vector3f& cubeCentreCameraSpace, float size)
	{
		const float halfSize = size * 0.5f;

		int x = 1;
		int y = 1;
		int z = 1;
		if (cubeCentreCameraSpace.x() < -halfSize)
		{
			x = 0;
		}
		if (cubeCentreCameraSpace.x() > halfSize)
		{
			x = 2;
		}

		if (cubeCentreCameraSpace.y() < -halfSize)
		{
			y = 0;
		}
		if (cubeCentreCameraSpace.y() > halfSize)
		{
			y = 2;
		}

		if (cubeCentreCameraSpace.z() < -halfSize)
		{
			z = 0;
		}
		if (cubeCentreCameraSpace.z() > halfSize)
		{
			z = 2;
		}

		int index = z * 9 + y * 3 + x;

		return index;
	}

	void computeBounds(const PolygonVertexArray& vertices, int32_t& min_x, int32_t& min_y, int32_t& max_x, int32_t& max_y, uint32_t /*width*/)
	{
		min_x = 10000;
		min_y = 10000;
		max_x = -10000;
		max_y = -10000;

		// Could choose to only look at the vertices which appear in the
		// index list, but the logic would be slightly more complicated.
		for (uint32_t ct = 0; ct < 8; ct++)
		{
			min_x = (std::min)(min_x, vertices[ct].x());
			max_x = (std::max)(max_x, vertices[ct].x());
			min_y = (std::min)(min_y, vertices[ct].y());
			max_y = (std::max)(max_y, vertices[ct].y());
		}
	}

	BoundingIndices computeBoundingIndices(Vector2f vertices[8])
	{
		BoundingIndices boundingIndices(0, 0, 0, 0);

		float minX = 10000;
		float minY = 10000;
		float maxX = -10000;
		float maxY = -10000;

		// Could choose to only look at the vertices which appear in the
		// index list, but the logic would be slightly more complicated.
		for (uint32_t ct = 0; ct < 8; ct++)
		{
			if (vertices[ct].x() < minX)
			{
				minX = vertices[ct].x();
				boundingIndices.minX = ct;
			}

			if (vertices[ct].y() < minY)
			{
				minY = vertices[ct].y();
				boundingIndices.minY = ct;
			}

			if (vertices[ct].x() > maxX)
			{
				maxX = vertices[ct].x();
				boundingIndices.maxX = ct;
			}

			if (vertices[ct].y() > maxY)
			{
				maxY = vertices[ct].y();
				boundingIndices.maxY = ct;
			}
		}

		return boundingIndices;
	}

	int determineView(const PolygonVertexArray& vertices)
	{
		int view = 4;

		if (vertices[0].x() > vertices[4].x()) view += 1;
		if (vertices[0].y() > vertices[4].y()) view += 3;
		if (vertices[3].x() < vertices[7].x()) view -= 1;
		if (vertices[3].y() < vertices[7].y()) view -= 3;

		return view;
	}

	// Was in VisibilityCalculator.cpp

	VisibilityCalculator::VisibilityCalculator()
	{
		mMaxLevelsToGenerate = 6;
		mMaxFootprintSize = 0.3;
		mMinDrawSize = 0.1;
		mMinTestSize = 0.1;
		mMaxNodeDistance = 200.0f;

		for (auto& face : mCubeFaces)
		{
			uint32_t faceSize = 1024;
			face = new VisibilityMask(faceSize, faceSize);
			face->clear();
		}
	}

	VisibilityCalculator::~VisibilityCalculator()
	{
		for (auto& face : mCubeFaces)
		{
			delete face;
			face = nullptr;
		}
	}


	// This function finds the material to use for a node. Only leaf nodes have a valid material, so for a given node it 
	// descends the tree to find the (approx) nearest non-zero leaf node to the camera, and then takes the material from that.
	// Note: This functions requres a camera position. How might a 'generic' version work without this? Just take the centre
	// leaf? Or the first non-zero one we find for a fixed traversal order? Or look at all the leaves and find the most common
	// (could be slow)? Might need a solution to this if we ever want to do it in a view-independant way.
	uint32_t getMaterialForNode(float centreX, float centreY, float centreZ, uint32_t nodeIndex, const Volume* volume, const Vector3f& cameraPos)
	{
		// When descending the tree I believe it would be more correct to compute the nearest child for every iteration.
		// If the camera is close to a node and near to the centre of one of it's faces then I think the nearest corner
		// of child nodes is not the same as the nearest corner of the start node. But in practice we are using this function
		// to get the material for distant nodes, so it probably doesn't matter and it seemed like it might be faster to
		// do it once at the start. We might need to come back to this in the future.
		uint8_t nearestChild = 0;
		if (cameraPos.x() > centreX) nearestChild |= 0x01;
		if (cameraPos.y() > centreY) nearestChild |= 0x02;
		if (cameraPos.z() > centreZ) nearestChild |= 0x04;

		const NodeStore& nodeData = getNodes(*volume).nodes();

		while (!isMaterialNode(nodeIndex))
		{
			// Based on Octree traversal method here: https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
			const uint8_t bitToggles[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };

			for (auto bt : bitToggles)
			{
				uint32_t childId = nearestChild ^ bt;
				uint32_t childIndex = nodeData[nodeIndex][childId];

				// Zero is a material used to denote empty space. But for the purpose of this function
				// we don't want to include it as it doesn't help provide a valid material for rendering.
				if (childIndex > 0)
				{
					nodeIndex = childIndex;
					break;
				}
			}
		}

		assert(isMaterialNode(nodeIndex));

		return nodeIndex;
	}

	Vector3f computeNodeNormal(Node node)
	{
		Vector3f normal(0.0f, 0.0f, 0.0f);
		for (uint32_t z = 0; z < 2; z++)
		{
			for (uint32_t y = 0; y < 2; y++)
			{
				for (uint32_t x = 0; x < 2; x++)
				{
					//assert(identifier.isOccupied());
					uint8_t childId = z << 2 | y << 1 | x;

					//bool hasChild = node.hasChild(childIndex);
					const uint32_t childIndex = node[childId];

					// If the child is not occupied then add a normal contribution in it's direction.
					if (childIndex == EmptyNodeIndex)
					{
						Vector3f normalContribution(int(x) * 2 - 1, int(y) * 2 - 1, int(z) * 2 - 1);
						normal += normalContribution;
					}
				}
			}
		}
		return normalize(normal);
	}

	Vector3f computeNormal(float x, float y, float z, uint32_t size, Volume* volume)
	{
		Vector3f centre(x, y, z);
		float fSize = size == 1 ? 1.0f : size * 0.5f + 0.5f;

		Vector3f normal(0.0f, 0.0f, 0.0f);

		for (int zOffset = -1; zOffset <= 1; zOffset += 2)
		{
			for (int yOffset = -1; yOffset <= 1; yOffset += 2)
			{
				for (int xOffset = -1; xOffset <= 1; xOffset += 2)
				{
					Vector3f offset(xOffset, yOffset, zOffset);
					offset *= fSize;

					Vector3f pos = centre + offset;

					bool occupied = volume->voxel(pos.x() + 0.5f, pos.y() + 0.5f, pos.z() + 0.5f);

					if (!occupied)
					{
						normal += offset;
					}
				}
			}
		}

		return normalize(normal);
	}

	Vector3f fakeNormalFromSize(int size)
	{
		unsigned int uHeight = logBase2(size);
		uHeight++; // Avoid black voxels
		unsigned int x = uHeight >> 0 & 0x01;
		unsigned int y = uHeight >> 1 & 0x01;
		unsigned int z = uHeight >> 2 & 0x01;
		return Vector3f(x, y, z) * 2.0f - 1.0f;
	}

	Glyph VisibilityCalculator::buildGlyphFromNode(float centreX, float centreY, float centreZ, uint32_t size, Node nodeForNormal, uint32_t nodeIndex, const Volume* volume, const Vector3f& cameraPos)
	{
		Glyph glyph;
		glyph.x = centreX;
		glyph.y = centreY;
		glyph.z = centreZ;
		glyph.size = size;

		//Vector3f norm = fakeNormalFromSize(size);
		Vector3f norm = computeNodeNormal(nodeForNormal);
		//Vector3f norm = computeNormal(centreX, centreY, centreZ, size, volume); // Normal from parent as we haven't updated 'node'

		/*uint32_t hash;
		Vector3f fCenter = glyphList.centreInLocalSpace();
		Vector3i centre(fCenter.x, fCenter.y, fCenter.z); // FIXME - Should round by adding 0.5 here?
		MurmurHash3_x86_32(&(centre), sizeof(centre), 42, &hash);
		norm.x = ((hash >> 0) & 0xFF);
		norm.y = ((hash >> 2) & 0xFF);
		norm.z = ((hash >> 4) & 0xFF);
		norm /= 127.5f; // Map 0.0 to 2.0
		norm -= 1.0f; // Map -1.0 to 1.0*/

		glyph.a = norm.x();
		glyph.b = norm.y();
		glyph.c = norm.z();
		glyph.d = getMaterialForNode(centreX, centreY, centreZ, nodeIndex, volume, cameraPos);

		return glyph;
	}

	uint32_t VisibilityCalculator::findVisibleOctreeNodesPerspective(CameraData* cameraData, const Volume* volume, Glyph* glyphs, uint32_t maxGlyphCount)
	{
		uint32_t glyphCount = 0;

		VisibilityMask* visMask = mCubeFaces[0];
		visMask->clear();

		//assert(rootGlyphList);
		NodeState nodeStateStack[33];

		const NodeStore& nodeData = getNodes(*volume).nodes();

		uint32_t nodeHeight = logBase2(VolumeSideLength);
		const uint32_t rootHeight = nodeHeight;

		// Initialise the root node state
		NodeState& rootNodeState = nodeStateStack[rootHeight];
		rootNodeState.mIndex = getRootNodeIndex(*volume);
		rootNodeState.mCentre = Vector3f(-0.5f);
		rootNodeState.mLowerCorner = Vector3i(std::numeric_limits<int32_t>::min());
		rootNodeState.mCentreX2 = (static_cast<Vector3i64>(rootNodeState.mLowerCorner) * int64_t(2)) + Vector3i64(VolumeSideLength - 1);

		// Note: We store our nodes in 'zyx' order, with 'x' in the LSB.
		const Vector3f& cameraPos = cameraData->position();
		const Vector3i64  cameraPosAsInt = static_cast<Vector3i64>(round(cameraPos));
		const Vector3i64  cameraPosAsIntX2 = cameraPosAsInt * int64_t(2);

		// Old floating point version below
		//if (cameraPos.x() > rootNodeState.mCentre.x()) rootNodeState.mNearestChild |= 0x01;
		//if (cameraPos.y() > rootNodeState.mCentre.y()) rootNodeState.mNearestChild |= 0x02;
		//if (cameraPos.z() > rootNodeState.mCentre.z()) rootNodeState.mNearestChild |= 0x04;

		if (cameraPosAsIntX2.x() > rootNodeState.mCentreX2.x()) rootNodeState.mNearestChild |= 0x01;
		if (cameraPosAsIntX2.y() > rootNodeState.mCentreX2.y()) rootNodeState.mNearestChild |= 0x02;
		if (cameraPosAsIntX2.z() > rootNodeState.mCentreX2.z()) rootNodeState.mNearestChild |= 0x04;

		Node node;

		// When we finish processing a node we increase 'nodeHeight' by one in order to move to its parent.
		// If the node we just finished was actually the root (i.e. we are done) then the new node height
		// will be greater than the root node height. Hence this is our termination condition.
		while (nodeHeight <= rootHeight)
		{
			// We sometimes need the parent node (e.g. when computing normals) and it is quick/easy
			// to get it here rather than putting it in the node state or finding it on demand.
			//Node parentNode = node;

			// We are at a given node height. The previous iteration will have filled
			// in the node state for the node we are supposed to process next.
			NodeState& nodeState = nodeStateStack[nodeHeight];
			node = nodeData[nodeState.mIndex];

			// It's trivial to compute the node size, so not much point storing it in the node state.
			uint64_t nodeSize = static_cast<uint64_t>(1) << nodeHeight;

			// Each node only needs to be processed once, but it gets entered multiple times (when first coming
			// in from its parent, and also when returning from processing it's children). Therefore we check if
			// this is the first time we have seen ths node (we won't have looked at it's children yet)
			//
			// NOTE - There is actually an argument for processing each (non-leaf) node more than once. I have
			// noticed that in the normal course of operation there are many more test calls than draw calls
			// (by about a factor of four) and initally this suprised me. However, if we imagine a surface slicing
			// diagonally through the volume then the visible voxels which make up this surface will often have
			// thier seven siblings hidden (lying just below the surface). In this case the parent node will pass
			// the visibility test and the first child should be drawn, but when we test the siblings they might
			// now be hidden (by the first voxel). If we were to test the parents visibility again we might
			// find it was now hidden (or perhaps we would actually need to draw the first 2-3 children before
			// this occured). At any rate, this provides a rationale for why we might want to re-test the
			// visibility of a parent after drawing a child. Alternatively, we could skip the parents visibility
			// check until the first child has been drawn.
			//
			// In practice this idea doesn't actually seem to work, presumably because we just end up testing the
			// parent more times and this does not make up for the times when we test less children (and the larger
			// parent is slower to test than the small child). But still, we should keep the idea in mind.
			if (nodeState.mProcessedChildCount == 0)
			{
				Vector3f cameraToNode = nodeState.mCentre - cameraPos;
				Vector3i64 cameraToNodeAsIntX2 = nodeState.mCentreX2 - cameraPosAsIntX2;
				//Vector3i64 cameraToNodeAsInt = cameraToNodeAsIntX2 / int64_t(2);

				const double nodeDistanceFromCameraSquared = dot(cameraToNode, cameraToNode);
				const int64_t nodeDistanceFromCameraSquaredAsIntX4 = dot(cameraToNodeAsIntX2, cameraToNodeAsIntX2);

				// 1.74 is sqrt(3) rounded up slightly. I'm not certain if or how much we really need to round up,
				// but if not then I think we need at least quite a lot of decimal places for when the node is 
				// large (2^32) and the camera is sat just inside the corner of it? Should be tested a bit more.
				double halfNodeDiagonal = nodeSize * (0.5f * 1.74f); // 1.732 = sqrt(3);
				double threshold = mMaxNodeDistance + halfNodeDiagonal;
				double thresholdSquared = threshold * threshold;
				double thresholdSquaredX4 = thresholdSquared * 4.0f;
				// Apply a limit as to how far away we render nodes.
				//if (nodeDistanceFromCameraSquared > thresholdSquared)
				if (nodeDistanceFromCameraSquaredAsIntX4 > thresholdSquaredX4)
				{
					// Go to parent
					nodeHeight += 1;
					continue;
				}

				//bool cameraInsideNode = isInside(cameraPos, nodeState.mCentre.x(), nodeState.mCentre.y(), nodeState.mCentre.z(), nodeSize);
				bool cameraInsideNode = isInsideInt(cameraPos, nodeState.mLowerCorner.x(), nodeState.mLowerCorner.y(), nodeState.mLowerCorner.z(), nodeSize);

				if (!cameraInsideNode)
				{
					bool visible = doCubeAgainstDirectional<Operations::Test>(visMask, *cameraData, nodeState.mLowerCorner, nodeSize);

					if (!visible)
					{
						// Go to parent
						nodeHeight += 1;
						continue;
					}

					// Could try a fast inv sqrt() here as we don't need much precision and actually need the inverse.
					// Also consider '_mm_rsqrt_ss', see https://www.sebastiansylvan.com/post/scalarsseintrinsics/
					// This footprint is not in any real units (e.g. pixels), it is just
					// a metric which increases with node size and decreses with distance.
					float footprintSize = static_cast<float>(nodeSize) / sqrt(nodeDistanceFromCameraSquared);

					// We stop traversing and simply render the node either if it is full (in which case we can't traverse further
					// anyway) or if the node is sufficiently small in screen-space that we have achived our desired level-of-detail.
					// FIXME - This call to isFull() was removed when restruturing the DAG code. Should we bring it back?
					if (/*node.isFull() ||*/ footprintSize < mMaxFootprintSize)
					{
						doCubeAgainstDirectional<Operations::Draw>(visMask, *cameraData, nodeState.mLowerCorner, nodeSize);

						Glyph glyph = buildGlyphFromNode(nodeState.mCentre.x(),
							nodeState.mCentre.y(), nodeState.mCentre.z(), nodeSize, node, nodeState.mIndex, volume, cameraPos);

						assert(glyphCount < maxGlyphCount);
						glyphs[glyphCount] = glyph;
						glyphCount++;

						if (glyphCount == maxGlyphCount)
						{
							return glyphCount;
						}

						// Go to parent
						nodeHeight += 1;
						continue;
					}
				}
			}

			// If we have processed all the children then get back to the parent.
			if (nodeState.mProcessedChildCount == 8)
			{
				nodeHeight += 1;
				continue;
			}

			// Based on Octree traversal method here: https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
			const uint8_t bitToggles[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };

			uint32_t childId = nodeState.mNearestChild ^ bitToggles[nodeState.mProcessedChildCount];

			nodeState.mProcessedChildCount++;

			//bool hasChild = node.hasChild(childIndex);
			const uint32_t childIndex = node[childId];

			if (childIndex == 0)
			{
				continue;
			}

			if (isMaterialNode(childIndex))
			{
				// Note that we don't bother starting a new glyph list in this case. A single voxel is unlikly to be bigger
				// than the max glyph list footprint, but even if it is we don't want a whole glyph list for a single voxel.

				// The node state is not used to 'communicate' between children, just to return to the parent's
				// state when moving back up the tree. Therefore we clear the state when entering a new child.
				nodeStateStack[nodeHeight - 1] = NodeState();

				NodeState& childNodeState = nodeStateStack[nodeHeight - 1];

				uint32_t childNodeSize = nodeSize >> 1;

				// Compute the centre of the child based on the parent and child index. Could stick the
				// '(float(childIndex >> 0 & 0x01) - 0.5f)' part in an eight-element LUT if we find that it is faster?
				childNodeState.mCentre.data[0] = nodeState.mCentre.x() + (float(childId >> 0 & 0x01) - 0.5f) * childNodeSize;
				childNodeState.mCentre.data[1] = nodeState.mCentre.y() + (float(childId >> 1 & 0x01) - 0.5f) * childNodeSize;
				childNodeState.mCentre.data[2] = nodeState.mCentre.z() + (float(childId >> 2 & 0x01) - 0.5f) * childNodeSize;

				childNodeState.mLowerCorner.data[0] = nodeState.mLowerCorner.x() + (childId >> 0 & 0x01) * childNodeSize;
				childNodeState.mLowerCorner.data[1] = nodeState.mLowerCorner.y() + (childId >> 1 & 0x01) * childNodeSize;
				childNodeState.mLowerCorner.data[2] = nodeState.mLowerCorner.z() + (childId >> 2 & 0x01) * childNodeSize;

				childNodeState.mCentreX2 = (static_cast<Vector3i64>(childNodeState.mLowerCorner) * int64_t(2)) + Vector3i64(childNodeSize - 1);

				//bool cameraInsideNode = isInside(cameraPos, childNodeState.mCentre.x(), childNodeState.mCentre.y(), childNodeState.mCentre.z(), childNodeSize);
				bool cameraInsideNode = isInsideInt(cameraPos, childNodeState.mLowerCorner.x(), childNodeState.mLowerCorner.y(), childNodeState.mLowerCorner.z(), nodeSize);

				if (!cameraInsideNode)
				{
					bool visible = doCubeAgainstDirectional<Operations::Test>(visMask, *cameraData, childNodeState.mLowerCorner, childNodeSize);

					if (visible)
					{
						doCubeAgainstDirectional<Operations::Draw>(visMask, *cameraData, childNodeState.mLowerCorner, childNodeSize);

						// Normal from parent as we haven't updated 'node'
						Glyph glyph = buildGlyphFromNode(childNodeState.mCentre.x(), childNodeState.mCentre.y(), childNodeState.mCentre.z(), childNodeSize, node, childIndex, volume, cameraPos);

						assert(glyphCount < maxGlyphCount);
						glyphs[glyphCount] = glyph;
						glyphCount++;

						if (glyphCount == maxGlyphCount)
						{
							return glyphCount;
						}
					}
				}
			}

			if (!isMaterialNode(childIndex))
			{
				// The node state is not used to 'communicate' between children, just to return to the parent's
				// state when moving back up the tree. Therefore we clear the state when entering a new child.
				nodeStateStack[nodeHeight - 1] = NodeState();

				NodeState& childNodeState = nodeStateStack[nodeHeight - 1];

				uint32_t childNodeSize = nodeSize >> 1;

				// Compute the centre of the child based on the parent and child index. Could stick the
				// '(float(childIndex >> 0 & 0x01) - 0.5f)' part in an eight-element LUT if we find that it is faster?
				childNodeState.mCentre.data[0] = nodeState.mCentre.x() + (float(childId >> 0 & 0x01) - 0.5f) * childNodeSize;
				childNodeState.mCentre.data[1] = nodeState.mCentre.y() + (float(childId >> 1 & 0x01) - 0.5f) * childNodeSize;
				childNodeState.mCentre.data[2] = nodeState.mCentre.z() + (float(childId >> 2 & 0x01) - 0.5f) * childNodeSize;

				childNodeState.mLowerCorner.data[0] = nodeState.mLowerCorner.x() + (childId >> 0 & 0x01) * childNodeSize;
				childNodeState.mLowerCorner.data[1] = nodeState.mLowerCorner.y() + (childId >> 1 & 0x01) * childNodeSize;
				childNodeState.mLowerCorner.data[2] = nodeState.mLowerCorner.z() + (childId >> 2 & 0x01) * childNodeSize;

				childNodeState.mCentreX2 = (static_cast<Vector3i64>(childNodeState.mLowerCorner) * int64_t(2)) + Vector3i64(childNodeSize - 1);

				childNodeState.mIndex = node[childId];

				// Note: We store our nodes in 'zyx' order, with 'x' in the LSB.
				// Old floating point version below.
				//if (cameraPos.x() > childNodeState.mCentre.x()) childNodeState.mNearestChild |= 0x01;
				//if (cameraPos.y() > childNodeState.mCentre.y()) childNodeState.mNearestChild |= 0x02;
				//if (cameraPos.z() > childNodeState.mCentre.z()) childNodeState.mNearestChild |= 0x04;
				if (cameraPosAsIntX2.x() > childNodeState.mCentreX2.x()) childNodeState.mNearestChild |= 0x01;
				if (cameraPosAsIntX2.y() > childNodeState.mCentreX2.y()) childNodeState.mNearestChild |= 0x02;
				if (cameraPosAsIntX2.z() > childNodeState.mCentreX2.z()) childNodeState.mNearestChild |= 0x04;

				nodeHeight -= 1;
			}
		}

		return glyphCount;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	//  ______              _                  _                                                  //
	//  | ___ \            | |                (_)                                                 //
	//  | |_/ /__ _ _   _  | |_ _ __ __ _  ___ _ _ __   __ _                                      //
	//  |    // _` | | | | | __| '__/ _` |/ __| | '_ \ / _` |                                     //
	//  | |\ \ (_| | |_| | | |_| | | (_| | (__| | | | | (_| |                                     //
	//  \_| \_\__,_|\__, |  \__|_|  \__,_|\___|_|_| |_|\__, |                                     //
	//               __/ |                              __/ |                                     //
	//              |___/                              |___/                                      //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// See "An Efficient Parametric Algorithm for Octree Traversal"
	///////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned char a;

	int first_node(double tx0, double ty0, double tz0, double txm, double tym, double tzm)
	{
		unsigned char answer = 0;   // initialize to 00000000
									// select the entry plane and set bits
		if (tx0 > ty0)
		{
			if (tx0 > tz0) // PLANE YZ
			{
				if (tym < tx0) answer |= (1 << 1);
				if (tzm < tx0) answer |= (1 << 2);
				return (int)answer;
			}
		}
		else {
			if (ty0 > tz0) // PLANE XZ
			{
				if (txm < ty0) answer |= (1 << 0);
				if (tzm < ty0) answer |= (1 << 2);
				return (int)answer;
			}
		}
		// PLANE XY
		if (txm < tz0) answer |= (1 << 0);
		if (tym < tz0) answer |= (1 << 1);
		return (int)answer;
	}

	int new_node(double txm, int x, double tym, int y, double tzm, int z)
	{
		if (txm < tym)
		{
			if (txm < tzm) { return x; }  // YZ plane
		}
		else
		{
			if (tym < tzm) { return y; } // XZ plane
		}
		return z; // XY plane;
	}

	std::string indent(uint level)
	{
		return std::string(level * 2, ' ');
	}

	void proc_subtree(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, const Internals::NodeStore& nodes, uint32 nodeIndex, RayVolumeIntersection& intersection, int level)
	{
		if (intersection) return;
		//childId = reverseBits(childId);

		//std::cout << indent(level) << childId << std::endl;

		double txm, tym, tzm;
		int currNode;

		if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
		{
			return;
		}

		if (isMaterialNode(nodeIndex))
		{
			if (nodeIndex > 0) // Occupied node
			{
				intersection.material = nodeIndex;
				intersection.distance = std::max(std::max(tx0, ty0), tz0);
				intersection.normal = Vector3d(0.0);
				if (tx0 > ty0 && tx0 > tz0)
				{
					intersection.normal[0] = -1.0;
				}
				if (ty0 > tx0 && ty0 > tz0)
				{
					intersection.normal[1] = -1.0;
				}
				if (tz0 > tx0 && tz0 > ty0)
				{
					intersection.normal[2] = -1.0;
				}

				// Flip normals if requred
				if (a & 1) { intersection.normal[0] *= -1.0f; }
				if (a & 2) { intersection.normal[1] *= -1.0f; }
				if (a & 4) { intersection.normal[2] *= -1.0f; }
			}

			return;
		}

		// FIXME - We need to handle infinite values here. Either as described in the paper, or
		// by adding a tiny offset to input vector components to make sure they are never zero.
		txm = 0.5*(tx0 + tx1);
		tym = 0.5*(ty0 + ty1);
		tzm = 0.5*(tz0 + tz1);

		currNode = first_node(tx0, ty0, tz0, txm, tym, tzm);

		do
		{
			// Note: The constants below are in the reverse order compared to the paper. Cubiquity uses the
			// LSBs in 'zyx' to index child nodes, but the paper *appears* to use 'xyz'? The paper actually
			// seems inconsistent, because Figure 1 imples 'xyz' order but Table 1 implies 'zyx'? Or they
			// are just numbering their bits differently? Maybe I am missunderstanding something.

			// FIXME: I think the calls to 'new_node' can probably be inlined and then simplified? It might
			// not be necessary to do quite so many comparisons in each case? E.g. for case '5' we could say
			// currNode = (tym < tx1 && tym < tz1) ? 7 : 8. Is that actually better? It would allow early-out
			// if the first comparison fails? Or invert the logic so that the statement will usually pass
			// (for predictive branching)?
			switch (currNode)
			{
			case 0:
				proc_subtree(tx0, ty0, tz0, txm, tym, tzm, nodes, nodes[nodeIndex][a], intersection, level + 1);
				currNode = new_node(txm, 1, tym, 2, tzm, 4);
				break;
			case 1:
				proc_subtree(txm, ty0, tz0, tx1, tym, tzm, nodes, nodes[nodeIndex][1 ^ a], intersection, level + 1);
				currNode = new_node(tx1, 8, tym, 3, tzm, 5);
				break;
			case 2:
				proc_subtree(tx0, tym, tz0, txm, ty1, tzm, nodes, nodes[nodeIndex][2 ^ a], intersection, level + 1);
				currNode = new_node(txm, 3, ty1, 8, tzm, 6);
				break;
			case 3:
				proc_subtree(txm, tym, tz0, tx1, ty1, tzm, nodes, nodes[nodeIndex][3 ^ a], intersection, level + 1);
				currNode = new_node(tx1, 8, ty1, 8, tzm, 7);
				break;
			case 4:
				proc_subtree(tx0, ty0, tzm, txm, tym, tz1, nodes, nodes[nodeIndex][4 ^ a], intersection, level + 1);
				currNode = new_node(txm, 5, tym, 6, tz1, 8);
				break;
			case 5:
				proc_subtree(txm, ty0, tzm, tx1, tym, tz1, nodes, nodes[nodeIndex][5 ^ a], intersection, level + 1);
				currNode = new_node(tx1, 8, tym, 7, tz1, 8);
				break;
			case 6:
				proc_subtree(tx0, tym, tzm, txm, ty1, tz1, nodes, nodes[nodeIndex][6 ^ a], intersection, level + 1);
				currNode = new_node(txm, 7, ty1, 8, tz1, 8);
				break;
			case 7:
				proc_subtree(txm, tym, tzm, tx1, ty1, tz1, nodes, nodes[nodeIndex][7 ^ a], intersection, level + 1);
				currNode = 8;
				break;
			}
		} while (currNode < 8);
	}

	void proc_subtree_iter(double tx0In, double ty0In, double tz0In, double tx1In, double ty1In, double tz1In, const Internals::NodeStore& nodes, uint32 nodeIndexIn, RayVolumeIntersection& intersection, int level)
	{
		struct State
		{
			void set(double tx0In, double ty0In, double tz0In, double tx1In, double ty1In, double tz1In, uint32 nodeIndexIn)
			{
				tx0 = tx0In; ty0 = ty0In; tz0 = tz0In; tx1 = tx1In; ty1 = ty1In; tz1 = tz1In; nodeIndex = nodeIndexIn; currNode = -1;
			}

			double tx0, ty0, tz0, tx1, ty1, tz1, txm, tym, tzm;
			uint32 nodeIndex;
			int currNode;
		};

		State stack[33]; // FIXME - How big should this be?
		State* pState = &(stack[0]);
		pState->set(tx0In, ty0In, tz0In, tx1In, ty1In, tz1In, nodeIndexIn);

		do
		{
			if (pState->currNode == -1)
			{
				//std::cout << indent(level) << state.childId << std::endl;

				if (pState->tx1 < 0.0 || pState->ty1 < 0.0 || pState->tz1 < 0.0)
				{
					level--;
					pState--;
					continue;
				}

				if (isMaterialNode(pState->nodeIndex))
				{
					if (pState->nodeIndex > 0) // Occupied node
					{
						intersection.material = pState->nodeIndex;
						intersection.distance = std::max(std::max(pState->tx0, pState->ty0), pState->tz0);
						intersection.normal = Vector3d(0.0);
						if (pState->tx0 > pState->ty0 && pState->tx0 > pState->tz0)
						{
							intersection.normal[0] = -1.0;
						}
						if (pState->ty0 > pState->tx0 && pState->ty0 > pState->tz0)
						{
							intersection.normal[1] = -1.0;
						}
						if (pState->tz0 > pState->tx0 && pState->tz0 > pState->ty0)
						{
							intersection.normal[2] = -1.0;
						}

						// Flip normals if requred
						if (a & 1) { intersection.normal[0] *= -1.0f; }
						if (a & 2) { intersection.normal[1] *= -1.0f; }
						if (a & 4) { intersection.normal[2] *= -1.0f; }

						return;
					}
					level--;
					pState--;
					continue;
				}

				// FIXME - We need to handle infinite values here. Either as described in the paper, or
				// by adding a tiny offset to input vector components to make sure they are never zero.
				pState->txm = 0.5*(pState->tx0 + pState->tx1);
				pState->tym = 0.5*(pState->ty0 + pState->ty1);
				pState->tzm = 0.5*(pState->tz0 + pState->tz1);

				pState->currNode = first_node(pState->tx0, pState->ty0, pState->tz0, pState->txm, pState->tym, pState->tzm);
			}

			// Note: The constants below are in the reverse order compared to the paper. Cubiquity uses the
			// LSBs in 'zyx' to index child nodes, but the paper *appears* to use 'xyz'? The paper actually
			// seems inconsistent, because Figure 1 imples 'xyz' order but Table 1 implies 'zyx'? Or they
			// are just numbering their bits differently? Maybe I am missunderstanding something.
			State* pNextState = pState + 1;
			switch (pState->currNode)
			{
			case 0:
				pNextState->set(pState->tx0, pState->ty0, pState->tz0, pState->txm, pState->tym, pState->tzm, nodes[pState->nodeIndex][a]);
				pState->currNode = new_node(pState->txm, 1, pState->tym, 2, pState->tzm, 4);
				break;
			case 1:
				//if (nodes[pState->nodeIndex].mChildren[1 ^ a] == 0) { level--; pState--; continue; }
				pNextState->set(pState->txm, pState->ty0, pState->tz0, pState->tx1, pState->tym, pState->tzm, nodes[pState->nodeIndex][1 ^ a]);
				pState->currNode = new_node(pState->tx1, 8, pState->tym, 3, pState->tzm, 5);
				break;
			case 2:
				pNextState->set(pState->tx0, pState->tym, pState->tz0, pState->txm, pState->ty1, pState->tzm, nodes[pState->nodeIndex][2 ^ a]);
				pState->currNode = new_node(pState->txm, 3, pState->ty1, 8, pState->tzm, 6);
				break;
			case 3:
				pNextState->set(pState->txm, pState->tym, pState->tz0, pState->tx1, pState->ty1, pState->tzm, nodes[pState->nodeIndex][3 ^ a]);
				pState->currNode = new_node(pState->tx1, 8, pState->ty1, 8, pState->tzm, 7);
				break;
			case 4:
				pNextState->set(pState->tx0, pState->ty0, pState->tzm, pState->txm, pState->tym, pState->tz1, nodes[pState->nodeIndex][4 ^ a]);
				pState->currNode = new_node(pState->txm, 5, pState->tym, 6, pState->tz1, 8);
				break;
			case 5:
				pNextState->set(pState->txm, pState->ty0, pState->tzm, pState->tx1, pState->tym, pState->tz1, nodes[pState->nodeIndex][5 ^ a]);
				pState->currNode = new_node(pState->tx1, 8, pState->tym, 7, pState->tz1, 8);
				break;
			case 6:
				pNextState->set(pState->tx0, pState->tym, pState->tzm, pState->txm, pState->ty1, pState->tz1, nodes[pState->nodeIndex][6 ^ a]);
				pState->currNode = new_node(pState->txm, 7, pState->ty1, 8, pState->tz1, 8);
				break;
			case 7:
				pNextState->set(pState->txm, pState->tym, pState->tzm, pState->tx1, pState->ty1, pState->tz1, nodes[pState->nodeIndex][7 ^ a]);
				pState->currNode = 8;
				break;
			case 8:
				level--;
				pState--;
				continue;
			}

			pState++;
			level++;

		} while (level >= 0);
	}


	RayVolumeIntersection ray_parameter(const Volume& volume, Ray3d ray)
	{
		// This algorithm is implmented with double precision. I have experimented with float
		// precision but found it is not sufficient. This is not surprising, considering we
		// (potentially) need sub-voxel precision covering the whole 2^32 space.
		//
		// It might be interesing to implement it using fixed point arithmetic at some point.
		// During the traversal I think the main operations are taking the average of two values
		// and some comparisons, which should all be easy enough with fixed point data encoded
		// in int64s.
		//
		// The current impementation always starts by intersecting with the root node and
		// traversing down fom there. I think there are two ways we can improve this:
		//
		//    - As large parts of the volume are often unoccupied it might make more sense to start
		//      from the 'effective root', i.e. the first node with more than one child. We don't
		//      curently have direct access to that (though we could quickly traverse to find it),
		//      but it could be considered in the future. This might also help the algorithm to
		//      work at floating point precision.
		//    - We could start traversal from the lowest node which encapsulates the ray (if we
		//      give the ray an end point or maximum length). This involves finding the bounds of 
		//      the ray and adjusting them to the ppropriate power-of-two. Details are still to
		//      be worked out.
		//
		// I think the two approaches are orthogonal and both useful. The first is probably quicker
		// to find (and can be cached for the volume?) while the second can be applied even the
		// effective root is the real root but the ray itself is much smaller.

		Ray3d inputRay = ray;

		RayVolumeIntersection intersection;
		intersection.material = 0;

		a = 0;

		const Internals::NodeStore& nodes = Internals::getNodes(volume).nodes();

		Vector3i rootLowerBound(std::numeric_limits<int32>::min());
		Vector3i rootUpperBound(std::numeric_limits<int32>::max());

		if (ray.mDir.x() < 0.0)
		{
			ray.mOrigin.data[0] += 0.5f;
			ray.mOrigin.data[0] = /*(std::numeric_limits<uint32>::max() + 1.0)*/ -ray.mOrigin.data[0];
			ray.mOrigin.data[0] -= 0.5f;
			ray.mDir.data[0] = -ray.mDir.data[0];
			a |= 1;
		}
		if (ray.mDir.y() < 0.0)
		{
			ray.mOrigin.data[1] += 0.5f;
			ray.mOrigin.data[1] = /*(std::numeric_limits<uint32>::max() + 1.0)*/ -ray.mOrigin.data[1];
			ray.mOrigin.data[1] -= 0.5f;
			ray.mDir.data[1] = -ray.mDir.data[1];
			a |= 2;
		}
		if (ray.mDir.z() < 0.0)
		{
			ray.mOrigin.data[2] += 0.5f;
			ray.mOrigin.data[2] = /*(std::numeric_limits<uint32>::max() + 1.0)*/ -ray.mOrigin.data[2];
			ray.mOrigin.data[2] -= 0.5f;
			ray.mDir.data[2] = -ray.mDir.data[2];
			a |= 4;
		}

		// FIXME - Do we need the +/-0.5 here? And the +1.0 above?
		double tx0 = ((rootLowerBound.x() - 0.5) - ray.mOrigin.x()) / ray.mDir.x();
		double tx1 = ((rootUpperBound.x() + 0.5) - ray.mOrigin.x()) / ray.mDir.x();
		double ty0 = ((rootLowerBound.y() - 0.5) - ray.mOrigin.y()) / ray.mDir.y();
		double ty1 = ((rootUpperBound.y() + 0.5) - ray.mOrigin.y()) / ray.mDir.y();
		double tz0 = ((rootLowerBound.z() - 0.5) - ray.mOrigin.z()) / ray.mDir.z();
		double tz1 = ((rootUpperBound.z() + 0.5) - ray.mOrigin.z()) / ray.mDir.z();

		if (std::max(std::max(tx0, ty0), tz0) < std::min(std::min(tx1, ty1), tz1))
		{
			proc_subtree_iter(tx0, ty0, tz0, tx1, ty1, tz1, nodes, getRootNodeIndex(volume), intersection, 0);
		}

		intersection.position = inputRay.mOrigin + (inputRay.mDir * intersection.distance);

		return intersection;
	}
}
