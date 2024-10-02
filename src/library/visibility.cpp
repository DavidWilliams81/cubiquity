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
#include "visibility.h"

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
	
	const int VisibilityMask::TileSize; // Should need this line as an int should be declarable in the header... but g++ gets upset.

	// Used for near to far octree traversal described here:
	// https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
	// Note: '4' comes before '3' in this array. This is not a mistake (see link above)
	const uint nearToFar[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };

	// Counter-clockwise winding in a right-handed (OpenGL-style) coordinate system.
	std::array<Vector4i, 6> cubeIndices = {
			Vector4i{4,6,2,0}, // min x
			Vector4i{1,3,7,5}, // max x
			Vector4i{4,0,1,5}, // min y
			Vector4i{6,7,3,2}, // max y
			Vector4i{0,2,3,1}, // min z
			Vector4i{4,5,7,6}  // max z
	};

	// Was in VisibilityMask.cpp

	VisibilityMask::VisibilityMask(uint32_t width, uint32_t height)
	{
		// The size of the rendered image does not need to match the size of the visibility
		// mask, therefore I'm not convinced it is worth the extra complexity of supporting
		// dimensions which are not a multiple of the tile size.
		if (width % TileSize || height % TileSize)
		{
			log(WARN, "Visibility mask dimensions should be a multiple of tile size (", TileSize, ")");
		}

		mWidth = width;
		mHeight = height;

		mWidthInTiles = (mWidth / TileSize);
		mHeightInTiles = (mHeight / TileSize);

		mTiles = new Tile[mWidthInTiles * mHeightInTiles];

		clear();
	}

	VisibilityMask::~VisibilityMask()
	{
		delete[] mTiles;
		mTiles = nullptr;

		std::cout << "Cached tiles count = " << mCachedTiles.size() << std::endl;
	}

	void VisibilityMask::clear()
	{
		memset(mTiles, 0, mWidthInTiles * mHeightInTiles * sizeof(*mTiles));
		mCachedTiles.clear(); // FIXME - May not need to clear this?
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
		uint32_t result = Internals::murmurHash3(mTiles, mWidthInTiles * mHeightInTiles * sizeof(*mTiles), 42);
		return result;
	}

	uint32_t VisibilityMask::getFaceSize()
	{
		assert(mWidth == mHeight);
		return mWidth;
	}

	bool VisibilityMask::pointInRect(const Vector2i& c, const Vector2i& clippedLowerLeft, const Vector2i& clippedUpperRight)
	{
		return c[0] >= clippedLowerLeft[0] && c[1] >= clippedLowerLeft[1] && c[0] <= clippedUpperRight[0] && c[1] <= clippedUpperRight[1];
	}

	// Determine which side of the edge v0->v1 the point p lies. Positive result means it
	// is to the left, negative result means it is to the right, and a result of zero means
	// that the point lies on the edge or that the two endpoints are the colocated. The result
	// also represents twice the signed are of the triangle formed by the three points.
	//
	// The actual expression is an expansion of the perp dot product of v0->v1 and v0->p,
	// which is a form of vector determinant. See http://geomalgorithms.com/vector_products.html
	int det(const Vector2i& v0, const Vector2i& v1, const Vector2i& p)
	{
		// Make sure the result doesn't overflow. See the section 'Integer Overflows':
		// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
		//assert(a.x() >= -16384 && a.x() <= 16383);
		//assert(b.x() >= -16384 && b.x() <= 16383);
		//assert(c.x() >= -16384 && c.x() <= 16383);

		return (v1[0] - v0[0]) * (p[1] - v0[1]) - (v1[1] - v0[1]) * (p[0] - v0[0]);
	}

	bool VisibilityMask::pointInQuad(const Vector2i& pointToTest, const QuadVertexArray& vertices)
	{
		// See https://fgiesen.wordpress.com/2013/02/10/optimizing-the-basic-rasterizer/
		int result = 0;
		for(int i = 0; i < 4; i++)
		{
			// If any determinant is negative then the sign bit of result gets set.
			result |= det(vertices[i], vertices[(i + 1) % 4], pointToTest);
		}
		return result >= 0; // Return true if sign bit not set (no determinants were negative).
	}

	VisibilityMask::Tile& VisibilityMask::getTile(int x, int y)
	{
		if (x >= 0 && x < mWidthInTiles && y >= 0 && y < mHeightInTiles)
		{
			Tile& tile = mTiles[y * mWidthInTiles + x];
			return tile;
		}
		else
		{
			return mBorderTile;
		}
	}

	bool VisibilityMask::drawPixel(uint32_t x, uint32_t y, bool writeEnabled)
	{
		assert(mTiles);
		assert(x < mWidth && y < mHeight);

		uint32_t tileX = x / TileSize;
		uint32_t tileY = y / TileSize;

		uint32_t offsetX = x % TileSize;
		uint32_t offsetY = y % TileSize;

		uint32_t offset = offsetX + offsetY * TileSize;

		Tile mask = Tile(0x01) << offset;

		Tile& tile = getTile(tileX, tileY);

		Tile oldResult = tile & mask;

		if (writeEnabled)
		{
			tile = tile | mask;
		}

		return oldResult == 0; // If pixel was not set before then it has now been written
	}

	bool VisibilityMask::testPixel(uint32_t x, uint32_t y)
	{
		assert(mTiles);
		assert(x < mWidth && y < mHeight);

		uint32_t tileX = x / TileSize;
		uint32_t tileY = y / TileSize;

		uint32_t offsetX = x % TileSize;
		uint32_t offsetY = y % TileSize;

		uint32_t offset = offsetX + offsetY * TileSize;

		Tile mask = Tile(0x01) << offset;

		Tile& tile = getTile(tileX, tileY);

		Tile result = tile & mask;

		return result != 0;
	}

	bool VisibilityMask::getPixel(uint x, uint y, const Tile& tile)
	{
		int shift = y * TileSize + x;
		Tile bitToTest = Tile(1);
		bitToTest <<= shift;
		return tile & bitToTest;
	}

	void VisibilityMask::setPixel(uint x, uint y, Tile& tile)
	{
		int shift = y * TileSize + x;
		Tile bitToSet = Tile(1);
		bitToSet <<= shift;
		tile |= bitToSet;
	}

	void VisibilityMask::setupQuad(const QuadVertexArray& vertices, const Vector2i& lowerCorner, Vector4i& w, Vector4i& A, Vector4i& B)
	{
		// Triangle setup. Note that the elemnts in arraya A and B are shifted by one position compared to what might be expected
		// when comparing to https://fgiesen.wordpress.com/2013/02/10/optimizing-the-basic-rasterizer/. This is because that article
		// applies a 'shift' when accessing A and B (e.g. w0 += A12 rather than w0 += A01) because they are using seperate variables
		// rather than an array. But because we add addays in our implementation it is easier to shift the position of the elements 
		// here at the start when building the array.
		A[0] = vertices[1][1] - vertices[2][1], B[0] = vertices[2][0] - vertices[1][0];
		A[1] = vertices[2][1] - vertices[3][1], B[1] = vertices[3][0] - vertices[2][0];
		A[2] = vertices[3][1] - vertices[0][1], B[2] = vertices[0][0] - vertices[3][0];
		A[3] = vertices[0][1] - vertices[1][1], B[3] = vertices[1][0] - vertices[0][0];

		// Barycentric coordinates at minX/minY corner
		w[0] = det(vertices[1], vertices[2], lowerCorner);
		w[1] = det(vertices[2], vertices[3], lowerCorner);
		w[2] = det(vertices[3], vertices[0], lowerCorner);
		w[3] = det(vertices[0], vertices[1], lowerCorner);
	}

	VisibilityMask::Tile VisibilityMask::rasteriseTile(const Vector4i& w_tile, const Vector4i& A, const Vector4i& B, const Bounds& boundsTileSpace)
	{
		const int minX = std::max(0, boundsTileSpace.lower.x());
		const int minY = std::max(0, boundsTileSpace.lower.y());

		const int maxX = std::min(TileSize - 1, boundsTileSpace.upper.x());
		const int maxY = std::min(TileSize - 1, boundsTileSpace.upper.y());

		VisibilityMask::Tile bitToTest = 0x0000000000000001;
		bitToTest <<= minY * TileSize + minX;

		Vector4i w_row = w_tile;

		w_row += B * minY;

		Tile rasterisedTile = 0;
		for (int y = minY; y <= maxY; y++)
		{
			// Barycentric coordinates at start of row
			Vector4i w = w_row;

			w += A * minX;

			for (int x = minX; x <= maxX; x++)
			{
				// If p is on or inside all edges, render pixel.
				if ((w[0] | w[1] | w[2] | w[3]) >= 0)
				{
					rasterisedTile |= bitToTest;
				}
				// One step to the right
				w += A;
				bitToTest <<= 1;
			}
			// One row step
			w_row += B;
			bitToTest <<= (TileSize - 1) - (maxX - minX);
		}

		return rasterisedTile;
	}

	bool VisibilityMask::drawNodeRef(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled)
	{
		bool drewPixel = false;
		for (int face = 0; face < 6; face++)
		{
			if (frontFaces[face])
			{
				const Vector4i& i = cubeIndices[face];
				const QuadVertexArray quadVertices{ vertices[i[0]], vertices[i[1]], vertices[i[2]], vertices[i[3]] };
				drewPixel = drawQuadRef(quadVertices, writeEnabled) || drewPixel;
			}
		}
		return drewPixel;
	}

	bool VisibilityMask::drawQuadRef(const QuadVertexArray& vertices, bool writeEnabled)
	{
		Bounds bounds = computeBounds(vertices);

		const Vector2i lowerLeft = bounds.lower;
		const Vector2i upperRight = bounds.upper;

		const Vector2i clippedLowerLeft = max(lowerLeft, Vector2i({ 0, 0 }));
		const Vector2i clippedUpperRight = min(upperRight, Vector2i({ int(mWidth) - 1, int(mHeight) - 1 }));

		Vector2i c;
		bool drewPixel = false;
		for (c[1] = clippedLowerLeft.y(); c[1] <= clippedUpperRight.y(); c[1]++)
		{
			for (c[0] = clippedLowerLeft.x(); c[0] <= clippedUpperRight.x(); c[0]++)
			{
				// Note: We have calls to both testPixel() and drawPixel(), which means the logic to find the tile within the image and
				// the pixel within the tile actually gets executed twice. I tried pulling this outside but it actually made performance
				// worse. I don't know why, note that drawPixel() is only sometimes called but that doesn't really explain it?
				if (testPixel(c[0], c[1]) == false) // Useful to avoid point-in-polygon tests for pixels which are already drawn.
				{
					if (pointInQuad(c, vertices))
					{
						drewPixel = drawPixel(c[0], c[1], writeEnabled) || drewPixel;
					}
				}
			}
		}

		return drewPixel;
	}

	bool VisibilityMask::drawNode(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled)
	{
		if (!writeEnabled)
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
		}

		//return drawNodeRef(vertices, frontFaces, writeEnabled);
		Bounds nodeBounds = computeBounds(vertices);
		int widthMinusOne = (nodeBounds.upper.x() - nodeBounds.lower.x());
		int heightMinusOne = (nodeBounds.upper.y() - nodeBounds.lower.y());
		
		if (widthMinusOne < TileSize && heightMinusOne < TileSize)
		{
			return drawNodeCached(vertices, frontFaces, nodeBounds, writeEnabled);
		}
		else
		{
			return drawNodeUncached(vertices, frontFaces, writeEnabled);
		}
	}

	int floordiv(int a, int b)
	{
		return (a < 0 ? a - (b - 1) : a) / b;
	}

	// Do a shift while allowing the shift amount to be negative
	VisibilityMask::Tile VisibilityMask::signedLeftShift(VisibilityMask::Tile value, int amount)
	{
		const int sizeOfTileInBits = static_cast<int>(sizeof(Tile)) * CHAR_BIT;
		assert(std::abs(amount) < sizeOfTileInBits && "Large shift is undefined behaviour");

		// Simple implementation using conditional
		amount >= 0 ? value <<= amount : value >>= -amount;

		// More complex implementation without conditions/branches. It shifts in both
		// directions with one masked out. I don't know if it's faster, but it might
		// be useful if we want to use SIMD, etc. Can probably be further refined.
		/*const uint sizeOfUIntInBits = sizeof(uint) * CHAR_BIT;
		uint uAmount = static_cast<uint>(amount);
		uint allSetIfPositive = (uAmount >> (sizeOfUIntInBits - 1)) - 1;
		uint allSetIfNegative = ~allSetIfPositive;
		// Should be branchless, otherwise see https://stackoverflow.com/a/12041701
		uAmount = std::abs(amount);
		value <<= (uAmount & allSetIfPositive);
		value >>= (uAmount & allSetIfNegative);*/

		return value;
	}

	bool VisibilityMask::blitTileRef(const Tile& tile, const Vector2i& position, bool writeEnabled)
	{
		bool drewPixel = false;
		for (uint y = 0; y < TileSize; y++)
		{
			for (uint x = 0; x < TileSize; x++)
			{
				if (getPixel(x, y, tile))
				{
					if (position.x() + x < mWidth && position.y() + y < mHeight)
					{
						if (testPixel(position.x() + x, position.y() + y) == false)
						{
							drawPixel(position.x() + x, position.y() + y, writeEnabled);
							drewPixel = true;
						}
					}
				}
			}
		}
		return drewPixel;
	}

	bool VisibilityMask::blitTile(const Tile& tile, const Vector2i& position, bool writeEnabled)
	{
		Tile drawnPixels = 0;
		Vector2i lowerLeftTilePos;
		// FIXME - If we clip the bounds then can we use normal int division instead of this floordiv() function?
		lowerLeftTilePos[0] = floordiv(position.x(), TileSize);
		lowerLeftTilePos[1] = floordiv(position.y(), TileSize);
		Vector2u offset;
		offset[0] = position[0] - (lowerLeftTilePos[0] * TileSize);
		offset[1] = position[1] - (lowerLeftTilePos[1] * TileSize);
		assert(offset[0] < TileSize && offset[1] < TileSize);

		// Later on we will position our tile using some bitshifts. Pixels on the left and right edges of
		// the tile will wrap around into the row above or below, and so we need to mask these out. We do
		// not have the same problem at the top and bottom of the tile because the bitshift always inserts
		// zeros at either end.
		Tile horzMask[2];
		horzMask[1] = 0x0101010101010101; // Tile full of zeros except for single column.
		horzMask[1] *= ((0x1 << offset[0]) - 1); // Expand the column to several columns.
		horzMask[0] = ~horzMask[1]; // The inverted set of columns.

		// Each potentially-cached tile overlaps at most four tiles in the visibility mask when correctly
		// positioned. It can be less than four (if either the x or y offset of zero) but I believe four
		// is the most common case. It might be desirable from a loop unrolling and/or SIMD point of view
		// to always process all four tiles, but currently we don't for the following reasons:
		//
		//   - If a tile visibility mask tile is not overlapped by the tile we are blitting, but we choose
		//     to process it anyway, then the magnitude of the bitshift can become larger than what is
		//     allowed by the standard (e.g.shifting a 64 - bit value more than 63 is undefined behaviour).
		//     We might be able to work around this by breaking the shift into two parts, using a mask,
		//     clamping the magnitude of the shift, or some other trick.
		//
		//   - We might have to take extra care to not write outside the bounds of the visibility image.
		//     I'm not too sure about this at the moment but it is something to consider.
		//
		//   - The size of the rendered node might be smaller than the size of the cached tile. This means
		//     that even if the cached tile overlaps four tiles in the visibility mask, we might be able
		//     to clip it and hence decrease the number of cases where we need to process the whole four
		//     tiles. I haven't tied exploiting this yet, and I think we would need to track the  size of
		//     the rendered node (or compute it from the tile via bitwise operations?).
		int maxTileX = offset.x() == 0 ? 0 : 1;
		int maxTileY = offset.y() == 0 ? 0 : 1;
		for (int tileY = 0; tileY <= maxTileY; tileY++)
		{
			for (int tileX = 0; tileX <= maxTileX; tileX++)
			{
				Vector2i tilePos;
				tilePos[0] = lowerLeftTilePos[0] + tileX;
				tilePos[1] = lowerLeftTilePos[1] + tileY;

				Tile tileCopy = tile;

				int shift = 0;
				shift += (int(offset.y()) - int(TileSize * tileY)) * int(TileSize);
				shift += int(offset.x()) - int(TileSize * tileX);

				tileCopy = signedLeftShift(tileCopy, shift);

				tileCopy &= horzMask[tileX];

				Tile& dstTile = getTile(tilePos[0], tilePos[1]);

				drawnPixels |= ((~dstTile) & tileCopy);

				if (writeEnabled)
				{
					dstTile |= tileCopy;
				}
			}
		}
		return drawnPixels != 0;
	}

	bool VisibilityMask::drawNodeCached(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, const Bounds& nodeBounds, bool writeEnabled)
	{
		bool drewPixel = false;

		PolygonVertexArray tileSpaceVertices;
		for (uint ct = 0; ct < 8; ct++)
		{
			tileSpaceVertices[ct] = vertices[ct] - nodeBounds.lower;
		}
		uint32 hash = murmurHash3(&tileSpaceVertices[0], sizeof(tileSpaceVertices[0]) * 8);

		Tile tile = 0;
		std::unordered_map<uint32, Tile>::iterator iter = mCachedTiles.find(hash);
		if (iter == mCachedTiles.end())
		{
			for (int face = 0; face < 6; face++)
			{
				if (frontFaces[face])
				{
					const Vector4i& i = cubeIndices[face];
					const QuadVertexArray quadVertices{ tileSpaceVertices[i[0]], tileSpaceVertices[i[1]], tileSpaceVertices[i[2]], tileSpaceVertices[i[3]] };
					drawQuadSmall(quadVertices, tile);
				}
			}
			mCachedTiles[hash] = tile;
		}
		else
		{
			tile = iter->second;
		}

		return blitTile(tile, nodeBounds.lower, writeEnabled);		
	}

	void VisibilityMask::drawQuadSmall(const QuadVertexArray& vertices, Tile& tile)
	{
		Bounds bounds = computeBounds(vertices);

		Vector4i A, B, w;
		Vector2i c = { 0, 0 };
		setupQuad(vertices, c, w, A, B);

		tile |= rasteriseTile(w, A, B, bounds);
	}

	bool VisibilityMask::drawNodeUncached(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled)
	{
		bool drewPixel = false;
		for (int face = 0; face < 6; face++)
		{
			if(frontFaces[face])
			{
				const Vector4i& i = cubeIndices[face];
				const QuadVertexArray quadVertices{ vertices[i[0]], vertices[i[1]], vertices[i[2]], vertices[i[3]] };
				drewPixel = drawQuadTiledNew(quadVertices, writeEnabled) || drewPixel;
			}
		}
		return drewPixel;
	}

	bool VisibilityMask::drawQuadTiledNew(const QuadVertexArray& vertices, bool writeEnabled)
	{
		bool drewPixel = false;

		const Bounds bounds = computeBounds(vertices);

		Bounds clippedBounds;
		clippedBounds.lower = max(bounds.lower, Vector2i({ 0, 0 }));
		clippedBounds.upper = min(bounds.upper, Vector2i({ int(mWidth) - 1, int(mHeight) - 1 }));

		int tileXBegin = clippedBounds.lower.x() / TileSize;
		int tileXEnd = clippedBounds.upper.x() / TileSize;
		int tileYBegin = clippedBounds.lower.y() / TileSize;
		int tileYEnd = clippedBounds.upper.y() / TileSize;

		Vector4i A, B;
		Vector2i c = { tileXBegin * TileSize, tileYBegin * TileSize };
		Vector4i w_tile_row;
		setupQuad(vertices, c, w_tile_row, A, B);

		for (int tileY = tileYBegin; tileY <= tileYEnd; tileY++)
		{
			Vector4i w_tile = w_tile_row;

			for (int tileX = tileXBegin; tileX <= tileXEnd; tileX++)
			{
				// Inversion means set bits indicate holes/gaps in mask.
				Tile& tile = mTiles[tileY * mWidthInTiles + tileX];
				Tile holes = ~tile;
				if (holes != 0)
				{
					Vector2i tilePos = { tileX * TileSize, tileY * TileSize };
					Bounds clippedBoundsTileSpace;
					clippedBoundsTileSpace.lower = clippedBounds.lower - tilePos;
					clippedBoundsTileSpace.upper = clippedBounds.upper - tilePos;

					Tile rasterisedTile = rasteriseTile(w_tile, A, B, clippedBoundsTileSpace);

					// If holes in mask line up with fragments then quad is visible.
					if (holes & rasterisedTile)
					{
						drewPixel = true;
						if (!writeEnabled)
						{
							return true;
						}
					}

					tile |= rasterisedTile;
				}
				w_tile += (A * TileSize);
			}
			w_tile_row += (B * TileSize);
		}

		return drewPixel;
	}

	Bounds computeBounds(const QuadVertexArray& vertices)
	{
		// See http://www.randygaul.net/2015/01/08/computing-aabb-trick/
		Bounds bounds;
		bounds.lower = vertices[0];
		bounds.upper = vertices[0];

		for (uint32_t ct = 1; ct < 4; ct++)
		{
			// FIXME - For some reason this version is slower (on Windows)
			// than the expanded version below. Should investigate why.
			//bounds.lower = (min)(bounds.lower, vertices[ct]);
			//bounds.upper = (max)(bounds.upper, vertices[ct]);

			bounds.lower[0] = (std::min)(bounds.lower.x(), vertices[ct][0]);
			bounds.lower[1] = (std::min)(bounds.lower.y(), vertices[ct][1]);
			bounds.upper[0] = (std::max)(bounds.upper.x(), vertices[ct][0]);
			bounds.upper[1] = (std::max)(bounds.upper.y(), vertices[ct][1]);
		}

		return bounds;
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

			bounds.lower[0] = (std::min)(bounds.lower.x(), vertices[ct].x());
			bounds.lower[1] = (std::min)(bounds.lower.y(), vertices[ct].y());
			bounds.upper[0] = (std::max)(bounds.upper.x(), vertices[ct].x());
			bounds.upper[1] = (std::max)(bounds.upper.y(), vertices[ct].y());
		}

		return bounds;
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

	// Was in VisibilityCalculator.cpp

	VisibilityCalculator::VisibilityCalculator()
	{
		mMaxFootprintSize = 0.3f;

		const uint32_t visMaskSize = 1024;
		mVisMask = new VisibilityMask(visMaskSize, visMaskSize);
		mVisMaskHalfFaceSize = static_cast<double>(visMaskSize) * 0.5;

		// Note: These cube vertices could be hard coded (they don't change)
		// but it probably takes less lines of code to generate them below.
		for (uint32 height = 0; height < 32; height++)
		{
			const double halfSize = static_cast<double>(uint32(1) << height) * 0.5;

			mCubeVertices[height][0] = Vector3d({ -halfSize, -halfSize, -halfSize });
			mCubeVertices[height][1] = Vector3d({ +halfSize, -halfSize, -halfSize });
			mCubeVertices[height][2] = Vector3d({ -halfSize, +halfSize, -halfSize });
			mCubeVertices[height][3] = Vector3d({ +halfSize, +halfSize, -halfSize });
			mCubeVertices[height][4] = Vector3d({ -halfSize, -halfSize, +halfSize });
			mCubeVertices[height][5] = Vector3d({ +halfSize, -halfSize, +halfSize });
			mCubeVertices[height][6] = Vector3d({ -halfSize, +halfSize, +halfSize });
			mCubeVertices[height][7] = Vector3d({ +halfSize, +halfSize, +halfSize });
		}
	}

	VisibilityCalculator::~VisibilityCalculator()
	{
		delete mVisMask;
		mVisMask = nullptr;
	}

	// This function finds the material to use for a node. Only leaf nodes have a valid material, so for a given node it 
	// descends the tree to find the (approx) nearest non-zero leaf node to the camera, and then takes the material from that.
	// Note: This functions requres a camera position. How might a 'generic' version work without this? Just take the centre
	// leaf? Or the first non-zero one we find for a fixed traversal order? Or look at all the leaves and find the most common
	// (could be slow)? Might need a solution to this if we ever want to do it in a view-independant way.
	uint32_t getMaterialForNode(float centreX, float centreY, float centreZ, uint32_t nodeIndex, const Volume* volume, const Vector3d& cameraPos)
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
			for (auto bt : nearToFar)
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

	// Note: We should probably make this operate on integers instead of floats.
	Vector3f computeNodeNormalRecursive(uint32 nodeIndex, const NodeStore& nodeData, int depth)
	{
		// Material nodes have no children, so we can't compute a normal for them.
		if (isMaterialNode(nodeIndex))
		{
			return { 0.0f, 0.0f, 0.0f };
		}

		Node node = nodeData[nodeIndex];

		Vector3f normal = { 0.0f, 0.0f, 0.0f };
		for (uint32_t z = 0; z < 2; z++)
		{
			for (int32_t y = 0; y < 2; y++)
			{
				for (uint32_t x = 0; x < 2; x++)
				{
					const uint8_t childId = z << 2 | y << 1 | x;
					const uint32_t childIndex = node[childId];

					Vector3f normalContribution = { 0.0f, 0.0f, 0.0f };
					if (childIndex == EmptyNodeIndex)
					{
						normalContribution = { float(x) * 2 - 1, float(y) * 2 - 1, float(z) * 2 - 1 };
					} else if (childIndex < MaterialCount)
					{
						normalContribution = { float(x) * 2 - 1, float(y) * 2 - 1, float(z) * 2 - 1 };
						normalContribution *= -1.0f;
					}
					else if (depth > 0)
					{
						// If a surface passes though a node then it might exactly cut the node in half,
						// meaning half the children are solid and half are empty. But it is much more
						// likely to pass through one or more of the children. It seems logical that these
						// children should have a greater contribution (arguably the only contribution?)
						// to the normal for the node. This is the logic behind the weighing factor below,
						// though it is still not clear what it should be set to.
						const float childWeighting = 100;
						normalContribution = computeNodeNormalRecursive(childIndex, nodeData, depth - 1);
						normalContribution *= childWeighting;
					}

					//normalContribution = normalize(normalContribution);
					normal += normalContribution;
				}
			}
		}
		float len = length(normal);
		if (len > 0.1)
		{
			normal = normalize(normal);
		}
		
		return normal;
	}

	Vector3f estimateNormalFromNeighbours(float x, float y, float z, uint32_t size, const Volume* volume)
	{
		Vector3f centre = { x, y, z };

		float fSize = size * 0.5f + 0.5f; // It should theorectcally be ok to sample just outside the node

		float fudgeFactor = 2.0f;    // In practice some extra scaling
		fSize *= fudgeFactor; // results in a smoother normal.

		Vector3f normal = {};

		for (int zOffset = -1; zOffset <= 1; zOffset += 1)
		{
			for (int yOffset = -1; yOffset <= 1; yOffset += 1)
			{
				for (int xOffset = -1; xOffset <= 1; xOffset += 1)
				{
					Vector3f offset({ static_cast<float>(xOffset), static_cast<float>(yOffset), static_cast<float>(zOffset) });
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

	uint32_t VisibilityCalculator::findVisibleOctreeNodes(const Volume* volume, CameraData* cameraData, NormalEstimation normalEstimation, bool subdivideMaterialNodes, Glyph* glyphs, uint32_t maxGlyphCount)
	{
		mNormalEstimation = normalEstimation;
		mSubdivideMaterialNodes = subdivideMaterialNodes;
		uint32_t glyphCount = 0;
		mVisMask->clear();


		for (uint32 height = 0; height < 32; height++)
		{
			Matrix4x4d viewMatrix = cameraData->viewMatrix();
			double halfSize = double(uint32(1) << height) * 0.5;
			Vector3d xAxis = { viewMatrix[0].x(), viewMatrix[0].y(), viewMatrix[0].z() };
			Vector3d yAxis = { viewMatrix[1].x(), viewMatrix[1].y(), viewMatrix[1].z() };
			Vector3d zAxis = { viewMatrix[2].x(), viewMatrix[2].y(), viewMatrix[2].z() };
			xAxis *= halfSize;
			yAxis *= halfSize;
			zAxis *= halfSize;

			mCubeVerticesViewSpace[height][0] = -xAxis - yAxis - zAxis;
			mCubeVerticesViewSpace[height][1] =  xAxis - yAxis - zAxis;
			mCubeVerticesViewSpace[height][2] = -xAxis + yAxis - zAxis;
			mCubeVerticesViewSpace[height][3] =  xAxis + yAxis - zAxis;
			mCubeVerticesViewSpace[height][4] = -xAxis - yAxis + zAxis;
			mCubeVerticesViewSpace[height][5] =  xAxis - yAxis + zAxis;
			mCubeVerticesViewSpace[height][6] = -xAxis + yAxis + zAxis;
			mCubeVerticesViewSpace[height][7] =  xAxis + yAxis + zAxis;

			/*for (int vertex = 0; vertex < 8; vertex++)
			{
				mCubeVerticesViewSpace[height][vertex] =
					mul(cameraData->viewMatrix(), mCubeVertices[height][vertex]);
			}*/
		}

		uint32_t rootHeight = logBase2(VolumeSideLength);
		Vector4d rootCentre = { 0.0, 0.0, 0.0, 1.0 };
		Vector4d rootCentreViewSpace = mul(cameraData->viewMatrix(), rootCentre);

		Vector3d rootCentre3 = { rootCentre[0], rootCentre[1], rootCentre[2]};
		Vector3d rootCentreViewSpace3 = { rootCentreViewSpace[0], rootCentreViewSpace[1], rootCentreViewSpace[2] };
		Vector3f rootNormal = { 0.0, 0.0, 0.0 };

		processNode(getRootNodeIndex(*volume), rootCentre3, rootCentreViewSpace3, rootHeight, rootNormal,
			volume, cameraData, glyphs, maxGlyphCount, glyphCount);

		return glyphCount;
	}

	void VisibilityCalculator::processNode(uint32 nodeIndex, const Vector3d& nodeCentre, const Vector3d& nodeCentreViewSpace, uint32 nodeHeight, const Vector3f& nodeNormal,
										   const Volume* volume, CameraData* cameraData, Glyph* glyphs, uint32_t maxGlyphCount, uint32_t& glyphCount)
	{
		const NodeStore& nodeData = getNodes(*volume).nodes();
		const Node& node = nodeData[nodeIndex];

		const uint32 childHeight = nodeHeight - 1;
		const double childSize = static_cast<double>(uint32(1) << childHeight);
		const double childHalfSize = childSize * 0.5;
		const double childHalfDiagonal = childSize * 1.73205080757 * 0.5;

		// Near to far octree traversal
		uint8 nearestChild = 0;
		const Vector3d& cameraPos = cameraData->position();
		if (cameraPos.x() > nodeCentre.x()) nearestChild |= 0x01;
		if (cameraPos.y() > nodeCentre.y()) nearestChild |= 0x02;
		if (cameraPos.z() > nodeCentre.z()) nearestChild |= 0x04;
		for(uint i = 0; i < 8; i++) // Iterate over the children
		{
			uint32_t childId = nearestChild ^ nearToFar[i]; // See octree traversal link above
			const uint32_t childIndex = isMaterialNode(nodeIndex) ? nodeIndex :node[childId];
			if (childIndex == 0) { continue; } // Empty child
			
			Vector3d childCentre = nodeCentre + mCubeVertices[childHeight][childId];
			Vector3d childCentreViewSpace = nodeCentreViewSpace + mCubeVerticesViewSpace[childHeight][childId];

			// Culling - note that testing against the four sides implicitly culls anything behind the
			// origin (as nothing can be in front of all planes behind the point where they intersect).
			// I don't think it's necessary to cull against the near clip plane, and it would 
			// probably cause a discontinuity as we later assume any node crossing z=0 is visible.
			bool insideFrustum = true;
			for (const Vector3d& planeNormal : cameraData->mNormalsViewSpace)
			{
				if (dot(childCentreViewSpace, planeNormal) < -childHalfDiagonal)
				{
					insideFrustum = false; break;
				}
			}
			if(!insideFrustum) { continue; }

			PolygonVertexArray corners2DInt;
			for (int ct = 0; ct < 8; ct++)
			{
				const Vector3d corner = childCentreViewSpace + mCubeVerticesViewSpace[childHeight][ct];

				// Much simplified version of applying the projection matrix, as we don't care about depth.
				const Matrix4x4d& projMatrix = cameraData->projMatrix();
				Vector2d corner2D = { corner[0] * projMatrix[0][0], corner[1] * projMatrix[1][1] };

				// I believe we can divide by 'z' instead of 'w' because we don't have a near clip
				// plane. It actually needs negating to match OpenGL but that doesn't matter for
				// occlusion testing purposes (I think this negation is where OpenGL switches from
				// RH to LH coordinate systems? See https://stackoverflow.com/a/12336360).
				const double invZ = 1.0 / -corner.z();
				corner2D *= invZ;

				// Map to window coordinates
				corner2D *= mVisMaskHalfFaceSize;
				corner2D += mVisMaskHalfFaceSize;

				// I'm not actually sure how OpenGL maps floatint-point viewport coordinates to ints, but
				// this seems reasonable (and it probably doesn't matter much for our purposes).
				// Note that the result may be outside the range [0, faceSize) as clipping is handled later.
				corners2DInt[ct] = static_cast<Vector2i>(round_to_int(corner2D));
			}

			// This footprint is not in any real units (e.g. pixels), it is just
			// a metric which increases with node size and decreses with distance.
			const double childFootprintSize = childSize / length(childCentreViewSpace);

			// Determine the set of visible faces based upon the camera position, and draw only those. Ideally this
			// would give exactly the same result as simply drawing all the faces because those facing away from the
			// camera should not generate fragments anyway. In practice, a face which is oriented only *slightly* away
			// from the camera can end up perpendicular to it after quantisation to the integer  pixel grid. Therefore
			// this back-face culling can introduce tiny differences, but I believe the result is more correct.
			//
			// Note: This approach should be extened to be hierarchal. If a given face of a node is
			// front-facing then the corresponding face of it's children must also be front-facing?
			FrontFaces frontFaces;
			frontFaces[0] = cameraData->position().x() < (childCentre.x() - childHalfSize);
			frontFaces[1] = cameraData->position().x() > (childCentre.x() + childHalfSize);
			frontFaces[2] = cameraData->position().y() < (childCentre.y() - childHalfSize);
			frontFaces[3] = cameraData->position().y() > (childCentre.y() + childHalfSize);
			frontFaces[4] = cameraData->position().z() < (childCentre.z() - childHalfSize);
			frontFaces[5] = cameraData->position().z() > (childCentre.z() + childHalfSize);

			const bool drawNotTraverse = // Attempt to draw (rather than traversing further) if:
				(childHeight == 0) || // We've reached a leaf
				(childFootprintSize <= mMaxFootprintSize) || // We are below the size threshold
				(isMaterialNode(childIndex) && !mSubdivideMaterialNodes); // We hit a material node.

			// Determine whether the node is visible, whilst also updating the visibility buffer if the node
			// is drawable. Nodes which can't be rasterised (due to straddling z=0) are assumed to be visible.
			const bool straddlesZeroPlane = childCentreViewSpace.z() >= -childHalfDiagonal;
			const bool isChildVisible = straddlesZeroPlane ||
										mVisMask->drawNode(corners2DInt, frontFaces, drawNotTraverse);

			Vector3f rawChildNormal = { 0.0f, 0.0f, 0.0f };
			Vector3f childNormal = { 0.0f, 0.0f, 0.0f };
			if (mNormalEstimation == NormalEstimation::FromChildren)
			{
				// There are some nodes for which we cannot compute normals from the chidren (e.g it could be a leaf
				// node, or it could be an internal node consisting of different materials but no empty space). In
				// this case we use the normal from the parent instead. Note that we only use the immediate parent
				// rather than propergating the normals right down the tree, as we have found this causes artifacts.
				rawChildNormal = computeNodeNormalRecursive(childIndex, nodeData, 3);
				if (length(rawChildNormal) > 0.01f)
				{
					childNormal = normalize(rawChildNormal);
				}
				else
				{
					childNormal = normalize(nodeNormal);
				}
			}

			if (isChildVisible)
			{
				if (drawNotTraverse)
				{
					if (mNormalEstimation == NormalEstimation::FromNeighbours)
					{
						// This generates high quality normals and also works for leaf 
						// nodes. There is no need for a contribution from the parent.
						childNormal = estimateNormalFromNeighbours(childCentre.x(), childCentre.y(), childCentre.z(), childSize, volume);
					}

					Glyph glyph;
					glyph.x = childCentre.x();
					glyph.y = childCentre.y();
					glyph.z = childCentre.z();
					glyph.size = childSize;

					glyph.a = childNormal.x();
					glyph.b = childNormal.y();
					glyph.c = childNormal.z();
					glyph.d = getMaterialForNode(childCentre.x(), childCentre.y(), childCentre.z(), childIndex, volume, cameraPos);

					assert(glyphCount < maxGlyphCount);
					glyphs[glyphCount] = glyph;
					glyphCount++;

					if (glyphCount == maxGlyphCount)
					{
						return;
					}
				}
				else
				{
					// Not drawable, descend further down the tree.
					processNode(childIndex, childCentre, childCentreViewSpace, childHeight, rawChildNormal,
						volume, cameraData, glyphs, maxGlyphCount, glyphCount);
				}
			}
		}
	}
}
