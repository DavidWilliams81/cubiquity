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

#include "geometry.h"

#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <algorithm>
#include <cstdint>
#include <numeric>

namespace Cubiquity
{
	using namespace Internals;

	enum class NormalEstimation { None, FromChildren, FromNeighbours };

	struct Bounds
	{
		vec2i lower;
		vec2i upper;
	};

	typedef std::array<vec2i, 8> PolygonVertexArray;

	typedef std::array<vec2i, 4> QuadVertexArray;

	typedef std::array<bool, 6> FrontFaces;

	class VisibilityMask
	{
	public:

		typedef uint64 Tile;
		static const int TileSize = 8;
		static_assert(sizeof(Tile) * CHAR_BIT == TileSize * TileSize);

		VisibilityMask(uint32_t width, uint32_t height);
		~VisibilityMask();

		uint32_t width() { return mWidth; }
		uint32_t height() { return mHeight; }

		void clear();
		void setOpaque();

		uint32_t hash();

		bool pointInRect(const vec2i& c, const vec2i& clippedLowerLeft, const vec2i& clippedUpperRight);
		bool pointInQuad(const vec2i& pointToTest, const QuadVertexArray& vertices);

		bool getPixel(uint x, uint y, const Tile& tile);
		void setPixel(uint x, uint y, Tile& tile);

		bool drawPixel(uint32_t x, uint32_t y, bool writeEnabled);
		bool testPixel(uint32_t x, uint32_t y);

		void setupQuad(const QuadVertexArray& vertices, const vec2i& lowerCorner, vec4i& w, vec4i& A, vec4i& B);
		Tile rasteriseTile(const vec4i& w_tile, const vec4i& A, const vec4i& B, const Bounds& boundsTileSpace);

		bool blitTileRef(const Tile& tile, const vec2i& position, bool writeEnabled);
		bool blitTile(const Tile& tile, const vec2i& position, bool writeEnabled);

		bool drawNodeRef(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled);
		bool drawQuadRef(const QuadVertexArray& vertices, bool writeEnabled);

		bool drawNode(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled);
		bool drawNodeCached(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, const Bounds& nodeBounds, bool writeEnabled);
		bool drawNodeUncached(const PolygonVertexArray& vertices, const FrontFaces& frontFaces, bool writeEnabled);
		void drawQuadSmall(const QuadVertexArray& vertices, Tile& tile);
		bool drawQuadTiledNew(const QuadVertexArray& vertices, bool writeEnabled);

		uint32_t getFaceSize();

	private:

		Tile& getTile(int x, int y);
		Tile signedLeftShift(VisibilityMask::Tile value, int amount);

		uint32_t mWidth = 0;
		uint32_t mWidthInTiles = 0;
		uint32_t mHeight = 0;
		uint32_t mHeightInTiles = 0;
		Tile* mTiles = nullptr;

		// Setting the border tile to be fully occupied means that fragments
		// tested against it will be considered hidden (which makes sense as
		// they are offscreen) and fragments drawn to it will have no effect.
		Tile mBorderTile = 0xffffffffffffffff;

		std::unordered_map<uint32, Tile> mCachedTiles;
	};

	Bounds computeBounds(const QuadVertexArray& vertices);
	Bounds computeBounds(const PolygonVertexArray& vertices);

	// Was in Renderer.h

	class CameraData
	{
	public:
		CameraData()
			:mPos({ 0.0, 0.0, 0.0 })
			, mViewMat()
			, mProjMat() {}
		CameraData(const vec3d& position, const vec3d& centre, const vec3d& up, double fovy, double aspect)
		{
			mPos = position;
			mProjMat = perspective_matrix(fovy, aspect, 0.1, 100.0);
			mViewMat = lookAtRH(position, centre, up);
			mInvViewMat = inverse(mViewMat);
			mViewProjMat = mul(mProjMat, mViewMat);

			// Based on the link below, but simplified for view space as camea is at origin and looking down z.
			// http://www.lighthouse3d.com/tutorials/view-frustum-culling/geometric-approach-extracting-the-planes/
			double yScale = tan(fovy / 2);
			double xScale = yScale * aspect;

			const vec3d upViewSpace = { 0.0, 1.0, 0.0 };
			const vec3d rightViewSpace = { 1.0, 0.0, 0.0 };
			const vec3d forwardViewSpace = { 0.0, 0.0, -1.0 };

			mLeftNormalViewSpace = cross(-upViewSpace,
				normalize(forwardViewSpace - rightViewSpace * xScale));
			mRightNormalViewSpace = cross(upViewSpace,
				normalize(forwardViewSpace + rightViewSpace * xScale));
			mBottomNormalViewSpace = cross(rightViewSpace,
				normalize(forwardViewSpace - upViewSpace * yScale));
			mTopNormalViewSpace = cross(-rightViewSpace,
				normalize(forwardViewSpace + upViewSpace * yScale));

			mNormalsViewSpace[0] = mLeftNormalViewSpace;
			mNormalsViewSpace[1] = mRightNormalViewSpace;
			mNormalsViewSpace[2] = mBottomNormalViewSpace;
			mNormalsViewSpace[3] = mTopNormalViewSpace;
		}

		const Matrix4x4d& viewMatrix() const { return mViewMat; }
		const Matrix4x4d& invViewMatrix() const { return mInvViewMat; }
		const Matrix4x4d& projMatrix() const { return mProjMat; }
		const Matrix4x4d& viewProjMatrix() const { return mViewProjMat; }
		const vec3d& position() const { return mPos; }

	private:
		vec3d mPos;
		Matrix4x4d mViewMat;
		Matrix4x4d mInvViewMat;
		Matrix4x4d mProjMat;
		Matrix4x4d mViewProjMat;

	public:
		// Normals point to interior of frustum
		std::array<vec3d, 4> mNormalsViewSpace;
		vec3d mLeftNormalViewSpace;
		vec3d mRightNormalViewSpace;
		vec3d mBottomNormalViewSpace;
		vec3d mTopNormalViewSpace;
	};

	// Was in Glyph.h

	struct Glyph
	{
		Glyph()
		{
		}

		Glyph(float x, float y, float z, float size, float a, float b, float c, float d)
		{
			this->x = x;
			this->y = y;
			this->z = z;
			this->size = size;

			this->a = a;
			this->b = b;
			this->c = c;
			this->d = d;
		}

		float x;
		float y;
		float z;
		float size;

		float a;
		float b;
		float c;
		float d;
	};

	// Was in VisibiltyCalculator.h

	typedef std::pair<vec3i, vec4f> CachedVertex;

	// Alternative name ideas: Vis(ibility)Query/Test(er)/Finder
	class VisibilityCalculator
	{

	public:

		VisibilityCalculator();
		~VisibilityCalculator();

		uint32_t findVisibleOctreeNodes(const Volume* volume, CameraData* cameraData, NormalEstimation normalEstimation,
			                            bool subdivideMaterialNodes, Glyph* glyphs, uint32_t maxGlyphCount);
		void processNode(uint32 nodeIndex, const vec3d& nodeCentre, const vec3d& nodeCentreViewSpace, uint32 nodeHeight, const vec3f& nodeNormal,
						 const Volume* volume, CameraData* cameraData, Glyph* glyphs, uint32_t maxGlyphCount, uint32_t& glyphCount);

		float mMaxFootprintSize;

		VisibilityMask* mVisMask;
		double mVisMaskHalfFaceSize;

		std::array<std::array<vec3d, 8>, 32> mCubeVertices;
		std::array<std::array<vec3d, 8>, 32> mCubeVerticesViewSpace;

	private:
		NormalEstimation mNormalEstimation;

		// It's not yet clear how useful this setting is. Subdividing material nodes greatly increases the number of
		// glyphs but may be necessary for point-based rendering. Also smaller nodes may draw into the vis mask faster.
		bool mSubdivideMaterialNodes = false;
	};

	uint32_t getMaterialForNode(float centreX, float centreY, float centreZ, uint32_t nodeIndex, Volume* volume, const vec3d& cameraPos);
	vec3f estimateNormalFromChildren(Node node);
	vec3f estimateNormalFromNeighbours(float x, float y, float z, uint32_t size, Volume* volume);

	void computeBounds(const PolygonVertexArray& vertices, int32_t& min_x, int32_t& min_y, int32_t& max_x, int32_t& max_y, uint32_t width);
}

#endif //CUBIQUITY_RENDERING_H
