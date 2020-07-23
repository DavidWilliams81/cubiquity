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
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <algorithm>
#include <cstdint>
#include <numeric>

namespace Cubiquity
{
	using namespace Internals;

	// Was in Cache.h

	template <typename T>
	struct Entry
	{
		Entry() : key(0) {}
		Entry(uint64_t k, const T& v)
			: key(k)
			, value(v) {}

		uint64_t key;
		T value;
	};

	template <typename T, int size>
	class Cache
	{
	public:
		Cache() : begin(0) {}
		Entry<T>& find(uint64_t key)
		{
			uint32_t i = begin;

			for (int ct = 0; ct < size; ct++)
			{
				Entry<T>& entry = entries[i % size];

				if (entry.key == key)
				{
					return entry;
				}

				i++;
			}

			i--;
			begin = i % size;

			return entries[begin];
		}

	private:
		uint32_t begin;
		std::array<Entry<T>, size> entries;
	};

	// Was in FixedCapacityVector.h

	// An alternative to std::vector (though many functions are not implemented) for
	// which the capacity cannot be changed. Hence it is much simpler and there are no
	// dynamic memory allocations, tracking of sizes, or extra pointer indirections.
	template <typename T, size_t Capacity>
	class FixedCapacityVector
	{
	public:

		typedef size_t size_type;

		FixedCapacityVector() : mEnd(mData.begin()) {} // Start empty

		FixedCapacityVector(const std::initializer_list<T>& dataIn)
			: mEnd(mData.begin()) // Start empty
		{
			// Add the elements from the initializer list
			for (const auto& element : dataIn)
			{
				push_back(element);
			}
		}

		void clear()
		{
			mEnd = mData.begin();
		}

		void push_back(const T& element)
		{
			assert(size() < Capacity);

			if (size() < Capacity)
			{
				*mEnd = element;
				mEnd++;
			}
		}

		size_type capacity() const
		{
			return Capacity;
		}

		T& operator[](size_t index) { assert(index < size()); return mData[index]; }
		const T& operator[](size_t index) const { assert(index < size()); return mData[index]; }

		typename std::array<T, Capacity>::iterator begin() noexcept { return mData.begin(); }
		typename std::array<T, Capacity>::const_iterator begin() const noexcept { return mData.begin(); }

		typename std::array<T, Capacity>::iterator end() noexcept { return mEnd; }
		typename std::array<T, Capacity>::const_iterator end() const noexcept { return mEnd; }

		size_type size() const noexcept { return end() - begin(); }

	private:
		std::array<T, Capacity> mData;
		typename std::array<T, Capacity>::iterator mEnd;
	};

	// Was in VisibilityMask.h
	namespace Operations
	{
		enum Operation
		{
			Draw,
			Test,
		};
	}
	typedef Operations::Operation Operation;

	struct Edge
	{
		// Doesn't intialise as this slows down creation of PolygonEdgeArray - but does that matter?
		Edge() {}

		Edge(int p0in, int p1in)
			:p0(p0in)
			, p1(p1in) {}

		bool operator==(const Edge& other) const {
			return
				p0 == other.p0 && p1 == other.p1;
		}

		bool operator!=(const Edge& other) const {
			return !(*this == other);
		}

		// Note: Could use bit field here as we only need 3 bits per value. Also having a
		// union with an int would then allow testing both members for zero at the same time.
		int p0;
		int p1;
	};

	struct Bounds
	{
		Vector2i lower;
		Vector2i upper;
	};

	typedef FixedCapacityVector<Edge, 12> PolygonEdgeArray;

	typedef std::array<Vector2i, 8> PolygonVertexArray;

	class VisibilityMask
	{
	public:
		VisibilityMask(uint32_t width, uint32_t height);
		~VisibilityMask();

		uint32_t width() { return mWidth; }
		uint32_t height() { return mHeight; }

		void clear();
		void setOpaque();

		uint32_t hash();

		bool pointInRect(const Vector2i& c, const Vector2i& clippedLowerLeft, const Vector2i& clippedUpperRight);
		bool pointInPolygon(const Vector2i& pointToTest, const PolygonVertexArray& vertices, const PolygonEdgeArray& edges);
		bool isConvex(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges);

		void drawPixel(uint32_t x, uint32_t y);
		bool testPixel(uint32_t x, uint32_t y);

		void drawConvexPolygon(const PolygonVertexArray& vertices);
		void drawConvexPolygon(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		void drawConvexPolygonTiny(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		void drawConvexPolygonSmall(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		void drawConvexPolygonLarge(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);

		bool testConvexPolygon(const PolygonVertexArray& vertices);
		bool testConvexPolygon(const PolygonVertexArray& vertices, const Bounds& bounds);
		bool testConvexPolygon(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);

		bool testConvexPolygonTiny(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		bool testConvexPolygonSmall(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		bool testConvexPolygonLarge(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);

		void drawConvexPolygonReference(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);
		bool testConvexPolygonReference(const PolygonVertexArray& vertices, const PolygonEdgeArray& segments, const Bounds& bounds);

		uint32_t getFaceSize();

	private:

		uint64_t hashCubeCorners(const PolygonVertexArray& vertices);

		bool allEdgesAreDiagonal(const PolygonVertexArray& vertices, const PolygonEdgeArray& edges);

		uint32_t mWidth = 0;
		uint32_t mWidthInTiles = 0;
		uint32_t mHeight = 0;
		uint32_t mHeightInTiles = 0;
		uint64_t* mTiles = nullptr;

		Cache<std::array<std::array<bool, 16>, 16>, 16> mSmallPolygonCache;
	};

	Bounds computeBounds(const PolygonVertexArray& vertices);

	int32_t perpDotProduct(const Vector2i& a, const Vector2i& b, const Vector2i& c);
	float perpDotProduct(const Vector4f& a, const Vector4f& b, const Vector4f& c);

	void segmentsFromCorners(const PolygonVertexArray& verticesIn, PolygonEdgeArray& edgesOut);
	void computeConvexHull(const PolygonVertexArray& verticesIn, PolygonEdgeArray& edgesOut);
	bool isInside(const Vector3f& pos, float x, float y, float z, float size);
	bool isInsideInt(const Vector3f& pos, int32_t x, int32_t y, int32_t z, int64_t size);

	// Was in Renderer.h

	struct BoundingIndices
	{
		BoundingIndices()
			// Initialised to invalid values
			:minX(Invalid), minY(Invalid)
			, maxX(Invalid), maxY(Invalid) {}

		BoundingIndices(uint8_t min_x, uint8_t min_y, uint8_t max_x, uint8_t max_y)
			:minX(min_x), minY(min_y)
			, maxX(max_x), maxY(max_y) {}

		bool operator==(const BoundingIndices& other) const {
			return
				minX == other.minX && minY == other.minY &&
				maxX == other.maxX && maxY == other.maxY;
		}

		bool operator!=(const BoundingIndices& other) const {
			return !(*this == other);
		}

		uint8_t minX;
		uint8_t minY;
		uint8_t maxX;
		uint8_t maxY;

		static const uint8_t Invalid = 255;
	};

	class CameraData
	{
	public:
		CameraData()
			:mPos(0.0f, 0.0f, 0.0f)
			, mViewMat()
			, mProjMat() {}
		CameraData(const Vector3f& position, const Vector3f& centre, const Vector3f& up, float fovy, float aspect)
		{
			mPos = position;
			mProjMat = perspective_matrix(fovy, aspect, 0.1f, 100.0f);
			mViewMat = lookAtRH(position, centre, up);
		}

		const Matrix4x4f& viewMatrix() const { return mViewMat; }
		const Matrix4x4f& projMatrix() const { return mProjMat; }
		const Vector3f& position() const { return mPos; }

	private:
		Vector3f mPos;
		Matrix4x4f mViewMat;
		Matrix4x4f mProjMat;
	};

	namespace ClipResults
	{
		enum ClipResult
		{
			NotClipped,
			Clipped,
			Culled
		};
	}
	typedef ClipResults::ClipResult ClipResult;

	namespace Modes
	{
		enum Mode
		{
			Normal,
			Cubemap,
		};
	}
	typedef Modes::Mode Mode;

	extern const PolygonEdgeArray seg[9];
	extern PolygonEdgeArray segmentCache[27];
	extern const BoundingIndices boundingIndices[9];

	int getSectionIndex(const Vector3f& cubeCentreCameraSpace, float size);

	int determineView(const PolygonVertexArray& vertices);

	inline double roundAwayFromZero(double x)
	{
		return x < 0 ? std::floor(x) : std::ceil(x);
	}

	//template<Operation op> bool doCubeAgainstDirectional(VisibilityMask* visMask, const CameraData& cameraData, const Vector3f& cubeCentreWorldSpace, float cubeSize)
	template<Operation op> bool doCubeAgainstDirectional(VisibilityMask* visMask, const CameraData& cameraData, const Vector3i& cubeLowerCornerWorldSpace, uint64_t cubeSizeInt)
	{
		// HACK - Might help by skipping large nodes which have precision problems when rendering? To be investigated...
		if (cubeSizeInt > 1000)
		{
			return true;
		}

		const float cubeSize = cubeSizeInt;
		const float halfCubeSize = cubeSize * 0.5f;
		const Vector3f& cubeCentreWorldSpace = static_cast<Vector3f>(cubeLowerCornerWorldSpace) + halfCubeSize;

		const Vector3f& cameraPosWorldSpace = cameraData.position();
		const Vector3f cubeCentreCameraSpace = cubeCentreWorldSpace - cameraPosWorldSpace;
		(void)(cubeCentreCameraSpace); // Suppress 'unused' warnings when not using asserts.

		// We do not allow a cube to be rendered if the camera is inside it. It would not really make sense to do so, because
		// when drawing we would be writing to every pixel and when testing we would be reading from every pixel. If we really
		// wanted to work with every pixel then we should have a seperate routine for this which just does a memcpy().
		assert((cubeCentreCameraSpace.x() > halfCubeSize) || (cubeCentreCameraSpace.y() > halfCubeSize) || (cubeCentreCameraSpace.z() > halfCubeSize) ||
			(cubeCentreCameraSpace.x() < -halfCubeSize) || (cubeCentreCameraSpace.y() < -halfCubeSize) || (cubeCentreCameraSpace.z() < -halfCubeSize));

		Vector4f cubeCentreViewSpace;
		const Matrix4x4f& viewMatrix = cameraData.viewMatrix();
		// Perform the view transform. The following lines of code are equililent to:
		//
		//     cubeCentreViewSpace = mul(viewMatrix, Vector4f(cubeCentreWorldSpace, 1.0f));
		//
		// But are optimised by expanding out the matrix multiplication and
		// skipping elements which are always zero or one.
		cubeCentreViewSpace[0] =
			cubeCentreWorldSpace.x() * viewMatrix[0][0] +
			cubeCentreWorldSpace.y() * viewMatrix[1][0] +
			cubeCentreWorldSpace.z() * viewMatrix[2][0] +
			/* 1.0 * */ viewMatrix[3][0];
		cubeCentreViewSpace[1] =
			cubeCentreWorldSpace.x() * viewMatrix[0][1] +
			cubeCentreWorldSpace.y() * viewMatrix[1][1] +
			cubeCentreWorldSpace.z() * viewMatrix[2][1] +
			/* 1.0 * */ viewMatrix[3][1];
		cubeCentreViewSpace[2] =
			cubeCentreWorldSpace.x() * viewMatrix[0][2] +
			cubeCentreWorldSpace.y() * viewMatrix[1][2] +
			cubeCentreWorldSpace.z() * viewMatrix[2][2] +
			/* 1.0 * */ viewMatrix[3][2];
		cubeCentreViewSpace[3] = 1.0f;

		// If the cube is completely behind the camera then it is not visible.
		// Note: Shouldn't we be checking against half the cube diagonal here?
		/*if (cubeCentreViewSpace.z() > halfCubeSize)
		{
		return false;
		}*/

		Vector4f xAxis = viewMatrix[0] * cubeSize;
		Vector4f yAxis = viewMatrix[1] * cubeSize;
		Vector4f zAxis = viewMatrix[2] * cubeSize;
		//Vector4f yAxis = viewMatrix[1] * cubeSize;
		//Vector4f zAxis = viewMatrix[2] * cubeSize;

		Vector4f cubeLowerCornerViewSpace = cubeCentreViewSpace;
		cubeLowerCornerViewSpace -= xAxis * 0.5f;
		cubeLowerCornerViewSpace -= yAxis * 0.5f;
		cubeLowerCornerViewSpace -= zAxis * 0.5f;

		PolygonVertexArray corners2DInt;
		//Vector2f corners2DFloat[8];
		Vector4f corners[8];

		corners[0] = cubeLowerCornerViewSpace;
		corners[1] = corners[0] + xAxis;
		corners[2] = corners[0] + yAxis;
		corners[3] = corners[1] + yAxis;
		corners[4] = corners[0] + zAxis;
		corners[5] = corners[1] + zAxis;
		corners[6] = corners[2] + zAxis;
		corners[7] = corners[3] + zAxis;

		bool allBehindCamera = true;
		for (int i = 0; i < 7; i++)
		{
			if (corners[i].z() < 0.0f) // Camera looks down negative z
			{
				allBehindCamera = false;
				break;
			}
		}
		if (allBehindCamera)
		{
			return false;
		}

		const Matrix4x4f& projMatrix = cameraData.projMatrix();

		assert(visMask->getFaceSize() % 2 == 0); // So we know that half the face size will be an int
		const int halfFaceSize = visMask->getFaceSize() / 2;
		const float halfFaceSizeAsFloat = static_cast<float>(halfFaceSize);

		int ct = 0;
		for (uint32_t zCt = 0; zCt < 2; zCt++)
		{
			for (uint32_t yCt = 0; yCt < 2; yCt++)
			{
				for (uint32_t xCt = 0; xCt < 2; xCt++)
				{
					// This is our rather crude handling of objects which intersect the near clip plane.
					// We just force all vertices to be slightly in front. This wil distort the cube, but
					// I don't think it affects the footprint (needs testing). If we find it does affect
					// the footprint, we could also try just discarding this vertex and computing the
					// convex polygon from those which remain (this means we need a way of flagging
					// vertices as valid though. Or perhaps convex calculation code could just ignore
					// vertices with positive z?).
					// FIXME - This has a significan effect on performance, at least on Windows. But why?
					// Could skip it if we know a given cube can't intersect the near clip plane?
					// Update, the performance hit seems to be gone following some code refactoring?
					corners[ct].data[2] = std::min(corners[ct].z(), -0.1f);


					// Apply the projection. We can also just divide by '-z' here (note the minus) which still does a projection but not
					// one which matches the OpenGL pipeline. Instead it gives a 90 degree field of view, which is usually more than our
					// real camera and so may be useful for catching stuff which is about to come onscreen (as well as being faster).
					// Actually, perhaps we can scale the 'z' value (or something?) to make the field of view larger than 90 degrees?

					/*std::cout << projMatrix << std::endl;
					exit(1);*/

					/*corners[ct] = mul(projMatrix, corners[ct]);
					corners[ct] = corners[ct] / corners[ct].w;*/

					corners[ct].data[0] *= projMatrix[0][0];
					corners[ct].data[1] *= projMatrix[1][1];

					const float invZ = 1.0f / -corners[ct].z();
					corners[ct].data[0] *= invZ;
					corners[ct].data[1] *= invZ;

					corners[ct] *= halfFaceSizeAsFloat;
					corners[ct] += halfFaceSizeAsFloat;

					// For some reason the vector opration below is much slower (on VS2017)
					// than doing the two seperate additions. Needs some investigation.
					//corners[ct] += 0.5f; // Redundant!
					corners[ct].data[0] += 0.5f;
					corners[ct].data[1] += 0.5f;

					//if (op == Operations::Draw)
					{
						corners2DInt[ct].data[0] = static_cast<int>(corners[ct].x() + 0.5f);
						corners2DInt[ct].data[1] = static_cast<int>(corners[ct].y() + 0.5f);
					}
					/*else if (op == Operations::Test)
					{
					corners2DInt[ct].x = roundAwayFromZero(corners[ct].x);
					corners2DInt[ct].y = roundAwayFromZero(corners[ct].y);
					}*/

					ct++;
				}
			}
		}

		// Rasterize
		if (op == Operations::Draw)
		{
			visMask->drawConvexPolygon(corners2DInt);
		}
		else if (op == Operations::Test)
		{
			bool visible = visMask->testConvexPolygon(corners2DInt);

			if (visible)
				return true;
		}

		return false;
	}

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

	// Alternative name ideas: Vis(ibility)Query/Test(er)/Finder
	class VisibilityCalculator
	{
		struct NodeState
		{
			NodeState()
				: mCentre(0.0f, 0.0f, 0.0f)
				, mLowerCorner(0, 0, 0)
				, mCentreX2(0, 0, 0)
				, mFootprintSize(0.0f)
				, mIndex(0xFFFFFFFF)
				, mProcessedChildCount(0)
				, mNearestChild(0) {}

			Vector3f mCentre;
			Vector3i mLowerCorner;
			Vector3i64 mCentreX2;
			float mFootprintSize;

			uint32_t mIndex;
			uint32_t mProcessedChildCount;

			uint8_t mNearestChild; // Index (0 - 7) of nearest child.
		};

	public:

		VisibilityCalculator();
		~VisibilityCalculator();

		Glyph buildGlyphFromNode(float centreX, float centreY, float centreZ, uint32_t size, Node nodeForNormal, uint32_t nodeIndex, const Volume* volume, const Vector3f& cameraPos);

		uint32_t findVisibleOctreeNodesPerspective(CameraData* cameraData, const Volume* volume, Glyph* glyphs, uint32_t maxGlyphCount);

		VisibilityMask* cubeFace(int i) { return mCubeFaces[i]; }

		uint32_t mMaxLevelsToGenerate;
		float mMaxFootprintSize;
		float mMinDrawSize;
		float mMinTestSize;
		float mMaxNodeDistance;

	protected:

		VisibilityMask* mCubeFaces[6];
	};

	uint32_t getMaterialForNode(float centreX, float centreY, float centreZ, uint32_t nodeIndex, Volume* volume, const Vector3f& cameraPos);
	Vector3f computeNodeNormal(Node node);
	Vector3f computeNormal(float x, float y, float z, uint32_t size, Volume* volume);

	Vector3f fakeNormalFromSize(int size);

	//template<Operation op> bool doCubeAgainstDirectional(VisibilityMask* visMask, const CameraData& cameraData, const Vector3f& cubeCentreWorldSpace, float cubeSize);
	template<Operation op> bool doCubeAgainstDirectional(VisibilityMask* visMask, const CameraData& cameraData, const Vector3i& cubeLowerCornerWorldSpace, uint64_t cubeSizeInt);

	void computeBounds(const PolygonVertexArray& vertices, int32_t& min_x, int32_t& min_y, int32_t& max_x, int32_t& max_y, uint32_t width);
	BoundingIndices computeBoundingIndices(Vector2f vertices[8]);

	// Raytracing
	double traceRay(Ray3f ray, Volume& volume);

	struct RayVolumeIntersection
	{
		explicit operator bool() { return material > 0; }

		Vector3d position;
		Vector3d normal;
		double distance;
		MaterialId material;
	};

	RayVolumeIntersection ray_parameter(const Volume& volume, Ray3d ray);
}

#endif //CUBIQUITY_RENDERING_H
