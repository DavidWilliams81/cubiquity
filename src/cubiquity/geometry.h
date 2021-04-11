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
#ifndef CUBIQUITY_GEOMETRY_H
#define CUBIQUITY_GEOMETRY_H

#include "base.h"

#include <algorithm>
#include <array>
#include <cassert>

#include <cstdint>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <random>

namespace Cubiquity
{
	const float Pi = 3.14159265358979f;

	////////////////////////////////////////////////////////////////////////////////
	//									Vector
	////////////////////////////////////////////////////////////////////////////////
	template <int Size, typename Type>
	class Vector
	{
	public:

        // Constructors
		Vector() {}
		Vector(Type t) { for(int i = 0; i < Size; ++i) { data[i] = t; } }
		Vector(Type x, Type y) : data{x,y} { static_assert(Size == 2, "Wrong no. of params"); }
		Vector(Type x, Type y, Type z) : data{x,y,z} { static_assert(Size == 3, "Wrong no. of params"); }
		Vector(Type x, Type y, Type z, Type w) : data{x,y,z,w} { static_assert(Size == 4, "Wrong no. of params"); }

		// For casting betwween Vector types of matching size.
		template <typename CastType> explicit Vector(const Vector<Size, CastType>& vector)
		{
			for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
		}

		// Indexed element access
		Type& operator[](int index) { assert(index < Size && "Index out of range"); return this->data[index]; }
		const Type& operator[](int index) const { assert(index < Size && "Index out of range"); return this->data[index]; }

		// Named element access. Could replace asserts with 'enable_if'?
		const Type& x() const { static_assert(Size > 0, "Vector too small to call x()"); return data[0]; }
		const Type& y() const { static_assert(Size > 1, "Vector too small to call y()"); return data[1]; }
		const Type& z() const { static_assert(Size > 2, "Vector too small to call z()"); return data[2]; }
		const Type& w() const { static_assert(Size > 3, "Vector too small to call w()"); return data[3]; }

		const Vector<3, Type> xxx() const { return Vector<3, Type>(x(), x(), x()); }
		const Vector<3, Type> xyz() const { return Vector<3, Type>(x(), y(), z()); }
		const Vector<3, Type> www() const { return Vector<3, Type>(w(), w(), w()); }

		// Find the biggest and smallest components in the vector
		int minComponentIndex() const
		{
			int minIndex = 0;
			for (int i = 1; i < Size; ++i)
			{
				if (data[i] < data[minIndex])
				{
					minIndex = i;
				}
			}
			return minIndex;
		}

		int maxComponentIndex() const
		{
			int maxIndex = 0;
			for (int i = 1; i < Size; ++i)
			{
				if (data[i] > data[maxIndex])
				{
					maxIndex = i;
				}
			}
			return maxIndex;
		}

		Type minComponentValue() const { Type t = data[0]; for(int i = 1; i < Size; ++i) { t = std::min(t, data[i]); } return t; }
		Type maxComponentValue() const { Type t = data[0]; for(int i = 1; i < Size; ++i) { t = std::max(t, data[i]); } return t; }

		// Unary arithmetic operators
		void operator+=(Type const& rhs) { for(int i = 0; i < Size; ++i) { data[i] += rhs; } }
		void operator+=(Vector<Size, Type> const& rhs) { for(int i = 0; i < Size; ++i) { data[i] += rhs.data[i]; } }
		void operator-=(Type const& rhs) { for(int i = 0; i < Size; ++i) { data[i] -= rhs; } }
		void operator-=(Vector<Size, Type> const& rhs) { for(int i = 0; i < Size; ++i) { data[i] -= rhs.data[i]; } }
		void operator*=(Type const& rhs) { for(int i = 0; i < Size; ++i) { data[i] *= rhs; } }
		void operator*=(Vector<Size, Type> const& rhs) { for(int i = 0; i < Size; ++i) { data[i] *= rhs.data[i]; } }
		void operator/=(Type const& rhs) { for(int i = 0; i < Size; ++i) { data[i] /= rhs; } }
		void operator/=(Vector<Size, Type> const& rhs) { for(int i = 0; i < Size; ++i) { data[i] /= rhs.data[i]; } }

		Type data[Size];
	};

	// Typedefs for commonly-used types
	typedef Vector<2, int> Vector2i;
	typedef Vector<2, float> Vector2f;
	typedef Vector<2, double> Vector2d;

    typedef Vector<3, int> Vector3i;
	typedef Vector<3, uint8_t> Vector3u8;
	typedef Vector<3, int64_t> Vector3i64;
	typedef Vector<3, float> Vector3f;
	typedef Vector<3, double> Vector3d;

	typedef Vector<4, int> Vector4i;
	typedef Vector<4, float> Vector4f;

    // Comparison operator
	template <int Size, typename Type>
	bool operator==(const Vector<Size, Type>& lhs, const Vector<Size, Type>& rhs)
	{
		// Note: The code below does seem sligtly slower than doing:
		//		'return memcmp(lhs.data, rhs.data, sizeof(lhs.data)) == 0;'
		// Is that expected? We probably shouldn't memcmp on floating point values,
		// but should we specialise for integers? More testing is needed here.
		for (int i = 0; i < Size; ++i) { if (lhs[i] != rhs[i]) return false; } return true;
	}

	template <int Size, typename Type>
	bool operator<(const Vector<Size, Type>& lhs, const Vector<Size, Type>& rhs)
	{
		for (int i = 0; i < Size; ++i)
		{
			if (lhs[i] < rhs[i]) return true;
			if (lhs[i] > rhs[i]) return false;
		}
		return false; // Equality
	}

	// Binary arithmetic operators
	template <int Size, typename Type>
	Vector<Size, Type> operator+(Vector<Size, Type> const& lhs, Type const& rhs)
	{ Vector<Size, Type> result( lhs ); result += rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator+(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{ Vector<Size, Type> result( lhs ); result += rhs; return result; }

    template <int Size, typename Type>
	Vector<Size, Type> operator-(Vector<Size, Type> const& lhs, Type const& rhs)
	{ Vector<Size, Type> result( lhs ); result -= rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator-(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{ Vector<Size, Type> result( lhs ); result -= rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator*(Vector<Size, Type> const& lhs, Type const& rhs)
	{ Vector<Size, Type> result( lhs ); result *= rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator*(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{ Vector<Size, Type> result( lhs ); result *= rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator/(Vector<Size, Type> const& lhs, Type const& rhs)
	{ Vector<Size, Type> result( lhs ); result /= rhs; return result; }

	template <int Size, typename Type>
	Vector<Size, Type> operator/(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{ Vector<Size, Type> result( lhs ); result /= rhs; return result; }

	/// Stream insertion operator.
	template <int Size, typename Type>
	std::ostream& operator<<(std::ostream& os, const Vector<Size, Type>& vector)
	{
		os << "(";
		for (uint32 i = 0; i < Size; ++i)
		{
			os << vector[i];
			if (i < (Size - 1))
			{
				os << ",";
			}
		}
		os << ")";
		return os;
	}

	// Component-wise absolute values (output is a Vector rather than a scaler).
	template <int Size, typename Type>
	Vector<Size, Type> abs(Vector<Size, Type> const& vec)
	{
		Vector<Size, Type> result;
		for(int i = 0; i < Size; ++i) { result[i] = std::abs(vec[i]); }
		return result;
	}

	// Component-wise minimum and maximum (output is a Vector rather than a scaler).
	template <int Size, typename Type>
	Vector<Size, Type> min(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{
		Vector<Size, Type> result;
		for(int i = 0; i < Size; ++i) { result[i] = std::min(lhs[i], rhs[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> max(Vector<Size, Type> const& lhs, Vector<Size, Type> const& rhs)
	{
		Vector<Size, Type> result;
		for(int i = 0; i < Size; ++i) { result[i] = std::max(lhs[i], rhs[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> floor(Vector<Size, Type> const& vec)
	{
		Vector<Size, Type> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::floor(vec[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> ceil(Vector<Size, Type> const& vec)
	{
		Vector<Size, Type> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::ceil(vec[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> round(Vector<Size, Type> const& vec)
	{
		Vector<Size, Type> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::round(vec[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, long> lround(Vector<Size, Type> const& vec)
	{
		Vector<Size, long> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::lround(vec[i]); }
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> fract(Vector<Size, Type> const& vec)
	{
		Vector<Size, Type> result;
		double intpart;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::modf(vec[i], &intpart);
		}
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> clamp(Vector<Size, Type> const& vec, Type minVal, Type maxVal)
	{
		Vector<Size, Type> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::max(minVal, std::min(vec[i], maxVal));
		}
		return result;
	}

	template <int Size, typename Type>
	Vector<Size, Type> mix(Vector<Size, Type> const& x, Vector<Size, Type> const& y, float a)
	{
		return x * (1.0f - a) + y * a;
	}

	// Common vector operations
	template <int Size, typename Type>
	Type dot(const Vector<Size, Type>& a, const Vector<Size, Type>& b)
	{
		Type result(0);
		for(int i = 0; i < Size; ++i) { result += a[i] * b[i]; }
		return result;
	}

	// Should only be valid for 3D vectors. Also, what should the return type be?
	template <int Size, typename Type>
	Vector3f cross(const Vector<Size, Type>& a, const Vector<Size, Type>& b)
	{
		return { a.y()*b.z() - a.z()*b.y(), a.z()*b.x() - a.x()*b.z(), a.x()*b.y() - a.y()*b.x() };
	}

	// FIXME - Should this always return float?
	template <int Size, typename Type>
	float length(const Vector<Size, Type>& vec)
	{
		return sqrt(dot(vec, vec));
	}

	template <int Size, typename Type>
	Vector<Size, Type> normalize(const Vector<Size, Type>& vec)
	{
		Vector<Size, Type> result;
		float len = length(vec);

		// Could do better error handling here but it's hard to define a correct
		// behaviour. I think I prefer to just not pass in invalid vectors?
		//assert(len >= 0.001f);

		result = vec * (Type(1.0) / len);
		return result;
	}

	////////////////////////////////////////////////////////////////////////////////
	//									Matrix
	////////////////////////////////////////////////////////////////////////////////
	class Matrix4x4f
	{
	public:
		// Set to identity matrix.
		Matrix4x4f()
		{
			for(int y = 0; y < 4; y++)
			{
				for(int x = 0; x < 4; x++)
				{
					data[x][y] = x == y ? 1.0f : 0.0f;
				}
			}
		}

		Vector4f& operator[](int index) { return data[index]; }

		const Vector4f& operator[](int index) const { return data[index]; }

	private:
		Vector4f data[4];
	};

	Matrix4x4f frustum_matrix(float x0, float x1, float y0, float y1, float n, float f);

	Matrix4x4f perspective_matrix(float fovy, float aspect, float n, float f);

	Matrix4x4f lookAtRH(const Vector3f& eye, const Vector3f& center, const Vector3f& up);

	Matrix4x4f translation_matrix(const Vector3f& pos);


	// Serialization
	/*template<class Type, int Size>
	std::ostream& operator<<(std::ostream& os, const linalg::vec<Type, Size>& vector)
	{
		os << "[";
		for (uint32_t i = 0; i < Size; i++)
		{
			os << vector[i];
			if (i < (Size - 1))
			{
				os << ",";
			}
		}
		os << "]";
		return os;
	}

	template<class Type, int Rows, int Cols>
	std::ostream& operator<<(std::ostream& os, const linalg::mat<Type, Rows, Cols>& matrix)
	{
		for (uint32_t row = 0; row < Rows; row++)
		{
			os << "[";
			for (uint32_t col = 0; col < Cols; col++)
			{
				// Note that 'setw()' only affects the next insertion, so we can
				// use it here for alignment without messing up the global state.
				// See https://stackoverflow.com/q/1532640/2337254
				os << std::setw(12);
				os << matrix[col][row];

				if (col < (Cols - 1))
				{
					os << ",";
				}
			}
			os << "]";
			if (row < (Rows - 1))
			{
				os << std::endl;
			}
		}

		return os;
	}*/

	template <int Size, typename Type>
	class Ray
	{
	public:
		Ray(const Vector<Size, Type>& origin, const Vector<Size, Type>& dir)
			: mOrigin(origin)
			, mDir(dir)	{}

		template <typename CastType> explicit Ray(const Ray<Size, CastType>& ray)
		{
			//for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
			mOrigin = static_cast<Vector<Size, Type>>(ray.mOrigin);
			mDir = static_cast<Vector<Size, Type>>(ray.mDir);
		}

	public:
		Vector<Size, Type> mOrigin;
		Vector<Size, Type> mDir; // FIXME - Would float always be sufficient for the direction?
	};

	typedef Ray<2, float> Ray2f;
	typedef Ray<2, double> Ray2d;
	typedef Ray<3, float> Ray3f;
	typedef Ray<3, double> Ray3d;

	template <int Size, typename Type>
	class Box
	{
	public:

		// Default constructed box starts off as invalid (min bigger than max))
		Box() { invalidate(); }

		Box(const Vector<Size, Type>& value)
		{
			mLower = value;
			mUpper = value;
		}

		Box(const Vector<Size, Type>& lower, const Vector<Size, Type>& upper)
		{
			mLower = lower;
			mUpper = upper;
		}

		// For casting betwween Box types of matching size.
		template <typename CastType> explicit Box(const Box<Size, CastType>& box)
		{
			mLower = static_cast<Vector<Size, Type>>(box.lower());
			mUpper = static_cast<Vector<Size, Type>>(box.upper());
		}

		const Vector<Size, Type>& lower() const { return mLower; }
		const Vector<Size, Type>& upper() const { return mUpper; }

		// Should try to remove these non-const versions...
		Vector<Size, Type>& lower() { return mLower; }
		Vector<Size, Type>& upper() { return mUpper; }

		Vector<Size, float> centre() const
		{
			return static_cast< Vector<Size, float> >(mLower + mUpper) * 0.5f;
		}


		void accumulate(const Vector<Size, Type>& value)
		{
			mLower = Cubiquity::min(mLower, value);
			mUpper = Cubiquity::max(mUpper, value);
		}

		void accumulate(const Box<Size, Type>& other)
		{
			mLower = Cubiquity::min(mLower, other.mLower);
			mUpper = Cubiquity::max(mUpper, other.mUpper);
		}

		bool contains(const Vector<Size, Type>& value) const
		{
			return (value.x() >= mLower.x()) && (value.y() >= mLower.y()) && (value.z() >= mLower.z()) &&
				(value.x() <= mUpper.x()) && (value.y() <= mUpper.y()) && (value.z() <= mUpper.z());
		}

		void dilate(Type amount)
		{
			mLower -= Vector<Size, Type>(amount);
			mUpper += Vector<Size, Type>(amount);
		}

		void invalidate()
		{
			mLower = Vector<Size, Type>(std::numeric_limits<Type>::max());
			mUpper = Vector<Size, Type>(std::numeric_limits<Type>::lowest());
		}

		Type sideLength(uint32_t sideIndex) const
		{
			// FIXME - The '+1' is to get the 'inclusive' length... what do we
			// really want here? Also how should negative side lengths be handled?
			return (mUpper[sideIndex] - mLower[sideIndex]) + 1;
		}
		
		float diagonalLength() const
		{
			return sqrtf(sideLength(0) * sideLength(0) + sideLength(1) *
				sideLength(1) + sideLength(2) * sideLength(2));
		}
		static Box<Size, Type> invalid()
		{
			return Box<Size, Type>(); // Default-constructed box is already invalid.
		}

		static Box<Size, Type> max()
		{
			return Box<Size, Type>(std::numeric_limits<Type>::lowest(), std::numeric_limits<Type>::max());
		}

	private:
		Vector<Size, Type> mLower;
		Vector<Size, Type> mUpper;
	};

	template <int Size, typename Type>
	bool operator==(const Box<Size, Type>& lhs, const Box<Size, Type>& rhs)
	{
		return lhs.lower() == rhs.lower() && lhs.upper() == rhs.upper();
	}

	template <int Size, typename Type>
	bool overlaps(const Box<Size, Type>& a, const Box<Size, Type>& b)
	{
		for (int i = 0; i < Size; i++)
		{
			if (a.upper()[i] < b.lower()[i] || a.lower()[i] > b.upper()[i])
			{
				return false;
			}
		}
		return true;
	}

	typedef Box<2, int> Box2i;
	typedef Box<2, float> Box2f;

	typedef Box<3, int> Box3i;
	typedef Box<3, float> Box3f;

	class Box3fSampler
	{
	public:
		Box3fSampler(const Box3f& bounds)
			: eng(0)
			, randX(bounds.lower().x(), bounds.upper().x())
			, randY(bounds.lower().y(), bounds.upper().y())
			, randZ(bounds.lower().z(), bounds.upper().z()) {}

		Vector3f next() { return Vector3f(randX(eng), randY(eng), randZ(eng)); }

	private:

		std::minstd_rand eng;
		std::uniform_real_distribution<> randX;
		std::uniform_real_distribution<> randY;
		std::uniform_real_distribution<> randZ;
	};

	class Box3iSampler
	{
	public:
		Box3iSampler(const Box3i& bounds)
			: eng(0)
			, randX(bounds.lower().x(), bounds.upper().x())
			, randY(bounds.lower().y(), bounds.upper().y())
			, randZ(bounds.lower().z(), bounds.upper().z()) {}

		Vector3i next() { return Vector3i(randX(eng), randY(eng), randZ(eng)); }

	private:

		std::minstd_rand eng;
		std::uniform_int_distribution<> randX;
		std::uniform_int_distribution<> randY;
		std::uniform_int_distribution<> randZ;
	};

	class Box3iSampler2
	{
	public:
		class iterator
		{
		public:
			iterator(uint32_t index, const Box3i& box) : eng(0), randX(box.lower().x(), box.upper().x()), randY(box.lower().y(), box.upper().y()), randZ(box.lower().z(), box.upper().z()), mIndex(index) { mValue = Vector3i(randX(eng), randY(eng), randZ(eng)); }
			iterator operator++() { ++mIndex; mValue = Vector3i(randX(eng), randY(eng), randZ(eng)); return *this; }
			bool operator!=(const iterator & other) const { return mIndex != other.mIndex; }
			const Vector3i& operator*() const { return mValue; }

		private:
			Vector3i mValue;
			std::minstd_rand eng;
			std::uniform_int_distribution<> randX;
			std::uniform_int_distribution<> randY;
			std::uniform_int_distribution<> randZ;
			uint32_t mIndex;
		};

		Box3iSampler2(uint32_t count, const Box3i& box) : mCount(count), mBox(box) {}

		iterator begin() const { return iterator(0, mBox); }
		iterator end() const { return iterator(mCount, mBox); }

	private:
		unsigned mCount;
		Box3i mBox;
	};

	class Triangle
	{
	public:
		Triangle();
		Triangle(const Cubiquity::Vector3f& vertex0, const Cubiquity::Vector3f& vertex1, const Cubiquity::Vector3f& vertex2)
			: vertices{ { vertex0, vertex1, vertex2 } } {}

		void flip();
		void translate(const Cubiquity::Vector3f& dir);
		void scale(float factor);

		float sideLength(uint32 index) const;
		Vector3f computeNormal() const;

		std::array<Cubiquity::Vector3f, 3> vertices;
	};

	float distance(const Cubiquity::Vector3f& point, const Cubiquity::Triangle& triangle);

	/// A (conceptual) list of non-indexed triangles, actually stored as an std::vector
	typedef std::vector<Cubiquity::Triangle> TriangleList;

	Box3f computeBounds(const TriangleList& triangles);
	void translate(TriangleList& triangles, const Cubiquity::Vector3f& dir);
	void scale(TriangleList& triangles, float factor);

	typedef std::pair<MaterialId, TriangleList> SubObject;
	typedef std::list<SubObject> SubObjectList;

	SubObject mergeSubObjects(const SubObjectList& subObjects, uint16 resultingMaterial);

	class Object
	{
	public:
		std::string name;
		SubObjectList subObjects;
	};

	TriangleList mergedTriangles(const Object& object);

	Box3f computeBounds(const Object& object);

	/// A simple geometry representation which holds a set of material identifiers and a corresponding list of triangles.
	typedef std::list<Object> Geometry;

	Object mergeObjects(const Geometry& geometry, std::string resultingObjectName, uint16 resultingMaterial);

	TriangleList mergedTriangles(const Geometry& geometry);

	Box3f computeBounds(const Geometry& geometry);
	void translate(Geometry& geometry, const Cubiquity::Vector3f& dir);
	void scale(Geometry& geometry, float factor);

	// Intersections
	struct RayBoxIntersection
	{
		explicit operator bool() { return exit >= entry; }

		float entry;
		float exit;
	};
	RayBoxIntersection intersect(const Ray3f& ray, const Box3f& box);
	
}

#endif // CUBIQUITY_GEOMETRY_H
