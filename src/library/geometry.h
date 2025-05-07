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
#include <iostream>
#include <random>
#include <span>

namespace Cubiquity
{

// Returns -1, 0, or 1. See https ://stackoverflow.com/a/4609795
template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

////////////////////////////////////////////////////////////////////////////////
//                                    Vector
////////////////////////////////////////////////////////////////////////////////

template <class Type>
struct Vec3
{
	typedef Type value_type;
	static constexpr int size() { return 3; }

	// Constructors
	Vec3() {}
	Vec3(const Type& val): x(val), y(val), z(val) {}
	Vec3(const Type& x, const Type& y, const Type& z) :x(x), y(y), z(z) {}

	// Converting constructor
	template <typename SrcType>
	explicit Vec3(const Vec3<SrcType>& v) :x(v.x), y(v.y), z(v.z) {}

	// Default comparison operators (requires C++20)
	friend auto operator<=>(const Vec3&, const Vec3&) = default;

	// Overloaded operators following the guidelines here:
	// https://en.cppreference.com/w/cpp/language/operators
	// https://learn.microsoft.com/en-us/cpp/cpp/operator-overloading
	Type& operator[](int index)
	{
		static_assert(sizeof(Vec3<Type>) == size() * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}
	const Type& operator[](int index) const
	{
		static_assert(sizeof(Vec3<Type>) == size() * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}

	friend std::ostream& operator<<(std::ostream& os, const Vec3& v)
	{
		os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
		return os;
	}

	// Unary negation returns by value - see https://stackoverflow.com/a/2155307
	Vec3 operator-() const { return { -x, -y, -z }; }

	// Binary operators with assignment and scalar on right-hand side
	Vec3& operator+= (const Type& rhs) { x +=  rhs; y +=  rhs; z +=  rhs; return *this; }
	Vec3& operator-= (const Type& rhs) { x -=  rhs; y -=  rhs; z -=  rhs; return *this; }
	Vec3& operator*= (const Type& rhs) { x *=  rhs; y *=  rhs; z *=  rhs; return *this; }
	Vec3& operator/= (const Type& rhs) { x /=  rhs; y /=  rhs; z /=  rhs; return *this; }
	Vec3& operator%= (const Type& rhs) { x %=  rhs; y %=  rhs; z %=  rhs; return *this; }

	Vec3& operator&= (const Type& rhs) { x &=  rhs; y &=  rhs; z &=  rhs; return *this; }
	Vec3& operator|= (const Type& rhs) { x |=  rhs; y |=  rhs; z |=  rhs; return *this; }
	Vec3& operator^= (const Type& rhs) { x ^=  rhs; y ^=  rhs; z ^=  rhs; return *this; }
	Vec3& operator<<=(const Type& rhs) { x <<= rhs; y <<= rhs; z <<= rhs; return *this; }
	Vec3& operator>>=(const Type& rhs) { x >>= rhs; y >>= rhs; z >>= rhs; return *this; }

	// Binary operators with assignment and vector on right-hand side
	Vec3& operator+= (const Vec3& rhs) { x +=  rhs.x; y +=  rhs.y; z +=  rhs.z; return *this; }
	Vec3& operator-= (const Vec3& rhs) { x -=  rhs.x; y -=  rhs.y; z -=  rhs.z; return *this; }
	Vec3& operator*= (const Vec3& rhs) { x *=  rhs.x; y *=  rhs.y; z *=  rhs.z; return *this; }
	Vec3& operator/= (const Vec3& rhs) { x /=  rhs.x; y /=  rhs.y; z /=  rhs.z; return *this; }
	Vec3& operator%= (const Vec3& rhs) { x %=  rhs.x; y %=  rhs.y; z %=  rhs.z; return *this; }

	Vec3& operator&= (const Vec3& rhs) { x &=  rhs.x; y &=  rhs.y; z &=  rhs.z; return *this; }
	Vec3& operator|= (const Vec3& rhs) { x |=  rhs.x; y |=  rhs.y; z |=  rhs.z; return *this; }
	Vec3& operator^= (const Vec3& rhs) { x ^=  rhs.x; y ^=  rhs.y; z ^=  rhs.z; return *this; }
	Vec3& operator<<=(const Vec3& rhs) { x <<= rhs.x; y <<= rhs.y; z <<= rhs.z; return *this; }
	Vec3& operator>>=(const Vec3& rhs) { x >>= rhs.x; y >>= rhs.y; z >>= rhs.z; return *this; }

	// Binary operators with scalar on right-hand side
	friend Vec3 operator+ (Vec3 lhs, const Type& rhs) { return lhs += rhs; }
	friend Vec3 operator- (Vec3 lhs, const Type& rhs) { return lhs -= rhs; }
	friend Vec3 operator* (Vec3 lhs, const Type& rhs) { return lhs *= rhs; }
	friend Vec3 operator/ (Vec3 lhs, const Type& rhs) { return lhs /= rhs; }
	friend Vec3 operator% (Vec3 lhs, const Type& rhs) { return lhs %= rhs; }

	friend Vec3 operator& (Vec3 lhs, const Type& rhs) { return lhs &= rhs; }
	friend Vec3 operator| (Vec3 lhs, const Type& rhs) { return lhs |= rhs; }
	friend Vec3 operator^ (Vec3 lhs, const Type& rhs) { return lhs ^= rhs; }
	friend Vec3 operator<<(Vec3 lhs, const Type& rhs) { return lhs <<= rhs; }
	friend Vec3 operator>>(Vec3 lhs, const Type& rhs) { return lhs >>= rhs; }

	// Binary operators with vector on right-hand side
	friend Vec3 operator+ (Vec3 lhs, const Vec3& rhs) { return lhs += rhs; }
	friend Vec3 operator- (Vec3 lhs, const Vec3& rhs) { return lhs -= rhs; }
	friend Vec3 operator* (Vec3 lhs, const Vec3& rhs) { return lhs *= rhs; }
	friend Vec3 operator/ (Vec3 lhs, const Vec3& rhs) { return lhs /= rhs; }
	friend Vec3 operator% (Vec3 lhs, const Vec3& rhs) { return lhs %= rhs; }

	friend Vec3 operator& (Vec3 lhs, const Vec3& rhs) { return lhs &= rhs; }
	friend Vec3 operator| (Vec3 lhs, const Vec3& rhs) { return lhs |= rhs; }
	friend Vec3 operator^ (Vec3 lhs, const Vec3& rhs) { return lhs ^= rhs; }
	friend Vec3 operator<<(Vec3 lhs, const Vec3& rhs) { return lhs <<= rhs; }
	friend Vec3 operator>>(Vec3 lhs, const Vec3& rhs) { return lhs >>= rhs; }

	// Other vector operations
	friend Vec3 abs(const Vec3& v) { return { std::abs(v.x), std::abs(v.y), std::abs(v.z) }; }
	friend bool all(const Vec3& v) { return v.x && v.y && v.z; }
	friend bool any(const Vec3& v) { return v.x || v.y || v.z; }
	friend Vec3 ceil(const Vec3& v) { return { std::ceil(v.x), std::ceil(v.y), std::ceil(v.z) }; }
	friend Vec3 cross(const Vec3& a, const Vec3& b) {
		return { a.y * b.z - a.z * b.y,
				 a.z * b.x - a.x * b.z,
				 a.x * b.y - a.y * b.x };
	}
	friend Type dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
	friend Vec3 floor(const Vec3& v) { return { std::floor(v.x), std::floor(v.y), std::floor(v.z) }; }
	friend Vec3 fract(const Vec3& v) { return v - floor(v); }
	friend Type length(const Vec3& v) {
		static_assert(std::is_floating_point_v<Type>);
		return sqrt(dot(v, v));
	}
	friend Vec3 max(const Vec3& v0, const Vec3& v1) {
		return { std::max(v0.x, v1.x), std::max(v0.y, v1.y), std::max(v0.z, v1.z) };
	}
	friend int max_index(const Vec3& v) {
		if (v.x >= v.y && v.x >= v.z) return 0;
		if (v.y >= v.x && v.y >= v.z) return 1;
		return 2;
	}
	friend Type max_value(const Vec3& v)
	{
		return std::max(std::max(v.x, v.y), v.z);
	}
	friend Vec3 min(const Vec3& v0, const Vec3& v1) {
		return { std::min(v0.x, v1.x), std::min(v0.y, v1.y), std::min(v0.z, v1.z) };
	}
	friend int min_index(const Vec3& v) {
		if (v.x <= v.y && v.x <= v.z) return 0;
		if (v.y <= v.x && v.y <= v.z) return 1;
		return 2;
	}
	friend Type min_value(const Vec3& v)
	{
		return std::min(std::min(v.x, v.y), v.z);
	}
	friend Vec3 normalize(const Vec3& v) {
		static_assert(std::is_floating_point_v<Type>);
		//assert(length(v) >= 0.001f);
		return v / length(v);
	}
	friend Vec3 pow(const Vec3& v, const Type& t) {
		return { std::pow(v.x, t), std::pow(v.y, t), std::pow(v.z, t) };
	}
	friend Vec3 sign(const Vec3& v) {
		return { sign(v.x), sign(v.y), sign(v.z) };
	}
	friend Vec3 step(float edge, const Vec3& v) {
		return { v.x < edge ? 0 : 1, v.y < edge ? 0 : 1, v.z < edge ? 0 : 1 };
	}

	// GLSL-compatible comparisons
	friend Vec3<bool> equal(const Vec3& v0, const Vec3& v1) {
		return { v0.x == v1.x, v0.y == v1.y, v0.z == v1.z };
	}
	friend Vec3<bool> lessThan(const Vec3& v0, const Vec3& v1) {
		return { v0.x < v1.x, v0.y < v1.y, v0.z < v1.z };
	}
	friend Vec3<bool> lessThanEqual(const Vec3& v0, const Vec3& v1) {
		return { v0.x <= v1.x, v0.y <= v1.y, v0.z <= v1.z };
	}
	friend Vec3<bool> greaterThan(const Vec3& v0, const Vec3& v1) {
		return { v0.x > v1.x, v0.y > v1.y, v0.z > v1.z };
	}
	friend Vec3<bool> greaterThanEqual(const Vec3& v0, const Vec3& v1) {
		return { v0.x >= v1.x, v0.y >= v1.y, v0.z >= v1.z };
	}

	Type x, y, z;
};

// Typedefs for basic vector types
typedef Vec3<i32> vec3i;
typedef Vec3<u32> vec3u;
typedef Vec3<bool> vec3b;
typedef Vec3<float> vec3f;
typedef Vec3<double> vec3d;

template <class VecType>
class Ray3
{
public:
	Ray3() {}

	Ray3(const VecType& origin, const VecType& dir)
		: mOrigin(origin)
		, mDir(dir)	{}

	template <typename CastType> explicit Ray3(const Ray3<CastType>& ray)
	{
		//for (u32 ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
		mOrigin = static_cast<VecType>(ray.mOrigin);
		mDir = static_cast<VecType>(ray.mDir);
	}

public:
	VecType mOrigin;
	VecType mDir; // FIXME - Would float always be sufficient for the direction?
};

typedef Ray3<vec3f> Ray3f;
typedef Ray3<vec3d> Ray3d;

template <class VecType>
class Box
{
public:

	// Default constructed box starts off as invalid (min bigger than max))
	Box() { invalidate(); }

	Box(const VecType& lower, const VecType& upper)
		:mExtents{ lower, upper } {}

	// For casting betwween Box types of matching size.
	template <typename CastType> explicit Box(const Box<CastType>& box)
		:mExtents { 
			static_cast<VecType>(box.lower()),
			static_cast<VecType>(box.upper())
		} {}

	const VecType& lower() const { return mExtents[0]; }
	const VecType& upper() const { return mExtents[1]; }

	// Should try to remove these non-const versions...
	VecType& lower() { return mExtents[0]; }
	VecType& upper() { return mExtents[1]; }

	void accumulate(const VecType& value)
	{
		lower() = min(lower(), value);
		upper() = max(upper(), value);
	}

	// Note: Could templatise on container
	void accumulate(const std::array<VecType, 3>& points)
	{
		for (auto const& point : points) {
			accumulate(point);
		}
	}

	void accumulate(const Box& other)
	{
		lower() = min(lower(), other.lower());
		upper() = max(upper(), other.upper());
	}

	bool contains(const VecType& value) const
	{
		return (value.x >= lower().x) && (value.y >= lower().y) && (value.z >= lower().z) &&
			(value.x <= upper().x) && (value.y <= upper().y) && (value.z <= upper().z);
	}

	bool contains(const Box& other) const
	{
		return contains(other.lower()) && contains(other.upper());
	}

	//void dilate(typename VecType::value_type amount)
	void dilate(VecType amount)
	{
		//VecType amountAsVec = { amount, amount, amount };
		VecType amountAsVec = amount;
		lower() -= amountAsVec;
		upper() += amountAsVec;
	}

	void invalidate()
	{
		lower() = VecType(std::numeric_limits<typename VecType::value_type>::max());
		upper() = VecType(std::numeric_limits<typename VecType::value_type>::lowest());
	}

	bool isValid()
	{
		return lower().x <= upper().x &&
			lower().y <= upper().y &&
			lower().z <= upper().z;
	}

	static Box invalid()
	{
		return Box(); // Default-constructed box is already invalid.
	}

public:
	VecType mExtents[2];
};

template <class VecType>
bool operator==(const Box<VecType>& lhs, const Box<VecType>& rhs)
{
	return lhs.lower() == rhs.lower() && lhs.upper() == rhs.upper();
}

template <class VecType>
bool overlaps(const Box<VecType>& a, const Box<VecType>& b)
{
	for (int i = 0; i < 3; i++)
	{
		if (a.upper()[i] < b.lower()[i] || a.lower()[i] > b.upper()[i])
		{
			return false;
		}
	}
	return true;
}

template <class VecType>
class Box3 : public Box<VecType>
{
public:
	Box3() : Box<VecType>() {}
	Box3(const VecType& lower, const VecType& upper) : Box<VecType>(lower, upper) {}
	explicit Box3(const Box<vec3i>& box) : Box<VecType>(box) {}

	vec3f centre() const { return (this->lower() + this->upper()) * 0.5f; }

	typename VecType::value_type volume() const
	{
		// FIXME - Handle negative volumes?
		VecType dims = this->upper() - this->lower();
		return dims.x * dims.y * dims.z;
	}
};

typedef Box3<vec3f> Box3f;
typedef Box3<vec3d> Box3d;

class Box3i : public Box<vec3i>
{
public:
	Box3i() : Box<vec3i>() {}
	Box3i(const vec3i& lower, const vec3i& upper) : Box<vec3i>(lower, upper) {}
	explicit Box3i(const Box<vec3f>& box) : Box<vec3i>(box) {}

	static Box3i max()
	{
		vec3i lower = vec3i(std::numeric_limits<int>::lowest());
		vec3i upper = vec3i(std::numeric_limits<int>::max());
		return Box3i(lower, upper);
	}

	// Note that these integer dimensions add one to include last voxel, and
	// return integers types which are wider than the internal representation.
	i64 width()  const { return static_cast<i64>(upper()[0]) - static_cast<i64>(lower()[0]) + 1; }
	i64 height() const { return static_cast<i64>(upper()[1]) - static_cast<i64>(lower()[1]) + 1; }
	i64 depth()  const { return static_cast<i64>(upper()[2]) - static_cast<i64>(lower()[2]) + 1; }

	// FIXME - Handle overflow.
	i64 voxelCount() const { return width() * height() * depth(); }
};

class Box3fSampler
{
public:
	Box3fSampler(const Box3f& bounds)
		: eng(0)
		, randX(bounds.lower().x, bounds.upper().x)
		, randY(bounds.lower().y, bounds.upper().y)
		, randZ(bounds.lower().z, bounds.upper().z) {}

	vec3f next() {
		vec3f result = { (float)randX(eng), (float)randY(eng), (float)randZ(eng) };
		return result;
	}

private:

	std::minstd_rand eng;
	std::uniform_real_distribution<> randX;
	std::uniform_real_distribution<> randY;
	std::uniform_real_distribution<> randZ;
};

class Triangle
{
public:
	Triangle();
	Triangle(float v0x, float v0y, float v0z, float v1x, float v1y, float v1z, float v2x, float v2y, float v2z)
		:Triangle(vec3f(v0x, v0y, v0z), vec3f(v1x, v1y, v1z), vec3f(v2x, v2y, v2z)) {}
	Triangle(const Cubiquity::vec3f& vertex0, const Cubiquity::vec3f& vertex1, const Cubiquity::vec3f& vertex2)
		: vertices{ { vertex0, vertex1, vertex2 } } {}

	void flip();
	void translate(const Cubiquity::vec3f& dir);
	void scale(float factor);

	float sideLength(int index) const;
	vec3f computeNormal() const;
	float area() const;
	Cubiquity::vec3f centre() const;

	std::array<Cubiquity::vec3f, 3> vertices;
};

float distance(const Cubiquity::vec3f& point, const Cubiquity::Triangle& triangle);

/// A (conceptual) list of non-indexed triangles, actually stored as an std::vector
typedef std::vector<Cubiquity::Triangle> TriangleList;
typedef std::span<Triangle> TriangleSpan;
typedef std::span<Triangle const> ConstTriangleSpan; // See https://stackoverflow.com/a/56895806

Box3f computeBounds(const std::array<Cubiquity::vec3f, 3>& points);
Box3f computeBounds(ConstTriangleSpan triangles);
void translate(TriangleList& triangles, const Cubiquity::vec3f& dir);
void scale(TriangleList& triangles, float factor);

// Intersections
struct RayBoxIntersection
{
	explicit operator bool() { return exit >= entry; }

	float entry;
	float exit;
};

// See https://tavianator.com/fast-branchless-raybounding-box-intersections/
template <class VecType>
RayBoxIntersection intersect(const Ray3<VecType>& ray, const Box3<VecType>& box)
{
	// Inverse direction could be precomputed and stored in the ray
	// if we find we often intersect the same ray with multiple boxes.
	const VecType invDir = VecType({ 1.0f, 1.0f, 1.0f }) / ray.mDir;

	const VecType lower = (box.lower() - ray.mOrigin) * invDir;
	const VecType upper = (box.upper() - ray.mOrigin) * invDir;

	const VecType minCorner = min(lower, upper);
	const VecType maxCorner = max(lower, upper);

	RayBoxIntersection intersection;
	intersection.entry = max_value(minCorner);
	intersection.exit = min_value(maxCorner);
	return intersection;
}

bool intersect(const Ray3f& ray, const Triangle& triangle, float& t);
	
}

#endif // CUBIQUITY_GEOMETRY_H
