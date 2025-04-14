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
const float Pi = 3.14159265358979f;

// See https ://stackoverflow.com/a/4609795
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
		return { std::copysign(1.0f, v.x),
			     std::copysign(1.0f, v.y),
			     std::copysign(1.0f, v.z) };
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

template <class Type>
struct Vec4
{
	typedef Type value_type;
	static constexpr int size() { return 4; }

	// Constructors
	Vec4() {}
	Vec4(const Type& val) : x(val), y(val), z(val), w(val) {}
	Vec4(const Type& x, const Type& y, const Type& z, const Type& w)
		:x(x), y(y), z(z), w(w) {}

	// Converting constructor
	template <typename SrcType>
	explicit Vec4(const Vec4<SrcType>& v) :x(v.x), y(v.y), z(v.z), w(v.w) {}

	// Default comparison operators (requires C++20)
	friend auto operator<=>(const Vec4&, const Vec4&) = default;

	// Overloaded operators following the guidelines here:
	// https://en.cppreference.com/w/cpp/language/operators
	// https://learn.microsoft.com/en-us/cpp/cpp/operator-overloading
	Type& operator[](int index)
	{
		static_assert(sizeof(Vec4<Type>) == size() * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}
	const Type& operator[](int index) const
	{
		static_assert(sizeof(Vec4<Type>) == size() * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}

	friend std::ostream& operator<<(std::ostream& os, const Vec4& v)
	{
		os << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
		return os;
	}

	// Unary negation - see https://stackoverflow.com/a/2155307
	Vec4 operator-() const { return {-x, -y, -z, -w}; }

	// Binary operators with assignment and scalar on right-hand side
	Vec4& operator+= (const Type& rhs) { x +=  rhs; y +=  rhs; z +=  rhs; w +=  rhs; return *this; }
	Vec4& operator-= (const Type& rhs) { x -=  rhs; y -=  rhs; z -=  rhs; w -=  rhs; return *this; }
	Vec4& operator*= (const Type& rhs) { x *=  rhs; y *=  rhs; z *=  rhs; w *=  rhs; return *this; }
	Vec4& operator/= (const Type& rhs) { x /=  rhs; y /=  rhs; z /=  rhs; w /=  rhs; return *this; }
	Vec4& operator%= (const Type& rhs) { x %=  rhs; y %=  rhs; z %=  rhs; w %=  rhs; return *this; }

	Vec4& operator&= (const Type& rhs) { x &=  rhs; y &=  rhs; z &=  rhs; w &=  rhs; return *this; }
	Vec4& operator|= (const Type& rhs) { x |=  rhs; y |=  rhs; z |=  rhs; w |=  rhs; return *this; }
	Vec4& operator^= (const Type& rhs) { x ^=  rhs; y ^=  rhs; z ^=  rhs; w ^=  rhs; return *this; }
	Vec4& operator<<=(const Type& rhs) { x <<= rhs; y <<= rhs; z <<= rhs; w <<= rhs; return *this; }
	Vec4& operator>>=(const Type& rhs) { x >>= rhs; y >>= rhs; z >>= rhs; w >>= rhs; return *this; }

	// Binary operators with assignment and vector on right-hand side
	Vec4& operator+= (const Vec4& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }
	Vec4& operator-= (const Vec4& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }
	Vec4& operator*= (const Vec4& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w; return *this; }
	Vec4& operator/= (const Vec4& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; w /= rhs.w; return *this; }
	Vec4& operator%= (const Vec4& rhs) { x %= rhs.x; y %= rhs.y; z %= rhs.z; w %= rhs.w; return *this; }

	Vec4& operator&= (const Vec4& rhs) { x &=  rhs.x; y &=  rhs.y; z &=  rhs.z; w &=  rhs.w; return *this; }
	Vec4& operator|= (const Vec4& rhs) { x |=  rhs.x; y |=  rhs.y; z |=  rhs.z; w |=  rhs.w; return *this; }
	Vec4& operator^= (const Vec4& rhs) { x ^=  rhs.x; y ^=  rhs.y; z ^=  rhs.z; w ^=  rhs.w; return *this; }
	Vec4& operator<<=(const Vec4& rhs) { x <<= rhs.x; y <<= rhs.y; z <<= rhs.z; w <<= rhs.w; return *this; }
	Vec4& operator>>=(const Vec4& rhs) { x >>= rhs.x; y >>= rhs.y; z >>= rhs.z; w >>= rhs.w; return *this; }

	// Binary operators with scalar on right-hand side
	friend Vec4 operator+ (Vec4 lhs, const Type& rhs) { return lhs += rhs; }
	friend Vec4 operator- (Vec4 lhs, const Type& rhs) { return lhs -= rhs; }
	friend Vec4 operator* (Vec4 lhs, const Type& rhs) { return lhs *= rhs; }
	friend Vec4 operator/ (Vec4 lhs, const Type& rhs) { return lhs /= rhs; }
	friend Vec4 operator% (Vec4 lhs, const Type& rhs) { return lhs %= rhs; }

	friend Vec4 operator& (Vec4 lhs, const Type& rhs) { return lhs &= rhs; }
	friend Vec4 operator| (Vec4 lhs, const Type& rhs) { return lhs |= rhs; }
	friend Vec4 operator^ (Vec4 lhs, const Type& rhs) { return lhs ^= rhs; }
	friend Vec4 operator<<(Vec4 lhs, const Type& rhs) { return lhs <<=rhs; }
	friend Vec4 operator>>(Vec4 lhs, const Type& rhs) { return lhs >>=rhs; }

	// Binary operators with vector on right-hand side
	friend Vec4 operator+ (Vec4 lhs, const Vec4& rhs) { return lhs += rhs; }
	friend Vec4 operator- (Vec4 lhs, const Vec4& rhs) { return lhs -= rhs; }
	friend Vec4 operator* (Vec4 lhs, const Vec4& rhs) { return lhs *= rhs; }
	friend Vec4 operator/ (Vec4 lhs, const Vec4& rhs) { return lhs /= rhs; }
	friend Vec4 operator% (Vec4 lhs, const Vec4& rhs) { return lhs %= rhs; }

	friend Vec4 operator& (Vec4 lhs, const Vec4& rhs) { return lhs &= rhs; }
	friend Vec4 operator| (Vec4 lhs, const Vec4& rhs) { return lhs |= rhs; }
	friend Vec4 operator^ (Vec4 lhs, const Vec4& rhs) { return lhs ^= rhs; }
	friend Vec4 operator<<(Vec4 lhs, const Vec4& rhs) { return lhs <<=rhs; }
	friend Vec4 operator>>(Vec4 lhs, const Vec4& rhs) { return lhs >>=rhs; }

	Type x, y, z, w;
};

// Typedefs for basic vector types
typedef Vec3<int32> vec3i;
typedef Vec3<uint32> vec3u;
typedef Vec3<bool> vec3b;
typedef Vec3<float> vec3f;
typedef Vec3<double> vec3d;

typedef Vec4<int32> vec4i;
typedef Vec4<uint32> vec4u;
typedef Vec4<bool> vec4b;
typedef Vec4<float> vec4f;
typedef Vec4<double> vec4d;

////////////////////////////////////////////////////////////////////////////////
//									Matrix
////////////////////////////////////////////////////////////////////////////////
template <typename Type>
class Matrix4x4
{
public:
	// Set to identity matrix.
	Matrix4x4()
	{
		for(int y = 0; y < 4; y++)
		{
			for(int x = 0; x < 4; x++)
			{
				data[x][y] = x == y ? 1.0f : 0.0f;
			}
		}
	}

	Matrix4x4(const Vec4<Type>& v0, const Vec4<Type>& v1, const Vec4<Type>& v2, const Vec4<Type>& v3)
	{
		data[0] = v0;
		data[1] = v1;
		data[2] = v2;
		data[3] = v3;
	}

	// For casting betwween vec types of matching size.
	template <typename CastType> explicit Matrix4x4(const Matrix4x4<CastType>& matrix)
	{
		for (uint32_t ct = 0; ct < 4; ++ct) { data[ct] = static_cast< Vec4<Type> >(matrix.data[ct]); }
	}

	Vec4<Type>& operator[](int index) { return data[index]; }

	const Vec4<Type>& operator[](int index) const { return data[index]; }

	void operator/=(Type const& rhs) { for (int i = 0; i < 4; ++i) { data[i] /= rhs; } }

public:
	Vec4<Type> data[4];
};

typedef Matrix4x4<float> Matrix4x4f;
typedef Matrix4x4<double> Matrix4x4d;

template <typename Type>
Matrix4x4<Type> frustum_matrix(Type x0, Type x1, Type y0, Type y1, Type n, Type f)
{
	const Type s = -1.0f, o = n;

	Matrix4x4<Type> result;
	result[0] = { 2 * n / (x1 - x0), 0, 0, 0 };
	result[1] = { 0, 2 * n / (y1 - y0), 0, 0 };
	result[2] = { -s * (x0 + x1) / (x1 - x0), -s * (y0 + y1) / (y1 - y0), s * (f + o) / (f - n), s };
	result[3] = { 0, 0, -(n + o) * f / (f - n), 0 };
	return result;
}

template <typename Type>
Matrix4x4<Type> perspective_matrix(Type fovy, Type aspect, Type n, Type f)
{
	Type y = n * std::tan(fovy / 2), x = y * aspect;

	return frustum_matrix(-x, x, -y, y, n, f);
}

template <typename Type>
Matrix4x4<Type> lookAtRH(const Vec3<Type>& eye, const Vec3<Type>& center, const Vec3<Type>& up)
{
	const Vec3<Type> f(normalize(center - eye));
	const Vec3<Type> s(normalize(cross(f, up)));
	const Vec3<Type> u(cross(s, f));

	Matrix4x4<Type> result;
	result[0][0] = s.x;
	result[1][0] = s.y;
	result[2][0] = s.z;
	result[0][1] = u.x;
	result[1][1] = u.y;
	result[2][1] = u.z;
	result[0][2] = -f.x;
	result[1][2] = -f.y;
	result[2][2] = -f.z;
	result[3][0] = -dot(s, eye);
	result[3][1] = -dot(u, eye);
	result[3][2] = dot(f, eye);
	return result;
}

template <typename Type>
Matrix4x4<Type> operator/(Matrix4x4<Type> const& lhs, Type const& rhs)
{
	Matrix4x4<Type> result(lhs); result /= rhs; return result;
}

template <typename Type>
Vec4<Type> mul(const Matrix4x4<Type>& a, const Vec4<Type>& b) { return a[0] * b.x + a[1] * b.y + a[2] * b.z + a[3] * b.w; }

template <typename Type>
Matrix4x4<Type> mul(const Matrix4x4<Type>& a, const Matrix4x4<Type>& b)
{
	return { mul(a,b[0]), mul(a,b[1]), mul(a,b[2]), mul(a,b[3]) };
}

template <typename Type>
Matrix4x4<Type> adjugate(const Matrix4x4<Type>& a) {
	return {
		{a[1][1] * a[2][2] * a[3][3] + a[3][1] * a[1][2] * a[2][3] + a[2][1] * a[3][2] * a[1][3] - a[1][1] * a[3][2] * a[2][3] - a[2][1] * a[1][2] * a[3][3] - a[3][1] * a[2][2] * a[1][3],
			a[0][1] * a[3][2] * a[2][3] + a[2][1] * a[0][2] * a[3][3] + a[3][1] * a[2][2] * a[0][3] - a[3][1] * a[0][2] * a[2][3] - a[2][1] * a[3][2] * a[0][3] - a[0][1] * a[2][2] * a[3][3],
			a[0][1] * a[1][2] * a[3][3] + a[3][1] * a[0][2] * a[1][3] + a[1][1] * a[3][2] * a[0][3] - a[0][1] * a[3][2] * a[1][3] - a[1][1] * a[0][2] * a[3][3] - a[3][1] * a[1][2] * a[0][3],
			a[0][1] * a[2][2] * a[1][3] + a[1][1] * a[0][2] * a[2][3] + a[2][1] * a[1][2] * a[0][3] - a[0][1] * a[1][2] * a[2][3] - a[2][1] * a[0][2] * a[1][3] - a[1][1] * a[2][2] * a[0][3]},
		{a[1][2] * a[3][3] * a[2][0] + a[2][2] * a[1][3] * a[3][0] + a[3][2] * a[2][3] * a[1][0] - a[1][2] * a[2][3] * a[3][0] - a[3][2] * a[1][3] * a[2][0] - a[2][2] * a[3][3] * a[1][0],
			a[0][2] * a[2][3] * a[3][0] + a[3][2] * a[0][3] * a[2][0] + a[2][2] * a[3][3] * a[0][0] - a[0][2] * a[3][3] * a[2][0] - a[2][2] * a[0][3] * a[3][0] - a[3][2] * a[2][3] * a[0][0],
			a[0][2] * a[3][3] * a[1][0] + a[1][2] * a[0][3] * a[3][0] + a[3][2] * a[1][3] * a[0][0] - a[0][2] * a[1][3] * a[3][0] - a[3][2] * a[0][3] * a[1][0] - a[1][2] * a[3][3] * a[0][0],
			a[0][2] * a[1][3] * a[2][0] + a[2][2] * a[0][3] * a[1][0] + a[1][2] * a[2][3] * a[0][0] - a[0][2] * a[2][3] * a[1][0] - a[1][2] * a[0][3] * a[2][0] - a[2][2] * a[1][3] * a[0][0]},
		{a[1][3] * a[2][0] * a[3][1] + a[3][3] * a[1][0] * a[2][1] + a[2][3] * a[3][0] * a[1][1] - a[1][3] * a[3][0] * a[2][1] - a[2][3] * a[1][0] * a[3][1] - a[3][3] * a[2][0] * a[1][1],
			a[0][3] * a[3][0] * a[2][1] + a[2][3] * a[0][0] * a[3][1] + a[3][3] * a[2][0] * a[0][1] - a[0][3] * a[2][0] * a[3][1] - a[3][3] * a[0][0] * a[2][1] - a[2][3] * a[3][0] * a[0][1],
			a[0][3] * a[1][0] * a[3][1] + a[3][3] * a[0][0] * a[1][1] + a[1][3] * a[3][0] * a[0][1] - a[0][3] * a[3][0] * a[1][1] - a[1][3] * a[0][0] * a[3][1] - a[3][3] * a[1][0] * a[0][1],
			a[0][3] * a[2][0] * a[1][1] + a[1][3] * a[0][0] * a[2][1] + a[2][3] * a[1][0] * a[0][1] - a[0][3] * a[1][0] * a[2][1] - a[2][3] * a[0][0] * a[1][1] - a[1][3] * a[2][0] * a[0][1]},
		{a[1][0] * a[3][1] * a[2][2] + a[2][0] * a[1][1] * a[3][2] + a[3][0] * a[2][1] * a[1][2] - a[1][0] * a[2][1] * a[3][2] - a[3][0] * a[1][1] * a[2][2] - a[2][0] * a[3][1] * a[1][2],
			a[0][0] * a[2][1] * a[3][2] + a[3][0] * a[0][1] * a[2][2] + a[2][0] * a[3][1] * a[0][2] - a[0][0] * a[3][1] * a[2][2] - a[2][0] * a[0][1] * a[3][2] - a[3][0] * a[2][1] * a[0][2],
			a[0][0] * a[3][1] * a[1][2] + a[1][0] * a[0][1] * a[3][2] + a[3][0] * a[1][1] * a[0][2] - a[0][0] * a[1][1] * a[3][2] - a[3][0] * a[0][1] * a[1][2] - a[1][0] * a[3][1] * a[0][2],
			a[0][0] * a[1][1] * a[2][2] + a[2][0] * a[0][1] * a[1][2] + a[1][0] * a[2][1] * a[0][2] - a[0][0] * a[2][1] * a[1][2] - a[1][0] * a[0][1] * a[2][2] - a[2][0] * a[1][1] * a[0][2]}
	};
}

template <typename Type>
Type determinant(const Matrix4x4<Type>& a) {
	return a[0][0] * (a[1][1] * a[2][2] * a[3][3] + a[3][1] * a[1][2] * a[2][3] + a[2][1] * a[3][2] * a[1][3] - a[1][1] * a[3][2] * a[2][3] - a[2][1] * a[1][2] * a[3][3] - a[3][1] * a[2][2] * a[1][3])
		+ a[0][1] * (a[1][2] * a[3][3] * a[2][0] + a[2][2] * a[1][3] * a[3][0] + a[3][2] * a[2][3] * a[1][0] - a[1][2] * a[2][3] * a[3][0] - a[3][2] * a[1][3] * a[2][0] - a[2][2] * a[3][3] * a[1][0])
		+ a[0][2] * (a[1][3] * a[2][0] * a[3][1] + a[3][3] * a[1][0] * a[2][1] + a[2][3] * a[3][0] * a[1][1] - a[1][3] * a[3][0] * a[2][1] - a[2][3] * a[1][0] * a[3][1] - a[3][3] * a[2][0] * a[1][1])
		+ a[0][3] * (a[1][0] * a[3][1] * a[2][2] + a[2][0] * a[1][1] * a[3][2] + a[3][0] * a[2][1] * a[1][2] - a[1][0] * a[2][1] * a[3][2] - a[3][0] * a[1][1] * a[2][2] - a[2][0] * a[3][1] * a[1][2]);
}

template <typename Type>
Matrix4x4<Type> inverse(const Matrix4x4<Type>& a) { return adjugate(a) / determinant(a); }

template <typename Type>
Matrix4x4<Type> translation_matrix(const Vec3<Type>& pos)
{
	Matrix4x4<Type> result;
	result[3] = Vec4<Type>({ pos.x, pos.y, pos.z, 1.0f });
	return result;
}

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
		//for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
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
	int64 width()  const { return static_cast<int64>(upper()[0]) - static_cast<int64>(lower()[0]) + 1; }
	int64 height() const { return static_cast<int64>(upper()[1]) - static_cast<int64>(lower()[1]) + 1; }
	int64 depth()  const { return static_cast<int64>(upper()[2]) - static_cast<int64>(lower()[2]) + 1; }

	// FIXME - Handle overflow.
	int64 voxelCount() const { return width() * height() * depth(); }
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

class Box3iSampler
{
public:
	Box3iSampler(const Box3i& bounds)
		: eng(0)
		, randX(bounds.lower().x, bounds.upper().x)
		, randY(bounds.lower().y, bounds.upper().y)
		, randZ(bounds.lower().z, bounds.upper().z) {}

	vec3i next() { return vec3i({ randX(eng), randY(eng), randZ(eng) }); }

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
		iterator(uint32_t index, const Box3i& box) : eng(0), randX(box.lower().x, box.upper().x), randY(box.lower().y, box.upper().y), randZ(box.lower().z, box.upper().z), mIndex(index) { mValue = vec3i({ randX(eng), randY(eng), randZ(eng) }); }
		iterator operator++() { ++mIndex; mValue = vec3i({ randX(eng), randY(eng), randZ(eng) }); return *this; }
		bool operator!=(const iterator & other) const { return mIndex != other.mIndex; }
		const vec3i& operator*() const { return mValue; }

	private:
		vec3i mValue;
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
