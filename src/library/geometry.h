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

template <class Type, int Size>
struct vec : std::array<Type, Size>
{
	// Constructors
	vec() {}
	vec(const Type& val)
	{
		// On Visual Studio 2022 unrolling the fill is slightly faster.
		//this->fill(val);
		if constexpr (Size > 0) { (*this)[0] = val; }
		if constexpr (Size > 1) { (*this)[1] = val; }
		if constexpr (Size > 2) { (*this)[2] = val; }
		if constexpr (Size > 3) { (*this)[3] = val; }
	}

	/*explicit vec(const vec& v)
	{
		std::copy(v.begin(), v.end(), this->begin());
	}*/
	
	// Using an std::initializer_list to initialise an underlying array is
	// suprisingly difficult: https://stackoverflow.com/q/5549524
	// I did not get the proposed variadic template solution to work (perhaps
	// because array is a base, not a member?), but the simple copy() works.
	vec(std::initializer_list<Type> l) // : std::array<Type, Size>(l)
	{
		// On Visual Studio 2022 unrolling the copy is slightly faster.
		//std::copy(l.begin(), l.end(), this->begin());
		if constexpr (Size > 0) { (*this)[0] = l.begin()[0]; }
		if constexpr (Size > 1) { (*this)[1] = l.begin()[1]; }
		if constexpr (Size > 2) { (*this)[2] = l.begin()[2]; }
		if constexpr (Size > 3) { (*this)[3] = l.begin()[3]; }
	}

	//vec(std::array<Type, Size> a) : std::array<Type, Size>(a) {}

	// For static_cast support
	template <typename CastType>
	explicit operator vec<CastType, Size>() const
	{
		vec<CastType, Size> result;
		for (uint32_t ct = 0; ct < Size; ++ct) { result[ct] = static_cast<CastType>((*this)[ct]); }
		return result;
	}

	// Named element access
	const Type& x() const { static_assert(Size > 0); return (*this)[0]; }
	const Type& y() const { static_assert(Size > 1); return (*this)[1]; }
	const Type& z() const { static_assert(Size > 2); return (*this)[2]; }
	const Type& w() const { static_assert(Size > 3); return (*this)[3]; }

	/**************************** Support functions ***************************/

	// Apply function in-place to each element of this vector. The explicit
	// unrolling seems to help performance, at least on Visual C++ 2022.
	template<class UnaryFunc>
	constexpr vec& apply(UnaryFunc func)
	{
		static_assert(Size <= 4);
		if constexpr (Size > 0) { func((*this)[0]); }
		if constexpr (Size > 1) { func((*this)[1]); }
		if constexpr (Size > 2) { func((*this)[2]); }
		if constexpr (Size > 3) { func((*this)[3]); }
		return *this;
	}

	// As above, but using the corresponding element of the params vector.
	template<class BinaryFunc>
	constexpr vec& apply(const vec& params, BinaryFunc func)
	{
		static_assert(Size <= 4);
		if constexpr (Size > 0) { func((*this)[0], params[0]); }
		if constexpr (Size > 1) { func((*this)[1], params[1]); }
		if constexpr (Size > 2) { func((*this)[2], params[2]); }
		if constexpr (Size > 3) { func((*this)[3], params[3]); }
		return *this;
	}

	// Similar to 'apply()', but as a non-member making from supplied vector.
	template<class UnaryFunc>
	friend constexpr vec make_from(const vec& v, UnaryFunc func)
	{
		static_assert(Size <= 4);
		vec<Type, Size> result;
		if constexpr (Size > 0) { result[0] = func(v[0]); }
		if constexpr (Size > 1) { result[1] = func(v[1]); }
		if constexpr (Size > 2) { result[2] = func(v[2]); }
		if constexpr (Size > 3) { result[3] = func(v[3]); }
		return result;
	}

	// As above, but taking two vectors as input.
	template<class BinaryFunc>
	friend constexpr auto make_from(const vec& v0,
		                            const vec& v1, BinaryFunc func)
	{
		static_assert(Size <= 4);
		using RetType = decltype(func(v0[0], v1[0])); // Can we avoid this?
		vec<RetType, Size> result;
		if constexpr (Size > 0) { result[0] = func(v0[0],v1[0]); }
		if constexpr (Size > 1) { result[1] = func(v0[1],v1[1]); }
		if constexpr (Size > 2) { result[2] = func(v0[2],v1[2]); }
		if constexpr (Size > 3) { result[3] = func(v0[3],v1[3]); }
		return result;
	}

	/************************** Overloaded operators **************************/

	// Overloaded operators following the guidelines here:
	// https://en.cppreference.com/w/cpp/language/operators
	// https://learn.microsoft.com/en-us/cpp/cpp/operator-overloading

	// Unary operators as non-members
	friend vec operator-(vec v) { return v.apply([](Type& t) {t = -t; }); }
	friend vec operator~(vec v) { return v.apply([](Type& t) {t = ~t; }); }

	// Binary operators with assignment and scalar on right-hand side
	vec& operator+= (const Type& rhs) {
		return apply([rhs](Type& t) {t += rhs; }); }
	vec& operator-= (const Type& rhs) {
		return apply([rhs](Type& t) {t -= rhs; }); }
	vec& operator*= (const Type& rhs) {
		return apply([rhs](Type& t) {t *= rhs; }); }
	vec& operator/= (const Type& rhs) {
		return apply([rhs](Type& t) {t /= rhs; }); }
	vec& operator%= (const Type& rhs) {
		return apply([rhs](Type& t) {t %= rhs; }); }

	vec& operator&= (const Type& rhs) {
		return apply([rhs](Type& t) {t &= rhs; }); }
	vec& operator|= (const Type& rhs) {
		return apply([rhs](Type& t) {t |= rhs; }); }
	vec& operator^= (const Type& rhs) {
		return apply([rhs](Type& t) {t ^= rhs; }); }
	vec& operator>>=(const Type& rhs) {
		return apply([rhs](Type& t) {t >>= rhs;}); }
	vec& operator<<=(const Type& rhs) {
		return apply([rhs](Type& t) {t <<= rhs;}); }

	// Binary operators with assignment and vector on right-hand side
	vec& operator+= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l += r; }); }
	vec& operator-= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l -= r; }); }
	vec& operator*= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l *= r; }); }
	vec& operator/= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l /= r; }); }
	vec& operator%= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l %= r; }); }

	vec& operator&= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l &= r; }); }
	vec& operator|= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l |= r; }); }
	vec& operator^= (const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l ^= r; }); }
	vec& operator>>=(const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l >>= r;}); }
	vec& operator<<=(const vec& rhs) {
		return apply(rhs, [](Type& l, const Type& r) {l <<= r;}); }

	// Binary operators as members with scalar on right-hand side
	friend vec operator+ (vec lhs, const Type& rhs) { return lhs += rhs; }
	friend vec operator- (vec lhs, const Type& rhs) { return lhs -= rhs; }
	friend vec operator* (vec lhs, const Type& rhs) { return lhs *= rhs; }
	friend vec operator/ (vec lhs, const Type& rhs) { return lhs /= rhs; }
	friend vec operator% (vec lhs, const Type& rhs) { return lhs %= rhs; }

	friend vec operator& (vec lhs, const Type& rhs) { return lhs &= rhs; }
	friend vec operator| (vec lhs, const Type& rhs) { return lhs |= rhs; }
	friend vec operator^ (vec lhs, const Type& rhs) { return lhs ^= rhs; }
	friend vec operator>>(vec lhs, const Type& rhs) { return lhs >>= rhs; }
	friend vec operator<<(vec lhs, const Type& rhs) { return lhs <<= rhs; }

	// Binary operators as members with vector on right-hand side
	friend vec operator+ (vec lhs, const vec& rhs) { return lhs += rhs; }
	friend vec operator- (vec lhs, const vec& rhs) { return lhs -= rhs; }
	friend vec operator* (vec lhs, const vec& rhs) { return lhs *= rhs; }
	friend vec operator/ (vec lhs, const vec& rhs) { return lhs /= rhs; }
	friend vec operator% (vec lhs, const vec& rhs) { return lhs %= rhs; }

	friend vec operator& (vec lhs, const vec& rhs) { return lhs &= rhs; }
	friend vec operator| (vec lhs, const vec& rhs) { return lhs |= rhs; }
	friend vec operator^ (vec lhs, const vec& rhs) { return lhs ^= rhs; }
	friend vec operator>>(vec lhs, const vec& rhs) { return lhs >>= rhs; }
	friend vec operator<<(vec lhs, const vec& rhs) { return lhs <<= rhs; }

	// Stream insertion operator.
	friend std::ostream& operator<<(std::ostream& os, const vec& vector) {
		os << "[";
		for (auto v : vector) { os << v << ","; }
		os << "\b]"; // Backspace to remove last comma.
		return os;
	}

	/************************* Other vector operations ************************/

	friend vec abs(vec v)
	{ 
		return make_from(v, [](const Type& t) { return std::abs(t); });
	}

	friend bool all(const vec& x)
	{
		for (int i = 0; i < Size; ++i) { if (!x[i]) return false; }
		return true;
	}

	friend bool any(const vec& x)
	{
		for (int i = 0; i < Size; ++i) { if (x[i]) return true; }
		return false;
	}

	friend vec ceil(vec v)
	{ 
		return make_from(v, [](const Type& t) { return std::ceil(t); });
	}

	friend vec cross(const vec& a, const vec& b)
	{
		static_assert(Size == 3, "Cross product only valid for 3D vectors");
		return { a.y() * b.z() - a.z() * b.y(),
				 a.z() * b.x() - a.x() * b.z(),
				 a.x() * b.y() - a.y() * b.x() };
	}

	friend Type dot(const vec& a, const vec& b)
	{
		Type result = 0;
		for (int i = 0; i < Size; ++i) { result += a[i] * b[i]; }
		return result;
	}

	friend vec floor(vec v)
	{
		return make_from(v, [](const Type& t) { return std::floor(t); });
	}

	friend vec fract(vec v)
	{
		return v - floor(v);
	}

	friend Type length(const vec& v)
	{
		// Floating point componenents only due to sqrt(). Integer vectors can
		// be cast to floating point vectors prior to calling this function.
		static_assert(std::is_floating_point<Type>::value);
		return sqrt(dot(v, v));
	}

	friend vec max(vec v0, const vec& v1)
	{
		return make_from(v0, v1, 
			[](const Type& x, const Type& y) { return std::max(x, y); });
	}

	friend vec min(vec v0, const vec& v1)
	{
		return make_from(v0, v1,
			[](const Type& x, const Type& y) { return std::min(x, y); });
	}

	friend vec mix(const vec& x, const vec& y, const vec& a)
	{
		return x * (vec(1) - a) + y * a;
	}

	friend vec normalize(const vec& v)
	{
		static_assert(std::is_floating_point<Type>::value);
		assert(length(v) >= 0.001f);
		return v / length(v);
	}

	friend vec pow(vec x, const Type& y)
	{ 
		return make_from(x, [y](const Type& t) { return std::pow(t, y); });
	}

	friend vec<long int, Size> round_to_int(const vec& v)
	{
		// Rounding in C++ is suprisingly complex (e.g. lround() vs lrint()) and 
		// built-in functions can be slow (https://stackoverflow.com/q/53962727).
		// Hence we use a simpler method here.
		return static_cast<vec<long int, Size>>(floor(v + vec(0.5)));
	}

	friend vec sign(vec v)
	{
		return v.apply([](Type& t) {t = std::copysign(1.0f, t); });
	}

	friend vec step(float edge, vec v)
	{
		return v.apply([edge](Type& t) {t = t < edge ? 0 : 1; });
	}

	// GLSL-compatible comparisons
	friend vec<bool, Size> equal(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l == r; });
	}

	friend vec<bool, Size> lessThan(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l < r; });
	}

	friend vec<bool, Size> lessThanEqual(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l <= r; });
	}

	friend vec<bool, Size> greaterThan(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l > r; });
	}

	friend vec<bool, Size> greaterThanEqual(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l >= r; });
	}
};

// Typedefs for basic vector types
typedef vec<int32, 2> vec2i;
typedef vec<uint32, 2> vec2u;
typedef vec<bool, 2> vec2b;
typedef vec<float, 2> vec2f;
typedef vec<double, 2> vec2d;

typedef vec<int32, 3> vec3i;
typedef vec<uint32, 3> vec3u;
typedef vec<bool, 3> vec3b;
typedef vec<float, 3> vec3f;
typedef vec<double, 3> vec3d;

typedef vec<int32, 4> vec4i;
typedef vec<uint32, 4> vec4u;
typedef vec<bool, 4> vec4b;
typedef vec<float, 4> vec4f;
typedef vec<double, 4> vec4d;

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

	Matrix4x4(const vec<Type, 4>& v0, const vec<Type, 4>& v1, const vec<Type, 4>& v2, const vec<Type, 4>& v3)
	{
		data[0] = v0;
		data[1] = v1;
		data[2] = v2;
		data[3] = v3;
	}

	// For casting betwween vec types of matching size.
	template <typename CastType> explicit Matrix4x4(const Matrix4x4<CastType>& matrix)
	{
		for (uint32_t ct = 0; ct < 4; ++ct) { data[ct] = static_cast< vec<Type, 4> >(matrix.data[ct]); }
	}

	vec<Type, 4>& operator[](int index) { return data[index]; }

	const vec<Type, 4>& operator[](int index) const { return data[index]; }

	void operator/=(Type const& rhs) { for (int i = 0; i < 4; ++i) { data[i] /= rhs; } }

public:
	vec<Type, 4> data[4];
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
Matrix4x4<Type> lookAtRH(const vec<Type, 3>& eye, const vec<Type, 3>& center, const vec<Type, 3>& up)
{
	const vec<Type, 3> f(normalize(center - eye));
	const vec<Type, 3> s(normalize(cross(f, up)));
	const vec<Type, 3> u(cross(s, f));

	Matrix4x4<Type> result;
	result[0][0] = s.x();
	result[1][0] = s.y();
	result[2][0] = s.z();
	result[0][1] = u.x();
	result[1][1] = u.y();
	result[2][1] = u.z();
	result[0][2] = -f.x();
	result[1][2] = -f.y();
	result[2][2] = -f.z();
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
vec<Type, 4> mul(const Matrix4x4<Type>& a, const vec<Type, 4>& b) { return a[0] * b.x() + a[1] * b.y() + a[2] * b.z() + a[3] * b.w(); }

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
Matrix4x4<Type> translation_matrix(const vec<Type, 3>& pos)
{
	Matrix4x4<Type> result;
	result[3] = vec<Type, 4>({ pos.x(), pos.y(), pos.z(), 1.0f });
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

template <class Type, int Size>
class Ray
{
public:
	Ray() {}

	Ray(const vec<Type, Size>& origin, const vec<Type, Size>& dir)
		: mOrigin(origin)
		, mDir(dir)	{}

	template <typename CastType> explicit Ray(const Ray<CastType, Size>& ray)
	{
		//for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
		mOrigin = static_cast<vec<Type, Size>>(ray.mOrigin);
		mDir = static_cast<vec<Type, Size>>(ray.mDir);
	}

public:
	vec<Type, Size> mOrigin;
	vec<Type, Size> mDir; // FIXME - Would float always be sufficient for the direction?
};

typedef Ray<float, 3> Ray3f;
typedef Ray<double, 3> Ray3d;

template <class Type, int Size>
class Box
{
public:

	// Default constructed box starts off as invalid (min bigger than max))
	Box() { invalidate(); }

	Box(const vec<Type, Size>& lower, const vec<Type, Size>& upper)
		:mExtents{ lower, upper } {}

	// For casting betwween Box types of matching size.
	template <typename CastType> explicit Box(const Box<CastType, Size>& box)
		:mExtents { 
			static_cast<vec<Type, Size>>(box.lower()),
			static_cast<vec<Type, Size>>(box.upper())
		} {}

	const vec<Type, Size>& lower() const { return mExtents[0]; }
	const vec<Type, Size>& upper() const { return mExtents[1]; }

	// Should try to remove these non-const versions...
	vec<Type, Size>& lower() { return mExtents[0]; }
	vec<Type, Size>& upper() { return mExtents[1]; }

	void accumulate(const vec<Type, Size>& value)
	{
		lower() = min(lower(), value);
		upper() = max(upper(), value);
	}

	// Note: Could templatise on container
	void accumulate(const std::array<vec<Type, Size>, 3>& points)
	{
		for (auto const& point : points) {
			accumulate(point);
		}
	}

	void accumulate(const Box<Type, Size>& other)
	{
		lower() = min(lower(), other.lower());
		upper() = max(upper(), other.upper());
	}

	bool contains(const vec<Type, Size>& value) const
	{
		return (value.x() >= lower().x()) && (value.y() >= lower().y()) && (value.z() >= lower().z()) &&
			(value.x() <= upper().x()) && (value.y() <= upper().y()) && (value.z() <= upper().z());
	}

	bool contains(const Box<Type, Size>& other) const
	{
		return contains(other.lower()) && contains(other.upper());
	}

	void dilate(Type amount)
	{
		vec<Type, Size> amountAsVec = { amount, amount, amount };
		lower() -= amountAsVec;
		upper() += amountAsVec;
	}

	void invalidate()
	{
		lower().fill(std::numeric_limits<Type>::max());
		upper().fill(std::numeric_limits<Type>::lowest());
	}

	bool isValid()
	{
		return lower().x() <= upper().x() &&
			lower().y() <= upper().y() &&
			lower().z() <= upper().z();
	}

	static Box<Type, Size> invalid()
	{
		return Box<Type, Size>(); // Default-constructed box is already invalid.
	}

public:
	vec<Type, Size> mExtents[2];
};

template <class Type, int Size>
bool operator==(const Box<Type, Size>& lhs, const Box<Type, Size>& rhs)
{
	return lhs.lower() == rhs.lower() && lhs.upper() == rhs.upper();
}

template <class Type, int Size>
bool overlaps(const Box<Type, Size>& a, const Box<Type, Size>& b)
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

template <class Type>
class Box3 : public Box<Type, 3>
{
public:
	Box3() : Box<Type, 3>() {}
	Box3(const vec<Type, 3>& lower, const vec<Type, 3>& upper) : Box<Type, 3>(lower, upper) {}
	explicit Box3(const Box<int32, 3>& box) : Box<Type, 3>(box) {}

	vec3f centre() const { return (this->lower() + this->upper()) * 0.5f; }

	Type volume() const
	{
		// FIXME - Handle negative volumes?
		vec<Type, 3> dims = this->upper() - this->lower();
		return dims.x() * dims.y() * dims.z();
	}
};

typedef Box3<float> Box3f;
typedef Box3<double> Box3d;

class Box3i : public Box<int32, 3>
{
public:
	Box3i() : Box<int32, 3>() {}
	Box3i(const vec<int32, 3>& lower, const vec<int32, 3>& upper) : Box<int32, 3>(lower, upper) {}
	explicit Box3i(const Box<float, 3>& box) : Box<int32, 3>(box) {}

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
		, randX(bounds.lower().x(), bounds.upper().x())
		, randY(bounds.lower().y(), bounds.upper().y())
		, randZ(bounds.lower().z(), bounds.upper().z()) {}

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
		, randX(bounds.lower().x(), bounds.upper().x())
		, randY(bounds.lower().y(), bounds.upper().y())
		, randZ(bounds.lower().z(), bounds.upper().z()) {}

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
		iterator(uint32_t index, const Box3i& box) : eng(0), randX(box.lower().x(), box.upper().x()), randY(box.lower().y(), box.upper().y()), randZ(box.lower().z(), box.upper().z()), mIndex(index) { mValue = vec3i({ randX(eng), randY(eng), randZ(eng) }); }
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
template <class Type>
RayBoxIntersection intersect(const Ray<Type, 3>& ray, const Box3<Type>& box)
{
	// Inverse direction could be precomputed and stored in the ray
	// if we find we often intersect the same ray with multiple boxes.
	const vec<Type, 3> invDir = vec<Type, 3>({ 1.0f, 1.0f, 1.0f }) / ray.mDir;

	const vec<Type, 3> lower = (box.lower() - ray.mOrigin) * invDir;
	const vec<Type, 3> upper = (box.upper() - ray.mOrigin) * invDir;

	const vec<Type, 3> minCorner = min(lower, upper);
	const vec<Type, 3> maxCorner = max(lower, upper);

	RayBoxIntersection intersection;
	intersection.entry = *(std::max_element(minCorner.begin(), minCorner.end()));
	intersection.exit = *(std::min_element(maxCorner.begin(), maxCorner.end()));
	return intersection;
}

bool intersect(const Ray3f& ray, const Triangle& triangle, float& t);
	
}

#endif // CUBIQUITY_GEOMETRY_H
