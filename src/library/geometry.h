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

template <class Type>
struct XY
{
	XY() {}
	XY(const Type& x, const Type& y)
		:x(x), y(y) {}

	static constexpr int sz() { return 2; }

	Type x, y;
};

template <class Type>
struct XYZ
{
	XYZ() {}
	XYZ(const Type& x, const Type& y, const Type& z)
		:x(x), y(y), z(z) {}

	static constexpr int sz() { return 3; }

	Type x, y, z;
};

template <class Type>
struct XYZW
{
	XYZW() {}
	XYZW(const Type& x, const Type& y, const Type& z, const Type& w)
		:x(x), y(y), z(z), w(w) {}

	static constexpr int sz() { return 4; }

	Type x, y, z, w;
};

////////////////////////////////////////////////////////////////////////////////
//                                    Vector
////////////////////////////////////////////////////////////////////////////////

//template <class Type, int Size>
template <template<class> class Base, class Type>
struct vec : public Base<Type>
{
	// Constructors
	vec() {}
	vec(const Type& val)
	{
		// On Visual Studio 2022 unrolling the fill is slightly faster.
		//this->fill(val);
		if constexpr (this->sz() > 0) { (*this)[0] = val; }
		if constexpr (this->sz() > 1) { (*this)[1] = val; }
		if constexpr (this->sz() > 2) { (*this)[2] = val; }
		if constexpr (this->sz() > 3) { (*this)[3] = val; }
	}

	vec(const Type& x, const Type& y) : Base<Type>(x, y) {}
	vec(const Type& x, const Type& y, const Type& z) : Base<Type>(x, y, z) {}
	vec(const Type& x, const Type& y, const Type& z, const Type& w) : Base<Type>(x, y, z, w) {}

	/*explicit vec(const vec& v)
	{
		std::copy(v.begin(), v.end(), this->begin());
	}*/

	// Using an std::initializer_list to initialise an underlying array is
	// suprisingly difficult: https://stackoverflow.com/q/5549524
	// I did not get the proposed variadic template solution to work (perhaps
	// because array is a base, not a member?), but the simple copy() works.
	//vec(std::initializer_list<Type> l) : Base<Type>(l)
	//{
		// On Visual Studio 2022 unrolling the copy is slightly faster.
		//std::copy(l.begin(), l.end(), this->begin());
		//if constexpr (Base<Type>::sz() > 0) { (*this)[0] = l.begin()[0]; }
		//if constexpr (Base<Type>::sz() > 1) { (*this)[1] = l.begin()[1]; }
		//if constexpr (Base<Type>::sz() > 2) { (*this)[2] = l.begin()[2]; }
		//if constexpr (Base<Type>::sz() > 3) { (*this)[3] = l.begin()[3]; }
	//}

	//vec(std::array<Type, Size> a) : std::array<Type, Size>(a) {}

	// For static_cast support
	template <typename CastType>
	explicit operator vec<Base, CastType>() const
	{
		vec<Base, CastType> result;
		for (uint32_t ct = 0; ct < Base<Type>::sz(); ++ct) { result[ct] = static_cast<CastType>((*this)[ct]); }
		return result;
	}

	/*template <typename SrcType>
	explicit vec(const vec<Base, SrcType>& v)
	{
		if constexpr (v.sz() > 0) { this->x = static_cast<SrcType>(v.x); }
		if constexpr (v.sz() > 1) { this->y = static_cast<SrcType>(v.y); }
		if constexpr (v.sz() > 2) { this->z = static_cast<SrcType>(v.z); }
		if constexpr (v.sz() > 3) { this->w = static_cast<SrcType>(v.w); }
	}*/

	// Named element access
	//const Type& x() const { static_assert(Base<Type>::sz() > 0); return (*this)[0]; }
	//const Type& y() const { static_assert(Base<Type>::sz() > 1); return (*this)[1]; }
	//const Type& z() const { static_assert(Base<Type>::sz() > 2); return (*this)[2]; }
	//const Type& w() const { static_assert(Base<Type>::sz() > 3); return (*this)[3]; }

	/**************************** Support functions ***************************/

	// Apply function in-place to each element of this vector. The explicit
	// unrolling seems to help performance, at least on Visual C++ 2022.
	template<class UnaryFunc>
	constexpr vec& apply(UnaryFunc func)
	{
		static_assert(this->sz() <= 4);
		if constexpr (this->sz() > 0) { func((*this)[0]); }
		if constexpr (this->sz() > 1) { func((*this)[1]); }
		if constexpr (this->sz() > 2) { func((*this)[2]); }
		if constexpr (this->sz() > 3) { func((*this)[3]); }
		return *this;
	}

	// As above, but using the corresponding element of the params vector.
	template<class BinaryFunc>
	constexpr vec& apply(const vec& params, BinaryFunc func)
	{
		static_assert(this->sz() <= 4);
		if constexpr (this->sz() > 0) { func((*this)[0], params[0]); }
		if constexpr (this->sz() > 1) { func((*this)[1], params[1]); }
		if constexpr (this->sz() > 2) { func((*this)[2], params[2]); }
		if constexpr (this->sz() > 3) { func((*this)[3], params[3]); }
		return *this;
	}

	// Similar to 'apply()', but as a non-member making from supplied vector.
	template<class UnaryFunc>
	friend constexpr vec make_from(const vec& v, UnaryFunc func)
	{
		static_assert(v.sz() <= 4);
		vec<Base, Type> result;
		if constexpr (v.sz() > 0) { result[0] = func(v[0]); }
		if constexpr (v.sz() > 1) { result[1] = func(v[1]); }
		if constexpr (v.sz() > 2) { result[2] = func(v[2]); }
		if constexpr (v.sz() > 3) { result[3] = func(v[3]); }
		return result;
	}

	// As above, but taking two vectors as input.
	template<class BinaryFunc>
	friend constexpr auto make_from(const vec& v0,
		                            const vec& v1, BinaryFunc func)
	{
		static_assert(v0.sz() <= 4);
		using RetType = decltype(func(v0[0], v1[0])); // Can we avoid this?
		vec<Base, RetType> result;
		if constexpr (v0.sz() > 0) { result[0] = func(v0[0], v1[0]); }
		if constexpr (v0.sz() > 1) { result[1] = func(v0[1], v1[1]); }
		if constexpr (v0.sz() > 2) { result[2] = func(v0[2], v1[2]); }
		if constexpr (v0.sz() > 3) { result[3] = func(v0[3], v1[3]); }
		return result;
	}

	/************************** Overloaded operators **************************/

	// Overloaded operators following the guidelines here:
	// https://en.cppreference.com/w/cpp/language/operators
	// https://learn.microsoft.com/en-us/cpp/cpp/operator-overloading

	Type& operator[](int index)
	{
		//static_assert(sizeof(vec) == 3 * sizeof(Type)); // Packed
		//assert(index < size && "Index out of range");
		//return *((&Base<Type>::_x) + index);
		Type* ptr = &(Base<Type>::x);
		ptr = ptr + index;
		return *ptr;
		//return Base<Type>::operator[](index);
	}
	const Type& operator[](int index) const
	{
		//static_assert(sizeof(vec) == 3 * sizeof(Type)); // Packed
		//assert(index < size && "Index out of range");
		const Type* ptr = &(Base<Type>::x);
		ptr = ptr + index;
		return *ptr;
		//return Base<Type>::operator[](index);
	}

	friend bool operator==(const vec& lhs, const vec& rhs)
	{
		for (uint32_t ct = 0; ct < Base<Type>::sz(); ++ct)
		{
			if (lhs[ct] != rhs[ct]) return false;
		}
		return true;
	}

	friend bool operator!=(const vec& lhs, const vec& rhs) { return !(lhs == rhs); }

	friend bool operator<(vec const& lhs, vec const& rhs)
	{
		for (uint32_t ct = 0; ct < Base<Type>::sz() - 1; ++ct)
		{
			if (lhs[ct] < rhs[ct]) return true;
			if (rhs[ct] < lhs[ct]) return false;
		}

		return lhs[Base<Type>::sz() - 1] < rhs[Base<Type>::sz() - 1];
	}

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
		for (int i = 0; i < Base<Type>::sz(); ++i) { os << vector[i] << ","; }
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
		for (int i = 0; i < Base<Type>::size(); ++i) { if (!x[i]) return false; }
		return true;
	}

	friend bool any(const vec& x)
	{
		for (int i = 0; i < Base<Type>::size(); ++i) { if (x[i]) return true; }
		return false;
	}

	friend vec ceil(vec v)
	{
		return make_from(v, [](const Type& t) { return std::ceil(t); });
	}

	friend vec cross(const vec& a, const vec& b)
	{
		static_assert(a.sz() == 3, "Cross product only valid for 3D vectors");
		return { a.y * b.z - a.z * b.y,
				 a.z * b.x - a.x * b.z,
				 a.x * b.y - a.y * b.x };
	}

	friend Type dot(const vec& a, const vec& b)
	{
		Type result = 0;
		for (int i = 0; i < Base<Type>::sz(); ++i) { result += a[i] * b[i]; }
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

	friend int max_index(vec v)
	{
		int result = 0;
		for (int i = 0; i < Base<Type>::sz(); ++i)
		{ 
			if (v[i] > v[result]) { result = i; }
		}
		return result;
	}

	friend Type max_value(vec v)
	{
		Type result = v[0];
		for (int i = 0; i < Base<Type>::sz(); ++i)
		{
			if (v[i] > result) { result = v[i]; }
		}
		return result;
	}

	friend vec min(vec v0, const vec& v1)
	{
		return make_from(v0, v1,
			[](const Type& x, const Type& y) { return std::min(x, y); });
	}

	friend int min_index(vec v)
	{
		int result = 0;
		for (int i = 0; i < Base<Type>::sz(); ++i)
		{
			if (v[i] < v[result]) { result = i; }
		}
		return result;
	}

	friend Type min_value(vec v)
	{
		Type result = v[0];
		for (int i = 0; i < Base<Type>::sz(); ++i)
		{
			if (v[i] < result) { result = v[i]; }
		}
		return result;
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

	friend vec<Base, long int> round_to_int(const vec& v)
	{
		// Rounding in C++ is suprisingly complex (e.g. lround() vs lrint()) and 
		// built-in functions can be slow (https://stackoverflow.com/q/53962727).
		// Hence we use a simpler method here.
		return static_cast<vec<Base, long int>>(floor(v + vec(0.5)));
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
	friend vec<Base, bool> equal(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l == r; });
	}

	friend vec<Base, bool> lessThan(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l < r; });
	}

	friend vec<Base, bool> lessThanEqual(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l <= r; });
	}

	friend vec<Base, bool> greaterThan(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l > r; });
	}

	friend vec<Base, bool> greaterThanEqual(vec x, vec y)
	{
		return make_from(
			x, y, [](const Type& l, const Type& r) {return l >= r; });
	}
};

template <class Type>
using Vec2 = vec<XY, Type>;

template <class Type>
using Vec3 = vec<XYZ, Type>;

template <class Type>
using Vec4 = vec<XYZW, Type>;

// Typedefs for basic vector types
using vec2i = Vec2<int32>;
using vec2u = Vec2<uint32>;
using vec2b = Vec2<bool>;
using vec2f = Vec2<float>;
using vec2d = Vec2<double>;

using vec3i = Vec3<int32>;
using vec3u = Vec3<uint32>;
using vec3b = Vec3<bool>;
using vec3f = Vec3<float>;
using vec3d = Vec3<double>;

using vec4i = Vec4<int32>;
using vec4u = Vec4<uint32>;
using vec4b = Vec4<bool>;
using vec4f = Vec4<float>;
using vec4d = Vec4<double>;

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
		for (int y = 0; y < 4; y++)
		{
			for (int x = 0; x < 4; x++)
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
		for (uint32_t ct = 0; ct < 4; ++ct) { data[ct] = static_cast<Vec4<Type>>(matrix.data[ct]); }
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

template <class Type>
class Ray3
{
public:
	Ray3() {}

	Ray3(const Vec3<Type>& origin, const Vec3<Type>& dir)
		: mOrigin(origin)
		, mDir(dir) {
	}

	template <typename CastType> explicit Ray3(const Ray3<CastType>& ray)
	{
		//for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
		mOrigin = static_cast<Vec3<Type>>(ray.mOrigin);
		mDir = static_cast<Vec3<Type>>(ray.mDir);
	}

public:
	Vec3<Type> mOrigin;
	Vec3<Type> mDir; // FIXME - Would float always be sufficient for the direction?
};

typedef Ray3<float> Ray3f;
typedef Ray3<double> Ray3d;

template <class Type>
class Box
{
public:

	// Default constructed box starts off as invalid (min bigger than max))
	Box() { invalidate(); }

	Box(const Vec3<Type>& lower, const Vec3<Type>& upper)
		:mExtents{ lower, upper } {
	}

	// For casting betwween Box types of matching size.
	template <typename CastType> explicit Box(const Box<CastType>& box)
		:mExtents{
			static_cast<Vec3<Type>>(box.lower()),
			static_cast<Vec3<Type>>(box.upper())
		} {
	}

	const Vec3<Type>& lower() const { return mExtents[0]; }
	const Vec3<Type>& upper() const { return mExtents[1]; }

	// Should try to remove these non-const versions...
	Vec3<Type>& lower() { return mExtents[0]; }
	Vec3<Type>& upper() { return mExtents[1]; }

	void accumulate(const Vec3<Type>& value)
	{
		lower() = min(lower(), value);
		upper() = max(upper(), value);
	}

	// Note: Could templatise on container
	void accumulate(const std::array<Vec3<Type>, 3>& points)
	{
		for (auto const& point : points) {
			accumulate(point);
		}
	}

	void accumulate(const Box<Type>& other)
	{
		lower() = min(lower(), other.lower());
		upper() = max(upper(), other.upper());
	}

	bool contains(const Vec3<Type>& value) const
	{
		return (value.x >= lower().x) && (value.y >= lower().y) && (value.z >= lower().z) &&
			(value.x <= upper().x) && (value.y <= upper().y) && (value.z <= upper().z);
	}

	bool contains(const Box<Type>& other) const
	{
		return contains(other.lower()) && contains(other.upper());
	}

	void dilate(Type amount)
	{
		Vec3<Type> amountAsVec = { amount, amount, amount };
		lower() -= amountAsVec;
		upper() += amountAsVec;
	}

	void invalidate()
	{
		lower() = Vec3<Type>(std::numeric_limits<Type>::max());
		upper() = Vec3<Type>(std::numeric_limits<Type>::lowest());
	}

	bool isValid()
	{
		return lower().x <= upper().x &&
			lower().y <= upper().y &&
			lower().z <= upper().z;
	}

	static Box<Type> invalid()
	{
		return Box<Type>(); // Default-constructed box is already invalid.
	}

public:
	Vec3<Type> mExtents[2];
};

template <class Type>
bool operator==(const Box<Type>& lhs, const Box<Type>& rhs)
{
	return lhs.lower() == rhs.lower() && lhs.upper() == rhs.upper();
}

template <class Type>
bool overlaps(const Box<Type>& a, const Box<Type>& b)
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

template <class Type>
class Box3 : public Box<Type>
{
public:
	Box3() : Box<Type>() {}
	Box3(const Vec3<Type>& lower, const Vec3<Type>& upper) : Box<Type>(lower, upper) {}
	explicit Box3(const Box<int32>& box) : Box<Type>(box) {}

	vec3f centre() const { return (this->lower() + this->upper()) * 0.5f; }

	Type volume() const
	{
		// FIXME - Handle negative volumes?
		Vec3<Type> dims = this->upper() - this->lower();
		return dims.x * dims.y * dims.z;
	}
};

typedef Box3<float> Box3f;
typedef Box3<double> Box3d;

class Box3i : public Box<int32>
{
public:
	Box3i() : Box<int32>() {}
	Box3i(const Vec3<int32>& lower, const Vec3<int32>& upper) : Box<int32>(lower, upper) {}
	explicit Box3i(const Box<float>& box) : Box<int32>(box) {}

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
		, randZ(bounds.lower().z, bounds.upper().z) {
	}

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
		, randZ(bounds.lower().z, bounds.upper().z) {
	}

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
		bool operator!=(const iterator& other) const { return mIndex != other.mIndex; }
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
		: vertices{ { vertex0, vertex1, vertex2 } } {
	}

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
RayBoxIntersection intersect(const Ray3<Type>& ray, const Box3<Type>& box)
{
	// Inverse direction could be precomputed and stored in the ray
	// if we find we often intersect the same ray with multiple boxes.
	const Vec3<Type> invDir = Vec3<Type>({ 1.0f, 1.0f, 1.0f }) / ray.mDir;

	const Vec3<Type> lower = (box.lower() - ray.mOrigin) * invDir;
	const Vec3<Type> upper = (box.upper() - ray.mOrigin) * invDir;

	const Vec3<Type> minCorner = min(lower, upper);
	const Vec3<Type> maxCorner = max(lower, upper);

	RayBoxIntersection intersection;
	intersection.entry = max_value(minCorner);
	intersection.exit = min_value(maxCorner);
	return intersection;
}

bool intersect(const Ray3f& ray, const Triangle& triangle, float& t);

}

#endif // CUBIQUITY_GEOMETRY_H
