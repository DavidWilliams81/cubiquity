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
struct Vec2
{
	typedef Type value_type;
	static constexpr bool is_vec() { return true; }
	static constexpr int size() { return 2; }

	// Constructors
	Vec2() {}
	Vec2(const Type& val) : x(val), y(val) {}
	Vec2(const Type& x, const Type& y) : x(x), y(y) {}

	// Converting constructor
	template <typename SrcType>
	explicit Vec2(const Vec2<SrcType>& v) :x(v.x), y(v.y) {}

	Type& operator[](int index)
	{
		assert(index < size() && "Index out of range");
		if (index == 0) return x;
		if (index == 1) return y;
	}
	const Type& operator[](int index) const
	{
		assert(index < size() && "Index out of range");
		if (index == 0) return x;
		if (index == 1) return y;
	}

	Type x, y;
};

template <class Type>
struct Vec3
{
	typedef Type value_type;
	static constexpr bool is_vec() { return true; }
	static constexpr int size() { return 3; }

	// Constructors
	Vec3() {}
	Vec3(const Type& val): x(val), y(val), z(val) {}
	Vec3(const Type& x, const Type& y, const Type& z) :x(x), y(y), z(z) {}

	// Converting constructor
	template <typename SrcType>
	explicit Vec3(const Vec3<SrcType>& v) :x(v.x), y(v.y), z(v.z) {}

	Type& operator[](int index)
	{
		static_assert(sizeof(Vec3<Type>) == 3 * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}
	const Type& operator[](int index) const
	{
		static_assert(sizeof(Vec3<Type>) == 3 * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}

	Type x, y, z;
};

template <class Type>
bool operator==(const Vec3<Type>& lhs, const Vec3<Type>& rhs)
{
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}
template <class Type>
bool operator!=(const Vec3<Type>& lhs, const Vec3<Type>& rhs) { return !(lhs == rhs); }

template <class Type>
bool operator<(Vec3<Type> const& lhs, Vec3<Type> const& rhs)
{
	if (lhs.x < rhs.x) return true;
	if (rhs.x < lhs.x) return false;

	if (lhs.y < rhs.y) return true;
	if (rhs.y < lhs.y) return false;

	return lhs.z < rhs.z;
}

template <class Type>
struct Vec4
{
	typedef Type value_type;
	static constexpr bool is_vec() { return true; }
	static constexpr int size() { return 4; }

	// Constructors
	Vec4() {}
	Vec4(const Type& val) : x(val), y(val), z(val), w(val) {}
	Vec4(const Type& x, const Type& y, const Type& z, const Type& w)
		:x(x), y(y), z(z), w(w) {}
	

	//vec(std::array<Type, Size> a) : std::array<Type, Size>(a) {}

	// For static_cast support
	template <typename CastType>
	explicit operator Vec4<CastType>() const
	{
		Vec4<CastType> result;
		for (uint32_t ct = 0; ct < 4; ++ct) { result[ct] = static_cast<CastType>((*this)[ct]); }
		return result;
	}

	Type& operator[](int index)
	{
		static_assert(sizeof(Vec4<Type>) == 4 * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}
	const Type& operator[](int index) const
	{
		static_assert(sizeof(Vec4<Type>) == 4 * sizeof(Type)); // Packed
		assert(index < size() && "Index out of range");
		return *((&x) + index);
	}

	Type x, y, z, w;
};

/****************************** Support functions *****************************/

// Apply function in-place to each element of this vector. The explicit
// unrolling seems to help performance, at least on Visual C++ 2022.
template<class VecType, class UnaryFunc>
constexpr VecType& apply(VecType& v, UnaryFunc func)
{
	if constexpr (v.size() > 0) func(v.x);
	if constexpr (v.size() > 1) func(v.y);
	if constexpr (v.size() > 2) func(v.z);
	if constexpr (v.size() > 3) func(v.w);
	return v;
}

// As above, but using the corresponding element of the params vector.
template<class VecType, class BinaryFunc>
constexpr VecType& apply(VecType& v, const VecType& params, BinaryFunc func)
{
	if constexpr (v.size() > 0) func(v.x, params.x);
	if constexpr (v.size() > 1) func(v.y, params.y);
	if constexpr (v.size() > 2) func(v.z, params.z);
	if constexpr (v.size() > 3) func(v.w, params.w);
	return v;
}

// Similar to 'apply()', but as a non-member making from supplied vector.
template <template<class> class Vec, class T, class UnaryFunc>
constexpr Vec<T> make_from(const Vec<T>& v, UnaryFunc func)
{
	Vec<T> result;
	if constexpr (v.size() > 0) result.x = func(v.x);
	if constexpr (v.size() > 1) result.y = func(v.y);
	if constexpr (v.size() > 2) result.z = func(v.z);
	if constexpr (v.size() > 3) result.w = func(v.w);
	return result;
}

// As above, but taking two vectors as input.
//template<class BinaryFunc>
template <template<class> class Vec, class T, class BinaryFunc>
constexpr auto make_from(const Vec<T>& v0,
	const Vec<T>& v1, BinaryFunc func)
{
	using RetType = decltype(func(v0.x, v1.x)); // Can we avoid this?
	Vec<RetType> result;
	if constexpr (v0.size() > 0) result.x = func(v0.x, v1.x);
	if constexpr (v0.size() > 1) result.y = func(v0.y, v1.y);
	if constexpr (v0.size() > 2) result.z = func(v0.z, v1.z);
	if constexpr (v0.size() > 3) result.w = func(v0.w, v1.w);
	return result;
}

// Stream insertion operator.
template <class VecType, typename = std::enable_if<VecType::is_vec()>>
std::ostream& operator<<(std::ostream& os, const VecType& vector) {
	os << "[";
	for (int i = 0; i < vector.size(); i++) { os << vector[i] << ","; }
	os << "\b]"; // Backspace to remove last comma.
	return os;
}

/**************************** Overloaded operators ****************************/

// Overloaded operators following the guidelines here:
// https://en.cppreference.com/w/cpp/language/operators
// https://learn.microsoft.com/en-us/cpp/cpp/operator-overloading

// Unary operators
template <template<class> class Vec, class T>
Vec<T> operator-(Vec<T> v) { return apply(v, [](T& t) {t = -t; }); }
template <template<class> class Vec, class T>
Vec<T> operator~(Vec<T> v) { return apply(v, [](T& t) {t = ~t; }); }

// Binary operators with assignment and scalar on right-hand side
template <template<class> class Vec, class T>
Vec<T>& operator+= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t += rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator-= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t -= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator*= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t *= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator/= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t /= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator%= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t %= rhs; }); }

template <template<class> class Vec, class T>
Vec<T>& operator&= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t &= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator|= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t |= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator^= (Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t ^= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator<<=(Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t <<= rhs; }); }
template <template<class> class Vec, class T>
Vec<T>& operator>>=(Vec<T>& lhs, const T& rhs) { return apply(lhs, [rhs](T& t) {t >>= rhs; }); }


// Binary operators with assignment and vector on right-hand side
template <template<class> class Vec, class T>
Vec<T>& operator+= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l += r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator-= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l -= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator*= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l *= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator/= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l /= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator%= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l %= r; }); }

template <template<class> class Vec, class T>
Vec<T>& operator&= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l &= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator|= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l |= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator^= (Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l ^= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator<<=(Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l <<= r; }); }
template <template<class> class Vec, class T>
Vec<T>& operator>>=(Vec<T>& lhs, const Vec<T>& rhs) { return apply(lhs, rhs, [](T& l, const T& r) {l >>= r; }); }


// Binary operators with scalar on right-hand side
template <template<class> class Vec, class T>
Vec<T> operator+ (Vec<T> lhs, const T& rhs) { return lhs += rhs; }
template <template<class> class Vec, class T>
Vec<T> operator- (Vec<T> lhs, const T& rhs) { return lhs -= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator* (Vec<T> lhs, const T& rhs) { return lhs *= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator/ (Vec<T> lhs, const T& rhs) { return lhs /= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator% (Vec<T> lhs, const T& rhs) { return lhs %= rhs; }

template <template<class> class Vec, class T>
Vec<T> operator& (Vec<T> lhs, const T& rhs) { return lhs &= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator| (Vec<T> lhs, const T& rhs) { return lhs |= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator^ (Vec<T> lhs, const T& rhs) { return lhs ^= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator<<(Vec<T> lhs, const T& rhs) { return lhs <<= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator>>(Vec<T> lhs, const T& rhs) { return lhs >>= rhs; }

// Binary operators with vector on right-hand side
template <template<class> class Vec, class T>
Vec<T> operator+ (Vec<T> lhs, const Vec<T>& rhs) { return lhs += rhs; }
// FIXME I don't understand why this enable_if is needed. Without it I get
// non-sensical errors in voxelization when working with std::span. These errors
// only occur on Vixual Studioo - GCC compiles fine. I'm hoping the issue gets
// fixed as I review the rest of the code.
template <template<class> class Vec, class T, typename = std::enable_if<Vec<T>::is_vec()>>
Vec<T> operator- (Vec<T> lhs, const Vec<T>& rhs) { return lhs -= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator* (Vec<T> lhs, const Vec<T>& rhs) { return lhs *= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator/ (Vec<T> lhs, const Vec<T>& rhs) { return lhs /= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator% (Vec<T> lhs, const Vec<T>& rhs) { return lhs %= rhs; }

template <template<class> class Vec, class T>
Vec<T> operator& (Vec<T> lhs, const Vec<T>& rhs) { return lhs &= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator| (Vec<T> lhs, const Vec<T>& rhs) { return lhs |= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator^ (Vec<T> lhs, const Vec<T>& rhs) { return lhs ^= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator<<(Vec<T> lhs, const Vec<T>& rhs) { return lhs <<= rhs; }
template <template<class> class Vec, class T>
Vec<T> operator>>(Vec<T> lhs, const Vec<T>& rhs) { return lhs >>= rhs; }


/*************************** Other vector operations **************************/

template <class VecType>
VecType abs(VecType v)
{
	return make_from(v, [](const VecType::value_type& t) { return std::abs(t); });
}

template <class VecType>
bool all(const VecType& x)
{
	for (int i = 0; i < 3; ++i) { if (!x[i]) return false; }
	return true;
}

template <class VecType>
bool any(const VecType& x)
{
	for (int i = 0; i < 3; ++i) { if (x[i]) return true; }
	return false;
}

template <class VecType>
VecType ceil(VecType v)
{
	return make_from(v, [](const VecType::value_type& t) { return std::ceil(t); });
}

template <class Type>
Vec3<Type> cross(const Vec3<Type>& a, const Vec3<Type>& b)
{
	return { a.y * b.z - a.z * b.y,
			 a.z * b.x - a.x * b.z,
			 a.x * b.y - a.y * b.x };
}

template <class VecType>
VecType::value_type dot(const VecType& a, const VecType& b)
{
	typename VecType::value_type result = 0;
	for (int i = 0; i < 3; ++i) { result += a[i] * b[i]; }
	return result;
}

template <class VecType>
VecType floor(VecType v)
{
	return make_from(v, [](const VecType::value_type& t) { return std::floor(t); });
}

template <class VecType>
VecType fract(VecType v)
{
	return v - floor(v);
}

template <class VecType>
VecType::value_type length(const VecType& v)
{
	// Floating point componenents only due to sqrt(). Integer vectors can
	// be cast to floating point vectors prior to calling this function.
	static_assert(std::is_floating_point<typename VecType::value_type>::value);
	return sqrt(dot(v, v));
}

template <class VecType>
VecType max(VecType v0, const VecType& v1)
{
	return make_from(v0, v1,
		[](const VecType::value_type& x, const VecType::value_type& y) { return std::max(x, y); });
}

template <class VecType>
int max_index(VecType v)
{
	if (v.x >= v.y && v.x >= v.z) return 0;
	if (v.y >= v.x && v.y >= v.z) return 1;
	return 2;
}

template <class VecType>
VecType::value_type max_value(VecType v)
{
	return std::max(std::max(v.x, v.y), v.z);
}

template <class VecType>
VecType min(VecType v0, const VecType& v1)
{
	return make_from(v0, v1,
		[](const VecType::value_type& x, const VecType::value_type& y) { return std::min(x, y); });
}

template <class VecType>
VecType::value_type min_value(VecType v)
{
	return std::min(std::min(v.x, v.y), v.z);
}

template <class VecType>
VecType mix(const VecType& x, const VecType& y, const VecType& a)
{
	return x * (Vec3(1) - a) + y * a;
}

template <class VecType>
VecType normalize(const VecType& v)
{
	static_assert(std::is_floating_point<typename VecType::value_type>::value);
	assert(length(v) >= 0.001f);
	return v / length(v);
}

template <class VecType>
VecType pow(VecType x, const typename VecType::value_type& y)
{
	return make_from(x, [y](const VecType::value_type& t) { return std::pow(t, y); });
}

template <class Type>
Vec2<long int> round_to_int(const Vec2<Type>& v)
{
	// Rounding in C++ is suprisingly complex (e.g. lround() vs lrint()) and 
	// built-in functions can be slow (https://stackoverflow.com/q/53962727).
	// Hence we use a simpler method here.
	return static_cast<Vec2<long int>>(floor(v + Vec2(0.5)));
}

template <class VecType>
VecType sign(VecType v)
{
	return apply(v, [](VecType::value_type& t) {t = std::copysign(1.0f, t); });
}

template <class VecType>
VecType step(float edge, VecType v)
{
	return apply(v, [edge](VecType::value_type& t) {t = t < edge ? 0 : 1; });
}

// GLSL-compatible comparisons
template <class VecType>
auto equal(VecType x, VecType y)
{
	return make_from(
		x, y, [](const VecType::value_type& l, const VecType::value_type& r) {return l == r; });
}

template <class VecType>
auto lessThan(VecType x, VecType y)
{
	return make_from(
		x, y, [](const VecType::value_type& l, const VecType::value_type& r) {return l < r; });
}

template <class VecType>
auto lessThanEqual(VecType x, VecType y)
{
	return make_from(
		x, y, [](const VecType::value_type& l, const VecType::value_type& r) {return l <= r; });
}

template <class VecType>
auto greaterThan(VecType x, VecType y)
{
	return make_from(
		x, y, [](const VecType::value_type& l, const VecType::value_type& r) {return l > r; });
}

template <class VecType>
auto greaterThanEqual(VecType x, VecType y)
{
	return make_from(
		x, y, [](const VecType::value_type& l, const VecType::value_type& r) {return l >= r; });
}

// Typedefs for basic vector types
typedef Vec2<int32> vec2i;
typedef Vec2<uint32> vec2u;
typedef Vec2<bool> vec2b;
typedef Vec2<float> vec2f;
typedef Vec2<double> vec2d;

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

template <class Type>
class Ray3
{
public:
	Ray3() {}

	Ray3(const Vec3<Type>& origin, const Vec3<Type>& dir)
		: mOrigin(origin)
		, mDir(dir)	{}

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
		:mExtents{ lower, upper } {}

	// For casting betwween Box types of matching size.
	template <typename CastType> explicit Box(const Box<CastType>& box)
		:mExtents { 
			static_cast<Vec3<Type>>(box.lower()),
			static_cast<Vec3<Type>>(box.upper())
		} {}

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
