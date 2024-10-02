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
#include <unordered_map>
#include <span>
#include <vector>
#include <random>
#include <type_traits>

namespace Cubiquity
{
	const float Pi = 3.14159265358979f;

	// See https ://stackoverflow.com/a/4609795
	template <typename T> int sign(T val) {
		return (T(0) < val) - (val < T(0));
	}

	////////////////////////////////////////////////////////////////////////////////
	//									Vector
	////////////////////////////////////////////////////////////////////////////////

	// Vector is a simple C++ aggregate. Could publicly inherit from std::array to
	// eliminate some boilerplate code (see https://stackoverflow.com/a/24281360),
	// but this only seems to works on  GCC 11 onwards (not Clang or MSVC).
	template <class Type, int Size>
	struct Vector
	{
		// Initialise array to single value (aggregates don't have constructors).
		static constexpr Vector<Type, Size> filled(const Type& value) {
			Vector<Type, Size> ret = { 0, 0, 0 };
			ret.fill(value);
			return ret;
		}

		// For static_cast support
		template <typename CastType>
		explicit operator Vector<CastType, Size>() const {
			Vector<CastType, Size> result;
			for (uint32_t ct = 0; ct < Size; ++ct) { result[ct] = static_cast<CastType>(mData[ct]); }
			return result;
		}

		// Nullary  arithmetic operators
		Vector<Type, Size> operator-() const {
			Vector<Type, Size> result; for (int i = 0; i < Size; ++i) { result[i] = -(*this)[i]; } return result; }
		Vector<Type, Size> operator~() const {
			Vector<Type, Size> result; for (int i = 0; i < Size; ++i) { result[i] = ~(*this)[i]; } return result;
		}

		// Unary arithmetic operators
		void operator+=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] += rhs; } }
		void operator+=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] += rhs[i]; } }
		void operator-=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] -= rhs; } }
		void operator-=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] -= rhs[i]; } }
		void operator*=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] *= rhs; } }
		void operator*=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] *= rhs[i]; } }
		void operator/=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] /= rhs; } }
		void operator/=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] /= rhs[i]; } }
		void operator%=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] %= rhs; } }
		void operator%=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] %= rhs[i]; } }

		void operator>>=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] >>= rhs; } }
		void operator>>=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] >>= rhs[i]; } }
		void operator<<=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] <<= rhs; } }
		void operator<<=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] <<= rhs[i]; } }

		void operator&=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] &= rhs; } }
		void operator&=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] &= rhs[i]; } }
		void operator|=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] |= rhs; } }
		void operator|=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] |= rhs[i]; } }
		void operator^=(Type const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] ^= rhs; } }
		void operator^=(Vector<Type, Size> const& rhs) { for (int i = 0; i < Size; ++i) { (*this)[i] ^= rhs[i]; } }

		// Binary arithmetic operators as friends (non-members)
		friend Vector<Type, Size> operator+(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result += rhs; return result; }
		friend Vector<Type, Size> operator+(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result += rhs; return result; }
		friend Vector<Type, Size> operator-(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result -= rhs; return result; }
		friend Vector<Type, Size> operator-(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result -= rhs; return result; }
		friend Vector<Type, Size> operator*(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result *= rhs; return result; }
		friend Vector<Type, Size> operator*(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result *= rhs; return result; }
		friend Vector<Type, Size> operator/(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result /= rhs; return result; }
		friend Vector<Type, Size> operator/(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result /= rhs; return result; }
		friend Vector<Type, Size> operator%(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result %= rhs; return result; }
		friend Vector<Type, Size> operator%(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result %= rhs; return result; }

		friend Vector<Type, Size> operator>>(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result >>= rhs; return result; }
		friend Vector<Type, Size> operator>>(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result >>= rhs; return result; }
		friend Vector<Type, Size> operator<<(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result <<= rhs; return result; }
		friend Vector<Type, Size> operator<<(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result <<= rhs; return result; }

		friend Vector<Type, Size> operator&(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result &= rhs; return result; }
		friend Vector<Type, Size> operator&(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result &= rhs; return result; }
		friend Vector<Type, Size> operator|(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result |= rhs; return result; }
		friend Vector<Type, Size> operator|(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result |= rhs; return result; }
		friend Vector<Type, Size> operator^(Vector<Type, Size> const& lhs, Type const& rhs) {
			Vector<Type, Size> result(lhs); result ^= rhs; return result; }
		friend Vector<Type, Size> operator^(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs) {
			Vector<Type, Size> result(lhs); result ^= rhs; return result; }

		// Comparison operators as friends (https://stackoverflow.com/a/1145635)
		friend bool operator==(const Vector<Type, Size>& lhs, const Vector<Type, Size>& rhs) { return lhs.mData == rhs.mData; }
		friend bool operator<(const Vector<Type, Size>& lhs, const Vector<Type, Size>& rhs) { return lhs.mData < rhs.mData; }

		// Stream insertion operator.
		friend std::ostream& operator<<(std::ostream& os, const Vector<Type, Size>& vector) {
			os << "[";
			for (auto v : vector) { os << v << ","; }
			os << "\b]"; // Backspace to remove last comma.
			return os;
		}

		// Indexed element access
		Type& operator[](int index) { assert(index < Size && "Index out of range"); return mData[index]; }
		const Type& operator[](int index) const { assert(index < Size && "Index out of range"); return mData[index]; }

		// Named element access
		const Type& x() const { static_assert(Size > 0, "Vector too small to call x()"); return (*this)[0]; }
		const Type& y() const { static_assert(Size > 1, "Vector too small to call y()"); return (*this)[1]; }
		const Type& z() const { static_assert(Size > 2, "Vector too small to call z()"); return (*this)[2]; }
		const Type& w() const { static_assert(Size > 3, "Vector too small to call w()"); return (*this)[3]; }

		// FIXME - Would be nice if we could avoid the copy here, and just return the curent vector cast to the correct type.
		Vector<Type, 3> xyz() const { static_assert(Size > 2, "Vector too small to call xyz()"); return Vector<Type, 3> { (*this)[0], (*this)[1], (*this)[2] }; }

		// Forwards for std::array interface.
		auto begin() noexcept { return mData.begin(); }
		auto begin() const noexcept { return mData.begin(); }
		auto cbegin() const noexcept { return mData.cbegin(); }

		auto end() noexcept { return mData.end(); }
		auto end() const noexcept { return mData.end(); }
		auto cend() const noexcept { return mData.cend(); }

		void fill(const Type& value) { mData.fill(value); }

		// Public (requirement for aggregate), but shouldn't be accessed directly.
		std::array<Type, Size> mData;
	};

	// Typedefs for commonly-used sizes
	template <typename Type> using Vector2 = Vector<Type, 4>;
	template <typename Type> using Vector3 = Vector<Type, 3>;
	template <typename Type> using Vector4 = Vector<Type, 4>;

	// Typedefs for vectors with fixed-size types
	typedef Vector<int32, 2> Vector2i32;
	typedef Vector<uint32, 2> Vector2u32;
	typedef Vector<float, 2> Vector2f;
	typedef Vector<double, 2> Vector2d;

	typedef Vector<int32, 3> Vector3i32;
	typedef Vector<uint32, 3> Vector3u32;
	typedef Vector<int64, 3> Vector3i64;
	typedef Vector<float, 3> Vector3f;
	typedef Vector<double, 3> Vector3d;
	typedef Vector<bool, 3> Vector3b;

	typedef Vector<int32, 4> Vector4i32;
	typedef Vector<uint32, 4> Vector4u32;
	typedef Vector<float, 4> Vector4f;
	typedef Vector<double, 4> Vector4d;
	typedef Vector<bool, 4> Vector4b;

	// Shorter typedefs where the user doesn't specify the size of the underlying type.
	// Internally we fix the size of the type to that which makes sense for the size of our
	// volume (rather than the machine native type) as this seems more useful for our application.
	typedef Vector2i32 Vector2i;
	typedef Vector2u32 Vector2u;

	typedef Vector3i32 Vector3i;
	typedef Vector3u32 Vector3u;

	typedef Vector4i32 Vector4i;
	typedef Vector4u32 Vector4u;

	// Applied in place so only works for operators which don't change type (e.g. not lround()).
	template<class Type, int Size, class UnaryOperation>
	void transform_in_place(Vector<Type, Size>& vec, UnaryOperation unary_op)
	{
		std::transform(vec.begin(), vec.end(), vec.begin(), unary_op);
	}

	// Component-wise minimum and maximum (output is a Vector rather than a scaler).
	template <class Type, int Size>
	Vector<Type, Size> min(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs)
	{
		Vector<Type, Size> result;
		for(int i = 0; i < Size; ++i) { result[i] = std::min(lhs[i], rhs[i]); }
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> max(Vector<Type, Size> const& lhs, Vector<Type, Size> const& rhs)
	{
		Vector<Type, Size> result;
		for(int i = 0; i < Size; ++i) { result[i] = std::max(lhs[i], rhs[i]); }
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> clamp(Vector<Type, Size> const& vec, Type minVal, Type maxVal)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::max(minVal, std::min(vec[i], maxVal));
		}
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> abs3(Vector<Type, Size> const& vec)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::abs(vec[i]); }
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> ceil(Vector<Type, Size> const& vec)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::ceil(vec[i]); }
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> floor(Vector<Type, Size> const& vec)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i) { result[i] = std::floor(vec[i]); }
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> fract(Vector<Type, Size> const& vec)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i) { result[i] = vec[i] - std::floor(vec[i]); }
		return result;
	}

	// Rounding in C++ is suprisingly complex (e.g. lround() vs lrint()) and built-in functions can
	// be very slow (see https://stackoverflow.com/q/53962727). Hence we use a simpler method here.
	template <class Type, int Size>
	Vector<long int, Size> round_to_int(Vector<Type, Size> const& vec)
	{
		Vector<long int, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::floor(vec[i] + Type(0.5));
		}
		return result;
	}

	// Common vector operations
	template <class Type, int Size>
	Type dot(const Vector<Type, Size>& a, const Vector<Type, Size>& b)
	{
		Type result(0);
		for(int i = 0; i < Size; ++i) { result += a[i] * b[i]; }
		return result;
	}

	// Should only be valid for 3D vectors. Also, what should the return type be?
	template <class Type, int Size>
	Vector<Type, Size> cross(const Vector<Type, Size>& a, const Vector<Type, Size>& b)
	{
		return { a.y()*b.z() - a.z()*b.y(), a.z()*b.x() - a.x()*b.z(), a.x()*b.y() - a.y()*b.x() };
	}

	template <class Type, int Size>
	Type length(const Vector<Type, Size>& vec)
	{
		static_assert(std::is_floating_point<Type>::value);
		return sqrt(dot(vec, vec));
	}

	template <class Type, int Size>
	Vector<Type, Size> normalize(const Vector<Type, Size>& vec)
	{
		Vector<Type, Size> result;
		float len = length(vec);

		// Could do better error handling here but it's hard to define a correct
		// behaviour. I think I prefer to just not pass in invalid vectors?
		//assert(len >= 0.001f);

		result = vec * (Type(1.0) / len);
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> mix(Vector<Type, Size> x, Vector<Type, Size> y, Vector<Type, Size> a)
	{
		Vector<Type, Size> ones = Vector<Type, Size>::filled(1);
		return x * (ones - a) + y * a;
	}

	template <class Type, int Size>
	Vector<Type, Size> sign(Vector<Type, Size> x)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::copysign(1.0f, x[i]);
		}
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> step(float edge, Vector<Type, Size> x)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] < edge ? 0 : 1;
		}
		return result;
	}

	template <class Type, int Size>
	Vector<Type, Size> pow(Vector<Type, Size> const& x, Vector<Type, Size> const& y)
	{
		Vector<Type, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = std::pow(x[i], y[i]);
		}
		return result;
	}

	template <class Type, int Size>
	bool all(Vector<Type, Size> x)
	{
		for (int i = 0; i < Size; ++i)
		{
			if (!x[i]) return false;
		}
		return true;
	}

	template <class Type, int Size>
	bool any(Vector<Type, Size> x)
	{
		for (int i = 0; i < Size; ++i)
		{
			if (x[i]) return true;
		}
		return false;
	}

	template <class Type, int Size>
	Vector<bool, Size> lessThan(Vector<Type, Size> x, Vector<Type, Size> y)
	{
		Vector<bool, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] < y[i];
		}
		return result;
	}

	template <class Type, int Size>
	Vector<bool, Size> lessThanEqual(Vector<Type, Size> x, Vector<Type, Size> y)
	{
		Vector<bool, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] <= y[i];
		}
		return result;
	}

	template <class Type, int Size>
	Vector<bool, Size> greaterThan(Vector<Type, Size> x, Vector<Type, Size> y)
	{
		Vector<bool, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] > y[i];
		}
		return result;
	}

	template <class Type, int Size>
	Vector<bool, Size> greaterThanEqual(Vector<Type, Size> x, Vector<Type, Size> y)
	{
		Vector<bool, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] >= y[i];
		}
		return result;
	}

	template <class Type, int Size>
	Vector<bool, Size> equal(Vector<Type, Size> x, Vector<Type, Size> y)
	{
		Vector<bool, Size> result;
		for (int i = 0; i < Size; ++i)
		{
			result[i] = x[i] == y[i];
		}
		return result;
	}

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

		Matrix4x4(const Vector4<Type>& v0, const Vector4<Type>& v1, const Vector4<Type>& v2, const Vector4<Type>& v3)
		{
			data[0] = v0;
			data[1] = v1;
			data[2] = v2;
			data[3] = v3;
		}

		// For casting betwween Vector types of matching size.
		template <typename CastType> explicit Matrix4x4(const Matrix4x4<CastType>& matrix)
		{
			for (uint32_t ct = 0; ct < 4; ++ct) { data[ct] = static_cast< Vector4<Type> >(matrix.data[ct]); }
		}

		Vector4<Type>& operator[](int index) { return data[index]; }

		const Vector4<Type>& operator[](int index) const { return data[index]; }

		void operator/=(Type const& rhs) { for (int i = 0; i < 4; ++i) { data[i] /= rhs; } }

	public:
		Vector4<Type> data[4];
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
	Matrix4x4<Type> lookAtRH(const Vector3<Type>& eye, const Vector3<Type>& center, const Vector3<Type>& up)
	{
		const Vector3<Type> f(normalize(center - eye));
		const Vector3<Type> s(normalize(cross(f, up)));
		const Vector3<Type> u(cross(s, f));

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
	Vector4<Type> mul(const Matrix4x4<Type>& a, const Vector4<Type>& b) { return a[0] * b.x() + a[1] * b.y() + a[2] * b.z() + a[3] * b.w(); }

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
	Matrix4x4<Type> translation_matrix(const Vector3<Type>& pos)
	{
		Matrix4x4<Type> result;
		result[3] = Vector4<Type>({ pos.x(), pos.y(), pos.z(), 1.0f });
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

		Ray(const Vector<Type, Size>& origin, const Vector<Type, Size>& dir)
			: mOrigin(origin)
			, mDir(dir)	{}

		template <typename CastType> explicit Ray(const Ray<CastType, Size>& ray)
		{
			//for (uint32_t ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
			mOrigin = static_cast<Vector<Type, Size>>(ray.mOrigin);
			mDir = static_cast<Vector<Type, Size>>(ray.mDir);
		}

	public:
		Vector<Type, Size> mOrigin;
		Vector<Type, Size> mDir; // FIXME - Would float always be sufficient for the direction?
	};

	typedef Ray<float, 2> Ray2f;
	typedef Ray<double, 2> Ray2d;
	typedef Ray<float, 3> Ray3f;
	typedef Ray<double, 3> Ray3d;

	template <class Type, int Size>
	class Box
	{
	public:

		// Default constructed box starts off as invalid (min bigger than max))
		Box() { invalidate(); }

		Box(const Vector<Type, Size>& lower, const Vector<Type, Size>& upper)
			:mExtents{ lower, upper } {}

		// For casting betwween Box types of matching size.
		template <typename CastType> explicit Box(const Box<CastType, Size>& box)
			:mExtents { 
				static_cast<Vector<Type, Size>>(box.lower()),
				static_cast<Vector<Type, Size>>(box.upper())
			} {}

		const Vector<Type, Size>& lower() const { return mExtents[0]; }
		const Vector<Type, Size>& upper() const { return mExtents[1]; }

		// Should try to remove these non-const versions...
		Vector<Type, Size>& lower() { return mExtents[0]; }
		Vector<Type, Size>& upper() { return mExtents[1]; }

		void accumulate(const Vector<Type, Size>& value)
		{
			lower() = Cubiquity::min(lower(), value);
			upper() = Cubiquity::max(upper(), value);
		}

		// Note: Could templatise on container
		void accumulate(const std::array<Vector<Type, Size>, 3>& points)
		{
			for (auto const& point : points) {
				accumulate(point);
			}
		}

		void accumulate(const Box<Type, Size>& other)
		{
			lower() = Cubiquity::min(lower(), other.lower());
			upper() = Cubiquity::max(upper(), other.upper());
		}

		bool contains(const Vector<Type, Size>& value) const
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
			Vector<Type, Size> amountAsVec = { amount, amount, amount };
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
		Vector<Type, Size> mExtents[2];
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
		Box3(const Vector<Type, 3>& lower, const Vector<Type, 3>& upper) : Box<Type, 3>(lower, upper) {}
		explicit Box3(const Box<int32, 3>& box) : Box<Type, 3>(box) {}

		Vector3f centre() const { return (this->lower() + this->upper()) * 0.5f; }

		Type volume() const
		{
			// FIXME - Handle negative volumes?
			Vector<Type, 3> dims = this->upper() - this->lower();
			return dims.x() * dims.y() * dims.z();
		}
	};

	typedef Box3<float> Box3f;
	typedef Box3<double> Box3d;

	class Box3i : public Box<int32, 3>
	{
	public:
		Box3i() : Box<int32, 3>() {}
		Box3i(const Vector<int32, 3>& lower, const Vector<int32, 3>& upper) : Box<int32, 3>(lower, upper) {}
		explicit Box3i(const Box<float, 3>& box) : Box<int32, 3>(box) {}

		static Box3i max()
		{
			Vector3i lower = Vector3i::filled(std::numeric_limits<int>::lowest());
			Vector3i upper = Vector3i::filled(std::numeric_limits<int>::max());
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

		Vector3f next() {
			Vector3f result = { randX(eng), randY(eng), randZ(eng) };
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

		Vector3i next() { return Vector3i({ randX(eng), randY(eng), randZ(eng) }); }

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
			iterator(uint32_t index, const Box3i& box) : eng(0), randX(box.lower().x(), box.upper().x()), randY(box.lower().y(), box.upper().y()), randZ(box.lower().z(), box.upper().z()), mIndex(index) { mValue = Vector3i({ randX(eng), randY(eng), randZ(eng) }); }
			iterator operator++() { ++mIndex; mValue = Vector3i({ randX(eng), randY(eng), randZ(eng) }); return *this; }
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

		float sideLength(int index) const;
		Vector3f computeNormal() const;
		float area() const;
		Cubiquity::Vector3f centre() const;

		std::array<Cubiquity::Vector3f, 3> vertices;
	};

	float distance(const Cubiquity::Vector3f& point, const Cubiquity::Triangle& triangle);

	/// A (conceptual) list of non-indexed triangles, actually stored as an std::vector
	typedef std::vector<Cubiquity::Triangle> TriangleList;
	typedef std::span<Triangle> TriangleSpan;
	typedef std::span<Triangle const> ConstTriangleSpan; // See https://stackoverflow.com/a/56895806

	Box3f computeBounds(const std::array<Cubiquity::Vector3f, 3>& points);
	Box3f computeBounds(ConstTriangleSpan triangles);
	void translate(TriangleList& triangles, const Cubiquity::Vector3f& dir);
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
		const Vector<Type, 3> invDir = Vector<Type, 3>({ 1.0f, 1.0f, 1.0f }) / ray.mDir;

		const Vector<Type, 3> lower = (box.lower() - ray.mOrigin) * invDir;
		const Vector<Type, 3> upper = (box.upper() - ray.mOrigin) * invDir;

		const Vector<Type, 3> minCorner = min(lower, upper);
		const Vector<Type, 3> maxCorner = max(lower, upper);

		RayBoxIntersection intersection;
		intersection.entry = *(std::max_element(minCorner.begin(), minCorner.end()));
		intersection.exit = *(std::min_element(maxCorner.begin(), maxCorner.end()));
		return intersection;
	}

	bool intersect(const Ray3f& ray, const Triangle& triangle, float& t);
	
}

#endif // CUBIQUITY_GEOMETRY_H
