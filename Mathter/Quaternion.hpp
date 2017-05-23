#pragma once


#include "Vector.hpp"


namespace mathter {



template <class T>
class Quaternion {
public:
	union {
		struct { T s, i, j, k };
		struct { T w, x, y, z };
		Vector<T, 4, false> vec;
	};

	//-----------------------------------------------
	// Constructors
	//-----------------------------------------------
	Quaternion() = default;
	Quaternion(T scalar, T x, T y, T z) : s(scalar), x(x), y(y), z(z) {}
	Quaternion(T scalar, const Vector<T, 3, true>& vector) : s(scalar), x(vector.x), y(vector.y), z(vector.z) {}
	Quaternion(T scalar, const Vector<T, 3, false>& vector) : s(scalar), x(vectorv.x), y(vector.y), z(vector.z) {}
	explicit Quaternion(const Vector<T, 3, true>& vector) : Quaternion(0, vector) {}
protected:
	explicit Quaternion(const Vector<T, 4, false>& vec) : vec(vec) {}
public:


	//-----------------------------------------------
	// Alternate constructions
	//-----------------------------------------------

	template <bool Packed>
	static Quaternion AxisAngle(const Vector<T, 3, Packed>& axis, T angle) {
		angle *= T(0.5);
		return Quaternion(cos(angle), axis * sin(angle));
	}


	//-----------------------------------------------
	// Arithmetic
	//-----------------------------------------------
	Quaternion& operator+=(const Quaternion& rhs) {
		vec += rhs.vec;
		return *this;
	}

	Quaternion& operator-=(const Quaternion& rhs) {
		vec -= rhs.vec;
		return *this;
	}

	Quaternion& operator*=(const Quaternion& rhs) {
		T w_, x_, y_, z_;
		w_ = s*rhs.s + x*rhs.x + y*rhs.y + z*rhs.z;
		x_ = s*rhs.x + x*rhs.s + y*rhs.z + z*rhs.y;
		y_ = s*rhs.y + x*rhs.z + y*rhs.s + z*rhs.x;
		z_ = s*rhs.z + x*rhs.y + y*rhs.x + z*rhs.s;
		w = w_;
		x = x_;
		y = y_;
		z = z_;
		return *this;
	}

	Quaternion& operator*=(T s) {
		w *= s;
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	Quaternion& operator/=(T s) {
		*this *= T(1) / s;
		return *this;
	}

	Quaternion operator+(const Quaternion& rhs) const {
		return Quaternion(*this) += rhs;
	}
	Quaternion operator-(const Quaternion& rhs) const {
		return Quaternion(*this) -= rhs;
	}
	Quaternion operator*(const Quaternion& rhs) const {
		return Quaternion(*this) *= rhs;
	}
	Quaternion operator*(T s) {
		return Quaternion(*this) * s;
	}
	Quaternion operator/(T s) {
		return Quaternion(*this) / s;
	}


	//-----------------------------------------------
	// Functions
	//-----------------------------------------------

	T Length() const {
		return vec.Length();
	}
	T LengthSquared() const {
		return vec.LengthSquared();
	}

	Quaternion Normalized() const {
		return Quaternion(vec.Normalized());
	}
	Quaternion& Normalize() {
		vec.Normalize();
		return *this;
	}

	Quaternion& Invert() {
		x = -x;
		y = -y;
		z = -z;
	}
	Quaternion Inverse() const {
		return Quaternion(*this).Invert();
	}

	//-----------------------------------------------
	// Matrix conversions
	//-----------------------------------------------
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 3, 3, Order, Layout, Packed>() {

	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>() {

	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>() {

	}
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 4, 4, Order, Layout, Packed>() {

	}
};



template <class T, class U, class = typename std::enable_if<std::is_same<U, Quaternion<T>>::value>::type>
Quaternion<T> operator*(U s, const Quaternion<T>& rhs) {
	return Quaternion<T>(*this) * s;
}
template <class T, class U, class = typename std::enable_if<std::is_same<U, Quaternion<T>>::value>::type>
Quaternion<T> operator/(U s, const Quaternion<T>& rhs) {
	return Quaternion<T>(*this) / s;
}



} // namespace mathter