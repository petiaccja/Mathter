//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma once


#include "Vector.hpp"


namespace mathter {



template <class T, bool Packed = false>
class Quaternion {
public:
	union {
		struct { T s, i, j, k; };
		struct { T w, x, y, z; };
		Vector<T, 4, Packed> vec;
	};

	//-----------------------------------------------
	// Constructors
	//-----------------------------------------------
	Quaternion() = default;
	Quaternion(const Quaternion& rhs) = default;
	Quaternion(T scalar, T x, T y, T z) : s(scalar), x(x), y(y), z(z) {}
	Quaternion(T scalar, const Vector<T, 3, true>& vector) : s(scalar), x(vector.x), y(vector.y), z(vector.z) {}
	Quaternion(T scalar, const Vector<T, 3, false>& vector) : s(scalar), x(vectorv.x), y(vector.y), z(vector.z) {}
	explicit Quaternion(const Vector<T, 3, true>& vector) : Quaternion(0, vector) {}
	Quaternion(const Quaternion<T, !Packed>& rhs) : vec(rhs.vec) {}

	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit Quaternion(const Matrix<U, 3, 3, Order, Layout, Packed>& rhs) {
		FromMatrix(rhs);
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit Quaternion(const Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>& rhs) {
		FromMatrix(rhs);
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit Quaternion(const Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>& rhs) {
		FromMatrix(rhs);
	}
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit Quaternion(const Matrix<U, 4, 4, Order, Layout, Packed>& rhs) {
		FromMatrix(rhs);
	}
protected:
	explicit Quaternion(const Vector<T, 4, false>& vec) : vec(vec) {}
public:
	//-----------------------------------------------
	// Assignment
	//-----------------------------------------------

	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	Quaternion& operator=(const Matrix<U, 3, 3, Order, Layout, Packed>& rhs) const {
		FromMatrix(rhs);
		return *this;
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	Quaternion& operator=(const Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>& rhs) const {
		FromMatrix(rhs);
		return *this;
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	Quaternion& operator=(const Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>& rhs) const {
		FromMatrix(rhs);
		return *this;
	}
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	Quaternion& operator=(const Matrix<U, 4, 4, Order, Layout, Packed>& rhs) const {
		FromMatrix(rhs);
		return *this;
	}



	//-----------------------------------------------
	// Alternate constructions
	//-----------------------------------------------

	template <bool VPacked>
	static Quaternion AxisAngle(const Vector<T, 3, VPacked>& axis, T angle) {
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
	// Comparison
	//-----------------------------------------------
	bool operator==(const Quaternion& rhs) const {
		return vec == rhs.vec;
	}
	bool operator!=(const Quaternion& rhs) const {
		return !(*this == rhs);
	}
	bool AlmostEqual(const Quaternion& rhs) const {
		return vec.AlmostEqual(rhs.vec);
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

	bool IsNormalized() const {
		T n = LengthSquared();
		return T(0.9999) <= n && n <= T(1.0001);
	}

	const T ScalarPart() const {
		return s;
	}
	const Vector<T, 3, Packed> VectorPart() const {
		return { x, y, z };
	}


	//-----------------------------------------------
	// Matrix conversions
	//-----------------------------------------------
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 3, 3, Order, Layout, Packed>() const {
		return ToMatrix<U, 3, 3, Order, Layout, Packed>();
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>() const {
		return ToMatrix<U, 3, 4, Order, Layout, Packed>();
	}
	template <class U, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>() const {
		return ToMatrix<U, 4, 3, Order, Layout, Packed>();
	}
	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	explicit operator Matrix<U, 4, 4, Order, Layout, Packed>() const {
		return ToMatrix<U, 4, 4, Order, Layout, Packed>();
	}

	//-----------------------------------------------
	// Truncate to vector
	//-----------------------------------------------
	template <class U, bool Packed>
	explicit operator Vector<U, 3, Packed>() const {
		return { x, y, z };
	}

	template <bool Packed>
	friend Vector<T, 3, Packed> operator*(const Vector<T, 3, Packed>& vec, const Quaternion<T, Packed>& q) {
		// sandwich product
		return (Vector<T, 3, Packed>)(q*Quaternion(vec)*q.Inverse());
	}

	template <bool Packed>
	friend Vector<T, 3, Packed>& operator*=(Vector<T, 3, Packed>& vec, const Quaternion<T, Packed>& q) {
		// sandwich product
		vec = (Vector<T, 3, Packed>)(q*Quaternion(vec)*q.Inverse());
		return vec;
	}


	//-----------------------------------------------
	// Apply to vector
	//-----------------------------------------------
	template <bool Packed>
	Vector<T, 3, Packed> operator*(const Vector<T, 3, Packed>& vec) const {
		// sandwich product
		return (*this)*Quaternion(vec)*Inverse();
	}

protected:
	//-----------------------------------------------
	// Matrix conversion helpers
	//-----------------------------------------------
	template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	Matrix<U, Rows, Columns, Order, Layout, Packed> ToMatrix() const {
		assert(IsNormalized());
		Matrix<U, Rows, Columns, Order, Layout, Packed> mat;
		auto elem = [&mat](int i, int j) -> U& {
			return Order == eMatrixOrder::PRECEDE_VECTOR ? mat(i, j) : mat(j, i);
		};
		elem(0, 0) = 1 - 2 * (j*j + k*k);	elem(0, 1) = 2 * (i*j - k*s);		elem(0, 2) = 2 * (i*k + j*s);
		elem(1, 0) = 2 * (i*j + k*s);		elem(1, 1) = 1 - 2 * (i*i + k*k);	elem(1, 2) = 2 * (j*k - i*s);
		elem(2, 0) = 2 * (i*k - j*s);		elem(2, 1) = 2 * (j*k + i*s);		elem(2, 2) = 1- 2 * (i*i + j*j);

		// Rest
		for (int j = 0; j < m.Width(); ++j) {
			for (int i = (j < 3 ? 3 : 0); i < m.Height(); ++i) {
				mat(i, j) = T(j == i);
			}
		}

		return mat;
	}

	template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	void FromMatrix(const Matrix<U, Rows, Columns, Order, Layout, Packed>& mat) const {
		assert(mat.IsRotationMatrix3D());
		auto elem = [&mat](int i, int j) -> U& {
			return Order == eMatrixOrder::PRECEDE_VECTOR ? mat(i, j) : mat(j, i);
		};
		w = std::sqrt(1 + elem(0, 0) + elem(1, 1) + elem(2, 2)) * T(0.5);
		T div = T(1) / (T(4) * w);
		x = (elem(2, 1) - elem(1, 2)) / div;
		y = (elem(0, 2) - elem(2, 0)) / div;
		z = (elem(1, 0) - elem(0, 1)) / div;
	}
};



template <class T, bool Packed, class U, class = typename std::enable_if<std::is_same<U, Quaternion<T, Packed>>::value>::type>
Quaternion<T, Packed> operator*(U s, const Quaternion<T, Packed>& rhs) {
	return Quaternion<T>(*this) * s;
}
template <class T, bool Packed, class U, class = typename std::enable_if<std::is_same<U, Quaternion<T, Packed>>::value>::type>
Quaternion<T, Packed> operator/(U s, const Quaternion<T, Packed>& rhs) {
	return Quaternion<T>(*this) / s;
}





} // namespace mathter