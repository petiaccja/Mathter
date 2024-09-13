// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/Types.hpp"
#include "../Matrix/Matrix.hpp"
#include "../Vector/Swizzle.hpp"
#include "../Vector/Vector.hpp"


namespace mathter {


struct CanonicalArg {};
constexpr CanonicalArg canonicalArg;


/// <summary> Allows you to do quaternion math and represent rotation in a compact way. </summary>
/// <typeparam name="T"> The scalar type of w, x, y and z. Use a builtin or custom floating or fixed point type. </typeparam>
/// <typeparam name="Packed"> If true, tightly packs quaternion members and disables padding due to overalignment in arrays.
///		Disables SIMD optimization. </typeparam>
/// <remarks>
/// These are plain mathematical quaternions, so expect the operations to work as mathematically defined.
/// There are helper functions to represent rotation with quaternions.
/// </remarks>
template <class T, eQuaternionLayout Layout = eQuaternionLayout::SCALAR_FIRST, bool Packed = false>
class Quaternion {
	static constexpr int scalarIdx = Layout == eQuaternionLayout::SCALAR_FIRST ? 0 : 3;
	static constexpr int vectorIdxI = Layout == eQuaternionLayout::SCALAR_FIRST ? 1 : 0;
	static constexpr int vectorIdxJ = Layout == eQuaternionLayout::SCALAR_FIRST ? 2 : 1;
	static constexpr int vectorIdxK = Layout == eQuaternionLayout::SCALAR_FIRST ? 3 : 2;

public:
	union {
		Elements<T, 4, Packed> elements;
		Swizzle<T, 4, Packed, scalarIdx> s; // The scalar imaginary element.
		Swizzle<T, 4, Packed, vectorIdxI> i; // The 1st imaginary element.
		Swizzle<T, 4, Packed, vectorIdxJ> j; // The 2nd imaginary element.
		Swizzle<T, 4, Packed, vectorIdxK> k; // The 3rd imaginary element.
		Swizzle<T, 4, Packed, scalarIdx> w; // The scalar imaginary element.
		Swizzle<T, 4, Packed, vectorIdxI> x; // The 1st imaginary element.
		Swizzle<T, 4, Packed, vectorIdxJ> y; // The 2nd imaginary element.
		Swizzle<T, 4, Packed, vectorIdxK> z; // The 3rd imaginary element.
		Swizzle<T, 4, Packed, scalarIdx> scalar; // The scalar part of the quaternion.
		Swizzle<T, 4, Packed, vectorIdxI, vectorIdxJ, vectorIdxK> vector; // The vector part of the quaternion.
		Swizzle<T, 4, Packed, scalarIdx, vectorIdxI, vectorIdxJ, vectorIdxK> canonical; // Canonical memory layout. Mainly for internal use.
	};

	//-----------------------------------------------
	// Constructors
	//-----------------------------------------------

	/// <summary> Does NOT zero-initialize values. </summary>
	Quaternion() : Quaternion(Vector<T, 4, Packed>{}) {}

	/// <summary> Convert from 4-dimensional vector. </summary>
	/// <remarks> The vector must share the layout of the quaternion. Use with care. </remarks>
	template <class TOther, bool PackedOther>
	explicit Quaternion(const Vector<TOther, 4, PackedOther>& vector);

	/// <summary> Convert from 4-dimensional vector. </summary>
	/// <remarks> The vector must have SIJK layout. Use with care. </remarks>
	template <class TOther, bool PackedOther>
	explicit Quaternion(CanonicalArg, const Vector<TOther, 4, PackedOther>& vector);

	/// <summary> Convert from a quaternion of different (or same) type. </summary>
	template <class TOther, eQuaternionLayout LayoutOther, bool PackedOther>
	Quaternion(const Quaternion<TOther, LayoutOther, PackedOther>& rhs);

	/// <summary> Set values directly. </summary>
	Quaternion(const T& scalar, const T& x = T(0), const T& y = T(0), const T& z = T(0));

	/// <summary> Sets the scalar part (w) and the vector part (xyz). This is not <see cref="AxisAngle"/> rotation. </summary>
	Quaternion(const T& scalar, const Vector<T, 3, Packed>& vector);

	/// <summary> Sets the scalar part to zero, and the vector part to given argument. </summary>
	explicit Quaternion(const Vector<T, 3, Packed>& vector);

	/// <summary> Convert a rotation matrix to equivalent quaternion. </summary>
	/// <remarks> Matrix must be in SO(3). </remarks>
	template <class TOther, int RowsOther, int ColumnsOther, eMatrixOrder OrderOther, eMatrixLayout LayoutOther, bool PackedOther>
	explicit Quaternion(const Matrix<TOther, RowsOther, ColumnsOther, OrderOther, LayoutOther, PackedOther>& matrix);


	//-----------------------------------------------
	// Conversion operators
	//-----------------------------------------------

	/// <summary> Converts the quaternion into a 3D rotation matrix. </summary>
	template <class TOther, int RowsOther, int ColumnsOther, eMatrixOrder OrderOther, eMatrixLayout LayoutOther, bool PackedOther>
	explicit operator Matrix<TOther, RowsOther, ColumnsOther, OrderOther, LayoutOther, PackedOther>() const;

	/// <summary> Converts the quaternion to a 4-dimensional vector. </summary>
	///	<remarks> The vector shares the layout of the quaternion. Use with care. </remarks>
	template <class TOther, bool PackedOther>
	explicit operator Vector<TOther, 4, PackedOther>() const;


	//-----------------------------------------------
	// Versor
	//-----------------------------------------------

	/// <summary> Rotate a vector by this unit quaternion. </summary>
	/// <remarks> If you get an undefined symbol for this, then the definition is in `RotationArithmetic.hpp`. </remarks>
	template <class TOther, bool PackedOther>
	auto operator()(const Vector<TOther, 3, PackedOther>& v) const;


	//-----------------------------------------------
	// Properties
	//-----------------------------------------------

	/// <summary> Returns the scalar part (w) of (w + xi + yj + zk). </summary>
	[[deprecated("use .scalar")]] T ScalarPart() const;

	/// <summary> Returns the vector part (x, y, z) of (w + xi + yj + zk). </summary>
	[[deprecated("use .vector")]] Vector<T, 3, Packed> VectorPart() const;
};


//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

namespace impl {

	template <class TM, int RowsM, int ColumnsM, eMatrixOrder OrderM, eMatrixLayout LayoutM, bool PackedM,
			  class TQ, eQuaternionLayout LayoutQ, bool PackedQ>
	void QuaternionToMatrix(Matrix<TM, RowsM, ColumnsM, OrderM, LayoutM, PackedM>& m,
							const Quaternion<TQ, LayoutQ, PackedQ>& q) {
		const auto elem = [&m](int i, int j) -> auto& {
			return OrderM == eMatrixOrder::PRECEDE_VECTOR ? m(i, j) : m(j, i);
		};
		elem(0, 0) = static_cast<TM>(TQ(1) - TQ(2) * (q.j * q.j + q.k * q.k));
		elem(0, 1) = static_cast<TM>(TQ(2) * (q.i * q.j - q.k * q.s));
		elem(0, 2) = static_cast<TM>(TQ(2) * (q.i * q.k + q.j * q.s));
		elem(1, 0) = static_cast<TM>(TQ(2) * (q.i * q.j + q.k * q.s));
		elem(1, 1) = static_cast<TM>(TQ(1) - TQ(2) * (q.i * q.i + q.k * q.k));
		elem(1, 2) = static_cast<TM>(TQ(2) * (q.j * q.k - q.i * q.s));
		elem(2, 0) = static_cast<TM>(TQ(2) * (q.i * q.k - q.j * q.s));
		elem(2, 1) = static_cast<TM>(TQ(2) * (q.j * q.k + q.i * q.s));
		elem(2, 2) = static_cast<TM>(TQ(1) - TQ(2) * (q.i * q.i + q.j * q.j));

		// Fill the rest as an identity matrix.
		for (int j = 0; j < m.Width(); ++j) {
			for (int i = (j < 3 ? 3 : 0); i < m.Height(); ++i) {
				m(i, j) = static_cast<TM>(j == i);
			}
		}
	}

	template <class TM, int RowsM, int ColumnsM, eMatrixOrder OrderM, eMatrixLayout LayoutM, bool PackedM,
			  class TQ, eQuaternionLayout LayoutQ, bool PackedQ>
	void MatrixToQuaternion(Quaternion<TQ, LayoutQ, PackedQ>& q,
							const Matrix<TM, RowsM, ColumnsM, OrderM, LayoutM, PackedM>& m) {
		auto elem = [&m](int i, int j) -> const auto& {
			return OrderM == eMatrixOrder::PRECEDE_VECTOR ? m(i, j) : m(j, i);
		};
		const auto s = std::sqrt(TM(1) + elem(0, 0) + elem(1, 1) + elem(2, 2)) * TM(0.5);
		const auto div = TM(1) / (TM(4) * s);
		const auto i = (elem(2, 1) - elem(1, 2)) * div;
		const auto j = (elem(0, 2) - elem(2, 0)) * div;
		const auto k = (elem(1, 0) - elem(0, 1)) * div;
		q.canonical = Vector(s, i, j, k);
	}

} // namespace impl


//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, bool PackedOther>
Quaternion<T, Layout, Packed>::Quaternion(const Vector<TOther, 4, PackedOther>& vector) {
	elements = Vector<T, 4, Packed>(vector).elements;
}


template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, bool PackedOther>
Quaternion<T, Layout, Packed>::Quaternion(CanonicalArg, const Vector<TOther, 4, PackedOther>& vector) {
	canonical = vector;
}


template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, eQuaternionLayout LayoutOther, bool PackedOther>
Quaternion<T, Layout, Packed>::Quaternion(const Quaternion<TOther, LayoutOther, PackedOther>& rhs) {
	canonical = rhs.canonical;
}


template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed>::Quaternion(const T& scalar, const T& x, const T& y, const T& z) {
	canonical = Vector{ scalar, x, y, z };
}


template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed>::Quaternion(const T& scalar, const Vector<T, 3, Packed>& vector) {
	canonical = Vector(scalar, vector);
}


template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed>::Quaternion(const Vector<T, 3, Packed>& vector) {
	canonical = Vector(static_cast<T>(0), vector);
}


template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, int RowsOther, int ColumnsOther, eMatrixOrder OrderOther, eMatrixLayout LayoutOther, bool PackedOther>
Quaternion<T, Layout, Packed>::Quaternion(const Matrix<TOther, RowsOther, ColumnsOther, OrderOther, LayoutOther, PackedOther>& matrix) {
	static_assert(RowsOther >= 3 && ColumnsOther >= 3, "Matrix does not represent a 3D rotation.");
	impl::MatrixToQuaternion(*this, matrix);
}


template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, int RowsOther, int ColumnsOther, eMatrixOrder OrderOther, eMatrixLayout LayoutOther, bool PackedOther>
Quaternion<T, Layout, Packed>::operator Matrix<TOther, RowsOther, ColumnsOther, OrderOther, LayoutOther, PackedOther>() const {
	static_assert(RowsOther >= 3 && ColumnsOther >= 3, "Matrix cannot represent a 3D rotation.");
	Matrix<TOther, RowsOther, ColumnsOther, OrderOther, LayoutOther, PackedOther> matrix;
	impl::QuaternionToMatrix(matrix, *this);
	return matrix;
}


template <class T, eQuaternionLayout Layout, bool Packed>
template <class TOther, bool PackedOther>
Quaternion<T, Layout, Packed>::operator Vector<TOther, 4, PackedOther>() const {
	using VecTarget = Vector<TOther, 4, PackedOther>;
	using VecNative = Vector<T, 4, Packed>;
	return VecTarget(VecNative(elements.array.data()));
}


template <class T, eQuaternionLayout Layout, bool Packed>
Vector(Quaternion<T, Layout, Packed>) -> Vector<T, 4, Packed>;


template <class T, eQuaternionLayout Layout, bool Packed>
T Quaternion<T, Layout, Packed>::ScalarPart() const {
	return scalar;
}


template <class T, eQuaternionLayout Layout, bool Packed>
Vector<T, 3, Packed> Quaternion<T, Layout, Packed>::VectorPart() const {
	return vector;
}

} // namespace mathter