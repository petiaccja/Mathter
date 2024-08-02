// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Vector/Math.hpp"
#include "../Vector/Vector.hpp"
#include "IdentityBuilder.hpp"

#include <algorithm>


namespace mathter {


template <class T, int Dim, bool Packed>
class ViewBuilder {
	using VectorT = Vector<T, Dim, Packed>;

public:
	ViewBuilder(const VectorT& eye, const VectorT& target, const std::array<VectorT, size_t(Dim - 2)>& bases, const std::array<bool, Dim>& flipAxes)
		: eye(eye), target(target), bases(bases), flipAxes(flipAxes) {}

	template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
	operator Matrix<U, Dim + 1, Dim + 1, Order, Layout, MPacked>() const {
		Matrix<U, Dim + 1, Dim + 1, Order, Layout, MPacked> m;
		Set(m);
		return m;
	}

	template <class U, eMatrixLayout Layout, bool MPacked>
	operator Matrix<U, Dim, Dim + 1, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
		Matrix<U, Dim, Dim + 1, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked> m;
		Set(m);
		return m;
	}

	template <class U, eMatrixLayout Layout, bool MPacked>
	operator Matrix<U, Dim + 1, Dim, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
		Matrix<U, Dim + 1, Dim, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked> m;
		Set(m);
		return m;
	}


private:
	template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
	void Set(Matrix<U, Rows, Columns, Order, Layout, MPacked>& m) const {
		Matrix<U, Dim, Dim, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, MPacked> rotation;

		std::decay_t<decltype(rotation.stripes)> extendedBases;
		auto& columns = rotation.stripes;

		// Arrange vectors like [side, look, bases...]
		const auto lookAxis = eye - target;
		extendedBases[Dim - 1] = lookAxis;
		std::copy(bases.begin(), bases.end(), extendedBases.begin() + 1);
		const auto sideAxis = Cross(extendedBases.begin() + 1, extendedBases.end());
		std::rotate(extendedBases.begin() + 1, extendedBases.end() - 1, extendedBases.end());
		extendedBases[0] = sideAxis;

		// Do Gram-Schmidt orthogonalization.
		GramSchmidtOrthogonalize(extendedBases.begin(), extendedBases.end(), rotation.stripes.begin());

		// Normalize columns.
		std::for_each(columns.begin(), columns.end(), [](auto& v) {
			v = NormalizePrecise(v);
		});

		// Reorder columns to [side, bases... look]
		std::rotate(columns.begin() + 1, columns.begin() + 2, columns.end());

		// Flip axes
		for (size_t i = 0; i < Dim; ++i) {
			columns[i] = flipAxes[i] ? -columns[i] : columns[i];
		}

		const auto translation = -eye * rotation;

		m = Identity();
		m.Insert(0, 0, Matrix<U, Dim, Dim, Order, Layout, MPacked>(rotation));
		if constexpr (Order == eMatrixOrder::FOLLOW_VECTOR) {
			m.Insert(Dim, 0, Matrix<U, 1, Dim, Order, Layout, MPacked>(translation));
		}
		else {
			m.Insert(0, Dim, Matrix<U, Dim, 1, Order, Layout, MPacked>(translation));
		}
	}

	const VectorT eye;
	const VectorT target;
	const std::array<VectorT, size_t(Dim - 2)> bases;
	const std::array<bool, Dim> flipAxes;
};


/// <summary> Creates a general, n-dimensional camera look-at matrix. </summary>
/// <param name="eye"> The camera's position. </param>
/// <param name="target"> The camera's target. </param>
/// <param name="bases"> Basis vectors fixing the camera's orientation. </param>
/// <param name="flipAxis"> Set any element to true to flip an axis in camera space. </param>
/// <remarks> The camera looks down the vector going from <paramref name="eye"/> to
///		<paramref name="target"/>, but it can still rotate around that vector. To fix the rotation,
///		an "up" vector must be provided in 3 dimensions. In higher dimensions,
///		we need multiple up vectors. Unfortunately I can't remember how these
///		basis vectors are used, but they are orthogonalized to each-other and to the look vector.
///		I can't remember the order of orthogonalization. </remarks>
template <class T, int Dim, bool Packed, size_t BaseDim, size_t FlipDim>
auto LookAt(const Vector<T, Dim, Packed>& eye,
			const Vector<T, Dim, Packed>& target,
			const std::array<Vector<T, Dim, Packed>, BaseDim>& bases,
			const std::array<bool, FlipDim>& flipAxes) {
	static_assert(BaseDim == Dim - 2, "You must provide 2 fewer bases than the dimension of the transform.");
	static_assert(Dim == FlipDim, "You must provide the same number of flips as the dimension of the transform.");
	return ViewBuilder<T, Dim, Packed>(eye, target, bases, flipAxes);
}


/// <summary> Creates a 2D look-at matrix. </summary>
/// <param name="eye"> The camera's position. </param>
/// <param name="target"> The camera's target. </param>
/// <param name="positiveYForward"> True if the camera looks towards +Y in camera space, false if -Y. </param>
/// <param name="flipX"> True to flip X in camera space. </param>
template <class T, bool Packed>
auto LookAt(const Vector<T, 2, Packed>& eye, const Vector<T, 2, Packed>& target, bool positiveYForward, bool flipX) {
	return LookAt(eye, target, std::array<Vector<T, 2, Packed>, 0>{}, std::array{ flipX, positiveYForward });
}


/// <summary> Creates a 3D look-at matrix. </summary>
/// <param name="eye"> The camera's position. </param>
/// <param name="target"> The camera's target. </param>
/// <param name="up"> Up direction in world space. </param>
/// <param name="positiveZForward"> True if the camera looks towards +Z in camera space, false if -Z. </param>
/// <param name="flipX"> True to flip X in camera space. </param>
/// <param name="flipY"> True to flip Y in camera space. </param>
/// <remarks> The camera space X is selected to be orthogonal to both the look direction and the <paramref name="up"/> vector.
///		Afterwards, the <paramref name="up"/> vector is re-orthogonalized to the camera-space Z and X vectors. </remarks>
template <class T, bool Packed>
auto LookAt(const Vector<T, 3, Packed>& eye, const Vector<T, 3, Packed>& target, const Vector<T, 3, Packed>& up, bool positiveZForward, bool flipX, bool flipY) {
	return LookAt(eye, target, std::array<Vector<T, 3, Packed>, 1>{ up }, std::array<bool, 3>{ flipX, flipY, positiveZForward });
}



} // namespace mathter