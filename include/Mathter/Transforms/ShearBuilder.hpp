// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Vector.hpp"
#include "IdentityBuilder.hpp"


namespace mathter {

namespace impl {

	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, int Dim>
	Matrix<T, Rows, Columns, Order, Layout, Packed> ShearMatrixGeneral(T slope, const Vector<T, Dim, Packed>& directionAxis, const Vector<T, Dim, Packed>& modulatorAxis) {
		using MNN = Matrix<T, Dim, Dim, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;
		using MN1 = Matrix<T, Dim, 1, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;
		using M1N = Matrix<T, 1, Dim, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;

		const MN1 d(directionAxis);
		const M1N mT(modulatorAxis);

		const auto shear = MNN(Identity()) + slope * d * mT;
		Matrix<T, Rows, Columns, Order, Layout, Packed> m = Identity();
		m.Insert(0, 0, Matrix<T, Dim, Dim, Order, Layout, Packed>(shear));
		return m;
	}


	template <class T>
	class ShearBuilderAligned {
	public:
		ShearBuilderAligned(T slope, int directionAxisIdx, int modulatorAxisIdx)
			: slope(slope), directionAxisIdx(directionAxisIdx), modulatorAxisIdx(modulatorAxisIdx) {}

		template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, Rows, Columns, Order, Layout, MPacked>() const {
			assert(directionAxisIdx != modulatorAxisIdx);

			Matrix<U, Rows, Columns, Order, Layout, MPacked> m = Identity();

			if constexpr (Order == eMatrixOrder::FOLLOW_VECTOR) {
				assert(modulatorAxisIdx < Rows);
				assert(directionAxisIdx < Columns);
				m(modulatorAxisIdx, directionAxisIdx) = slope;
			}
			else {
				assert(directionAxisIdx < Rows);
				assert(modulatorAxisIdx < Columns);
				m(directionAxisIdx, modulatorAxisIdx) = slope;
			}

			return m;
		}

	private:
		T slope;
		int directionAxisIdx;
		int modulatorAxisIdx;
	};


	template <class T, int Dim, bool Packed>
	class ShearBuilderGeneral {
	public:
		ShearBuilderGeneral(T slope, const Vector<T, Dim, Packed>& directionAxis, const Vector<T, Dim, Packed>& modulatorAxis)
			: slope(slope), directionAxis(directionAxis), modulatorAxis(modulatorAxis) {}


		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 4, Order, Layout, MPacked>() const {
			return ShearMatrixGeneral<U, 4, 4, Order, Layout, MPacked>(U(slope), Vector<U, Dim, MPacked>(directionAxis), Vector<U, Dim, MPacked>(modulatorAxis));
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 3, Order, Layout, MPacked>() const {
			return ShearMatrixGeneral<U, 3, 3, Order, Layout, MPacked>(U(slope), Vector<U, Dim, MPacked>(directionAxis), Vector<U, Dim, MPacked>(modulatorAxis));
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
			return ShearMatrixGeneral<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>(U(slope), Vector<U, Dim, MPacked>(directionAxis), Vector<U, Dim, MPacked>(modulatorAxis));
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
			return ShearMatrixGeneral<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>(U(slope), Vector<U, Dim, MPacked>(directionAxis), Vector<U, Dim, MPacked>(modulatorAxis));
		}

	private:
		T slope;
		Vector<T, Dim, Packed> directionAxis;
		Vector<T, Dim, Packed> modulatorAxis;
	};

} // namespace impl


/// <summary> Creates a shear matrix aligned to a certain axis. </summary>
/// <param name="slope"> Strength of the shear. </param>
/// <param name="directionAxisIdx"> Points are moved along this axis. </param>
/// <param name="modulatorAxisIdx"> The displacement of points is proportional to this coordinate's value. </param>
/// <remarks> The formula for displacement along the pricipal axis is
///		<paramref name="slope"/>&ast;pos[<paramref name="modulatorAxisIdx"/>]. </remarks>
template <class T>
auto Shear(T slope, int directionAxisIdx, int modulatorAxisIdx) {
	return impl::ShearBuilderAligned(slope, directionAxisIdx, modulatorAxisIdx);
}


/// <summary> Creates a general shear matrix. </summary>
/// <param name="slope"> Strength of the shear. </param>
/// <param name="directionAxis"> Points are moved along this axis. Must be a unit vector. </param>
/// <param name="modulatorAxis"> The displacement of points is proportional to distance from origin projected to this axis. Must be a unit vector. </param>
/// <remarks> <paramref name="directionAxis"/> and <paramref name="modulatorAxis"/> must be perpendicular. </remarks>
template <class T, int Dim, bool Packed>
auto Shear(T slope, const Vector<T, Dim, Packed>& directionAxis, const Vector<T, Dim, Packed>& modulatorAxis) {
	return impl::ShearBuilderGeneral(slope, directionAxis, modulatorAxis);
}

} // namespace mathter