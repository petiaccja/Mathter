// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once


#include "../Common/OptimizationUtil.hpp"
#include "../Matrix/Matrix.hpp"
#include "../Quaternion/Quaternion.hpp"
#include "../Vector/Vector.hpp"
#include "ZeroBuilder.hpp"


namespace mathter {

namespace impl {

	class ZeroBuilder {
	public:
		template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
		operator Matrix<T, Rows, Columns, Order, Layout, Packed>() const {
			using Mat = Matrix<T, Rows, Columns, Order, Layout, Packed>;
			Mat m;
			for (auto& stripe : m.stripes) {
				stripe = Vector<T, Mat::stripeDim, Packed>(static_cast<T>(0));
			}
			return m;
		}

		template <class T, int Dim, bool Packed>
		operator Vector<T, Dim, Packed>() const {
			return Vector<T, Dim, Packed>(T(0));
		}

		template <class T, eQuaternionLayout Layout, bool Packed>
		operator Quaternion<T, Layout, Packed>() const {
			return Quaternion<T, Layout, Packed>(T(0), T(0), T(0), T(0));
		}
	};

} // namespace impl


/// <summary> Creates a matrix/vector/quaternion with all elements zero. </summary>
inline auto Zero() {
	return impl::ZeroBuilder{};
}


} // namespace mathter