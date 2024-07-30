// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Quaternion/Quaternion.hpp"
#include "ZeroBuilder.hpp"


namespace mathter {

namespace impl {

	class IdentityBuilder {
	public:
		template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
		operator Matrix<T, Rows, Columns, Order, Layout, Packed>() const {
			Matrix<T, Rows, Columns, Order, Layout, Packed> m = Zero();
			for (int i = 0; i < std::min(Rows, Columns); ++i) {
				m(i, i) = T(1);
			}
			return m;
		}

		template <class T, eQuaternionLayout Layout, bool Packed>
		operator Quaternion<T, Layout, Packed>() const {
			return Quaternion<T, Layout, Packed>{ T(1), T(0), T(0), T(0) };
		}
	};

} // namespace impl


/// <summary> Creates an identity matrix or identity quaternion. </summary>
/// <remarks> If the matrix is not square, it will look like a truncated larger square identity matrix. </remarks>
inline auto Identity() {
	return impl::IdentityBuilder{};
}

} // namespace mathter