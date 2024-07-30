// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Vector.hpp"
#include "IdentityBuilder.hpp"

#include <cmath>


namespace mathter {

namespace impl {

	template <class T>
	class Rotation2DBuilder {
	public:
		Rotation2DBuilder(const T& angle) : angle(angle) {}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 3, Order, Layout, MPacked>() const {
			Matrix<U, 3, 3, Order, Layout, MPacked> m;
			Set(m);
			return m;
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 2, 2, Order, Layout, MPacked>() const {
			Matrix<U, 2, 2, Order, Layout, MPacked> m;
			Set(m);
			return m;
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 2, 3, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
			Matrix<U, 2, 3, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked> m;
			Set(m);
			return m;
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 2, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
			Matrix<U, 3, 2, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked> m;
			Set(m);
			return m;
		}

	private:
		template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		void Set(Matrix<U, Rows, Columns, Order, Layout, MPacked>& m) const {
			m = Identity();

			const auto C = std::cos(angle);
			const auto S = std::sin(angle);

			auto elem = [&m](int i, int j) -> U& {
				return Order == eMatrixOrder::FOLLOW_VECTOR ? m(i, j) : m(j, i);
			};

			// Indices according to follow vector order
			elem(0, 0) = U(C);
			elem(0, 1) = U(S);
			elem(1, 0) = U(-S);
			elem(1, 1) = U(C);
		}

		T angle;
	};

} // namespace impl


/// <summary> Creates a 2D rotation matrix. </summary>
/// <param name="angle"> Counter-clockwise angle in radians. </param>
template <class T>
auto Rotation(const T& angle) {
	return impl::Rotation2DBuilder{ angle };
}

} // namespace mathter