// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Vector/Vector.hpp"
#include "ZeroBuilder.hpp"

#include <array>
#include <utility>


namespace mathter {

namespace impl {

	template <class T, int Dim>
	class ScaleBuilder {
	public:
		explicit ScaleBuilder(const std::array<T, Dim>& scale) : scale(scale) {}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, Dim + 1, Dim + 1, Order, Layout, MPacked>() const {
			Matrix<U, Dim + 1, Dim + 1, Order, Layout, MPacked> m;
			Set(m);
			return m;
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, Dim, Dim, Order, Layout, MPacked>() const {
			Matrix<U, Dim, Dim, Order, Layout, MPacked> m;
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
			m = Zero();
			size_t i;
			for (i = 0; i < scale.size(); ++i) {
				m(i, i) = static_cast<U>(scale[i]);
			}
			for (; i < std::min(Rows, Columns); ++i) {
				m(i, i) = static_cast<U>(1);
			}
		}

		std::array<T, Dim> scale;
	};

	template <class T, size_t N>
	ScaleBuilder(const std::array<T, N>&) -> ScaleBuilder<T, int(N)>;

} // namespace impl


/// <summary> Creates a scaling matrix. </summary>
/// <param name="scale"> A vector containing the scales of respective axes. </summary>
template <class T, int Dim, bool Packed>
auto Scale(const Vector<T, Dim, Packed>& scale) {
	std::array<T, Dim> arr;
	std::copy(scale.begin(), scale.end(), arr.begin());
	return impl::ScaleBuilder{ arr };
}


/// <summary> Creates a scaling matrix. </summary>
/// <param name="scales"> A list of scalars corresponding to scaling on respective axes. </summary>
template <class... Scales, std::enable_if_t<(... && is_scalar_v<std::decay_t<Scales>>), int> = 0>
auto Scale(const Scales&... scales) {
	using T = common_arithmetic_type_t<std::decay_t<Scales>...>;
	return impl::ScaleBuilder{ std::array{ static_cast<T>(scales)... } };
}


/// <summary> Creates a scaling matrix. </summary>
/// <param name="scale"> An array containing the scales of respective axes. </summary>
template <class T, int Dim>
auto Scale(const std::array<T, Dim>& scale) {
	return impl::ScaleBuilder{ scale };
}

} // namespace mathter