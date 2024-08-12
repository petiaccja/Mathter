// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Vector.hpp"
#include "IdentityBuilder.hpp"

#include <array>
#include <utility>


namespace mathter {

namespace impl {

	template <class T, int Dim>
	class TranslationBuilder {
	public:
		explicit TranslationBuilder(const std::array<T, Dim>& translation) : translation(translation) {}

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
			m = Identity();
			if constexpr (Order == eMatrixOrder::FOLLOW_VECTOR) {
				for (size_t i = 0; i < translation.size(); ++i) {
					m(Rows - 1, i) = static_cast<U>(translation[i]);
				}
			}
			else {
				for (size_t i = 0; i < translation.size(); ++i) {
					m(i, Columns - 1) = static_cast<U>(translation[i]);
				}
			}
		}

		std::array<T, Dim> translation;
	};

	template <class T, size_t N>
	TranslationBuilder(const std::array<T, N>&) -> TranslationBuilder<T, int(N)>;

} // namespace impl


/// <summary> Creates a translation matrix. </summary>
/// <param name="scale"> A vector containing the translation along respective axes. </summary>
template <class T, int Dim, bool Packed>
auto Translation(const Vector<T, Dim, Packed>& scale) {
	std::array<T, Dim> arr;
	std::copy(scale.begin(), scale.end(), arr.begin());
	return impl::TranslationBuilder{ arr };
}


/// <summary> Creates a translation matrix. </summary>
/// <param name="scales"> A list of scalars corresponding to translation along respective axes. </summary>
template <class... Scales, std::enable_if_t<(... && is_scalar_v<std::decay_t<Scales>>), int> = 0>
auto Translation(const Scales&... scales) {
	using T = common_arithmetic_type_t<std::decay_t<Scales>...>;
	return impl::TranslationBuilder{ std::array{ static_cast<T>(scales)... } };
}


/// <summary> Creates a translation matrix. </summary>
/// <param name="scale"> An array containing the translation along respective axes. </summary>
template <class T, int Dim>
auto Translation(const std::array<T, Dim>& scale) {
	return impl::TranslationBuilder{ scale };
}

} // namespace mathter