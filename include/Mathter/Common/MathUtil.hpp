// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include <tuple>


namespace mathter {

template <class T>
std::tuple<T, T> Fast2Sum(const T& a, const T& b) {
	const auto s = a + b;
	const auto z = s - a;
	const auto c = b - z;
	return { s, c };
}

} // namespace mathter