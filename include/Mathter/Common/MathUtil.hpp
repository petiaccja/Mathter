// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include <cmath>
#include <complex>
#include <tuple>
#include <type_traits>


namespace mathter {

template <class T>
std::tuple<T, T> Fast2Sum(const T& a, const T& b) {
	const auto s = a + b;
	const auto z = s - a;
	const auto c = b - z;
	return { s, c };
}


template <class FloatingPoint>
auto Normalize(FloatingPoint x) -> std::enable_if_t<std::is_floating_point_v<FloatingPoint>, FloatingPoint> {
	return std::copysign(FloatingPoint(1), x);
}


template <class FloatingPoint>
auto Normalize(const std::complex<FloatingPoint>& z) -> std::complex<FloatingPoint> {
	const auto absReal = std::abs(std::real(z));
	const auto absImag = std::abs(std::imag(z));
	const auto prescaler = FloatingPoint(1) / std::max(absReal, absImag);
	const auto zp = z * prescaler; // Crude prescaling necessary because std::abs(z) may be inaccurate for denormal inputs.
	const auto mag = std::sqrt(std::norm(zp)); // Don't use std::abs, std::sqrt cannot overflow due to scaling.
	return std::isfinite(mag) ? zp / mag : std::complex<FloatingPoint>(1);
}

} // namespace mathter