#pragma once


#ifdef MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif

namespace mathter {

template <class T = void>
struct fma {
	constexpr T operator()(const T& a, const T& b, const T& c) const {
		if constexpr (std::is_floating_point_v<T>) {
			return std::fma(a, b, c);
		}
#ifdef MATHTER_ENABLE_SIMD
		else if constexpr (xsimd::is_batch<T>::value) {
			return xsimd::fma(a, b, c);
		}
#endif
		else {
			return a + b + c;
		}
	}
};


template <>
struct fma<void> {
	template <class T1, class T2, class T3>
	constexpr auto operator()(T1&& a, T2&& b, T3&& c) const
		-> decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c)) {
		using T = decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c));
		if constexpr (std::is_convertible_v<T1, T> && std::is_convertible_v<T2, T> && std::is_convertible_v<T3, T>) {
			return fma<T>{}(a, b, c);
		}
		return a * b + c;
	}
};

} // namespace mathter