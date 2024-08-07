#pragma once


#ifdef MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif

#include <cmath>
#include <type_traits>


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
		using T = std::remove_reference_t<decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c))>;
		if constexpr (std::is_convertible_v<T1, T> && std::is_convertible_v<T2, T> && std::is_convertible_v<T3, T>) {
			return fma<T>{}(a, b, c);
		}
		return a * b + c;
	}
};


template <class T = void>
struct max {
	constexpr T operator()(const T& lhs, const T& rhs) const {
#ifdef MATHTER_ENABLE_SIMD
		if constexpr (xsimd::is_batch<T>::value) {
			return xsimd::max(lhs, rhs);
		}
		else {
#endif
			return std::max(lhs, rhs);
#ifdef MATHTER_ENABLE_SIMD
		}
#endif
	}
};


template <>
struct max<void> {
	template <class T1, class T2>
	constexpr auto operator()(T1&& lhs, T2&& rhs) const
		-> decltype(std::forward<T1>(lhs) + std::forward<T2>(rhs)) {
		using T = std::remove_reference_t<decltype(std::forward<T1>(lhs) + std::forward<T2>(rhs))>;
		return max<T>{}(static_cast<T>(lhs), static_cast<T>(rhs));
	}
};


template <class T = void>
struct min {
	constexpr T operator()(const T& lhs, const T& rhs) const {
#ifdef MATHTER_ENABLE_SIMD
		if constexpr (xsimd::is_batch<T>::value) {
			return xsimd::min(lhs, rhs);
		}
		else {
#endif
			return std::min(lhs, rhs);
#ifdef MATHTER_ENABLE_SIMD
		}
#endif
	}
};


template <>
struct min<void> {
	template <class T1, class T2>
	constexpr auto operator()(T1&& lhs, T2&& rhs) const
		-> decltype(std::forward<T1>(lhs) + std::forward<T2>(rhs)) {
		using T = std::remove_reference_t<decltype(std::forward<T1>(lhs) + std::forward<T2>(rhs))>;
		return min<T>{}(static_cast<T>(lhs), static_cast<T>(rhs));
	}
};


template <class T = void>
struct abs {
	constexpr T operator()(const T& arg) const {
#ifdef MATHTER_ENABLE_SIMD
		if constexpr (xsimd::is_batch<T>::value) {
			return xsimd::abs(arg);
		}
		else {
#endif
			return std::abs(arg);
#ifdef MATHTER_ENABLE_SIMD
		}
#endif
	}
};


template <>
struct abs<void> {
	template <class T>
	constexpr std::remove_reference_t<T> operator()(T&& arg) const {
		return abs<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};

} // namespace mathter