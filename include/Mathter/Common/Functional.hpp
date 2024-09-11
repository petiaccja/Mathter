#pragma once


#ifdef MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif

#include "TypeTraits.hpp"

#include <cmath>
#include <type_traits>


namespace mathter {

template <class T, class A>
struct supports_fma {
#ifdef MATHTER_ENABLE_SIMD
	template <class T_, class A_>
	static constexpr auto get(std::nullptr_t) -> std::enable_if_t<xsimd::fma3<A_>::supported(), bool> {
		return true;
	}
#endif

	template <class T_, class A_>
	static constexpr bool get(const void*) {
		return false;
	}

	static constexpr bool value = get<T, A>(nullptr);
};


template <class T = void>
struct madd {
	constexpr T operator()(const T& a, const T& b, const T& c) const {
		if constexpr (std::is_floating_point_v<T>) {
#if defined(FP_FAST_FMA) && defined(FP_FAST_FMAF) && defined(FP_FAST_FMAL)
			return std::fma(a, b, c);
#else
			return a * b + c;
#endif
		}
		else {
			return a * b + c;
		}
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct madd<xsimd::batch<T, A>> {
	xsimd::batch<T, A> operator()(const xsimd::batch<T, A>& a, const xsimd::batch<T, A>& b, const xsimd::batch<T, A>& c) const {
		return xsimd::fma(a, b, c);
	}
};
#endif


template <>
struct madd<void> {
	template <class T1, class T2, class T3>
	constexpr auto operator()(T1&& a, T2&& b, T3&& c) const
		-> decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c)) {
		using T = std::remove_reference_t<decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c))>;
		return madd<T>{}(static_cast<T>(std::forward<T1>(a)),
						 static_cast<T>(std::forward<T2>(b)),
						 static_cast<T>(std::forward<T3>(c)));
	}
};


template <class T = void>
struct fma {
	constexpr T operator()(const T& a, const T& b, const T& c) const {
		if constexpr (std::is_floating_point_v<T>) {
			return std::fma(a, b, c);
		}
		else {
			return a * b + c;
		}
	}
};


#ifdef MATHTER_ENABLE_SIMD

template <class T, class A, bool>
struct fma_base_simd {};


template <class T, class A>
struct fma_base_simd<T, A, true> {
	auto operator()(const xsimd::batch<T, A>& a, const xsimd::batch<T, A>& b, const xsimd::batch<T, A>& c) const {
		return xsimd::fma(a, b, c);
	}
};


template <class T, class A>
struct fma<xsimd::batch<T, A>> : fma_base_simd<T, A, supports_fma<T, A>::value> {};

#endif


template <>
struct fma<void> {
	template <class T1, class T2, class T3>
	constexpr auto operator()(T1&& a, T2&& b, T3&& c) const
		-> std::invoke_result_t<fma<common_arithmetic_type_t<T1, T2, T3>>,
								common_arithmetic_type_t<T1, T2, T3>,
								common_arithmetic_type_t<T1, T2, T3>,
								common_arithmetic_type_t<T1, T2, T3>> {
		using T = std::remove_reference_t<decltype(std::forward<T1>(a) * std::forward<T2>(b) + std::forward<T3>(c))>;
		return fma<T>{}(static_cast<T>(std::forward<T1>(a)),
						static_cast<T>(std::forward<T2>(b)),
						static_cast<T>(std::forward<T3>(c)));
	}
};


template <class T = void>
struct max {
	constexpr T operator()(const T& lhs, const T& rhs) const {
		return std::max(lhs, rhs);
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct max<xsimd::batch<T, A>> {
	xsimd::batch<T, A> operator()(const xsimd::batch<T, A>& lhs, const xsimd::batch<T, A>& rhs) const {
		return xsimd::max(lhs, rhs);
	}
};
#endif


template <>
struct max<void> {
	template <class T1, class T2>
	constexpr auto operator()(T1&& lhs, T2&& rhs) const
		-> common_arithmetic_type_t<std::decay_t<T1>, std::decay_t<T2>> {
		using T = common_arithmetic_type_t<std::decay_t<T1>, std::decay_t<T2>>;
		return max<T>{}(static_cast<T>(std::forward<T1>(lhs)), static_cast<T>(std::forward<T2>(rhs)));
	}
};


template <class T = void>
struct min {
	constexpr T operator()(const T& lhs, const T& rhs) const {
		return std::min(lhs, rhs);
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct min<xsimd::batch<T, A>> {
	xsimd::batch<T, A> operator()(const xsimd::batch<T, A>& lhs, const xsimd::batch<T, A>& rhs) const {
		return xsimd::min(lhs, rhs);
	}
};
#endif


template <>
struct min<void> {
	template <class T1, class T2>
	constexpr auto operator()(T1&& lhs, T2&& rhs) const
		-> common_arithmetic_type_t<std::decay_t<T1>, std::decay_t<T2>> {
		using T = common_arithmetic_type_t<std::decay_t<T1>, std::decay_t<T2>>;
		return min<T>{}(static_cast<T>(std::forward<T1>(lhs)), static_cast<T>(std::forward<T2>(rhs)));
	}
};


template <class T = void>
struct abs {
	constexpr auto operator()(const T& arg) const {
		return std::abs(arg);
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct abs<xsimd::batch<T, A>> {
	xsimd::batch<remove_complex_t<T>, A> operator()(const xsimd::batch<T, A>& arg) const {
		return xsimd::abs(arg);
	}
};
#endif


template <>
struct abs<void> {
	template <class T>
	constexpr auto operator()(T&& arg) const {
		return abs<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};


template <class T = void>
struct real {
	constexpr auto operator()(const T& arg) const {
		if constexpr (is_complex_v<std::decay_t<T>>) {
			return std::real(arg);
		}
		else {
			return arg;
		}
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct real<xsimd::batch<T, A>> {
	xsimd::batch<remove_complex_t<T>, A> operator()(const xsimd::batch<T, A>& arg) const {
		return xsimd::real(arg);
	}
};
#endif


template <>
struct real<void> {
	template <class T>
	constexpr auto operator()(T&& arg) const {
		return real<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};


template <class T = void>
struct imag {
	constexpr auto operator()(const T& arg) const {
		if constexpr (is_complex_v<std::decay_t<T>>) {
			return std::imag(arg);
		}
		else {
			return static_cast<T>(0);
		}
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct imag<xsimd::batch<T, A>> {
	xsimd::batch<remove_complex_t<T>, A> operator()(const xsimd::batch<T, A>& arg) const {
		return xsimd::imag(arg);
	}
};
#endif


template <>
struct imag<void> {
	template <class T>
	constexpr auto operator()(T&& arg) const {
		return imag<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};


template <class T = void>
struct conj {
	constexpr auto operator()(const T& arg) const {
		if constexpr (is_complex_v<std::decay_t<T>>) {
			return std::conj(arg);
		}
		else {
			return arg;
		}
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct conj<xsimd::batch<T, A>> {
	xsimd::batch<remove_complex_t<T>, A> operator()(const xsimd::batch<T, A>& arg) const {
		if constexpr (is_complex_v<std::decay_t<T>>) {
			return xsimd::conj(arg);
		}
		else {
			return arg;
		}
	}
};
#endif


template <>
struct conj<void> {
	template <class T>
	constexpr auto operator()(T&& arg) const {
		return conj<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};


template <class T = void>
struct sqrt {
	constexpr auto operator()(const T& arg) const {
		return std::sqrt(arg);
	}
};


#ifdef MATHTER_ENABLE_SIMD
template <class T, class A>
struct sqrt<xsimd::batch<T, A>> {
	xsimd::batch<remove_complex_t<T>, A> operator()(const xsimd::batch<T, A>& arg) const {
		return xsimd::sqrt(arg);
	}
};
#endif


template <>
struct sqrt<void> {
	template <class T>
	constexpr auto operator()(T&& arg) const {
		return sqrt<std::decay_t<T>>{}(std::forward<T>(arg));
	}
};

} // namespace mathter