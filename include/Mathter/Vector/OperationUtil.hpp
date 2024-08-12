// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/TypeTraits.hpp"
#include "SIMDUtil.hpp"
#include "Vector.hpp"

#include <cstddef>
#include <type_traits>



namespace mathter {

template <int NumNonMasked, class Batch, class Element>
Batch FillMasked(Batch batch, Element value) {
#ifdef MATHTER_ENABLE_SIMD
	struct MaskGenerator {
		static constexpr bool get(unsigned idx, [[maybe_unused]] unsigned size) noexcept {
			return idx < NumNonMasked;
		}
	};
	if constexpr (NumNonMasked == Batch::size) {
		return batch;
	}
	else {
		const auto fillers = Batch{ value };
		const auto mask = xsimd::make_batch_bool_constant<Batch, MaskGenerator>();
		return xsimd::select(mask, batch, fillers);
	}
#else
	return batch;
#endif
}


template <int Dim, bool Packed, class Fun, class... Args>
struct is_batched_invocable {
	using ScalarResult = std::invoke_result_t<Fun, Args...>;
	using BatchResult = MakeBatch<ScalarResult, Dim, Packed>;

	template <class... Batches>
	static constexpr auto Get(std::nullptr_t)
		-> std::enable_if_t<std::is_invocable_v<Fun, Batches...>, bool> {
		return false;
	}

	template <class...>
	static constexpr auto Get(...) {
		return false;
	}

	static constexpr bool value = Get<MakeBatch<Args, Dim, Packed>...>(nullptr);
};


template <int Dim, bool Packed, class Fun, class... Args>
inline constexpr auto is_batched_invocable_v = is_batched_invocable<Dim, Packed, Fun, Args...>::value;


template <class T, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoUnaryOpScalar(const Vector<T, Dim, Packed>& arg, Fun&& fun, std::index_sequence<Indices...>) {
	using R = std::invoke_result_t<Fun, T>;
	return Vector<R, Dim, Packed>{ fun(arg[Indices])... };
}


template <class T, int Dim, bool Packed, class Fun>
auto DoUnaryOp(const Vector<T, Dim, Packed>& arg, Fun&& fun) {
	using R = std::invoke_result_t<Fun, T>;
	if constexpr (is_batched_invocable_v<Dim, Packed, Fun, T>) {
		const auto batch = fun(arg.elements.Load());
		return Vector<R, Dim, Packed>(batch);
	}
	else {
		return DoUnaryOpScalar(arg, fun, std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2, class Fun, size_t... Indices>
auto DoBinaryOpScalar(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs, Fun&& fun, std::index_sequence<Indices...>) {
	using T = std::invoke_result_t<Fun, T1, T2>;
	constexpr auto Packed = Packed1 && Packed2;
	return Vector<T, Dim, Packed>{ T(fun(lhs[Indices], rhs[Indices]))... };
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2, class Fun>
auto DoBinaryOp(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs, Fun&& fun) {
	using T = std::invoke_result_t<Fun, T1, T2>;
	constexpr auto Packed = Packed1 && Packed2;
	if constexpr (is_batched_invocable_v<Dim, Packed, Fun, T1, T2>) {
		return Vector<T, Dim, Packed>{ fun(lhs.elements.Load(), rhs.elements.Load()) };
	}
	else {
		return DoBinaryOpScalar(lhs, rhs, fun, std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, class T3, int Dim, bool Packed1, bool Packed2, bool Packed3, class Fun, size_t... Indices>
auto DoTernaryOpScalar(const Vector<T1, Dim, Packed1>& a, const Vector<T2, Dim, Packed2>& b, const Vector<T3, Dim, Packed3>& c, Fun&& fun, std::index_sequence<Indices...>) {
	using T = std::invoke_result_t<Fun, T1, T2, T3>;
	constexpr auto Packed = Packed1 && Packed2 && Packed3;
	return Vector<T, Dim, Packed>{ T(fun(a[Indices], b[Indices], c[Indices]))... };
}


template <class T1, class T2, class T3, int Dim, bool Packed1, bool Packed2, bool Packed3, class Fun>
auto DoTernaryOp(const Vector<T1, Dim, Packed1>& a, const Vector<T2, Dim, Packed2>& b, const Vector<T3, Dim, Packed3>& c, Fun&& fun) {
	using T = std::invoke_result_t<Fun, T1, T2, T3>;
	constexpr auto Packed = Packed1 && Packed2 && Packed3;
	if constexpr (is_batched_invocable_v<Dim, Packed, Fun, T1, T2, T3>) {
		return Vector<T, Dim, Packed>{ fun(a.elements.Load(), b.elements.Load(), c.elements.Load()) };
	}
	else {
		return DoTernaryOpScalar(a, b, c, fun, std::make_index_sequence<Dim>{});
	}
}

} // namespace mathter