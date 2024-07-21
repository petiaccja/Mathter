// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/TypeTraits.hpp"
#include "Vector.hpp"



namespace mathter {

template <int NumNonMasked, class Batch, class Element>
Batch FillMasked(Batch batch, Element value) {
#ifdef MATHTER_ENABLE_SIMD
	struct MaskGenerator {
		static constexpr bool get(unsigned idx, unsigned size) noexcept {
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


template <class T, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoUnaryOpScalar(const Vector<T, Dim, Packed>& arg, Fun&& fun, std::index_sequence<Indices...>) {
	using R = std::invoke_result_t<Fun, T>;
	return Vector<R, Dim, Packed>{ fun(arg[Indices])... };
}


template <class T, int Dim, bool Packed, class Fun>
auto DoUnaryOp(const Vector<T, Dim, Packed>& arg, Fun&& fun) {
	using R = std::invoke_result_t<Fun, T>;
	using BatchT = typename Vector<T, Dim, Packed>::Batch;
	if constexpr (arg.isBatched && std::is_invocable_v<Fun, BatchT>) {
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
	using BatchT1 = typename Vector<T1, Dim, Packed1>::Batch;
	using BatchT2 = typename Vector<T2, Dim, Packed2>::Batch;
	if constexpr (IsBatched<T, Dim, Packed1>() && std::is_invocable_v<Fun, BatchT1, BatchT2>) {
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
	using BatchT1 = typename Vector<T1, Dim, Packed>::Batch;
	using BatchT2 = typename Vector<T2, Dim, Packed>::Batch;
	using BatchT3 = typename Vector<T3, Dim, Packed>::Batch;
	if constexpr (IsBatched<T, Dim, Packed>() && std::is_invocable_v<Fun, BatchT1, BatchT2, BatchT3>) {
		return Vector<T, Dim, Packed>{ fun(a.elements.Load(), b.elements.Load(), c.elements.Load()) };
	}
	else {
		return DoTernaryOpScalar(a, b, c, fun, std::make_index_sequence<Dim>{});
	}
}

} // namespace mathter