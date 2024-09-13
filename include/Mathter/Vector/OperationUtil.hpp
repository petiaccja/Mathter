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
		const auto mask = xsimd::make_batch_bool_constant<typename Batch::value_type, typename Batch::arch_type, MaskGenerator>();
		return xsimd::select(mask, batch, fillers);
	}
#else
	return batch;
#endif
}


template <int NumNonMasked, class Batch>
Batch FillMaskedWithFirst(Batch batch) {
#ifdef MATHTER_ENABLE_SIMD
	using Scalar = xsimd::scalar_type_t<Batch>;
	using Real = remove_complex_t<Scalar>;
	using Uint = make_sized_integer_t<sizeof(Real), false>;
	struct Generator {
		static constexpr Uint get(unsigned idx, [[maybe_unused]] unsigned size) noexcept {
			return idx < NumNonMasked ? Uint(idx) : Uint(0);
		}
	};

	using UintBatch = xsimd::batch<Uint, typename Batch::arch_type>;

	const auto mask = xsimd::make_batch_constant<typename UintBatch::value_type, typename UintBatch::arch_type, Generator>();
	return xsimd::swizzle(batch, mask);
#else
	return batch;
#endif
}


template <class Fun, class... Vectors>
struct is_invocable_batched {
	using Scalar = std::invoke_result_t<Fun, scalar_type_t<Vectors>...>;
	static constexpr auto Dim = (... + dimension_v<Vectors>) / sizeof...(Vectors);
	static constexpr auto Packed = (... && is_packed_v<Vectors>);
	using Batch = typename Vector<Scalar, Dim, Packed>::Batch;

	template <class... Batches>
	static constexpr auto Get(std::nullptr_t)
		-> std::enable_if_t<(... && std::is_convertible_v<Batches, Batch>)
								&& std::is_invocable_v<Fun, std::conditional_t<false, Batches, Batch>...>,
							bool> {
		return true;
	}

	template <class...>
	static constexpr auto Get(...) {
		return false;
	}

	static constexpr bool value = Get<typename Vectors::Batch...>(nullptr);
};


template <class Fun, class... Vectors>
inline constexpr auto is_invocable_batched_v = is_invocable_batched<Fun, Vectors...>::value;


#if defined(MATHTER_ENABLE_SIMD) && _M_AMD64
static_assert(is_invocable_batched_v<std::plus<>, Vector<float, 4, false>, Vector<float, 4, false>>);
#endif


template <class Fun, class... Batches>
auto invoke_batched(Fun&& fun, Batches&&... batches) {
	using R = std::invoke_result_t<Fun, typename std::decay_t<Batches>::value_type...>;
	const auto Dim = (... + std::decay_t<Batches>::size) / sizeof...(Batches);
	using Batch = MakeBatch<R, Dim, false>;
	return std::invoke(std::forward<Fun>(fun), Batch(std::forward<Batches>(batches))...);
}


template <class T, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoUnaryOpScalar(const Vector<T, Dim, Packed>& arg, Fun&& fun, std::index_sequence<Indices...>) {
	using R = std::invoke_result_t<Fun, T>;
	return Vector<R, Dim, Packed>{ fun(arg[Indices])... };
}


template <class T, int Dim, bool Packed, class Fun>
auto DoUnaryOp(const Vector<T, Dim, Packed>& arg, Fun&& fun) {
	using R = std::invoke_result_t<Fun, T>;
	if constexpr (is_invocable_batched_v<Fun, std::decay_t<decltype(arg)>>) {
		return Vector<R, Dim, Packed>{ invoke_batched(std::forward<Fun>(fun), arg.elements.Load()) };
	}
	else {
		return DoUnaryOpScalar(arg, std::forward<Fun>(fun), std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2, class Fun, size_t... Indices>
auto DoBinaryOpScalar(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs, Fun&& fun, std::index_sequence<Indices...>) {
	using T = std::invoke_result_t<Fun, T1, T2>;
	constexpr auto Packed = Packed1 && Packed2;
	return Vector<T, Dim, Packed>{ fun(lhs[Indices], rhs[Indices])... };
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2, class Fun>
auto DoBinaryOp(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs, Fun&& fun) {
	using T = std::invoke_result_t<Fun, T1, T2>;
	constexpr auto Packed = Packed1 && Packed2;
	if constexpr (is_invocable_batched_v<Fun, std::decay_t<decltype(lhs)>, std::decay_t<decltype(rhs)>>) {
		return Vector<T, Dim, Packed>{ invoke_batched(fun, lhs.elements.Load(), rhs.elements.Load()) };
	}
	else {
		return DoBinaryOpScalar(lhs, rhs, std::forward<Fun>(fun), std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, class T3, int Dim, bool Packed1, bool Packed2, bool Packed3, class Fun, size_t... Indices>
auto DoTernaryOpScalar(const Vector<T1, Dim, Packed1>& a, const Vector<T2, Dim, Packed2>& b, const Vector<T3, Dim, Packed3>& c, Fun&& fun, std::index_sequence<Indices...>) {
	using T = std::invoke_result_t<Fun, T1, T2, T3>;
	constexpr auto Packed = Packed1 && Packed2 && Packed3;
	return Vector<T, Dim, Packed>{ fun(a[Indices], b[Indices], c[Indices])... };
}


template <class T1, class T2, class T3, int Dim, bool Packed1, bool Packed2, bool Packed3, class Fun>
auto DoTernaryOp(const Vector<T1, Dim, Packed1>& a, const Vector<T2, Dim, Packed2>& b, const Vector<T3, Dim, Packed3>& c, Fun&& fun) {
	using T = std::invoke_result_t<Fun, T1, T2, T3>;
	constexpr auto Packed = Packed1 && Packed2 && Packed3;
	if constexpr (is_invocable_batched_v<Fun, std::decay_t<decltype(a)>, std::decay_t<decltype(b)>, std::decay_t<decltype(c)>>) {
		return Vector<T, Dim, Packed>{ invoke_batched(std::forward<Fun>(fun), a.elements.Load(), b.elements.Load(), c.elements.Load()) };
	}
	else {
		return DoTernaryOpScalar(a, b, c, std::forward<Fun>(fun), std::make_index_sequence<Dim>{});
	}
}

} // namespace mathter