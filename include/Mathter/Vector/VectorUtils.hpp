// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/TypeTraits.hpp"
#include "VectorImpl.hpp"



namespace mathter {

template <class T, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoUnaryOpScalar(const Vector<T, Dim, Packed>& arg, Fun&& fun, std::index_sequence<Indices...>) {
	return Vector<T, Dim, Packed>{ T(fun(arg[Indices]))... };
}


template <class T, int Dim, bool Packed, class Fun>
auto DoUnaryOp(const Vector<T, Dim, Packed>& arg, Fun&& fun) {
	using MyBatchT = Batch<T, Dim, Packed>;
	if constexpr (IsBatched<T, Dim, Packed>()) {
		return Vector<T, Dim, Packed>{ fun(MyBatchT::load_unaligned(arg.data())) };
	}
	else {
		return DoUnaryOpScalar(arg, fun, std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoBinaryOpScalar(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs, Fun&& fun, std::index_sequence<Indices...>) {
	using T = common_arithmetic_type_t<T1, T2>;
	return Vector<T, Dim, Packed>{ T(fun(lhs[Indices], rhs[Indices]))... };
}


template <class T1, class T2, int Dim, bool Packed, class Fun>
auto DoBinaryOp(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs, Fun&& fun) {
	using T = common_arithmetic_type_t<T1, T2>;
	using MyBatchT = Batch<T, Dim, Packed>;
	using MyBatchT1 = Batch<T1, Dim, Packed>;
	using MyBatchT2 = Batch<T2, Dim, Packed>;
	if constexpr (IsBatched<T, Dim, Packed>()
				  && std::is_convertible_v<MyBatchT1, MyBatchT>
				  && std::is_convertible_v<MyBatchT2, MyBatchT>) {
		return Vector<T, Dim, Packed>{ fun(MyBatchT1::load_unaligned(lhs.data()), MyBatchT2::load_unaligned(rhs.data())) };
	}
	else {
		return DoBinaryOpScalar(lhs, rhs, fun, std::make_index_sequence<Dim>{});
	}
}


template <class T1, class T2, class T3, int Dim, bool Packed, class Fun, size_t... Indices>
auto DoTernaryOpScalar(const Vector<T1, Dim, Packed>& a, const Vector<T2, Dim, Packed>& b, const Vector<T3, Dim, Packed>& c, Fun&& fun, std::index_sequence<Indices...>) {
	using T = common_arithmetic_type_t<T1, T2, T3>;
	return Vector<T, Dim, Packed>{ T(fun(a[Indices], b[Indices], c[Indices]))... };
}


template <class T1, class T2, class T3, int Dim, bool Packed, class Fun>
auto DoTernaryOp(const Vector<T1, Dim, Packed>& a, const Vector<T2, Dim, Packed>& b, const Vector<T3, Dim, Packed>& c, Fun&& fun) {
	using T = common_arithmetic_type_t<T1, T2>;
	using MyBatchT = Batch<T, Dim, Packed>;
	using MyBatchT1 = Batch<T1, Dim, Packed>;
	using MyBatchT2 = Batch<T2, Dim, Packed>;
	using MyBatchT3 = Batch<T3, Dim, Packed>;
	if constexpr (IsBatched<T, Dim, Packed>()
				  && std::is_convertible_v<MyBatchT1, MyBatchT>
				  && std::is_convertible_v<MyBatchT2, MyBatchT>
				  && std::is_convertible_v<MyBatchT3, MyBatchT>) {
		return Vector<T, Dim, Packed>{ fun(MyBatchT1::load_unaligned(a.data()), MyBatchT2::load_unaligned(b.data()), MyBatchT3::load_unaligned(c.data())) };
	}
	else {
		return DoTernaryOpScalar(a, b, c, fun, std::make_index_sequence<Dim>{});
	}
}


template <int NumNonMasked, class Batch, class Element>
inline Batch FillMasked(Batch batch, Element value) {
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

} // namespace mathter