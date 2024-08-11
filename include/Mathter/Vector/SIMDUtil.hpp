#pragma once


#if MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif

#include <array>


namespace mathter {

constexpr int GetBatchSize(int Dim, bool Packed) {
	return Packed	? Dim :
		   Dim == 3 ? 4 :
		   Dim == 5 ? 8 :
		   Dim == 6 ? 8 :
		   Dim == 7 ? 8 :
					  Dim;
}


template <class T, int Dim, bool Packed>
constexpr int IsBatched() {
#if MATHTER_ENABLE_SIMD
	if constexpr (!Packed) {
		using BatchT = typename xsimd::make_sized_batch<T, GetBatchSize(Dim, Packed)>::type;
		return !std::is_void_v<BatchT>;
	}
#endif
	return false;
}


namespace impl {
	template <class T, int Dim, bool Packed, bool>
	struct MakeBatchHelper {
		using type = void;
	};

#if MATHTER_ENABLE_SIMD
	template <class T, int Dim, bool Packed>
	struct MakeBatchHelper<T, Dim, Packed, true> {
		using type = typename xsimd::make_sized_batch<T, GetBatchSize(Dim, Packed)>::type;
	};
#endif

} // namespace impl


template <class T, int Dim, bool Packed>
using MakeBatch = typename impl::MakeBatchHelper<T, Dim, Packed, IsBatched<T, Dim, Packed>()>::type;


template <class T, int Dim, bool Packed>
constexpr int GetStorageSize() {
	return IsBatched<T, Dim, Packed>() ? GetBatchSize(Dim, Packed) : Dim;
}


template <class T, int Dim, bool Packed>
constexpr int GetStorageAlignment() {
	if constexpr (IsBatched<T, Dim, Packed>()) {
		using B = MakeBatch<T, Dim, Packed>;
		static_assert(!std::is_void_v<B>, "IsBatched should prevent this case from ever happening.");
		return alignof(B);
	}
	return alignof(T);
}


template <class T, int Dim, bool Packed>
using Storage = std::array<T, GetStorageSize<T, Dim, Packed>()>;

static_assert(std::is_standard_layout_v<Storage<float, 3, false>>);



} // namespace mathter