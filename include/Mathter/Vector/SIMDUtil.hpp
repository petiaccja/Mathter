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
	template <class T, int Dim, bool Packed>
	struct MakeBatchHelper {
		static auto GetType() {
			if constexpr (IsBatched<T, Dim, Packed>()) {
#if MATHTER_ENABLE_SIMD
				using B = typename xsimd::make_sized_batch<T, GetBatchSize(Dim, Packed)>::type;
				return static_cast<B*>(nullptr);
#else
				return static_cast<void*>(nullptr);
#endif
			}
			else {
				return static_cast<void*>(nullptr);
			}
		}
	};
} // namespace impl


template <class T, int Dim, bool Packed>
using MakeBatch = std::decay_t<std::remove_pointer_t<decltype(impl::MakeBatchHelper<T, Dim, Packed>::GetType())>>;


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