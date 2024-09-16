// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/TypeTraits.hpp"
#include "../Common/Types.hpp"
#include "SIMDUtil.hpp"

#include <array>
#include <cassert>
#include <type_traits>


namespace mathter {

/// <summary> Enables element swizzling (reordering elements) for vectors. </summary>
/// <remarks>
/// This class is not supposed to be used on its own. Use the dedicated xx, xy, etc.
///	members of Vector.
/// </remarks>
template <class T, int Dim, bool Packed, int... Indices>
struct Swizzle {
	using SourceStorage = Storage<T, Dim, Packed>;
	using TargetStorage = Storage<T, sizeof...(Indices), Packed>;

	static constexpr auto isSourceBatched = IsBatched<T, Dim, Packed>();
	static constexpr auto isTargetBatched = IsBatched<T, sizeof...(Indices), Packed>();
	using SourceBatch = MakeBatch<T, Dim, Packed>;
	using TargetBatch = MakeBatch<T, sizeof...(Indices), Packed>;

	alignas(GetStorageAlignment<T, Dim, Packed>()) SourceStorage array;


	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy[2] returns y. </summary>
	T& operator[](size_t idx);

	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy[2] returns y. </summary>
	const T& operator[](size_t idx) const;

	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy(2) returns y. </summary>
	T& operator()(size_t idx);

	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy(2) returns y. </summary>
	const T& operator()(size_t idx) const;

	/// <summary> Assigns the swizzled elements of the other to the un-swizzled slots of this. </summary>
	template <class TOther, int DimOther, bool PackedOther, int... IndicesOther>
	Swizzle& operator=(const Swizzle<TOther, DimOther, PackedOther, IndicesOther...>& rhs);

	/// <summary> Assigns the elements of a vector to the un-swizzled slots of this. </summary>
	template <class TOther, int DimOther, bool PackedOther>
	Swizzle& operator=(const Vector<TOther, DimOther, PackedOther>& rhs);

	/// <summary> Assign a scalar if the swizzle is only a single element of the vector. </summary>
	template <class TOther, std::enable_if_t<is_scalar_v<TOther> && sizeof...(Indices) == 1, int> = 0>
	Swizzle& operator=(const TOther& value);


	template <class T2, class = std::enable_if_t<std::is_convertible_v<T, T2> && sizeof...(Indices) == 1, T>>
	operator T2() const;


	/// <summary> Do not use. Used in Vector internally.</summary>
	template <class ContexprIfHack = T>
	TargetStorage Linearize() const;

private:
	static constexpr size_t GetSourceIndex(size_t idx) {
		constexpr std::array table = { static_cast<size_t>(Indices)... };
		return table[idx];
	}

	template <class ContexprIfHack = T>
	static SourceStorage Delinearize(const TargetStorage& target);
	template <class ContexprIfHack = T>
	static SourceStorage Blend(const SourceStorage& old, const SourceStorage& fresh);
};


namespace check_standard_layout {

	// This is extremely important to allow safe access to unions of Swizzles.
	static_assert(std::is_standard_layout_v<Swizzle<float, 3, false, 1, 2>>);

} // namespace check_standard_layout


namespace impl {

	template <auto... Indices>
	struct LinearizationGenerator {
		static constexpr unsigned get(unsigned idx, unsigned size) {
			constexpr std::array<unsigned, sizeof...(Indices)> table{ static_cast<unsigned>(Indices)... };
			return idx < table.size() ? table[idx] : 0u;
		}
	};

	template <auto... Indices>
	struct DelinearizationGenerator {
		static constexpr unsigned get_occurence(unsigned value, std::array<unsigned, sizeof...(Indices)> values, unsigned location = 0) {
			if (values.size() <= location) {
				return 0;
			}
			return value == values[location] ? location : get_occurence(value, values, location + 1);
		}

		static constexpr unsigned get(unsigned idx, unsigned size) {
			constexpr std::array<unsigned, sizeof...(Indices)> table{ static_cast<unsigned>(Indices)... };
			return get_occurence(idx, table);
		}
	};

	template <auto... Indices>
	struct BlendGenerator {
		static constexpr unsigned sentinel = 1000'000;

		static constexpr unsigned get_occurence(unsigned value, std::array<unsigned, sizeof...(Indices)> values, unsigned location = 0) {
			if (values.size() <= location) {
				return sentinel;
			}
			return value == values[location] ? location : get_occurence(value, values, location + 1);
		}

		static constexpr unsigned get(unsigned idx, unsigned size) {
			constexpr std::array<unsigned, sizeof...(Indices)> table{ static_cast<unsigned>(Indices)... };
			return get_occurence(idx, table) != sentinel;
		}
	};

} // namespace impl


template <class T, int Dim, bool Packed, int... Indices>
T& Swizzle<T, Dim, Packed, Indices...>::operator[](size_t idx) {
	assert(idx < Dim);
	return array[GetSourceIndex(idx)];
}


template <class T, int Dim, bool Packed, int... Indices>
const T& Swizzle<T, Dim, Packed, Indices...>::operator[](size_t idx) const {
	assert(idx < Dim);
	return array[GetSourceIndex(idx)];
}


template <class T, int Dim, bool Packed, int... Indices>
T& Swizzle<T, Dim, Packed, Indices...>::operator()(size_t idx) {
	assert(idx < Dim);
	return array[GetSourceIndex(idx)];
}


template <class T, int Dim, bool Packed, int... Indices>
const T& Swizzle<T, Dim, Packed, Indices...>::operator()(size_t idx) const {
	assert(idx < Dim);
	return array[GetSourceIndex(idx)];
}


template <class T, int Dim, bool Packed, int... Indices>
template <class TOther, int DimOther, bool PackedOther, int... IndicesOther>
Swizzle<T, Dim, Packed, Indices...>& Swizzle<T, Dim, Packed, Indices...>::operator=(const Swizzle<TOther, DimOther, PackedOther, IndicesOther...>& rhs) {
	static_assert(std::is_convertible_v<TOther, T> && sizeof...(Indices) == sizeof...(IndicesOther), "incompatible swizzle operators");

	const auto rhsLin = rhs.Linearize();
	auto lhsLin = TargetStorage{};
	std::copy(rhsLin.begin(), rhsLin.end(), lhsLin.begin());
	array = Blend(array, Delinearize(lhsLin));
	return *this;
}


template <class T, int Dim, bool Packed, int... Indices>
template <class TOther, std::enable_if_t<is_scalar_v<TOther> && sizeof...(Indices) == 1, int>>
Swizzle<T, Dim, Packed, Indices...>& Swizzle<T, Dim, Packed, Indices...>::operator=(const TOther& value) {
	static_assert(sizeof...(Indices) == 1);
	(..., (array[Indices] = value)); // This expands to a single assignment if sizeof...(Indices) == 1.
	return *this;
}


template <class T, int Dim, bool Packed, int... Indices>
template <class T2, class>
Swizzle<T, Dim, Packed, Indices...>::operator T2() const {
	return static_cast<T2>((*this)[0]);
}


template <class T, int Dim, bool Packed, int... Indices>
template <class>
typename Swizzle<T, Dim, Packed, Indices...>::TargetStorage Swizzle<T, Dim, Packed, Indices...>::Linearize() const {
#if MATHTER_ENABLE_SIMD
	if constexpr (isSourceBatched && sizeof...(Indices) <= Dim) {
		// Lambda trickery: if constexpr evaluated unless template, lambda is template.
		using Integer = make_sized_integer_t<sizeof(remove_complex_t<T>), false>;
		using SourceIntBatch = xsimd::batch<Integer, typename SourceBatch::arch_type>;

		const auto sourceBatch = SourceBatch::load_aligned(array.data());
		const auto mask = xsimd::make_batch_constant<typename SourceIntBatch::value_type, typename SourceIntBatch::arch_type, impl::LinearizationGenerator<Indices...>>();
		const auto linearizedBatch = xsimd::swizzle(sourceBatch, mask);

		alignas(GetStorageAlignment<T, Dim, Packed>()) SourceStorage linearized;
		linearizedBatch.store_aligned(linearized.data());

		TargetStorage result;
		std::copy(linearized.begin(), linearized.begin() + sizeof...(Indices), result.begin());
		return result;
	}
#endif
	return TargetStorage{ array[Indices]... };
}


template <class T, int Dim, bool Packed, int... Indices>
template <class>
auto Swizzle<T, Dim, Packed, Indices...>::Delinearize(const TargetStorage& target) -> SourceStorage {
#if MATHTER_ENABLE_SIMD
	if constexpr (isSourceBatched && sizeof...(Indices) <= Dim) {
		using Integer = make_sized_integer_t<sizeof(remove_complex_t<T>), false>;
		using SourceIntBatch = xsimd::batch<Integer, typename SourceBatch::arch_type>;

		alignas(GetStorageAlignment<T, Dim, Packed>()) SourceStorage targetPadded;
		std::copy(target.begin(), target.end(), targetPadded.begin());
		const auto targetBatch = SourceBatch::load_aligned(targetPadded.data());
		const auto mask = xsimd::make_batch_constant<typename SourceIntBatch::value_type, typename SourceIntBatch::arch_type, impl::DelinearizationGenerator<Indices...>>();
		const auto delinearizedBatch = xsimd::swizzle(targetBatch, mask);

		alignas(GetStorageAlignment<T, Dim, Packed>()) SourceStorage delinearized;
		delinearizedBatch.store_aligned(delinearized.data());

		return delinearized;
	}
#endif
	SourceStorage delinearized;
	for (size_t i = 0; i < target.size(); ++i) {
		delinearized[GetSourceIndex(int(i))] = target[i];
	}
	return delinearized;
}


template <class T, int Dim, bool Packed, int... Indices>
template <class>
auto Swizzle<T, Dim, Packed, Indices...>::Blend(const SourceStorage& old, const SourceStorage& fresh) -> SourceStorage {
#if MATHTER_ENABLE_SIMD
	if constexpr (isSourceBatched && sizeof...(Indices) <= Dim) {
		const auto oldBatch = SourceBatch::load_unaligned(old.data());
		const auto freshBatch = SourceBatch::load_unaligned(fresh.data());
		const auto mask = xsimd::make_batch_bool_constant<typename SourceBatch::value_type, typename SourceBatch::arch_type, impl::BlendGenerator<Indices...>>();
		const auto blendedBatch = xsimd::select(mask, freshBatch, oldBatch);

		alignas(GetStorageAlignment<T, Dim, Packed>()) SourceStorage blended;
		blendedBatch.store_aligned(blended.data());

		return blended;
	}
#endif
	SourceStorage blended;
	for (size_t i = 0; i < old.size(); ++i) {
		blended[i] = impl::BlendGenerator<Indices...>().get(i, old.size()) ? fresh[i] : old[i];
	}
	return blended;
}

} // namespace mathter


#include "Vector.hpp"


namespace mathter {

template <class T, int Dim, bool Packed, int... Indices>
template <class TOther, int DimOther, bool PackedOther>
Swizzle<T, Dim, Packed, Indices...>& Swizzle<T, Dim, Packed, Indices...>::operator=(const Vector<TOther, DimOther, PackedOther>& rhs) {
	static_assert(std::is_convertible_v<TOther, T> && sizeof...(Indices) == DimOther, "incompatible vector");

	auto lhsLin = TargetStorage{};
	std::copy(rhs.begin(), rhs.end(), lhsLin.begin());
	array = Blend(array, Delinearize(lhsLin));
	return *this;
}

} // namespace mathter