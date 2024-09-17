// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Vector.hpp"

#include <type_traits>


namespace mathter {

/// <summary> Concatenates the arguments, and returns the concatenated parts as a vector. </summary>
template <class Left,
		  class Right,
		  std::enable_if_t<is_vector_v<std::decay_t<Left>>
							   || is_vector_v<std::decay_t<Right>>
							   || is_swizzle_v<std::decay_t<Left>>
							   || is_swizzle_v<std::decay_t<Right>>,
						   int> = 0>
[[deprecated("use concatenating ctor with CTAD")]]
auto operator|(Left&& lhs, Right&& rhs) {
	return Vector(std::forward<Left>(lhs), std::forward<Right>(rhs));
}

} // namespace mathter