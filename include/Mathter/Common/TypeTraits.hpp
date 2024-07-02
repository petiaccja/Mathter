#pragma once

#include "TypeTraits.hpp"

#include <type_traits>


namespace mathter {

//------------------------------------------------------------------------------
// Common arithmetic type
//------------------------------------------------------------------------------

namespace impl {

	template <class... T>
	struct common_arithmetic_type {
		template <class... U>
		static constexpr auto eval(std::nullptr_t) -> std::enable_if<true, decltype((... + std::declval<U>()))> {
			return {};
		}

		template <class... U>
		static constexpr auto eval(const void*) -> std::enable_if<false> {
			return {};
		}

		using enable_if_type = decltype(eval<T...>(nullptr));
	};

} // namespace impl


template <class... T>
struct common_arithmetic_type : impl::common_arithmetic_type<T...>::enable_if_type {};


template <class... T>
using common_arithmetic_type_t = typename common_arithmetic_type<T...>::type;

// Tests
static_assert(std::is_same_v<common_arithmetic_type_t<float, double>, double>);
static_assert(std::is_same_v<common_arithmetic_type_t<double, float>, double>);


} // namespace mathter