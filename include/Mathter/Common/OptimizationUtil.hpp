#pragma once

#include <cstddef>
#include <functional>
#include <utility>


#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-attributes"
#endif


namespace mathter {


#ifdef _MSC_VER
#define MATHTER_NOINLINE [[msvc::noinline]]
#define MATHTER_FORCEINLINE [[msvc::forceinline]] inline
#define MATHTER_FLATTEN_BLOCK [[msvc::flatten]]
#else
#define MATHTER_NOINLINE [[gnu::noinline]]
#define MATHTER_FORCEINLINE [[gnu::always_inline]] inline
#define MATHTER_FLATTEN_FUNCTION [[gnu::flatten]]
#endif


template <class Fun, size_t... Indices>
MATHTER_FORCEINLINE decltype(auto) Apply(Fun&& fun, std::integer_sequence<size_t, Indices...>) {
	return std::invoke(std::forward<Fun>(fun), Indices...);
}


template <ptrdiff_t First, ptrdiff_t Last, ptrdiff_t Step, ptrdiff_t Limit, class Fun, class... Args>
MATHTER_FORCEINLINE void ForUnrolled(std::integral_constant<ptrdiff_t, First>,
									 std::integral_constant<ptrdiff_t, Last> last,
									 std::integral_constant<ptrdiff_t, Step> step,
									 std::integral_constant<ptrdiff_t, Limit> limit,
									 Fun&& fun,
									 Args&&... args) {
	constexpr auto Count = (Last - First) / Step;
	if constexpr (Count <= Limit) {
		if constexpr (0 < Step ? First < Last : First > Last) {
			fun(First, args...);
			ForUnrolled(std::integral_constant<ptrdiff_t, First + Step>{}, last, step, limit, std::forward<Fun>(fun), std::forward<Args>(args)...);
		}
	}
	else {
		for (auto i = First; i < Last; i += Step) {
			fun(i, args...);
		}
	}
}


template <ptrdiff_t First, ptrdiff_t Last, ptrdiff_t Step, ptrdiff_t Limit, class Fun, class... Args>
MATHTER_FORCEINLINE auto ForUnrolled(Fun&& fun, Args&&... args) -> std::enable_if_t<!std::is_convertible_v<Fun, std::integral_constant<ptrdiff_t, First>>> {
	ForUnrolled(std::integral_constant<ptrdiff_t, First>{},
				std::integral_constant<ptrdiff_t, Last>{},
				std::integral_constant<ptrdiff_t, Step>{},
				std::integral_constant<ptrdiff_t, Limit>{},
				std::forward<Fun>(fun),
				std::forward<Args>(args)...);
}

} // namespace mathter


#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif