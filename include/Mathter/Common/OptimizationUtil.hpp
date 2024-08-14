#pragma once

#include <cstddef>
#include <functional>
#include <utility>


namespace mathter {


#ifdef _MSC_VER
#define MATHTER_NOINLINE [[msvc::noinline]]
#define MATHTER_FORCEINLINE [[msvc::forceinline]]
#define MATHTER_FLATTEN_BLOCK [[msvc::flatten]]
#else
#define MATHTER_NOINLINE [[gnu::noinline]]
#define MATHTER_FORCEINLINE [[gnu::always_inline]]
#define MATHTER_FLATTEN_FUNCTION [[gnu::flatten]]
#endif


namespace impl {

	template <size_t... Indices, class Func, class... Args>
	MATHTER_FORCEINLINE auto LoopUnrollHelper(std::index_sequence<Indices...>, Func&& func, Args&&... args) {
		return func(std::forward<Args>(args)..., Indices...);
	}

} // namespace impl


template <size_t Iterations, class Func, class... Args>
MATHTER_FORCEINLINE auto LoopUnroll(Func&& func, Args&&... args) {
	return impl::LoopUnrollHelper(std::make_index_sequence<Iterations>{}, std::forward<Func>(func), std::forward<Args>(args)...);
}


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
	constexpr auto Count = Last - First / Step;
	if constexpr (Count <= Limit) {
		if constexpr (0 < Step ? First < Last : First > Last) {
			fun(First, args...);
			ForUnrolled(std::integral_constant<ptrdiff_t, First + Step>{}, last, step, limit, std::move(fun), std::forward<Args>(args)...);
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