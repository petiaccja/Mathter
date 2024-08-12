#pragma once

#include <cstddef>
#include <utility>


namespace mathter {


#ifdef _MSC_VER
#define MATHTER_NOINLINE __declspec(noinline)
#define MATHTER_FORCEINLINE __forceinline
#define MATHTER_FLATTEN [[msvc::flatten]]
#else
#define MATHTER_NOINLINE __attribute__((noinline))
#define MATHTER_FORCEINLINE __attribute__((always_inline))
#define MATHTER_FLATTEN __attribute__((flatten))
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


template <ptrdiff_t First, ptrdiff_t Last, ptrdiff_t Step, ptrdiff_t Limit, class Fun, class... Args>
void ForUnrolled(std::integral_constant<ptrdiff_t, First>,
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
auto ForUnrolled(Fun&& fun, Args&&... args) -> std::enable_if_t<!std::is_convertible_v<Fun, std::integral_constant<ptrdiff_t, First>>> {
	ForUnrolled(std::integral_constant<ptrdiff_t, First>{},
				std::integral_constant<ptrdiff_t, Last>{},
				std::integral_constant<ptrdiff_t, Step>{},
				std::integral_constant<ptrdiff_t, Limit>{},
				std::forward<Fun>(fun),
				std::forward<Args>(args)...);
}

} // namespace mathter