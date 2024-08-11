#pragma once

#include <cstddef>
#include <utility>


namespace mathter {

namespace impl {

	template <class Func, size_t... Indices>
	auto LoopUnroll(Func&& func, std::index_sequence<Indices...>) {
		return func(Indices...);
	}

} // namespace impl


template <size_t Iterations, class Func>
auto LoopUnroll(Func&& func) {
	return impl::LoopUnroll(std::forward<Func>(func), std::make_index_sequence<Iterations>{});
}

} // namespace mathter