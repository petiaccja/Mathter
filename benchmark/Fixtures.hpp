#include "Benchmark.hpp"

#include <variant>


template <class Op>
struct GenericNAryFixture {
	template <class Lhs, class... Args, size_t Count>
	MATHTER_FORCEINLINE auto Latency(const Lhs& lhs, const std::array<Args, Count>&... args) const {
		return std::tuple(DependentLoop(op, lhs, args...), Count);
	}

	template <class Lhs, class... Args, size_t Lanes, size_t Count>
	MATHTER_FORCEINLINE auto Throughput(const std::array<Lhs, Lanes>& lhs, const std::array<Args, Count>&... args) const {
		return std::tuple(IndependentLoop(op, lhs, args...), Count);
	}

	Op op;
};

template <class Op>
GenericNAryFixture(const Op&) -> GenericNAryFixture<Op>;


template <class Op>
struct GenericUnaryFixture {
	template <class Lhs, size_t Count>
	MATHTER_FORCEINLINE auto Latency(const Lhs& arg, const std::array<std::monostate, Count>& counter) const {
		return std::tuple(DependentLoop(op, arg, counter), Count);
	}

	template <class Lhs, size_t Lanes, size_t Count>
	MATHTER_FORCEINLINE auto Throughput(const std::array<Lhs, Lanes>& arg, const std::array<std::monostate, Count>& counter) const {
		return std::tuple(IndependentLoop(op, arg, counter), Count);
	}

	Op op;
};

template <class Op>
GenericUnaryFixture(const Op&) -> GenericUnaryFixture<Op>;