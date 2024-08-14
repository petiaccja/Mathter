#include "Benchmark.hpp"

#include <variant>


template <class Op>
struct GenericBinaryFixture {
	template <class Lhs, class Rhs, size_t Count>
	MATHTER_FORCEINLINE auto Latency(const Lhs& lhs, const std::array<Rhs, Count>& rhs) const {
		return std::tuple(DependentLoop(op, lhs, rhs), Count);
	}

	template <class Lhs, class Rhs, size_t Lanes, size_t Count>
	MATHTER_FORCEINLINE auto Throughput(const std::array<Lhs, Lanes>& lhs, const std::array<Rhs, Count>& rhs) const {
		return std::tuple(IndependentLoop(op, lhs, rhs), Count);
	}

	Op op;
};

template <class Op>
GenericBinaryFixture(const Op&) -> GenericBinaryFixture<Op>;


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