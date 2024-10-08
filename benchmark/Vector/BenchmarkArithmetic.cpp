#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "../Input.hpp"

#include <Mathter/Transforms/RandomBuilder.hpp>
#include <Mathter/Vector.hpp>

#include <random>


using namespace mathter;

namespace {

template <class Vec, int Count>
std::array<Vec, Count> MakeInput() {
	using Scalar = scalar_type_t<Vec>;
	using Real = remove_complex_t<Scalar>;
	using namespace std::complex_literals;

	std::mt19937_64 rne;
	std::uniform_real_distribution<Real> rng(Real(-0.5), Real(0.5));

	std::array<Vec, Count> r;
	for (auto& v : r) {
		if constexpr (is_complex_v<Scalar>) {
			const Vec real = Random(rng, rne);
			const Vec imag = Random(rng, rne);
			v = real + Scalar(1if) * imag;
		}
		else {
			v = Random(rng, rne);
		}
	}
	return r;
};

#define VECTOR_BINOP_BENCHMARK_CASE(TYPE, DIM, PACKED, OP, OPTEXT)                  \
	BENCHMARK_CASE(#TYPE "." #DIM " " OPTEXT " " #TYPE "." #DIM " (P=" #PACKED ")", \
				   "[Vector][Arithmetic]",                                          \
				   50,                                                              \
				   64,                                                              \
				   GenericNAryFixture{ OP{} },                                      \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 1>()[0],              \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 4>(),                 \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 64>());


VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, std::plus<>, "+");

VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, std::multiplies<>, "*");

VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, std::divides<>, "/");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, std::plus<>, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, std::plus<>, "+");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, std::multiplies<>, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, std::multiplies<>, "*");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, std::divides<>, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, std::divides<>, "/");

} // namespace