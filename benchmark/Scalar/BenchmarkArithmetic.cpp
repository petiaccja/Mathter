#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "Input.hpp"

#include <Mathter/Common/TypeTraits.hpp>

#include <complex>
#include <random>


using namespace mathter;

namespace {

#define SCALAR_BINOP_BENCHMARK_CASE(TYPE, OP, OPTEXT) \
	BENCHMARK_CASE(#TYPE " " OPTEXT " " #TYPE,        \
				   "[Scalar][Arithmetic]",            \
				   50,                                \
				   64,                                \
				   GenericBinaryFixture{ OP },        \
				   MakeInput<TYPE, 1>()[0],           \
				   MakeInput<TYPE, 8>(),              \
				   MakeInput<TYPE, 128>());


SCALAR_BINOP_BENCHMARK_CASE(float, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(float, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(float, std::divides<>{}, "/");

SCALAR_BINOP_BENCHMARK_CASE(double, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(double, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(double, std::divides<>{}, "/");


} // namespace