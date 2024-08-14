#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "Input.hpp"

#include <Mathter/Common/TypeTraits.hpp>


using namespace mathter;

namespace {

#define SCALAR_BINOP_BENCHMARK_CASE(TYPE, OP, OPTEXT)      \
	BENCHMARK_CASE(#TYPE " " OPTEXT " " #TYPE,             \
				   "[Scalar][Arithmetic]",                 \
				   50,                                     \
				   64,                                     \
				   GenericBinaryFixture{ OP },             \
				   MakeConstantInput<TYPE, 1>(TYPE(1))[0], \
				   MakeConstantInput<TYPE, 16>(TYPE(1)),    \
				   MakeConstantInput<TYPE, 256>(TYPE(1)));


SCALAR_BINOP_BENCHMARK_CASE(int32_t, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(int32_t, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(int32_t, std::divides<>{}, "/");

SCALAR_BINOP_BENCHMARK_CASE(int64_t, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(int64_t, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(int64_t, std::divides<>{}, "/");

SCALAR_BINOP_BENCHMARK_CASE(float, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(float, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(float, std::divides<>{}, "/");

SCALAR_BINOP_BENCHMARK_CASE(double, std::multiplies<>{}, "*");
SCALAR_BINOP_BENCHMARK_CASE(double, std::plus<>{}, "+");
SCALAR_BINOP_BENCHMARK_CASE(double, std::divides<>{}, "/");


} // namespace