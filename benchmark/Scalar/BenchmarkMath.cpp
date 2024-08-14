#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "Input.hpp"

#include <Mathter/Common/TypeTraits.hpp>


using namespace mathter;

#define WRAP_FUNC(FUNC) [](const auto& arg, auto&&...) { return FUNC((arg)); }

namespace {

#define SCALAR_MATH_BENCHMARK_CASE_MINMAX(TYPE, OP, OPTEXT, MIN, MAX) \
	BENCHMARK_CASE(#TYPE ": " OPTEXT,                                 \
				   "[Scalar][Math]",                                  \
				   50,                                                \
				   64,                                                \
				   GenericUnaryFixture{ OP },                         \
				   MakeInput<TYPE, 1>(MIN, MAX)[0],                   \
				   MakeInput<TYPE, 8>(MIN, MAX),                      \
				   std::array<std::monostate, 64>());


#define SCALAR_MATH_BENCHMARK_CASE(TYPE, OP, OPTEXT) SCALAR_MATH_BENCHMARK_CASE_MINMAX(TYPE, OP, OPTEXT, 0.1f, 0.9f)


SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::abs), "abs");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::sqrt), "sqrt");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::exp), "exp");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::log), "log");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::sin), "sin");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::cos), "cos");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::tan), "tan");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::asin), "asin");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::acos), "acos");
SCALAR_MATH_BENCHMARK_CASE(float, WRAP_FUNC(std::atan), "atan");


SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::abs), "abs");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::sqrt), "sqrt");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::exp), "exp");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::log), "log");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::sin), "sin");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::cos), "cos");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::tan), "tan");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::asin), "asin");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::acos), "acos");
SCALAR_MATH_BENCHMARK_CASE(double, WRAP_FUNC(std::atan), "atan");

} // namespace