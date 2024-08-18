#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "../Input.hpp"

#include <Mathter/Vector.hpp>


using namespace mathter;

namespace {

struct SumOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& vec) const {
		const auto result = Sum(vec);
		auto copy = vec;
		copy[0] = -result;
		return copy;
	}
};


struct DotOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) const {
		const auto result = Dot(lhs, rhs);
		auto copy = lhs;
		copy[0] = -result;
		return copy;
	}
};


struct LengthOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& vec) const {
		const auto result = Length(vec);
		auto copy = vec;
		copy[0] = -result;
		return copy;
	}
};


struct LengthPreciseOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& vec) const {
		const auto result = LengthPrecise(vec);
		auto copy = vec;
		copy[0] = -result;
		return copy;
	}
};


struct NormalizehOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& vec) const {
		auto result = Normalize(vec);
		result[0] = static_cast<T>(1);
		return result;
	}
};


struct NormalizePreciseOp {
	template <class T, int Dim, bool Packed>
	auto operator()(const Vector<T, Dim, Packed>& vec) const {
		auto result = NormalizePrecise(vec);
		result[0] = static_cast<T>(1);
		return result;
	}
};


struct CrossOp {
	template <class... Vec>
	auto operator()(const Vec&... vectors) const {
		return Cross(vectors...);
	}
};


#define VECTOR_MATH_BINARY_BENCHMARK_CASE(TYPE, DIM, PACKED, OP, OPTEXT) \
	BENCHMARK_CASE(#TYPE "." #DIM ": " OPTEXT " (P=" #PACKED ")",        \
				   "[Vector][Math]",                                     \
				   50,                                                   \
				   64,                                                   \
				   GenericNAryFixture{ OP },                             \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 1>()[0],   \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 4>(),      \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 64>());


#define VECTOR_MATH_TERNARY_BENCHMARK_CASE(TYPE, DIM, PACKED, OP, OPTEXT) \
	BENCHMARK_CASE(#TYPE "." #DIM ": " OPTEXT " (P=" #PACKED ")",         \
				   "[Vector][Math]",                                      \
				   50,                                                    \
				   64,                                                    \
				   GenericNAryFixture{ OP },                              \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 1>()[0],    \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 4>(),       \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 64>(),      \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 64>());


#define VECTOR_MATH_UNARY_BENCHMARK_CASE(TYPE, DIM, PACKED, OP, OPTEXT) \
	BENCHMARK_CASE(#TYPE "." #DIM ": " OPTEXT " (P=" #PACKED ")",       \
				   "[Vector][Math]",                                    \
				   50,                                                  \
				   64,                                                  \
				   GenericUnaryFixture{ OP },                           \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 1>()[0],  \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 4>(),     \
				   std::array<std::monostate, 64>());


VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, SumOp{}, "Sum");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 3, false, SumOp{}, "Sum");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 4, false, SumOp{}, "Sum");
VECTOR_MATH_BINARY_BENCHMARK_CASE(float, 2, false, DotOp{}, "Dot");
VECTOR_MATH_BINARY_BENCHMARK_CASE(float, 3, false, DotOp{}, "Dot");
VECTOR_MATH_BINARY_BENCHMARK_CASE(float, 4, false, DotOp{}, "Dot");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 3, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 4, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 3, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 4, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 3, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 4, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 3, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 4, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(float, 2, false, CrossOp{}, "Cross");
VECTOR_MATH_BINARY_BENCHMARK_CASE(float, 3, false, CrossOp{}, "Cross");
VECTOR_MATH_TERNARY_BENCHMARK_CASE(float, 4, false, CrossOp{}, "Cross");

VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, SumOp{}, "Sum");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 3, false, SumOp{}, "Sum");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 4, false, SumOp{}, "Sum");
VECTOR_MATH_BINARY_BENCHMARK_CASE(double, 2, false, DotOp{}, "Dot");
VECTOR_MATH_BINARY_BENCHMARK_CASE(double, 3, false, DotOp{}, "Dot");
VECTOR_MATH_BINARY_BENCHMARK_CASE(double, 4, false, DotOp{}, "Dot");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 3, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 4, false, LengthOp{}, "Length");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 3, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 4, false, LengthPreciseOp{}, "LengthPrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 3, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 4, false, NormalizehOp{}, "Normalize");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 3, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 4, false, NormalizePreciseOp{}, "NormalizePrecise");
VECTOR_MATH_UNARY_BENCHMARK_CASE(double, 2, false, CrossOp{}, "Cross");
VECTOR_MATH_BINARY_BENCHMARK_CASE(double, 3, false, CrossOp{}, "Cross");
VECTOR_MATH_TERNARY_BENCHMARK_CASE(double, 4, false, CrossOp{}, "Cross");

} // namespace