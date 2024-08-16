#include "../Benchmark.hpp"
#include "../Fixtures.hpp"

#include <array>
#include <cstddef>
#include <functional>

namespace {
namespace msvc_repro_matrix2x2 {

	template <class T, int Dim>
	struct Vector {
		std::array<T, Dim> elements;

		auto& operator[](size_t idx) {
			return elements[idx];
		}

		const auto& operator[](size_t idx) const {
			return elements[idx];
		}

		static Vector All(T value) {
			Vector v;
			std::fill(v.elements.begin(), v.elements.end(), value);
			return v;
		}
	};


	template <class T, int Dim, class Op, size_t... Indices>
	auto BinaryOp(const Vector<T, Dim>& lhs, const Vector<T, Dim>& rhs, Op&& op, std::index_sequence<Indices...>) {
		return Vector<T, Dim>{ op(lhs[Indices], rhs[Indices])... };
	}


	template <class T, int Dim>
	auto operator+(const Vector<T, Dim>& lhs, const Vector<T, Dim>& rhs) {
		return BinaryOp(lhs, rhs, std::plus<>{}, std::make_index_sequence<Dim>{});
	}


	template <class T, int Dim>
	auto operator*(const Vector<T, Dim>& lhs, const Vector<T, Dim>& rhs) {
		return BinaryOp(lhs, rhs, std::multiplies<>{}, std::make_index_sequence<Dim>{});
	}


	template <class T, int Rows, int Columns>
	struct Matrix {
		std::array<Vector<T, Columns>, Rows> stripes;

		auto& operator()(size_t row, size_t col) {
			return stripes[row][col];
		}

		const auto& operator()(size_t row, size_t col) const {
			return stripes[row][col];
		}
	};


	template <class T, int Rows1, int Match, int Columns2>
	auto operator*(const Matrix<T, Rows1, Match>& lhs, const Matrix<T, Match, Columns2>& rhs) {
		using Vec = Vector<T, Columns2>;
		Matrix<T, Rows1, Columns2> m;
		for (size_t rowIdx = 0; rowIdx < Rows1; ++rowIdx) {
			m.stripes[rowIdx] = Vec::All(lhs(rowIdx, 0)) * rhs.stripes[0];
			for (size_t runIdx = 1; runIdx < Match; ++runIdx) {
				m.stripes[rowIdx] = Vec::All(lhs(rowIdx, runIdx)) * rhs.stripes[runIdx] + m.stripes[rowIdx];
			}
		}
		return m;
	}


	template <class T, int Rows1, int Match, int Columns2>
	auto Unrolled(const Matrix<T, Rows1, Match>& lhs, const Matrix<T, Match, Columns2>& rhs) {
		using Vec = Vector<T, Columns2>;
		return Matrix<T, Rows1, Columns2>{
			Vector<T, Columns2>{ { lhs(0, 0) * rhs(0, 0) + lhs(0, 1) * rhs(1, 0), lhs(0, 0) * rhs(0, 1) + lhs(0, 1) * rhs(1, 1) } },
			Vector<T, Columns2>{ { lhs(1, 0) * rhs(0, 0) + lhs(1, 1) * rhs(1, 0), lhs(1, 0) * rhs(0, 1) + lhs(1, 1) * rhs(1, 1) } }
		};
	}

	using Mat = Matrix<float, 2, 2>;


	template <size_t Size>
	auto MakeInput() {
		std::array<Mat, Size> out;
		for (auto& m : out) {
			for (size_t i = 0; i < 2; ++i) {
				for (size_t j = 0; j < 2; ++j) {
					m(i, j) = float(i == j);
				}
			}
		}
		return out;
	}


	BENCHMARK_CASE("MSVC Repro: float.22 * float.22 (loop)",
				   "[MSVC Repro]",
				   50,
				   64,
				   GenericNAryFixture{ std::multiplies<>{} },
				   MakeInput<1>()[0],
				   MakeInput<4>(),
				   MakeInput<64>());


	BENCHMARK_CASE("MSVC Repro: float.22 * float.22 (unrolled)",
				   "[MSVC Repro]",
				   50,
				   64,
				   GenericNAryFixture{ [](const auto& lhs, const auto& rhs) { return Unrolled(lhs, rhs); } },
				   MakeInput<1>()[0],
				   MakeInput<4>(),
				   MakeInput<64>());

} // namespace msvc_repro_matrix2x2
} // namespace