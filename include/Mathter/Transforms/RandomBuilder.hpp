// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Matrix.hpp"
#include "../Quaternion/Quaternion.hpp"
#include "../Vector/Vector.hpp"
#include "ZeroBuilder.hpp"


namespace mathter {

namespace impl {

	template <class Distribution, class Engine>
	class RandomBuilder {
	public:
		RandomBuilder(Distribution& distribution, Engine& engine) : distribution(distribution), engine(engine) {}

		template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
		operator Matrix<T, Rows, Columns, Order, Layout, Packed>() const {
			Matrix<T, Rows, Columns, Order, Layout, Packed> m;
			for (size_t i = 0; i < Rows; ++i) {
				for (size_t j = 0; j < Columns; ++j) {
					m(i, j) = static_cast<T>(distribution(engine));
				}
			}
			return m;
		}

		template <class T, int Dim, bool Packed>
		operator Vector<T, Dim, Packed>() const {
			Vector<T, Dim, Packed> v;
			for (size_t i = 0; i < Dim; ++i) {
				v[i] = static_cast<T>(distribution(engine));
			}
			return v;
		}

		template <class T, eQuaternionLayout Layout, bool Packed>
		operator Quaternion<T, Layout, Packed>() const {
			return Quaternion<T, Layout, Packed>{
				static_cast<T>(distribution(engine)),
				static_cast<T>(distribution(engine)),
				static_cast<T>(distribution(engine)),
				static_cast<T>(distribution(engine))
			};
		}

		Distribution& distribution;
		Engine& engine;
	};

} // namespace impl


/// <summary> Creates an identity matrix or identity quaternion. </summary>
/// <remarks> If the matrix is not square, it will look like a truncated larger square identity matrix. </remarks>
template <class Distribution, class Engine>
inline auto Random(Distribution& distribution, Engine& engine) {
	return impl::RandomBuilder{ distribution, engine };
}

} // namespace mathter