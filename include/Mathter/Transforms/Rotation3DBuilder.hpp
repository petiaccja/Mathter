// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Matrix/Math.hpp"
#include "../Matrix/Matrix.hpp"
#include "../Quaternion/Quaternion.hpp"
#include "../Vector.hpp"
#include "IdentityBuilder.hpp"

#include <array>
#include <cmath>
#include <stdexcept>


namespace mathter {

namespace impl {

	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class A>
	Matrix<T, Rows, Columns, Order, Layout, Packed> RotationMatrixSingleAxis(A angle, int axis) {
		assert(0 <= axis && axis < 3);

		Matrix<T, Rows, Columns, Order, Layout, Packed> m = Identity();

		// Indices according to follow vector order
		auto elem = [&m](int i, int j) -> auto& {
			return Order == eMatrixOrder::FOLLOW_VECTOR ? m(i, j) : m(j, i);
		};

		const auto c = static_cast<T>(std::cos(angle));
		const auto s = static_cast<T>(std::sin(angle));
		if (axis == 0) {
			// Rotate around X
			elem(1, 1) = c;
			elem(1, 2) = s;
			elem(2, 1) = -s;
			elem(2, 2) = c;
		}
		else if (axis == 1) {
			// Rotate around Y
			elem(0, 0) = c;
			elem(0, 2) = -s;
			elem(2, 0) = s;
			elem(2, 2) = c;
		}
		else {
			// Rotate around Z
			elem(0, 0) = c;
			elem(0, 1) = s;
			elem(1, 0) = -s;
			elem(1, 1) = c;
		}

		return m;
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	Matrix<T, Rows, Columns, Order, Layout, Packed> RotationMatrixAxisAngle(const Vector<T, 3, Packed>& axis, T angle) {
		using M33 = Matrix<T, 3, 3, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;
		using M31 = Matrix<T, 3, 1, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;
		using M13 = Matrix<T, 1, 3, eMatrixOrder::PRECEDE_VECTOR, Layout, Packed>;

		const auto c = std::cos(angle);
		const auto s = std::sin(angle);

		const M31 u(axis);
		const M13 uT(axis);
		const M33 cross = {
			T(0), -axis.z, axis.y,
			axis.z, T(0), -axis.x,
			-axis.y, axis.x, T(0)
		};

		const auto r = c * M33(Identity()) + s * cross + (T(1) - c) * u * uT;
		Matrix<T, Rows, Columns, Order, Layout, Packed> m = Identity();
		m.Insert(0, 0, Matrix<T, 3, 3, Order, Layout, Packed>(r));
		return m;
	}


	template <class T, eQuaternionLayout Layout, bool Packed>
	Quaternion<T, Layout, Packed> RotationQuaternionAxisAngle(const Vector<T, 3, Packed>& axis, const T& angle) {
		const auto c = std::cos(angle / T(2));
		const auto s = std::sin(angle / T(2));

		return Quaternion<T, Layout, Packed>(c, s * axis);
	}

	template <class T>
	class Rotation3DAxisBuilder {
	public:
		Rotation3DAxisBuilder(const T& angle, int axis) : angle(angle), axis(axis) {}
		Rotation3DAxisBuilder& operator=(const Rotation3DAxisBuilder&) = delete;

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 4, Order, Layout, MPacked>() const {
			return RotationMatrixSingleAxis<U, 4, 4, Order, Layout, MPacked>(U(angle), axis);
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 3, Order, Layout, MPacked>() const {
			return RotationMatrixSingleAxis<U, 3, 3, Order, Layout, MPacked>(U(angle), axis);
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
			return RotationMatrixSingleAxis<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>(U(angle), axis);
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
			return RotationMatrixSingleAxis<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>(U(angle), axis);
		}

		template <class U, eQuaternionLayout Layout, bool QPacked>
		operator Quaternion<U, Layout, QPacked>() const {
			Vector<U, 3, QPacked> vec = Zero();
			vec[axis] = U(1);
			return RotationQuaternionAxisAngle<U, Layout, QPacked>(vec, U(angle));
		}

	private:
		T angle;
		int axis;
	};


	template <class T>
	class Rotation3DTriAxisBuilder {
	public:
		Rotation3DTriAxisBuilder(const std::array<T, 3>& angles, std::array<int, 3> axes) : angles(angles), axes(axes) {}
		Rotation3DTriAxisBuilder& operator=(const Rotation3DTriAxisBuilder&) = delete;

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 4, Order, Layout, MPacked>() const {
			return Make<U, 4, 4, Order, Layout, MPacked>();
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 3, Order, Layout, MPacked>() const {
			return Make<U, 3, 3, Order, Layout, MPacked>();
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
			return Make<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>();
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
			return Make<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>();
		}

		template <class U, eQuaternionLayout Layout, bool QPacked>
		operator Quaternion<U, Layout, QPacked>() const {
			Vector<U, 3, QPacked> vec0 = Zero();
			vec0[axes[0]] = U(1);
			Vector<U, 3, QPacked> vec1 = Zero();
			vec1[axes[1]] = U(1);
			Vector<U, 3, QPacked> vec2 = Zero();
			vec2[axes[2]] = U(1);
			const auto r0 = RotationQuaternionAxisAngle<U, Layout, QPacked>(vec0, U(angles[0]));
			const auto r1 = RotationQuaternionAxisAngle<U, Layout, QPacked>(vec1, U(angles[1]));
			const auto r2 = RotationQuaternionAxisAngle<U, Layout, QPacked>(vec2, U(angles[2]));
			return r2 * r1 * r0;
		}

	private:
		template <class U, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		auto Make() const {
			const auto r0 = RotationMatrixSingleAxis<U, 3, 3, Order, Layout, MPacked>(U(angles[0]), axes[0]);
			const auto r1 = RotationMatrixSingleAxis<U, 3, 3, Order, Layout, MPacked>(U(angles[1]), axes[1]);
			const auto r2 = RotationMatrixSingleAxis<U, 3, 3, Order, Layout, MPacked>(U(angles[2]), axes[2]);
			const auto r = Order == eMatrixOrder::FOLLOW_VECTOR ? r0 * r1 * r2 : r2 * r1 * r0;
			Matrix<U, Rows, Columns, Order, Layout, MPacked> m = Identity();
			m.Insert(0, 0, r);
			return m;
		}

		std::array<T, 3> angles;
		std::array<int, 3> axes;
	};


	template <class T, bool Packed>
	class Rotation3DAxisAngleBuilder {
	public:
		Rotation3DAxisAngleBuilder(const Vector<T, 3, Packed>& axis, T angle) : axis(axis), angle(angle) {}
		Rotation3DAxisAngleBuilder& operator=(const Rotation3DAxisAngleBuilder&) = delete;

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 4, Order, Layout, MPacked>() const {
			return RotationMatrixAxisAngle<U, 4, 4, Order, Layout, MPacked>(axis, U(angle));
		}

		template <class U, eMatrixOrder Order, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 3, Order, Layout, MPacked>() const {
			return RotationMatrixAxisAngle<U, 3, 3, Order, Layout, MPacked>(axis, U(angle));
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>() const {
			return RotationMatrixAxisAngle<U, 3, 4, eMatrixOrder::PRECEDE_VECTOR, Layout, MPacked>(axis, U(angle));
		}

		template <class U, eMatrixLayout Layout, bool MPacked>
		operator Matrix<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>() const {
			return RotationMatrixAxisAngle<U, 4, 3, eMatrixOrder::FOLLOW_VECTOR, Layout, MPacked>(axis, U(angle));
		}

		template <class U, eQuaternionLayout Layout, bool QPacked>
		operator Quaternion<U, Layout, QPacked>() const {
			return RotationQuaternionAxisAngle<U, Layout, QPacked>(axis, U(angle));
		}

	private:
		const Vector<T, 3, Packed> axis;
		const T angle;
	};

} // namespace impl


/// <summary> Rotates around coordinate axis. </summary>
/// <param name="axis"> 0 for X, 1 for Y, 2 for Z and so on... </param>
/// <param name="angle"> Angle of rotation in radians. </param>
/// <remarks> Positive angles rotate according to the right-hand rule in right-handed
///		coordinate systems (left-handed according to left-hand rule).
template <class T>
auto RotationAxis(T angle, int axis) {
	return impl::Rotation3DAxisBuilder(angle, axis);
}


/// <summary> Rotates around coordinate axis. </summary>
/// <typeparam name="Axis"> 0 for X, 1 for Y, 2 for Z and so on... </typeparam>
/// <param name="angle"> Angle of rotation in radians. </param>
/// <remarks> Positive angles rotate according to the right-hand rule in right-handed
///		coordinate systems (left-handed according to left-hand rule).
template <int Axis, class T>
auto RotationAxis(T angle) {
	return RotationAxis(angle, Axis);
}


/// <summary> Rotates around the X axis according to the right (left) hand rule in right (left) handed systems. </summary>
/// <param name="angle"> Angle of rotation in radians. </param>
template <class T>
auto RotationX(T angle) {
	return RotationAxis<0>(angle);
}


/// <summary> Rotates around the Y axis according to the right (left) hand rule in right (left) handed systems. </summary>
/// <param name="angle"> Angle of rotation in radians. </param>
template <class T>
auto RotationY(T angle) {
	return RotationAxis<1>(angle);
}


/// <summary> Rotates around the Z axis according to the right (left) hand rule in right (left) handed systems. </summary>
/// <param name="angle"> Angle of rotation in radians. </param>
template <class T>
auto RotationZ(T angle) {
	return RotationAxis<2>(angle);
}


/// <summary> Rotates around three axes in succession. </summary>
/// <remarks> Axes: 0 for X, 1 for Y and 2 for Z.
///		Angles in radians. Each rotation according to the right (and left) hand rule in right (and left) handed systems. </remarks>
template <int FirstAxis, int SecondAxis, int ThirdAxis, class T>
auto RotationAxis3(T angle0, T angle1, T angle2) {
	return impl::Rotation3DTriAxisBuilder(std::array<T, 3>{ angle0, angle1, angle2 }, std::array<int, 3>{ FirstAxis, SecondAxis, ThirdAxis });
}


/// <summary> Rotation matrix from Euler angles. Rotations are Z-X-Z. </summary>
/// <param name="z1"> Angle of the first rotation around Z in radians. </param>
/// <param name="x2"> Angle of the second rotation around X in radians. </param>
/// <param name="z3"> Angle of the third rotation around Z in radians. </param>
/// <remarks> Each rotation according to the right (and left) hand rule in right (and left) handed systems. </remarks>
template <class T>
auto RotationEuler(T z1, T x2, T z3) {
	return RotationAxis3<2, 0, 2>(z1, x2, z3);
}


/// <summary> Rotation matrix from roll-pitch-yaw angles. Rotations are X-Y-Z. </summary>
/// <param name="z1"> Angle of the first rotation around X in radians. </param>
/// <param name="x2"> Angle of the second rotation around Y in radians. </param>
/// <param name="z3"> Angle of the third rotation around Z in radians. </param>
/// /// <remarks> Each rotation according to the right (and left) hand rule in right (and left) handed systems. </remarks>
template <class T>
auto RotationRPY(T x1, T y2, T z3) {
	return RotationAxis3<0, 1, 2>(x1, y2, z3);
}


/// <summary> Rotates around an arbitrary axis. </summary>
/// <param name="axis"> Axis of rotation, must be normalized. </param>
/// <param name="angle"> Angle of rotation in radians. </param>
/// <remarks> Right-hand (left-hand) rule is followed in right-handed (left-handed) systems. </remarks>
template <class T, bool Vpacked, class U>
auto RotationAxisAngle(const Vector<T, 3, Vpacked>& axis, U angle) {
	assert(std::abs(T(1) - Length(axis)) < T(0.05) && "The axis must be normalized.");
	return impl::Rotation3DAxisAngleBuilder(axis, T(angle));
}

} // namespace mathter