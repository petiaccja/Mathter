//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#include "Mathter/Matrix.hpp"
#include "Mathter/Quaternion.hpp"
#include "Mathter/Vector.hpp"


using namespace mathter;

void Example() {
	// Declare some types.
	using Vec2 = Vector<float, 2, false>;
	using Vec3 = Vector<float, 3, false>;
	using Vec4 = Vector<float, 4, false>;
	using Mat33 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
	using Mat43 = Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
	using Mat34 = Matrix<float, 3, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
	using Mat44 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
	using Quat = Quaternion<float, false>;

	// Create some vectors.
	Vec2 a = { 1, 2 };
	Vec3 v1 = { 1, 2, 3 };
	Vec3 v2 = { a, 3 };
	Vec4 v3 = a.xxyy;

	// Use the transform builder functions.
	Mat44 m = Translation(3, 2, 1);
	Mat43 m43 = Translation(3, 2, 1);
	// Mat34 m34 = Translation(3, 2, 1); // Doesn't compile. This should be PRECEDE_VECTOR.

	// Expansion and reduction of vectors.
	Vec3 r11 = v1 * m;
	Vec3 r12 = Transpose(m)*v1;
	Vec3 r2 = Vec3((v2 | 1) * m);

	// Shared builder functions.
	Quat q = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), Deg2Rad(60.f));
	Mat33 mr = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), Deg2Rad(60.f));
	mr = Identity();
	q = Identity();

	// Conversion between rotators.
	mr = Mat33(q);
	q = Quat(mr);

	// Rotate vectors via quaternions.
	Vec3 vr = v1 * q;

	// Solve a linear equation.
	Vector<float, 6> b = { 1, 2, 3, 4, 5, 6 };
	Matrix<float, 6, 6> M;
	for (int row = 0; row < M.RowCount(); ++row) {
		for (int col = 0; col < M.RowCount(); ++col) {
			M(row, col) = rand() % 1000 / 1000.f;
		}
	}
	Vector<float, 6> x = DecomposeLUP(M).Solve(b); // Mx = b
}