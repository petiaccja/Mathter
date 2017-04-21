#include "Mathter/Vector.hpp" // vectors
#include "Mathter/Matrix.hpp" // matrices
#include "Mathter/Geometry.hpp" // lines, planes, intersection

using namespace mathter;

void foo() {

	Vector<float, 3, false> u; // SIMD accelerated 3-dimensional vector
	Vector<float, 3, false> v(1,2,3);

	u[0] = 4;
	u(1) = 5;
	u.z = 6;

	Vector<float, 6> u6(7,8,9, u); // 7,8,9,4,5,6
	Vector<float, 6> v6 = u | v; // 4,5,6,1,2,3
	Vector<double, 6> w6 = { u,v }; // 4,5,6,1,2,3

	
	Matrix<float, 3, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false> M =
	{
		1,	2,	3,
		4,	5,	6,
		7,	8,	9,
		10,	11,	12,
	};


	M.SetTranslation(10, 10, 10);
	Vector<float, 3> v_transformed1 = (v | 1)*M;
	
	Matrix<float, 4, 3, eMatrixOrder::PRECEDE_VECTOR> M_T = M.Transposed();
	Vector<float, 3> v_transformed2 = M_T*(v | 1);




}