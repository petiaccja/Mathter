Mathter
===
**It's not really production ready, but keep reading if interested.**

What is this?
---
Mathter is a linear algebra library focused on 2D and 3D game development. It uses the latest C++ features to provide a painless and highy configurable API. As of now, only fixed-size vectors and matrices are supported.

Features
---
- Header-only: just include Mathter/*
- Configureable
  - Float, double, int, or roll your own (only float and double is tested as of now)
  - Dimensions as template parameters
  - Post- and pre-multiplication style for geometrical transforms
  - Row-major and column-major memory layout
  - Tightly pack elements (no SIMD alignment and padding)
- Fast: SIMD accelerated, loops unrolled for small matrices
- Vectors and Matrix arithmetic
- Common matrix operations
- Geometrical transforms (projective coming soon)
- Hyperplanes, parametric lines, intersection
- Matrix decompositions (LU (SVD coming soon)
- Quaternions (coming soon)

Examples
---
Include stuff:
```c++
#include <Mathter/Vector.hpp> // vectors
#include <Mathter/Matrix.hpp> // matrices
#include <Mathter/Geometry.hpp> // lines, planes, intersection
using namespace mathter;
```

Vectors have three parameters:
- Scalar type (T)
- Dimension (Dim)
- Whether to pack tightly (true) or use SIMD (false, default)
```c++
Vector<float, 3, false> u; // SIMD accelerated 3-dimensional vector
Vector<float, 3, false> v(1,2,3);
```

You can take advantage of automatic concatenation:
```c++
Vector<float, 6> u6(7,8,9, u); // 7,8,9,4,5,6
Vector<float, 6> v6 = u | v; // 4,5,6,1,2,3
Vector<double, 6> w6 = { u,v }; // 4,5,6,1,2,3
```
Be careful, different scalar types are implicitly cast.
You might want to disable compiler warnings for this.

Matrices have many parameters:
- Scalar type (T)
- Number of columns, or **width** (Columns)
- Number of rows, or **height** (Rows)
- Geometrical transformation matrices written after (row-vector)(default) or before (column-vector) the vector (Order)
- Row-major (default) or column-major storage (Layout)
- Enable tight packing (for uploading to GPU or serialization)
```c++
Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false> M =
{
	1,	2,	3,
	4,	5,	6,
	7,	8,	9,
	10,	11,	12,
};
```
Elements layed out as on paper. This is not affected by memory layout or order of multiplication.

Let's transform a vector:
```c++
M.SetTranslation(10, 10, 10);
Vector<float, 3> v_transformed1 = (v | 1)*M;
```
Notice we had to append a one to the vector to apply the translation. Since this is not necesserily efficient and nice, I might add implicit extension to vector-matrix multiplication.

The other way around:
```c++
Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> M_T = M.Transposed();
Vector<float, 3> v_transformed2 = M_T*(v | 1);
```


So how do I use it?
---
Grab the files from Mathter/, add to you include directory, you are done.

While the interface is clean, the library uses extensive template magic to achieve stuff with less code, and it's insides are not for the faint-hearted. It applies to compilers as well, so do not expect much from things older than the VS2015 toolset. There will be GCC and Clang support, but not tested yet.

**TL;DR; You don't. It's not ready yet.**



















