Mathter
===
**As the code is fresh, it is subject to interface changes and may contain bugs.**

What is this?
---
Mathter is a linear algebra library focused on 3D game development.
Most other math libraries use a certain set of conventions for vectors and matrices.
Even if you are lucky enough to have these documented, they are most certainly
not *your* conventions. Getting fed up with this was the primary motive to
write another library. You can tell Mathter to use *your* own conventions
by setting a few template parameters, and you get it all generated for you
thanks to C++ template (dark) magic.

Features
---
- (To Be Added (TBA) features in parentheses)
- Header-only: just include Mathter/*
- Vectors, Matrices, (Quaternions TBA), Hyperplanes, Parametric lines
- In any dimension
- Configurable:
  - Base data type: float, double, (int, std::complex, custom TBA)
  - Dimension
  - Preferred multiplication order: vector*matrix or matrix*vector
  - Memory layout: row-major and column-major matrices
  - Tight packing: removes alignment requirements and padding from between elements (disables SIMD as well)
- SIMD acceleration, optimized for small matrices
- Vectors:
  - Arithmetic operators
  - Dot product, cross product (in higher dimensions as well)
  - Concatenation
- Matrices:
  - Arithmetic operators
  - Common: trace, determinant, inverse, transpose
  - Decomposition: LU decomposition, (SVD TBA)
- Geometry:
  - Hyperplane-line and line segment intersection
  - Line-line, line segment-line segment intersection in 2D
  - Distance from plane (line and line segment TBA)

Shortcomings:
- Dimensions specified only at compile time
- Slow compilation: explicit template specialization in your project solves this


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
- Number of rows, or **width** (Columns)
- Number of columns, or **height** (Rows)
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
Notice we had to append a one to the vector to apply the translation. If you leave it off, it is automatically appended during multiplication.

The other way around:
```c++
Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> M_T = M.Transposed();
Vector<float, 3> v_transformed2 = M_T*(v | 1);
```


So how do I use it?
---
Grab the files from Mathter/, add to you include directory, you are done.

A new compiler that supports C++14 is required. It works with MSVC 2015 and 2017.
Clang and GCC support will be added, but not officially supported at the moment (it may work though).

Be aware that the code is fresh, meaning it is subject to bugs and interface changes.


License
---
This code is distrubuted under The Unlicense.
It translates as "do whatever the hell you want with the code, without getting you ass kicked".
You are encouraged (but not required) to indicate the source of this software.


















