Mathter
===

Clang5, GCC7: 
[![Build Status](https://travis-ci.org/petiaccja/Mathter.svg?branch=master)](https://travis-ci.org/petiaccja/Mathter)

MSVC 2017:
[![Build status](https://ci.appveyor.com/api/projects/status/6uvfnfgp5paha8kw?svg=true)](https://ci.appveyor.com/project/petiaccja/mathter)

Introduction
---
Mathter is a **linear algebra** library with focus on **game development**, however it may be useful for other applications where **small-matrix** linear algebra or **3D** coordinate calculations are needed.

Yet another 3D math library...? I wrote Mathter because non of the libraries I tried used the conventions I wanted, and each of them had something about the syntax that bugged me. So will Mathter support your convetions? Yes, it supports everything through templates. Will Mathter's syntax bug you? Yes, it will, but it was designed to do it as little as possible. Read on for details.

Features
---
- General:
  - SIMD (for important single-precision math only, for now)
  - C++17
  - Compile-time vector size (NO runtime resize on vectors/matrices)
- Template parametrization & convention configurations:
  - T: any scalar type can be specified for Vectors, Matrices, etc. (int, float, double, complex). Of course, some of these make no sense for all operations.
  - Dimensions: vector length, matrix row and column count.
  - Order: whether transform matrices (i.e. scale, rotate) should act as v\*M or M\*v.
  - Layout: row-major or column major layout for matrices.
  - Packed: disable SIMD, forced alignment, and pack vector/matrix elements tightly in memory. Use Packed=true to handle GPU uploads.
- Mathematical primitives:
  - Vectors [v[i], v(i), v.x, swizzling: v.zyx, v.xy]
  - Matrices [m(row,col)]
  - Quaternions [q.w, q.x, q.y, q.z]
- Arithmetic:
  - Vector\*Vector
  - Vector\*Matrix
  - Vector\*Quat
  - Matrix\*Matrix
  - Quat\*Quat
- Transformations (they work in higher dimensions, too):
  - Rotation (only 2D, 3D), scale, translation, orthographic projection, perspective projection, camera look-at, shear
- Common functions:
  - Vectors: length, dot product, cross product, normalization, ...
  - Matrices: trace, determinant, inverse, norm, transpose, ...
  - Quaternions: length/abs, normalization, conjugate/inverse, exp, log, pow, ...
- Matrix decompositions:
  - LU & LUP
  - QR
  - SVD
- Geometry:
  - Lines, line segments, rays
  - Planes (hyperplanes in N-dimensions)
  - Triangles in 3D
  - Bezier curves
  - Intersections: line-hyperplane, segment-hyperplane, ray-triangle (3D), line/segment-line/segment (2D)
- Utility:
  - Math constants
  - Radians <=> degrees
  - Standard stream I/O of vectors and matrices (this is kinda hacked together TBH)
- Usage:
  - If it makes sense, it compiles (if it doesn't make sense, it doesn't compile)
  - Implicit conversions at many places, for convenience
  - Concise and powerful syntax via some C++ magic
  - Prefers mathematical notation (i.e. quaternions are chained qr = q3*q2*q1)
  - Many things are generalized to higher dimensions (i.e. hyperplanes, N-dimensional transforms & functions), but interfaces are not compromised

For details, browse [Mathter's wiki](https://github.com/petiaccja/Mathter/wiki).

Code example
---

Before you do anything, you want to declare types to you own taste:
```c++
using namespace mathter;
using Vec2 = Vector<float, 2, false>;
using Vec3 = Vector<float, 3, false>;
using Vec4 = Vector<float, 4, false>;
using Mat33 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat43 = Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat34 = Matrix<float, 3, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat44 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Quat = Quaternion<float, false>;
```
Remember, all stuff is in the namespace mathter.

You have multiple ways to initialize your variables. Beware, **default constructors will garbage-initialize** (NOT null-initialize). This saves 2 clock cycles when you initialize later.
```c++
Vec2 a = { 1, 2 };
Vec3 v1 = { 1, 2, 3 };
Vec3 v2 = { a, 3 };
Vec4 v3 = a.xxyy;
Mat33 mdecl = {
  1, 0, 0,
  0, 1, 0,
  0, 0, 1,
};
```

Transforms are for the most part represented by matrices. You can use the builder functions to construct them. It is important that the builder functions don't create any matrices, you have to assign them to a matrix with a well-defined type. As long as the destination matrix can represent the transform, you'll be fine:
```c++
Mat44 m = Translation(3, 2, 1); // Always fine, but the 3 values are placed differently with FOLLOW vs. PRECEDE.
Mat43 m43 = Translation(3, 2, 1); // This is FOLLOW_VECTOR, fine.
// Mat34 m34 = Translation(3, 2, 1); // Doesn't compile. This should be PRECEDE_VECTOR.
```

You can multiply vectors and matrices fairly liberally. The first two operations extend the 3-vector to a 4-vector by appending a 1 at the end, then multiply by the 4x4 matrix, do the perspective division with the last element of the resulting 4-vector, then truncate the result to a 3-vector. The last line does this manually, without perspective division.
```c++
Vec3 r11 = v1 * m;
Vec3 r12 = Transpose(m)*v1;
Vec3 r2 = Vec3((v2 | 1) * m);
```

The rotation builders work much the same for matrices and quaternions:
```c++
Quat q = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), Deg2Rad(60.f));
Mat33 mr = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), Deg2Rad(60.f));
mr = Identity();
q = Identity();
```

And you can also convert between the 2 representations:
```c++
mr = Mat33(q);
q = Quat(mr);
```

Applying rotations to vectors with quaternions is straightforward, but in general it's faster with matrices:
```c++
Vec3 vr = v1 * q;
```

You can also do other things like solving linear equations via an LUP decomposition of the matrix. Apart from rotations, things work in higher dimensions, too:
```c++
Vector<float, 6> b = { 1, 2, 3, 4, 5, 6 };
Matrix<float, 6, 6> M;
for (int row = 0; row < M.RowCount(); ++row) {
	for (int col = 0; col < M.RowCount(); ++col) {
		M(row, col) = rand() % 1000 / 1000.f;
	}
}
Vector<float, 6> x = DecomposeLUP(M).Solve(b); // Mx = b
```

Installation
---
Mathter is **header-only**, you only need to add the Mathter folder to your include path. You will **need a C++17 compliant compiler** to use the library. Builds on Clang5, GCC7 and MSVC 2017.


License
---
This code is distrubuted under **The Unlicense**. (I may change it MIT or something, but it will always be just as permissive.)
It translates as "**do whatever the hell you want with the code**, without getting your ass kicked". You are not required to credit me (but I'm happy if you do 🙂). In any case, don't be a cunt and claim it's your work.