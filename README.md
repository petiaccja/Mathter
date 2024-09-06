Mathter
===

![Language](https://img.shields.io/badge/Language-C++17-blue)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)
[![Build & test](https://github.com/petiaccja/Mathter/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/petiaccja/Mathter/actions/workflows/build_and_test.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_Mathter&metric=alert_status)](https://sonarcloud.io/dashboard?id=petiaccja_Mathter)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_Mathter&metric=coverage)](https://sonarcloud.io/dashboard?id=petiaccja_Mathter)


Introduction
---
Mathter is a header-only linear algebra library for game development and scientific applications.

Why yet another 3D math library?
- Existing libraries often have fixed conventions & notation, but Mathter is fully configurable:
  - Scalar types: floating point, integer, or complex
  - Dimensions: arbitrary vector and matrix sizes, algorithms generalized to 2D, 3D, and up, where possible
  - Multiplication order: matrices before vectors, vectors before matrices
  - Memory layout: row-major & column-major matrices, sijk & ijks quaternions
  - Packing: tightly packed objects help with memory transfer to GPU, aligned and padded objects help with SIMD
  - Customized projections: view transforms and projections are fully configurable to accomodate any normalized device coordinates
- Intuitive & safe API:
  - Conventions & notation configured via template arguments
  - Full vector swizzling
  - Several small features, like arbitrary vector concatenation, using translation matrices without augmented vectors, or automatic perspective division
  - Multiplication order is enforced at compile time
- Extendable:
  - Geometric transforms are free functions, independent of the underlying linalg objects
  - You can add your own within your project, no need to modify Mathter


Example code:
```c++
#include <Mathter/Decompositions/DecomposeSVD.hpp>
#include <Mathter/Matrix.hpp>
#include <Mathter/Transforms.hpp>
#include <Mathter/Vector.hpp>

using namespace mathter;

using Vec3 = Vector<float, 3, false>
using Mat44 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;

const Mat44 preRotation = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), 0.4f);
const Mat44 scale = Scale(4, 5, 6);
const Mat44 postRotation = RotationAxisAngle(Normalize(Vec3(1, 2, 3)), -0.4f);
const Mat44 translation = Translation(3, 2, 1);

const auto transform = preRotation * scale * postRotation * translation;
const auto local = transform.template Extract<3, 3>(0, 0);
const auto [u, s, v] = DecomposeSVD(local);

const Vec3 original = { 5, 4, 6 };
const Vec3 transformed = original * transform;
const Vec3 transformed = transform * original; // Compilation error due to matrix order.
```


For more detailed information about using Mathter, **read the** [**guide**](https://github.com/petiaccja/Mathter/blob/master/docs/Guide.md).

Features
---
- General:
  - Requires C++17 and above
  - SIMD acceleration (using the [XSimd](https://github.com/xtensor-stack/xsimd) library)
    - Optional, and Mathter works without dependencies as well
  - Header-only
  - Arbitrary compile-time dimensions
    - Mathter does not support dynamic sizes, like [Eigen](https://eigen.tuxfamily.org)
    - Larger matrices, like 5x5, are supported, unlike in other 3D libraries
    - Algorithms (e.g. shear transform, decompositions) generalize to higher dimensions
  - Prefers mathematical notation
    - E.g. quaternions are chained like q3\*q2\*q1
- Linear algebraic objects:
  - Vectors
    - Swizzlers
  - Matrices
  - Quaternions
- Geometric primitives:
  - Bezier curves
  - Hyperplanes
  - Lines
  - Line segments
  - Rays
  - Triangles
- Arithmetc & algorithms:
  - Arithmetic operations on vectors, matrices, and quaternions
  - Mathematical functions:
    - Dot product
    - Cross product
    - Inverse
    - Conjugate transpose
    - \+ many more
  - Intersections of geometric primitives
- Coordinate transformations:
  - Identity
  - Orthographic projection
  - Perspective projection
  - Rotations (in 2D and 3D)
  - Scaling
  - Shear
  - Translation
  - Camera look-at
- Matrix decompositions & systems of linear equations
  - LU & LUP, QR & LQ, SVD
  - Each have solvers for:
    - Systems of equations
    - Matrix inverse
    - Psudoinverse (where applicable)
    - Linear least squares (where applicable)

Installation
---

**Getting Mathter:**
1. Using conan: https://conan.io/center/mathter
2. Using vcpkg: https://vcpkg.io/en/package/mathter
3. Using CMake: Mathter uses CMake as its build system; configure and install Mathter and use the installed `MathterConfig.cmake`.
4. Manually: Mathter is header-only, so you can just copy the files and get going

**Compiler support:**
- Requires **C++17** or above (tested with '17 and '20)
- Supports major compilers:
  - **GCC** (tested on v13)
  - **Clang** (tested on v17)
  - **MSVC** (tested on v19.3x)
    - Use the included `Mathter.natvis` file to get pretty printed types in the VS debugger

**Build flags:**
- Set `MATHTER_BUILD_TESTS:BOOL=OFF` and `MATHTER_BUILD_BENCHMARKS:BOOL=OFF` if you don't need them
- Set `MATHTER_ENABLE_SIMD:BOOL` according to whether you have [XSimd](https://github.com/xtensor-stack/xsimd) installed


Building the tests
---

To build Mathter and run tests, you'll need [conan (version 2)](https://docs.conan.io/2/installation.html):

```
pip install conan
```

You then have to create a `conan` profile:

```
conan profile detect
```

It's recommended to verify and edit your conan profile as needed.

You then have to install the dependencies using `conan`, configure with `CMake`, and build:

```
conan install . --build=missing -pr:h=default -pr:b=default
cmake . --preset conan-debug
cmake --build build/debug
```

Instead of the `default` profile, you can also use the profiles in `.github/build_profiles`.

Once the build is finished, you can run `./bin/UnitTest`, which is a [Catch2](https://github.com/catchorg/Catch2) application.

If in doubt, you can always look at the CI workflows for an example of building and running the tests.


License
---
The code is using the **MIT license**, which is a very permissive license suitable for non-commercial and commercial uses alike. However, you have to include the copyright notice in your code. Read the full license for the exact terms.