Guide for Mathter
===

1. [Contents of Mathter](#structure)
2. [Installation](#installation)
3. [Configuration](#configuration)
4. [Vectors](#vectors)
5. [Matrices](#matrices)
6. [Quaternions](#quaternions)
7. [Transforms](#transforms)


Contents of Mathter
<a name="structure"></a>
---

The library is split into multiple pieces that address independent domains of 3D math calculations.

#### Math primitives

Mathter provides an implementation for the following mathematical primitives:
- Vector
- Matrix
- Quaternion

<sup>Vectors and matrices are a distinct type to make the coding of 3D math easier. You can use matrices for vectors, but special vector functions are not provided for those.</sup>

The primitive types support multiple ways to initialize them and access their elements. Their size and scalar type is specificied as template parameters, among other options.


#### Arithmetic operators

The library provides overloaded arithmetic operators to multiply, divide, add and subtract the primitives. You can do all the sensible operations, such as multiply matrices with matrices, matrices with vectors, vectors with vectors, quaternions with quaternions or quaternions with vectors, etc., using operators \*, /, +, and - in C++.


#### Transform builders

In Mathter, you can build transform matrices and quaternions via free functions such as ```mathter::Perspective(...)```. The functions are evaluated lazily, assign the result of ```mathter::Perspective(...)``` to a matrix that can hold the transform to produce a transform matrix.


#### Vector swizzling

Vectors can be accessed like ```v.xxyy``` just like in GLSL or HLSL. You can also do arithmetic on the swizzlers.


#### Common vector functions

The library provides common functions on vectors, such as length, normalization or dot product. These are free functions in the mathter namespace.


#### Common vector functions

Similarly to the vectors, you can calculate the determinant, trace or inverse of a matrix. These are also provided as free functions in the mathter namespace.


#### Matrix decompositions

Several common matrix decompositions are included. Use the ```Decompose*``` functions to analyze a matrix.


#### Utility

The most important piece is probably numerical constants and conversion between degrees and radians.


#### Geometry

Additional primitives for lines, line segments, rays, triangles and hyperplanes are provided. You can do intersection testing on most reasonable combination of these.



Installation
<a name="installation"></a>
---

Mathter is a header-only library, you simply have to add the *Mathter* to your include path. You can compile the tests with CMake if you want to test your compiler or run the tests yourself. You will need a C++17 compliant compiler for both the code and the tests.

When using Mathter, ignore the files in the subfolders, and include these
```c++
#include <Mathter/Vector.hpp>
#include <Mathter/Matrix.hpp>
#include <Mathter/Quaternion.hpp>
#include <Mathter/Utility.hpp>
#include <Mathter/Geometry.hpp>
#include <Mathter/IoStream.hpp>
```
as needed.

There is also an extra ```.natvis``` file. If you add this file in a Microsoft Visual Studio project, it will display the primitives in pretty text in the debugger. (Column-major matrices appear transposed!)



Configuration
<a name="configuration"></a>
---

Mathter was designed to support any notation or convention that you or your graphics API wants (even when these two are different). You are no longer locked to the notation you hate, neither do you have to suffer with DirectX conventions when using OpenGL.

Here you find a list of important parameters that you first want to set to make Mathter behave as you want it to. It is recommended that you typedef Mathter types and wrap Mathter functions to give yourself a simpler interface without the parameters.


#### Template parameters

##### Post- and pre-multiplication of vectors by matrices

Matrices accept a template parameter ```Order```, which can be either of:
- ```eMatrixOreder::FOLLOW_VECTOR```: you will have multiply matrices and vector like ```v*M```
- ```eMatrixOreder::PRECEDE_VECTOR```: you will have multiply matrices and vector like ```M*v```

To avoid the post/pre naming confusion, these tell whether your matrices *follow* or *precede* the vectors you multiply with them. All transformation matrices you construct will adhere to the option you specify here. Your rules of vector-matrix multiplication are enforced at compile-time. Keep in mind that this also affects the order in which you multiply matrices when chaining transforms: the matrix first acting on the vector needs to be right next to the vector, for example ```M3*M2*M1*v``` or ```v*M1*M2*M3```.

##### Row-major vs column-major matrix layout

Matrices takes another template parameter called ```Layout```, which can be either of:
- ```eMatrixLayout::ROW_MAJOR```
- ```eMatrixLayout::COLUMN_MAJOR```

These don't affect the syntax of the code you write, only change the way the two dimensional matrix is laid out in the sequential memory. Note however, that **performance can be affected** with vector-matrix multiplications. It is faster when you have row-vectors with row-major layout, and column-vectors with column-major layout, because it can better use SIMD. Nevertheless, the correctness of the calculations is never affected.

##### Packing

All types accept an additional parameter called ```Packed```. Its default value is false, meaning the type will have an arbitrary size and arbitrary memory alignment which is optimal for the speed of calculations. When set to true, however, Mathter will impose no extra alignment and will pack objects tightly. Take a 3x3 float matrix as an example. Without packing, the matrix will consume the space of 3x4 float values, because its rows (or columns) are packed to 4 floats to leverage SIMD. When packed, however, the matrix will consume the space of exactly 3x3 floats, as one would expect. Packing is useful when uploading to the GPU, but on the CPU, you should always use unpacked types for performance. Like layout, packing does not change arithmetic or behaviour.


#### Runtime parameters

##### Camera look-at matrices

Mathter has no concept of handedness. When making look-at transforms, you instead have the option to specify the look direction and the up-vector of the camera to any vector. Additionally, you have the option to invert any of the coordinate axes during the transform. With these parameters, you can produce any coordinate system, that you would like.

##### Perspective projection matrices

Perpsective matrices work in the same coordinate system as the look-at matrix. (You may still flip axes, but it is not adviced.) You can specify arbitrary near and far planes (as long as they have the same sign), and you can specify arbitrary NDC Z coordinates. This way, you can do projections regardless of what NDC coordinates your graphics API uses or if the camera looks towards the positive or negative Z axis.

#### Global macros

##### Uninitalized object behaviour

In the following code sample, ```v``` is not initialized:
```c++
Vector<float, 3> v;
```
By default, in **debug builds** Mathter initializes all members of ```v``` to a NaN, or the maximum value the type can represent when it has no NaN. In **release builds** (when ```NDEBUG``` defined), Mathter does not initialize the member of ```v``` at all, meaning they will contain memory garbage. This behaviour is to help with debugging while maintaining maximum performance in release mode. You may want to enable floating point exceptions as well to catch signaling NaNs.

If you, however, want to **change this behviour**, you have the option to define exactly one of:

- ```MATHTER_NULL_INITIALIZE```: initializes vector, matrix and quaternion types with all zeros
- ```MATHTER_INVALID_INITIALIZE```: initializes with all signaling NaNs (default in debug builds)
- ```MATHTER_DONT_INITIALIZE```: does not initialize at all (default in release builds)

When any of these are defined, the default behaviour is overridden regardless of the build type. (The defined value does not matter, just defined them to ```1```.)

It is **important** to note that this feature is not meant for production, but only to aid debugging. Don't ever rely on your variables being initialized to zero, because it is very likely that the flag will be forgotten when someone else builds your application. Always fully initialize Mathter variables, and use these flags to help you catch uninitialized variables.

The reason **why** Mathter does not initialize objects to a default value is because these types have no logical restrictions on what data they should hold, neither do they have a single sensible default value. It is not possible to decide whether matrices should be initialized to identity or all-zero -- both are equally good and equally bad.

Vectors
<a name="vectors"></a>
---

#### The Vector class 

Vectors have three template parameters:
```c++
template <class T, int Dim, bool Packed>
class Vector;
```

With ```T```, you can specify which underlying type the vector uses for its elements. You can set it to ```float```, ```double```, ```int```, ```std::complex``` or any other type or your preference. Floating point types are well-tested, integers also work, there is no reason others wouldn't work but they are not tested. Note, however, that some operations will fail to compile or will provide incorrect results when they don't make sense. For example, you cannot normalize an integer vector because the result is not representable, you cannot ```Min/Max``` a complex vector because complex numbers have no ordering, and the length of an integer vector will get rounded to an integer, which is not what you generally want.

The dimension specifies the number of elements in the vector. It can be any positive number, including 1. Dynamically sized vectors (and matrices) are not supported yet.

The ```Packed``` flag specifies whether the vector should ditch SIMD and forced alignment to tightly pack its elements. Read more in [configuration](#configuration).


#### Creating vectors

You can use one of the several constructors:
```c++
using Vec3 = Vector<float, 3, false>;
using Vec2 = Vector<float, 2, false>;

float data[3] = {...};

Vec2 v1 = {1, 2}; // the elements one by one
Vec3 v2(42.f); // all elements to the same value
Vec3 v3 = { v1, 1 }; // by concatenating multiple entities
Vec3 v4 = { data }; // from a pointer

// these are useful when working with homogeneous coordinates
Vec3 v5 = Vec3(v1); // by appending a one to a shorter vector
Vec2 v6 = Vec2(v4); // by trimming the last element of a longer vector
```

Additionally, you can use explicit conversions between vectors of the same size (but different type or packing).

#### Element access

You can access the individual elements of the vectors by indices or coordinate axis names:
```c++
Vec3 v;
v.x;
v(0);
v[0]; // same as the above two
```

#### Arithmetic on vectors

The elementwise operators are naturally provided by overloads in C++:
```c++
Vec3 v1 = { 1, 2, 3 };
Vec3 v2 = { 1, 2, 3 };
Vec3 v3 = v1 * v2;
```

You can also use scalars to act on vectors:
```c++
Vec3 v1 = { 1, 2, 3 };
Vec3 r1 = 2.0f * v1; // { 2, 4, 6 }
Vec3 r2 = 2.0f / v1; // { 2, 1, 0.666 }
Vec3 r3 = 1.0f + v1; // { 2, 3, 4 }
```

#### Swizzle operators

Swizzle operators can be used to select or modify subsets of a vector:
```c++
Vec3 v;
v.xy = v.yx;
```

You can also perform arithmetic on swizzle operators:
```c++
Vec3 v;
Vec3 r = 0.26f * v.xxx + 0.68f * v.yyy + 0.06f * v.zzz;
```


#### Vector functions

Common vector functions are provided by the library as free functions:
```c++
Vec3 v1 = ...;
Vec3 v2 = ...;
float l = Length(v1);
float angle = std::acos(std::clamp(Dot(v1, v2), -1.0f, 1.0f));
```

For the complete list, use your code completion or browse the code.


#### Limitations

If you have a keen eye, you noticed that the example codes above initialize float vectors with integers. This is not entirely safe, but the library forces the conversion so you don't get a compiler warning.

A similar issue arises when multiplying vectors with scalars of different types. This is allowed, but again, an explicit type conversion is forced by the library, giving no warnings.

You can only perform operations on vectors of the same type, and you must use explicit casts to convert vectors.




Matrices
<a name="matrices"></a>
---

#### The matrix class

Matrices have 6 template parameters:
```c++
class Matrix<class T,
             int Rows,
             int Columns,
             eMatrixOrder Order,
             eMatrixLayout Layout,
             bool Packed>;
```

```T``` specifies what type the elements of the matrix are. Generally, the same applies as for [vectors](#vectors).

```Rows``` and ```Columns``` specify the number of rows and columns of the matrix. Like vectors, dynamically sized matrices are not supported yet.

The ```Order```, ```Layout``` and ```Packed``` are explained in detail in the [configuration section](#configuration).


#### Creating matrices

Matrices can be created by specifying each element of the matrix:

```c++
using Mat44 = Matrix<float, 4, 4>; // Using default template arguments.

// Elements not initialized. See #configuration for more info.
Mat44 m1; 

// First row of the matrix is {1, 2, 3, 4}, regardless of Order or Layout.
Mat44 m2 = {
     1,  2,  3,  4,
     5,  6,  7,  8,
     9, 10, 11, 12,
    13, 14, 15, 16
};
```

Currently, matrices rows/columns be constructed from vectors, and there are no concatenation operations for matrices.


#### Converting matrices

Matrices can be readily converted via an explicit cast, when the multiplication order is the same:
```c++
// Note that both are FOLLOW_VECTOR
using Mat44FRd = Matrix<double, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat44FC_Packed = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;

Mat44FRd calculations = ...;
Mat44FC_Packed upload = Mat44FC_Packed(calculations);
```

When you want to convert between different orderings, you have to be more specific:
```c++
// Note that this one is PRECEDE_VECTOR
using Mat44PRd = Matrix<double, 4, 4, eMatrixOrder::PRECEDE_VECTOR, eMatrixLayout::ROW_MAJOR, false>;

// This one will correctly represent the same transform, but with different notation:
Mat44PRd calculations_p = matrix_representation_cast<Mat44PRd>(calculations);

// This one represents an incorrect transform, but you can do it nevertheless:
Mat44PRd wrong_p = matrix_reinterpret_cast<Mat44PRd>(calculations);
```

#### Element access

You need to index the matrix via rows and columns:
```c++
Mat33 m = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
};

// Use m(row, column)
float m01 = m(0, 1); // 2
float m10 = m(1, 0); // 4
```

Unfortunately, this is one of the few things that you cannot change in Mathter.

#### Arithmetic

You can use the operators \*, + and - to multiply, add and subtract matrices where the sizes of the operands are correct. You can also use \* and / to multiply and divide the matrices by scalars.

```c++
Matrix<float, 4, 1> u; // A column-vector.
Matrix<float, 1, 4> v; // A row-vector.
Matrix<float, 4, 4> o = u * v; // The outer product.
Matrix<float, 1, 1> i = v * u; // The inner-product.
```

You can also use compound operators when applicable:
```c++
Mat44 m, n;
m *= n;
m += n;
m *= 2.0f;
```


#### Matrix functions

Similarly to vector functions, matrix functions are also provided as free functions:

```c++
Mat44 m;

// All of U,S,V are 4x4 matrices.
auto [U,S,V] = DecomposeSVD(Determinant(m) * Inverse(Transpose(m)));
```


#### Submatrices (WARNING: deprecated)

There is an implementation of submatrices in the library, however, it awaits a major overhaul, and as such is deprecated and will be removed. If this doesn't bother you, it can be used.

An example:
```c++
Mat22 m = {
    1, 2,
    3, 4
};
Mat44 M = Zero(); // All elements zeroed.

// M will now look like
// 1, 2, 0, 0,
// 3, 4, 0, 0,
// 0, 0, 0, 0,
// 0, 0, 0, 0
M.Submatrix<2,2>(0,0) = m;
```


Quaternions
<a name="quaternions"></a>
---

The Quaternion class:
```c++
template <class T, bool Packed>
class Quaternion;
```

Quaternions only have the type and the packing parameters, which should by now be familiar to you. If not, then read the section about [vectors](#vectors) and [configuration](#configuration)


#### Creating quaternions

In general, when creating quaternions, you have to specify the scalar real and the imaginary vector components:
```c++
using Quat = Quaternion<float, false>;

Quat q1 = { 1, 0, 0, 0 }; // Unit quaternion, with the real component 1 and the imaginaries 0.
Quat q2 = { 1, Vec3(0,0,0) }; // Same as above.
Quat q3 = { Vec3(2,3,4) }; // The real components is zero, while the imaginaries are 2i + 3j + 4k.

using namespace quat_literals;
Quat q4 = 1 + 2_i + 3_j + 4_k; // This is not constexpr, so slow, but it's fun.
```

#### Element access

You can use the set {w,x,y,z} or {s,i,j,k}:
```c++
Quat q;
q.w = 1;
q.i = 0;
q.j = 0;
q.k = 0;
```

#### Arithmetic

The most important regarding quaternions is multiplication:
```c++
Quat rotation1;
Quat rotation2;
Quat rotationCombined = rotation2 * rotation1;
```

As you can see, the quaternions are multiplied in reverse order to chain rotations. Unlike the order of matrices, this cannot be configured in Mathter. The reason for this is that Mathter treats quaternions as mathematical objects in the first place, and rotators in the second. Quaternions are simply complex numbers with 3 imaginary elements as opposed to one, and the distributive law determines the multiplication's result.

There is not much to say about the other operations:
```c++
Quat q1, q2;
Quat q3 = (q1 + q2) * 3;
```

#### Quaternions functions

Similarly to Matrix and Vector functions, Quaternion functions are also provided as free functions. You can find things important for 3D, such as ```Normalize```, ```Length```, and ```Inverse```. These are also found by their mathematical names like ```Abs``` and ```Conjugate```. Additionally, there are mathematical functions like ```Exp``` and ```Log```.

#### Conversion to and from matrices

Since quaternions and matrices can both represent rotations, it makes sense to allow conversions between them. You can use typecasts for this purpose:

```c++
Quat q = RotationAxisAngle(Vec3(1, 0, 0), Deg2Rad(45.f));
Mat33 m = Mat33(q);
q = Quat(m);
```


Transforms
<a name="transforms"></a>
---

Arguably, the most important part of a 3D game math library is to do 3D transforms. In Mathter, you have the most common transform available to you as builder functions.

An example translation matrix:
```c++
Mat44 translation = Translation(1, 2, 3);
```

The transforms are lazily constructed. The builder functions returns a proxy object that can then be implicitly converted to any type that is suitable for representing said transform. Here is an example:
```c++
Mat33 r1 = RotationX(1.0f); // Okay, just a linear mapping.
Mat44 r2 = RotationX(1.0f); // Okay, may be an affine transform.
Mat43F r3 = RotationX(1.0f); // Only follow-vector!
Mat34P r4 = RotationX(1.0f); // Only precede-vector!
Quat q = RotationX(1.0f); // Equality is important, no discrimination against quaterions.
```

You can easily find the full list of transformations by browsing the sources.