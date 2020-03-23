Introduction
===


Structure of the library
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


Vectors
---

#### The Vector class 

Vectors have three template parameters:
```c++
template <class T, int Dim, bool Packed>
class Vector;
```

With ```T```, you can specify which underlying type the vector uses for its elements. You can set it to ```float```, ```double```, ```int```, ```std::complex``` or any other type or your preference. Floating point types are well-tested, integers also work, there is no reason others wouldn't work but they are not tested. Note, however, that some operations will fail to compile or provide incorrect results when they don't make sense. For example, you cannot normalize an integer vector because the result is not representable, you cannot ```Min/Max``` a complex vector because complex numbers have no ordering, and the length of an integer vector will get rounded to an integer, which is not what you generally want.

The dimensions specifies the number of elements in the vector. It can be any positive number, including 1. Dynamically sized vectors (and matrices) are not supported yet.

The ```Packed``` flag specifies whether the vector should ditch SIMD and forced alignment to tightly pack its elements. This is useful when you want to upload to the GPU and you needs strict layout requirements.


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

#### Element access

You can access the individual elements of the vectors by indices or coordinate axis names:
```c++
Vec3 v;
v.x;
v(0);
v[0]; // same as the above two
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

A similar issue arises when multi

