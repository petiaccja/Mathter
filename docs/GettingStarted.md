# Gettings started

## Library structure

Mathter's code is grouped by its features:
- **Common**: this is internal to Mathter, you shouldn't use it, though you can if you want to.
- **Decompositions**: methods for matrix decompositions. Include headers as needed.
- **Geometry**: geometric primitives (i.e. Line) and intersections. Include them via the top-level `Geometry.hpp`.
- `IoStream.hpp`: writing and reading math types to and from standard I/O streams.
- **Matrix**: matrices, arithmetic, and functions. Include the top-level `Matrix.hpp`.
- **Quaternion**: quaternions, arithmetic, and functions. Include the top-level `Quaternion.hpp`.
- **Transforms**: all the transforms, such as translation and rotation. Include the top-level `Transforms.hpp`.
- `Utility.hpp`: basic utilities, like conversion between angles and radians.
- **Vector**: vectors, arithmetic, and functions. Include the top-level `Vector.hpp`.

To get a glance of the available features, you can browse through the folder structure. You should be familiar with the contents based on the file names.

I suggest that you use the modules of Mathter through the top-level header files, and not include files deeper in the folder structure.

## Template parameters explained

Mathter is highly configurable via template parameters. For `Vector`s, `Matrix`es, and `Quaternion`s, the following template parameters are available to you:
- `T`: the type of the element of a vector, matrix, or quaternion. For vectors and matrices, integers, floating point, and standard complex types are supported. For quaternions, only floating point types are supported. Note that while you can instantiate the math types with many scalar types, some operations won't make sense. For example, you cannot normalize a vector of integers, because the result cannot be represented.
- `Dim`, `Rows`, `Columns`: specifies the size of vectors or matrices. The sizes are fixed at compilation time as Mathter does not support dynamic sizes. There are no restrictions on the size at compilation time, so you can create a 16x4 matrix or a vector of length 11.
- `Order`: defines whether matrices multiply vectors from the left or right. The order also affects how you concatenate matrix transforms. There are two options:
    - `PRECEDE_VECTOR`:
        - Multiplying vectors from the left: `M*v`.
        - Combining transforms right-to-left: `M3*M2*M1` applies `M1` first.
    - `FOLLOW_VECTOR`:
        - Multiplying vectors from the right: `v*M`.
        - Combining transforms left-to-right: `M1*M2*M3` applies `M1` first.
- `Layout`: 
    - For matrices: selects either row-major or column-major memory layout. The memory layout does not affect calculations and the interface in any way, but it might affect performance.
    - For quaternions: stores either the vector part or the scalar part first. This does not affect calculations and the interface, and is also irrelevant for performance.
- `Packed`: when packing is turned off, Mathter will use SIMD and may pad or overalign objects to achieve maximum performance. When packing is off, the elements of math types are tightly packed and Mathter does not explicitely use SIMD (though the compiler might). Use packed types when memcpy'ing to the GPU, otherwise use non-packed types to improve performance.


## Integrating Mathter into your project

Mathter's templates are meant to give you full flexibility in choosing your math conventions, however, they can be verbose to type every time. It's best if you define your own types via simple typedefs in a header file accessible throughout your entire project:

```c++
#pragma once

#include <Mathter/*.hpp>

using namespace mathter;

using Vec3 = Vector<float, 3, false>;
using Mat44 = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Quat = Quaternion<float, eQuatLayout::SCALAR_FIRST, false>;

using Vec3P = Vector<float, 3, true>;
using Mat44P = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;
using QuatP = Quaternion<float, eQuatLayout::VECTOR_FIRST, true>;
```

If you need to interface with other applications or hardware, you can define multiple sets of primitives. In the example code above, there is an SIMD set for CPU calculations, and there is a packed set for GPU uploads and downloads. Keep in mind that the multiplication order or memory layout does not have to match, Mathter will convert it seamlessly while preserving the meaning of the transforms.

For methods with a lot of dynamic configuration options, such as perspective transforms, it may be a good idea to define wrapper functions:

```c++
auto MyPerspective(float fovX,
                   float aspectRation,
                   float nearPlane,
                   float farPlane) {
    // You always configure NDC space like this in your graphics API:
    constexpr float projNearPlane = 0.0f;
    constexpr float projFarPlane = 1.0f;
    return Perspective(fovX,
                       aspectRatio,
                       nearPlane,
                       farPlane,
                       projNearPlane,
                       projFarPlane);
}
```