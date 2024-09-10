# Vectors, matrices, and quaternions

## Vectors

While it's possible to use row or column matrices, Mathter comes with a dedicated vector type. This makes it easier to implement and use graphics-oriented features, such as element swizzling.

### Constructing vectors

Vectors have many different constructors to cover a wide range of use cases:

```c++
// Elements set to sNaN in debug, uninitialized in release.
Vector<float, 3> v1;

// [1, 2, 3, 1] Appends a 1 to the 3D vector to get 4D homogeneous coordinates.
Vector<float, 4> v2(Vector(1, 2, 3));

// [1, 2, 3] Truncates 4D homogeneous coordinates to 3D vectors by doing the perspective division.
Vector<float, 3> v3(Vector(2, 4, 6, 2));

// [1, 1, 1] All elements the same.
Vector<float, 3> v4(1.0f);

// [3, 2, 1, 1] Pass any swizzles, scalars, or vectors, they are concatenated, and CTAD figures out the element type, size, and packing mode.
Vector v(v2.zyx, 1.0f);
```

You can also convert between vectors of the same size but different scalar type or packing mode.

**Tip:** you can enable floating point exceptions to trap uninitialized vectors that have been initialized to signaling NaNs.

### Element access

You can use both the `[]` and `()` operators to access single elements:

```c++
Vector v(1, 2, 3);
v[1]; // == 2
v(2); // == 3
```

Additionally, vectors support iterators and the `data` method the same way as the STL:

```c++
Vector v(3, 1, 2);
std::sort(v.begin(), v.end());
std::ranges::sort(v);
```

Swizzlers are available for vectors with 4 or fewer elements, up to 4 elements:

```c++
Vector v(3, 1, 2);
auto a = v.xxxx; // [3, 3, 3, 3] The type is Swizzle<..., 0, 0, 0, 0>!
Vector b = v.xxxx; // [3, 3, 3, 3]
auto& b = v.z; // [2] The type is Swizzle<..., 0>&, not float&!
```

Keep it in mind that the single elements like `.x` are also swizzlers, not `float&`s as you would expect. This has a minimal effect in practice as they are implicitly convertible to `float`, they can participate in arithmetic as usual, and they can be assigned to.

**Good to know:** vectors are implemented as a union of swizzlers. Each swizzler has all elements of the vector, thus they have the exact same memory layout, which prevents UB, and makes it possible to carelessly copy and pass swizzlers around by value even after destroying the underlying vector object.

### Mathematical functions

Mathter implements the most common (`Length`, `Normalize`, `Dot`) and some of the less common (Gram-Schmidt orthogonalization) mathematical functions that work on vectors.

These are all implemented as free functions as opposed to being member functions, and they are found in `Vector/Math.hpp`. For a full list, you can browse that file and read the doc comments.

### Vector concatenation

**Possible deprecation**: this feature has become much less relevant now with CTAD. I recommend you use the vector constructor with CTAD as I might retire this in favor of bitwise operations.

It's possible to concatenate vector using the pipe operator:

```c++
Vector a(1, 2);
Vector b(3, 4);
auto c = a | b; // [1, 2, 3, 4]
Vector c(a, b); // Prefer this instead of the pipe operator.
```

## Matrices

Aside from the many template parameters that are explained in the gettings started part, matrices in Mathter are pretty standard.

### Constructing matrices

Matrices can be created from their elements, or for row and column matrices, from vectors:

```c++
// This will always give you the matrix as you see it below, regardless
// or the matrix's multiplication order or memory layout.
Matrix<float, 2, 3> m = {
    1, 2, 3,
    4, 5, 6
};

// You can convert a vector to a column (or row) matrix.
Matrix<float, 3, 1> r = Vector(1.0f, 2.0f, 3.0f);
```

Additionally, matrices of different scalar type, order, and layout can be converted to one another via the constructors. This will always preserve the transform, so, for example, converting a left-multiplying matrix to a right-multiplying matrix will transpose its elements:

```c++
// Affine transform with rotation, scaling, and translation.
Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> p = {...};

// Elements transposed, the shape is also transposed.
Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR> f(p);
```

### Element access

You can access individual elements, rows, columns, and submatrices:

```c++
Matrix<float, 3, 4> m;

// Get/set individual element.
m(1, 2) = 3.0f;

// Get/set row (or column).
Vector<float, 3> row0 = m.Row(0);
m.Row(1, row0);

// Get/set submatrix.
auto sm = m.Extract<3, 3>(0, 0);
sm *= Matrix<float, 3, 3>(Scale(1, 2, 1));
m.Insert(0, 0, sm);
```

### Mathematical functions

Mathter implements common math functions for matrices such as `Determinant` or `Inverse`. For the full list, you can check out `Matrix/Math.hpp`. Just like vectors, these are implemented as free functions.


## Quaternions

### Constructing quaternions

You can construct quaternions from its elements or from a scalar and vector part:

```c++
const auto [x, y, z] = ...;

// Specify elements one by one.
Quaternion<float> q(s, x, y, z);

// Specify scalar part and vector part.
Quaternion<float> q(s, Vector(x, y, z));

// Construct from vector part only. Useful when representing a regular 3D vector as a quaternion.
Quaternion<float> q(Vector(x, y, z));
```

Quaternions of different types and layout can be converted to each other. Additionally, quaternions can also be converted to and from 3x3 rotation matrices. Be careful with this, the behaviour is undefined if the matrix is not actually a rotation matrix.

```c++
Matrix<float, 3, 3> m = RotationX(1.0f);
auto q = Qauternion<float>(m);
```

### Element access

Quaternion elements are accessed via swizzlers, similarly to vectors:

```c++
Quaternion<float> q = RotationAxisAngle(Normalize(Vector(1.0f, 2.0f, 3.0f)), 0.7f);
auto vector = Vector(q.vector);
auto scalar = float(q.scalar);
```

Note that the quaternion's data layout is irrelevant, `q.x` returns the first element of the vector part regardless of whether the quaternion is stored in vector first or scalar first memory layout.

### Mathematical functions

Like vectors and matrices, quaternion also come with a useful set of mathematical functions, defined in the header file `Quaternion/Math.hpp`.

### Literals

You can also use literals to define quaternions:

```c++
using namespace mathter::quat_literals;
Quaternion<float> q = 1.0f + 2.0_if + 3.0_jf + 4.0_kf;
```

I'm unsure about the performance of this. The expressions are not `constexpr`, so compile time evaluation is not guaranteed, but if the compiler can do constant propagation through the SIMD intrinsics then this should be basically noop.

## Arithmetic

Vectors, matrices, and quaternions overload the arithmetic operators so that you can write concise code. It's **important** to remember that Mathter always uses mathematically correct notation, and does not redefine the meaning or order of operations to better align with graphics programming.

### Arithmetic on vectors

There is nothing notable for vector arithmetic, as they implement the usual elementwise operations. When doing arithmetic on mixed vectors and scalars, the scalars are implicitly extended to full vector length.

### Arithmetic on swizzlers

The vector arithmetic operators are overload for swizzlers as well, thus, you can immediately use them in arithmetic formulas:

```c++
auto v = a.xxx + a.yyy + a.zzz; // Returns a Vector<T, 3>
```

However, vector mathematical functions and vector-matrix multiplication are **not** overloaded for swizzlers, therefore you have to first convert them to a vector:

```c++
auto n = Normalize(Vector(v.xyz));
const auto w = M * Vector(v.xyz);
```

### Arithmetic on matrices

- Must be the same `Order`: you cannot mix matrices with follow-vector and precede-vector semantics in arithmetic code.
- Shapes must match: you can multiply a 2x3 and a 3x4 matrix, the result will be a 2x4 matrix.
- Vectors must be on the correct side: matrices with precede-vector semantics must multiply vectors from the right, as in `M*v`. The operators are not overloaded for the incorrect order. (Same for follow-vector order.)
- Matrix multiplication is as per mathematics, regardless of `Order`: you concatenate transforms by matrix multiplication, but the multiplication order is the opposite for follow-vector and precede-vector semantics.

For vector-matrix multiplication, there are some extras to aid graphics programming:
- You can multiply a 3-vector by a 4x4 matrix: this implicitly augments the vector to 4D homogeneous coordinates, does the multiplication, the perspective division, and returns the result in 3D space.
- This also applies for affine transforms represented by 3x4 or 4x3 matrices, but these don't require perspective divisions

While matrix layout does not affect the mathematics, it does affect performance in some cases. Column major matrices generally go better with preceding the vector, and row major matrices go better with following the vector. Multiplying matrices of the same layout is usually faster. Doing it the other way is still pretty fast, but may not be optimal.

### Arithmetic on quaternions

Quaternions are just complex numbers on steroids, and they are implemented as such in Mathter as opposed to being implemented as a special rotation class. While most things are as expected, two peculiarities are worth mentioning:

- Combining quaternion rotations: quaternions are always combined right-to-left, that is, `q = q3*q2*q1` will apply the rotation `q1` first. Rotations are not commutative, so the order is important. This is in contrast to matrices that may be combined as `R1*R2*R3` or `R3*R2*R1` depending on their multiplication order.
- Applying qauternions to vectors: currently, the `*` operator is overloaded for this purpose, and both `v2 = q*v` and `v2 = v*q` apply the quaternion rotation to the vector. This syntax may be confusing, so I'm considering **deprecating** this in favour of `q(v)` or `q.Apply(v)`.



