# Transforms

Mathter implements common geometric transforms as a set of lazy-evaluated free functions that are found in the `Transforms/` folder. I recommend browsing the sources to discover all available transforms. If you need to figure out how a particular transform is created, the docstrings will help you, so I won't repeat that in this guide.

## Anatomy of a transform

In Mathter, transforms don't return a matrix or a quaternion, but a proxy object which I refer to as a *builder*:

```c++
// Rotate around the X axis by 1.7 radians.
const auto builder = RotationX(1.7f);
```

You shouldn't interact with this builder object, nor can you use it in code to rotate a vector. You must convert it to an actual mathematical primitive:

```c++
// The quaternion and the matrix represent the same rotation.
Quaternion<float> q = builder;
Matrix<float, 3, 3> m = builder;
```

## Using transforms

Although you can, typically, you will not explicitly hold the builder object around, but use one of the following patterns:

```c++
auto q = Quaternion<float>(RotationX(1.7f));
Matrix<float, 3, 3> m = RotationX(1.7f);
```

As you can see, quaternions and rotation matrices are treated with the same syntax. In general, transforms can be converted to any mathematical primitive that can represent that transform. To give you an example, a 3x4 matrix with precede-vector order can hold a 3D translation, and any vector, matrix, or quaternion can hold a `Zero()` transform.

