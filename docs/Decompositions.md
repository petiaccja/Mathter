# Matrix decompositions, equation systems

## Matrix decompositions in Mathter

Mathter implements the following matrix decompositions:
- LU and LUP: with and without pivoting
- QR and LQ
- SVD

You can use structured bindings to get the matrix factors:

```c++
const auto [Q, R] = DecomposeQR(A);
```

## Solving linear systems of equations

In addition to returning the matrix factors, all decompositions have a `Solve` method to solve linear systems of equations.

Let's consider the following equation system:
```
2x + 3y + 1z =  8
1x + 4y + 1z = -5
6x + 1y - 4z =  3
```

One could rewrite this in matrix form in the following way:
```
[[2, 3,  1],    [[x],    [[ 8],
 [1, 4,  1],  *  [y],  =  [-5],
 [6, 1, -4]]     [z]]     [ 3]]
```

The above form is very easy to express in Mathter. Note that the matrix must have precede-vector semantics as that's how the above formula is written:

```c++
const Matrix<float, 3, 3, eMatrixOrder::PRECEDE_VECTOR> A = {
    2, 3,  1,
    1, 4,  1,
    6, 1, -4
};
const Vector<float, 3> b = {8, -5, 3}
```

Regardless of the decomposition, you can solve it by this code:

```c++
const Vector<float, 3> x = DecomposeQR(A).Solve(b);
```

This should give you the vector `[x, y, z]^T` that satisfies the above equation.

## Solving least squares problems

In linear least squares problems, you have more equations than you have unknowns, so a least squares problem in matrix form would look like this:

```
[[2, 3],    [[x],    [[ 8],
 [1, 4],  *  [y]]  =  [-5],
 [6, 1]               [ 4],
 [9, 2]]              [ 3]]
```

You can follow the same method as above to solve the least squares problem for `x` and `y`. Note that this is only possible with the QR decomposition and the SVD, but not with LU/LUP.

## Inverse and pseudoinverse

In addition to `Solve`, all decompositions have an `Inverse` method as well. For square matrices, this returns the inverse, whereas for rectangular matrices it returns the Moore-Penrose pseudoinverse. Unless you want the pseudoinverse, you are best of relying on the usual `Inverse` method in the matrix mathematical functions header.

## To invert or not to invert

#### Small matrices

Raise your hands if you've heard a thousand times that you should not invert matrices. This is generally good advice, however it's usually given in the context of massive sparse matrices with millions of elements. For tiny matrices of 4x4 or smaller, which you most often use, Mathter uses explicit formulas to calculate the inverse. With these explicit formulas, inverting the matrix and multiplying by the inverse is much faster that even doing a lightweight LUP decomposition.

#### Mid-sized matrices 

For larger matrices, 5x5 and above, Mathter uses the LUP decomposition to compute the inverse. In this case, if you solve an equation system only once, you're probably better off with a direct solution as opposed to inverting.

#### Solving multiple equation systems at once

There are basically two options:
- Invert the matrix, and multiply all the load vectors by the inverted matrix
- Concatenate the load vectors into a large matrix, and pass the matrix to `Solve`

## Limitations

The major limitation is is that Mathter does not support dynamically sized matrices, so these methods are only applicable to fixed size problems. The algorithms are also tuned for small matrix sizes, and would be slow on much larger problems.




