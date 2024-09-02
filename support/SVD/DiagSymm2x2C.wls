A = {{a11, aoff}, {Conjugate[aoff], a22}} (* The matrix of which the SVD we want to compute. *)

(* V diagonalizes A such that V * A^T * A * V^T = D is diagonal. *)
LA = L * IdentityMatrix[2] - A;
v = {{cv}, {sv}}

(* Calculate eigenvalues. *)
det = Det[LA]
Print[Simplify[Solve[det == 0, {L}]]]

(* Calculate the elements of D. *)
x = (LA) . v

solSv1 = Solve[x[[1, 1]] == 0, {sv}]
solSv2 = Solve[x[[2, 1]] == 0, {sv}]

Print[solSv1]
Print[solSv2]


