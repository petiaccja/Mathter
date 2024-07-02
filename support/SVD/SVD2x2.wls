A = {{a11, a12}, {a21, a22}} (* The matrix of which the SVD we want to compute. *)

(* V diagonalizes A such that V * A^T * A * V^T = D is diagonal. *)
cv = Cos[v] (* v is the rotation angle of the rotation matrix V. *)
sv = Sin[v]
tv = Tan[v]
V = {{cv, -sv}, {sv, cv}}

(* Calculate the elements of D. *)
D1 = V . Transpose[A] . A . Transpose[V]

(* Solve off-diagonals to be zero. *)
trigRelation := cv^2 + sv^2 == 1
offdiagonalRelation := D1[[1, 2]] == 0
cvSv = FullSimplify[Solve[trigRelation && offdiagonalRelation, {cv, sv}]]

(* Simplified case for RQ preconditioning: a21 = 0. *)
cvSvRq = FullSimplify[Replace[cvSv, a21 -> 0, All]]

Print["Solutions without preconditioning:"]
Print[cvSv]
Print["Solutions with RQ preconditioning:"]
Print[cvSvRq]