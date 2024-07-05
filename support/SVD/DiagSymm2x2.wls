A = {{a11, aoff}, {aoff, a22}} (* The matrix of which the SVD we want to compute. *)

(* V diagonalizes A such that V * A^T * A * V^T = D is diagonal. *)
cv = Cos[v] (* v is the rotation angle of the rotation matrix V. *)
sv = Sin[v]
tv = Tan[v]
V = {{cv, -sv}, {sv, cv}}

(* Calculate the elements of D. *)
D1 = Transpose[V] . A . V

(* Solve off-diagonals to be zero. *)
trigRelation := cv^2 + sv^2 == 1
offdiagonalRelation := D1[[1, 2]] == 0
cvSv = FullSimplify[Solve[trigRelation && offdiagonalRelation, {cv, sv}]]

Print["Solutions:"]
Print[cvSv]