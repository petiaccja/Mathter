R = {{r11, r12}, {0, r22}} (* The matrix of which the SVD we want to compute. *)

(* V diagonalizes A such that V * A^T * A * V^T = D is diagonal. *)
V = {{cv, -Conjugate[sv]}, {sv, Conjugate[cv]}}

(* Calculate the elements of D. *)
D1 = Simplify[V . Conjugate[Transpose[R]] . R . Conjugate[Transpose[V]]]

(* Solve off-diagonals to be zero. *)

(* u = cv / Conjugate[sv] *)
eq1 = Simplify[Replace[D1[[1, 2]], {sv -> Conjugate[cv / u]}, All]]
eq2 = Simplify[Replace[D1[[1, 2]], {cv -> u * Conjugate[sv]}, All]]
Print[Simplify[Replace[Solve[eq1 == 0, u, Complexes], {cv -> 1}, All]]]
Print[Simplify[Replace[Solve[eq2 == 0, u, Complexes], {sv -> 1}, All]]]

(* v = Conjugate[sv] / cv *)
eq3 = Simplify[Replace[D1[[1, 2]], {sv -> Conjugate[v * cv]}, All]]
eq4 = Simplify[Replace[D1[[1, 2]], {cv -> Conjugate[sv] / v}, All]]
Print[Simplify[Replace[Solve[eq3 == 0, v, Complexes], {cv -> 1}, All]]]
Print[Simplify[Replace[Solve[eq4 == 0, v, Complexes], {sv -> 1}, All]]]
