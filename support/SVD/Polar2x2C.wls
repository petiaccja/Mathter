A = {{a11, a12}, {a21, a22}} (* The matrix of which the polar decomposition we want to compute. *)

U = {{u11, u12}, {u21, u22}}
P = {{p11, p12}, {Conjugate[p12], p22}}

(* Calculate the elements of UP. *)
UP = Simplify[U . P]

Print[UP]