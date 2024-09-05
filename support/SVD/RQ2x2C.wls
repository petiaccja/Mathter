A = {{a11, a12}, {a21, a22}} (* The matrix of which the RQ decomposition we want to compute. *)
R = {{r11, r12}, {0, r22}} (* Define R as an upper triangular matrix. *)
Q = {{cq, -Conjugate[sq]}, {sq, Conjugate[cq]}} (* Define Q as a rotation matrix, where q is the rotation angle. *)

(* Expand the elements of R*Q. *)
RQ = Simplify[R . Q]

Print[RQ]