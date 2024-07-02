A = {{a11, a12}, {a21, a22}} (* The matrix of which the RQ decomposition we want to compute. *)

R = {{r11, r12}, {0, r22}} (* Define R as an upper triangular matrix. *)
Q = {{cq, -sq}, {sq, cq}} (* Define Q as a rotation matrix, where q is the rotation angle. *)

(* Expand the elements of R*Q. *)
RQ = R.Q

(* Using the 4 equations given by RQ=A and cos(q)^2 + sin(q)^2 == 1, solve for
the 5 unknown elements or R and Q. *)
solution = Simplify[Solve[RQ == A && sq^2 + cq^2 == 1, {cq, sq, r11, r12, r22}]]

(* Can choose either solution. There seems to be no numerical difference. *)
Print[solution[[1]]]
Print[solution[[2]]]
