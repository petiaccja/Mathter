(* Example for applying a givens rotation. *)

A = {
    {a11, a12, a13, a14},
    {a21, a22, a23, a24},
    {a31, a32, a33, a34},
    {a41, a42, a43, a44}
}

V = {
    {cv, 0, -sv, 0},
    {0, 1, 0, 0},
    {sv, 0, cv, 0},
    {0, 0, 0, 1}
}

Print[Simplify[A.V]]
Print[Transpose[A].A]