Define $A$ as:

$$
\begin{equation}
A = 
\begin{bmatrix}
a_{11} & a_{off} \\
a_{off} & a_{22}
\end{bmatrix}
\end{equation}
$$

Define $V$ as:

$$
\begin{equation}
A = 
\begin{bmatrix}
c_v & -s_v \\
s_v & c_v
\end{bmatrix}
\end{equation}
$$


Solutions to $V^T A V$:

$$
\begin{equation}
\begin{cases}
    z = a_{11} - a_{22} \\
    p_1 = \sqrt{z^2 + 2 a_{off}^2} \\
    p_a = \sqrt{(z + p_1)/p_1} \\
    p_b = \sqrt{(-z + p_1)/p_1} \\
    \begin{cases}
        c_v = -\frac{p_a}{\sqrt{2}} \\
        s_v = -\frac{\sqrt{2} a_{off}}{p_1 p_a}
    \end{cases} \\
    \begin{cases}
        c_v = \frac{p_a}{\sqrt{2}} \\
        s_v = \frac{\sqrt{2} a_{off}}{p_1 p_a}
    \end{cases} \\
    \begin{cases}
        c_v = -\frac{p_b}{\sqrt{2}} \\
        s_v = \frac{\sqrt{2} a_{off}}{p_1 p_b}
    \end{cases} \\
    \begin{cases}
        c_v = \frac{p_b}{\sqrt{2}} \\
        s_v = -\frac{\sqrt{2} a_{off}}{p_1 p_b}
    \end{cases}
\end{cases}
\end{equation}
$$

Numerically stable solution:

$$
\begin{equation}
\begin{cases}
    z = a_{11} - a_{22} \\
    p_1 = \sqrt{z^2 + 2 a_{off}^2} \\
    p_{a/b} = \sqrt{(|z| + p_1)/p_1} \\
    \begin{cases}
        c_v = -\frac{p_{a/b}}{\sqrt{2}} \\
        s_v = -\text{sign}(z)\frac{\sqrt{2} a_{off}}{p_1 p_{a/b}}
    \end{cases}
\end{cases}
\end{equation}
$$