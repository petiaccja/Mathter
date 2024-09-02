## Computing the 2x2 singular value decomposition

This document takes a practical approach. The theory and proofs for the existence of solutions and their properties can be found on [Wikipedia's article on SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition#Based_on_the_spectral_theorem). This document assumes real matrices.


### Diagonalising the 2x2 matrix

The first step is to compute a rotation matrix $V$ that diagonalizes $A$:
$$
VA^TAV^T = D,\\
D\text{ is a diagonal matrix, and}\\
V\text{ is a rotation matrix}.
$$

Let's define $A$ and $V$ in terms of their elements:
$$
\begin{equation}
A = 
\begin{bmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{bmatrix}
V = 
\begin{bmatrix}
c_v & -s_v \\
s_v & c_v
\end{bmatrix}
\end{equation}
$$

Note:
- $c_v$ and $s_v$ refer to $cos(v)$ and $sin(v)$, where $v$ is the rotation angle of $V$. This solution does not use any trigonometry, hence the treatment of the elements as simple variables.
- $V$ is a rotation matrix (SO(2)) with a determinant of $+1$. Alternatively, $V$ could also be a non-special orthogonal matrix (O(2)) with a determinant of $\pm1$, which allows the negation of the singular values.

After performing the matrix multiplications, the elements of $VA^TAV^T$ are as follows:
$$
\begin{equation}
D = VA^TAV^T = 
\begin{bmatrix}
d_{11}(a_{ij}, c_v, s_v) & d_{12}(a_{ij}, c_v, s_v)\\
d_{21}(a_{ij}, c_v, s_v) & d_{22}(a_{ij}, c_v, s_v)
\end{bmatrix}
\end{equation}
$$

(Each $d_{ij}$ is a function of $a_{ij}$, $c_v$, and $s_v$. The expanded formulas for $d_{ij}$ are too long to fit on this page.)

Based on the expanded formulas, the off-diagonal elements of $D$, $d_{12}$ and $d_{21}$, are equal. If $D$ is diagonal, these off-diagonals are zero, so we have to find $c_v$ and $s_v$ that ensures this. The single equation for $0 = d_{12} = d_{21}$ can be supplemented by the trigonometric relation between the cosine and sine, giving two equations for two unknowns:

$$
\begin{equation}
\begin{cases}
0 = d_{12}\\
1 = c_v^2 + s_v^2 = cos(v)^2 + sin(v)^2
\end{cases}
\end{equation}
$$

Solving this system of equations for $c_v$ and $s_v$ gives the following four solutions:

$$
\begin{equation}
	\begin{cases}
		z = a_{11}^2 - a_{12}^2 + a_{21}^2 - a_{22}^2 \\
		p_1 = \sqrt{((a_{12} + a_{21})^2 + (a_{11} - a_{22})^2) \cdot ((a_{12} - a_{21})^2 + (a_{11} + a_{22})^2)} \\
		p_a = \sqrt{(z + p_1)/p_1} \\
		p_b = \sqrt{(-z + p_1)/p_1} \\
		\begin{cases}
			c_v = -\frac{p_a}{\sqrt{2}} \\
			s_v = \frac{(-z + p_1) \cdot p_a}{2\sqrt{2}(a_{11}a_{12} + a_{21}a_{22})}
		\end{cases} \\
		\begin{cases}
			c_v = \frac{p_a}{\sqrt{2}} \\
			s_v = \frac{(z - p_1) \cdot p_a}{2\sqrt{2}(a_{11}a_{12} + a_{21}a_{22})}
		\end{cases} \\
		\begin{cases}
			c_v = -\frac{p_b}{\sqrt{2}}, \\
			s_v = \frac{(-z - p_1) \cdot p_b}{2\sqrt{2}(a_{11}a_{12} + a_{21}a_{22})}
		\end{cases} \\
		\begin{cases}
			c_v = \frac{p_b}{\sqrt{2}}, \\
			s_v = \frac{(z + p_1) \cdot p_b}{2\sqrt{2}(a_{11}a_{12} + a_{21}a_{22})}
		\end{cases}
	\end{cases}
\end{equation}
$$

### Fixing the divisor with RQ preconditioning

Consider the following part of the divisor for $s_v$:
$$
\begin{equation}
a_{11}a_{12} + a_{21}a_{22}
\end{equation}
$$

The numerical evaluation with 23-bit numbers (such as FP32) yields a 46-bit results for both $a_{11}a_{12}$ and $a_{21}a_{22}$, however, the least significant 23 bits are discarded due to floating point rounding. In the scenario where, for example, the first 21 bits of the 46 bit results are the same, and the results are of opposite sign, the addition will cancel out those 21 bits, leaving us with only 2 bits of precision. As this imprecise factor contributes multiplicatively to $s_v$, it can render the results unusable.

To fix this, we can use RQ decomposition to split $A$ into an upper-triangular matrix and a rotation matrix:

$$
\begin{equation}
A = RQ = 
\begin{bmatrix}
r_{11} & r_{12} \\
r_{21}=0 & r_{22}
\end{bmatrix}
\begin{bmatrix}
c_q & -s_q \\
s_q & c_q
\end{bmatrix}
\end{equation}
$$

After the RQ decomposition, the SVD of $R$ can be calculated, and we can restore $A$ using $Q$:

$$
\begin{equation}
\begin{split}
A = USV \\
RQ = USV \\
RQ = USV'Q \\
V = V'Q
\end{split}
\end{equation}
$$

The original $V$ matrix in the SVD of $A$ can also be calculated using $Q$, as shown in the equation above.

The divisor after the RQ decomposition reads:

$$
\begin{equation}
\begin{split}
r_{21} = 0 \\
r_{11}r_{12} + r_{21}r_{22} = r_{11}r_{12}
\end{split}
\end{equation}
$$

The catastrophic cancellation of significant bits is gone, plus we can also simplify all the formulas above by replacing $a_{21}$ with zeroes.

*It may be possible to use only the cosine ($c_v$) in which the problematic divisor is not present, and calculate the sine and its sign with other methods. My original code uses the tangent, in which the problematic divisor is present, so RQ decomposition seemed the best option. I'll continue with RQ decomposition as it's simple, accurate, and substantially simplifies all the above formulas.*

**In the rest of the document, I'm assuming an upper-triangular $A$**, for which the four simplified solutions are as follows:

$$
\begin{equation}
\begin{cases}
z = a_{11}^2 - a_{12}^2 - a_{22}^2 \\
p_1 = \sqrt{(a_{12}^2 + (a_{11} - a_{22})^2) \cdot (a_{12}^2 + (a_{11} + a_{22})^2)} \\
p_a = \sqrt{(z + p_1)/p_1} \\
p_b = \sqrt{(-z + p_1)/p_1} \\
	\begin{cases}
		c_v = -\frac{p_a}{\sqrt{2}} \\
		s_v = \frac{\sqrt{2} a_{11} a_{12}}{p_1 p_a}
	\end{cases} \\
	\begin{cases}
		c_v = \frac{p_a}{\sqrt{2}} \\
		s_v = -\frac{\sqrt{2} a_{11} a_{12}}{p_1 p_a}
	\end{cases} \\
	\begin{cases}
		c_v = -\frac{p_b}{\sqrt{2}} \\
		s_v = -\frac{\sqrt{2} a_{11} a_{12}}{p_1 p_b}
	\end{cases} \\
	\begin{cases}
		c_v = \frac{p_b}{\sqrt{2}} \\
		s_v = \frac{\sqrt{2} a_{11} a_{12}}{p_1 p_b}
	\end{cases}
\end{cases}
\end{equation}
$$

### Fixing $p_a$ and $p_b$ by selecting the better of the four solutions

Both $p_a$ and $p_b$ contain the term $\pm z + p_1$. This term can eventually be brought to a form $\pm z + \sqrt{z^2 + 4g^2}$ by rearranging the terms. (See the next section for explanation of $g$.) This is prone to the same cancellation effects that made us introduce the RQ decomposition. Think about the case when $+z$ is negative and $g$ is much smaller than $z$. In this case, the value of $z$ and $\sqrt{z^2 + 4g^2}$ will be very close, but the opposite sign, effectively cancelling out most of the most significant bits. As this term also multiplicatively contributes to the final results, the algorithm will fail for some inputs.

To avoid this issue, we can dynamically choose between employing $p_a$ (thus selecting either of the first two solutions), or employing $p_b$ (thus selecting either of the latter two solutions). Since the formulas for the first two and latter two solutions are very similar, we can easily switch between them using $\text{sign}(z)$.

### Fixing singularity when $p_1$ is close to zero

Look at the following equation for $p_1$:
$$
\begin{equation}
p_1 = \sqrt{(a_{12}^2 + (a_{11} - a_{22})^2) \cdot (a_{12}^2 + (a_{11} + a_{22})^2)}
\end{equation}
$$

We can extract the terms under the square root like so:
$$
\begin{equation}
\begin{split}
f_1 = a_{12}^2 + (a_{11} - a_{22})^2 \\
f_2 = a_{12}^2 + (a_{11} + a_{22})^2 \\
p_1 = \sqrt{f_1 \cdot f_2}
\end{split}
\end{equation}
$$

Consider the scenario when $a_{11} = a_{22} = 0.5$, and $a_{12}$ is tiny ($\simeq 10^{-19}$). In this case, $f_1 = a_{12}^2$, which is a denormal floating point number at 32 bits of precision, while $f_2$ will simply be 1 as the denormal part in the sum vanishes. The low precision denormal is propagated multiplicatively into $s_v$, corrupting the results. While $c_v$ seems to be immune to this problem as $p_1 >> z$, it's probably better to fix this as opposed to calculating $s_v$ differently.

Let's rearrange the terms of $p_1$ in a more floating-point-friendly way:
$$
\begin{equation}
\begin{split}
z = a_{11}^2 - a_{12}^2 - a_{22}^2 \\
g = a_{11} a_{12} \\
p_1 = \sqrt{z^2 + 4g^2}
\end{split}
\end{equation}
$$

Repeating the short calculation above with the same numbers, $z$ is still denormal, however, $g$ has a value of `0.5e-19f`. When we feed this into the square root, $g >> z$, meaning the denormal's digits have no significance in $p_1$. To avoid getting a denormal for $g^2$, `std::hypot` can be used.

### Fixing inaccuracies of $z$

The formula $1 - x^2$, when evaluated using floating point numbers, has the largest error when $1 - x \simeq \sqrt{\epsilon}$. Consider the value $x = 0.9999$, in which case the accurate value after extending to FP64 is 0.00020002318454714896, however, using 32-bit precision the result is 0.000200033188, meaning the error is around $420 \cdot \epsilon$. If we rearrange the formula as $(1 - x)(1 + x)$, we get the results almost as accurately as FP32 can represent it: 0.000200023191.

Using this knowledge, we can rearrange $z$ as follows:

$$
\begin{equation}
z = (a_{11} - a_{22})(a_{11} + a_{22}) - a_{12}^2
\end{equation}\\
\begin{equation}
z = (a_{11} - a_{12})(a_{11} + a_{12}) - a_{22}^2
\end{equation}
$$

When we choose the first rearrangement, either $g = a_{11} a_{12}$ or $z$ will be much larger in magnitude than the error we introduce when $a_{11}$ and $a_{12}$ are close, therefore the error won't show up in $p_1 = \sqrt{z^2 + 4g^2}$. On the other hand, with the second rearrangement, choosing $a_{12} = 0$ will return us to the degenerate case of $a_{11}^2 - a_{22}^2$ while $g = 0$, meaning the entire error of the calculation flows into $p_1 = \sqrt{z^2 + 4g^2}$.

In light of this, we can safely use the first rearrangement at all times, because its inaccuracies will be shadowed by much larger terms in the calculations.

### Putting all fixes together

After applying the fixes above, the solution looks like this:

$$
\begin{equation}
\begin{cases}
z = (a_{11} - a_{22})(a_{11} - a_{22}) - a_{12}^2 \\
g = a_{11} a_{12} \\
p_1 = \sqrt{z^2 + 4g^2} \\
p_{a/b} = \sqrt{(|z| + p_1)/p_1} \\
c_v = -\frac{p_{a/b}}{\sqrt{2}} \\
s_v = \text{sign}(z)\cdot\frac{\sqrt{2} g}{p_1 p_{a/b}}
\end{cases}
\end{equation}
$$

### Computing $U$ and $S$ given $A$ and $V$

Using the definition of the SVD, we can compute $US$ as follows:

$$
\begin{equation}
\begin{split}
A = USV \\
AV^{-1} = USVV^{-1} \\
AV^{-1} = US \\
\text{Since }V\text{ is a rotation matrix:}\\
AV^T = US \\
US = 
\begin{bmatrix}
\mu_{11} & \mu_{12} \\
\mu_{21} & \mu_{22}
\end{bmatrix} =
\begin{bmatrix}
a_{11} & a_{12} \\
0 & a{22}
\end{bmatrix}
\begin{bmatrix}
c_v & s_v \\
-s_v & c_v
\end{bmatrix} = \\ =
\begin{bmatrix}
c_v a_{11} - s_v a_{12} & s_v a_{11} + c_v a_{12} \\
-s_v a_{22} & c_v a_{22}
\end{bmatrix}
\end{split}
\end{equation}
$$

Since we know that $S$ is a diagonal matrix, we can express $US$ like this as well:

$$
\begin{equation}
\begin{split}
US = 
\begin{bmatrix}
\mu_{11} & \mu_{12} \\
\mu_{21} & \mu_{22}
\end{bmatrix} =
\begin{bmatrix}
c_u & -s_u \\
s_u & c_u
\end{bmatrix}
\begin{bmatrix}
s_{11} & 0 \\
0 & s_{22}
\end{bmatrix} = \\ =
\begin{bmatrix}
 c_u s_{11} & -s_u s_{22} \\
s_u s_{11} & c_u s_{22}
\end{bmatrix}
\end{split}
\end{equation}
$$

We can derive the formula for $s_{ii}$ by starting with the well-known trigonometric identity:
$$
\begin{equation}
\begin{split}
c_u^2 + s_u^2 = 1 \\
(c_u s_{11})^2 + (s_u s_{11})^2 = s_{11}^2 \\
\mu_{11}^2 + \mu_{21}^2 = s_{11}^2 \\
\pm\sqrt{\mu_{11}^2 + \mu_{21}^2} = s_{11}
\end{split}
\end{equation}
$$

Calculating $c_u$ and $s_u$ follows trivially:

$$
\begin{equation}
\begin{split}
c_u = \mu_{11} / s_{11} \\
s_u = \mu_{21} / s_{11}
\end{split}
\end{equation}
$$

It's **important** to note that $U$ has a determinant of $+1$ as we defined it. This means that if $US$ has a negative determinant, one has to make one of $s_{ii}$ negative before calculating $c_u$ and $s_u$. Alternatively, if positive singular values are required, the proper row/column of $U$ can be sign-flipped so that $U$ has a determinant of $-1$. This would mean $U$ is no longer a pure rotation matrix (SO(2)), but a rotoreflection matrix (O(2)).

### Computing the RQ decomposition

The RQ decomposition rewrites the matrix $A$ as a product of an upper-triangular matrix $R$ and a rotation matrix $Q$:

$$
\begin{equation}
A = RQ
\end{equation}
$$

*It's the same as the well-known QR decomposition, just flipped so that it can be used for the SVD. $R$ could also be lower triangular, it makes no difference.*

We can find the solution by expanding the matrix product:
$$
\begin{equation}
\begin{split}
A = 
\begin{bmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{bmatrix} = 
\begin{bmatrix}
r_{11} & r_{12} \\
0 & r_{22}
\end{bmatrix}
\begin{bmatrix}
c_q & -s_q \\
s_q & c_q
\end{bmatrix} = 
RQ \\
\begin{bmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{bmatrix} = 
\begin{bmatrix}
c_q r_{11} + s_q r_{12} & -s_q r_{11} + c_q r_{12}\\
s_q r_{22} & c_q r_{22}
\end{bmatrix}
\end{split}
\end{equation}
$$

This and the trigonometric relation between the cosine and the sine gives us 5 equations for the 5 unknowns:

$$
\begin{equation}
\begin{cases}
1 = c_q^2 + s_q^2 \\
a_{11} = c_q r_{11} + s_q r_{12} \\
a_{12} = -s_q r_{11} + c_q r_{12}\\
a_{21} = s_q r_{22} \\
a_{22} = c_q r_{22}
\end{cases}
\end{equation}
$$

Solving for the unknowns gives us two solutions for the elements of $R$ and $Q$:

$$
\begin{equation}
	\begin{cases}
		\begin{cases}
			c_q = -\frac{a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			s_q = -\frac{a_{21}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{11} = \frac{a_{12}*a_{21} - a_{11}*a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{12} = \frac{-a_{11}*a_{21} - a_{12}*a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{22} = -\sqrt{a_{21}^2 + a_{22}^2} \\
		\end{cases} \\
		\begin{cases}
			c_q = \frac{a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			s_q = \frac{a_{21}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{11} = \frac{-a_{12}*a_{21} + a_{11}*a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{12} = \frac{a_{11}*a_{21} + a_{12}*a_{22}}{\sqrt{a_{21}^2 + a_{22}^2}} \\
			r_{22} = \sqrt{a_{21}^2 + a_{22}^2}
		\end{cases}
	\end{cases}
\end{equation}
$$

To me, the solutions seem to be equally stable numerically, so picking either works. The code is already stable numerically, so aside from scaling the input $a_{ij}$ to avoid overflow and underflow, no special treatment is necessary.