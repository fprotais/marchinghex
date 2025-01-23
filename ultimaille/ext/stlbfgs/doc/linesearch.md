# Short Moré-Thuente paper recap
## Definitions
It is assumed that the objective function \(\phi: \mathbb{R} \rightarrow \mathbb{R}\) defined on \([0, \infty]\) is smooth and continuously differentiable with \(\phi'(0) < 0\) we search for such a step \(\alpha > 0\) that that so-called *Strong Wolfe Conditions* hold (Equations 1.1 and 1.2 in the paper):
\[  \phi(\alpha) \le \phi(0) + \mu \phi'(0) \alpha \\
   |\phi'(\alpha)| \le \eta |\phi'(0)|             \]

Usually, \(\mu\) is a small value below \(1/2\) and \(\eta\) is a value close to \(1\). Note that \(\mu \le \eta\).

During the procedure we make use of an auxiliary function \(\psi(\alpha)\) defined as follows (just before Equation. 2.1 in the paper):
\[ \psi(\alpha) := \phi(\alpha) - \mu \phi'(0) \alpha \]

## Iterative search for step length
The algorithm can be summarized as follows (follows the *Search Algorithm* in Section 2 of the paper).

**N.B.** This algorithm uses function \(\psi\) in steps (2) and (3) until the following holds:
\[ \psi(\alpha_t) \le 0, \phi'(\alpha_t) > 0 \]
After this statement becomes true the algorithm above starts using function \(\phi\) in the steps (2) and (3) instead.

**Algorithm:** For a given step \(\alpha_t\) and an interval of values \([\alpha_l, \alpha_u]\):
1. Check if the Strong Wolfe Conditions hold for \(\alpha_t\). If they do - terminate the procedure with \(\alpha_t\) as the result.
2. Generate next step length from the interval \([\alpha_l, \alpha_u]\) and the current step using either function \(\psi\) or \(\phi\). This part is shown in the Section 4 "TRIAL VALUE SELECTION" in the paper.
3. Update the interval \([\alpha_l, \alpha_u]\) using either function \(\psi\) or \(\phi\) and the current step \(\alpha_t\). This part is covered in two algorithms: *Updating Algorithm* (right after theorem 2.1 in the paper) when the \(\psi\) function is still used and *Modified Updating Algorithm* (shown after theorem 3.2 in the paper) used after we switch to using function \(\phi\). These algorithms differ solely in the function used within them.

