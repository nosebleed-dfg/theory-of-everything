"""
Identify the remaining unidentified eigenvalues of the 120-cell.
These are: mult 25, mult 36, mult 25, mult 36 (and their conjugates).

From the data:
  1.0745770822 x 25  -- unknown
  1.4818013007 x 36  -- unknown
  2.8218059200 x 36  -- unknown
  3.4480705698 x 25  -- unknown
  6.4773523480 x 25  -- unknown (conjugate of 3.4480705698?)
  6.6963927793 x 36  -- unknown (conjugate of ?)

Wait - let me recheck conjugate pairs. The sum of identified conjugate pairs:
  (0.1459, 6.8541) -> sum = 7
  (0.3820, 2.6180) -> sum = 3
  (0.6972, 4.3028) -> sum = 5
  (1.7639, 6.2361) -> sum = 8
  (2.2087, 6.7913) -> sum = 9
  (2.3820, 4.6180) -> sum = 7
  (3.5858, 6.4142) -> sum = 10
  (4.3820, 6.6180) -> sum = 11

Remaining unidentified:
  1.0746 x 25
  1.4818 x 36
  2.8218 x 36
  3.4481 x 25
  6.4774 x 25
  6.6964 x 36

Check potential pairs:
  1.0746 + 6.4774 = 7.5520 (not integer)
  1.0746 + 3.4481 = 4.5227 (not integer)
  1.4818 + 6.6964 = 8.1782 (not integer)
  1.4818 + 2.8218 = 4.3036 (not integer)

Hmm, none of these pair up with integer sums.
Maybe some involve BOTH sqrt(2) and sqrt(5)?
Or maybe sqrt(5) pairs but with fractional a?

Let me try: (1.0746 + 3.4481)/2 = 2.2613...
  (3.4481 - 1.0746)/2 = 1.1868... = 1.1868 / sqrt(5) = 0.5308... no

Let me try: what if they are a + b*sqrt(5) + c*sqrt(2)?
Or what if they satisfy a degree-4 minimal polynomial?
"""

import numpy as np
from fractions import Fraction
from sympy import isprime, factorint, nsimplify, sqrt, Rational, solve, Symbol, Poly

s2 = np.sqrt(2)
s5 = np.sqrt(5)
phi = (1 + s5) / 2

# The unidentified eigenvalues
unid = [
    (1.0745770822, 25),
    (1.4818013007, 36),
    (2.8218059200, 36),
    (3.4480705698, 25),
    (6.4773523480, 25),
    (6.6963927793, 36),
]

print("Trying to identify unidentified eigenvalues")
print()

# Try: a/2 + b/2 * sqrt(5) + c/2 * sqrt(2) + d/2 * sqrt(10)
# The minimal polynomial over Q should be degree 4 at most.

# Let me try to find the minimal polynomial numerically.
# If lambda satisfies a polynomial of degree 2: lambda^2 + p*lambda + q = 0
# then p = -(lambda + conjugate), q = lambda * conjugate

# For degree 4 (if over Q(sqrt(5), sqrt(2))):
# Look for the set of 4 conjugates: a +- b*s5 +- c*s2 +- ...

# Let me check if pairs sum to a half-integer:
print("Checking half-integer sums:")
for i in range(len(unid)):
    for j in range(i+1, len(unid)):
        s = unid[i][0] + unid[j][0]
        for denom in [1, 2, 3, 4, 5, 6]:
            if abs(s * denom - round(s * denom)) < 1e-5:
                print(f"   {unid[i][0]:.10f} + {unid[j][0]:.10f} = {s:.10f} ~= {round(s*denom)}/{denom}")

print()

# Let me try: maybe these come in groups of 3 (one eigenvalue space decomposed under H4)?
# Or groups related by different Galois actions.

# Try nsimplify with specific extensions
for ev, mult in unid:
    print(f"\n   lambda = {ev:.10f} (mult {mult})")

    # Try a + b*sqrt(5) where a,b are rationals with denominator up to 12
    best_err = 1e-4
    best_form = None
    for a_num in range(-20, 80):
        for a_den in [1, 2, 3, 4, 5, 6, 10, 12]:
            a = a_num / a_den
            for b_num in range(-20, 20):
                for b_den in [1, 2, 3, 4, 5, 6, 10, 12]:
                    b = b_num / b_den
                    # Try sqrt(5)
                    val = a + b * s5
                    if abs(val - ev) < best_err:
                        best_err = abs(val - ev)
                        best_form = f"{a_num}/{a_den} + {b_num}/{b_den}*sqrt(5)"
                    # Try sqrt(2)
                    val = a + b * s2
                    if abs(val - ev) < best_err:
                        best_err = abs(val - ev)
                        best_form = f"{a_num}/{a_den} + {b_num}/{b_den}*sqrt(2)"
                    # Try sqrt(10)
                    val = a + b * np.sqrt(10)
                    if abs(val - ev) < best_err:
                        best_err = abs(val - ev)
                        best_form = f"{a_num}/{a_den} + {b_num}/{b_den}*sqrt(10)"
                    # Try sqrt(3)
                    val = a + b * np.sqrt(3)
                    if abs(val - ev) < best_err:
                        best_err = abs(val - ev)
                        best_form = f"{a_num}/{a_den} + {b_num}/{b_den}*sqrt(3)"
                    # Try sqrt(6)
                    val = a + b * np.sqrt(6)
                    if abs(val - ev) < best_err:
                        best_err = abs(val - ev)
                        best_form = f"{a_num}/{a_den} + {b_num}/{b_den}*sqrt(6)"

    if best_form and best_err < 1e-6:
        print(f"      FOUND: {best_form} (err={best_err:.2e})")
    else:
        print(f"      Best so far: {best_form} (err={best_err:.2e})")

        # Try two-root forms: a + b*sqrt(5) + c*sqrt(2)
        best_err2 = 1e-4
        best_form2 = None
        for a_num in range(-20, 80):
            for a_den in [1, 2, 4]:
                a = a_num / a_den
                for b_num in range(-10, 10):
                    for b_den in [1, 2, 4]:
                        b = b_num / b_den
                        for c_num in range(-10, 10):
                            for c_den in [1, 2, 4]:
                                c = c_num / c_den
                                val = a + b * s5 + c * s2
                                if abs(val - ev) < best_err2:
                                    best_err2 = abs(val - ev)
                                    best_form2 = f"{a_num}/{a_den} + {b_num}/{b_den}*s5 + {c_num}/{c_den}*s2"

        if best_form2 and best_err2 < 1e-6:
            print(f"      TWO-ROOT: {best_form2} (err={best_err2:.2e})")
        else:
            print(f"      TWO-ROOT best: {best_form2} (err={best_err2:.2e})")

print()
print()

# Actually, let me try to compute the MINIMAL POLYNOMIAL of each eigenvalue
# If x is an eigenvalue of the Laplacian, it satisfies the characteristic polynomial
# which has integer coefficients. So the minimal polynomial divides it.
#
# For a numerical eigenvalue, find the minimal poly by checking:
# is it rational? Is it root of degree 2? degree 3? degree 4?

print("=" * 60)
print("MINIMAL POLYNOMIAL SEARCH")
print("=" * 60)

for ev, mult in unid:
    print(f"\n   lambda = {ev:.15f} (mult {mult})")

    # Check degree 2: x^2 + bx + c = 0  =>  c = -(x^2 + bx)
    # For integer b: c should be integer
    # b = -trace, c = product of roots (Vieta's)

    # Try: minimal poly with integer coefficients up to degree 4
    for deg in [2, 3, 4]:
        # Build Vandermonde-like system
        powers = [ev**k for k in range(deg + 1)]
        # We want: a0 + a1*x + ... + a_{deg-1}*x^{deg-1} + x^deg = 0
        # So: x^deg = -(a0 + a1*x + ... + a_{deg-1}*x^{deg-1})

        # For degree 2: x^2 + a1*x + a0 = 0
        # a1 = -(x + conjugate), a0 = x*conjugate
        # Since x is algebraic of degree <= 4, try integer and half-integer coefficients

        if deg == 2:
            for a1_num in range(-20, 20):
                for a1_den in [1]:
                    a1 = a1_num / a1_den
                    for a0_num in range(-50, 50):
                        for a0_den in [1]:
                            a0 = a0_num / a0_den
                            val = ev**2 + a1*ev + a0
                            if abs(val) < 1e-8:
                                print(f"      Deg 2: x^2 + {a1_num}*x + {a0_num} = 0 (residual={abs(val):.2e})")

        elif deg == 4:
            # x^4 + a*x^3 + b*x^2 + c*x + d = 0
            # Try small integer coefficients
            found = False
            for a3 in range(-20, 5):
                if found:
                    break
                for a2 in range(-5, 50):
                    if found:
                        break
                    for a1 in range(-50, 50):
                        if found:
                            break
                        for a0 in range(-50, 50):
                            val = ev**4 + a3*ev**3 + a2*ev**2 + a1*ev + a0
                            if abs(val) < 1e-6:
                                print(f"      Deg 4: x^4 + {a3}*x^3 + {a2}*x^2 + {a1}*x + {a0} = 0 (residual={abs(val):.2e})")
                                found = True
                                break

print()

# Now let me try a COMPLETELY different approach.
# The characteristic polynomial of the 120-cell Laplacian has integer coefficients.
# I can compute it using sympy with exact arithmetic.
# But 600x600 is too big for sympy.
#
# Instead: use the fact that the 120-cell is highly symmetric.
# Its eigenvalues can be computed from the spectrum of the H4 Coxeter group.
# The adjacency spectrum of the 120-cell is known in the literature.
# Let me look for it.

# Actually, let me just compute Tr(L+) using a different method:
# Tr(L+) = (1/n) * sum_{i<j} R_{ij} where R_{ij} is the effective resistance between i and j
# And for a vertex-transitive graph, R_{ij} depends only on the distance class.
# Tr(L+) = (1/n) * sum_i n_i * R_i where n_i is the number of pairs at distance i
# and R_i is the resistance at distance i.
#
# Or more directly: Tr(L+) = sum of 1/lambda for nonzero lambda (with multiplicity)
# I already have this numerically: 253.813053...
# I need to find the EXACT rational value.

# The trace of L+ for a vertex-transitive d-regular graph on n vertices:
# Tr(L+) is always rational (irrationals cancel in conjugate pairs).
# So let me find the exact fraction by computing with higher precision.

# Use mpmath for higher precision
from mpmath import mp, mpf, matrix, eighe, fsum

mp.dps = 50  # 50 decimal places

print("=" * 60)
print("HIGH PRECISION COMPUTATION")
print("=" * 60)

# Rebuild the Laplacian with mpmath
# Actually, since A120 is integer-valued, I can use numpy's integer matrix
# and compute eigenvalues with mpmath

# Convert A120 to mpmath matrix - but this is 600x600, might be slow
# Let me try a different approach: compute Tr(L+) = Tr((L + J/n)^{-1}) - 1/0
# Actually: Tr(L+) = Tr((L + (1/n)*J)^{-1}) - n/0... no.
#
# Tr(L+) = sum_{i=1}^{n-1} 1/lambda_i where lambda_i are nonzero eigenvalues
# For a connected graph, this equals: Tr(L^{-1} restricted to complement of null space)
#
# Alternative formula: Tr(L+) = sum_i sum_j (L+)_{ij} delta_{ij} = sum_i (L+)_{ii}
# And for vertex-transitive: all (L+)_{ii} are equal = Tr(L+)/n
# So Tr(L+) = n * (L+)_{11}
#
# And (L+)_{11} = (1/n) * sum_{j} R_{1j}/2... no.
#
# Actually: (L+)_{ii} = (1/n^2) * sum_{j} Omega_{ij} where Omega is the resistance matrix
# Hmm this is getting circular.
#
# Let me just use the Kirchhoff matrix tree theorem approach.
# Tr(L+) = (1/n) * sum_{i<j} R_{ij}  [no, that's different]
#
# For vertex-transitive: Tr(L+) = (n-1)/n * Kirchhoff_index / n  [not right either]
#
# The cleanest: Tr(L+) = sum of 1/lambda_i for nonzero eigenvalues.
# I need exact eigenvalues. Let me work with the characteristic polynomial.

# For the 120-cell, the adjacency eigenvalues are known from group theory.
# The graph is the Cayley graph of... no, it's not a Cayley graph.
# But its spectrum can be computed from the H4 representation theory.

# Let me try yet another approach: compute Tr(L+) = Tr((L + J/n)^{-1} - J/n)
# = Tr((L + J/n)^{-1}) - 1  [since Tr(J/n) = 1]
# And (L + J/n) is invertible since we've filled in the zero eigenvalue.

# With numpy this gives us a good numerical answer. The question is making it exact.
# Let me try with very aggressive fraction search.

# The eigenvalues involve at most sqrt(2), sqrt(5), sqrt(13), sqrt(21)
# Wait - sqrt(13) and sqrt(21)? Let me check.
# The nsimplify said: 5/2 - sqrt(13)/2 for the mult-16 eigenvalue.
# But the PAIR had sum 5 and product 3.
# x = 5/2 - sqrt(13)/2 means x^2 - 5x + 3 = 0
# Check: (5/2 - sqrt(13)/2)^2 - 5*(5/2-sqrt(13)/2) + 3
#       = 25/4 - 5*sqrt(13)/2 + 13/4 - 25/2 + 5*sqrt(13)/2 + 3
#       = 38/4 - 25/2 + 3 = 9.5 - 12.5 + 3 = 0 YES!

# And the other pair: 9/2 - sqrt(21)/2 has x^2 - 9x + 15 = 0
# Check: (9/2)^2 - 21/4 - 9*(9/2-sqrt(21)/2) + 15
#       = 81/4 - 21/4 - 81/2 + 9*sqrt(21)/2 + 15 = 60/4 - 81/2 + 15 = 15 - 40.5 + 15 = -10.5
# Hmm that's wrong. Let me recheck.
# x = 9/2 - sqrt(21)/2, x^2 = 81/4 - 9*sqrt(21)/2 + 21/4 = 102/4 - 9*sqrt(21)/2
# -9x = -81/2 + 9*sqrt(21)/2
# x^2 - 9x + 15 = 102/4 - 9*sqrt(21)/2 - 81/2 + 9*sqrt(21)/2 + 15
#                = 102/4 - 162/4 + 60/4 = 0/4 = 0 YES!

# OK so what about the unidentified ones? Let me check if they satisfy quartic polynomials.
# The minimal polynomial search found degree-4 results. Let me look harder.

print("\nDegree-4 minimal polynomial search (wider range):")
for ev, mult in unid:
    print(f"\n   lambda = {ev:.15f} (mult {mult})")
    found = False
    for a3 in range(-30, 5):
        if found:
            break
        for a2 in range(0, 100):
            if found:
                break
            for a1 in range(-200, 200):
                if found:
                    break
                for a0 in range(-100, 200):
                    val = ev**4 + a3*ev**3 + a2*ev**2 + a1*ev + a0
                    if abs(val) < 1e-5:
                        # Verify: also check that the other roots make sense
                        x = Symbol('x')
                        poly = x**4 + a3*x**3 + a2*x**2 + a1*x + a0
                        roots = solve(poly, x)
                        root_vals = [complex(r.evalf()) for r in roots]
                        real_roots = [r.real for r in root_vals if abs(r.imag) < 1e-8 and r.real > -0.1]
                        print(f"      x^4 + {a3}x^3 + {a2}x^2 + {a1}x + {a0} = 0 (res={abs(val):.2e})")
                        print(f"      Real roots: {[round(r, 6) for r in sorted(real_roots)]}")
                        found = True
                        break
