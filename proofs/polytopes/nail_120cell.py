"""
Nail down ALL 120-cell eigenvalues exactly and compute Tr(L+).

The key discovery from the minimal polynomial search:
- The mult-25 unidentified eigenvalues satisfy x^3 - 11x^2 + 33x - 24 = 0
- The mult-36 unidentified eigenvalues satisfy x^3 - 11x^2 + 33x - 28 = 0

For sum of reciprocals: if roots satisfy x^3 + ax^2 + bx + c = 0,
then sum(1/x_i) = -b/c (from Vieta's formulas applied to 1/x).
"""

import numpy as np
from sympy import (Symbol, solve, Rational, sqrt, simplify, Poly,
                   isprime, factorint, factor, nroots, N, re, im)
from fractions import Fraction

x = Symbol('x')

# Cubic 1: x^3 - 11x^2 + 33x - 24 = 0 (roots: eigenvalues with mult 25)
poly1 = x**3 - 11*x**2 + 33*x - 24
roots1_sym = solve(poly1, x)
roots1_num = nroots(poly1)
print("Cubic 1: x^3 - 11x^2 + 33x - 24 = 0")
print("Numerical roots:")
for r in sorted(roots1_num, key=lambda r: float(re(r))):
    print(f"   {float(re(N(r))):.10f} + {float(im(N(r))):.2e}i")

# Verify: By Vieta's, sum of roots = 11, sum of products of pairs = 33, product = 24
print(f"   Vieta check: sum={11}, sum_pairs={33}, product={24}")
print(f"   Sum of reciprocals = -b/c = -(33)/(-24) = 33/24 = {Fraction(33,24)}")
print()

# Actually: for x^3 + px^2 + qx + r = 0, sum(1/xi) = -q/r
# Our poly: x^3 - 11x^2 + 33x - 24 = 0 -> p=-11, q=33, r=-24
# sum(1/xi) = -q/r = -33/(-24) = 33/24 = 11/8
print(f"   sum(1/xi) = 33/24 = {Fraction(33,24)}")
# Verify numerically
sum_recip_1_num = sum(1/float(re(N(r))) for r in roots1_num)
print(f"   Numerical check: {sum_recip_1_num:.10f} vs {float(Fraction(33,24)):.10f}")
print()

# Cubic 2: x^3 - 11x^2 + 33x - 28 = 0 (roots: eigenvalues with mult 36)
poly2 = x**3 - 11*x**2 + 33*x - 28
roots2_num = nroots(poly2)
print("Cubic 2: x^3 - 11x^2 + 33x - 28 = 0")
print("Numerical roots:")
for r in sorted(roots2_num, key=lambda r: float(re(r))):
    print(f"   {float(re(N(r))):.10f} + {float(im(N(r))):.2e}i")
print(f"   Vieta: sum={11}, sum_pairs={33}, product={28}")
print(f"   sum(1/xi) = 33/28 = {Fraction(33,28)}")
sum_recip_2_num = sum(1/float(re(N(r))) for r in roots2_num)
print(f"   Numerical check: {sum_recip_2_num:.10f} vs {float(Fraction(33,28)):.10f}")
print()

# Verify the numerical eigenvalues against these cubics
target_25 = [1.0745770822, 3.4480705698, 6.4773523480]
target_36 = [1.4818013007, 2.8218059200, 6.6963927793]

print("Verification - plugging numerical eigenvalues into cubics:")
for t in target_25:
    val = t**3 - 11*t**2 + 33*t - 24
    print(f"   Cubic1({t:.10f}) = {val:.2e}")
for t in target_36:
    val = t**3 - 11*t**2 + 33*t - 28
    print(f"   Cubic2({t:.10f}) = {val:.2e}")
print()

# ============================================================
# EXACT Tr(L+) COMPUTATION
# ============================================================
print("=" * 70)
print("EXACT Tr(L+) FOR THE 120-CELL")
print("=" * 70)
print()

# All eigenvalue groups and their contributions to Tr(L+):
# For conjugate pairs a +- b*sqrt(d), the sum of reciprocals is:
#   1/(a-b*sqrt(d)) + 1/(a+b*sqrt(d)) = 2a/(a^2 - b^2*d)

contributions = []

# Pair 1: 7/2 +- 3*sqrt(5)/2, mult 4
# Product of conjugates = 49/4 - 45/4 = 1
# Sum of reciprocals = 2*(7/2)/1 = 7
c = 4 * Fraction(7, 1)
contributions.append(("7/2 +- 3sqrt(5)/2", "4 each", Fraction(7,1), c))

# Pair 2: 3/2 +- sqrt(5)/2, mult 9
# Product = 9/4 - 5/4 = 1
# Sum recip = 2*(3/2)/1 = 3
c = 9 * Fraction(3, 1)
contributions.append(("3/2 +- sqrt(5)/2", "9 each", Fraction(3,1), c))

# Pair 3: 5/2 +- sqrt(13)/2, mult 16
# Product = 25/4 - 13/4 = 3
# Sum recip = 2*(5/2)/3 = 5/3
c = 16 * Fraction(5, 3)
contributions.append(("5/2 +- sqrt(13)/2", "16 each", Fraction(5,3), c))

# Pair 4: 4 +- sqrt(5), mult 24
# Product = 16 - 5 = 11
# Sum recip = 2*4/11 = 8/11
c = 24 * Fraction(8, 11)
contributions.append(("4 +- sqrt(5)", "24 each", Fraction(8,11), c))

# Pair 5: 9/2 +- sqrt(21)/2, mult 16
# Product = 81/4 - 21/4 = 15
# Sum recip = 2*(9/2)/15 = 9/15 = 3/5
c = 16 * Fraction(3, 5)
contributions.append(("9/2 +- sqrt(21)/2", "16 each", Fraction(3,5), c))

# Pair 6: 7/2 +- sqrt(5)/2, mult 24
# Product = 49/4 - 5/4 = 11
# Sum recip = 2*(7/2)/11 = 7/11
c = 24 * Fraction(7, 11)
contributions.append(("7/2 +- sqrt(5)/2", "24 each", Fraction(7,11), c))

# Pair 7: 5 +- sqrt(2), mult 48
# Product = 25 - 2 = 23
# Sum recip = 2*5/23 = 10/23
c = 48 * Fraction(10, 23)
contributions.append(("5 +- sqrt(2)", "48 each", Fraction(10,23), c))

# Pair 8: 11/2 +- sqrt(5)/2, mult 30
# Product = 121/4 - 5/4 = 29
# Sum recip = 2*(11/2)/29 = 11/29
c = 30 * Fraction(11, 29)
contributions.append(("11/2 +- sqrt(5)/2", "30 each", Fraction(11,29), c))

# Rational eigenvalues
contributions.append(("3", "40", Fraction(1,3), 40 * Fraction(1,3)))
contributions.append(("4", "18", Fraction(1,4), 18 * Fraction(1,4)))
contributions.append(("5", "8", Fraction(1,5), 8 * Fraction(1,5)))
contributions.append(("6", "8", Fraction(1,6), 8 * Fraction(1,6)))

# Cubic roots (using Vieta's formulas)
contributions.append(("Cubic1 (x^3-11x^2+33x-24=0)", "25 each", Fraction(11,8), 25 * Fraction(11,8)))
contributions.append(("Cubic2 (x^3-11x^2+33x-28=0)", "36 each", Fraction(33,28), 36 * Fraction(33,28)))

print(f"{'Eigenvalue':<40} {'Mult':<10} {'sum 1/lam':<12} {'Contribution'}")
print("-" * 80)
total = Fraction(0)
for name, mult, sum_recip, contrib in contributions:
    print(f"   {name:<38} {str(mult):<10} {str(sum_recip):<12} {contrib}")
    total += contrib

print("-" * 80)
print(f"   {'TOTAL':<38} {'':<10} {'':<12} {total}")
print()
print(f"Tr(L^+) = {total}")
print(f"         = {float(total):.15f}")
print()

# Numerical check from eigenvalue computation
print(f"Numerical value from eigenvalues: 253.813053051396")
print(f"Error: {abs(float(total) - 253.813053051396):.2e}")
print()

print(f"Numerator = {total.numerator}")
print(f"Denominator = {total.denominator}")
print()

if isprime(total.numerator):
    print(f"*** NUMERATOR {total.numerator} IS PRIME! ***")
else:
    print(f"Numerator factorization: {factorint(total.numerator)}")
    print(f"*** NUMERATOR IS NOT PRIME ***")
print()

# Multiplicity check
mults_check = (
    1 +      # zero eigenvalue
    4*2 +    # pair 1
    9*2 +    # pair 2
    16*2 +   # pair 3
    24*2 +   # pair 4
    16*2 +   # pair 5
    24*2 +   # pair 6
    48*2 +   # pair 7
    30*2 +   # pair 8
    40 +     # eigenvalue 3
    18 +     # eigenvalue 4
    8 +      # eigenvalue 5
    8 +      # eigenvalue 6
    25*3 +   # cubic 1 (3 roots)
    36*3     # cubic 2 (3 roots)
)
print(f"Total multiplicity: {mults_check} (should be 600)")
print()

# ============================================================
# COMPLETE SUMMARY
# ============================================================
print("=" * 70)
print("COMPLETE RESULTS: ALL 6 REGULAR 4-POLYTOPES")
print("=" * 70)
print()

results = []

# 5-cell
tr = Fraction(4, 5)
results.append(("5-cell (simplex)", 5, 4, tr, tr.numerator, bool(isprime(tr.numerator))))

# 8-cell
tr = Fraction(4,2) + Fraction(6,4) + Fraction(4,6) + Fraction(1,8)
results.append(("8-cell (tesseract)", 16, 4, tr, tr.numerator, bool(isprime(tr.numerator))))

# 16-cell
tr = Fraction(4,6) + Fraction(3,8)
results.append(("16-cell", 8, 6, tr, tr.numerator, bool(isprime(tr.numerator))))

# 24-cell
tr = Fraction(4,4) + Fraction(9,8) + Fraction(8,10) + Fraction(2,12)
results.append(("24-cell", 24, 8, tr, tr.numerator, bool(isprime(tr.numerator))))

# 600-cell
tr = Fraction(16,9) + Fraction(25,12) + Fraction(18,7) + Fraction(16,15) + Fraction(72,36) + Fraction(180,80)
results.append(("600-cell", 120, 12, tr, tr.numerator, bool(isprime(tr.numerator))))

# 120-cell
results.append(("120-cell", 600, 4, total, total.numerator, bool(isprime(total.numerator))))

print(f"{'Polytope':<22} {'V':<6} {'d':<4} {'Tr(L+)':<25} {'Num':<15} {'Prime?'}")
print("-" * 85)
for name, V, d, tr, num, prime in results:
    prime_str = "YES!" if prime else "no"
    print(f"   {name:<20} {V:<6} {d:<4} {str(tr):<25} {num:<15} {prime_str}")

print()
print("For comparison, the 3D Platonic solids:")
print()
results_3d = [
    ("Tetrahedron", 4, 3, Fraction(3,4)),
    ("Cube", 8, 3, Fraction(29,12)),
    ("Octahedron", 6, 4, Fraction(13,12)),
    ("Icosahedron", 12, 5, Fraction(7,3)),
    ("Dodecahedron", 20, 3, Fraction(137,15)),
]
for name, V, d, tr in results_3d:
    prime = bool(isprime(tr.numerator))
    prime_str = "YES!" if prime else "no"
    print(f"   {name:<20} {V:<6} {d:<4} {str(tr):<25} {tr.numerator:<15} {prime_str}")

# Count total primes
total_3d_primes = sum(1 for _,_,_,tr in results_3d if isprime(tr.numerator))
total_4d_primes = sum(1 for _,_,_,_,_,p in results if p)
print()
print(f"3D: {total_3d_primes}/5 have prime numerators (ALL FIVE)")
print(f"4D: {total_4d_primes}/6 have prime numerators")

if not isprime(total.numerator):
    print(f"\n120-cell numerator factorization: {total.numerator} = {factorint(total.numerator)}")
