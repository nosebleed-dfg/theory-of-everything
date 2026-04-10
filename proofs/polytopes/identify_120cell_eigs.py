"""
Identify the exact eigenvalues of the 120-cell Laplacian.
The 120-cell is the dual of the 600-cell in 4D.
Its symmetry group is H4 (order 14400).

The eigenvalues must be expressible in terms of:
- Integers
- sqrt(2), sqrt(5), phi = (1+sqrt(5))/2
- Possibly sqrt(2)*sqrt(5) = sqrt(10)

Strategy: use the minimal polynomials. Each eigenvalue satisfies a polynomial
with integer coefficients, and we can identify them by looking at conjugate pairs.
"""

import numpy as np
from fractions import Fraction
from sympy import (isprime, factorint, sqrt as ssqrt, Rational, nsimplify,
                   minimal_polynomial, Symbol, solve, factor, simplify)
from itertools import permutations, product
from collections import Counter

phi = (1 + np.sqrt(5)) / 2
inv_phi = 1 / phi
s5 = np.sqrt(5)
s2 = np.sqrt(2)

def even_permutations_4():
    all_perms = list(permutations([0, 1, 2, 3]))
    even = []
    for p in all_perms:
        inv = 0
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    inv += 1
        if inv % 2 == 0:
            even.append(p)
    return even

def deduplicate(verts, tol=1e-10):
    unique = []
    for v in verts:
        found = False
        for u in unique:
            if all(abs(a - b) < tol for a, b in zip(v, u)):
                found = True
                break
        if not found:
            unique.append(v)
    return unique

def get_eigenvalue_multiplicities(eigs, tol=1e-6):
    sorted_eigs = sorted(eigs)
    groups = []
    current_group = [sorted_eigs[0]]
    for e in sorted_eigs[1:]:
        if abs(e - current_group[-1]) < tol:
            current_group.append(e)
        else:
            groups.append(current_group)
            current_group = [e]
    groups.append(current_group)
    return [(np.mean(g), len(g)) for g in groups]


print("Rebuilding 600-cell and 120-cell...")

# Build 600-cell
verts_600 = []
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4
        v[i] = s
        verts_600.append(tuple(v))
for signs in product([0.5, -0.5], repeat=4):
    verts_600.append(tuple(signs))
even_perms = even_permutations_4()
for perm in even_perms:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                v = [0.0]*4
                v[perm[0]] = 0
                v[perm[1]] = s1 * 0.5
                v[perm[2]] = s2 * phi / 2
                v[perm[3]] = s3 * inv_phi / 2
                verts_600.append(tuple(v))

verts_600 = deduplicate(verts_600)
assert len(verts_600) == 120, f"Expected 120 vertices, got {len(verts_600)}"

# Build 600-cell adjacency
target_ip = phi / 2
n = 120
A600 = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(i+1, n):
        ip = sum(a*b for a, b in zip(verts_600[i], verts_600[j]))
        if abs(ip - target_ip) < 1e-8:
            A600[i][j] = 1
            A600[j][i] = 1

# Find K4 subgraphs (tetrahedral cells of 600-cell)
adj_list = {i: set() for i in range(n)}
for i in range(n):
    for j in range(n):
        if A600[i][j] == 1:
            adj_list[i].add(j)

cells = []
for i in range(n):
    ni = adj_list[i]
    for j in ni:
        if j > i:
            nij = ni & adj_list[j]
            for k in nij:
                if k > j:
                    nijk = nij & adj_list[k]
                    for l in nijk:
                        if l > k:
                            cells.append((i, j, k, l))

assert len(cells) == 600, f"Expected 600 cells, got {len(cells)}"

# Build 120-cell adjacency
cell_sets = [frozenset(c) for c in cells]
n120 = 600
A120 = np.zeros((n120, n120), dtype=np.int8)
for i in range(n120):
    for j in range(i+1, n120):
        if len(cell_sets[i] & cell_sets[j]) == 3:
            A120[i][j] = 1
            A120[j][i] = 1

deg120 = A120.sum(axis=1)
assert min(deg120) == max(deg120) == 4

print("Computing 120-cell eigenvalues with high precision...")
L120 = 4 * np.eye(n120) - A120.astype(np.float64)
eigs_120 = sorted(np.linalg.eigvalsh(L120))

groups = get_eigenvalue_multiplicities(eigs_120)
print(f"\nFound {len(groups)} distinct eigenvalues.")
print()

# The eigenvalues and multiplicities:
eig_mults = [(g[0], g[1]) for g in groups]

# Let me try to identify ALL eigenvalues systematically.
# The 120-cell graph has symmetry group H4 (order 14400).
# Its irreps have dimensions: 1, 4, 5, 6, 8, 9, 10, 16, 18, 20, 24, 25, 30, 36, 40, 48
# The multiplicities we see are: 1, 4, 9, 16, 25, 36, 24, 16, 24, 9, 36, 40, 25, 48, 18, 16, 30, 24, 8, 8, 24, 48, 25, 30, 36, 16, 4
# These are squares and other numbers that match H4 irrep dimensions.

# For eigenvalues involving sqrt(5): look for conjugate pairs under sqrt(5) -> -sqrt(5)
# For eigenvalues involving sqrt(2): look for conjugate pairs under sqrt(2) -> -sqrt(2)

# Group eigenvalues by their sum with potential conjugate:
print("Looking for conjugate pairs (sum should be rational):")
print()

identified = {}
remaining = list(range(len(eig_mults)))

# First, identify the obvious rational ones
for idx in remaining[:]:
    ev, mult = eig_mults[idx]
    rounded = round(ev)
    if abs(ev - rounded) < 1e-8:
        identified[idx] = (rounded, 1, mult, f"{rounded}")  # (value, denom, mult, name)
        remaining.remove(idx)
        print(f"   Rational: {ev:.10f} = {rounded} (mult {mult})")

print()

# Now look for sqrt(5) conjugate pairs
# If lambda = a + b*sqrt(5), then conjugate is a - b*sqrt(5)
# Sum = 2a, difference = 2b*sqrt(5), product = a^2 - 5b^2

print("Looking for sqrt(5) conjugate pairs:")
found_pairs = []
for i in remaining:
    for j in remaining:
        if i >= j:
            continue
        ev_i, mult_i = eig_mults[i]
        ev_j, mult_j = eig_mults[j]
        if mult_i != mult_j:
            continue
        s = ev_i + ev_j
        # Check if sum is rational (integer or simple fraction)
        s_round = round(s * 2) / 2  # try half-integers
        if abs(s - round(s)) < 1e-6:
            # Sum is integer. Check if product is rational.
            p = ev_i * ev_j
            p_round = round(p * 4) / 4  # try quarter-integers
            if abs(p - round(p * 4) / 4) < 1e-5:
                # This looks like a conjugate pair a +- b*sqrt(5)
                a = s / 2
                b_sq = (a**2 - p) / 5
                if b_sq > 0:
                    b = np.sqrt(b_sq)
                    # Verify
                    v1 = a - b * s5
                    v2 = a + b * s5
                    if abs(v1 - ev_i) < 1e-6 and abs(v2 - ev_j) < 1e-6:
                        found_pairs.append((i, j, a, b, mult_i))
                        print(f"   {ev_i:.10f} = {a:.6f} - {b:.6f}*sqrt(5) (mult {mult_i})")
                        print(f"   {ev_j:.10f} = {a:.6f} + {b:.6f}*sqrt(5) (mult {mult_i})")
                        print(f"      sum = {s:.6f}, product = {p:.6f}")
                        # Try to identify a and b as rationals
                        a_frac = Fraction(round(a * 12), 12)
                        b_frac_sq = Fraction(round(b_sq * 400), 400)
                        print(f"      a = {a:.10f}, a_frac = {a_frac}")
                        print(f"      b^2 = {b_sq:.10f}, b^2_frac = {b_frac_sq}")
                        print()

# For remaining eigenvalues, try sqrt(2) pairs
print()
print("Looking for sqrt(2) conjugate pairs:")
remaining_after_s5 = [i for i in remaining if not any(i in (p[0], p[1]) for p in found_pairs)]

for i in remaining_after_s5:
    for j in remaining_after_s5:
        if i >= j:
            continue
        ev_i, mult_i = eig_mults[i]
        ev_j, mult_j = eig_mults[j]
        if mult_i != mult_j:
            continue
        s = ev_i + ev_j
        if abs(s - round(s)) < 1e-6:
            p = ev_i * ev_j
            a = s / 2
            b_sq = (a**2 - p) / 2
            if b_sq > 0:
                b = np.sqrt(b_sq)
                v1 = a - b * s2
                v2 = a + b * s2
                if abs(v1 - ev_i) < 1e-6 and abs(v2 - ev_j) < 1e-6:
                    print(f"   {ev_i:.10f} = {a:.6f} - {b:.6f}*sqrt(2) (mult {mult_i})")
                    print(f"   {ev_j:.10f} = {a:.6f} + {b:.6f}*sqrt(2) (mult {mult_i})")
                    a_frac = Fraction(round(a * 12), 12)
                    b_frac_sq = Fraction(round(b_sq * 100), 100)
                    print(f"      a = {a:.10f}, b^2 = {b_sq:.10f}")
                    print()

# Now let me try a completely different approach:
# Use sympy nsimplify to identify each eigenvalue
print()
print("=" * 60)
print("SYMPY NSIMPLIFY IDENTIFICATION")
print("=" * 60)
print()

from sympy import nsimplify as ns, sqrt as sq, Rational as R, pi

for ev, mult in eig_mults:
    if abs(ev) < 1e-8:
        print(f"   0 (mult {mult})")
        continue
    # Try nsimplify with various irrationals
    try:
        exact = ns(ev, rational=False, tolerance=1e-8)
        print(f"   {ev:.10f} -> {exact} (mult {mult})")
    except:
        print(f"   {ev:.10f} -> FAILED (mult {mult})")

print()
print("=" * 60)
print("COMPUTING EXACT Tr(L^+) USING SYMPY NSIMPLIFY")
print("=" * 60)
print()

# Compute the trace numerically with high precision
tr_numerical = sum(mult / ev for ev, mult in eig_mults if abs(ev) > 1e-8)
print(f"Numerical Tr(L^+) = {tr_numerical:.15f}")

# Try to identify as a fraction
# The denominator should divide lcm of all eigenvalue denominators
# For the 120-cell with degree 4, eigenvalues are roots of char poly with integer coefficients
# The trace of pseudoinverse should be rational

# Use continued fraction expansion
from sympy import nsimplify
exact_tr = nsimplify(tr_numerical, rational=True, tolerance=1e-8)
print(f"nsimplify(Tr) = {exact_tr}")
print(f"As float: {float(exact_tr):.15f}")
print(f"Error: {abs(float(exact_tr) - tr_numerical):.2e}")

# Also try fraction identification more carefully
print()
print("Searching for exact fraction p/q...")
best_err = 1
best_frac = None
for q in range(1, 100000):
    p = round(tr_numerical * q)
    err = abs(p/q - tr_numerical)
    if err < best_err:
        best_err = err
        best_frac = Fraction(p, q)
        if err < 1e-10:
            break

print(f"Best fraction: {best_frac}")
print(f"Error: {best_err:.2e}")
print(f"Numerator = {best_frac.numerator}")
print(f"Denominator = {best_frac.denominator}")
print(f"Prime? {isprime(best_frac.numerator)}")
if not isprime(best_frac.numerator):
    print(f"Factorization: {factorint(best_frac.numerator)}")
