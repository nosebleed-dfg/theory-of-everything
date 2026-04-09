"""
Verification of all results. Cross-check exact traces against numerical eigenvalues.
Also verify the 120-cell using a completely independent method:
compute Tr(L+) = Tr((L + J/n)^{-1}) - 1/(eigenvalue that was 0, now 1)

For vertex-transitive graph: Tr(L+) = n * (L+)_{00}
And (L+)_{00} can be computed from the matrix directly.
"""

import numpy as np
from fractions import Fraction
from itertools import permutations, product
from collections import Counter
from sympy import isprime, factorint

phi = (1 + np.sqrt(5)) / 2
inv_phi = 1 / phi

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


print("=" * 70)
print("INDEPENDENT VERIFICATION USING MATRIX PSEUDOINVERSE")
print("=" * 70)
print()

# Method: compute L+ directly as a matrix and take its trace.
# For a graph Laplacian L of a connected graph:
# L+ = (L - J/n)^{-1} + J/n  ... no that's not right.
# L+ = (L + J/n)^{-1} - J/n  ...
# Actually: L+ is the pseudoinverse. For L = sum lambda_i v_i v_i^T,
# L+ = sum (1/lambda_i) v_i v_i^T (over nonzero lambda_i).
#
# Numerically: use numpy.linalg.pinv

# === 5-cell ===
n = 5
A = np.ones((n, n)) - np.eye(n)
L = np.diag(A.sum(axis=1)) - A
Lp = np.linalg.pinv(L)
tr = np.trace(Lp)
print(f"5-cell:     Tr(L+) = {tr:.15f}  (exact: {float(Fraction(4,5)):.15f})")

# === 8-cell ===
n = 16
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if bin(i ^ j).count('1') == 1:
            A[i][j] = 1
L = np.diag(A.sum(axis=1)) - A
Lp = np.linalg.pinv(L)
tr = np.trace(Lp)
print(f"8-cell:     Tr(L+) = {tr:.15f}  (exact: {float(Fraction(103,24)):.15f})")

# === 16-cell ===
verts_16cell = []
for i in range(4):
    v_pos = [0]*4; v_pos[i] = 1; verts_16cell.append(tuple(v_pos))
    v_neg = [0]*4; v_neg[i] = -1; verts_16cell.append(tuple(v_neg))
n = 8
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            ip = sum(a*b for a, b in zip(verts_16cell[i], verts_16cell[j]))
            if ip != -1:
                A[i][j] = 1
L = np.diag(A.sum(axis=1)) - A
Lp = np.linalg.pinv(L)
tr = np.trace(Lp)
print(f"16-cell:    Tr(L+) = {tr:.15f}  (exact: {float(Fraction(25,24)):.15f})")

# === 24-cell ===
verts_24cell = set()
for perm in set(permutations([1, 1, 0, 0])):
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            v = list(perm)
            nz = [i for i, x in enumerate(v) if x != 0]
            if len(nz) == 2:
                v[nz[0]] *= s1
                v[nz[1]] *= s2
                verts_24cell.add(tuple(v))
verts_24cell = sorted(verts_24cell)
n = len(verts_24cell)
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            dist_sq = sum((a-b)**2 for a, b in zip(verts_24cell[i], verts_24cell[j]))
            if abs(dist_sq - 2) < 1e-10:
                A[i][j] = 1
L = np.diag(A.sum(axis=1)) - A
Lp = np.linalg.pinv(L)
tr = np.trace(Lp)
exact_24 = Fraction(4,4) + Fraction(9,8) + Fraction(8,10) + Fraction(2,12)
print(f"24-cell:    Tr(L+) = {tr:.15f}  (exact: {float(exact_24):.15f})")

# === 600-cell ===
verts_600 = []
for i in range(4):
    for s in [1, -1]:
        v = [0.0]*4; v[i] = s; verts_600.append(tuple(v))
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
assert len(verts_600) == 120
n = 120
target_ip = phi / 2
A = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        ip = sum(a*b for a, b in zip(verts_600[i], verts_600[j]))
        if abs(ip - target_ip) < 1e-8:
            A[i][j] = 1; A[j][i] = 1
L = np.diag(A.sum(axis=1)) - A
Lp = np.linalg.pinv(L)
tr_600 = np.trace(Lp)
exact_600 = Fraction(16,9) + Fraction(25,12) + Fraction(18,7) + Fraction(16,15) + Fraction(72,36) + Fraction(180,80)
print(f"600-cell:   Tr(L+) = {tr_600:.15f}  (exact: {float(exact_600):.15f})")

# === 120-cell ===
print("\nBuilding 120-cell (this takes a moment)...")
adj_list = {i: set() for i in range(n)}
for i in range(n):
    for j in range(n):
        if A[i][j] == 1:
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
assert len(cells) == 600
cell_sets = [frozenset(c) for c in cells]
n120 = 600
A120 = np.zeros((n120, n120), dtype=np.int8)
for i in range(n120):
    for j in range(i+1, n120):
        if len(cell_sets[i] & cell_sets[j]) == 3:
            A120[i][j] = 1; A120[j][i] = 1

L120 = 4 * np.eye(n120) - A120.astype(float)
print("Computing pseudoinverse of 600x600 matrix...")
Lp120 = np.linalg.pinv(L120)
tr_120 = np.trace(Lp120)
print(f"120-cell:   Tr(L+) = {tr_120:.15f}  (claimed: {float(Fraction(1564270151,6163080)):.15f})")
print()

# Now let me try harder to identify the exact fraction.
# The denominator should be related to the LCM of eigenvalue denominators.
# From the eigenvalue analysis:
# Denominators involved: 1 (for integers), products of conjugate pairs
# Products: 1, 1, 3, 11, 15, 11, 23, 29
# Plus: 3, 4, 5, 6, 24, 28

# The LCM of all these: lcm(1, 3, 11, 15, 23, 29, 4, 5, 6, 24, 28)
from math import lcm
from functools import reduce
denoms = [1, 3, 11, 15, 23, 29, 4, 5, 6, 24, 28, 8]  # 8 from 11/8 (cubic1)
l = reduce(lcm, denoms)
print(f"LCM of denominators: {l}")
print(f"Expected denominator divides: {l}")
print(f"Actual denominator: {Fraction(1564270151, 6163080).denominator}")

# Let me check: is 6163080 divisible by the LCM?
print(f"LCM = {l}")
print(f"6163080 / {l} = {6163080 / l}")
print(f"6163080 = {factorint(6163080)}")
print(f"LCM factors = {factorint(l)}")

# Hmm, the denominator should be lcm of {1, 3, 5, 7, 8, 11, 15, 23, 24, 28, 29}
# Wait - the contributions have denominators:
# 1 (pair 1: 28)
# 1 (pair 2: 27)
# 3 (pair 3: 80/3)
# 11 (pair 4: 192/11)
# 5 (pair 5: 48/5)
# 11 (pair 6: 168/11)
# 23 (pair 7: 480/23)
# 29 (pair 8: 330/29)
# 3 (eigval 3: 40/3)
# 2 (eigval 4: 9/2)
# 5 (eigval 5: 8/5)
# 3 (eigval 6: 4/3)
# 8 (cubic1: 275/8)
# 7 (cubic2: 297/7)

contribs_exact = [
    Fraction(28),
    Fraction(27),
    Fraction(80, 3),
    Fraction(192, 11),
    Fraction(48, 5),
    Fraction(168, 11),
    Fraction(480, 23),
    Fraction(330, 29),
    Fraction(40, 3),
    Fraction(9, 2),
    Fraction(8, 5),
    Fraction(4, 3),
    Fraction(275, 8),
    Fraction(297, 7),
]

total_exact = sum(contribs_exact)
print(f"\nExact sum: {total_exact}")
print(f"Decimal: {float(total_exact):.15f}")
print(f"Numerator: {total_exact.numerator}")
print(f"Denominator: {total_exact.denominator}")
print(f"pinv trace: {tr_120:.15f}")
print(f"Match? {abs(float(total_exact) - tr_120) < 1e-8}")
print()

# Let me verify each individual contribution
print("Verifying contributions individually using eigenvalue data...")
eigs_120 = sorted(np.linalg.eigvalsh(L120))
groups = get_eigenvalue_multiplicities(eigs_120)

# Compute numerical trace from groups
tr_from_groups = sum(m/e for e, m in groups if abs(e) > 1e-6)
print(f"Trace from eigenvalue groups: {tr_from_groups:.15f}")
print(f"Trace from pinv: {tr_120:.15f}")
print(f"Match: {abs(tr_from_groups - tr_120) < 1e-8}")
print()

# Print ALL eigenvalue groups with more precision
print("All eigenvalue groups:")
for e, m in groups:
    if abs(e) > 1e-6:
        print(f"   {e:.15f} x {m}  ->  contrib = {m/e:.15f}")
print(f"   Sum of contributions = {tr_from_groups:.15f}")
print()

# Try to match numerical to exact
print(f"Exact total as fraction: {total_exact}")
print(f"As decimal: {float(total_exact):.15f}")
print(f"Numerical: {tr_from_groups:.15f}")
print(f"Difference: {float(total_exact) - tr_from_groups:.2e}")
print()

if abs(float(total_exact) - tr_from_groups) < 1e-8:
    print("*** EXACT RESULT CONFIRMED ***")
    print(f"Tr(L+) = {total_exact}")
    print(f"Numerator = {total_exact.numerator}")
    print(f"Denominator = {total_exact.denominator}")
    if isprime(total_exact.numerator):
        print(f"*** NUMERATOR IS PRIME! ***")
    else:
        print(f"Numerator factorization: {factorint(total_exact.numerator)}")
else:
    print("*** MISMATCH - need to recheck eigenvalue identifications ***")
    print(f"   Exact: {float(total_exact):.15f}")
    print(f"   Numerical: {tr_from_groups:.15f}")
    print(f"   Error: {abs(float(total_exact) - tr_from_groups):.2e}")
