"""
Verify 600-cell exact computation more carefully.
Also double-check: are there any multiplicity issues?
"""

import numpy as np
from fractions import Fraction
from sympy import isprime, factorint
from itertools import permutations, product

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

# Build 600-cell
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
eigs = sorted(np.linalg.eigvalsh(L))

# Verify via pinv
Lp = np.linalg.pinv(L)
tr_pinv = np.trace(Lp)
print(f"600-cell Tr(L+) via pinv: {tr_pinv:.15f}")
print()

# Eigenvalue groups
print("Eigenvalue groups:")
from collections import Counter

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

groups = get_eigenvalue_multiplicities(eigs)
for e, m in groups:
    print(f"   {e:.15f} x {m}")

print()
print("Verification of exact eigenvalue identifications:")
s5 = np.sqrt(5)

# 0 x 1
# 9 - 3*sqrt(5) x 4 = 2.291796...
# 10 - 2*sqrt(5) x 9 = 5.527864...
# 9 x 16
# 12 x 25
# 14 x 36
# 10 + 2*sqrt(5) x 9 = 14.472136...
# 15 x 16
# 9 + 3*sqrt(5) x 4 = 15.708204...

exact_eigs = [
    (0, 1),
    (9 - 3*s5, 4),
    (10 - 2*s5, 9),
    (9, 16),
    (12, 25),
    (14, 36),
    (10 + 2*s5, 9),
    (15, 16),
    (9 + 3*s5, 4),
]

print(f"{'Exact':<25} {'Numerical':<25} {'Mult':<6} {'Match?'}")
for (exact_val, exact_mult), (num_val, num_mult) in zip(exact_eigs, groups):
    match = abs(exact_val - num_val) < 1e-8 and exact_mult == num_mult
    print(f"   {exact_val:<23.15f} {num_val:<23.15f} {num_mult:<6} {'YES' if match else 'NO!'}")

print()

# Exact trace computation
# Pair 1: 9 +- 3*sqrt(5), mult 4
# Sum of reciprocals = 2*9 / (81 - 45) = 18/36 = 1/2
# Contribution: 4 * 1/2 = 2
pair1 = Fraction(4) * Fraction(18, 36)
print(f"Pair (9+-3sqrt5): 4 * 18/36 = {pair1}")

# Pair 2: 10 +- 2*sqrt(5), mult 9
# Sum of reciprocals = 2*10 / (100 - 20) = 20/80 = 1/4
# Contribution: 9 * 1/4 = 9/4
pair2 = Fraction(9) * Fraction(20, 80)
print(f"Pair (10+-2sqrt5): 9 * 20/80 = {pair2}")

# Rational eigenvalues
r1 = Fraction(16, 9)    # 9 x 16
r2 = Fraction(25, 12)   # 12 x 25
r3 = Fraction(36, 14)   # 14 x 36
r4 = Fraction(16, 15)   # 15 x 16

total = pair1 + pair2 + r1 + r2 + r3 + r4
print(f"\nTotal = {pair1} + {pair2} + {r1} + {r2} + {r3} + {r4}")
print(f"      = {total}")
print(f"      = {float(total):.15f}")
print(f"pinv:   {tr_pinv:.15f}")
print(f"Error: {abs(float(total) - tr_pinv):.2e}")
print()
print(f"Numerator: {total.numerator}")
print(f"Denominator: {total.denominator}")
print(f"Prime? {isprime(total.numerator)}")
print()

# Also verify r3 = 36/14 = 18/7
print(f"Note: 36/14 = {Fraction(36,14)}")

# Let me also look at whether 3701 has any special structure
print(f"\n3701 = {factorint(3701)}")
print(f"3701 is prime: {isprime(3701)}")

# Does 3701 relate to dodecahedral/icosahedral constants?
# phi = (1+sqrt(5))/2 ~ 1.618
# For the dodecahedron (3D): 137 = 5^3 + 12 = 125 + 12
# For the 600-cell (4D icosahedron): 3701 = ?
print(f"\n3701 / 137 = {3701/137:.6f}")  # 27.007... interesting
print(f"3701 - 137*27 = {3701 - 137*27}")  # = 3701 - 3699 = 2
print(f"3701 = 137 * 27 + 2")
print(f"     = 137 * 27 + 2")
print(f"     = 137 * 3^3 + 2")
print()

# The denominator 315 = 5 * 7 * 9 = 5 * 63
print(f"315 = {factorint(315)}")
# Compare with 3D dodecahedron denominator 15 = 3 * 5
# 315 / 15 = 21
print(f"315 / 15 = {315/15}")
