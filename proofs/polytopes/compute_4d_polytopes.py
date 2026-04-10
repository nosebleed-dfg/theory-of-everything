"""
Compute Tr(L^+) for all 6 regular 4D polytopes.
Uses numpy for eigenvalue computation.

The 6 regular 4-polytopes:
1. 5-cell (simplex): 5 vertices, each connected to all others (K5), degree 4
2. 8-cell (tesseract): 16 vertices, degree 4
3. 16-cell: 8 vertices, degree 6
4. 24-cell: 24 vertices, degree 8
5. 600-cell: 120 vertices, degree 12
6. 120-cell: 600 vertices, degree 4
"""

import numpy as np
from fractions import Fraction
from itertools import permutations, product
from collections import Counter
from sympy import isprime, factorint

phi = (1 + np.sqrt(5)) / 2
inv_phi = 1 / phi


def even_permutations_4():
    """Return all 12 even permutations of (0,1,2,3)"""
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
    """Group eigenvalues by rounding and count multiplicities."""
    # Use a smarter clustering approach
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


def compute_trace_numerical(eig_groups):
    """Compute Tr(L+) from eigenvalue groups, skipping zero eigenvalue."""
    total = 0.0
    for eigval, mult in eig_groups:
        if abs(eigval) > 1e-6:
            total += mult / eigval
    return total


# ============================================================
print("=" * 60)
print("1. 5-CELL (4-simplex, K5)")
print("   5 vertices, degree 4")
print()

n = 5
A = np.ones((n, n), dtype=int) - np.eye(n, dtype=int)
L = np.diag(A.sum(axis=1)) - A
eigs = sorted(np.linalg.eigvalsh(L.astype(float)))
groups = get_eigenvalue_multiplicities(eigs)
print(f"   Eigenvalue groups: {[(round(e,6), m) for e, m in groups]}")

# Exact: eigenvalue 0 (mult 1), eigenvalue 5 (mult 4)
tr = Fraction(4, 5)
print(f"   Tr(L^+) = {tr} = {float(tr):.10f}")
print(f"   Numerical check: {compute_trace_numerical(groups):.10f}")
print(f"   Numerator = {tr.numerator}, Denominator = {tr.denominator}")
print(f"   Numerator prime? {isprime(tr.numerator)}")
if not isprime(tr.numerator):
    print(f"   Factorization: {factorint(tr.numerator)}")
print()

# ============================================================
print("=" * 60)
print("2. 8-CELL (Tesseract, Q4)")
print("   16 vertices, degree 4")
print()

n = 16
A = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if bin(i ^ j).count('1') == 1:
            A[i][j] = 1
L = np.diag(A.sum(axis=1)) - A
eigs = sorted(np.linalg.eigvalsh(L.astype(float)))
groups = get_eigenvalue_multiplicities(eigs)
print(f"   Eigenvalue groups: {[(round(e,6), m) for e, m in groups]}")

# Exact: Laplacian eigenvalues of Q4: 0(1), 2(4), 4(6), 6(4), 8(1)
tr = Fraction(4, 2) + Fraction(6, 4) + Fraction(4, 6) + Fraction(1, 8)
print(f"   Exact: 4/2 + 6/4 + 4/6 + 1/8 = {tr}")
print(f"   Numerical check: {compute_trace_numerical(groups):.10f}")
print(f"   Numerator = {tr.numerator}, Denominator = {tr.denominator}")
print(f"   Numerator prime? {isprime(tr.numerator)}")
if not isprime(tr.numerator):
    print(f"   Factorization: {factorint(tr.numerator)}")
print()

# ============================================================
print("=" * 60)
print("3. 16-CELL")
print("   8 vertices, degree 6")
print()

# Vertices: +-e_i in R^4
verts_16cell = []
for i in range(4):
    v_pos = [0, 0, 0, 0]
    v_pos[i] = 1
    verts_16cell.append(tuple(v_pos))
    v_neg = [0, 0, 0, 0]
    v_neg[i] = -1
    verts_16cell.append(tuple(v_neg))

n = 8
A = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j:
            ip = sum(a * b for a, b in zip(verts_16cell[i], verts_16cell[j]))
            # Adjacent iff NOT antipodal (ip != -1)
            if ip != -1:
                A[i][j] = 1

L = np.diag(A.sum(axis=1)) - A
deg = A.sum(axis=1)[0]
print(f"   Degree: {deg}")
eigs = sorted(np.linalg.eigvalsh(L.astype(float)))
groups = get_eigenvalue_multiplicities(eigs)
print(f"   Eigenvalue groups: {[(round(e,6), m) for e, m in groups]}")

# Compute exact trace from groups
tr_num = compute_trace_numerical(groups)
print(f"   Numerical Tr(L^+) = {tr_num:.10f}")

# The 16-cell graph is the cocktail party graph CP(4) = K_{2,2,2,2}
# Complete multipartite graph K(2,2,2,2) with 4 parts of size 2
# Adjacency eigenvalues: for K(n1,...,nk) with all ni=m, k parts:
# Eigenvalue km=N (once, but actually eigenvalue of adj = (k-1)*m for the uniform vector)
# Wait, let me just read off from numerical.
# Eigenvalues should be: 0, 6, 6, 6, 8, 8, 8, 8? Let me check.
# CP(4) = complement of 4*K2
# Adjacency eigenvalues of 4*K2: 1 (4 times), -1 (4 times)
# Complement: A_bar = J - I - A
# Eigenvalue for all-ones vector of A: sum of row = 1 for each K2, so 1
# A_bar all-ones eigenvalue: N-1-1 = 6
# For other eigenvectors of A with eigenvalue mu: A_bar eigenvalue = -1-mu
# Eigenvalue 1 (but NOT the all-ones eigenvector): A_bar gives -1-1 = -2
# Eigenvalue -1: A_bar gives -1-(-1) = 0
# So adj eigenvalues of CP(4): 6 (once), -2 (3 times), 0 (4 times)
# Wait that's only 8, and multiplicities: 1+3+4 = 8.
# But we need: the eigenvectors of 4*K2.
# 4*K2 has 4 components each K2. For each K2: eigvals 1, -1.
# So the 8 eigenvalues of adj(4*K2): 1,1,1,1,-1,-1,-1,-1
# The all-ones eigenvector of 4*K2 has eigenvalue 1.
# The other three eigenvalue-1 eigenvectors are "inter-component" combinations.
# The four eigenvalue -1 eigenvectors are "within-component" differences.
# Complement eigenvalues:
# - all-ones: N-1-mu = 7-1 = 6 (once)
# - other eigenvalue 1: -1-1 = -2 (three times)
# - eigenvalue -1: -1-(-1) = 0 (four times)
# Adjacency eigenvalues of CP(4): 6(1), -2(3), 0(4)
# Laplacian eigenvalues: degree - adj = 6 - adj
# 6 - 6 = 0 (once)
# 6 - (-2) = 8 (three times)
# 6 - 0 = 6 (four times)
# Laplacian: 0(1), 6(4), 8(3)

tr_exact = Fraction(4, 6) + Fraction(3, 8)
print(f"   Exact: 4/6 + 3/8 = {Fraction(4,6)} + {Fraction(3,8)} = {tr_exact}")
print(f"   Numerator = {tr_exact.numerator}, Denominator = {tr_exact.denominator}")
print(f"   Numerator prime? {isprime(tr_exact.numerator)}")
if not isprime(tr_exact.numerator):
    print(f"   Factorization: {factorint(tr_exact.numerator)}")
print()

# ============================================================
print("=" * 60)
print("4. 24-CELL")
print("   24 vertices, degree 8")
print()

# Generate 24-cell vertices: all permutations of (+-1, +-1, 0, 0)
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
print(f"   Number of vertices: {len(verts_24cell)}")

n = len(verts_24cell)
A = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j:
            dist_sq = sum((a - b) ** 2 for a, b in zip(verts_24cell[i], verts_24cell[j]))
            if abs(dist_sq - 2) < 1e-10:
                A[i][j] = 1

deg = A.sum(axis=1)
print(f"   Degree: {deg[0]} (all same: {len(set(deg)) == 1})")

L = np.diag(deg) - A
eigs = sorted(np.linalg.eigvalsh(L.astype(float)))
groups = get_eigenvalue_multiplicities(eigs)
print(f"   Eigenvalue groups: {[(round(e,6), m) for e, m in groups]}")

tr_num = compute_trace_numerical(groups)
print(f"   Numerical Tr(L^+) = {tr_num:.10f}")

# Try to identify exact eigenvalues and compute exact trace
print("   Attempting exact identification...")
for eigval, mult in groups:
    print(f"      {eigval:.10f} x {mult}")
print()

# ============================================================
print("=" * 60)
print("5. 600-CELL")
print("   120 vertices, degree 12")
print()

# Build 600-cell vertices
verts_600 = []

# Type 1: 8 vertices - permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1, -1]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s * 1.0
        verts_600.append(tuple(v))

# Type 2: 16 vertices - (+-1/2, +-1/2, +-1/2, +-1/2)
for signs in product([0.5, -0.5], repeat=4):
    verts_600.append(tuple(signs))

# Type 3: 96 vertices - all EVEN permutations of (0, +-1/2, +-phi/2, +-1/(2*phi))
even_perms = even_permutations_4()

for perm in even_perms:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                v = [0.0, 0.0, 0.0, 0.0]
                v[perm[0]] = 0
                v[perm[1]] = s1 * 0.5
                v[perm[2]] = s2 * phi / 2
                v[perm[3]] = s3 * inv_phi / 2
                verts_600.append(tuple(v))

verts_600 = deduplicate(verts_600)
print(f"   Vertices: {len(verts_600)}")

# Check norms
norms = [np.sqrt(sum(x ** 2 for x in v)) for v in verts_600]
print(f"   Norm range: [{min(norms):.10f}, {max(norms):.10f}]")

# Build adjacency: inner product = phi/2 for adjacent vertices
target_ip = phi / 2
n = len(verts_600)
A = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(i + 1, n):
        ip = sum(a * b for a, b in zip(verts_600[i], verts_600[j]))
        if abs(ip - target_ip) < 1e-8:
            A[i][j] = 1
            A[j][i] = 1

deg = A.sum(axis=1)
print(f"   Degree range: [{min(deg)}, {max(deg)}]")

if min(deg) == max(deg) == 12:
    print("   Degree check PASSED!")
    L = 12 * np.eye(n) - A.astype(float)
    print("   Computing eigenvalues of 120x120 matrix...")
    eigs = sorted(np.linalg.eigvalsh(L))
    groups = get_eigenvalue_multiplicities(eigs)
    print(f"   Eigenvalue groups:")
    for eigval, mult in groups:
        print(f"      {eigval:.10f} x {mult}")

    tr_num = compute_trace_numerical(groups)
    print(f"   Numerical Tr(L^+) = {tr_num:.10f}")
else:
    print(f"   DEGREE CHECK FAILED!")
    print(f"   Degree distribution: {Counter(deg.astype(int))}")
    # Debug: show all distinct inner products
    all_ips = set()
    for i in range(min(n, 50)):
        for j in range(i + 1, min(n, 50)):
            ip = sum(a * b for a, b in zip(verts_600[i], verts_600[j]))
            all_ips.add(round(ip, 6))
    print(f"   Distinct inner products (sample): {sorted(all_ips)}")

print()

# ============================================================
print("=" * 60)
print("6. 120-CELL")
print("   600 vertices, degree 4")
print()

# Construct from 600-cell dual: find tetrahedral cells of 600-cell
if min(deg) == max(deg) == 12:
    print("   Building adjacency list for 600-cell...")
    adj_list = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj_list[i].add(j)

    print("   Finding tetrahedral cells (K4 subgraphs)...")
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

    print(f"   Found {len(cells)} tetrahedral cells (expected 600)")

    if len(cells) == 600:
        # Cell centers = 120-cell vertices
        cell_centers = []
        for cell in cells:
            center = [0.0, 0.0, 0.0, 0.0]
            for v_idx in cell:
                for d in range(4):
                    center[d] += verts_600[v_idx][d]
            center = tuple(c / 4 for c in center)
            cell_centers.append(center)

        # Normalize
        norms_cc = [np.sqrt(sum(x ** 2 for x in c)) for c in cell_centers]
        print(f"   Center norms: [{min(norms_cc):.8f}, {max(norms_cc):.8f}]")

        # Two cells adjacent iff they share a face (3 vertices)
        print("   Building 120-cell adjacency (600x600)...")
        cell_sets = [frozenset(c) for c in cells]

        n120 = 600
        A120 = np.zeros((n120, n120), dtype=np.int8)
        for i in range(n120):
            for j in range(i + 1, n120):
                if len(cell_sets[i] & cell_sets[j]) == 3:
                    A120[i][j] = 1
                    A120[j][i] = 1

        deg120 = A120.sum(axis=1)
        print(f"   120-cell degree range: [{min(deg120)}, {max(deg120)}]")

        if min(deg120) == max(deg120) == 4:
            print("   Degree 4 check PASSED!")
            print("   Computing eigenvalues of 600x600 matrix...")
            L120 = 4 * np.eye(n120) - A120.astype(float)
            eigs_120 = sorted(np.linalg.eigvalsh(L120))
            groups_120 = get_eigenvalue_multiplicities(eigs_120)
            print(f"   Eigenvalue groups:")
            for eigval, mult in groups_120:
                print(f"      {eigval:.10f} x {mult}")

            tr_num_120 = compute_trace_numerical(groups_120)
            print(f"   Numerical Tr(L^+) = {tr_num_120:.10f}")
        else:
            print(f"   DEGREE CHECK FAILED!")
            print(f"   Degree dist: {Counter(deg120.astype(int))}")
    else:
        print(f"   Wrong cell count. Cannot build 120-cell.")
else:
    print("   Skipped (600-cell failed)")

print()
print("=" * 60)
print("COMPUTATION COMPLETE")
print("=" * 60)
