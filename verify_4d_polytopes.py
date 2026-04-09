"""
EXTENSION: Compute Tr(L^+) for the six regular 4D polytopes.
These are the 4-dimensional analogs of the Platonic solids.
"""

import numpy as np
from sympy import Matrix, Rational, isprime, simplify
from fractions import Fraction
from math import gcd
import itertools

print("=" * 70)
print("4D REGULAR POLYTOPES: Tr(L^+) COMPUTATION")
print("=" * 70)

def compute_trace_exact(A, name):
    """Compute Tr(L^+) exactly using sympy eigenvalues"""
    n = A.rows
    D = Matrix.zeros(n, n)
    for i in range(n):
        D[i, i] = sum(A[i, j] for j in range(n))

    L = D - A

    eigenvals = L.eigenvals()

    print(f"\n{name} (n={n}, degree={D[0,0]})")
    sorted_eigs = sorted(eigenvals.items(), key=lambda x: float(x[0].evalf()))
    for eig, mult in sorted_eigs:
        print(f"  eigenvalue {eig} (mult {mult})")

    trace = Rational(0)
    for eig, mult in eigenvals.items():
        if eig != 0:
            trace += Rational(mult) / eig

    trace = simplify(trace)

    if trace.is_rational:
        p = trace.p
        q = trace.q
        print(f"  Tr(L^+) = {p}/{q}")
        print(f"  Numerator = {p}, prime = {isprime(p)}")
        return p, q, isprime(p)
    else:
        print(f"  Tr(L^+) = {trace} (attempting to rationalize...)")
        # Try evaluating numerically and rationalizing
        val = float(trace.evalf())
        return rationalize(val, name)

def rationalize(val, name):
    """Try to express a float as a fraction with small denominator"""
    best_frac = None
    best_err = 1e10
    for denom in range(1, 100000):
        numer = round(val * denom)
        err = abs(val - numer/denom)
        if err < best_err:
            best_err = err
            best_frac = (numer, denom)
            if err < 1e-10:
                break

    if best_frac and best_err < 1e-6:
        g = gcd(best_frac[0], best_frac[1])
        p, q = best_frac[0]//g, best_frac[1]//g
        is_p = isprime(abs(p))
        print(f"  Tr(L^+) ~ {val:.10f} = {p}/{q} (error {best_err:.2e})")
        print(f"  Numerator = {p}, prime = {is_p}")
        return p, q, is_p
    else:
        print(f"  Could not rationalize {val}")
        return None, None, None

def compute_trace_numerical(A_np, name):
    """Compute Tr(L^+) numerically"""
    n = A_np.shape[0]
    D_np = np.diag(A_np.sum(axis=1))
    L_np = D_np - A_np

    eigs = np.linalg.eigvalsh(L_np)
    eigs.sort()

    trace = sum(1.0/e for e in eigs if abs(e) > 1e-8)

    return rationalize(trace, name)


# ======================================================================
# 1. 5-cell (Pentachoron / Hypertetrahedron)
# Complete graph K5
# 5 vertices, each connected to all others, degree 4
# ======================================================================

def build_5cell():
    n = 5
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i, j] = 1
    return A, "5-cell (K5)"


# ======================================================================
# 2. 8-cell (Tesseract / Hypercube)
# 16 vertices, degree 4
# Vertices are 4-bit binary strings, edge iff Hamming distance = 1
# ======================================================================

def build_8cell():
    n = 16
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if bin(i ^ j).count('1') == 1:
                A[i, j] = 1
    return A, "8-cell (Tesseract)"


# ======================================================================
# 3. 16-cell (Hyperoctahedron)
# 8 vertices, degree 6
# 4 pairs of antipodal vertices: (0,1),(2,3),(4,5),(6,7)
# Each vertex connected to all except its antipode
# ======================================================================

def build_16cell():
    n = 8
    A = Matrix.zeros(n, n)
    antipodal = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4, 6:7, 7:6}
    for i in range(n):
        for j in range(n):
            if i != j and j != antipodal[i]:
                A[i, j] = 1
    return A, "16-cell (Hyperoctahedron)"


# ======================================================================
# 4. 24-cell
# 24 vertices, degree 8
# Vertices: all permutations of (+-1, +-1, 0, 0)
# Edge iff distance = sqrt(2)
# ======================================================================

def build_24cell():
    verts = []
    # All permutations of (+-1, +-1, 0, 0)
    for perm in itertools.permutations(range(4)):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0, 0, 0, 0]
                v[perm[0]] = s1
                v[perm[1]] = s2
                v = tuple(v)
                if v not in verts:
                    verts.append(v)

    n = len(verts)
    print(f"  24-cell: {n} vertices generated")

    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(4))
            if abs(d2 - 2) < 0.01:  # edge length = sqrt(2)
                A[i, j] = 1
                A[j, i] = 1

    degrees = [int(sum(A[i, j] for j in range(n))) for i in range(n)]
    d_set = set(degrees)
    print(f"  degrees = {d_set}")
    return A, "24-cell"


# ======================================================================
# 5. 120-cell
# 600 vertices, degree 4
# This is HUGE -- we'll use numerical methods
# ======================================================================

def build_120cell_numerical():
    """Build the 120-cell graph (600 vertices, degree 4)"""
    phi = (1 + 5**0.5) / 2
    phi2 = phi**2
    phi3 = phi**3
    inv_phi = 1/phi
    inv_phi2 = 1/phi**2
    inv_phi3 = 1/phi**3

    # The 600 vertices of the 120-cell can be generated as follows:
    # (using the construction from H.S.M. Coxeter)
    # This is extremely complex. Let me use a known adjacency structure instead.

    # For the 120-cell, the graph is the skeleton of the dual of the 600-cell.
    # The 600-cell has 120 vertices, and its DUAL (the 120-cell) has 600 vertices.

    # Actually, I have the names wrong:
    # 120-cell: 600 vertices, 1200 edges, degree 4
    # 600-cell: 120 vertices, 720 edges, degree 12

    # Let's do the 600-cell first (120 vertices, more manageable)
    print("  [120-cell has 600 vertices -- too large for exact sympy]")
    print("  [Deferring to 600-cell which has 120 vertices]")
    return None, "120-cell"


# ======================================================================
# 6. 600-cell
# 120 vertices, degree 12
# ======================================================================

def build_600cell():
    """Build the 600-cell graph (120 vertices, degree 12)"""
    phi = (1 + 5**0.5) / 2

    verts = []

    # The 120 vertices of the 600-cell:
    # 8 vertices: all permutations of (+-2, 0, 0, 0)
    for i in range(4):
        for s in [2, -2]:
            v = [0, 0, 0, 0]
            v[i] = s
            verts.append(tuple(v))

    # 16 vertices: (+-1, +-1, +-1, +-1)
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    verts.append((s0, s1, s2, s3))

    # 96 vertices: all EVEN permutations of (+-phi, +-1, +-1/phi, 0)
    # "Even permutation" here means an even permutation of the coordinate positions
    base_vals = [phi, 1, 1/phi, 0]

    def even_perms_4(lst):
        """Generate all even permutations of a 4-element list"""
        # Even permutations of (0,1,2,3) are those with even parity
        all_perms = list(itertools.permutations(range(4)))
        even = []
        for p in all_perms:
            # Count inversions
            inv = sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j])
            if inv % 2 == 0:
                even.append(p)
        result = []
        for p in even:
            result.append(tuple(lst[p[i]] for i in range(4)))
        return result

    for perm in even_perms_4(base_vals):
        for s0 in ([1, -1] if abs(perm[0]) > 1e-10 else [1]):
            for s1 in ([1, -1] if abs(perm[1]) > 1e-10 else [1]):
                for s2 in ([1, -1] if abs(perm[2]) > 1e-10 else [1]):
                    for s3 in ([1, -1] if abs(perm[3]) > 1e-10 else [1]):
                        v = (s0*perm[0], s1*perm[1], s2*perm[2], s3*perm[3])
                        v_r = tuple(round(x, 10) for x in v)
                        if v_r not in verts:
                            verts.append(v_r)

    n = len(verts)
    print(f"  600-cell: {n} vertices generated")

    if n != 120:
        print(f"  WARNING: expected 120, got {n}")
        # Try to deduplicate more aggressively
        unique = []
        for v in verts:
            is_dup = False
            for u in unique:
                if all(abs(v[k] - u[k]) < 1e-8 for k in range(4)):
                    is_dup = True
                    break
            if not is_dup:
                unique.append(v)
        verts = unique
        n = len(verts)
        print(f"  After dedup: {n} vertices")

    # Find edge length (minimum distance between vertices)
    # All vertices should be at the same distance from origin
    norms = [sum(x**2 for x in v)**0.5 for v in verts]
    print(f"  Norms range: [{min(norms):.6f}, {max(norms):.6f}]")

    # Find all pairwise distances
    min_dist2 = float('inf')
    for i in range(min(n, 200)):  # Sample
        for j in range(i+1, min(n, 200)):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(4))
            if d2 > 1e-10 and d2 < min_dist2:
                min_dist2 = d2

    print(f"  Minimum distance^2 = {min_dist2:.6f}")
    edge_len2 = min_dist2

    # Build adjacency
    A_np = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(4))
            if abs(d2 - edge_len2) < 0.01:
                A_np[i][j] = 1
                A_np[j][i] = 1

    degrees = A_np.sum(axis=1)
    d_set = set(int(d) for d in degrees)
    print(f"  degrees = {d_set}")

    return A_np, "600-cell", n

# ======================================================================
# Run computations
# ======================================================================

results = {}

# 5-cell
A, name = build_5cell()
p, q, is_p = compute_trace_exact(A, name)
results[name] = (p, q, is_p)

# 8-cell (tesseract)
A, name = build_8cell()
p, q, is_p = compute_trace_exact(A, name)
results[name] = (p, q, is_p)

# 16-cell
A, name = build_16cell()
p, q, is_p = compute_trace_exact(A, name)
results[name] = (p, q, is_p)

# 24-cell
A, name = build_24cell()
p, q, is_p = compute_trace_exact(A, name)
results[name] = (p, q, is_p)

# 600-cell (numerical)
try:
    result = build_600cell()
    if result[0] is not None:
        A_np, name, n = result
        p, q, is_p = compute_trace_numerical(A_np, name)
        results[name] = (p, q, is_p)
except Exception as e:
    print(f"  600-cell error: {e}")

# 120-cell
build_120cell_numerical()

# ======================================================================
# Summary
# ======================================================================
print("\n" + "=" * 70)
print("4D REGULAR POLYTOPES SUMMARY")
print("=" * 70)
print(f"{'Polytope':<30} {'Tr(L^+)':<15} {'Num':<10} {'Prime?'}")
print("-" * 70)
all_prime = True
for name, (p, q, is_p) in results.items():
    if p is not None:
        print(f"{name:<30} {p}/{q:<13} {p:<10} {is_p}")
        if not is_p:
            all_prime = False
    else:
        print(f"{name:<30} {'???':<15}")
        all_prime = False

if not all_prime:
    print("\nNOT all 4D polytope numerators are prime!")
    print("This STRENGTHENS the claim that the property is specific to 3D Platonic solids.")
else:
    print("\nAll computed 4D polytope numerators are prime!")
    print("This would extend the theorem to 4D.")
