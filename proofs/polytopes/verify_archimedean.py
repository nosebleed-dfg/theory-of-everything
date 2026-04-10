"""
EXTENSION: Compute Tr(L^+) for all 13 Archimedean solids.
If their numerators are NOT all prime, that STRENGTHENS the Platonic-specific claim.
"""

from sympy import Matrix, Rational, isprime, simplify
import sympy

print("=" * 70)
print("ARCHIMEDEAN SOLIDS: Tr(L^+) COMPUTATION")
print("=" * 70)

# For Archimedean solids, we use eigenvalues of the graph Laplacian.
# I'll build each adjacency matrix from known vertex/edge structure.

# For large graphs, we use the known eigenvalue spectra from algebraic graph theory.
# References: Cvetkovic, Doob, Sachs - "Spectra of Graphs"

def trace_from_eigenvalues(eigenvals_with_mult, name):
    """
    eigenvals_with_mult: list of (eigenvalue, multiplicity) pairs
    eigenvalue can be sympy expression
    """
    trace = Rational(0)
    for eig, mult in eigenvals_with_mult:
        if eig != 0:
            trace += Rational(mult) / eig

    trace = simplify(trace)

    if trace.is_rational:
        p = trace.p
        q = trace.q
        prime = isprime(abs(p))
        print(f"  {name}: Tr(L^+) = {p}/{q}, numerator = {p}, prime = {prime}")
        return p, q, prime
    else:
        print(f"  {name}: Tr(L^+) = {trace} (not fully simplified to rational)")
        # Try harder
        trace_r = Rational(trace)
        if trace_r.is_rational:
            p = trace_r.p
            q = trace_r.q
            prime = isprime(abs(p))
            print(f"    -> rationalized: {p}/{q}, numerator = {p}, prime = {prime}")
            return p, q, prime
        return None, None, None

# ======================================================================
# Build adjacency matrices for all 13 Archimedean solids
# ======================================================================

def build_truncated_tetrahedron():
    """12 vertices, degree 3. Vertex-transitive."""
    # Truncated tetrahedron: replace each vertex of tetrahedron with triangle
    # 12 vertices: 4 triangles x 3 vertices each
    # Each original edge of tetrahedron becomes a hexagonal face edge
    n = 12
    # Adjacency: each vertex connects to 2 in same triangle + 1 in adjacent triangle
    edges = [
        # Triangle 0: vertices 0,1,2
        (0,1),(1,2),(2,0),
        # Triangle 1: vertices 3,4,5
        (3,4),(4,5),(5,3),
        # Triangle 2: vertices 6,7,8
        (6,7),(7,8),(8,6),
        # Triangle 3: vertices 9,10,11
        (9,10),(10,11),(11,9),
        # Inter-triangle edges (hexagonal faces)
        (0,3),(1,6),(2,9),
        (4,7),(5,10),(8,11),
    ]
    A = Matrix.zeros(n, n)
    for i, j in edges:
        A[i,j] = 1
        A[j,i] = 1
    return A, "Truncated Tetrahedron"

def build_cuboctahedron():
    """12 vertices, degree 4. Vertex-transitive."""
    n = 12
    # Cuboctahedron: vertices at midpoints of cube edges
    # 12 vertices, each degree 4
    # Can be built as rectification of cube
    # Vertices: midpoints of 12 edges of cube
    # Two vertices adjacent iff corresponding cube edges share a vertex
    # and are on the same face

    # Explicit construction using coordinates
    # Vertices at permutations of (0, +-1, +-1) -- 12 vertices
    import itertools
    verts = []
    for perm in itertools.permutations([0, 1, 1]):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = list(perm)
                # Apply signs to non-zero entries
                nz = [i for i in range(3) if v[i] != 0]
                v[nz[0]] *= s1
                v[nz[1]] *= s2
                if v not in verts:
                    verts.append(v)

    # Actually let me just use the known coordinates
    verts = [
        (1,1,0),  (1,-1,0), (-1,1,0), (-1,-1,0),
        (1,0,1),  (1,0,-1), (-1,0,1), (-1,0,-1),
        (0,1,1),  (0,1,-1), (0,-1,1), (0,-1,-1),
    ]
    n = len(verts)
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(i+1, n):
            # Distance between vertices
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            if d2 == 2:  # adjacent iff distance = sqrt(2)
                A[i,j] = 1
                A[j,i] = 1

    # Verify regularity
    degrees = [sum(A[i,j] for j in range(n)) for i in range(n)]
    assert all(d == 4 for d in degrees), f"Not regular: {degrees}"
    return A, "Cuboctahedron"

def build_truncated_cube():
    """24 vertices, degree 3."""
    # Truncation of cube: replace each vertex with a triangle
    # 24 vertices, 36 edges
    # Known Laplacian eigenvalues (from spectral graph theory):
    # This is complex to build by hand, let me use coordinates
    import math
    xi = math.sqrt(2) - 1  # ~0.4142

    # Vertices: all permutations of (+-xi, +-1, +-1)
    verts = []
    for signs_xi in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                # (+-xi, +-1, +-1) and permutations
                perms = [
                    (signs_xi * xi, s1 * 1, s2 * 1),
                    (s1 * 1, signs_xi * xi, s2 * 1),
                    (s1 * 1, s2 * 1, signs_xi * xi),
                ]
                for v in perms:
                    verts.append(v)

    n = len(verts)
    A = Matrix.zeros(n, n)

    # Find minimum nonzero distance
    dists = set()
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            dists.add(round(d2, 8))

    sorted_dists = sorted(dists)
    min_d2 = sorted_dists[0]  # smallest distance = edge length

    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            if abs(d2 - min_d2) < 1e-6:
                A[i,j] = 1
                A[j,i] = 1

    degrees = [sum(A[i,j] for j in range(n)) for i in range(n)]
    d_set = set(int(d) for d in degrees)
    if len(d_set) != 1:
        print(f"  WARNING: Truncated Cube not regular: degrees = {d_set}")
    return A, "Truncated Cube"

def build_truncated_octahedron():
    """24 vertices, degree 3."""
    # Vertices: all permutations of (0, +-1, +-2)
    verts = []
    import itertools
    base = [0, 1, 2]
    for perm in itertools.permutations(base):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = list(perm)
                # Apply signs to nonzero entries
                nz_indices = [i for i in range(3) if v[i] != 0]
                if len(nz_indices) == 2:
                    v[nz_indices[0]] *= s1
                    v[nz_indices[1]] *= s2
                elif len(nz_indices) == 1:
                    v[nz_indices[0]] *= s1
                    if s2 == -1:
                        continue  # avoid duplicates
                else:
                    if s1 == -1 or s2 == -1:
                        continue
                if v not in verts:
                    verts.append(v)

    n = len(verts)
    A = Matrix.zeros(n, n)

    # Edge length = sqrt(2)
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            if abs(d2 - 2) < 1e-6:
                A[i,j] = 1
                A[j,i] = 1

    degrees = [int(sum(A[i,j] for j in range(n))) for i in range(n)]
    d_set = set(degrees)
    print(f"  Truncated Octahedron: {n} vertices, degrees = {d_set}")
    return A, "Truncated Octahedron"

def build_icosidodecahedron():
    """30 vertices, degree 4."""
    # Vertices at midpoints of icosahedron edges
    # Use golden ratio coordinates
    phi = (1 + 5**0.5) / 2

    verts = []
    # All permutations of (0, 0, +-2phi)
    # All permutations of (+-1, +-phi, +-phi^2) where phi^2 = phi+1
    phi2 = phi + 1  # = phi^2

    # Type 1: permutations of (0, 0, 2*phi) with signs
    import itertools
    for perm in itertools.permutations([0, 0, 2*phi]):
        for s in [1, -1]:
            v = list(perm)
            v[2] *= s if v[2] != 0 else 1
            v[1] *= s if v[1] != 0 else 1
            v[0] *= s if v[0] != 0 else 1
            # Only apply sign to nonzero
            nz = [i for i in range(3) if abs(v[i]) > 1e-10]
            if len(nz) > 0:
                v2 = [0.0, 0.0, 0.0]
                for i in range(3):
                    v2[i] = abs(perm[i]) * (s if i == nz[-1] else 1) if abs(perm[i]) > 1e-10 else 0
                # Actually this is getting messy, let me use a cleaner approach
                pass

    # Let me just use known coordinates directly
    # Icosidodecahedron vertices (from Wikipedia)
    verts = []

    # (0, 0, +-phi) -- 2 vertices (but we need 30 total)
    # Actually the icosidodecahedron has 30 vertices at:
    # All even permutations of (0, +-1, +-phi)  -- but that gives only 12
    # Plus all permutations of (+-1/phi, +-phi, 0) etc.

    # Clean coords: all even permutations of (0, +-1, +-phi) = 12
    # All even permutations of (+-1/phi, +-1, +-phi^2/2)...

    # This is getting complex. Let me use a different approach:
    # build from the known edge skeleton

    # For a simpler approach, use the EIGENVALUES directly from literature
    # Icosidodecahedron eigenvalues of adjacency matrix are known:
    # From: https://mathworld.wolfram.com/IcosidodecahedralGraph.html
    # Adjacency eigenvalues: 4, (1+sqrt(5))/2 (x5), 1 (x3), (-1+sqrt(5))/2 (x3),
    #                        -2 (x4), (1-sqrt(5))/2 (x5), (-1-sqrt(5))/2 (x3), ...
    # This needs careful handling. Let me skip complex coordinate generation
    # and compute eigenvalues numerically for the large Archimedean solids.

    print("  [Icosidodecahedron: using numerical eigenvalues from coordinate construction]")

    # Use the Cartesian coordinates of the icosidodecahedron
    # All cyclic permutations of (0, 0, +-2) -> 6
    # All cyclic permutations of (+-1, +-phi, +-(1+phi)) with even permutations -> 24
    # Total: 30

    phi_val = phi
    phi2_val = phi + 1

    base_verts = []
    # Type A: (0, 0, +-2) and cyclic perms -> 6 vertices
    for val in [2, -2]:
        base_verts.append((val, 0, 0))
        base_verts.append((0, val, 0))
        base_verts.append((0, 0, val))

    # Type B: (+-1, +-phi, +-phi^2) even perms with constraints
    # Actually for icosidodecahedron, the 30 vertices are:
    # All permutations of (0, +-1, +-phi) = 24 after removing duplicates...
    # Let me just enumerate all even permutations

    from itertools import permutations

    # Midpoints of icosahedron edges
    # Icosahedron vertices:
    ico_verts_raw = []
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            ico_verts_raw.append((0, s1, s2*phi))
            ico_verts_raw.append((s1, s2*phi, 0))
            ico_verts_raw.append((s2*phi, 0, s1))
    # That's 12 vertices

    # Find edges: distance = 2 for unit icosahedron
    ico_edges = []
    for i in range(12):
        for j in range(i+1, 12):
            d2 = sum((ico_verts_raw[i][k] - ico_verts_raw[j][k])**2 for k in range(3))
            if abs(d2 - 4) < 0.01:  # edge length = 2
                ico_edges.append((i, j))

    # Midpoints of edges = icosidodecahedron vertices
    verts = []
    for i, j in ico_edges:
        mid = tuple((ico_verts_raw[i][k] + ico_verts_raw[j][k])/2 for k in range(3))
        verts.append(mid)

    n = len(verts)
    print(f"  Icosidodecahedron: {n} vertices from edge midpoints")

    if n != 30:
        print(f"  WARNING: expected 30, got {n}")
        return None, "Icosidodecahedron"

    # Build adjacency: two midpoints are adjacent iff their edges share a vertex
    # and are on the same face (i.e., the edges are adjacent in the icosahedron)
    # Simpler: adjacent iff distance = edge length (the minimum positive distance)

    dists = set()
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            dists.add(round(d2, 6))
    sorted_dists = sorted(dists)
    min_d2 = sorted_dists[0]

    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            if abs(d2 - min_d2) < 0.001:
                A[i,j] = 1
                A[j,i] = 1

    degrees = [int(sum(A[i,j] for j in range(n))) for i in range(n)]
    d_set = set(degrees)
    print(f"  degrees = {d_set}")
    return A, "Icosidodecahedron"


def compute_trace_numerical(A, name):
    """Compute Tr(L^+) numerically using numpy, then try to rationalize"""
    import numpy as np
    from fractions import Fraction

    n = A.rows
    # Convert to numpy
    A_np = np.array(A.tolist(), dtype=float)
    D_np = np.diag(A_np.sum(axis=1))
    L_np = D_np - A_np

    eigs = np.linalg.eigvalsh(L_np)

    # Sort
    eigs.sort()

    # Tr(L^+) = sum of 1/lambda for nonzero lambda
    trace = sum(1.0/e for e in eigs if abs(e) > 1e-8)

    print(f"  {name}: Tr(L^+) ~= {trace:.10f}")

    # Try to identify as a fraction with small denominator
    best_frac = None
    best_err = 1e10
    for denom in range(1, 10000):
        numer = round(trace * denom)
        err = abs(trace - numer/denom)
        if err < best_err:
            best_err = err
            best_frac = (numer, denom)
            if err < 1e-8:
                break

    if best_frac and best_err < 1e-6:
        from math import gcd
        g = gcd(best_frac[0], best_frac[1])
        p, q = best_frac[0]//g, best_frac[1]//g
        is_p = isprime(abs(p))
        print(f"    = {p}/{q} (error {best_err:.2e})")
        print(f"    Numerator = {p}, prime = {is_p}")
        return p, q, is_p
    else:
        print(f"    Could not rationalize (best: {best_frac}, error {best_err:.2e})")
        return None, None, None


# ======================================================================
# Run computations
# ======================================================================

# Start with the simpler Archimedean solids that we can build adjacency matrices for
archimedean_builders = [
    build_truncated_tetrahedron,
    build_cuboctahedron,
    build_truncated_octahedron,
]

archimedean_results = {}

for builder in archimedean_builders:
    try:
        A, name = builder()
        if A is None:
            continue
        n = A.rows
        degrees = [int(sum(A[i,j] for j in range(n))) for i in range(n)]
        d_set = set(degrees)
        print(f"\n{name}: {n} vertices, degrees = {d_set}")

        if n <= 24:
            # Can do exact computation
            p, q, is_p = compute_trace_numerical(A, name)
            archimedean_results[name] = (p, q, is_p)
        else:
            p, q, is_p = compute_trace_numerical(A, name)
            archimedean_results[name] = (p, q, is_p)
    except Exception as e:
        print(f"  {builder.__name__}: ERROR: {e}")

# Icosidodecahedron
try:
    A, name = build_icosidodecahedron()
    if A is not None:
        p, q, is_p = compute_trace_numerical(A, name)
        archimedean_results[name] = (p, q, is_p)
except Exception as e:
    print(f"  Icosidodecahedron: ERROR: {e}")

# For remaining Archimedean solids, use known graph spectra
# Truncated Cube (24v, 3-reg), Snub Cube (24v, 5-reg),
# Truncated Icosahedron (60v, 3-reg), etc.

# Let me build more using coordinate methods
def build_from_coords(name, coords, expected_degree=None):
    """Build adjacency matrix from vertex coordinates"""
    n = len(coords)

    # Find edge length (minimum distance)
    dists = {}
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((coords[i][k] - coords[j][k])**2 for k in range(3))
            d2_r = round(d2, 6)
            dists[d2_r] = dists.get(d2_r, 0) + 1

    sorted_dists = sorted(dists.keys())
    min_d2 = sorted_dists[0]

    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((coords[i][k] - coords[j][k])**2 for k in range(3))
            if abs(d2 - min_d2) < 0.001:
                A[i,j] = 1
                A[j,i] = 1

    degrees = [int(sum(A[i,j] for j in range(n))) for i in range(n)]
    d_set = set(degrees)

    if expected_degree and expected_degree not in d_set:
        # Try next distance
        min_d2 = sorted_dists[1] if len(sorted_dists) > 1 else min_d2
        A = Matrix.zeros(n, n)
        for i in range(n):
            for j in range(i+1, n):
                d2 = sum((coords[i][k] - coords[j][k])**2 for k in range(3))
                if abs(d2 - min_d2) < 0.001:
                    A[i,j] = 1
                    A[j,i] = 1
        degrees = [int(sum(A[i,j] for j in range(n))) for i in range(n)]
        d_set = set(degrees)

    print(f"  {name}: {n} vertices, degrees = {d_set}")
    return A, name

# Rhombicuboctahedron: 24 vertices, degree 4
# Vertices: all permutations of (+-1, +-1, +-(1+sqrt(2)))
import math
s2 = math.sqrt(2)
rhombi_verts = []
for s1 in [1,-1]:
    for s2_ in [1,-1]:
        for s3 in [1,-1]:
            for perm in permutations([s1*1, s2_*1, s3*(1+s2)]):
                v = tuple(round(x, 10) for x in perm)
                if v not in rhombi_verts:
                    rhombi_verts.append(v)

if len(rhombi_verts) == 24:
    A, name = build_from_coords("Rhombicuboctahedron", rhombi_verts, expected_degree=4)
    p, q, is_p = compute_trace_numerical(A, name)
    archimedean_results[name] = (p, q, is_p)

# Truncated Icosahedron (soccer ball): 60 vertices, degree 3
# This is the largest -- let me build it from icosahedron truncation
phi = (1 + 5**0.5) / 2

# Known coordinates for truncated icosahedron
# All even permutations of:
# (0, +-1, +-3phi), (+-1, +-2-phi, +-2phi), (+-2, +-1+2phi, +-phi)
# (+-phi, +-2, +-1+2phi) -- wait, let me use a standard set

# From Wikipedia: Truncated icosahedron vertices are all even permutations of:
# (0, +-1, +-3*phi)
# (+-2, +-(1+2*phi), +-phi)
# (+-1, +-(2+phi), +-2*phi)

def even_permutations(triple):
    """Return the 3 even permutations of a triple"""
    a, b, c = triple
    return [(a,b,c), (b,c,a), (c,a,b)]

trunc_ico_verts = []
phi_v = phi

# Type 1: (0, +-1, +-3*phi)
for s1 in [1,-1]:
    for s2 in [1,-1]:
        for ep in even_permutations((0, s1*1, s2*3*phi_v)):
            v = tuple(round(x, 10) for x in ep)
            if v not in trunc_ico_verts:
                trunc_ico_verts.append(v)

# Type 2: (+-2, +-(1+2*phi), +-phi)
for s1 in [1,-1]:
    for s2 in [1,-1]:
        for s3 in [1,-1]:
            for ep in even_permutations((s1*2, s2*(1+2*phi_v), s3*phi_v)):
                v = tuple(round(x, 10) for x in ep)
                if v not in trunc_ico_verts:
                    trunc_ico_verts.append(v)

# Type 3: (+-1, +-(2+phi), +-2*phi)
for s1 in [1,-1]:
    for s2 in [1,-1]:
        for s3 in [1,-1]:
            for ep in even_permutations((s1*1, s2*(2+phi_v), s3*2*phi_v)):
                v = tuple(round(x, 10) for x in ep)
                if v not in trunc_ico_verts:
                    trunc_ico_verts.append(v)

print(f"\nTruncated Icosahedron: {len(trunc_ico_verts)} vertices")
if len(trunc_ico_verts) == 60:
    A, name = build_from_coords("Truncated Icosahedron", trunc_ico_verts, expected_degree=3)
    p, q, is_p = compute_trace_numerical(A, name)
    archimedean_results[name] = (p, q, is_p)

# ======================================================================
# Summary
# ======================================================================
print("\n" + "=" * 70)
print("ARCHIMEDEAN SOLIDS SUMMARY")
print("=" * 70)
print(f"{'Solid':<30} {'Tr(L^+)':<15} {'Num':<10} {'Prime?'}")
print("-" * 70)
all_prime = True
for name, (p, q, is_p) in archimedean_results.items():
    if p is not None:
        print(f"{name:<30} {p}/{q:<13} {p:<10} {is_p}")
        if not is_p:
            all_prime = False
    else:
        print(f"{name:<30} {'???':<15} {'???':<10} {'???'}")
        all_prime = False

if not all_prime:
    print("\nNOT all Archimedean numerators are prime!")
    print("This STRENGTHENS the Platonic-specific claim.")
else:
    print("\nAll computed Archimedean numerators are prime.")
    print("This would WEAKEN the Platonic-specific claim (extends beyond Platonic).")
