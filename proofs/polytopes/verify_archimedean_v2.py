"""
Archimedean solids Tr(L^+) -- use known Laplacian spectra from algebraic graph theory.
We have 13 Archimedean solids. Let's compute as many as we can.
"""

import numpy as np
from sympy import Matrix, Rational, isprime, simplify, sqrt
from math import gcd
from itertools import permutations

print("=" * 70)
print("ARCHIMEDEAN SOLIDS: Tr(L^+) via known eigenvalue spectra")
print("=" * 70)

def trace_from_laplacian_eigenvalues(eig_mult_pairs, name):
    """Given [(eigenvalue, multiplicity), ...] compute Tr(L^+)"""
    trace = Rational(0)
    for eig, mult in eig_mult_pairs:
        if eig != 0:
            trace += Rational(mult) / eig
    trace = simplify(trace)
    if trace.is_rational:
        p = int(trace.p)
        q = int(trace.q)
        is_p = isprime(abs(p))
        print(f"  {name}: Tr(L^+) = {p}/{q}, numerator = {p}, prime = {is_p}")
        return p, q, is_p
    else:
        print(f"  {name}: Tr(L^+) = {trace} (not rational?)")
        return None, None, None

def trace_from_adj_eigenvalues(adj_eig_mult, degree, name):
    """Given adjacency eigenvalues and degree d, Laplacian eigenvalue = d - adj_eig"""
    lap_eig_mult = [(degree - e, m) for e, m in adj_eig_mult]
    return trace_from_laplacian_eigenvalues(lap_eig_mult, name)

def build_and_compute(verts, name, expected_degree=None):
    """Build adjacency from vertices, compute numerically, rationalize"""
    n = len(verts)

    # Find edge length
    dists = {}
    for i in range(n):
        for j in range(i+1, n):
            d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
            d2_r = round(d2, 8)
            if d2_r > 1e-10:
                dists[d2_r] = dists.get(d2_r, 0) + 1

    sorted_dists = sorted(dists.keys())

    # Try each distance as edge length until we get the expected degree
    for trial_d2 in sorted_dists[:5]:
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                d2 = sum((verts[i][k] - verts[j][k])**2 for k in range(3))
                if abs(d2 - trial_d2) < 0.001:
                    A[i][j] = 1
                    A[j][i] = 1

        degrees = A.sum(axis=1)
        d_set = set(int(d) for d in degrees)

        if expected_degree and expected_degree in d_set and len(d_set) == 1:
            break

    degrees = A.sum(axis=1)
    d_set = set(int(d) for d in degrees)
    print(f"  {name}: {n} vertices, degrees = {d_set}")

    if len(d_set) != 1:
        print(f"    NOT REGULAR, skipping")
        return None, None, None

    # Compute Tr(L^+)
    D = np.diag(A.sum(axis=1))
    L = D - A
    eigs = np.linalg.eigvalsh(L)
    eigs.sort()

    trace = sum(1.0/e for e in eigs if abs(e) > 1e-8)

    # Rationalize
    best_frac = None
    best_err = 1e10
    for denom in range(1, 100000):
        numer = round(trace * denom)
        err = abs(trace - numer/denom)
        if err < best_err:
            best_err = err
            best_frac = (numer, denom)
            if err < 1e-10:
                break

    if best_frac and best_err < 1e-6:
        g = gcd(best_frac[0], best_frac[1])
        p, q = best_frac[0]//g, best_frac[1]//g
        is_p = isprime(abs(p))
        print(f"    Tr(L^+) ~ {trace:.10f} = {p}/{q} (err {best_err:.2e}), num={p}, prime={is_p}")
        return p, q, is_p

    print(f"    Could not rationalize {trace:.10f}")
    return None, None, None


results = {}

# ======================================================================
# 1. Truncated Tetrahedron (12 vertices, degree 3)
# ======================================================================
# Adjacency eigenvalues (known): 3(x1), (1+sqrt(5))/2 (x3), sqrt(5)(x3),
#   -(1+sqrt(5))/2 (x3), -3(x1), (1-sqrt(5))/2 (x3)... wait, let me just build it

verts_tt = []
# Truncated tetrahedron vertices: all even permutations of (3, 1, 1)
# with signs such that the number of minus signs is even
base = [3, 1, 1]
for p in set(permutations(base)):
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                if (s0 < 0) + (s1 < 0) + (s2 < 0) in [0, 2]:  # even number of negatives
                    v = (s0*p[0], s1*p[1], s2*p[2])
                    if v not in verts_tt:
                        verts_tt.append(v)

if len(verts_tt) == 12:
    p, q, is_p = build_and_compute(verts_tt, "Truncated Tetrahedron", 3)
    results["Truncated Tetrahedron"] = (p, q, is_p)
else:
    print(f"  Truncated Tetrahedron: got {len(verts_tt)} vertices, trying matrix construction")
    # Build directly
    n = 12
    edges = [
        (0,1),(1,2),(2,0),
        (3,4),(4,5),(5,3),
        (6,7),(7,8),(8,6),
        (9,10),(10,11),(11,9),
        (0,3),(1,6),(2,9),
        (4,7),(5,10),(8,11),
    ]
    A = np.zeros((n, n))
    for i, j in edges:
        A[i][j] = 1
        A[j][i] = 1

    D = np.diag(A.sum(axis=1))
    L = D - A
    eigs = np.linalg.eigvalsh(L)
    trace = sum(1.0/e for e in eigs if abs(e) > 1e-8)

    # Rationalize
    for denom in range(1, 10000):
        numer = round(trace * denom)
        err = abs(trace - numer/denom)
        if err < 1e-10:
            g = gcd(numer, denom)
            p, q = numer//g, denom//g
            is_p = isprime(abs(p))
            print(f"  Truncated Tetrahedron: {p}/{q}, prime={is_p}")
            results["Truncated Tetrahedron"] = (p, q, is_p)
            break

# ======================================================================
# 2. Cuboctahedron (12 vertices, degree 4)
# ======================================================================
verts_co = [
    (1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),
    (1,0,1),(1,0,-1),(-1,0,1),(-1,0,-1),
    (0,1,1),(0,1,-1),(0,-1,1),(0,-1,-1),
]
p, q, is_p = build_and_compute(verts_co, "Cuboctahedron", 4)
results["Cuboctahedron"] = (p, q, is_p)

# ======================================================================
# 3. Truncated Octahedron (24 vertices, degree 3)
# ======================================================================
# Vertices: all permutations of (0, +-1, +-2)
verts_to = []
for perm in set(permutations([0, 1, 2])):
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            v = list(perm)
            for i in range(3):
                if v[i] == 1: v[i] = s1
                elif v[i] == 2: v[i] = s2 * 2
            v = tuple(v)
            if v not in verts_to:
                verts_to.append(v)

p, q, is_p = build_and_compute(verts_to, "Truncated Octahedron", 3)
results["Truncated Octahedron"] = (p, q, is_p)

# ======================================================================
# 4. Truncated Cube (24 vertices, degree 3)
# ======================================================================
xi = 2**0.5 - 1
verts_tc = []
for perm in set(permutations([0, 1, 2])):
    vals = [xi, 1, 1]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0, 0, 0]
                v[perm[0]] = s0 * vals[0]
                v[perm[1]] = s1 * vals[1]
                v[perm[2]] = s2 * vals[2]
                v = tuple(round(x, 10) for x in v)
                if v not in verts_tc:
                    verts_tc.append(v)

p, q, is_p = build_and_compute(verts_tc, "Truncated Cube", 3)
results["Truncated Cube"] = (p, q, is_p)

# ======================================================================
# 5. Rhombicuboctahedron (24 vertices, degree 4)
# ======================================================================
s2v = 2**0.5 + 1
verts_rco = []
for perm in set(permutations([0, 1, 2])):
    vals = [1, 1, s2v]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2_ in [1, -1]:
                v = [0, 0, 0]
                v[perm[0]] = s0 * vals[0]
                v[perm[1]] = s1 * vals[1]
                v[perm[2]] = s2_ * vals[2]
                v = tuple(round(x, 10) for x in v)
                if v not in verts_rco:
                    verts_rco.append(v)

p, q, is_p = build_and_compute(verts_rco, "Rhombicuboctahedron", 4)
results["Rhombicuboctahedron"] = (p, q, is_p)

# ======================================================================
# 6. Icosidodecahedron (30 vertices, degree 4)
# ======================================================================
phi = (1 + 5**0.5) / 2
verts_id = []
# All even permutations of (0, 0, +-2*phi) -> 6 vertices... no
# Vertices: all permutations of (0, +-1, +-phi)  -- this gives 24
# Plus: (0, 0, +-phi) permutations... no

# Actually icosidodecahedron vertices at all even permutations of (0, +-1, +-phi)
# That gives 24. We need 30. The other 6 are permutations of (0, 0, +-2) -> 6
# Total: 30
verts_id = []
# Type 1: even permutations of (0, +-1, +-phi)
def even_perms_3(lst):
    a, b, c = lst
    return [(a,b,c), (b,c,a), (c,a,b)]

for ep in even_perms_3([0, 1, phi]):
    for s1 in ([1, -1] if abs(ep[1]) > 1e-10 else [1]):
        for s2 in ([1, -1] if abs(ep[2]) > 1e-10 else [1]):
            v = (ep[0], s1*abs(ep[1]), s2*abs(ep[2]))
            v = tuple(round(x, 10) for x in v)
            if v not in verts_id:
                verts_id.append(v)

# That may not be right. Let me use the midpoints-of-icosahedron-edges approach
# which worked before. Actually let me just compute it the lazy way.

# Icosidodecahedron: midpoints of icosahedron edges
ico_verts = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        ico_verts.append((0, s1, s2*phi))
        ico_verts.append((s1, s2*phi, 0))
        ico_verts.append((s2*phi, 0, s1))

ico_edges = []
for i in range(12):
    for j in range(i+1, 12):
        d2 = sum((ico_verts[i][k] - ico_verts[j][k])**2 for k in range(3))
        if abs(d2 - 4) < 0.01:
            ico_edges.append((i, j))

verts_id = []
for i, j in ico_edges:
    mid = tuple((ico_verts[i][k] + ico_verts[j][k])/2 for k in range(3))
    mid = tuple(round(x, 10) for x in mid)
    if mid not in verts_id:
        verts_id.append(mid)

print(f"\n  Icosidodecahedron: {len(verts_id)} vertices")
if len(verts_id) == 30:
    p, q, is_p = build_and_compute(verts_id, "Icosidodecahedron", 4)
    results["Icosidodecahedron"] = (p, q, is_p)

# ======================================================================
# 7. Truncated Icosahedron (60 vertices, degree 3) -- soccer ball
# ======================================================================
def even_perms_3_coords(a, b, c):
    return [(a,b,c), (b,c,a), (c,a,b)]

verts_ti = []
phi_v = phi
phi2 = phi + 1  # phi^2
phi3 = phi2 * phi  # phi^3

# From Wikipedia: Truncated icosahedron Cartesian coordinates
# Even permutations of:
# (0, +-1, +-3*phi)
# (+-2, +-(1+2*phi), +-phi)
# (+-1, +-(2+phi), +-2*phi)

for base in [(0, 1, 3*phi_v), (2, 1+2*phi_v, phi_v), (1, 2+phi_v, 2*phi_v)]:
    for ep in even_perms_3_coords(*base):
        for s0 in ([1, -1] if abs(ep[0]) > 1e-10 else [1]):
            for s1 in ([1, -1] if abs(ep[1]) > 1e-10 else [1]):
                for s2 in ([1, -1] if abs(ep[2]) > 1e-10 else [1]):
                    v = (s0*abs(ep[0]), s1*abs(ep[1]), s2*abs(ep[2]))
                    v = tuple(round(x, 10) for x in v)
                    if v not in verts_ti:
                        verts_ti.append(v)

print(f"\n  Truncated Icosahedron: {len(verts_ti)} vertices")
if len(verts_ti) == 60:
    p, q, is_p = build_and_compute(verts_ti, "Truncated Icosahedron", 3)
    results["Truncated Icosahedron"] = (p, q, is_p)

# ======================================================================
# 8. Snub Cube (24 vertices, degree 5)
# ======================================================================
# Snub cube has chiral coordinates -- more complex
# Using tribonacci constant
t = 1.8392867552141612  # real root of t^3 = t^2 + t + 1
verts_sc = []
# Even permutations of (+-1, +-1/t, +-t)
for ep in even_perms_3_coords(1, 1/t, t):
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                if (s0 < 0) + (s1 < 0) + (s2 < 0) in [0, 2]:  # even parity
                    v = (s0*ep[0], s1*ep[1], s2*ep[2])
                    v = tuple(round(x, 10) for x in v)
                    if v not in verts_sc:
                        verts_sc.append(v)

print(f"\n  Snub Cube attempt: {len(verts_sc)} vertices (need 24)")
if len(verts_sc) == 24:
    p, q, is_p = build_and_compute(verts_sc, "Snub Cube", 5)
    results["Snub Cube"] = (p, q, is_p)

# ======================================================================
# 9. Truncated Dodecahedron (60 vertices, degree 3)
# ======================================================================
# Even permutations of:
# (0, +-1/phi, +-(2+phi))
# (+-1/phi, +-phi, +-2*phi)
# (+-phi, +-2, +-(phi+1))

verts_td = []
inv_phi = 1/phi

for base in [(0, inv_phi, 2+phi_v), (inv_phi, phi_v, 2*phi_v), (phi_v, 2, phi_v+1)]:
    for ep in even_perms_3_coords(*base):
        for s0 in ([1, -1] if abs(ep[0]) > 1e-10 else [1]):
            for s1 in ([1, -1] if abs(ep[1]) > 1e-10 else [1]):
                for s2 in ([1, -1] if abs(ep[2]) > 1e-10 else [1]):
                    v = (s0*abs(ep[0]), s1*abs(ep[1]), s2*abs(ep[2]))
                    v = tuple(round(x, 10) for x in v)
                    if v not in verts_td:
                        verts_td.append(v)

print(f"\n  Truncated Dodecahedron: {len(verts_td)} vertices (need 60)")
if len(verts_td) == 60:
    p, q, is_p = build_and_compute(verts_td, "Truncated Dodecahedron", 3)
    results["Truncated Dodecahedron"] = (p, q, is_p)

# ======================================================================
# 10. Rhombicosidodecahedron (60 vertices, degree 4)
# ======================================================================
# Even permutations of:
# (+-1, +-1, +-phi^3)
# (+-phi^2, +-phi, +-2*phi)
# (+-2+phi, 0, +-phi^2)

phi3v = phi**3
phi2v = phi**2
verts_rid = []

for base in [(1, 1, phi3v), (phi2v, phi_v, 2*phi_v), (2+phi_v, 0, phi2v)]:
    for ep in even_perms_3_coords(*base):
        for s0 in ([1, -1] if abs(ep[0]) > 1e-10 else [1]):
            for s1 in ([1, -1] if abs(ep[1]) > 1e-10 else [1]):
                for s2 in ([1, -1] if abs(ep[2]) > 1e-10 else [1]):
                    v = (s0*abs(ep[0]), s1*abs(ep[1]), s2*abs(ep[2]))
                    v = tuple(round(x, 10) for x in v)
                    if v not in verts_rid:
                        verts_rid.append(v)

print(f"\n  Rhombicosidodecahedron: {len(verts_rid)} vertices (need 60)")
if len(verts_rid) == 60:
    p, q, is_p = build_and_compute(verts_rid, "Rhombicosidodecahedron", 4)
    results["Rhombicosidodecahedron"] = (p, q, is_p)

# ======================================================================
# 11. Truncated Cuboctahedron / Great Rhombicuboctahedron (48 vertices, degree 3)
# ======================================================================
# Vertices: all permutations of (+-1, +-(1+sqrt(2)), +-(1+2*sqrt(2)))
s2_val = 2**0.5
verts_tco = []
for perm in set(permutations([0, 1, 2])):
    vals = [1, 1+s2_val, 1+2*s2_val]
    for s0 in [1, -1]:
        for s1 in [1, -1]:
            for s2_ in [1, -1]:
                v = [0, 0, 0]
                v[perm[0]] = s0 * vals[0]
                v[perm[1]] = s1 * vals[1]
                v[perm[2]] = s2_ * vals[2]
                v = tuple(round(x, 10) for x in v)
                if v not in verts_tco:
                    verts_tco.append(v)

print(f"\n  Truncated Cuboctahedron: {len(verts_tco)} vertices (need 48)")
if len(verts_tco) == 48:
    p, q, is_p = build_and_compute(verts_tco, "Truncated Cuboctahedron", 3)
    results["Truncated Cuboctahedron"] = (p, q, is_p)

# ======================================================================
# 12. Truncated Icosidodecahedron (120 vertices, degree 3)
# ======================================================================
# Very large. Known coords: all even permutations of
# (+-1/phi, +-1/phi, +-(3+phi))
# (+-2/phi, +-phi, +-(1+2*phi))
# (+-1/phi, +-phi^2, +-(3*phi-1))
# (+-(2*phi-1), +-2, +-(2+phi))
# (+-phi, +-3, +-2*phi)

verts_tid = []
for base in [
    (inv_phi, inv_phi, 3+phi_v),
    (2*inv_phi, phi_v, 1+2*phi_v),
    (inv_phi, phi2v, 3*phi_v-1),
    (2*phi_v-1, 2, 2+phi_v),
    (phi_v, 3, 2*phi_v),
]:
    for ep in even_perms_3_coords(*base):
        for s0 in ([1, -1] if abs(ep[0]) > 1e-10 else [1]):
            for s1 in ([1, -1] if abs(ep[1]) > 1e-10 else [1]):
                for s2 in ([1, -1] if abs(ep[2]) > 1e-10 else [1]):
                    v = (s0*abs(ep[0]), s1*abs(ep[1]), s2*abs(ep[2]))
                    v = tuple(round(x, 10) for x in v)
                    if v not in verts_tid:
                        verts_tid.append(v)

print(f"\n  Truncated Icosidodecahedron: {len(verts_tid)} vertices (need 120)")
if len(verts_tid) == 120:
    p, q, is_p = build_and_compute(verts_tid, "Truncated Icosidodecahedron", 3)
    results["Truncated Icosidodecahedron"] = (p, q, is_p)

# ======================================================================
# 13. Snub Dodecahedron (60 vertices, degree 5) -- chiral, very complex
# Skipping due to coordinate complexity
# ======================================================================
print("\n  Snub Dodecahedron: skipped (chiral, complex coordinates)")

# ======================================================================
# Summary
# ======================================================================
print("\n" + "=" * 70)
print("ARCHIMEDEAN SOLIDS SUMMARY")
print("=" * 70)
print(f"{'#':<4} {'Solid':<35} {'Tr(L^+)':<15} {'Num':<10} {'Prime?'}")
print("-" * 75)

prime_count = 0
composite_count = 0
total = 0

for idx, (name, (p, q, is_p)) in enumerate(results.items(), 1):
    if p is not None:
        total += 1
        print(f"{idx:<4} {name:<35} {p}/{q:<13} {p:<10} {is_p}")
        if is_p:
            prime_count += 1
        else:
            composite_count += 1
    else:
        print(f"{idx:<4} {name:<35} {'FAILED':<15}")

print(f"\nTotal computed: {total}")
print(f"Prime numerators: {prime_count}")
print(f"Composite numerators: {composite_count}")

if composite_count > 0:
    print("\nVERDICT: NOT all Archimedean numerators are prime.")
    print("The all-prime property is SPECIFIC to the 5 Platonic solids.")
