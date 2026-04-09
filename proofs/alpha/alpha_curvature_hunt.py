#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ALPHA_CURVATURE_HUNT — hunts the 1.89 ppb gap between algebraic and spectral roads to 1/alpha
nos3bl33d

Curvature, volume, conformal coupling, heat kernel, winding, and Regge corrections.
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from mpmath import mp, mpf, pi as mpi, phi as mphi, sqrt as msqrt, sin as msin
from mpmath import cos as mcos, exp as mexp, gamma as mgamma, log as mlog
from mpmath import fsum, power, nstr, acos, atan
from scipy import sparse
from scipy.sparse.linalg import eigsh
import itertools
import time

# Set very high precision
mp.dps = 80  # 80 decimal places

# ============================================================================
# FUNDAMENTAL CONSTANTS (high precision)
# ============================================================================
phi = mphi  # golden ratio (1+sqrt(5))/2
PI = mpi
NIST_INV_ALPHA = mpf('137.035999084')  # NIST 2018 CODATA
NIST_UNCERTAINTY = mpf('0.000000021')

# Road 1 algebraic formula
def road1_inv_alpha():
    """(20*phi^6 - 30/(2*pi)^3) / phi^2 * (1 + 1/(2*phi^27))"""
    base = (20 * phi**6 - mpf(30) / (2*PI)**3) / phi**2
    correction = 1 + 1 / (2 * phi**27)
    return base * correction

# Road 1 base (without correction)
def road1_base():
    return (20 * phi**6 - mpf(30) / (2*PI)**3) / phi**2

ROAD1 = road1_inv_alpha()
ROAD1_BASE = road1_base()

print("=" * 80)
print("FINE STRUCTURE CONSTANT GAP HUNTER")
print("=" * 80)
print()
print(f"phi           = {nstr(phi, 30)}")
print(f"Road 1 base   = {nstr(ROAD1_BASE, 20)}")
print(f"Road 1 (full) = {nstr(ROAD1, 15)}")
print(f"NIST          = {nstr(NIST_INV_ALPHA, 15)}")
print(f"Road 1 - NIST = {nstr(ROAD1 - NIST_INV_ALPHA, 6)} ({nstr((ROAD1 - NIST_INV_ALPHA)/NIST_INV_ALPHA * 1e9, 4)} ppb)")
print()

# ============================================================================
# SECTION 0: BUILD THE 120-CELL (CORRECT CONSTRUCTION)
# ============================================================================
print("=" * 80)
print("SECTION 0: BUILDING THE 120-CELL")
print("=" * 80)

def build_120_cell_vertices():
    """
    Build 600 vertices of the 120-cell in R^4, using the standard construction.

    The 600 vertices of the 120-cell with circumradius 2*sqrt(2) are:
    All permutations and sign changes of:
      (0, 0, 2, 2)                    -- 24
      (1, 1, 1, sqrt(5))              -- 32 (with even # of minus on non-sqrt5)
      (phi^-2, phi, phi, phi)          -- 32
      (phi^-1, phi^-1, phi^-1, phi^2) -- 32
    All EVEN permutations with ALL sign changes of:
      (0, phi^-1, 1, phi^2)           -- 192
      (0, phi, phi^2, 1/phi^2...

    Actually, let me use the known coordinate list.
    The 120-cell vertices are the 600 points obtained from:

    24 from all permutations of (0,0,+/-2,+/-2)
    64 from all sign changes of (+/-1,+/-1,+/-1,+/-sqrt(5))
    64 from all sign changes of (+/-phi, +/-phi, +/-phi, +/-phi^-2)
    64 from all sign changes of (+/-phi^2, +/-phi^-1, +/-phi^-1, +/-phi^-1)

    Plus 384 from even permutations of all sign changes of:
    (0, +/-phi^-1, +/-phi, +/-sqrt(5))
    (0, +/-phi^-2, +/-1, +/-phi^2)
    (+/-phi^-1, +/-1, +/-phi, +/-2)
    """
    p = (1 + np.sqrt(5)) / 2  # phi
    q = 1 / p               # 1/phi = phi-1
    p2 = p * p               # phi^2
    q2 = q * q               # 1/phi^2
    s5 = np.sqrt(5)

    verts = set()

    def roundv(v, decimals=10):
        return tuple(round(x, decimals) for x in v)

    # Helper: all permutations of a tuple (handling duplicates from repeated values)
    def all_perms(t):
        if len(t) <= 1:
            yield t
            return
        seen = set()
        for i, elem in enumerate(t):
            if elem in seen:
                continue
            seen.add(elem)
            rest = t[:i] + t[i+1:]
            for perm_rest in all_perms(rest):
                yield (elem,) + perm_rest

    # Helper: all sign changes of a tuple
    def all_signs(t):
        for signs in itertools.product([1, -1], repeat=len(t)):
            yield tuple(s * x for s, x in zip(signs, t))

    # Helper: all permutations with all sign changes
    def all_perms_signs(t):
        for perm in all_perms(t):
            for sv in all_signs(perm):
                verts.add(roundv(sv))

    # Helper: even permutations of 4 elements
    def even_perms_4():
        """Return even permutation indices for S_4."""
        # A_4 has 12 elements
        return [
            (0,1,2,3), (0,2,3,1), (0,3,1,2),
            (1,0,3,2), (1,2,0,3), (1,3,2,0),
            (2,0,1,3), (2,1,3,0), (2,3,0,1),
            (3,0,2,1), (3,1,0,2), (3,2,1,0),
        ]

    def even_perms_all_signs(t):
        for perm_idx in even_perms_4():
            perm = tuple(t[i] for i in perm_idx)
            for sv in all_signs(perm):
                verts.add(roundv(sv))

    # GROUP 1: All permutations of (0, 0, 2, 2) with all sign changes -> 24
    all_perms_signs((0.0, 0.0, 2.0, 2.0))

    # GROUP 2: All sign changes of all permutations of (1, 1, 1, sqrt(5)) -> 64
    # Actually: all permutations (4 of them since 3 repeated) times 16 sign combos = 64
    all_perms_signs((1.0, 1.0, 1.0, s5))

    # GROUP 3: All sign changes of all permutations of (phi^-2, phi, phi, phi) -> 64
    all_perms_signs((q2, p, p, p))

    # GROUP 4: All sign changes of all permutations of (phi^-1, phi^-1, phi^-1, phi^2) -> 64
    all_perms_signs((q, q, q, p2))

    # GROUP 5: Even permutations with all sign changes -> 384 total
    # These come from 8 base types, each giving 192, but with overlaps...
    # Actually: three types, each giving 12 even perms * 16 signs = 192, minus overlaps
    # Let me just use the correct three orbits:

    # (0, phi^-1, phi, sqrt(5)) -- 192 vertices
    even_perms_all_signs((0.0, q, p, s5))

    # (0, phi^-2, 1, phi^2) -- 192 vertices
    even_perms_all_signs((0.0, q2, 1.0, p2))

    # (phi^-1, 1, phi, 2) -- 192 vertices
    even_perms_all_signs((q, 1.0, p, 2.0))

    result = np.array(sorted(verts))
    return result


def build_120_cell_graph(vertices, tol=0.005):
    """
    Build adjacency for the 120-cell. Each vertex has degree 4.
    With circumradius 2*sqrt(2), edge length = 3 - sqrt(5) = 2/phi^2 ~ 0.764.
    """
    n = len(vertices)
    p = (1 + np.sqrt(5)) / 2
    expected_edge = 3.0 - np.sqrt(5)  # 2/phi^2 = 3 - sqrt(5) ~ 0.764
    expected_sq = expected_edge ** 2

    # Use vectorized distance computation
    from scipy.spatial.distance import cdist
    dist_sq_matrix = cdist(vertices, vertices, 'sqeuclidean')

    adjacency = [[] for _ in range(n)]
    edge_count = 0

    for i in range(n):
        for j in range(i + 1, n):
            if abs(dist_sq_matrix[i, j] - expected_sq) < tol:
                adjacency[i].append(j)
                adjacency[j].append(i)
                edge_count += 1

    return adjacency, edge_count


print("Building 120-cell vertices...")
t0 = time.time()
vertices_120 = build_120_cell_vertices()
print(f"  Vertices: {len(vertices_120)} (expect 600)")

# Verify: all vertices should have the same norm (circumradius = 2*sqrt(2))
norms = np.linalg.norm(vertices_120, axis=1)
print(f"  Norm range: {norms.min():.6f} to {norms.max():.6f} (expect {2*np.sqrt(2):.6f})")

print("Building adjacency graph...")
adj_120, n_edges = build_120_cell_graph(vertices_120)
degrees = [len(a) for a in adj_120]
print(f"  Edges: {n_edges} (expect 1200)")
print(f"  Degree range: {min(degrees)}-{max(degrees)} (expect 4-4)")
print(f"  Build time: {time.time() - t0:.1f}s")

N_VERTS = len(vertices_120)

if N_VERTS != 600:
    print(f"  WARNING: Got {N_VERTS} vertices instead of 600!")
    print(f"  This affects spectral calculations. Proceeding anyway...")

# Build Laplacian matrix
print("Building graph Laplacian...")
row, col, data = [], [], []
for i in range(N_VERTS):
    row.append(i)
    col.append(i)
    data.append(float(len(adj_120[i])))  # degree
    for j in adj_120[i]:
        row.append(i)
        col.append(j)
        data.append(-1.0)

L = sparse.csr_matrix((data, (row, col)), shape=(N_VERTS, N_VERTS))

# Get ALL eigenvalues (manageable size)
print("Computing eigenvalues of Laplacian...")
L_dense = L.toarray()
eigenvalues_all = np.linalg.eigvalsh(L_dense)
eigenvalues_all = np.sort(eigenvalues_all)

# Group eigenvalues by multiplicity
def group_eigenvalues(evals, tol=1e-6):
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i + count] - val) < tol:
            count += 1
        groups.append((val, count))
        i += count
    return groups

eigen_groups = group_eigenvalues(eigenvalues_all)
print(f"  Distinct eigenvalues: {len(eigen_groups)}")
print(f"  Eigenvalue spectrum (lambda, multiplicity):")
for lam, mult in eigen_groups[:20]:
    print(f"    lambda = {lam:10.6f}, mult = {mult}")
if len(eigen_groups) > 20:
    print(f"    ... ({len(eigen_groups) - 20} more)")

# ============================================================================
# KEY SPECTRAL QUANTITIES
# ============================================================================
print()
print("=" * 80)
print("KEY SPECTRAL QUANTITIES")
print("=" * 80)

def green_function_mp(m_sq, eigen_groups_list, N=None):
    """Compute lattice Green's function G(m^2) = (1/N) sum_k mult_k/(lambda_k + m^2)"""
    if N is None:
        N = sum(mult for _, mult in eigen_groups_list)
    m2 = mpf(str(m_sq))
    total = mpf(0)
    for lam, mult in eigen_groups_list:
        total += mpf(mult) / (mpf(str(lam)) + m2)
    return total / N

# Convert eigen_groups to high precision
eigen_groups_mp = [(float(lam), mult) for lam, mult in eigen_groups]

# The conformal mass from xi = 1/5
xi_current = mpf('1') / mpf('5')
R_sq = mpf(8)  # R = 2*sqrt(2), R^2 = 8 for unit edge 120-cell
R_ricci = mpf(6) / R_sq  # = 3/4
m_sq_current = xi_current * R_ricci  # = (1/5)(3/4) = 3/20

print(f"Current conformal coupling: xi = {nstr(xi_current, 10)}")
print(f"S^3 radius: R = 2*sqrt(2), R^2 = {nstr(R_sq, 10)}")
print(f"Ricci scalar: R_ricci = 6/R^2 = {nstr(R_ricci, 10)}")
print(f"Current mass^2: m^2 = xi * R_ricci = {nstr(m_sq_current, 10)}")

G1_lattice = green_function_mp(m_sq_current, eigen_groups_mp, N=N_VERTS)
G1_continuum = mpf(3) / (PI * phi**2)

print(f"G1 (lattice, m^2={nstr(m_sq_current,4)}) = {nstr(G1_lattice, 20)}")
print(f"G1 (continuum, 3/(pi*phi^2))            = {nstr(G1_continuum, 20)}")
print(f"Lattice - Continuum deficit              = {nstr(G1_lattice - G1_continuum, 10)}")

# Z_spectral via 120-cell: 1/alpha = 20*phi^4 * Z
Z_NIST = NIST_INV_ALPHA / (20 * phi**4)
print(f"\nZ_NIST = NIST / (20*phi^4) = {nstr(Z_NIST, 20)}")

# ============================================================================
# SECTION 1: SCALAR CURVATURE COUPLING
# ============================================================================
print()
print("=" * 80)
print("SECTION 1: SCALAR CURVATURE COUPLING -- ANGULAR DEFICIT OF DODECAHEDRON")
print("=" * 80)

dihedral_120cell = 4 * PI / 5
dihedral_dodecahedron = mpi - atan(mpf(2))  # arctan(2) complement

print(f"Dihedral angle of 120-cell (between cells): {nstr(dihedral_120cell * 180 / PI, 10)} degrees")
print(f"Dihedral angle of regular dodecahedron: {nstr(dihedral_dodecahedron * 180 / PI, 10)} degrees")

dihedral_inflated = 2 * PI / 3  # 120 deg -- exact value on S^3
angular_inflation = dihedral_inflated - dihedral_dodecahedron
angular_deficit_flat = 2 * PI - 3 * dihedral_dodecahedron
angular_deficit_curved = 2 * PI - 3 * dihedral_inflated

print(f"\nFlat dodecahedron dihedral: {nstr(dihedral_dodecahedron * 180 / PI, 10)} deg")
print(f"S^3 inflated dihedral: {nstr(dihedral_inflated * 180 / PI, 10)} deg")
print(f"Angular inflation per face: {nstr(angular_inflation * 180 / PI, 10)} deg")
print(f"Angular deficit at edge (flat):  {nstr(angular_deficit_flat * 180 / PI, 10)} deg")
print(f"Angular deficit at edge (S^3):   {nstr(angular_deficit_curved * 180 / PI, 10)} deg (should be 0)")

# Angular deficit per vertex of dodecahedron
vertex_angle = 3 * PI / 5  # interior angle of regular pentagon = 108 deg
angular_deficit_vertex = 2 * PI - 3 * vertex_angle  # pi/5

print(f"\nPentagonal interior angle: {nstr(vertex_angle * 180 / PI, 10)} deg")
print(f"Angular deficit per dodecahedral vertex: {nstr(angular_deficit_vertex * 180 / PI, 10)} deg = pi/5")
print(f"Total Gaussian curvature (Descartes) = 20 * pi/5 = 4*pi = {nstr(20 * angular_deficit_vertex, 10)}")

R_val = 2 * msqrt(2)
vol_S3 = 2 * PI**2 * R_val**3
integrated_curvature_S3 = R_ricci * vol_S3

print(f"\nS^3 radius R = {nstr(R_val, 10)}")
print(f"Volume of S^3 = 2*pi^2*R^3 = {nstr(vol_S3, 10)}")
print(f"Integrated curvature = R_ricci * Vol = {nstr(integrated_curvature_S3, 10)}")
print(f"Curvature per cell (120 cells) = {nstr(integrated_curvature_S3 / 120, 10)}")
print(f"Curvature per vertex (600 verts) = {nstr(integrated_curvature_S3 / 600, 10)}")

curvature_correction_ratio = angular_deficit_vertex / (2 * PI)  # = 1/10
print(f"\nAngular deficit / (2*pi) = pi/5 / (2*pi) = {nstr(curvature_correction_ratio, 10)} = 1/10")

# ============================================================================
# SECTION 2: VOLUME ELEMENT CORRECTION
# ============================================================================
print()
print("=" * 80)
print("SECTION 2: VOLUME ELEMENT CORRECTION (S^3 vs FLAT)")
print("=" * 80)

# Edge length a = 2/phi^2 = 3 - sqrt(5), R = 2*sqrt(2)
a_edge_val = 2 / phi**2  # actual edge length in these coordinates
a_over_R = a_edge_val / R_val  # edge/radius
vol_ratio_edge = msin(a_over_R)**2 / a_over_R**2

print(f"Edge length a = 2/phi^2 = {nstr(a_edge_val, 10)}")
print(f"a/R = {nstr(a_over_R, 10)}")
print(f"Volume element ratio at edge scale: sin^2(a/R)/(a/R)^2 = {nstr(vol_ratio_edge, 15)}")
print(f"Fractional correction: {nstr(1 - vol_ratio_edge, 10)} = {nstr((1 - vol_ratio_edge) * 100, 4)}%")

# Propagator ratio
nn_angle = a_over_R
G_ratio_nn = nn_angle / msin(nn_angle)

print(f"\nGreen's function ratio at nearest neighbor:")
print(f"G_S3/G_flat = (r/R)/sin(r/R) = {nstr(G_ratio_nn, 15)}")
print(f"Fractional excess: {nstr(G_ratio_nn - 1, 10)} = {nstr((G_ratio_nn - 1) * 1e6, 6)} ppm")

# ============================================================================
# ROAD 2 SETUP AND GAP ANALYSIS
# ============================================================================
print()
print("=" * 80)
print("GAP ANALYSIS")
print("=" * 80)

ROAD2_CURRENT = mpf('137.035998858')
gap_current = ROAD1 - ROAD2_CURRENT
fractional_gap = gap_current / NIST_INV_ALPHA

print(f"Road 1          = {nstr(ROAD1, 15)}")
print(f"Road 2 (xi=1/5) = {nstr(ROAD2_CURRENT, 15)}")
print(f"Gap             = {nstr(gap_current, 6)}")
print(f"Fractional gap  = {nstr(fractional_gap, 6)} = {nstr(fractional_gap * 1e9, 4)} ppb")

G1_deficit = G1_lattice - G1_continuum
deficit_1loop = mpf('-8.343e-7')
deficit_2loop = mpf('+8.362e-7')

print(f"\nG1 deficit (lattice - continuum) = {nstr(G1_deficit, 10)}")
print(f"1-loop deficit = {nstr(deficit_1loop, 6)}")
print(f"2-loop deficit = {nstr(deficit_2loop, 6)}")
print(f"Sum of deficits = {nstr(deficit_1loop + deficit_2loop, 6)}")

if abs(G1_deficit) > 1e-15:
    sensitivity = deficit_1loop * ROAD1_BASE / G1_deficit
    print(f"Sensitivity d(1/alpha)/dG1 ~ {nstr(sensitivity, 6)}")

# ============================================================================
# SECTION 3: EXACT CONFORMAL COUPLING FOR NIST
# ============================================================================
print()
print("=" * 80)
print("SECTION 3: FIND EXACT xi THAT GIVES NIST VALUE")
print("=" * 80)

t_NIST = (NIST_INV_ALPHA - ROAD2_CURRENT) / (ROAD1 - ROAD2_CURRENT)
print(f"Interpolation parameter t = (NIST - Road2)/(Road1 - Road2) = {nstr(t_NIST, 10)}")
print(f"  NIST is {nstr(t_NIST * 100, 4)}% of the way from Road 2 to Road 1")

G1_NIST = G1_lattice + t_NIST * (G1_continuum - G1_lattice)
print(f"\nG1 needed for NIST = {nstr(G1_NIST, 20)}")
print(f"G1_lattice (current) = {nstr(G1_lattice, 20)}")
print(f"G1_continuum = {nstr(G1_continuum, 20)}")

# Binary search for xi_exact
def G1_of_xi(xi_val):
    m2 = xi_val * R_ricci
    return green_function_mp(float(m2), eigen_groups_mp, N=N_VERTS)

print("\nBinary search for xi_exact...")
xi_lo = mpf('0.001')
xi_hi = mpf('2.0')
for _ in range(200):
    xi_mid = (xi_lo + xi_hi) / 2
    G1_mid = G1_of_xi(xi_mid)
    if G1_mid > G1_NIST:
        xi_lo = xi_mid
    else:
        xi_hi = xi_mid

xi_exact = (xi_lo + xi_hi) / 2
G1_at_exact = G1_of_xi(xi_exact)
m2_exact = xi_exact * R_ricci

print(f"xi_exact = {nstr(xi_exact, 30)}")
print(f"m^2_exact = xi_exact * 3/4 = {nstr(m2_exact, 30)}")
print(f"G1(xi_exact) = {nstr(G1_at_exact, 20)}")
print(f"G1_NIST target = {nstr(G1_NIST, 20)}")

print(f"\nxi_exact analysis:")
print(f"  xi = {nstr(xi_exact, 20)}")
print(f"  1/xi = {nstr(1/xi_exact, 20)}")
print(f"  xi * 5 = {nstr(xi_exact * 5, 20)}")
print(f"  xi * 6 = {nstr(xi_exact * 6, 20)}")
print(f"  xi * 8 = {nstr(xi_exact * 8, 20)}")
print(f"  xi * 10 = {nstr(xi_exact * 10, 20)}")
print(f"  xi * 20 = {nstr(xi_exact * 20, 20)}")
print(f"  xi * phi = {nstr(xi_exact * phi, 20)}")
print(f"  xi * phi^2 = {nstr(xi_exact * phi**2, 20)}")
print(f"  xi * pi = {nstr(xi_exact * PI, 20)}")
print(f"  xi - 1/5 = {nstr(xi_exact - mpf('0.2'), 15)}")
print(f"  xi - 1/6 = {nstr(xi_exact - mpf(1)/6, 15)}")
print(f"  xi - 1/8 = {nstr(xi_exact - mpf('0.125'), 15)}")

# Standard conformal couplings
for d_eff in range(2, 20):
    xi_d = mpf(d_eff - 2) / (4 * (d_eff - 1))
    if abs(xi_d - xi_exact) / max(abs(xi_exact), mpf('1e-20')) < 0.05:
        print(f"  CLOSE to d={d_eff} conformal: xi = (d-2)/(4(d-1)) = {nstr(xi_d, 15)}, diff = {nstr(xi_d - xi_exact, 10)}")

# Check phi powers
print(f"\n  Checking xi_exact against phi powers:")
for n in range(-5, 6):
    ratio = xi_exact / phi**n
    # Check if this is a simple fraction
    for q in range(1, 50):
        p = int(float(ratio * q + mpf('0.5')))
        if p > 0 and abs(ratio - mpf(p)/q) < mpf('1e-5'):
            residual_ppb = float(abs(ratio - mpf(p)/q) / xi_exact * 1e9)
            if residual_ppb < 100:
                print(f"    xi / phi^{n:+d} ~ {p}/{q} (residual: {residual_ppb:.1f} ppb)")

# ============================================================================
# SECTION 4: HEAT KERNEL ON S^3
# ============================================================================
print()
print("=" * 80)
print("SECTION 4: HEAT KERNEL ON S^3 -- EXACT GREEN'S FUNCTION")
print("=" * 80)

print("Continuum S^3 spectrum:")
print(f"  Eigenvalues: lambda_l = l(l+2)/R^2, l = 0, 1, 2, ...")
print(f"  Multiplicity: (l+1)^2")
print(f"  R^2 = {nstr(R_sq, 10)}")

def G_S3_continuum(m_sq, R_squared, l_max=5000):
    """Green's function on S^3 using continuum eigenvalues."""
    m2 = mpf(str(m_sq))
    R2 = mpf(str(R_squared))
    Vol = 2 * PI**2 * R2**mpf('1.5')
    total = mpf(0)
    for l in range(l_max + 1):
        lam_l = mpf(l) * (mpf(l) + 2) / R2
        mult = (mpf(l) + 1)**2
        total += mult / (lam_l + m2)
    return total / Vol

print(f"\nComputing continuum S^3 Green's function (l_max=5000)...")
G_S3_exact = G_S3_continuum(float(m_sq_current), float(R_sq), l_max=5000)
print(f"G_S3(m^2=3/20, R^2=8) = {nstr(G_S3_exact, 20)}")
print(f"G1_lattice (120-cell)  = {nstr(G1_lattice, 20)}")
print(f"G1_continuum (flat)    = {nstr(G1_continuum, 20)}")
print(f"G_S3 - G1_lattice      = {nstr(G_S3_exact - G1_lattice, 10)}")
print(f"G_S3 - G1_continuum    = {nstr(G_S3_exact - G1_continuum, 10)}")

# Convergence check
G_S3_1000 = G_S3_continuum(float(m_sq_current), float(R_sq), l_max=1000)
G_S3_2000 = G_S3_continuum(float(m_sq_current), float(R_sq), l_max=2000)
print(f"\nConvergence check:")
print(f"  l_max=1000: {nstr(G_S3_1000, 20)}")
print(f"  l_max=2000: {nstr(G_S3_2000, 20)}")
print(f"  l_max=5000: {nstr(G_S3_exact, 20)}")

# ============================================================================
# SECTION 5: WINDING NUMBER CORRECTIONS
# ============================================================================
print()
print("=" * 80)
print("SECTION 5: WINDING NUMBER / TOPOLOGY CORRECTIONS")
print("=" * 80)

m_val = msqrt(m_sq_current)
circumference = 2 * PI * R_val

print(f"Mass: m = sqrt(3/20) = {nstr(m_val, 10)}")
print(f"S^3 circumference: 2*pi*R = {nstr(circumference, 10)}")
print(f"m * circumference = {nstr(m_val * circumference, 10)}")

print(f"\nWinding number contributions to G(0, m^2):")
total_winding = mpf(0)
for n in range(1, 11):
    dist_n = circumference * n
    G_n = mexp(-m_val * dist_n) / (4 * PI * dist_n)
    total_winding += 2 * G_n
    print(f"  n={n}: distance = {nstr(dist_n, 6)}, G_n = {nstr(G_n, 10)}, cumulative = {nstr(total_winding, 10)}")

print(f"\nTotal winding correction = {nstr(total_winding, 10)}")
print(f"Ratio to G1_lattice = {nstr(total_winding / G1_lattice, 10)}")
print(f"In ppb: {nstr(total_winding / G1_lattice * 1e9, 4)}")

graph_diameter = 15
half_circum = PI * R_val
print(f"\nGraph diameter: {graph_diameter} edges")
print(f"S^3 half-circumference: pi*R = {nstr(half_circum, 6)}")
print(f"Graph reaches {nstr(mpf(graph_diameter)/circumference, 4)} full circumferences")

# ============================================================================
# SECTION 6: REGGE CALCULUS
# ============================================================================
print()
print("=" * 80)
print("SECTION 6: REGGE CALCULUS -- DEFICIT ANGLES IN THE 120-CELL")
print("=" * 80)

flat_deficit_per_edge = 2 * PI - 3 * dihedral_dodecahedron
print(f"Dihedral of flat dodecahedron: {nstr(dihedral_dodecahedron * 180 / PI, 10)} deg")
print(f"Deficit per edge (flat): {nstr(flat_deficit_per_edge * 180 / PI, 10)} deg")
print(f"Deficit per edge (S^3): 0 deg (exact tessellation)")

a_edge = mpf(1)
A_pentagon = a_edge**2 / 4 * msqrt(25 + 10 * msqrt(5))
total_Regge_action_flat = 1200 * a_edge * flat_deficit_per_edge

print(f"\nPentagonal face area: {nstr(A_pentagon, 10)}")
print(f"Total Regge action (edge hinges): 1200 * a * delta = {nstr(total_Regge_action_flat, 10)}")
print(f"Regge action per vertex: {nstr(total_Regge_action_flat / 600, 10)}")
print(f"Regge action per cell: {nstr(total_Regge_action_flat / 120, 10)}")

curvature_per_edge = flat_deficit_per_edge / (2 * PI)
print(f"\nCurvature fraction per edge: delta/(2*pi) = {nstr(curvature_per_edge, 10)}")

# ============================================================================
# SECTION 7: DECOMPOSE ALGEBRAIC CORRECTION
# ============================================================================
print()
print("=" * 80)
print("SECTION 7: DECOMPOSE ALGEBRAIC CORRECTION")
print("=" * 80)

C_exact = NIST_INV_ALPHA / ROAD1_BASE - 1
C_original = 1 / (2 * phi**27)

print(f"Original correction: C = 1/(2*phi^27) = {nstr(C_original, 20)}")
print(f"Exact correction for NIST: C_exact     = {nstr(C_exact, 20)}")
print(f"Residual: C_exact - C_original         = {nstr(C_exact - C_original, 15)}")
print(f"Fractional residual: {nstr((C_exact - C_original)/C_original, 10)}")

residual = C_exact - C_original
print(f"\nResidual = {nstr(residual, 15)}")
print(f"Residual / C_original = {nstr(residual / C_original, 15)}")
print(f"Residual / C_original^2 = {nstr(residual / C_original**2, 15)}")

k_coeff = -residual / C_original**2
print(f"If C_exact = C_orig - k*C_orig^2, then k = {nstr(k_coeff, 15)}")

# Geometric series
C_geom = C_original / (1 - C_original)
inv_alpha_geom = ROAD1_BASE * (1 + C_geom)
ppb_geom = float((inv_alpha_geom - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"\nGeometric series: C = 1/(2*phi^27 - 1)")
print(f"  1/alpha_geom = {nstr(inv_alpha_geom, 15)} ({ppb_geom:.4f} ppb)")

# Check: 1/(2*phi^27 + k) for various k
print(f"\n1/C_exact = {nstr(1/C_exact, 20)}")
print(f"2*phi^27  = {nstr(2*phi**27, 20)}")
curv_addition = 1/C_exact - 2*phi**27
print(f"Difference: 1/C_exact - 2*phi^27 = {nstr(curv_addition, 15)}")

# ============================================================================
# SECTION 8: NUMERICAL GAP DECOMPOSITION
# ============================================================================
print()
print("=" * 80)
print("SECTION 8: NUMERICAL GAP DECOMPOSITION")
print("=" * 80)

gap_road1_NIST = ROAD1 - NIST_INV_ALPHA
gap_NIST_road2 = NIST_INV_ALPHA - ROAD2_CURRENT
gap_road1_road2 = ROAD1 - ROAD2_CURRENT

print(f"Road 1 - NIST   = {nstr(gap_road1_NIST, 10)} ({nstr(gap_road1_NIST/NIST_INV_ALPHA * 1e9, 4)} ppb)")
print(f"NIST - Road 2   = {nstr(gap_NIST_road2, 10)} ({nstr(gap_NIST_road2/NIST_INV_ALPHA * 1e9, 4)} ppb)")
print(f"Road 1 - Road 2 = {nstr(gap_road1_road2, 10)} ({nstr(gap_road1_road2/NIST_INV_ALPHA * 1e9, 4)} ppb)")

# Decompose Road 1 overshoot
ov = gap_road1_NIST / NIST_INV_ALPHA
print(f"\n--- Road 1 overshoot = {nstr(ov, 10)} ---")
print(f"  vs C_original^2 = {nstr(C_original**2, 10)}")
print(f"  ov / C_original^2 = {nstr(ov / C_original**2, 10)}")
print(f"  ov * 20 * phi^4 = {nstr(ov * 20 * phi**4, 10)}")

# Decompose Road 2 undershoot
un = gap_NIST_road2 / NIST_INV_ALPHA
print(f"\n--- Road 2 undershoot = {nstr(un, 10)} ---")
print(f"  un * phi^27 = {nstr(un * phi**27, 10)}")
print(f"  un * phi^27 * 2 = {nstr(un * phi**27 * 2, 10)}")
print(f"  un / C_original = {nstr(un / C_original, 6)}")

# Try products of two known quantities
target_gap = gap_road1_road2 / NIST_INV_ALPHA
print(f"\n--- Scanning products for total gap (fractional) = {nstr(target_gap, 10)} ---")

known = {
    'phi': phi, '1/phi': 1/phi, 'phi^2': phi**2, 'phi^3': phi**3,
    'pi': PI, 'pi^2': PI**2, '1/pi': 1/PI,
    'G1_cont': G1_continuum, 'G1_lat': G1_lattice,
    '1/120': mpf(1)/120, '1/600': mpf(1)/600, '1/720': mpf(1)/720,
    '1/1200': mpf(1)/1200,
    'R_ricci': R_ricci, '1/R^2': 1/R_sq,
    'delta_v': angular_deficit_vertex,
    'vol_corr': 1 - vol_ratio_edge,
    'C_orig': C_original,
}

close_matches = []
for name1, v1 in known.items():
    for name2, v2 in known.items():
        prod = v1 * v2
        if abs(prod) > 1e-15 and abs(prod) < 1:
            ratio = target_gap / prod
            if abs(ratio) > 0.5 and abs(ratio) < 20:
                for q in range(1, 13):
                    p = int(float(abs(ratio) * q + mpf('0.5')))
                    if p > 0 and abs(abs(ratio) - mpf(p)/q) < abs(ratio) * 0.01:
                        ppb_err = float(abs(target_gap - prod * mpf(p)/q * (1 if ratio > 0 else -1)) / abs(target_gap) * 1e3)
                        if ppb_err < 100:
                            sign = "+" if ratio > 0 else "-"
                            close_matches.append((ppb_err, f"gap ~ {sign}{p}/{q} * {name1} * {name2} (err: {ppb_err:.1f} per mille)"))

close_matches.sort()
for err, desc in close_matches[:20]:
    print(f"  {desc}")

# ============================================================================
# SECTION 9: DIMENSIONAL ANALYSIS
# ============================================================================
print()
print("=" * 80)
print("SECTION 9: DIMENSIONAL ANALYSIS")
print("=" * 80)

print("Dimension-dependent quantities:")
for d in [3, 4, 5, 6]:
    xi_d = mpf(d - 2) / (4 * (d - 1))
    R_d = mpf(d) * (d - 1) / R_sq
    m2_d = xi_d * R_d
    vol_d = 2 * PI**((d + 1) / mpf(2)) / mgamma((d + 1) / mpf(2)) * R_val**d
    G1_d = green_function_mp(float(m2_d), eigen_groups_mp, N=N_VERTS)

    print(f"\n  d = {d}:")
    print(f"    xi = (d-2)/(4(d-1)) = {nstr(xi_d, 10)}")
    print(f"    R_scalar = d(d-1)/R^2 = {nstr(R_d, 10)}")
    print(f"    m^2 = xi * R_scalar = {nstr(m2_d, 10)}")
    print(f"    G1(m^2) = {nstr(G1_d, 15)}")

# ============================================================================
# SECTION 10: COMPREHENSIVE xi SCAN
# ============================================================================
print()
print("=" * 80)
print("SECTION 10: COMPREHENSIVE xi SCAN FOR NIST MATCH")
print("=" * 80)

# We need the linear model relating G1 to 1/alpha
f_road1 = C_original  # 1/(2*phi^27)
f_road2 = ROAD2_CURRENT / ROAD1_BASE - 1

if abs(G1_continuum - G1_lattice) > 1e-20:
    b_coeff = (f_road1 - f_road2) / (G1_continuum - G1_lattice)
    a_coeff = f_road1 - b_coeff * G1_continuum

    print(f"Linear model: f(G1) = a + b*G1")
    print(f"  a = {nstr(a_coeff, 15)}")
    print(f"  b = {nstr(b_coeff, 15)}")

    f_NIST = NIST_INV_ALPHA / ROAD1_BASE - 1
    G1_for_NIST = (f_NIST - a_coeff) / b_coeff
    print(f"  G1 needed for NIST: {nstr(G1_for_NIST, 15)}")

    print(f"\nxi scan (m^2 = xi * 3/4):")
    print(f"{'xi':>15s} {'m^2':>10s} {'G1':>18s} {'1/alpha':>18s} {'ppb':>10s} {'note':>20s}")
    print("-" * 95)

    best_ppb = 1e10
    best_xi = mpf(0)
    best_inv_alpha = mpf(0)

    xi_candidates = [
        (mpf(1)/mpf(8), "1/8 (d=3)"),
        (mpf(1)/mpf(6), "1/6 (d=4)"),
        (mpf(3)/mpf(16), "3/16 (d=5)"),
        (mpf(1)/mpf(5), "1/5 (current)"),
        (mpf(5)/mpf(24), "5/24 (d=7)"),
        (mpf(1)/mpf(4), "1/4"),
        (mpf(3)/mpf(10), "3/10"),
        (mpf(1)/mpf(3), "1/3"),
        (1 / (2 * phi**2), "1/(2*phi^2)"),
        (mpf(2)/mpf(9), "2/9"),
        (mpf(3)/mpf(14), "3/14"),
        (mpf(7)/mpf(36), "7/36"),
        (mpf(4)/mpf(21), "4/21"),
        (mpf(5)/mpf(26), "5/26"),
        (mpf(11)/mpf(60), "11/60"),
        (xi_exact, "xi_exact (fit)"),
    ]

    for xi_val, note in xi_candidates:
        if xi_val <= 0:
            continue
        m2_val = xi_val * R_ricci
        G1_val = green_function_mp(float(m2_val), eigen_groups_mp, N=N_VERTS)
        f_val = a_coeff + b_coeff * G1_val
        inv_alpha_val = ROAD1_BASE * (1 + f_val)
        ppb = float((inv_alpha_val - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)

        if abs(ppb) < abs(best_ppb):
            best_ppb = ppb
            best_xi = xi_val
            best_inv_alpha = inv_alpha_val

        print(f"{nstr(xi_val, 10):>15s} {nstr(m2_val, 6):>10s} {nstr(G1_val, 14):>18s} {nstr(inv_alpha_val, 12):>18s} {ppb:>10.4f} {note:>20s}")

    print(f"\nBest: xi = {nstr(best_xi, 20)}, 1/alpha = {nstr(best_inv_alpha, 15)}, {best_ppb:.4f} ppb")

# ============================================================================
# SECTION 11: COMBINED CORRECTIONS TO ROAD 1
# ============================================================================
print()
print("=" * 80)
print("SECTION 11: COMBINED CORRECTIONS TO ROAD 1")
print("=" * 80)

target = NIST_INV_ALPHA
overshoot = ROAD1 - target

print(f"Road 1 overshoots NIST by: {nstr(overshoot, 10)}")
corr_needed = 1 - NIST_INV_ALPHA / ROAD1
print(f"Multiplicative correction needed: 1 - {nstr(corr_needed, 15)}")

corrections = {
    "R_ricci/(8*pi^2)": R_ricci / (8 * PI**2),
    "1/(4*phi^54)": 1 / (4 * phi**54),
    "pi/(2*phi^54)": PI / (2 * phi**54),
    "C_orig^2": C_original**2,
    "3/(20*phi^54)": 3 / (20 * phi**54),
    "1/(phi^27*pi)": 1 / (phi**27 * PI),
    "deficit^2/(4*pi^2)": (angular_deficit_vertex)**2 / (4 * PI**2),
    "1/(120*phi^27)": 1 / (120 * phi**27),
    "1/(600*phi^27)": 1 / (600 * phi**27),
    "1/(40*phi^31)": 1 / (40 * phi**31),
}

print(f"\n{'Formula':>30s} {'Value':>18s} {'Needed':>18s} {'Ratio':>12s} {'ppb res':>10s}")
print("-" * 92)

for name, val in corrections.items():
    ratio = corr_needed / val if abs(val) > 1e-50 else mpf(0)
    inv_alpha_corrected = ROAD1 * (1 - val)
    ppb_residual = float((inv_alpha_corrected - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
    print(f"{name:>30s} {nstr(val, 12):>18s} {nstr(corr_needed, 12):>18s} {nstr(ratio, 8):>12s} {ppb_residual:>10.4f}")

# ============================================================================
# SECTION 12: LATTICE vs S^3 EIGENVALUE COMPARISON
# ============================================================================
print()
print("=" * 80)
print("SECTION 12: LATTICE vs S^3 CONTINUUM EIGENVALUE COMPARISON")
print("=" * 80)

print("S^3 continuum eigenvalues (l(l+2)/8) vs lattice:")
print(f"{'l':>4s} {'lambda_S3':>12s} {'mult':>6s} {'cumul':>8s} {'near_lat':>12s} {'lat_mult':>10s}")
print("-" * 60)

cum_modes = 0
for l in range(16):
    lam_S3 = l * (l + 2) / 8.0
    mult_S3 = (l + 1) ** 2
    cum_modes += mult_S3

    best_match = None
    best_dist = 999
    for lam_lat, mult_lat in eigen_groups:
        dist = abs(lam_lat - lam_S3)
        if dist < best_dist:
            best_dist = dist
            best_match = (lam_lat, mult_lat)

    print(f"{l:>4d} {lam_S3:>12.6f} {mult_S3:>6d} {cum_modes:>8d} {best_match[0]:>12.6f} {best_match[1]:>10d}")

# ============================================================================
# SECTION 13: SYSTEMATIC CORRECTION SCAN
# ============================================================================
print()
print("=" * 80)
print("SECTION 13: SYSTEMATIC CORRECTION COMBINATIONS")
print("=" * 80)

delta_needed = NIST_INV_ALPHA - ROAD1  # negative
small = delta_needed / ROAD1_BASE
print(f"Additive correction needed: delta = {nstr(delta_needed, 10)}")
print(f"Small number (delta/base) = {nstr(small, 15)}")
print(f"  * phi^27 = {nstr(small * phi**27, 15)}")
print(f"  * phi^54 = {nstr(small * phi**54, 15)}")
print(f"  / C_original = {nstr(small / C_original, 15)}")
print(f"  / C_original^2 = {nstr(small / C_original**2, 15)}")

# The extra correction
extra = C_exact - C_original
print(f"\nExtra correction (C_exact - C_original) = {nstr(extra, 15)}")

# Scan: extra = a * phi^b * pi^c / d
print("\n--- Scanning a*phi^b*pi^c/d for the extra correction ---")
matches = []
for a_num in range(-10, 11):
    if a_num == 0:
        continue
    for b_pow in range(-60, 5):
        for c_pow in range(-3, 4):
            for d_denom in [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 40, 60, 120, 600, 1200, 720]:
                candidate = mpf(a_num) * phi**b_pow * PI**c_pow / d_denom
                if abs(candidate) > 1e-20 and abs(candidate) < 1e-3:
                    ratio = extra / candidate
                    if abs(ratio - 1) < 1e-3:
                        ppb_err = float(abs(ratio - 1) * 1e9)
                        if ppb_err < 50:
                            matches.append((ppb_err, f"extra ~ {a_num}*phi^{b_pow}*pi^{c_pow}/{d_denom} (err: {ppb_err:.2f} ppb)"))

matches.sort()
print(f"Found {len(matches)} close matches:")
for ppb, desc in matches[:30]:
    print(f"  {desc}")

# ============================================================================
# SECTION 14: EXACT CLOSED FORM SEARCH
# ============================================================================
print()
print("=" * 80)
print("SECTION 14: EXACT CLOSED FORM FOR C_exact")
print("=" * 80)

K_factor = 1 - C_exact / C_original
print(f"If C_exact = C_original * (1 - K):")
print(f"  K = {nstr(K_factor, 15)}")
print(f"  K * phi^27 = {nstr(K_factor * phi**27, 15)}")

epsilon_needed = 1 - C_exact * 2 * phi**27
print(f"\nIf C_exact = 1/(2*phi^27) * (1 - eps):")
print(f"  eps = {nstr(epsilon_needed, 15)}")
print(f"  eps * phi^27 = {nstr(epsilon_needed * phi**27, 15)}")
print(f"  eps * 600 = {nstr(epsilon_needed * 600, 15)}")
print(f"  eps * 1200 = {nstr(epsilon_needed * 1200, 15)}")
print(f"  eps * 120 = {nstr(epsilon_needed * 120, 15)}")

# SECTION: Search 1/(2*phi^27 + k) for integer k
print(f"\nSearching 1/(2*phi^27 + k) for integer k:")
for k in range(-100, 101):
    C_trial = 1 / (2 * phi**27 + k)
    trial = ROAD1_BASE * (1 + C_trial)
    ppb_trial = float((trial - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
    if abs(ppb_trial) < 1:
        print(f"  k = {k}: 1/alpha = {nstr(trial, 15)} ({ppb_trial:.4f} ppb)")

# Fractional k
print(f"\nSearching 1/(2*phi^27 + p/q):")
for q in range(1, 31):
    for p in range(-300, 301):
        k_frac = mpf(p) / q
        if abs(k_frac) > 300:
            continue
        C_trial = 1 / (2 * phi**27 + k_frac)
        trial = ROAD1_BASE * (1 + C_trial)
        ppb_trial = float((trial - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
        if abs(ppb_trial) < 0.02:
            print(f"  p/q = {p}/{q}: 1/(2*phi^27 + {p}/{q}) -> {nstr(trial, 15)} ({ppb_trial:.6f} ppb)")

# ============================================================================
# SECTION 15: TARGETED SCAN WITH CORRECTION = n/d * phi^p * pi^q
# ============================================================================
print()
print("=" * 80)
print("SECTION 15: TARGETED CORRECTION TO 1/(2*phi^27)")
print("=" * 80)

# Search for: C = 1/(2*phi^27 + correction)
# where correction is a simple expression
print(f"\nDenominator offset needed: {nstr(curv_addition, 15)}")

# Fine scan
print("\nSearching: 1/(2*phi^27 + a*phi^n*pi^m/d) for NIST match:")
best_total = []

for a_num in range(-5, 6):
    if a_num == 0:
        continue
    for n_pow in range(-10, 28):
        for m_pow in range(-2, 3):
            for d_denom in [1, 2, 3, 4, 5, 6, 8, 10, 12, 20, 24, 40, 60, 120]:
                offset = mpf(a_num) * phi**n_pow * PI**m_pow / d_denom
                denom = 2 * phi**27 + offset
                if abs(denom) < 1:  # avoid division by zero or near-zero
                    continue
                C_trial = 1 / denom
                trial = ROAD1_BASE * (1 + C_trial)
                ppb_trial = float((trial - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
                if abs(ppb_trial) < 0.15:
                    best_total.append((abs(ppb_trial), ppb_trial,
                        f"1/(2*phi^27 + {a_num}*phi^{n_pow}*pi^{m_pow}/{d_denom})", nstr(trial, 15)))

best_total.sort()
print(f"Found {len(best_total)} candidates within 0.15 ppb:")
for absp, ppb, formula, val in best_total[:40]:
    print(f"  {formula:50s} -> {val} ({ppb:+.6f} ppb)")

# ============================================================================
# SECTION 15b: DEEP ANALYSIS OF THE NUMBER 188
# ============================================================================
print()
print("=" * 80)
print("SECTION 15b: DEEP ANALYSIS OF DENOMINATOR OFFSET 188.126...")
print("=" * 80)

# The NIST-exact formula is: 1/alpha = base * (1 + 1/(2*phi^27 + D))
# where D = 188.126... (the denominator offset)
D_exact = curv_addition  # 1/C_exact - 2*phi^27
print(f"D_exact = {nstr(D_exact, 20)}")
print(f"D_exact = {nstr(D_exact, 6)} (rounded)")

# What is 188.126 in terms of geometric quantities of the 120-cell?
print(f"\n--- Decomposing D = {nstr(D_exact, 8)} ---")

# Key numbers of the 120-cell:
# V=600, E=1200, F=720, C=120
# Euler char = V - E + F - C = 0 (for 4D polytope boundary = S^3)
# Vertex degree = 4, face sides = 5, cells per edge = 3
# Dihedral angle on S^3 = 120 deg
# Angular deficit (flat) per edge = 10.305 deg = pi/5 * (3/pi) * ... etc.

print(f"  D / phi = {nstr(D_exact / phi, 10)}")
print(f"  D / phi^2 = {nstr(D_exact / phi**2, 10)}")
print(f"  D / phi^3 = {nstr(D_exact / phi**3, 10)}")
print(f"  D / phi^4 = {nstr(D_exact / phi**4, 10)}")
print(f"  D / phi^5 = {nstr(D_exact / phi**5, 10)}")
print(f"  D / pi = {nstr(D_exact / PI, 10)}")
print(f"  D / pi^2 = {nstr(D_exact / PI**2, 10)}")
print(f"  D * pi = {nstr(D_exact * PI, 10)}")
print(f"  D / (4*pi) = {nstr(D_exact / (4*PI), 10)}")
print(f"  D / (2*pi) = {nstr(D_exact / (2*PI), 10)}")
print(f"  D / 120 = {nstr(D_exact / 120, 10)}")
print(f"  D / 60 = {nstr(D_exact / 60, 10)}")
print(f"  D / 600 = {nstr(D_exact / 600, 10)}")
print(f"  D / 1200 = {nstr(D_exact / 1200, 10)}")
print(f"  D / 720 = {nstr(D_exact / 720, 10)}")

# Check: is D close to any product of small integers and phi/pi?
print(f"\n--- Checking D against simple expressions ---")
candidates_188 = []
for a in range(1, 30):
    for b in range(1, 30):
        for phi_pow in range(-5, 8):
            for pi_pow in range(-3, 4):
                val = mpf(a) / mpf(b) * phi**phi_pow * PI**pi_pow
                if abs(val - D_exact) / D_exact < 0.001:
                    err_ppb = float(abs(val - D_exact) / D_exact * 1e6)
                    if err_ppb < 100:
                        candidates_188.append((err_ppb, f"{a}/{b} * phi^{phi_pow} * pi^{pi_pow} = {nstr(val, 10)}"))

candidates_188.sort()
for err, desc in candidates_188[:25]:
    print(f"  D ~ {desc} (err: {err:.1f} ppm)")

# Check combinations with 120-cell numbers
print(f"\n--- D vs 120-cell number combinations ---")
cell_nums = {
    'V': 600, 'E': 1200, 'F': 720, 'C': 120,
    'V-E+F': 120, 'E/V': 2, 'F/C': 6, 'V/C': 5,
    'E/C': 10, 'E-V': 600, 'V-F': -120, 'F-C': 600,
    'deg': 4, 'face_sides': 5, 'cells_per_edge': 3,
}
for name, n in cell_nums.items():
    if n != 0:
        ratio = D_exact / n
        print(f"  D / {name}({n}) = {nstr(ratio, 10)}")
        # Is this ratio close to a phi power times pi power?
        for pp in range(-4, 6):
            for qp in range(-2, 3):
                val = phi**pp * PI**qp
                if abs(ratio / val - 1) < 0.01:
                    r2 = ratio / val
                    print(f"    = phi^{pp} * pi^{qp} * {nstr(r2, 8)}")

# The REALLY interesting check: 188 = 4*47, or 2^2 * 47
# But also: V - E + F - C = 600 - 1200 + 720 - 120 = 0 (Euler)
# V + E + F + C = 600 + 1200 + 720 + 120 = 2640
# 188 could be (V*something + E*something + ...)

# More direct: is D = a + b*phi for integers a,b?
# Since D = 188.126..., check:
# phi = 1.618..., so D/phi = 116.28, 188/phi = 116.17
# D - 188 = 0.126... Is that phi/5/phi^2/something?
D_frac = D_exact - 188
print(f"\n  D - 188 = {nstr(D_frac, 15)}")
print(f"  D - 188 in terms of phi: {nstr(D_frac / phi, 10)}")
print(f"  D - 188 in terms of 1/phi: {nstr(D_frac * phi, 10)}")
print(f"  D - 188 in terms of pi: {nstr(D_frac / PI, 10)}")
print(f"  D - 188 in terms of 1/pi: {nstr(D_frac * PI, 10)}")
print(f"  (D-188)*8 = {nstr(D_frac * 8, 10)}")
print(f"  (D-188)*10 = {nstr(D_frac * 10, 10)}")
print(f"  (D-188)*120 = {nstr(D_frac * 120, 10)}")
print(f"  (D-188)*600 = {nstr(D_frac * 600, 10)}")

# Is D perhaps the integral of the Regge curvature over some subset?
# Regge action per edge = a * delta_flat = 0.764 * 0.1799 = 0.1375
# 188 / 0.1375 ~ 1367 ~ 1200 + 167
print(f"\n  total_Regge / (2*pi) = {nstr(total_Regge_action_flat / (2*PI), 10)}")
print(f"  total_Regge / pi = {nstr(total_Regge_action_flat / PI, 10)}")
print(f"  D / (total_Regge/(2*pi)) = {nstr(D_exact / (total_Regge_action_flat / (2*PI)), 10)}")

# Very important: D = 188.126 ~ 30*2*pi = 188.496 (close!)
print(f"\n  30 * 2*pi = {nstr(30 * 2 * PI, 10)}")
print(f"  D / (30*2*pi) = {nstr(D_exact / (30 * 2 * PI), 10)}")
# Not exact. How about:
print(f"  60*pi = {nstr(60*PI, 10)}")
print(f"  D / (60*pi) = {nstr(D_exact / (60*PI), 10)}")

# 188.126 ~ 188 + 1/8 = 188.125 (VERY close!)
print(f"\n  188 + 1/8 = {nstr(mpf(188) + mpf(1)/8, 15)}")
print(f"  D - (188 + 1/8) = {nstr(D_exact - 188 - mpf(1)/8, 15)}")
inv_alpha_188_125 = ROAD1_BASE * (1 + 1/(2*phi**27 + mpf(1505)/8))  # 188+1/8 = 1505/8
ppb_188_125 = float((inv_alpha_188_125 - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"  1/(2*phi^27 + 1505/8) -> {nstr(inv_alpha_188_125, 15)} ({ppb_188_125:.6f} ppb)")

# Or 188 + 1/phi^8?
print(f"  188 + 1/phi^8 = {nstr(188 + 1/phi**8, 15)}")
ppb_test = float((ROAD1_BASE * (1 + 1/(2*phi**27 + 188 + 1/phi**8)) - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"  -> ppb = {ppb_test:.6f}")

# 188 + phi/phi^2 = 188 + 1/phi = 188.618... no too big
# 188 + 1/(2*phi^2) = 188.3820... too big
# 188 + 1/8 = 188.125 vs D = 188.126...
# 188 + 1/phi^7 = 188 + 0.03444 = 188.034 no
# How about 188 + (3-sqrt(5))/6? That's 188 + 0.1273 = 188.127 VERY close
print(f"  188 + (3-sqrt(5))/6 = {nstr(188 + (3 - msqrt(5))/6, 15)}")
print(f"  = 188 + edge/6 = 188 + a/(6) where a = 2/phi^2")
ppb_edge6 = float((ROAD1_BASE * (1 + 1/(2*phi**27 + 188 + 2/(6*phi**2))) - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"  -> ppb = {ppb_edge6:.6f}")

# 188 + 2/(6*phi^2) = 188 + 1/(3*phi^2) = 188 + 1/(3*2.618) = 188 + 0.12732
print(f"  188 + 1/(3*phi^2) = {nstr(188 + 1/(3*phi**2), 15)}")
ppb_3phi2 = float((ROAD1_BASE * (1 + 1/(2*phi**27 + 188 + 1/(3*phi**2))) - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"  -> ppb = {ppb_3phi2:.6f}")

# Wait: 188 + 1/(2*pi*phi^2)?
# 1/(2*pi*phi^2) = 1/(2*3.14159*2.618) = 1/16.45 = 0.0608 no

# Let me try brute force: D = p + q/r * phi^s * pi^t
print(f"\n--- Fine-grained search near D = 188.126 ---")
best_D = []
for p_int in [187, 188, 189]:
    for q_num in range(-20, 21):
        if q_num == 0:
            continue
        for r_den in [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30]:
            for s_pow in range(-8, 8):
                for t_pow in range(-2, 3):
                    val = p_int + mpf(q_num) / r_den * phi**s_pow * PI**t_pow
                    if abs(val - D_exact) < 0.0001:
                        ppb_err = float(abs(val - D_exact) / D_exact * 1e6)
                        trial = ROAD1_BASE * (1 + 1/(2*phi**27 + val))
                        ppb_nist = float((trial - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
                        best_D.append((abs(ppb_nist), ppb_nist,
                            f"D = {p_int} + {q_num}/{r_den}*phi^{s_pow}*pi^{t_pow}", nstr(val, 12)))

best_D.sort()
print(f"Found {len(best_D)} candidates within 0.0001 of D:")
for absp, ppb_n, formula, val_str in best_D[:30]:
    print(f"  {formula:45s} = {val_str} ({ppb_n:+.6f} ppb from NIST)")

# Also try: D = n * phi^k for integer n, small k
print(f"\n--- D as n * phi^k ---")
for k in range(-5, 6):
    n_approx = D_exact / phi**k
    n_int = int(float(n_approx + mpf('0.5')))
    val = n_int * phi**k
    err = float(abs(val - D_exact))
    if err < 0.5:
        trial = ROAD1_BASE * (1 + 1/(2*phi**27 + val))
        ppb_n = float((trial - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
        print(f"  D ~ {n_int} * phi^{k} = {nstr(val, 10)} (err={err:.4f}, {ppb_n:+.4f} ppb)")

# ============================================================================
# SECTION 16: HIGHER-ORDER CORRECTION FORMULAS
# ============================================================================
print()
print("=" * 80)
print("SECTION 16: HIGHER-ORDER CORRECTIONS")
print("=" * 80)

# exp(1/(2*phi^27))
inv_alpha_exp = ROAD1_BASE * mexp(C_original)
ppb_exp = float((inv_alpha_exp - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base * exp(1/(2*phi^27)): {nstr(inv_alpha_exp, 15)} ({ppb_exp:.4f} ppb)")

# 1 + x - x^2/2 (log series)
inv_alpha_log = ROAD1_BASE * (1 + C_original - C_original**2 / 2)
ppb_log = float((inv_alpha_log - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base * (1 + x - x^2/2): {nstr(inv_alpha_log, 15)} ({ppb_log:.4f} ppb)")

# 1 + x - x^2
inv_alpha_sq = ROAD1_BASE * (1 + C_original - C_original**2)
ppb_sq = float((inv_alpha_sq - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base * (1 + x - x^2): {nstr(inv_alpha_sq, 15)} ({ppb_sq:.4f} ppb)")

# cos(pi * x)
inv_alpha_cos = ROAD1_BASE * mcos(PI * C_original) / mcos(PI * C_original) * (1 + C_original)  # trivial
# Actually try: base / (1 - x)
inv_alpha_div = ROAD1_BASE / (1 - C_original)
ppb_div = float((inv_alpha_div - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base / (1 - x): {nstr(inv_alpha_div, 15)} ({ppb_div:.4f} ppb)")

# sqrt(1 + 2x) * base
inv_alpha_sqrt = ROAD1_BASE * msqrt(1 + 2 * C_original)
ppb_sqrt = float((inv_alpha_sqrt - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base * sqrt(1+2x): {nstr(inv_alpha_sqrt, 15)} ({ppb_sqrt:.4f} ppb)")

# (1+x)^(1-x) * base
inv_alpha_pow = ROAD1_BASE * (1 + C_original)**(1 - C_original)
ppb_pow = float((inv_alpha_pow - NIST_INV_ALPHA) / NIST_INV_ALPHA * 1e9)
print(f"base * (1+x)^(1-x): {nstr(inv_alpha_pow, 15)} ({ppb_pow:.4f} ppb)")

# ============================================================================
# SECTION 17: SMOKING GUN SEARCH
# ============================================================================
print()
print("=" * 80)
print("SECTION 17: SMOKING GUN -- DEFICIT CANCELLATION ANALYSIS")
print("=" * 80)

if abs(deficit_1loop) > 0:
    ratio_deficits = -deficit_2loop / deficit_1loop
    excess = ratio_deficits - 1
    print(f"Ratio |deficit_2loop/deficit_1loop| = {nstr(ratio_deficits, 10)}")
    print(f"Excess over exact cancellation = {nstr(excess, 10)}")
    print(f"Excess = {nstr(excess * 1000, 6)} per mille")
    print(f"Excess * pi = {nstr(excess * PI, 10)}")
    print(f"Excess * phi = {nstr(excess * phi, 10)}")
    print(f"Excess * phi^2 = {nstr(excess * phi**2, 10)}")
    print(f"R_ricci/(4*pi) = {nstr(R_ricci / (4*PI), 10)}")
    print(f"1/(4*phi^2) = {nstr(1/(4*phi**2), 10)}")
    print(f"Excess / R_ricci = {nstr(excess / R_ricci, 10)}")
    print(f"Excess * 600 = {nstr(excess * 600, 10)}")

# Is the gap precisely a curvature term times a topological invariant?
print(f"\n--- Does gap = (angular_deficit/2pi)^n * C_original * base? ---")
for n in range(1, 8):
    test = (angular_deficit_vertex / (2*PI))**n * C_original * ROAD1_BASE
    ppb_test = float(test / NIST_INV_ALPHA * 1e9)
    print(f"  n={n}: (pi/5/(2pi))^{n} * C * base = {nstr(test, 10)} = {ppb_test:.4f} ppb (gap = {nstr(gap_road1_road2, 6)})")

# What about the VOLUME correction to the propagator at the right scale?
print(f"\n--- Volume / curvature corrections at various scales ---")
for scale_name, r_scale in [("edge", a_edge_val), ("nn_geodesic", a_edge_val),
                              ("cell_radius", msqrt(3)*phi*a_edge_val), ("diameter", 15*a_edge_val)]:
    chi = r_scale / R_val if r_scale < PI * R_val else PI
    vol_r = msin(chi)**2 / chi**2 if chi > 0.001 else mpf(1)
    prop_r = chi / msin(chi) if chi > 0.001 else mpf(1)
    print(f"  scale={scale_name}: r={nstr(r_scale,4)}, chi={nstr(chi,4)}, vol_ratio={nstr(vol_r,10)}, prop_ratio={nstr(prop_r,10)}")

# ============================================================================
# SECTION 18: CURVATURE CORRECTION TO LOOP INTEGRAL
# ============================================================================
print()
print("=" * 80)
print("SECTION 18: CURVATURE CORRECTION TO LOOP INTEGRAL")
print("=" * 80)

# Compute I_flat vs I_lattice (sum of 1/(lambda+m^2)^2)
I_lattice = mpf(0)
for lam, mult in eigen_groups_mp:
    I_lattice += mpf(mult) / (mpf(str(lam)) + m_sq_current)**2
I_lattice /= N_VERTS

I_flat = 1 / (8 * PI * m_val)

print(f"I_flat (1/(8*pi*m)) = {nstr(I_flat, 15)}")
print(f"I_lattice (sum 1/(lam+m^2)^2 / N) = {nstr(I_lattice, 15)}")
print(f"I_lattice / I_flat = {nstr(I_lattice / I_flat, 15)}")
print(f"I_lattice - I_flat = {nstr(I_lattice - I_flat, 10)}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print()
print("=" * 80)
print("FINAL SUMMARY -- SORTED BY DISTANCE TO NIST")
print("=" * 80)

results = [
    ("Road 1 (algebraic)", ROAD1, float((ROAD1 - NIST_INV_ALPHA)/NIST_INV_ALPHA * 1e9)),
    ("Road 2 (spectral, xi=1/5)", ROAD2_CURRENT, float((ROAD2_CURRENT - NIST_INV_ALPHA)/NIST_INV_ALPHA * 1e9)),
    ("Road 1 * exp(x)", inv_alpha_exp, ppb_exp),
    ("base * (1+x-x^2/2)", inv_alpha_log, ppb_log),
    ("base * (1+x-x^2)", inv_alpha_sq, ppb_sq),
    ("base / (1-x)", inv_alpha_div, ppb_div),
    ("base * sqrt(1+2x)", inv_alpha_sqrt, ppb_sqrt),
    ("base * (1+x)^(1-x)", inv_alpha_pow, ppb_pow),
    ("Geometric series 1/(2phi^27-1)", inv_alpha_geom, ppb_geom),
]

# Add best correction formulas
for absp, ppb, formula, val_str in best_total[:10]:
    try:
        val_mpf = mpf(val_str)
        results.append((formula, val_mpf, ppb))
    except:
        pass

results.sort(key=lambda x: abs(x[2]))

print(f"{'Method':>55s} {'1/alpha':>18s} {'ppb':>10s}")
print("-" * 88)
for name, val, ppb in results[:30]:
    marker = " <--" if abs(ppb) < 0.02 else ""
    print(f"{name:>55s} {nstr(val, 15):>18s} {ppb:>+10.4f}{marker}")

print(f"\nNIST value:       {nstr(NIST_INV_ALPHA, 15)}")
print(f"NIST uncertainty: +/- {nstr(NIST_UNCERTAINTY, 6)} ({nstr(NIST_UNCERTAINTY/NIST_INV_ALPHA * 1e9, 3)} ppb)")

print("\n" + "=" * 80)
print("COMPUTATION COMPLETE")
print("=" * 80)
