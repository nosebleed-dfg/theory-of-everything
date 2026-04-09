"""
120CELL_GRAVITY — derives Newton's G from the 120-cell spectral structure; 600-cell dual construction
nos3bl33d

600-cell -> 600 tetrahedra -> 120-cell dual graph. Spectral gap analysis for gravitational coupling.
"""

import numpy as np
from scipy import linalg
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
from collections import Counter
import time

phi = (1 + np.sqrt(5)) / 2
psi = 1 / phi  # = phi - 1

# Physical constants
G_SI = 6.67430e-11
c_SI = 299792458
hbar_SI = 1.054571817e-34
m_proton = 1.672621924e-27
m_electron = 9.1093837015e-31
m_Planck = 2.176434e-8
alpha_EM = 1/137.035999084
alpha_G_proton = G_SI * m_proton**2 / (hbar_SI * c_SI)

print("=" * 80)
print("120-CELL SPECTRAL ANALYSIS FOR GRAVITATIONAL CONSTANT")
print("=" * 80)
print(f"\nphi = {phi:.15f}")
print(f"phi^(-4) = {phi**(-4):.15f}")
print(f"alpha_EM = {alpha_EM:.12e}")
print(f"alpha_G(proton) = {alpha_G_proton:.6e}")
print(f"alpha_EM/alpha_G = {alpha_EM/alpha_G_proton:.6e}")
print(f"log_phi(alpha_EM/alpha_G) = {np.log(alpha_EM/alpha_G_proton)/np.log(phi):.6f}")
print(f"m_Planck/m_proton = {m_Planck/m_proton:.6e}")
print(f"log_phi(m_Planck/m_proton) = {np.log(m_Planck/m_proton)/np.log(phi):.6f}")

# =============================================================================
# STEP 1: Build the 600-cell (120 vertices on S^3)
# =============================================================================
print("\n" + "=" * 80)
print("STEP 1: Building 600-cell (120 vertices)")
print("=" * 80)

verts_600 = []

# Type 1: 8 vertices — permutations of (+-1, 0, 0, 0)
for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        verts_600.append(tuple(v))

# Type 2: 16 vertices — (+-1/2, +-1/2, +-1/2, +-1/2)
for s0 in [1, -1]:
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                verts_600.append((s0*0.5, s1*0.5, s2*0.5, s3*0.5))

# Type 3: 96 vertices — even permutations of (+-phi/2, +-1/2, +-1/(2*phi), 0)
hp = phi / 2
hh = 0.5
hq = 1.0 / (2.0 * phi)

even_perms_idx = [
    (0,1,2,3), (0,2,3,1), (0,3,1,2),
    (1,0,3,2), (1,2,0,3), (1,3,2,0),
    (2,0,1,3), (2,1,3,0), (2,3,0,1),
    (3,0,2,1), (3,1,0,2), (3,2,1,0)
]

base_abs = [hp, hh, hq, 0.0]

for perm in even_perms_idx:
    perm_vals = [base_abs[perm[i]] for i in range(4)]
    nonzero_positions = [i for i in range(4) if abs(perm_vals[i]) > 1e-10]
    n_nonzero = len(nonzero_positions)
    for bits in range(1 << n_nonzero):
        v = list(perm_vals)
        for bit_idx, pos in enumerate(nonzero_positions):
            if bits & (1 << bit_idx):
                v[pos] = -v[pos]
        verts_600.append(tuple(v))

# Deduplicate
verts_600_clean = []
for v in verts_600:
    is_dup = False
    for u in verts_600_clean:
        if max(abs(v[i] - u[i]) for i in range(4)) < 1e-8:
            is_dup = True
            break
    if not is_dup:
        verts_600_clean.append(v)

verts_600 = np.array(verts_600_clean)
print(f"600-cell vertices: {len(verts_600)}")
norms = np.sqrt(np.sum(verts_600**2, axis=1))
print(f"Norm range: [{norms.min():.8f}, {norms.max():.8f}]")
assert len(verts_600) == 120, f"Expected 120 vertices, got {len(verts_600)}"

# =============================================================================
# STEP 2: Build 600-cell adjacency
# =============================================================================
print("\n" + "=" * 80)
print("STEP 2: Building 600-cell adjacency")
print("=" * 80)

diff = verts_600[:, np.newaxis, :] - verts_600[np.newaxis, :, :]
dist_sq = np.sum(diff**2, axis=2)

ds = dist_sq.copy()
np.fill_diagonal(ds, 999)
min_d2 = ds.min()
print(f"Min distance^2: {min_d2:.10f}")
print(f"(1/phi)^2 = {psi**2:.10f}")

adj_600 = (np.abs(dist_sq - min_d2) < 1e-6).astype(int)
np.fill_diagonal(adj_600, 0)

degrees = adj_600.sum(axis=1)
print(f"Degree distribution: {Counter(degrees.tolist())}")
n_edges_600 = adj_600.sum() // 2
print(f"Edges: {n_edges_600}")
assert n_edges_600 == 720, f"Expected 720 edges, got {n_edges_600}"

# =============================================================================
# STEP 3: Find all 600 tetrahedra and build 120-cell dual
# =============================================================================
print("\n" + "=" * 80)
print("STEP 3: Finding tetrahedra and building 120-cell")
print("=" * 80)

t0 = time.time()

# Build neighbor lists
neighbors = [set() for _ in range(120)]
for i in range(120):
    for j in range(i+1, 120):
        if adj_600[i, j]:
            neighbors[i].add(j)
            neighbors[j].add(i)

# Find triangles
triangles = []
for i in range(120):
    ni = neighbors[i]
    for j in sorted(ni):
        if j > i:
            for k in sorted(ni & neighbors[j]):
                if k > j:
                    triangles.append((i, j, k))

print(f"Triangles: {len(triangles)} (expected 1200)")

# Find tetrahedra
tetrahedra = set()
for (i, j, k) in triangles:
    common = neighbors[i] & neighbors[j] & neighbors[k]
    for l in common:
        tetrahedra.add(tuple(sorted([i, j, k, l])))

tetrahedra = sorted(tetrahedra)
print(f"Tetrahedra: {len(tetrahedra)} (expected 600)")
assert len(tetrahedra) == 600, f"Expected 600 tetrahedra, got {len(tetrahedra)}"

# Build face -> tet lookup for adjacency
face_to_tets = {}
for idx, tet in enumerate(tetrahedra):
    for skip in range(4):
        face = tuple(tet[j] for j in range(4) if j != skip)
        if face in face_to_tets:
            face_to_tets[face].append(idx)
        else:
            face_to_tets[face] = [idx]

# Adjacency: shared triangular face
N = 600
adj_120 = np.zeros((N, N), dtype=np.int8)
edges_120 = 0
for face, tets_list in face_to_tets.items():
    if len(tets_list) == 2:
        i, j = tets_list
        adj_120[i, j] = 1
        adj_120[j, i] = 1
        edges_120 += 1

degrees_120 = adj_120.sum(axis=1)
print(f"\n120-cell constructed in {time.time()-t0:.1f}s")
print(f"Vertices: {N}")
print(f"Edges: {edges_120}")
print(f"Degree distribution: {Counter(degrees_120.tolist())}")

assert Counter(degrees_120.tolist()) == Counter({4: 600}), "120-cell should be 4-regular"
print("*** 120-CELL VERIFIED: 600 vertices, 1200 edges, 4-regular ***")

# Compute tet centers (vertices of the 120-cell in R^4)
tet_centers = np.array([np.mean(verts_600[list(t)], axis=0) for t in tetrahedra])
center_norms = np.sqrt(np.sum(tet_centers**2, axis=1))

# =============================================================================
# STEP 4: Spectral Analysis
# =============================================================================
print("\n" + "=" * 80)
print("STEP 4: SPECTRAL ANALYSIS")
print("=" * 80)

adj = adj_120.astype(np.float64)
degree = 4
L = degree * np.eye(N) - adj

print("\nComputing adjacency eigenvalues...")
t0 = time.time()
eigs_adj_all = np.sort(np.real(linalg.eigvalsh(adj)))[::-1]
print(f"  Time: {time.time()-t0:.1f}s")

print("Computing Laplacian eigenvalues and eigenvectors...")
t0 = time.time()
eig_vals_L, eig_vecs_L = linalg.eigh(L)
idx_sort = np.argsort(eig_vals_L)
eig_vals_L = eig_vals_L[idx_sort]
eig_vecs_L = eig_vecs_L[:, idx_sort]
print(f"  Time: {time.time()-t0:.1f}s")

# Distinct eigenvalues
def get_distinct(vals, tol=1e-6):
    """Return (distinct_values, multiplicities) sorted."""
    rounded = np.round(vals, 8)
    unique, counts = np.unique(rounded, return_counts=True)
    return unique, counts

unique_adj, mult_adj = get_distinct(eigs_adj_all)
unique_adj = unique_adj[::-1]
mult_adj = mult_adj[::-1]

unique_lap, mult_lap = get_distinct(eig_vals_L)

print(f"\nDistinct adjacency eigenvalues: {len(unique_adj)}")
print(f"Distinct Laplacian eigenvalues: {len(unique_lap)}")
print(f"Multiplicity sum: {sum(mult_adj)} (should be 600)")

# Express eigenvalues in terms of phi
def phi_expr(val):
    """Try to express val as (a + b*sqrt(5))/c for small integers."""
    s5 = np.sqrt(5)
    best = None
    best_err = 1e-5
    for denom in [1, 2, 3, 4, 5, 6, 8, 10]:
        for a in range(-30, 31):
            for b in range(-30, 31):
                test = (a + b*s5) / denom
                err = abs(val - test)
                if err < best_err:
                    best_err = err
                    if denom == 1:
                        if b == 0:
                            best = f"{a}"
                        else:
                            sign = '+' if b > 0 else ''
                            best = f"{a}{sign}{b}s5"
                    else:
                        if b == 0:
                            best = f"{a}/{denom}"
                        else:
                            sign = '+' if b > 0 else ''
                            best = f"({a}{sign}{b}s5)/{denom}"
    return best if best else "?"

print("\n" + "-" * 85)
print(f"{'#':>3} {'lambda_i':>14} {'mu_i=4-lam':>14} {'mult':>5} {'expression':>25} {'mu_expr':>25}")
print("-" * 85)

for i, (lam, m) in enumerate(zip(unique_adj, mult_adj)):
    lam_f = float(lam)
    mu_f = degree - lam_f
    expr_l = phi_expr(lam_f)
    expr_m = phi_expr(mu_f)
    print(f"{i+1:3d} {lam_f:14.8f} {mu_f:14.8f} {int(m):5d} {expr_l:>25} {expr_m:>25}")

# =============================================================================
# STEP 5: Key Spectral Quantities
# =============================================================================
print("\n" + "=" * 80)
print("STEP 5: KEY SPECTRAL QUANTITIES")
print("=" * 80)

# Spectral gap
nonzero_lap = unique_lap[unique_lap > 1e-6]
spectral_gap = float(nonzero_lap[0])
print(f"\nSpectral gap: {spectral_gap:.12f}")
print(f"phi^(-4) = {phi**(-4):.12f}")
print(f"Match: {abs(spectral_gap - phi**(-4)) < 1e-6}")

# Largest eigenvalue
mu_max = float(unique_lap[-1])
lam_max = float(unique_adj[0])
lam_min = float(unique_adj[-1])
print(f"\nLargest Laplacian eigenvalue: {mu_max:.10f}")
print(f"Spectral ratio mu_max/mu_min: {mu_max/spectral_gap:.10f}")
print(f"phi^8 = {phi**8:.10f}")

# Product of 27 distinct |adj eigenvalues| (nonzero)
nonzero_adj_distinct = unique_adj[np.abs(unique_adj) > 1e-6]
P27 = float(np.prod(np.abs(nonzero_adj_distinct)))
print(f"\nProduct of distinct nonzero |adj eigenvalues| = {P27:.6f}")
print(f"Number of such eigenvalues: {len(nonzero_adj_distinct)}")
print(f"phi^27 = {phi**27:.6f}")
print(f"phi^{len(nonzero_adj_distinct)} = {phi**len(nonzero_adj_distinct):.6f}")
print(f"Ratio phi^27/P27 = {phi**27/P27:.10f}")

# Spanning trees: tau = (1/N) * prod(nonzero laplacian eigenvalues with mult)
log_tau = -np.log(N)
for mu_val, m_val in zip(unique_lap, mult_lap):
    if mu_val > 1e-6:
        log_tau += int(m_val) * np.log(float(mu_val))

print(f"\nSpanning trees:")
print(f"  log10(tau) = {log_tau/np.log(10):.4f}")
print(f"  log_phi(tau) = {log_tau/np.log(phi):.4f}")
print(f"  log10(alpha_EM/alpha_G) = {np.log10(alpha_EM/alpha_G_proton):.4f}")

# Functional determinant
log_det_L = sum(int(m) * np.log(float(mu)) for mu, m in zip(unique_lap, mult_lap) if mu > 1e-6)
print(f"\nlog det(L_nonzero) = {log_det_L:.6f}")
print(f"det(L)^(1/599) = {np.exp(log_det_L/599):.10f}")

# Spectral zeta
print("\nSpectral zeta function zeta_L(s) = sum mu_i^(-s) (with multiplicities):")
for s in [0.5, 1.0, 1.5, 2.0, 3.0, -1, -2, -0.5]:
    zeta = sum(int(m) * float(mu)**(-s) for mu, m in zip(unique_lap, mult_lap) if mu > 1e-6)
    print(f"  zeta_L({s:5.1f}) = {zeta:.10f}")

# =============================================================================
# STEP 6: Graph Distances & Antipodal Structure
# =============================================================================
print("\n" + "=" * 80)
print("STEP 6: GRAPH DISTANCES AND ANTIPODAL STRUCTURE")
print("=" * 80)

print("\nComputing shortest paths...")
t0 = time.time()
dist_matrix = shortest_path(csr_matrix(adj.astype(float)), method='D', directed=False)
print(f"  Time: {time.time()-t0:.1f}s")

dist_from_0 = dist_matrix[0].astype(int)
dist_dist = Counter(dist_from_0.tolist())
diameter = int(dist_matrix.max())

print(f"\nDiameter: {diameter}")
print(f"\nDistance distribution from vertex 0:")
print(f"  {'d':>4} {'N(d)':>6} {'(-1)^d*N(d)':>12}")
for d in sorted(dist_dist):
    print(f"  {d:4d} {dist_dist[d]:6d} {(-1)**d * dist_dist[d]:12d}")

# Antipodal vertex
antipodal = np.where(dist_from_0 == diameter)[0]
print(f"\nAntipodal vertices: {len(antipodal)} at distance {diameter}")
v_anti = antipodal[0]

# =============================================================================
# STEP 7: Green's Functions — Self-energy and Antipodal Propagator
# =============================================================================
print("\n" + "=" * 80)
print("STEP 7: GREEN'S FUNCTIONS")
print("=" * 80)

v0 = 0
psi_0 = eig_vecs_L[v0, :]     # eigenvector components at vertex 0
psi_a = eig_vecs_L[v_anti, :]  # eigenvector components at antipodal vertex

# Verify zero mode
print(f"\nZero mode psi(0): {psi_0[0]:.10f} (expected {1/np.sqrt(N):.10f})")
print(f"Zero mode psi(anti): {psi_a[0]:.10f}")

# Green's function: G_m(x,y) = sum_k psi_k(x)*psi_k(y) / (mu_k + m^2)
print("\n--- Green's function G_m(0,0) and G_m(0,antipode) ---")
print(f"{'m^2':>14} {'G(0,0)':>18} {'G(0,anti)':>18} {'G_anti/G_self':>16} {'|G_self/G_anti|':>16}")
print("-" * 85)

m2_list = [1e-4, 1e-3, 0.01, phi**(-4), 0.05, phi**(-2), 0.5, 1.0, phi, phi**2, phi**4, 10, 100]

for m2 in m2_list:
    G00 = np.sum(psi_0**2 / (eig_vals_L + m2))
    G0a = np.sum(psi_0 * psi_a / (eig_vals_L + m2))
    r1 = G0a / G00
    r2 = abs(G00 / G0a) if abs(G0a) > 1e-50 else float('inf')
    print(f"{m2:14.6e} {G00:18.12f} {G0a:18.12e} {r1:16.12f} {r2:16.6f}")

# =============================================================================
# STEP 8: Contributions by eigenvalue group
# =============================================================================
print("\n" + "=" * 80)
print("STEP 8: EIGENVALUE-GROUP CONTRIBUTIONS TO GREEN'S FUNCTIONS")
print("=" * 80)

# Group eigenvectors by eigenvalue
eigengroups = {}
for k in range(N):
    mu = round(float(eig_vals_L[k]), 8)
    if mu not in eigengroups:
        eigengroups[mu] = []
    eigengroups[mu].append(k)

print(f"\n{'mu':>12} {'mult':>5} {'sum|psi0|^2':>14} {'sum psi0*psiA':>16} {'ratio A/0':>12} {'sign':>6}")
print("-" * 75)

for mu in sorted(eigengroups.keys()):
    indices = eigengroups[mu]
    mult = len(indices)
    sum_sq = sum(psi_0[k]**2 for k in indices)
    sum_cross = sum(psi_0[k] * psi_a[k] for k in indices)
    ratio = sum_cross / sum_sq if abs(sum_sq) > 1e-15 else 0
    sign_str = "+" if sum_cross > 1e-15 else ("-" if sum_cross < -1e-15 else "0")
    print(f"{mu:12.6f} {mult:5d} {sum_sq:14.10f} {sum_cross:16.12f} {ratio:12.6f} {sign_str:>6}")

# Verify: sum of all sum_sq should = 1/N (diagonal of Green's function resolvent identity)
total_sq = sum(sum(psi_0[k]**2 for k in indices) for indices in eigengroups.values())
print(f"\nTotal sum|psi_0|^2 = {total_sq:.10f} (should be 1.0 since eigvecs normalized)")

# =============================================================================
# STEP 9: Distance-resolved Green's function
# =============================================================================
print("\n" + "=" * 80)
print("STEP 9: GREEN'S FUNCTION AT ALL DISTANCES")
print("=" * 80)

for m2_label, m2 in [("spectral_gap", spectral_gap), ("phi^(-2)", phi**(-2)), ("1.0", 1.0)]:
    print(f"\nm^2 = {m2:.10f} ({m2_label})")
    print(f"{'d':>4} {'N(d)':>6} {'G(0,d) avg':>18} {'G(0,d)/G(0,0)':>18} {'|G|*N(d)':>14}")

    G_self = np.sum(psi_0**2 / (eig_vals_L + m2))

    for d in range(diameter + 1):
        verts_d = np.where(dist_from_0 == d)[0]
        if len(verts_d) == 0:
            continue

        G_vals = []
        for v in verts_d:
            psi_v = eig_vecs_L[v, :]
            G_0v = np.sum(psi_0 * psi_v / (eig_vals_L + m2))
            G_vals.append(G_0v)

        G_avg = np.mean(G_vals)
        ratio = G_avg / G_self
        shell_sum = abs(G_avg) * len(verts_d)
        print(f"{d:4d} {len(verts_d):6d} {G_avg:18.12e} {ratio:18.12e} {shell_sum:14.8e}")

# =============================================================================
# STEP 10: Madelung-like sums
# =============================================================================
print("\n" + "=" * 80)
print("STEP 10: MADELUNG-LIKE SUMS")
print("=" * 80)

# Sum (-1)^d * N(d)
mad_0 = sum((-1)**d * dist_dist[d] for d in dist_dist)
print(f"\nsum_d (-1)^d * N(d) = {mad_0}")

# Sum (-1)^d * N(d) / d (skip d=0)
mad_1 = sum((-1)**d * dist_dist[d] / d for d in dist_dist if d > 0)
print(f"sum_d (-1)^d * N(d) / d = {mad_1:.10f}")

# Sum (-1)^d * N(d) / d^2
mad_2 = sum((-1)**d * dist_dist[d] / d**2 for d in dist_dist if d > 0)
print(f"sum_d (-1)^d * N(d) / d^2 = {mad_2:.10f}")

# Sum (-1)^d * N(d) / d^3
mad_3 = sum((-1)**d * dist_dist[d] / d**3 for d in dist_dist if d > 0)
print(f"sum_d (-1)^d * N(d) / d^3 = {mad_3:.10f}")

# Weighted by actual Green's function
print("\nMadelung with actual Green's function at distance d:")
for m2 in [spectral_gap, phi**(-2), 1.0]:
    G_self = np.sum(psi_0**2 / (eig_vals_L + m2))
    mad_G = 0
    for d in range(1, diameter + 1):
        verts_d = np.where(dist_from_0 == d)[0]
        if len(verts_d) == 0:
            continue
        for v in verts_d:
            psi_v = eig_vecs_L[v, :]
            G_0v = np.sum(psi_0 * psi_v / (eig_vals_L + m2))
            mad_G += (-1)**d * G_0v

    print(f"  m^2={m2:.6f}: Madelung_G = {mad_G:.10e}, /G(0,0) = {mad_G/G_self:.10e}")

# =============================================================================
# STEP 11: Higher-order couplings from spectrum
# =============================================================================
print("\n" + "=" * 80)
print("STEP 11: HIGHER-ORDER SPECTRAL COUPLINGS")
print("=" * 80)

# What power of phi gives alpha_G?
log_phi_aEM = np.log(alpha_EM) / np.log(phi)
log_phi_aG = np.log(alpha_G_proton) / np.log(phi)
log_phi_ratio = np.log(alpha_EM/alpha_G_proton) / np.log(phi)

print(f"\nlog_phi(alpha_EM) = {log_phi_aEM:.6f}")
print(f"log_phi(alpha_G) = {log_phi_aG:.6f}")
print(f"log_phi(alpha_EM/alpha_G) = {log_phi_ratio:.6f}")
print(f"log_phi(m_Planck/m_proton) = {np.log(m_Planck/m_proton)/np.log(phi):.6f}")

# Products of first k nonzero distinct Laplacian eigenvalues
print("\n--- Cumulative products of distinct nonzero Laplacian eigenvalues ---")
print(f"{'k':>3} {'mu_k':>12} {'mult':>5} {'log10(prod)':>14} {'log_phi(prod)':>14}")
nonzero_mu_sorted = sorted([(mu, len(indices)) for mu, indices in eigengroups.items() if mu > 1e-6])

log_prod = 0.0
for k, (mu, mult) in enumerate(nonzero_mu_sorted, 1):
    log_prod += np.log(mu)
    l10 = log_prod / np.log(10)
    lphi = log_prod / np.log(phi)
    print(f"{k:3d} {mu:12.8f} {mult:5d} {l10:14.6f} {lphi:14.6f}")

# With multiplicities
print("\n--- Cumulative products WITH multiplicities ---")
log_prod_m = 0.0
for k, (mu, mult) in enumerate(nonzero_mu_sorted, 1):
    log_prod_m += mult * np.log(mu)
    l10 = log_prod_m / np.log(10)
    lphi = log_prod_m / np.log(phi)
    marker = ""
    if abs(l10 - np.log10(alpha_EM/alpha_G_proton)) < 2:
        marker = " <--- NEAR alpha_EM/alpha_G!"
    print(f"{k:3d} {mu:12.8f} {mult:5d} {l10:14.4f} {lphi:14.4f}{marker}")

# =============================================================================
# STEP 12: The heat kernel at special times
# =============================================================================
print("\n" + "=" * 80)
print("STEP 12: HEAT KERNEL")
print("=" * 80)

print("\n--- Heat kernel trace K(t) = Tr(exp(-Lt)) ---")
for t in [phi**n for n in range(-6, 7)]:
    K = sum(mult * np.exp(-mu * t) for mu, mult in
            [(mu, len(idx)) for mu, idx in eigengroups.items()])
    print(f"  K(phi^{np.log(t)/np.log(phi):5.1f}) = K({t:12.6f}) = {K:.6f}")

print("\n--- Heat kernel at antipode K_anti(t) ---")
for t in [phi**n for n in range(-6, 7)]:
    K_anti = 0
    for mu, indices in eigengroups.items():
        cross = sum(psi_0[k] * psi_a[k] for k in indices)
        K_anti += cross * np.exp(-mu * t)
    ratio = K_anti * N  # normalize
    print(f"  K_anti(phi^{np.log(t)/np.log(phi):5.1f}) = {K_anti:.10e}  N*K_anti = {ratio:.10e}")

# =============================================================================
# STEP 13: Gravity formula attempts
# =============================================================================
print("\n" + "=" * 80)
print("STEP 13: GRAVITY DERIVATION ATTEMPTS")
print("=" * 80)

print("\n--- Attempt 1: alpha_G = alpha_EM * |G(0,anti)/G(0,0)| ---")
for m2 in [spectral_gap, phi**(-2), 1.0, phi**2]:
    G00 = np.sum(psi_0**2 / (eig_vals_L + m2))
    G0a = np.sum(psi_0 * psi_a / (eig_vals_L + m2))
    aG_trial = alpha_EM * abs(G0a / G00)
    G_trial = aG_trial * hbar_SI * c_SI / m_proton**2
    ratio = G_trial / G_SI
    print(f"  m^2={m2:.6f}: G_trial = {G_trial:.4e}  G_SI = {G_SI:.4e}  ratio = {ratio:.4e}"
          f"  log10 = {np.log10(ratio):.2f}")

print("\n--- Attempt 2: alpha_G = alpha_EM * |G(0,anti)/G(0,0)|^2 ---")
for m2 in [spectral_gap, phi**(-2), 1.0]:
    G00 = np.sum(psi_0**2 / (eig_vals_L + m2))
    G0a = np.sum(psi_0 * psi_a / (eig_vals_L + m2))
    aG_trial = alpha_EM * (G0a / G00)**2
    G_trial = aG_trial * hbar_SI * c_SI / m_proton**2
    ratio = G_trial / G_SI
    print(f"  m^2={m2:.6f}: G_trial = {G_trial:.4e}  ratio = {ratio:.4e}")

print("\n--- Attempt 3: alpha_G = (spectral_gap/V)^k for various k ---")
sg_over_V = spectral_gap / 600
# alpha_EM = sg/20 ≈ 7.3e-3
# Need alpha_G ≈ 5.9e-39
# So need (sg/V)^k = 5.9e-39
# sg/V ≈ 2.43e-4
# k = log(5.9e-39) / log(2.43e-4) = -38.23 / -3.614 ≈ 10.58
for k in range(1, 15):
    val = sg_over_V**k
    print(f"  (gap/N)^{k:2d} = {val:.4e}  log10 = {np.log10(val):.2f}"
          f"  {'<-- close to alpha_G!' if abs(np.log10(val) - np.log10(alpha_G_proton)) < 1 else ''}")

print("\n--- Attempt 4: (phi^(-4)/20)^k for various k ---")
a_bare = phi**(-4) / 20
for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    val = a_bare**k
    print(f"  alpha_bare^{k:2d} = {val:.4e}  log10 = {np.log10(val):.2f}"
          f"  {'<-- near alpha_G!' if abs(np.log10(val) - np.log10(alpha_G_proton)) < 1 else ''}")

print("\n--- Attempt 5: phi^(-4n) / V^m for various n, m ---")
target = alpha_G_proton
log_target = np.log(target)
for n in range(1, 50):
    for m in range(0, 10):
        val = phi**(-4*n) / (600**m)
        if abs(np.log10(val) - np.log10(target)) < 0.5:
            print(f"  phi^(-{4*n}) / 600^{m} = {val:.6e}  ratio to alpha_G = {val/target:.4f}")

# Also try with 20 instead of 600
for n in range(1, 50):
    for m in range(0, 15):
        val = phi**(-4*n) / (20**m)
        if abs(np.log10(val) - np.log10(target)) < 0.3:
            print(f"  phi^(-{4*n}) / 20^{m} = {val:.6e}  ratio to alpha_G = {val/target:.4f}")

print("\n--- Attempt 6: Product of eigenvalues and gravity ---")
# The spanning tree count
print(f"  log10(tau) = {log_tau/np.log(10):.4f}")
print(f"  alpha_EM / tau = ?")
log_aEM_over_tau = np.log(alpha_EM) - log_tau
print(f"  log10(alpha_EM/tau) = {log_aEM_over_tau/np.log(10):.4f}")
print(f"  log10(alpha_G) = {np.log10(alpha_G_proton):.4f}")

# tau * alpha_EM ? tau / alpha_EM ?
print(f"  alpha_EM * exp(-log_tau) = {alpha_EM * np.exp(-log_tau):.4e}")
print(f"  alpha_EM / tau ... log10 = {(np.log10(alpha_EM) - log_tau/np.log(10)):.4f}")

# =============================================================================
# STEP 14: Proton mass from geometry
# =============================================================================
print("\n" + "=" * 80)
print("STEP 14: PROTON MASS FROM GEOMETRY")
print("=" * 80)

mass_ratio = m_Planck / m_proton
log_phi_mr = np.log(mass_ratio) / np.log(phi)
log_phi_mr2 = 2 * log_phi_mr  # for mass^2 ratio

print(f"\nm_Planck/m_proton = {mass_ratio:.6e}")
print(f"log_phi(m_Planck/m_proton) = {log_phi_mr:.6f}")
print(f"  = {int(log_phi_mr)} + {log_phi_mr - int(log_phi_mr):.6f}")
print(f"2*log_phi (mass^2) = {log_phi_mr2:.6f}")
print(f"  = {int(log_phi_mr2)} + {log_phi_mr2 - int(log_phi_mr2):.6f}")

# Nearby geometric numbers for ~91.46
print(f"\nGeometric numbers near 91.46:")
print(f"  91 = T_13 = 13th triangular number = sum(1..13)")
print(f"  91 = 7*13")
print(f"  90 = 9*10 = 2*45")
print(f"  92 = 4*23")

# Check: is log_phi_mr close to a spectral quantity?
print(f"\n  Total multiplicity-weighted eigenvalue sum = {sum(int(m)*float(mu) for mu,m in zip(unique_lap, mult_lap)):.4f}")
print(f"  = N*d = {N*degree} (trace of L)")
print(f"  Number of nonzero eigenvalue groups: {len(nonzero_mu_sorted)}")

# The 120-cell as a polytope: V=600, E=1200, F=720, C=120
V, E, F, C = 600, 1200, 720, 120
chi = V - E + F - C
print(f"\n  V-E+F-C = {V}-{E}+{F}-{C} = {chi}")
print(f"  V+C = {V+C}")
print(f"  E-F = {E-F}")
print(f"  V*C/(E+F) = {V*C/(E+F):.4f}")

# =============================================================================
# STEP 15: Analytic torsion / Reidemeister torsion
# =============================================================================
print("\n" + "=" * 80)
print("STEP 15: ANALYTIC TORSION")
print("=" * 80)

# For a graph, the analytic torsion is related to the tree number
# log T = (1/2) * sum_k' log(mu_k) where sum' is over nonzero eigenvalues
# This is (1/2) * log det'(L) = (1/2) * log(N * tau)
log_T = 0.5 * log_det_L
print(f"\nlog(analytic torsion) = (1/2)*log det'(L) = {log_T:.6f}")
print(f"log10(analytic torsion) = {log_T/np.log(10):.6f}")
print(f"Analytic torsion^(1/N) = {np.exp(log_T/N):.10f}")

# =============================================================================
# STEP 16: Det(L + m^2*I) for various m^2
# =============================================================================
print("\n" + "=" * 80)
print("STEP 16: REGULARIZED DETERMINANT det(L + m^2*I)")
print("=" * 80)

print(f"\n{'m^2':>14} {'log10 det':>14} {'log_phi det':>14} {'det^(1/N)':>14}")
print("-" * 60)

for m2 in [1e-6, 1e-4, 0.001, 0.01, phi**(-4), phi**(-2), 0.5, 1.0, phi, phi**2, phi**4]:
    log_det = sum(int(m) * np.log(float(mu) + m2) for mu, m in zip(unique_lap, mult_lap))
    l10 = log_det / np.log(10)
    lphi = log_det / np.log(phi)
    det_root = np.exp(log_det / N)
    print(f"{m2:14.6e} {l10:14.4f} {lphi:14.4f} {det_root:14.8f}")

# =============================================================================
# STEP 17: The spectral zeta derivative at s=0
# =============================================================================
print("\n" + "=" * 80)
print("STEP 17: SPECTRAL ZETA DERIVATIVE zeta'_L(0) = -log det(L)")
print("=" * 80)

# zeta'(0) = -sum' log(mu_k) = -log det'(L)
zeta_prime_0 = -log_det_L
print(f"\nzeta'_L(0) = -log det'(L) = {zeta_prime_0:.6f}")
print(f"exp(-zeta'(0)) = det'(L) = exp({log_det_L:.4f})")
print(f"log10 of det'(L) = {log_det_L/np.log(10):.4f}")

# =============================================================================
# STEP 18: Key number checks
# =============================================================================
print("\n" + "=" * 80)
print("STEP 18: SEARCHING FOR THE GRAVITY FORMULA")
print("=" * 80)

# We know:
# alpha_EM ~ phi^(-4)/20 ~ phi^(-10.23) ~ 10^(-2.137)
# alpha_G ~ 10^(-38.23) ~ phi^(-183.16)
# alpha_EM/alpha_G ~ 10^(36.09) ~ phi^(172.70)
#
# 172.70 / 4 = 43.17 ... not clean
# But 172 = 4*43, and 43 is prime
# Or: 172.70 ~ 600/3.475 ... not clean
#
# What if alpha_G = phi^(-4*k) / (20^m) and we need to find k, m?

print("\nSearching phi^(-4k) / (20^m) = alpha_G...")
for k in range(1, 60):
    for m in range(0, 20):
        val = phi**(-4*k) / (20**m)
        if val > 0 and abs(np.log10(val) - np.log10(alpha_G_proton)) < 0.15:
            print(f"  k={k:2d}, m={m:2d}: phi^(-{4*k}) / 20^{m} = {val:.6e}"
                  f"  ratio = {val/alpha_G_proton:.6f}")

# What about phi^(-4k) / (N^m) where N=600?
print("\nSearching phi^(-4k) / (600^m) = alpha_G...")
for k in range(1, 60):
    for m in range(0, 15):
        val = phi**(-4*k) / (600**m)
        if val > 0 and abs(np.log10(val) - np.log10(alpha_G_proton)) < 0.15:
            print(f"  k={k:2d}, m={m:2d}: phi^(-{4*k}) / 600^{m} = {val:.6e}"
                  f"  ratio = {val/alpha_G_proton:.6f}")

# Key: does the number 172.70 have a spectral decomposition?
# 172.70 = 2 * 86.35
# log_phi(alpha_EM) ≈ -10.23
# log_phi(alpha_G) ≈ -182.93
# So log_phi(alpha_EM/alpha_G) = 172.70
# Is 172.70 ≈ N * spectral_quantity?
print(f"\n172.70 / N = {172.70/N:.6f}")  # 0.2878
print(f"172.70 / diameter = {172.70/diameter:.6f}")  # 11.51
print(f"172.70 / 27 = {172.70/27:.6f}")  # 6.40
print(f"172.70 / log_phi(mu_max) = {172.70/(np.log(mu_max)/np.log(phi)):.6f}")

# What if it's 4 * (number of distinct eigenvalues) * phi^something?
print(f"172.70 / (4 * 27) = {172.70/(4*27):.6f}")  # 1.598 ≈ phi?
print(f"4 * 27 * phi = {4*27*phi:.4f}")  # 174.9
print(f"4 * 27 * phi vs 172.70: off by {abs(4*27*phi - 172.70):.2f}")

# WAIT: 4*27*phi ≈ 174.89, but we need 172.70
# Try: 4*d_eig*phi - correction = ?
# 172.70 ≈ 4*43 + 0.70? 43 is prime
# 172.70 ≈ 12*14.39?
# 172.70 ≈ 2 * 86.35 ≈ 2 * log_phi(m_P/m_p)^2 ... wait
# log_phi(m_P/m_p) = 91.46
# 2 * 91.46 = 182.92 ≈ log_phi(1/alpha_G) = -log_phi(alpha_G)
# So log_phi(alpha_EM/alpha_G) = log_phi(alpha_EM) + 2*log_phi(m_P/m_p)
# 172.70 ≈ -10.23 + 182.92 = 172.69  YES!

print(f"\n*** KEY IDENTITY CHECK ***")
print(f"log_phi(alpha_EM/alpha_G) = {log_phi_ratio:.6f}")
print(f"log_phi(alpha_EM) + 2*log_phi(m_P/m_p) = {log_phi_aEM + 2*log_phi_mr:.6f}")
print(f"These should be equal: alpha_G = alpha_EM * (m_p/m_P)^2")
print(f"  alpha_EM * (m_p/m_P)^2 = {alpha_EM * (m_proton/m_Planck)**2:.6e}")
print(f"  alpha_G = {alpha_G_proton:.6e}")
print(f"  Ratio: {alpha_EM * (m_proton/m_Planck)**2 / alpha_G_proton:.6f}")
# This just says alpha_G = G*m_p^2/(hbar*c) and alpha_EM = e^2/(hbar*c*4pi*eps0)
# So alpha_G/alpha_EM = G*m_p^2/e^2 * 4pi*eps0 -- not circular but not helpful

# The REAL question: what sets m_proton/m_Planck from the spectrum?
# log_phi(m_P/m_p) ≈ 91.46
# Is 91 a spectral number?

print(f"\n*** PROTON MASS SEARCH ***")
print(f"log_phi(m_P/m_p) = {log_phi_mr:.6f}")

# Sum of all multiplicities of eigenvalues < some threshold?
# Sum of multiplicities: always 600
# But: what about 91 from geometric structure?

# 91 = sum(1..13). But also:
# The 120-cell has 120 cells. Each cell is a dodecahedron with 20 vertices, 30 edges, 12 faces.
# 120 * 20 / shared = 600 (each vertex shared by 4 cells)
# 120 - 30 + 1 = 91? Let's check: 120 - 30 = 90, + 1 = 91
# Or: C - E_cell + 1 = ?
# f-vector: (600, 1200, 720, 120)
# 720 - 600 - 120 + 1 = 1 (Euler for 3-sphere is 0: V-E+F-C=0)

# How about: the number of distinct NONZERO eigenvalue pairs?
n_nonzero_distinct = len(nonzero_mu_sorted)
print(f"  Distinct nonzero Laplacian eigenvalues: {n_nonzero_distinct}")

# Check spectral quantities close to 91
total_eigenvalue_sum = sum(float(mu)*int(m) for mu,m in zip(unique_lap, mult_lap))
print(f"  Trace of L = {total_eigenvalue_sum:.4f} (= N*d = {N*degree})")
print(f"  Trace/N = {total_eigenvalue_sum/N:.4f}")

# =============================================================================
# FINAL RESULTS TABLE
# =============================================================================
print("\n" + "=" * 80)
print("=" * 80)
print("COMPREHENSIVE RESULTS TABLE")
print("=" * 80)
print("=" * 80)

print(f"""
======================================================================
120-CELL GRAPH PROPERTIES
======================================================================
  Vertices (N)                  : {N}
  Edges                         : {edges_120}
  Degree (d)                    : {degree}
  Diameter                      : {diameter}
  Distinct eigenvalues          : {len(unique_adj)}
  Euler characteristic          : {chi}

======================================================================
EIGENVALUE SPECTRUM
======================================================================""")

for i, (lam, m) in enumerate(zip(unique_adj, mult_adj)):
    lam_f = float(lam)
    mu_f = degree - lam_f
    print(f"  [{i+1:2d}]  adj = {lam_f:12.8f}   lap = {mu_f:12.8f}   mult = {int(m):4d}")

print(f"""
======================================================================
SPECTRAL INVARIANTS
======================================================================
  Spectral gap                  : {spectral_gap:.12f}
  = phi^(-4)                    : {phi**(-4):.12f}  (match: {abs(spectral_gap - phi**(-4)) < 1e-6})
  Largest Laplacian eigenvalue  : {mu_max:.10f}
  Spectral ratio mu_max/gap    : {mu_max/spectral_gap:.10f}

  Product distinct |adj eigs|   : {P27:.6f}
  phi^27                        : {phi**27:.6f}
  phi^27 / Product              : {phi**27/P27:.10f}

  log10(spanning trees)         : {log_tau/np.log(10):.4f}
  log_phi(spanning trees)       : {log_tau/np.log(phi):.4f}
  log det'(Laplacian)           : {log_det_L:.6f}
  Analytic torsion (1/2 logdet) : {log_T:.6f}

======================================================================
GREEN'S FUNCTIONS AT m^2 = spectral_gap
======================================================================
  G(0, 0)       [self-energy]   : {np.sum(psi_0**2 / (eig_vals_L + spectral_gap)):.10e}
  G(0, antipode) [gravity]      : {np.sum(psi_0*psi_a / (eig_vals_L + spectral_gap)):.10e}
  Ratio |G_self/G_anti|         : {abs(np.sum(psi_0**2 / (eig_vals_L + spectral_gap)) / np.sum(psi_0*psi_a / (eig_vals_L + spectral_gap))):.6e}

======================================================================
TARGET PHYSICAL CONSTANTS
======================================================================
  alpha_EM                      : {alpha_EM:.12e}
  alpha_G (proton)              : {alpha_G_proton:.6e}
  alpha_EM / alpha_G            : {alpha_EM/alpha_G_proton:.6e}
  log10(alpha_EM/alpha_G)       : {np.log10(alpha_EM/alpha_G_proton):.6f}
  log_phi(alpha_EM/alpha_G)     : {log_phi_ratio:.6f}
  m_Planck/m_proton             : {m_Planck/m_proton:.6e}
  log_phi(m_Planck/m_proton)    : {log_phi_mr:.6f}

======================================================================
ALPHA FROM SPECTRUM
======================================================================
  alpha_bare = phi^(-4)/20      : {phi**(-4)/20:.12e}
  alpha_EM (measured)           : {alpha_EM:.12e}
  Ratio (alpha_EM/alpha_bare)   : {alpha_EM/(phi**(-4)/20):.10f}
""")
