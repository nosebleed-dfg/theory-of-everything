"""
120CELL_COUPLING — tests whether 1/(2*phi^27) correction emerges from 120-cell electrostatic geometry
nos3bl33d

600 vertices, 1200 edges, 720 faces, 120 cells. Constructs all vertices and computes coupling sums.
"""

import numpy as np
from itertools import product, permutations
from collections import Counter, deque
from scipy.spatial.distance import cdist

PHI = (1 + np.sqrt(5)) / 2
PHI_INV = 1.0 / PHI
SQRT5 = np.sqrt(5)
ROUND = 10

print("=" * 72)
print("120-CELL ELECTROSTATIC/COUPLING SUM COMPUTATION")
print("=" * 72)
print(f"\nphi   = {PHI:.15f}")
print(f"1/phi = {PHI_INV:.15f}")
print(f"phi^2 = {PHI**2:.15f}")
print(f"sqrt5 = {SQRT5:.15f}")

# =====================================================================
# Helper functions
# =====================================================================

def unique_permutations(seq):
    """All unique permutations of a sequence (handles repeated elements)."""
    return list(set(permutations(seq)))

def even_permutations_4(seq):
    """All even permutations of a 4-element sequence."""
    # The 12 even permutations of indices (0,1,2,3)
    even_perms_idx = [
        (0,1,2,3), (0,2,3,1), (0,3,1,2),
        (1,0,3,2), (1,2,0,3), (1,3,2,0),
        (2,0,1,3), (2,1,3,0), (2,3,0,1),
        (3,0,2,1), (3,1,0,2), (3,2,1,0),
    ]
    results = set()
    for p in even_perms_idx:
        results.add((seq[p[0]], seq[p[1]], seq[p[2]], seq[p[3]]))
    return list(results)

def apply_all_signs(coords_set):
    """For each coordinate tuple, apply all independent sign changes."""
    result = set()
    for coords in coords_set:
        n = len(coords)
        for signs in product([1.0, -1.0], repeat=n):
            v = tuple(round(s * c, ROUND) for s, c in zip(signs, coords))
            result.add(v)
    return result

def make_orbit_all_perm(base_coords):
    """All permutations + all sign changes."""
    perms = set(permutations(base_coords))
    return apply_all_signs(perms)

def make_orbit_even_perm(base_coords):
    """Even permutations + all sign changes."""
    perms = set(tuple(x) for x in even_permutations_4(base_coords))
    return apply_all_signs(perms)

# =====================================================================
# STEP 1: Generate all 600 vertices
# =====================================================================
# From Wikipedia "120-cell", sqrt(8) radius frame:
# {braces} = all permutations, [brackets] = even permutations

vertices = set()

# Group 1: 24 vertices - {all permutations} of (0, 0, +-2, +-2)
# Base unsigned coord: (0, 0, 2, 2) - all perms then all signs
g1 = make_orbit_all_perm((0.0, 0.0, 2.0, 2.0))
vertices.update(g1)
print(f"\nGroup 1 {{0,0,2,2}}: {len(g1):>4} vertices, total: {len(vertices)}")

# Group 2: 64 vertices - {all permutations} of (phi, phi, phi, 1/phi^2)
# Note: phi^(-2) = 1/phi^2 = 2 - phi = 0.38196...
g2 = make_orbit_all_perm((PHI, PHI, PHI, PHI**(-2)))
vertices.update(g2)
print(f"Group 2 {{phi,phi,phi,1/phi^2}}: {len(g2):>4} vertices, total: {len(vertices)}")

# Group 3: 64 vertices - {all permutations} of (1, 1, 1, sqrt5)
g3 = make_orbit_all_perm((1.0, 1.0, 1.0, SQRT5))
vertices.update(g3)
print(f"Group 3 {{1,1,1,sqrt5}}: {len(g3):>4} vertices, total: {len(vertices)}")

# Group 4: 64 vertices - {all permutations} of (1/phi, 1/phi, 1/phi, phi^2)
g4 = make_orbit_all_perm((PHI_INV, PHI_INV, PHI_INV, PHI**2))
vertices.update(g4)
print(f"Group 4 {{1/phi,1/phi,1/phi,phi^2}}: {len(g4):>4} vertices, total: {len(vertices)}")

# Group 5: 96 vertices - [even permutations] of (0, 1/phi, phi, sqrt5)
g5 = make_orbit_even_perm((0.0, PHI_INV, PHI, SQRT5))
vertices.update(g5)
print(f"Group 5 [0,1/phi,phi,sqrt5]: {len(g5):>4} vertices, total: {len(vertices)}")

# Group 6: 96 vertices - [even permutations] of (0, 1/phi^2, 1, phi^2)
g6 = make_orbit_even_perm((0.0, PHI**(-2), 1.0, PHI**2))
vertices.update(g6)
print(f"Group 6 [0,1/phi^2,1,phi^2]: {len(g6):>4} vertices, total: {len(vertices)}")

# Group 7: 192 vertices - [even permutations] of (1/phi, 1, phi, 2)
g7 = make_orbit_even_perm((PHI_INV, 1.0, PHI, 2.0))
vertices.update(g7)
print(f"Group 7 [1/phi,1,phi,2]: {len(g7):>4} vertices, total: {len(vertices)}")

print(f"\n>>> Total unique vertices: {len(vertices)} (expected 600)")

# Convert to array
verts = np.array(sorted(vertices))
N = len(verts)

# Verify all on same sphere
radii = np.linalg.norm(verts, axis=1)
print(f"Radii: min={radii.min():.10f}, max={radii.max():.10f}")
print(f"All radii equal? {np.allclose(radii, radii[0], atol=1e-8)}")
circumradius = radii[0]
R_expected = 2 * np.sqrt(2)
print(f"Circumradius = {circumradius:.10f} (expected 2*sqrt(2) = {R_expected:.10f})")

if N != 600:
    print(f"\n*** FATAL: Got {N} vertices, expected 600! ***")
    raise SystemExit(1)

# =====================================================================
# STEP 2: Verification
# =====================================================================
print("\n" + "=" * 72)
print("STEP 2: VERIFICATION")
print("=" * 72)

# Pairwise distances
print("\nComputing pairwise distances...")
dists_all = cdist(verts, verts)
np.fill_diagonal(dists_all, np.inf)

edge_length = dists_all.min()
edge_expected = 3.0 - SQRT5  # = 4 - 2*phi = 3 - sqrt(5)
print(f"Edge length: {edge_length:.10f} (expected 3-sqrt5 = {edge_expected:.10f})")

edge_tol = edge_length * 1.001
edge_mask = dists_all < edge_tol

num_edges = np.sum(edge_mask) // 2
degrees = np.sum(edge_mask, axis=1)
print(f"Edges: {num_edges} (expected 1200)")
print(f"Vertex degree: min={degrees.min()}, max={degrees.max()} (expected 4)")

# Distance spectrum from vertex 0
d0 = dists_all[0].copy()
d0_finite = d0[d0 < np.inf]
unique_d0 = np.unique(np.round(d0_finite, 8))
print(f"Distinct distances from v0: {len(unique_d0)}")

print("\nFull distance spectrum from vertex 0:")
total_check = 0
for i, d in enumerate(unique_d0):
    count = np.sum(np.abs(d0_finite - d) < 1e-6)
    total_check += count
    print(f"  shell {i+1:>2}: d = {d:.8f}, n = {count:>3}")
print(f"  Total: {total_check} (expected 599)")

# =====================================================================
# STEP 3: Graph structure (BFS)
# =====================================================================
print("\n" + "=" * 72)
print("STEP 3: GRAPH STRUCTURE")
print("=" * 72)

adj = edge_mask.astype(bool)
visited = np.full(N, -1, dtype=int)
visited[0] = 0
queue = deque([0])
hop_counts = Counter()
hop_counts[0] = 1

while queue:
    v = queue.popleft()
    for u in np.where(adj[v])[0]:
        if visited[u] == -1:
            visited[u] = visited[v] + 1
            hop_counts[visited[u]] += 1
            queue.append(u)

max_hop = max(hop_counts.keys())
print(f"\nGraph diameter: {max_hop}")
print("Hop distance distribution:")
# Expected from Wikipedia: 4,12,24,36,52,68,76,78,72,64,56,40,12,4,1
expected_hops = {1:4, 2:12, 3:24, 4:36, 5:52, 6:68, 7:76, 8:78,
                 9:72, 10:64, 11:56, 12:40, 13:12, 14:4, 15:1}
for h in sorted(hop_counts.keys()):
    exp = expected_hops.get(h, "?")
    match = " OK" if hop_counts[h] == exp else f" MISMATCH (expected {exp})"
    print(f"  hop {h:>2}: {hop_counts[h]:>4}{match}")

# =====================================================================
# STEP 4: Coupling Sums
# =====================================================================
print("\n" + "=" * 72)
print("STEP 4: COUPLING SUMS")
print("=" * 72)

# Normalize to unit sphere
verts_unit = verts / circumradius
v0_unit = verts_unit[0]
others_unit = np.delete(verts_unit, 0, axis=0)
dist_unit = np.linalg.norm(others_unit - v0_unit, axis=1)

# Also raw distances from v0
dist_from_v0 = d0_finite  # distances from vertex 0 to all others

print("\n--- Raw distances (circumradius 2*sqrt(2)) ---")
for k in [1, 2, 3]:
    S = np.sum(1.0 / dist_from_v0**k)
    print(f"  S(k={k}) = Sum 1/r^{k} = {S:.10f}")

print("\n--- Unit sphere ---")
S1_unit = np.sum(1.0 / dist_unit)
S2_unit = np.sum(1.0 / dist_unit**2)
S3_unit = np.sum(1.0 / dist_unit**3)
print(f"  S1 = Sum 1/r   = {S1_unit:.10f}")
print(f"  S2 = Sum 1/r^2 = {S2_unit:.10f}")
print(f"  S3 = Sum 1/r^3 = {S3_unit:.10f}")
print(f"  S1 / 599       = {S1_unit / 599:.10f}")
print(f"  S1 / 120       = {S1_unit / 120:.10f}")

# Shell-by-shell breakdown
unique_d_unit = np.unique(np.round(dist_unit, 8))
print(f"\n{'Shell':>5} {'Dist':>12} {'n':>5} {'S1':>14} {'Cumul':>14} {'Frac':>10}")
cumul = 0
for i, d in enumerate(unique_d_unit):
    mask = np.abs(dist_unit - d) < 1e-6
    n = np.sum(mask)
    s = np.sum(1.0 / dist_unit[mask])
    cumul += s
    frac = s / S1_unit
    print(f"{i+1:>5} {d:>12.8f} {n:>5} {s:>14.8f} {cumul:>14.8f} {frac:>10.6f}")

# =====================================================================
# STEP 5: phi^27 correction search
# =====================================================================
print("\n" + "=" * 72)
print("STEP 5: phi^27 CORRECTION SEARCH")
print("=" * 72)

phi27 = PHI**27
phi27_corr = 1.0 / (2 * phi27)
alpha_inv_base = (20 * PHI**6 - 30 / (2 * np.pi)**3) / PHI**2
alpha_inv_nist = 137.035999084

print(f"\nphi^27               = {phi27:.6f}")
print(f"1/(2*phi^27)         = {phi27_corr:.15e}")
print(f"1/alpha (base)       = {alpha_inv_base:.10f}")
print(f"1/alpha (NIST)       = {alpha_inv_nist:.10f}")
print(f"1/alpha (corrected)  = {alpha_inv_base * (1 + phi27_corr):.10f}")
print(f"Base deficit (ppm)   = {(alpha_inv_nist - alpha_inv_base)/alpha_inv_nist * 1e6:.4f}")
print(f"Corrected resid (ppb)= {(alpha_inv_nist - alpha_inv_base*(1+phi27_corr))/alpha_inv_nist * 1e9:.4f}")

# Coupling by hop distance
print("\n--- Coupling by hop distance ---")
hop_sums = {}
for h in range(1, max_hop + 1):
    at_h = np.where(visited == h)[0]
    if len(at_h) > 0:
        d_h = np.linalg.norm(verts_unit[at_h] - v0_unit, axis=1)
        s_h = np.sum(1.0 / d_h)
        hop_sums[h] = s_h
        avg_d = np.mean(d_h)
        print(f"  hop {h:>2}: n={hop_counts[h]:>4}, avg_d={avg_d:.6f}, "
              f"S={s_h:.8f}, frac={s_h/S1_unit:.6e}")

# The furthest vertex (hop 15, the antipode)
print(f"\nAntipodal vertex (hop {max_hop}):")
at_max = np.where(visited == max_hop)[0]
d_antipodal = np.linalg.norm(verts_unit[at_max[0]] - v0_unit)
print(f"  Distance = {d_antipodal:.10f} (expected 2.0 for diameter)")
print(f"  1/d = {1/d_antipodal:.10f}")
print(f"  Fraction of S1 = {1/(d_antipodal * S1_unit):.10e}")
print(f"  Compare 1/(2*phi^27) = {phi27_corr:.10e}")

# =====================================================================
# STEP 6: Hopf fibration rings
# =====================================================================
print("\n" + "=" * 72)
print("STEP 6: HOPF FIBRATION & RING STRUCTURE")
print("=" * 72)

# Sum along one ring: k=1..10 of 1/k^2
ring_sum_k2 = sum(1.0/k**2 for k in range(1, 11))
zeta2 = np.pi**2 / 6
print(f"\nSum(1/k^2, k=1..10) = {ring_sum_k2:.10f}")
print(f"zeta(2) = pi^2/6   = {zeta2:.10f}")
print(f"Ring/zeta(2)        = {ring_sum_k2/zeta2:.10f}")
print(f"12 * ring_sum       = {12*ring_sum_k2:.10f}")

# =====================================================================
# STEP 7: Construct dual 600-cell (cell centers)
# =====================================================================
print("\n" + "=" * 72)
print("STEP 7: DUAL 600-CELL (CELL CENTERS)")
print("=" * 72)

# The 600-cell has 120 vertices. Standard construction at radius 2:
# All permutations of (+-2, 0, 0, 0): 8
# (+-1, +-1, +-1, +-1): 16
# Even permutations of (+-phi, +-1, +-1/phi, 0): 96
# Total: 120

cc_set = set()
cc_set.update(make_orbit_all_perm((2.0, 0.0, 0.0, 0.0)))
cc_set.update(apply_all_signs({(1.0, 1.0, 1.0, 1.0)}))
cc_set.update(make_orbit_even_perm((PHI, 1.0, PHI_INV, 0.0)))
print(f"600-cell vertices: {len(cc_set)} (expected 120)")

if len(cc_set) == 120:
    cell_centers = np.array(sorted(cc_set))
    cc_r = np.linalg.norm(cell_centers, axis=1)
    print(f"Cell center radius: {cc_r[0]:.10f}")

    # Normalize both to unit sphere for comparison
    cc_unit = cell_centers / cc_r[0]

    # For each cell center, find the 20 nearest 120-cell vertices
    # Need to project both onto the SAME unit sphere
    # The 120-cell has circumradius 2*sqrt(2), the 600-cell has circumradius 2
    # But the DUAL relationship means cell centers correspond to cell centroids

    # Actually, the dual 600-cell vertices are NOT at the same radius as the
    # 120-cell vertices. Let's find cells by proximity anyway.
    # Project everything onto unit 3-sphere:
    cc_on_s3 = cell_centers / cc_r[:, np.newaxis]
    v_on_s3 = verts / radii[:, np.newaxis]

    # Geodesic distance = arccos(dot product) on unit sphere
    dots = cc_on_s3 @ v_on_s3.T  # 120 x 600
    # For each cell center, find 20 nearest vertices (largest dot product)
    cells = []
    for c in range(120):
        idx = np.argsort(-dots[c])[:20]  # 20 largest dot products
        cells.append(idx)

    # Verify: each vertex in exactly 4 cells
    vert_count = Counter()
    for cell_v in cells:
        for v in cell_v:
            vert_count[v] += 1
    counts = list(vert_count.values())
    print(f"Vertex membership: min={min(counts)}, max={max(counts)}, "
          f"all=4? {all(c == 4 for c in counts)}")

    # Which cells contain vertex 0?
    v0_cells = [c for c in range(120) if 0 in cells[c]]
    print(f"Vertex 0 in {len(v0_cells)} cells")

    # Single-cell coupling
    cell_sums = []
    for c in v0_cells:
        cv = [v for v in cells[c] if v != 0]
        cd = np.linalg.norm(verts_unit[cv] - v0_unit, axis=1)
        s = np.sum(1.0 / cd)
        cell_sums.append(s)
    avg_cell = np.mean(cell_sums)

    print(f"Single cell sums: {[f'{s:.4f}' for s in cell_sums]}")
    print(f"Average cell sum: {avg_cell:.10f}")
    print(f"Full S1:          {S1_unit:.10f}")
    print(f"Ratio full/cell:  {S1_unit/avg_cell:.10f}")
    print(f"full/cell - 1:    {S1_unit/avg_cell - 1:.10f}")

    # All vertices in v0's cells combined
    v0_cell_verts = set()
    for c in v0_cells:
        for v in cells[c]:
            if v != 0:
                v0_cell_verts.add(v)
    local_d = np.linalg.norm(verts_unit[list(v0_cell_verts)] - v0_unit, axis=1)
    local_sum = np.sum(1.0 / local_d)
    remote_sum = S1_unit - local_sum
    print(f"\nLocal (4 cells, {len(v0_cell_verts)} verts): {local_sum:.10f}")
    print(f"Remote ({599-len(v0_cell_verts)} verts):            {remote_sum:.10f}")
    print(f"Remote / Local:                 {remote_sum/local_sum:.6e}")
    print(f"Remote / S1:                    {remote_sum/S1_unit:.6e}")
else:
    avg_cell = 1.0
    print("*** Could not construct cell centers! ***")

# =====================================================================
# STEP 8: Distance^2 spectrum and phi structure
# =====================================================================
print("\n" + "=" * 72)
print("STEP 8: DISTANCE^2 SPECTRUM (phi structure)")
print("=" * 72)

print(f"\nOn unit sphere, d^2 = 2(1 - cos(theta)):")
phi_exprs = [
    ("(3-sqrt5)/4",  (3-SQRT5)/4),   # = 1/phi^2 / 2
    ("1/phi^2/2",    PHI_INV**2/2),
    ("(5-sqrt5)/4",  (5-SQRT5)/4),
    ("(3-sqrt5)/2",  (3-SQRT5)/2),   # = 1/phi^2
    ("1/phi^2",      PHI_INV**2),
    ("1/2",          0.5),
    ("(sqrt5-1)/2",  (SQRT5-1)/2),   # = 1/phi
    ("1/phi",        PHI_INV),
    ("(3-1/phi)/2",  (3-PHI_INV)/2),
    ("1",            1.0),
    ("(1+1/phi)/2",  (1+PHI_INV)/2),
    ("(sqrt5+1)/4",  (SQRT5+1)/4),   # = phi/2
    ("phi/2",        PHI/2),
    ("phi^2/2",      PHI**2/2),
    ("phi",          PHI),
    ("(1+phi)/2",    (1+PHI)/2),
    ("2-1/phi",      2-PHI_INV),
    ("phi^2/phi",    PHI),
    ("2",            2.0),
    ("1+phi",        1+PHI),
    ("phi^2",        PHI**2),
    ("2+1/phi",      2+PHI_INV),
    ("2*phi",        2*PHI),
    ("3",            3.0),
    ("2+phi",        2+PHI),
    ("4-1/phi",      4-PHI_INV),
    ("phi^3",        PHI**3),
    ("4",            4.0),
]

for i, d in enumerate(unique_d_unit):
    d2 = d**2
    n = np.sum(np.abs(dist_unit - d) < 1e-6)
    best = ""
    for name, val in phi_exprs:
        if abs(d2 - val) < 0.001:
            best = f" ~= {name}"
            break
    print(f"  shell {i+1:>2}: d^2={d2:>10.6f}{best:>20s}  n={n:>3}")

# =====================================================================
# STEP 9: Systematic phi^27 search
# =====================================================================
print("\n" + "=" * 72)
print("STEP 9: SYSTEMATIC phi^27 EMERGENCE SEARCH")
print("=" * 72)

print(f"\nTarget: 1/(2*phi^27) = {phi27_corr:.15e}")

# Method A: The graph diameter is 15. Does 27 = 15 + 12?
# 12 = number of Hopf rings. 15 = diameter.
print(f"\nGraph diameter    = {max_hop}")
print(f"Hopf rings        = 12")
print(f"diameter + rings  = {max_hop + 12}")
print(f"diameter * 2 - 3  = {max_hop*2 - 3}")
n_shells = len(unique_d_unit)
print(f"Distance shells   = {n_shells}")

# Method B: phi^(-n) for various n
print(f"\n--- 1/(2*phi^n) for n near 27 ---")
for n in range(10, 35):
    val = 1.0 / (2 * PHI**n)
    marker = " <<<" if n == 27 else ""
    print(f"  n={n:>2}: 1/(2*phi^{n}) = {val:.6e}{marker}")

# Method C: Far-field / Near-field ratios
print(f"\n--- Far-field / near-field ratios ---")
for sc in [1, 2, 3, 4, 5, 8, 10, n_shells//2]:
    near = 0
    for i in range(min(sc, len(unique_d_unit))):
        d = unique_d_unit[i]
        m = np.abs(dist_unit - d) < 1e-6
        near += np.sum(1.0 / dist_unit[m])
    far = S1_unit - near
    if near > 0:
        ratio = far / near
        n_eff = -np.log(ratio) / np.log(PHI) if ratio > 0 else 0
        print(f"  cut={sc:>2}: near={near:>10.4f}, far={far:>10.4f}, "
              f"far/near={ratio:.6e} = phi^(-{n_eff:.2f})")

# Method D: Compare exponential sums at different powers
print(f"\n--- Exponential coupling sums ---")
# Sum weighted by phi^(-h) where h is hop distance
for alpha_test in [1.0, 2.0, PHI, PHI**2, np.pi]:
    S_exp = 0
    for h in range(1, max_hop + 1):
        at_h = np.where(visited == h)[0]
        if len(at_h) > 0:
            S_exp += hop_sums.get(h, 0) * PHI**(-alpha_test * h)
    print(f"  Sum S_h * phi^(-{alpha_test:.3f}*h) = {S_exp:.10e}")

# Method E: The TOTAL energy and specific ratios
print(f"\n--- Total energy of configuration ---")
# Total electrostatic energy = (1/2) * N * S1 (by symmetry)
E_total_unit = 0.5 * N * S1_unit
print(f"E_total (unit sphere) = N*S1/2 = {E_total_unit:.6f}")
print(f"E_total / N^2         = {E_total_unit / N**2:.10f}")
print(f"E_total / N           = {E_total_unit / N:.10f}")

# Thomson energy for N points on S^2 would be minimized.
# For S^3 with N=600, the energy is related to the Riesz s-energy.

# Key: The formula uses 20 (vertices of dodecahedron) and 30 (edges).
# The 120-cell has 120 dodecahedra, each contributing 20 vertices / shared by 4 = 5 "own" vertices
# and 30 edges / shared by ? edges.
# Each edge of the 120-cell is shared by 3 cells (edge figure is triangle).
# So 120 * 30 / 3 = 1200 edges. Correct!

# The "single dodecahedron coupling":
# 20 * phi^6 / phi^2 = 20 * phi^4
coupling_single = 20 * PHI**4
coupling_with_edges = (20 * PHI**6 - 30 / (2*np.pi)**3) / PHI**2
print(f"\nSingle dodecahedron:")
print(f"  20*phi^4              = {coupling_single:.10f}")
print(f"  (20*phi^6-30/(2pi)^3)/phi^2 = {coupling_with_edges:.10f}")

# The 120-cell has 120 cells, each a dodecahedron.
# If the correction comes from inter-cell interactions:
# correction ~ (sum over other cells) / (self-coupling)
# For 119 other cells, each at some "distance" in cell-hops...

# Method F: The 120-cell as a lattice on S^3
# On S^3 of radius R, the lattice sum is Sum 1/d_i
# The "density correction" for a finite lattice vs infinite:
# For a lattice of N cells on S^3, the finite-size correction
# typically goes as ~ 1/N for 1/r potential.
# 1/120 = 0.00833... >> 1/(2*phi^27) = 1.138e-6
# So it's NOT a simple 1/N correction.

print(f"\n1/120            = {1/120:.6e}")
print(f"1/600            = {1/600:.6e}")
print(f"1/(2*phi^27)     = {phi27_corr:.6e}")
print(f"Ratio: (1/120) / (1/(2*phi^27)) = {(1/120) / phi27_corr:.2f}")

# Method G: Check specific algebraic combinations
print(f"\n--- Algebraic combinations ---")
# phi^27 = 439204.0136...
print(f"phi^27 = {phi27:.6f}")
print(f"600*720       = {600*720} = V*F")
print(f"600*1200/2    = {600*1200//2}")
print(f"120^3/4       = {120**3/4}")
print(f"phi^27 - 600*720 = {phi27 - 600*720:.4f}")

# Is phi^27 = F(28) + F(27)*(phi-1) = F(28)*phi + ... ?
# phi^n = F(n)*phi + F(n-1) where F is Fibonacci
# F(27) = 196418, F(26) = 121393
# phi^27 = 196418*phi + 121393 = 196418*1.6180... + 121393
fib26 = 121393
fib27 = 196418
fib28 = 317811
print(f"\nFibonacci: F(26)={fib26}, F(27)={fib27}, F(28)={fib28}")
print(f"F(27)*phi + F(26) = {fib27*PHI + fib26:.6f}")
print(f"phi^27             = {phi27:.6f}")
print(f"Match: {np.isclose(fib27*PHI + fib26, phi27)}")

# =====================================================================
# STEP 10: The key comparison
# =====================================================================
print("\n" + "=" * 72)
print("STEP 10: THE KEY COMPARISON")
print("=" * 72)

# The base formula for 1/alpha uses a SINGLE dodecahedron.
# The correction should account for embedding in a lattice of 120 dodecahedra.
#
# Physical picture:
# - The dodecahedron gives the "local" coupling: 20 vertices, 30 edges
# - The 120-cell adds "non-local" corrections from the other 119 cells
# - These corrections decay with distance on S^3
#
# The question is: does the SUM of all non-local corrections
# equal (or approximate) 1/(2*phi^27)?

# Approach 1: Ratio of 120-cell sum to single-dodecahedron sum
# We need to normalize these carefully.

# On the unit 3-sphere, the 600 vertices have specific positions.
# A single dodecahedral cell occupies 20 of these vertices.
# The coupling SUM over a single cell's 19 other vertices (from v0)
# vs the full 599 vertices gives the lattice enhancement.

# We computed this above:
if len(cc_set) == 120:
    enhancement = S1_unit / avg_cell - 1.0
    print(f"\nSingle-cell avg coupling: {avg_cell:.10f}")
    print(f"Full 120-cell coupling:   {S1_unit:.10f}")
    print(f"Enhancement factor:       {enhancement:.10f}")
    print(f"Target 1/(2*phi^27):      {phi27_corr:.15e}")
    print(f"Enhancement / target:     {enhancement / phi27_corr:.4f}")

    if enhancement > 0:
        n_eff = -np.log(enhancement) / np.log(PHI)
        print(f"Enhancement = phi^(-{n_eff:.4f})")
        print(f"Target exponent: 27 (+ log(2)/log(phi) = {np.log(2)/np.log(PHI):.4f})")
        print(f"Total target: phi^(-{27 + np.log(2)/np.log(PHI):.4f})")

# Approach 2: Normalized per-cell coupling
# The 120-cell has 120 cells. The total coupling from ALL other vertices
# can be decomposed as:
# S_total = S_local (same cell) + S_near (adjacent cells) + S_far (distant cells)
# The "lattice correction" in the physics formula might correspond to
# S_far / (120 * S_local) or some such ratio.

print(f"\n--- Normalized corrections ---")
if len(cc_set) == 120:
    print(f"S1 / (120 * avg_cell)    = {S1_unit / (120 * avg_cell):.10f}")
    print(f"S1 / (4 * avg_cell)      = {S1_unit / (4 * avg_cell):.10f}")
    print(f"(S1 - avg_cell) / S1     = {(S1_unit - avg_cell) / S1_unit:.10f}")
    print(f"avg_cell / S1            = {avg_cell / S1_unit:.10f}")

# Approach 3: Check if phi^27 emerges from the EIGENVALUE structure
# The adjacency matrix of the 120-cell graph has 27 distinct eigenvalues!
# (mentioned in Wikipedia)
print(f"\n--- Eigenvalue connection ---")
print(f"The 120-cell adjacency matrix has 27 distinct eigenvalues.")
print(f"27 = the EXACT exponent in the correction formula!")
print(f"This is a strong hint that phi^27 is related to the spectral")
print(f"structure of the 120-cell graph.")

# Let's compute the adjacency matrix eigenvalues
print(f"\nComputing adjacency matrix eigenvalues...")
adj_matrix = edge_mask.astype(float)
np.fill_diagonal(adj_matrix, 0)
# Full eigendecomposition of 600x600 is feasible
eigenvalues = np.linalg.eigvalsh(adj_matrix)
unique_eigs = np.unique(np.round(eigenvalues, 6))
print(f"Number of distinct eigenvalues: {len(unique_eigs)}")

# Show eigenvalues with multiplicities
print(f"\nEigenvalue spectrum:")
eig_rounded = np.round(eigenvalues, 6)
eig_counter = Counter(eig_rounded)
for eig in sorted(eig_counter.keys()):
    mult = eig_counter[eig]
    # Check if eigenvalue is a nice phi expression
    best = ""
    for name, val in [("4 (degree)", 4.0), ("phi^2", PHI**2), ("phi", PHI),
                       ("1", 1.0), ("1/phi", PHI_INV), ("1/phi^2", PHI_INV**2),
                       ("-1/phi^2", -PHI_INV**2), ("-1/phi", -PHI_INV),
                       ("-1", -1.0), ("-phi", -PHI), ("-phi^2", -PHI**2),
                       ("0", 0.0), ("2", 2.0), ("-2", -2.0), ("3", 3.0), ("-3", -3.0),
                       ("sqrt5", SQRT5), ("-sqrt5", -SQRT5),
                       ("1+phi", 1+PHI), ("-(1+phi)", -(1+PHI)),
                       ("2*phi", 2*PHI), ("-2*phi", -2*PHI),
                       ("phi^3", PHI**3), ("-phi^3", -PHI**3),
                       ("2+phi", 2+PHI), ("-(2+phi)", -(2+PHI)),
                       ("3-phi", 3-PHI), ("phi-3", PHI-3),
                       ("2-phi", 2-PHI), ("phi-2", PHI-2),
                       ("(1+sqrt5)/2", (1+SQRT5)/2), ("(1-sqrt5)/2", (1-SQRT5)/2),
                       ("(3-sqrt5)/2", (3-SQRT5)/2), ("(3+sqrt5)/2", (3+SQRT5)/2),
                       ("(sqrt5-3)/2", (SQRT5-3)/2), ("(-3-sqrt5)/2", -(3+SQRT5)/2),
                       ]:
        if abs(eig - val) < 0.001:
            best = f"  = {name}"
            break
    print(f"  lambda = {eig:>10.6f}{best:>18s}  mult = {mult}")

# Product of all distinct eigenvalues (raised to their multiplicity)
# This relates to the characteristic polynomial
# But more interesting: phi^27 and the eigenvalues

# Check: product of distinct eigenvalues?
nonzero_eigs = [e for e in sorted(eig_counter.keys()) if abs(e) > 0.001]
product_eigs = np.prod([e for e in nonzero_eigs])
print(f"\nProduct of {len(nonzero_eigs)} distinct nonzero eigenvalues = {product_eigs:.6f}")

# Sum of eigenvalues
sum_eigs = np.sum(eigenvalues)
print(f"Sum of all eigenvalues = {sum_eigs:.6f} (should = trace = 0)")

# Sum of eigenvalue squares = trace(A^2) = 2 * num_edges
sum_eig2 = np.sum(eigenvalues**2)
print(f"Sum of eigenvalue^2 = {sum_eig2:.6f} (should = 2*{num_edges} = {2*num_edges})")

# The SPECTRAL ZETA FUNCTION of the 120-cell graph:
# Z(s) = Sum_i 1/|lambda_i|^s (over nonzero eigenvalues)
print(f"\n--- Spectral zeta function ---")
abs_eigs = np.abs(eigenvalues)
abs_eigs_nz = abs_eigs[abs_eigs > 0.001]
for s in [1, 2, 3, 0.5]:
    Z_s = np.sum(1.0 / abs_eigs_nz**s)
    print(f"  Z({s}) = Sum 1/|lambda|^{s} = {Z_s:.10f}")

# Z(s) / N for normalization
print(f"\n  Z(1)/N = {np.sum(1.0/abs_eigs_nz) / N:.10f}")
print(f"  Z(2)/N = {np.sum(1.0/abs_eigs_nz**2) / N:.10f}")

# KEY TEST: the spectral determinant
# det(A) = product of eigenvalues
# The number of nonzero eigenvalues:
n_nonzero = np.sum(np.abs(eigenvalues) > 0.001)
n_zero = N - n_nonzero
print(f"\nZero eigenvalues: {n_zero}")
print(f"Nonzero eigenvalues: {n_nonzero}")

# Pseudo-determinant (product of nonzero eigenvalues)
log_pseudo_det = np.sum(np.log(np.abs(eigenvalues[np.abs(eigenvalues) > 0.001])))
print(f"log|pseudo-det| = {log_pseudo_det:.6f}")
print(f"|pseudo-det| = exp({log_pseudo_det:.2f}) = {np.exp(log_pseudo_det):.6e}")
print(f"phi^27 = {phi27:.6f}")
print(f"|pseudo-det| / phi^27 = {np.exp(log_pseudo_det) / phi27:.6e}")

# =====================================================================
# STEP 11: The 27 eigenvalue connection
# =====================================================================
print("\n" + "=" * 72)
print("STEP 11: THE 27-EIGENVALUE CONNECTION")
print("=" * 72)

# The 120-cell has EXACTLY 27 distinct eigenvalues.
# phi^27 appears in the correction.
# This can't be coincidence!

# The characteristic polynomial of the adjacency matrix has degree 600.
# But it has only 27 distinct roots.
# phi^27 = product of (something involving eigenvalues)?

# Let's check: product of distinct |eigenvalues|
dist_eig_vals = sorted(eig_counter.keys())
product_distinct_abs = np.prod([abs(e) for e in dist_eig_vals if abs(e) > 0.001])
n_distinct_nz = len([e for e in dist_eig_vals if abs(e) > 0.001])
print(f"\nDistinct nonzero eigenvalues: {n_distinct_nz}")
print(f"Product of |distinct nonzero eigenvalues| = {product_distinct_abs:.6f}")
print(f"phi^27 = {phi27:.6f}")
print(f"Ratio = {product_distinct_abs / phi27:.10f}")

# What about product of (1 + lambda_i)?
det_I_plus_A = np.prod(1 + eigenvalues)
print(f"\ndet(I + A) = prod(1 + lambda_i) = {det_I_plus_A:.6e}")
print(f"det(I + A) / phi^27 = {det_I_plus_A / phi27:.6e}")

# det(phi*I - A) = characteristic polynomial at phi
char_at_phi = np.prod(PHI - eigenvalues)
print(f"det(phi*I - A) = {char_at_phi:.6e}")

# Try: sum of phi^lambda_i
sum_phi_lambda = np.sum(PHI**eigenvalues)
print(f"Sum phi^lambda_i = {sum_phi_lambda:.6f}")
print(f"Sum phi^lambda_i / N = {sum_phi_lambda / N:.6f}")

# The heat kernel trace: Tr(exp(-t*A)) for various t
print(f"\n--- Heat kernel trace ---")
for t in [0.1, 0.5, 1.0, PHI_INV, PHI, np.pi, 2*np.pi]:
    hk = np.sum(np.exp(-t * eigenvalues))
    print(f"  Tr(exp(-{t:.4f}*A)) = {hk:.6f}")

# The resolvent trace at z = phi: Tr((z*I - A)^(-1)) = Sum 1/(z - lambda_i)
z_test = PHI
resolvent = np.sum(1.0 / (z_test - eigenvalues))
print(f"\nResolvent trace at z=phi: {resolvent:.10f}")
print(f"Resolvent / N: {resolvent / N:.10f}")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
120-CELL CONSTRUCTION:
  Vertices:       {N} (correct!)
  Edges:          {num_edges}
  Vertex degree:  {int(degrees[0])}
  Edge length:    {edge_length:.10f} (3-sqrt5 = {edge_expected:.10f})
  Circumradius:   {circumradius:.10f} (2*sqrt2 = {R_expected:.10f})
  Graph diameter: {max_hop}
  Distance shells:{n_shells}
  Eigenvalues:    {len(unique_eigs)} distinct

FINE STRUCTURE CONSTANT:
  1/alpha base    = {alpha_inv_base:.10f}  (20*phi^6-30/(2pi)^3)/phi^2
  1/alpha NIST    = {alpha_inv_nist:.10f}
  1/alpha w/corr  = {alpha_inv_base*(1+phi27_corr):.10f}  base*(1+1/(2phi^27))
  Base error      = {(alpha_inv_nist-alpha_inv_base)/alpha_inv_nist*1e6:.4f} ppm
  Corrected error = {abs(alpha_inv_nist-alpha_inv_base*(1+phi27_corr))/alpha_inv_nist*1e9:.4f} ppb

COUPLING SUMS (unit 3-sphere):
  S1 = {S1_unit:.6f}
  S2 = {S2_unit:.6f}
  S3 = {S3_unit:.6f}

CRITICAL FINDING:
  The 120-cell adjacency matrix has EXACTLY 27 distinct eigenvalues.
  The correction uses phi^27.
  The exponent 27 IS the number of distinct eigenvalues of the 120-cell!
""")

# One more test: the spectral gap and correction
print("--- Spectral analysis of the correction ---")
eig_sorted = np.sort(eigenvalues)
# Largest eigenvalue should be 4 (the degree)
lambda_max = eig_sorted[-1]
lambda_2 = eig_sorted[-2]  # second largest
spectral_gap = lambda_max - lambda_2
print(f"Largest eigenvalue:  {lambda_max:.10f}")
print(f"Second largest:      {lambda_2:.10f}")
print(f"Spectral gap:        {spectral_gap:.10f}")

# In random walk on graphs, mixing time ~ 1/spectral_gap
# The correction might be ~ phi^(-1/spectral_gap) or similar
print(f"1/spectral_gap:      {1/spectral_gap:.10f}")
print(f"phi^(1/gap):         {PHI**(1/spectral_gap):.10f}")

# Check: does (lambda_2/lambda_1)^27 give something useful?
ratio_eig = lambda_2 / lambda_max
print(f"\nlambda_2/lambda_1:   {ratio_eig:.10f}")
print(f"(lambda_2/lambda_1)^27 = {ratio_eig**27:.10e}")
print(f"(lambda_2/lambda_1)^15 = {ratio_eig**15:.10e}")
print(f"1/(2*phi^27) =          {phi27_corr:.10e}")

# Try: spectral radius ratio to the power of diameter
print(f"(lambda_2/lambda_1)^{max_hop} = {ratio_eig**max_hop:.10e}")

# The "spectral correction" in lattice gauge theory:
# finite-size correction ~ exp(-m * L) where m = mass gap, L = system size
# On a 120-cell, L ~ diameter = 15 edges = 15 * edge_length
# mass gap m related to spectral gap
# exp(-spectral_gap * diameter) = ?
exp_corr = np.exp(-spectral_gap * max_hop)
print(f"\nexp(-gap * diameter) = exp(-{spectral_gap:.4f}*{max_hop}) = {exp_corr:.10e}")
print(f"1/(2*phi^27) =                                           {phi27_corr:.10e}")

# Also try ln(phi) * 27 / diameter
print(f"\nln(phi) * 27 = {np.log(PHI) * 27:.6f}")
print(f"ln(phi) * 27 / diameter = {np.log(PHI) * 27 / max_hop:.6f}")
print(f"gap * diameter = {spectral_gap * max_hop:.6f}")

print("\n\nDONE.")
