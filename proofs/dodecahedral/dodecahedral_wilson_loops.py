"""
DODECAHEDRAL_WILSON_LOOPS — non-abelian Wilson loop spectrum for the 16-band metallic group on T^11
nos3bl33d

16x16 Berry phase eigenvalue curves: Z_2 invariants, winding numbers, band inversions, A5 irrep mixing.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from numpy.linalg import eigh, svd, det, eigvals
from itertools import combinations
import time

PHI = (1 + np.sqrt(5)) / 2

# ============================================================
# Dodecahedron construction (from existing code)
# ============================================================

def build_dodecahedron():
    adj = [
        [1, 4, 5], [0, 2, 6], [1, 3, 7], [2, 4, 8], [3, 0, 9],
        [0, 14, 10], [1, 10, 11], [2, 11, 12], [3, 12, 13], [4, 13, 14],
        [5, 6, 15], [10, 16, 19], [6, 7, 16], [7, 8, 17], [8, 9, 18],
        [9, 5, 19], [11, 15, 17], [12, 16, 18], [13, 17, 19], [14, 18, 15],
    ]
    edges = []
    for u, nbrs in enumerate(adj):
        for v in nbrs:
            if u < v:
                edges.append((u, v))
    assert len(edges) == 30, f"Expected 30 edges, got {len(edges)}"
    return adj, edges


def find_nontree_edges(adj, edges, nv=20):
    visited = [False]*nv
    tree = set()
    queue = [0]
    visited[0] = True
    while queue:
        u = queue.pop(0)
        for v in adj[u]:
            if not visited[v]:
                visited[v] = True
                tree.add((min(u,v), max(u,v)))
                queue.append(v)
    nontree = [e for e in edges if e not in tree]
    assert len(nontree) == 11, f"Expected 11 non-tree edges, got {len(nontree)}"
    return nontree


def build_H(k, adj, nt_idx, nv=20):
    """Build the 20x20 Bloch Hamiltonian at quasi-momentum k in R^11."""
    H = np.zeros((nv, nv), dtype=complex)
    for u, nbrs in enumerate(adj):
        for v in nbrs:
            if u < v:
                e = (u, v)
                if e in nt_idx:
                    ph = np.exp(1j * k[nt_idx[e]])
                    H[u,v] += ph
                    H[v,u] += np.conj(ph)
                else:
                    H[u,v] += 1.0
                    H[v,u] += 1.0
    return H


# ============================================================
# Wilson loop computation
# ============================================================

def compute_wilson_loop(k_base, dir_loop, dir_perp, k_perp_val,
                        adj, nt_idx, band_indices, N_loop, nv=20):
    """
    Compute the Wilson loop matrix for a set of bands.

    Parameters:
        k_base: base k-vector (11D, zeros typically)
        dir_loop: integer, the BZ direction to integrate around (0-10)
        dir_perp: integer, the BZ direction treated as parameter (0-10)
        k_perp_val: float, the value of k[dir_perp]
        adj: adjacency list
        nt_idx: dict mapping non-tree edges to indices
        band_indices: list of band indices to track (e.g., range(16))
        N_loop: number of discretization points around the loop
        nv: number of vertices (20)

    Returns:
        W: the Wilson loop matrix (m x m unitary matrix)
        unitarity_error: ||W^dag W - I||
    """
    m = len(band_indices)
    dk = 2 * np.pi / N_loop

    # Collect eigenvectors at each point around the loop
    all_vecs = []
    for j in range(N_loop):
        k = k_base.copy()
        k[dir_perp] = k_perp_val
        k[dir_loop] = dk * j
        _, evecs = eigh(build_H(k, adj, nt_idx, nv))
        # Extract the bands we care about: nv x m matrix
        all_vecs.append(evecs[:, band_indices])

    # Build Wilson loop as product of overlap matrices
    # W = M_0 * M_1 * ... * M_{N-1}
    # where M_j = <u(k_j)|u(k_{j+1})> is m x m
    W = np.eye(m, dtype=complex)
    for j in range(N_loop):
        j_next = (j + 1) % N_loop
        # Overlap matrix: m x m
        M = all_vecs[j].conj().T @ all_vecs[j_next]

        # SVD-regularize M to handle gauge ambiguities at degeneracies
        # This replaces M with the closest unitary matrix
        U_svd, S_svd, Vh_svd = svd(M)
        M_reg = U_svd @ Vh_svd  # Closest unitary (polar decomposition)

        W = W @ M_reg

    # Measure unitarity
    unitarity_err = np.linalg.norm(W.conj().T @ W - np.eye(m))

    return W, unitarity_err


def wilson_loop_spectrum(dir_loop, dir_perp, adj, nt_idx,
                         band_indices, N_loop=30, N_perp=60, nv=20):
    """
    Compute the Wilson loop spectrum as a function of k_perp.

    Returns:
        k_perp_vals: array of k_perp values [0, 2*pi)
        theta_spectrum: array of shape (N_perp, m) with Wilson loop
                        eigenvalue phases theta_n in [-pi, pi]
        unitarity_errors: array of unitarity errors
    """
    m = len(band_indices)
    k_base = np.zeros(11)
    dk_perp = 2 * np.pi / N_perp

    k_perp_vals = np.array([dk_perp * i for i in range(N_perp)])
    theta_spectrum = np.zeros((N_perp, m))
    unitarity_errors = np.zeros(N_perp)

    for i in range(N_perp):
        W, err = compute_wilson_loop(
            k_base, dir_loop, dir_perp, k_perp_vals[i],
            adj, nt_idx, band_indices, N_loop, nv
        )
        unitarity_errors[i] = err

        # Eigenvalues of the Wilson loop
        evals = eigvals(W)
        # Extract phases and sort
        phases = np.angle(evals)
        phases.sort()
        theta_spectrum[i] = phases

    return k_perp_vals, theta_spectrum, unitarity_errors


# ============================================================
# Topology extraction from Wilson loop spectrum
# ============================================================

def count_theta_crossings(theta_spectrum, theta_target, tolerance=0.15):
    """
    Count how many times eigenvalue curves cross theta = theta_target.

    Uses a simple crossing detection: look for sign changes in
    (theta - theta_target) between consecutive k_perp points,
    handling the branch cut at +/-pi.
    """
    N_perp, m = theta_spectrum.shape
    crossings = 0

    for band_idx in range(m):
        for i in range(N_perp):
            i_next = (i + 1) % N_perp
            t1 = theta_spectrum[i, band_idx] - theta_target
            t2 = theta_spectrum[i_next, band_idx] - theta_target

            # Wrap to [-pi, pi]
            t1 = (t1 + np.pi) % (2*np.pi) - np.pi
            t2 = (t2 + np.pi) % (2*np.pi) - np.pi

            # Crossing = sign change and both values small
            if abs(t1) < tolerance and abs(t2) < tolerance and t1 * t2 < 0:
                crossings += 1

    return crossings


def compute_winding_numbers(theta_spectrum):
    """
    Compute the total winding number of the Wilson loop eigenvalues.

    The winding number is how many times each eigenvalue curve winds
    around the unit circle as k_perp goes from 0 to 2*pi.
    """
    N_perp, m = theta_spectrum.shape
    windings = np.zeros(m)

    for band_idx in range(m):
        total_wind = 0.0
        for i in range(N_perp):
            i_next = (i + 1) % N_perp
            dtheta = theta_spectrum[i_next, band_idx] - theta_spectrum[i, band_idx]
            # Unwrap: if jump > pi, it wrapped around
            if dtheta > np.pi:
                dtheta -= 2*np.pi
            elif dtheta < -np.pi:
                dtheta += 2*np.pi
            total_wind += dtheta
        windings[band_idx] = total_wind / (2*np.pi)

    return windings


def track_eigenvalue_continuity(theta_spectrum):
    """
    Re-sort eigenvalue curves for continuity across k_perp steps.

    The naive sorting at each k_perp can cause artificial jumps when
    eigenvalue curves cross. This does a nearest-neighbor tracking
    to produce continuous curves.

    Returns:
        tracked: array of shape (N_perp, m) with continuously tracked phases
    """
    N_perp, m = theta_spectrum.shape
    tracked = np.zeros_like(theta_spectrum)
    tracked[0] = theta_spectrum[0]

    for i in range(1, N_perp):
        prev = tracked[i-1]
        curr = theta_spectrum[i].copy()

        # Hungarian-style assignment: match each current eigenvalue
        # to the nearest previous one
        used = [False] * m
        assignment = [-1] * m

        # Distance matrix (accounting for periodic boundary in theta)
        dist = np.zeros((m, m))
        for a in range(m):
            for b in range(m):
                d = abs(curr[b] - prev[a])
                d = min(d, 2*np.pi - d)  # periodic distance
                dist[a, b] = d

        # Greedy assignment (good enough for smooth curves)
        for _ in range(m):
            best_a, best_b = -1, -1
            best_d = np.inf
            for a in range(m):
                if assignment[a] >= 0:
                    continue
                for b in range(m):
                    if used[b]:
                        continue
                    if dist[a, b] < best_d:
                        best_d = dist[a, b]
                        best_a, best_b = a, b
            assignment[best_a] = best_b
            used[best_b] = True

        for a in range(m):
            tracked[i, a] = curr[assignment[a]]

    return tracked


def detect_degeneracies(theta_spectrum, threshold=0.05):
    """
    Find (k_perp_index, band_pair) where Wilson loop eigenvalues
    are nearly degenerate. These are topologically significant.
    """
    N_perp, m = theta_spectrum.shape
    degeneracies = []

    for i in range(N_perp):
        phases = np.sort(theta_spectrum[i])
        for a in range(m):
            for b in range(a+1, m):
                d = abs(phases[b] - phases[a])
                d = min(d, 2*np.pi - d)
                if d < threshold:
                    degeneracies.append((i, a, b, d))

    return degeneracies


# ============================================================
# A5 irrep analysis
# ============================================================

def analyze_a5_clustering(theta_spectrum, k_perp_vals):
    """
    At k=0 the 16 bands decompose as chi_3'(3) + chi_4(4) + chi_4'(4) + chi_5(5).
    Check if the Wilson loop eigenvalues cluster into groups of these sizes.
    """
    N_perp, m = theta_spectrum.shape

    results = []
    for i in range(N_perp):
        phases = np.sort(theta_spectrum[i])

        # Hierarchical clustering: find natural gaps
        gaps = np.diff(phases)
        # Handle wrap-around gap
        wrap_gap = (phases[0] + 2*np.pi) - phases[-1]
        all_gaps = np.append(gaps, wrap_gap)

        # Find the 3 largest gaps (to split into 4 groups if they exist)
        gap_indices = np.argsort(all_gaps)[::-1]

        # Take top 3 gaps
        split_points = sorted(gap_indices[:3])

        # Build groups
        groups = []
        start = 0
        for sp in split_points:
            if sp < m - 1:  # internal gap
                groups.append(list(range(start, sp+1)))
                start = sp + 1
        groups.append(list(range(start, m)))

        group_sizes = sorted([len(g) for g in groups])
        results.append({
            'k_perp': k_perp_vals[i],
            'group_sizes': group_sizes,
            'phases': phases,
            'max_gap': all_gaps.max(),
            'groups': groups
        })

    return results


# ============================================================
# ASCII spectrum visualization
# ============================================================

def ascii_wilson_spectrum(k_perp_vals, theta_spectrum, title="", width=78, height=30):
    """
    Render the Wilson loop spectrum as ASCII art.

    Horizontal axis: k_perp from 0 to 2*pi
    Vertical axis: theta from -pi to pi
    Each eigenvalue curve is plotted as a character
    """
    N_perp, m = theta_spectrum.shape

    # Create grid
    grid = [[' '] * width for _ in range(height)]

    # Map k_perp to x, theta to y
    def to_x(k):
        return int(k / (2*np.pi) * (width - 1))

    def to_y(theta):
        # theta in [-pi, pi] -> row in [0, height-1], with -pi at bottom
        return height - 1 - int((theta + np.pi) / (2*np.pi) * (height - 1))

    # Symbols for different eigenvalue curves
    symbols = '0123456789abcdef'

    # Plot each eigenvalue curve
    for band_idx in range(min(m, 16)):
        sym = symbols[band_idx % len(symbols)]
        for i in range(N_perp):
            x = to_x(k_perp_vals[i])
            y = to_y(theta_spectrum[i, band_idx])
            x = max(0, min(width-1, x))
            y = max(0, min(height-1, y))
            if grid[y][x] == ' ':
                grid[y][x] = sym
            elif grid[y][x] != sym:
                grid[y][x] = '*'  # Overlap/crossing marker

    # Add axis labels
    lines = []
    if title:
        lines.append(title)
    lines.append(f"  +pi {'=' * (width-8)} ")
    for row in range(height):
        if row == 0:
            label = '+pi|'
        elif row == height // 2:
            label = '  0|'
        elif row == height - 1:
            label = '-pi|'
        else:
            label = '   |'
        lines.append(f"{label}{''.join(grid[row])}")
    lines.append(f"     {'=' * (width-4)}")
    lines.append(f"     0{' ' * (width//2 - 4)}pi{' ' * (width//2 - 3)}2*pi")
    lines.append(f"     {'k_perp (direction perpendicular to loop)':^{width-4}}")
    lines.append(f"     Symbols: 0-f = 16 eigenvalue curves, * = crossing")

    return '\n'.join(lines)


def ascii_spectrum_density(k_perp_vals, theta_spectrum, title="", width=78, height=30):
    """
    Alternative visualization: density plot showing how many eigenvalues
    fall in each region. Better for seeing clustering.
    """
    N_perp, m = theta_spectrum.shape
    grid = np.zeros((height, width))

    for band_idx in range(m):
        for i in range(N_perp):
            x = int(k_perp_vals[i] / (2*np.pi) * (width - 1))
            y = height - 1 - int((theta_spectrum[i, band_idx] + np.pi) / (2*np.pi) * (height - 1))
            x = max(0, min(width-1, x))
            y = max(0, min(height-1, y))
            grid[y, x] += 1

    density_chars = ' .,:;+*#@'
    max_val = max(grid.max(), 1)

    lines = []
    if title:
        lines.append(title)
    for row in range(height):
        if row == 0:
            label = '+pi|'
        elif row == height // 2:
            label = '  0|'
        elif row == height - 1:
            label = '-pi|'
        else:
            label = '   |'
        row_str = ''
        for col in range(width):
            idx = int(grid[row, col] / max_val * (len(density_chars) - 1))
            row_str += density_chars[min(idx, len(density_chars)-1)]
        lines.append(f"{label}{row_str}")
    lines.append(f"     {'=' * width}")
    lines.append(f"     0{' ' * (width//2 - 4)}pi{' ' * (width//2 - 3)}2*pi")
    lines.append(f"     Density: {density_chars} (low to high)")

    return '\n'.join(lines)


# ============================================================
# Main computation
# ============================================================

def main():
    print()
    print("*" * 72)
    print("*  DODECAHEDRAL BLOCH HAMILTONIAN: NON-ABELIAN WILSON LOOPS")
    print("*  16-Band Metallic Group Topological Analysis")
    print(f"*  phi = {PHI:.10f}")
    print("*" * 72)

    adj, edges = build_dodecahedron()
    nontree = find_nontree_edges(adj, edges)
    nt_idx = {e: i for i, e in enumerate(nontree)}

    print(f"\nDodecahedron: 20 vertices, 30 edges, b_1 = 11")
    print(f"Non-tree edges: {nontree}")

    # k=0 spectrum
    H0 = build_H(np.zeros(11), adj, nt_idx)
    e0, v0 = eigh(H0)
    print(f"\nSpectrum at k=0:")
    print(f"  Eigenvalues: {np.round(e0, 6)}")
    print(f"  16-band group [0-15]: {np.round(e0[:16], 6)}")
    print(f"  3-band group [16-18]: {np.round(e0[16:19], 6)}")
    print(f"  Isolated band [19]:   {np.round(e0[19], 6)}")

    # A5 irrep decomposition of the 16-band group at k=0
    print(f"\n  A5 irrep decomposition at k=0:")
    print(f"    -sqrt(5) = {-np.sqrt(5):.6f} x3  (chi_3', dim 3)")
    print(f"    -2       = -2.000000        x4  (chi_4,  dim 4)")
    print(f"     0       =  0.000000        x4  (chi_4', dim 4)")
    print(f"     1       =  1.000000        x5  (chi_5,  dim 5)")
    print(f"    Total: 3 + 4 + 4 + 5 = 16")

    band_indices = list(range(16))
    m = len(band_indices)

    # Parameters
    N_loop = 30   # points around the Wilson loop
    N_perp = 60   # points in the perpendicular direction

    # Representative direction pairs
    # From previous Chern analysis:
    #   (0,7) had Chern = +1 for bands [16-18]
    #   We pick diverse pairs to probe different topological sectors
    pairs = [
        (0, 1),   # first two non-tree edges
        (0, 2),   # first + third
        (0, 3),   # different pairing
        (0, 7),   # the topologically active Chern slice
        (1, 4),   # another Chern-active pair
        (2, 5),   # symmetry-related
    ]

    print(f"\n{'='*72}")
    print(f"WILSON LOOP COMPUTATION")
    print(f"{'='*72}")
    print(f"  Band group: [0-15] (16 bands)")
    print(f"  N_loop = {N_loop} (discretization around loop)")
    print(f"  N_perp = {N_perp} (k_perp sampling)")
    print(f"  Direction pairs: {pairs}")

    all_spectra = {}
    all_windings = {}
    all_crossings = {}

    for pair_idx, (dir_perp, dir_loop) in enumerate(pairs):
        t0 = time.time()
        print(f"\n--- Pair ({dir_perp},{dir_loop}): "
              f"loop in dir {dir_loop}, k_perp in dir {dir_perp} ---")

        k_vals, theta_raw, unitarity_errs = wilson_loop_spectrum(
            dir_loop, dir_perp, adj, nt_idx, band_indices,
            N_loop=N_loop, N_perp=N_perp
        )

        dt = time.time() - t0
        print(f"  Computed in {dt:.1f}s")
        print(f"  Max unitarity error: {unitarity_errs.max():.2e}")
        print(f"  Mean unitarity error: {unitarity_errs.mean():.2e}")

        # Track eigenvalues continuously
        theta_tracked = track_eigenvalue_continuity(theta_raw)

        # Compute winding numbers
        windings = compute_winding_numbers(theta_tracked)
        total_winding = sum(windings)

        print(f"  Winding numbers: {np.round(windings, 3)}")
        print(f"  Total winding: {total_winding:.6f}")
        print(f"  Integer windings: {np.round(windings).astype(int)}")

        # Count crossings at theta = 0 and theta = pi
        cross_0 = count_theta_crossings(theta_tracked, 0.0)
        cross_pi = count_theta_crossings(theta_tracked, np.pi)
        cross_neg_pi = count_theta_crossings(theta_tracked, -np.pi)

        print(f"  Crossings at theta=0:   {cross_0}")
        print(f"  Crossings at theta=pi:  {cross_pi}")
        print(f"  Crossings at theta=-pi: {cross_neg_pi}")

        # Z_2 invariant: parity of crossings at pi
        z2 = (cross_pi + cross_neg_pi) % 2
        print(f"  Z_2 invariant (from pi crossings): {z2}")

        # Detect degeneracies
        degens = detect_degeneracies(theta_tracked, threshold=0.03)
        if degens:
            print(f"  Near-degeneracies: {len(degens)}")
            # Group by k_perp index
            degen_k = set(d[0] for d in degens)
            print(f"  At k_perp indices: {sorted(degen_k)[:10]}{'...' if len(degen_k)>10 else ''}")
        else:
            print(f"  No near-degeneracies found")

        # A5 clustering
        clustering = analyze_a5_clustering(theta_tracked, k_vals)
        k0_cluster = clustering[0]
        print(f"  A5 clustering at k_perp=0: groups {k0_cluster['group_sizes']}")

        # Check how clustering evolves
        size_evolution = [c['group_sizes'] for c in clustering]
        unique_patterns = list(set(tuple(s) for s in size_evolution))
        print(f"  Distinct clustering patterns: {len(unique_patterns)}")
        for pat in sorted(unique_patterns)[:5]:
            count = sum(1 for s in size_evolution if tuple(s) == pat)
            print(f"    {list(pat)}: {count}/{N_perp} samples")

        # Eigenvalue phase distribution at k_perp=0
        print(f"  Wilson eigenvalue phases at k_perp=0:")
        phases_at_0 = np.sort(theta_tracked[0])
        for n in range(m):
            print(f"    theta_{n:2d} = {phases_at_0[n]:+.6f}  "
                  f"(e^{{i*theta}} = {np.cos(phases_at_0[n]):+.4f} "
                  f"+ {np.sin(phases_at_0[n]):+.4f}i)")

        # ASCII visualization
        print()
        print(ascii_wilson_spectrum(
            k_vals, theta_tracked,
            title=f"  Wilson Loop Spectrum: loop dir={dir_loop}, k_perp dir={dir_perp}",
            width=72, height=24
        ))
        print()
        print(ascii_spectrum_density(
            k_vals, theta_tracked,
            title=f"  Density Plot: loop dir={dir_loop}, k_perp dir={dir_perp}",
            width=72, height=24
        ))

        # Store results
        all_spectra[(dir_perp, dir_loop)] = theta_tracked
        all_windings[(dir_perp, dir_loop)] = windings
        all_crossings[(dir_perp, dir_loop)] = (cross_0, cross_pi)

    # ================================================================
    # SYNTHESIS
    # ================================================================
    print(f"\n{'='*72}")
    print("SYNTHESIS: NON-ABELIAN TOPOLOGY OF THE 16-BAND GROUP")
    print(f"{'='*72}")

    # Winding number summary
    print(f"\n--- WINDING NUMBERS ---")
    print(f"{'Pair':>10} | {'Total':>8} | {'Non-zero':>10} | {'Int windings':>30}")
    print(f"{'-'*10}-+-{'-'*8}-+-{'-'*10}-+-{'-'*30}")
    for (dp, dl), winds in all_windings.items():
        int_winds = np.round(winds).astype(int)
        nonzero = np.count_nonzero(int_winds)
        total = int(np.round(sum(winds)))
        nz_vals = int_winds[int_winds != 0]
        print(f"  ({dp},{dl})   | {total:+6d}   | {nonzero:6d}     | {nz_vals}")

    # Crossing summary
    print(f"\n--- Z_2 CROSSINGS ---")
    print(f"{'Pair':>10} | {'theta=0':>8} | {'theta=pi':>10} | {'Z_2':>5}")
    print(f"{'-'*10}-+-{'-'*8}-+-{'-'*10}-+-{'-'*5}")
    for (dp, dl), (c0, cpi) in all_crossings.items():
        z2 = cpi % 2
        print(f"  ({dp},{dl})   | {c0:6d}   | {cpi:8d}   | {z2:3d}")

    # Eigenvalue spread analysis
    print(f"\n--- EIGENVALUE SPREAD ---")
    for (dp, dl), theta in all_spectra.items():
        spread = theta.max(axis=1) - theta.min(axis=1)
        full_range = theta.max() - theta.min()
        print(f"  ({dp},{dl}): range [{theta.min():.4f}, {theta.max():.4f}], "
              f"mean spread = {spread.mean():.4f}")

    # Look for universal structure
    print(f"\n--- UNIVERSAL STRUCTURE ---")
    # Check if the same eigenvalue pattern appears in all pairs
    phase_at_0 = {}
    for (dp, dl), theta in all_spectra.items():
        phase_at_0[(dp, dl)] = np.sort(theta[0])

    # Compare phase patterns between pairs
    print(f"  Eigenvalue phases at k_perp = 0:")
    for (dp, dl), phases in phase_at_0.items():
        print(f"    ({dp},{dl}): {np.round(phases, 3)}")

    # Check if eigenvalue phases at k_perp=0 are direction-independent
    # (they should be if the Wilson loop is trivial at k=0)
    ref_phases = list(phase_at_0.values())[0]
    all_same = True
    for phases in phase_at_0.values():
        if not np.allclose(phases, ref_phases, atol=0.1):
            all_same = False
            break
    if all_same:
        print(f"  -> Wilson loop phases at k_perp=0 are direction-INDEPENDENT")
        print(f"     This means the topology is encoded in how they EVOLVE")
    else:
        print(f"  -> Wilson loop phases at k_perp=0 are direction-DEPENDENT")
        print(f"     The topology depends on which loop we integrate around")

    # Topology classification
    print(f"\n--- TOPOLOGICAL CLASSIFICATION ---")

    total_nonzero_windings = 0
    total_z2 = 0
    for (dp, dl), winds in all_windings.items():
        int_winds = np.round(winds).astype(int)
        total_nonzero_windings += np.count_nonzero(int_winds)
    for (dp, dl), (c0, cpi) in all_crossings.items():
        total_z2 += cpi % 2

    if total_nonzero_windings == 0 and total_z2 == 0:
        print(f"  The 16-band group appears TOPOLOGICALLY TRIVIAL")
        print(f"  (no non-zero windings, no Z_2 invariants)")
        print(f"  BUT: it may have FRAGILE topology detectable only by")
        print(f"  higher-dimensional invariants (C_2 on 4D sub-tori)")
    elif total_nonzero_windings > 0:
        print(f"  NON-TRIVIAL WINDING detected!")
        print(f"  Total non-zero winding eigenvalues: {total_nonzero_windings}")
        print(f"  This indicates STABLE non-abelian topology")
    else:
        print(f"  No winding but Z_2 = {total_z2} non-trivial")
        print(f"  This indicates Z_2 TOPOLOGICAL INSULATOR structure")

    # Eigenvalue grouping interpretation
    print(f"\n--- A5 SYMMETRY INTERPRETATION ---")
    print(f"  At k=0, 16 bands = chi_3'(3) + chi_4(4) + chi_4'(4) + chi_5(5)")
    print(f"  Wilson loop eigenvalues should cluster as (3,4,4,5) near k_perp=0")
    print(f"  and evolve/mix as k_perp increases, breaking the A5 structure")
    print()
    print(f"  The mixing pattern encodes how the A5 irreps recombine")
    print(f"  under the symmetry-breaking effect of non-zero momentum.")
    print(f"  Eigenvalue crossings = irrep mixing = non-abelian Berry phase")

    # Final particle physics interpretation
    print(f"\n--- PHYSICAL INTERPRETATION ---")
    print(f"  The non-abelian Wilson loop of the 16-band group on T^11")
    print(f"  captures the gauge-invariant topological content that")
    print(f"  individual Chern numbers cannot access (bands cross).")
    print(f"")
    print(f"  Key structural features:")
    print(f"    - 16 = 3 + 4 + 4 + 5 (A5 irrep dimensions)")
    print(f"    - The Berry holonomy mixes these irreps non-trivially")
    print(f"    - Eigenvalue winding = protected edge modes")
    print(f"    - Eigenvalue degeneracies = symmetry-enforced crossings")
    print(f"")
    print(f"  For Standard Model analogy:")
    print(f"    16 = 16 of SO(10) (one generation of fermions!)")
    print(f"    3+4+4+5 could map to quark/lepton families under A5")
    print(f"    The Wilson loop spectrum encodes the 'flavor mixing'")

    # Convergence check on one pair
    print(f"\n{'='*72}")
    print("CONVERGENCE CHECK")
    print(f"{'='*72}")

    dp, dl = pairs[0]
    print(f"  Testing pair ({dp},{dl}) at k_perp = pi/2")
    for N_test in [15, 20, 30, 45, 60]:
        W, err = compute_wilson_loop(
            np.zeros(11), dl, dp, np.pi/2,
            adj, nt_idx, band_indices, N_test
        )
        evals = eigvals(W)
        phases = np.sort(np.angle(evals))
        det_W = det(W)
        print(f"  N={N_test:2d}: det(W)={abs(det_W):.6f}, "
              f"unitarity={err:.2e}, "
              f"phases[0:3]={np.round(phases[:3], 4)}")

    print(f"\n{'='*72}")
    print("COMPUTATION COMPLETE")
    print(f"{'='*72}")


if __name__ == "__main__":
    main()
