"""
DODECAHEDRAL_KNOT_SPECTRUM — classifies all cycles on the dodecahedron; eigenmode energies vs particle mass ratios
nos3bl33d

DFS cycle classification up to length 15, Laplacian eigenspace projection, combination energies, mass ratios.
"""

import numpy as np
from numpy.linalg import eigh
import networkx as nx
from itertools import combinations
from collections import defaultdict
import time

# --------------------------------------------------------------------------- #
# CONSTANTS
# --------------------------------------------------------------------------- #
phi = (1 + np.sqrt(5)) / 2
phi_sq = phi**2            # phi + 1
phi_4 = phi**4             # 3*phi + 2 = 6.854101966...
V = 20                     # vertices
E_edges = 30               # edges
F = 12                     # faces
d = 3                      # degree (regular graph)

# Known particle mass ratios (CODATA 2022 / PDG 2024)
KNOWN_RATIOS = {
    "m_p / m_e":        1836.15267343,
    "m_muon / m_e":     206.7682830,
    "m_tau / m_e":      3477.48,
    "m_tau / m_muon":   16.8170,
    "m_W / m_Z":        0.88145,
    "m_p / m_pion":     6.7224,
    "m_pion / m_e":     273.13,
    "m_pion / m_muon":  1.3210,
    "m_kaon / m_pion":  3.540,
    "m_proton / m_neutron": 0.99862,
    "m_W / m_e":        157299.0,
    "m_Z / m_e":        178450.0,
    "m_higgs / m_e":    245070.0,
}

# --------------------------------------------------------------------------- #
# 1. BUILD THE DODECAHEDRON
# --------------------------------------------------------------------------- #
def build_dodecahedron():
    """Build dodecahedron graph with adjacency list and adjacency matrix."""
    edges = [
        (0, 1), (0, 4), (0, 5),
        (1, 2), (1, 6),
        (2, 3), (2, 7),
        (3, 4), (3, 8),
        (4, 9),
        (5, 10), (5, 14),
        (6, 10), (6, 11),
        (7, 11), (7, 12),
        (8, 12), (8, 13),
        (9, 13), (9, 14),
        (10, 15),
        (11, 16),
        (12, 17),
        (13, 18),
        (14, 19),
        (15, 16), (15, 19),
        (16, 17),
        (17, 18),
        (18, 19),
    ]

    adj = [[] for _ in range(V)]
    A = np.zeros((V, V), dtype=float)

    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        A[u, v] = 1.0
        A[v, u] = 1.0

    # Validate
    assert len(edges) == E_edges
    for i in range(V):
        assert len(adj[i]) == d, f"Vertex {i} has degree {len(adj[i])}, expected {d}"

    # Degree matrix
    D = np.diag(np.full(V, float(d)))
    # Laplacian L = D - A
    L = D - A

    return adj, A, L, edges


def build_networkx_graph(edges):
    """Build NetworkX graph for cycle finding."""
    G = nx.Graph()
    G.add_nodes_from(range(V))
    G.add_edges_from(edges)
    return G


# --------------------------------------------------------------------------- #
# 2. COMPUTE LAPLACIAN SPECTRUM
# --------------------------------------------------------------------------- #
def compute_spectrum(L):
    """
    Compute eigenvalues and eigenvectors of the Laplacian.

    Known eigenvalues of the dodecahedron Laplacian:
      {0, 3-sqrt(5), 2, 3, 5, 3+sqrt(5)}
    with multiplicities {1, 3, 5, 4, 4, 3}
    """
    eigenvalues, eigenvectors = eigh(L)

    # Round to remove numerical noise, then find distinct values
    rounded = np.round(eigenvalues, 10)
    distinct = sorted(set(rounded))

    # Group by eigenvalue
    eigenspaces = {}
    for mu in distinct:
        mask = np.abs(rounded - mu) < 1e-8
        indices = np.where(mask)[0]
        eigenspaces[mu] = {
            'value': mu,
            'multiplicity': len(indices),
            'vectors': eigenvectors[:, indices],  # columns are eigenvectors
            'indices': indices,
        }

    return eigenvalues, eigenvectors, eigenspaces


# --------------------------------------------------------------------------- #
# 3. FIND ALL CYCLES ON THE DODECAHEDRON (up to max_length)
# --------------------------------------------------------------------------- #
def find_all_cycles_dfs(adj, max_length=15):
    """
    Find all simple cycles on the dodecahedron up to a given length.
    Uses DFS with backtracking. Returns cycles as lists of vertices.

    A cycle is a closed simple path: v0 -> v1 -> ... -> vk -> v0 where all
    vertices are distinct (except v0=vk).

    To avoid duplicates: we fix the minimum-index vertex as the start,
    and require the second vertex to be less than the last vertex
    (picks one of the two orientations).
    """
    all_cycles = []
    n = len(adj)

    def dfs(path, visited):
        current = path[-1]
        length = len(path)

        if length > max_length:
            return

        for neighbor in adj[current]:
            if neighbor == path[0] and length >= 3:
                # Found a cycle -- record it if canonical
                # Canonical: start is the minimum vertex, and path[1] < path[-1]
                cycle = list(path)
                min_idx = cycle.index(min(cycle))
                # Rotate so minimum is first
                cycle = cycle[min_idx:] + cycle[:min_idx]
                # Pick orientation: path[1] < path[-1]
                if len(cycle) >= 3 and cycle[1] > cycle[-1]:
                    cycle = [cycle[0]] + cycle[1:][::-1]

                cycle_tuple = tuple(cycle)
                all_cycles.append(cycle_tuple)

            elif neighbor not in visited and neighbor > path[0]:
                # Only explore neighbors with index > start vertex
                # (otherwise we'd double-count)
                # Actually: we need neighbor > start to avoid finding the same
                # cycle from multiple starting vertices. But we also need to
                # allow going to vertices < current. Let me fix this.
                pass

            if neighbor not in visited and length < max_length:
                if neighbor != path[0]:  # don't revisit start except to close
                    visited.add(neighbor)
                    path.append(neighbor)
                    dfs(path, visited)
                    path.pop()
                    visited.discard(neighbor)

    # Start DFS from each vertex (but canonicalize to avoid duplicates)
    for start in range(n):
        visited = {start}
        dfs([start], visited)

    # Deduplicate
    unique_cycles = list(set(all_cycles))
    unique_cycles.sort(key=lambda c: (len(c), c))

    return unique_cycles


def find_all_cycles_networkx(G, max_length=15):
    """
    Find all simple cycles using NetworkX, filtered by max length.
    More reliable than hand-rolled DFS for a small graph.
    """
    all_cycles = []
    for cycle in nx.simple_cycles(G, length_bound=max_length):
        if len(cycle) >= 3:
            # Canonicalize: rotate so minimum vertex is first
            min_idx = cycle.index(min(cycle))
            cycle = cycle[min_idx:] + cycle[:min_idx]
            # Pick orientation: second element < last element
            if len(cycle) >= 3 and cycle[1] > cycle[-1]:
                cycle = [cycle[0]] + cycle[1:][::-1]
            all_cycles.append(tuple(cycle))

    unique = list(set(all_cycles))
    unique.sort(key=lambda c: (len(c), c))
    return unique


# --------------------------------------------------------------------------- #
# 4. COMPUTE CYCLE ENERGIES
# --------------------------------------------------------------------------- #
def cycle_indicator_vector(cycle, n=V):
    """
    Create the indicator vector for a cycle.
    v[i] = 1 if vertex i is in the cycle, 0 otherwise.
    Normalized to unit length.
    """
    v = np.zeros(n)
    for vertex in cycle:
        v[vertex] = 1.0
    return v / np.linalg.norm(v)


def cycle_edge_indicator(cycle, n=V):
    """
    Create a vector weighted by the cycle's edge structure.
    For each vertex in the cycle, weight = number of cycle-edges incident to it / 2.
    For a simple cycle, every vertex has exactly 2 incident cycle-edges.
    So this is just 2/len(cycle) * indicator. But for vertex weighting we use
    the number of cycle edges touching that vertex.
    """
    v = np.zeros(n)
    L = len(cycle)
    for i in range(L):
        v[cycle[i]] += 1.0  # Each vertex in cycle gets weight 1
    # Every cycle vertex has exactly 2 cycle-edges, so this IS the indicator
    return v / np.linalg.norm(v)


def compute_cycle_energy(cycle, eigenspaces):
    """
    Compute the energy of a cycle by projecting onto eigenspaces.

    E_cycle = sum_i |<cycle | eigenspace_i>|^2 * mu_i

    where the sum is over eigenspaces, and |<cycle | eigenspace_i>|^2 is the
    total squared projection of the cycle indicator onto eigenspace i.
    """
    psi = cycle_indicator_vector(cycle)

    energy = 0.0
    projections = {}

    for mu, space in eigenspaces.items():
        # Project psi onto this eigenspace
        vectors = space['vectors']  # shape (V, multiplicity)
        # Total projection = sum of squared inner products with basis vectors
        proj_sq = 0.0
        for j in range(vectors.shape[1]):
            coeff = np.dot(psi, vectors[:, j])
            proj_sq += coeff**2

        projections[mu] = proj_sq
        energy += proj_sq * mu

    return energy, projections


def compute_cycle_laplacian_energy(cycle, L_matrix):
    """
    Alternative energy: E = <psi | L | psi> where psi is the cycle indicator.
    This is the Rayleigh quotient -- measures how much the cycle "oscillates."

    For a cycle of length k: this counts the number of edges in the cycle
    that connect vertices INSIDE the cycle to vertices OUTSIDE it,
    weighted appropriately.
    """
    psi = cycle_indicator_vector(cycle)
    return float(psi @ L_matrix @ psi)


def compute_cycle_edge_energy(cycle, adj):
    """
    Edge-based energy: count edges in the cycle vs edges leaving the cycle.

    Internal edges (within cycle): e_in
    Boundary edges (one end in cycle, one out): e_out

    Energy = e_out / e_in (higher = more "surface tension" = higher mass)
    """
    cycle_set = set(cycle)
    L = len(cycle)
    e_in = L  # The cycle itself has exactly L edges
    e_out = 0
    for v in cycle:
        for w in adj[v]:
            if w not in cycle_set:
                e_out += 1
    return e_out / e_in if e_in > 0 else 0.0


# --------------------------------------------------------------------------- #
# 5. TWO-MODE COMBINATIONS
# --------------------------------------------------------------------------- #
def compute_overlap_integral(eigvec_i, eigvec_j):
    """
    Compute overlap integral of two eigenvectors on the vertices.
    For orthogonal eigenvectors from different eigenspaces, this is 0.
    For eigenvectors in the same eigenspace, use their inner product.

    The "interaction energy" for a two-mode combination uses the
    vertex-wise product integral: sum_v psi_i(v) * psi_j(v)
    """
    return np.dot(eigvec_i, eigvec_j)


def compute_mode_pair_energies(eigenvalues, eigenvectors, eigenspaces):
    """
    For all pairs of eigenmodes, compute:
    E_combined = E_1 + E_2 + V_interaction

    where V_interaction = lambda * sum_v |psi_1(v)|^2 * |psi_2(v)|^2
    (density-density interaction, like a Hubbard-type term).

    The coupling lambda comes from the lattice: lambda = 1/V (mean field).
    """
    n_modes = len(eigenvalues)
    pairs = []

    for i in range(n_modes):
        for j in range(i, n_modes):
            mu_i = eigenvalues[i]
            mu_j = eigenvalues[j]
            psi_i = eigenvectors[:, i]
            psi_j = eigenvectors[:, j]

            # Direct sum
            E_sum = mu_i + mu_j

            # Density-density interaction
            rho_i = psi_i**2
            rho_j = psi_j**2
            V_int = np.sum(rho_i * rho_j) / V  # normalized by vertex count

            # Combined energy
            E_combined = E_sum + V_int

            pairs.append({
                'modes': (i, j),
                'eigenvalues': (mu_i, mu_j),
                'E_sum': E_sum,
                'V_interaction': V_int,
                'E_combined': E_combined,
            })

    return pairs


# --------------------------------------------------------------------------- #
# 6. RATIO MATCHING
# --------------------------------------------------------------------------- #
def find_ratio_matches(energies, known_ratios, tolerance=0.05):
    """
    Check if any ratio of energies matches known particle mass ratios.

    For each pair (E_i, E_j) with E_i > E_j > 0, compute E_i/E_j
    and check against known ratios.

    Also check sqrt(E_i/E_j) (harmonic oscillator dispersion).
    Also check (E_i/E_j)^2 (quadratic dispersion).
    Also check E_i^2 / E_j^2 (for E~m*c^2 with relativistic mode).
    """
    # Filter positive energies and sort
    pos_energies = sorted([e for e in energies if e > 1e-12])

    matches = []

    for i in range(len(pos_energies)):
        for j in range(len(pos_energies)):
            if i == j:
                continue
            E_i = pos_energies[i]
            E_j = pos_energies[j]
            if E_j < 1e-12:
                continue

            ratio = E_i / E_j
            sqrt_ratio = np.sqrt(ratio)
            sq_ratio = ratio**2

            for name, target in known_ratios.items():
                for test_ratio, rtype in [(ratio, "linear"),
                                           (sqrt_ratio, "sqrt"),
                                           (sq_ratio, "squared")]:
                    if target > 0:
                        rel_err = abs(test_ratio - target) / target
                        if rel_err < tolerance:
                            matches.append({
                                'name': name,
                                'target': target,
                                'computed': test_ratio,
                                'rel_error': rel_err,
                                'E_num': E_i,
                                'E_den': E_j,
                                'dispersion': rtype,
                            })

    return matches


def find_ratio_matches_from_list(energy_list, known_ratios, tolerance=0.05):
    """
    Same as above but takes a list of (label, energy) tuples.
    Returns matches with labels for identification.
    """
    # Filter positive
    pos = [(label, e) for label, e in energy_list if e > 1e-12]

    matches = []
    for i, (label_i, E_i) in enumerate(pos):
        for j, (label_j, E_j) in enumerate(pos):
            if i == j or E_j < 1e-12:
                continue

            ratio = E_i / E_j
            sqrt_ratio = np.sqrt(ratio) if ratio > 0 else 0
            sq_ratio = ratio**2

            for name, target in known_ratios.items():
                for test_ratio, rtype in [(ratio, "linear"),
                                           (sqrt_ratio, "sqrt"),
                                           (sq_ratio, "squared")]:
                    if target > 0 and test_ratio > 0:
                        rel_err = abs(test_ratio - target) / target
                        if rel_err < tolerance:
                            matches.append({
                                'particle_ratio': name,
                                'target': target,
                                'computed': test_ratio,
                                'rel_error': rel_err,
                                'E_numerator': (label_i, E_i),
                                'E_denominator': (label_j, E_j),
                                'dispersion': rtype,
                            })

    return matches


# --------------------------------------------------------------------------- #
# 7. EXTENDED CYCLE ENERGY MODELS
# --------------------------------------------------------------------------- #
def compute_winding_energy(cycle, L_matrix, eigenspaces):
    """
    Compute an energy that accounts for the topological "winding" of the cycle.

    For a cycle C of length L on the dodecahedron:
    E_winding = L^2 / V * <psi_C | L | psi_C>

    The L^2/V factor comes from the winding: on the universal cover, a cycle
    of length L wraps around ~L/girth times. The girth of the dodecahedron is 5.
    """
    psi = cycle_indicator_vector(cycle)
    laplacian_energy = float(psi @ L_matrix @ psi)
    L = len(cycle)
    winding_number = L / 5.0  # approximate winding in units of girth
    return winding_number**2 * laplacian_energy


def compute_golden_energy(cycle, eigenspaces):
    """
    Golden energy model: each eigenvalue contributes weighted by phi.

    E_golden = sum_i |c_i|^2 * phi^(2*rank_i) * mu_i

    where rank_i is the position of eigenvalue mu_i in the sorted list (0-indexed),
    and c_i = projection of cycle onto eigenspace i.
    """
    psi = cycle_indicator_vector(cycle)
    sorted_mus = sorted(eigenspaces.keys())

    energy = 0.0
    for rank, mu in enumerate(sorted_mus):
        space = eigenspaces[mu]
        vectors = space['vectors']
        proj_sq = sum(np.dot(psi, vectors[:, j])**2 for j in range(vectors.shape[1]))
        energy += proj_sq * phi**(2 * rank) * mu

    return energy


def compute_phi_weighted_energy(cycle, eigenspaces):
    """
    Phi-weighted energy: eigenvalues scaled by their phi-connection.

    In the framework, eigenvalue mu relates to phi via the characteristic
    polynomial. The energy uses the phi-scaling directly:

    E_phi = sum_i |c_i|^2 * mu_i * (mu_i / phi)^2

    This weights higher eigenvalues more strongly (massive modes).
    """
    psi = cycle_indicator_vector(cycle)

    energy = 0.0
    for mu, space in eigenspaces.items():
        vectors = space['vectors']
        proj_sq = sum(np.dot(psi, vectors[:, j])**2 for j in range(vectors.shape[1]))
        if mu > 1e-12:
            energy += proj_sq * mu * (mu / phi)**2

    return energy


# --------------------------------------------------------------------------- #
# 8. SYMMETRY CLASSIFICATION
# --------------------------------------------------------------------------- #
def classify_cycle_symmetry(cycle, G):
    """
    Classify cycle by its automorphism group.

    Returns:
    - n_automorphisms: number of graph automorphisms that map cycle to itself
    - is_face: whether the cycle is a pentagonal face
    - is_hamiltonian: whether the cycle visits all 20 vertices
    """
    cycle_set = frozenset(cycle)
    L = len(cycle)

    is_face = (L == 5)  # girth = 5, faces are pentagons
    is_hamiltonian = (L == V)

    # Check if cycle forms a face by verifying all edges exist
    if is_face:
        cycle_edges = set()
        for i in range(L):
            e = tuple(sorted([cycle[i], cycle[(i+1) % L]]))
            cycle_edges.add(e)
        # It's a face if all 5 edges are graph edges
        all_edges_exist = all(G.has_edge(u, v) for u, v in cycle_edges)
        is_face = all_edges_exist

    return {
        'length': L,
        'n_vertices': len(set(cycle)),
        'is_face': is_face,
        'is_hamiltonian': is_hamiltonian,
    }


# --------------------------------------------------------------------------- #
# MAIN COMPUTATION
# --------------------------------------------------------------------------- #
def main():
    print("=" * 80)
    print("DODECAHEDRAL KNOT SPECTRUM: Cycle Classification & Eigenmode Energies")
    print("=" * 80)
    print()

    # ------------------------------------------------------------------- #
    # STEP 1: Build the dodecahedron
    # ------------------------------------------------------------------- #
    print("--- STEP 1: Building dodecahedron ---")
    adj, A, L_matrix, edges = build_dodecahedron()
    G = build_networkx_graph(edges)

    print(f"  Vertices: {V}")
    print(f"  Edges:    {E_edges}")
    print(f"  Faces:    {F}")
    print(f"  Degree:   {d}")
    print(f"  Girth:    {nx.girth(G)}")
    print(f"  Diameter: {nx.diameter(G)}")
    print()

    # ------------------------------------------------------------------- #
    # STEP 2: Compute Laplacian spectrum
    # ------------------------------------------------------------------- #
    print("--- STEP 2: Laplacian Spectrum ---")
    eigenvalues, eigenvectors, eigenspaces = compute_spectrum(L_matrix)

    print(f"  {'Eigenvalue':>12}  {'Exact':>14}  {'Multiplicity':>12}  {'phi-relation':>20}")
    print(f"  {'-'*12}  {'-'*14}  {'-'*12}  {'-'*20}")

    exact_names = {
        0.0: "0",
        round(3 - np.sqrt(5), 10): "3 - sqrt(5)",
        2.0: "2",
        3.0: "3",
        5.0: "5",
        round(3 + np.sqrt(5), 10): "3 + sqrt(5)",
    }

    phi_relations = {
        0.0: "trivial",
        round(3 - np.sqrt(5), 10): "= 3 - sqrt(5) = 2/phi^2",
        2.0: "= 2",
        3.0: "= 3 = phi^2 + 1/phi^2",
        5.0: "= 5 = phi^4 - phi^2 + 1",
        round(3 + np.sqrt(5), 10): "= 3 + sqrt(5) = 2*phi^2",
    }

    sorted_mus = sorted(eigenspaces.keys())
    for mu in sorted_mus:
        space = eigenspaces[mu]
        exact = exact_names.get(round(mu, 10), "?")
        phi_rel = phi_relations.get(round(mu, 10), "?")
        print(f"  {mu:12.8f}  {exact:>14}  {space['multiplicity']:>12}  {phi_rel}")

    print(f"\n  Total eigenvalues: {sum(s['multiplicity'] for s in eigenspaces.values())} (should be {V})")
    print(f"  Trace of L (= sum of eigenvalues): {sum(eigenvalues):.6f} (should be {d * V} = {d*V})")
    print()

    # ------------------------------------------------------------------- #
    # STEP 3: Single eigenmode energies and ratios
    # ------------------------------------------------------------------- #
    print("--- STEP 3: Single Eigenmode Energies ---")
    print()
    print("  If E ~ mu (linear dispersion):")
    nonzero_mus = [mu for mu in sorted_mus if mu > 1e-12]
    mu_min = min(nonzero_mus)

    single_energies = []
    for mu in nonzero_mus:
        ratio = mu / mu_min
        sqrt_ratio = np.sqrt(ratio)
        single_energies.append(("mu=" + exact_names.get(round(mu, 10), f"{mu:.4f}"), mu))
        print(f"    mu = {mu:.8f}  E/E_min = {ratio:.6f}  sqrt(E/E_min) = {sqrt_ratio:.6f}")

    print()
    print("  Eigenvalue ratios (all pairs):")
    for i, mu_i in enumerate(nonzero_mus):
        for j, mu_j in enumerate(nonzero_mus):
            if i < j:
                r = mu_i / mu_j
                name_i = exact_names.get(round(mu_i, 10), f"{mu_i:.4f}")
                name_j = exact_names.get(round(mu_j, 10), f"{mu_j:.4f}")
                print(f"    ({name_i}) / ({name_j}) = {r:.8f},  sqrt = {np.sqrt(r):.8f}")

    print()
    print("  phi-based eigenvalue ratios:")
    for mu in nonzero_mus:
        for k in range(-4, 10):
            target = phi**k
            if abs(mu / mu_min - target) / target < 0.01:
                print(f"    mu={mu:.6f}/mu_min={mu_min:.6f} = {mu/mu_min:.8f} ~= phi^{k} = {target:.8f}")
            if abs(np.sqrt(mu / mu_min) - target) / target < 0.01:
                print(f"    sqrt(mu={mu:.6f}/mu_min={mu_min:.6f}) = {np.sqrt(mu/mu_min):.8f} ~= phi^{k} = {target:.8f}")

    print()

    # ------------------------------------------------------------------- #
    # STEP 4: Find ALL cycles
    # ------------------------------------------------------------------- #
    print("--- STEP 4: Finding all cycles (up to length 15) ---")
    t0 = time.time()

    cycles = find_all_cycles_networkx(G, max_length=15)

    t1 = time.time()
    print(f"  Found {len(cycles)} distinct cycles in {t1-t0:.2f}s")

    # Count by length
    length_counts = defaultdict(int)
    for c in cycles:
        length_counts[len(c)] += 1

    print(f"\n  {'Length':>6}  {'Count':>8}")
    print(f"  {'-'*6}  {'-'*8}")
    for L in sorted(length_counts.keys()):
        print(f"  {L:6d}  {length_counts[L]:8d}")

    print()

    # ------------------------------------------------------------------- #
    # STEP 5: Classify cycles and compute energies
    # ------------------------------------------------------------------- #
    print("--- STEP 5: Cycle Energy Spectrum ---")
    print()

    cycle_data = []
    for cycle in cycles:
        sym = classify_cycle_symmetry(cycle, G)
        E_proj, projections = compute_cycle_energy(cycle, eigenspaces)
        E_lap = compute_cycle_laplacian_energy(cycle, L_matrix)
        E_edge = compute_cycle_edge_energy(cycle, adj)
        E_wind = compute_winding_energy(cycle, L_matrix, eigenspaces)
        E_gold = compute_golden_energy(cycle, eigenspaces)
        E_phi = compute_phi_weighted_energy(cycle, eigenspaces)

        cycle_data.append({
            'cycle': cycle,
            'length': len(cycle),
            'symmetry': sym,
            'E_projection': E_proj,
            'E_laplacian': E_lap,
            'E_edge': E_edge,
            'E_winding': E_wind,
            'E_golden': E_gold,
            'E_phi_weighted': E_phi,
            'projections': projections,
        })

    # Group by energy (E_laplacian, rounded)
    print("  ENERGY MODEL 1: Rayleigh quotient E = <psi|L|psi>")
    print()
    energy_groups_lap = defaultdict(list)
    for cd in cycle_data:
        key = round(cd['E_laplacian'], 6)
        energy_groups_lap[key].append(cd)

    print(f"  {'E_laplacian':>12}  {'Count':>6}  {'Lengths':>20}  {'Has face?':>10}")
    print(f"  {'-'*12}  {'-'*6}  {'-'*20}  {'-'*10}")
    for E in sorted(energy_groups_lap.keys()):
        group = energy_groups_lap[E]
        lengths = sorted(set(cd['length'] for cd in group))
        has_face = any(cd['symmetry']['is_face'] for cd in group)
        has_ham = any(cd['symmetry']['is_hamiltonian'] for cd in group)
        label = ""
        if has_face:
            label += " [FACE]"
        if has_ham:
            label += " [HAMILTONIAN]"
        print(f"  {E:12.6f}  {len(group):6d}  {str(lengths):>20}  {label}")

    print()

    # Distinct Laplacian energies
    distinct_E_lap = sorted(energy_groups_lap.keys())
    print(f"  Distinct Laplacian energies: {len(distinct_E_lap)}")
    print()

    # ------------------------------------------------------------------- #
    # STEP 5b: Projection energy model
    # ------------------------------------------------------------------- #
    print("  ENERGY MODEL 2: Eigenspace projection E = sum |c_i|^2 * mu_i")
    print()
    energy_groups_proj = defaultdict(list)
    for cd in cycle_data:
        key = round(cd['E_projection'], 6)
        energy_groups_proj[key].append(cd)

    distinct_E_proj = sorted(energy_groups_proj.keys())
    print(f"  Distinct projection energies: {len(distinct_E_proj)}")
    print(f"  Min: {distinct_E_proj[0]:.6f}  Max: {distinct_E_proj[-1]:.6f}")
    print()

    # Show first 30 distinct energies
    print(f"  {'#':>4}  {'E_proj':>12}  {'Count':>6}  {'Lengths':>15}")
    print(f"  {'-'*4}  {'-'*12}  {'-'*6}  {'-'*15}")
    for idx, E in enumerate(distinct_E_proj[:30]):
        group = energy_groups_proj[E]
        lengths = sorted(set(cd['length'] for cd in group))
        print(f"  {idx:4d}  {E:12.8f}  {len(group):6d}  {str(lengths):>15}")
    if len(distinct_E_proj) > 30:
        print(f"  ... ({len(distinct_E_proj) - 30} more)")
    print()

    # ------------------------------------------------------------------- #
    # STEP 5c: Golden energy model
    # ------------------------------------------------------------------- #
    print("  ENERGY MODEL 3: Golden energy E_golden = sum |c_i|^2 * phi^(2*rank) * mu_i")
    print()
    energy_groups_gold = defaultdict(list)
    for cd in cycle_data:
        key = round(cd['E_golden'], 4)
        energy_groups_gold[key].append(cd)

    distinct_E_gold = sorted(energy_groups_gold.keys())
    print(f"  Distinct golden energies: {len(distinct_E_gold)}")
    print(f"  Min: {distinct_E_gold[0]:.6f}  Max: {distinct_E_gold[-1]:.6f}")
    print()

    # ------------------------------------------------------------------- #
    # STEP 6: Eigenspace projection profiles
    # ------------------------------------------------------------------- #
    print("--- STEP 6: Eigenspace Projection Profiles ---")
    print()
    print("  How cycles project onto each eigenspace:")
    print()

    # For each cycle length, show the average projection profile
    for L in sorted(length_counts.keys())[:8]:
        cycles_of_length = [cd for cd in cycle_data if cd['length'] == L]
        n_cycles = len(cycles_of_length)

        avg_proj = {}
        for mu in sorted_mus:
            vals = [cd['projections'][mu] for cd in cycles_of_length]
            avg_proj[mu] = np.mean(vals)

        print(f"  Length {L:2d} ({n_cycles:4d} cycles):")
        for mu in sorted_mus:
            bar = "#" * int(avg_proj[mu] * 50)
            exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
            print(f"    mu={exact:>14}: |c|^2 = {avg_proj[mu]:.6f}  {bar}")
        print()

    # ------------------------------------------------------------------- #
    # STEP 7: Two-mode combinations
    # ------------------------------------------------------------------- #
    print("--- STEP 7: Two-Mode Combination Energies ---")
    print()

    pairs = compute_mode_pair_energies(eigenvalues, eigenvectors, eigenspaces)

    # Get distinct combined energies
    pair_energies = []
    for p in pairs:
        if p['E_combined'] > 1e-12:
            label = f"modes({p['modes'][0]},{p['modes'][1]})"
            pair_energies.append((label, p['E_combined']))

    distinct_pair_E = sorted(set(round(e, 6) for _, e in pair_energies))
    print(f"  Total mode pairs: {len(pairs)}")
    print(f"  Distinct combined energies: {len(distinct_pair_E)}")
    print(f"  Range: [{distinct_pair_E[0]:.6f}, {distinct_pair_E[-1]:.6f}]")
    print()

    # Show some
    print(f"  {'E_combined':>12}  {'E_sum':>10}  {'V_int':>10}  {'Modes':>15}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*15}")
    shown = set()
    for p in sorted(pairs, key=lambda x: x['E_combined']):
        key = round(p['E_combined'], 6)
        if key not in shown and p['E_combined'] > 1e-12:
            shown.add(key)
            print(f"  {p['E_combined']:12.6f}  {p['E_sum']:10.6f}  {p['V_interaction']:10.6f}  {str(p['modes']):>15}")
        if len(shown) >= 20:
            break
    print()

    # ------------------------------------------------------------------- #
    # STEP 8: Ratio matching against known particle masses
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 8: RATIO MATCHING AGAINST KNOWN PARTICLE MASSES")
    print("=" * 80)
    print()

    # Collect ALL energy values from all models
    all_energies = []

    # Single eigenvalues
    for mu in nonzero_mus:
        exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
        all_energies.append((f"eigenvalue {exact}", mu))

    # Cycle Laplacian energies (distinct)
    for E in distinct_E_lap:
        if E > 1e-12:
            all_energies.append((f"cycle_lap E={E:.4f}", E))

    # Cycle projection energies (distinct)
    for E in distinct_E_proj:
        if E > 1e-12:
            all_energies.append((f"cycle_proj E={E:.4f}", E))

    # Cycle golden energies (distinct)
    for E in distinct_E_gold:
        if E > 1e-12:
            all_energies.append((f"cycle_gold E={E:.4f}", E))

    # Two-mode pair energies (distinct)
    for E in distinct_pair_E:
        if E > 1e-12:
            all_energies.append((f"pair E={E:.4f}", E))

    print(f"  Total distinct energy values to check: {len(all_energies)}")
    print()

    # Find matches at 5% tolerance
    print("  --- Matches at 5% tolerance ---")
    matches = find_ratio_matches_from_list(all_energies, KNOWN_RATIOS, tolerance=0.05)

    if matches:
        # Sort by error
        matches.sort(key=lambda m: m['rel_error'])

        # Deduplicate by (particle_ratio, dispersion)
        seen = set()
        unique_matches = []
        for m in matches:
            key = (m['particle_ratio'], m['dispersion'], round(m['computed'], 4))
            if key not in seen:
                seen.add(key)
                unique_matches.append(m)

        print(f"  Found {len(unique_matches)} unique matches:")
        print()
        print(f"  {'Particle Ratio':>25}  {'Target':>12}  {'Computed':>12}  {'Error%':>8}  {'Dispersion':>10}  {'Numerator':>30}  {'Denominator':>30}")
        print(f"  {'-'*25}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*10}  {'-'*30}  {'-'*30}")

        for m in unique_matches[:50]:
            print(f"  {m['particle_ratio']:>25}  {m['target']:12.4f}  {m['computed']:12.4f}  {m['rel_error']*100:7.2f}%  {m['dispersion']:>10}  {m['E_numerator'][0]:>30}  {m['E_denominator'][0]:>30}")
    else:
        print("  No matches found at 5% tolerance.")

    print()

    # Also check at 10% tolerance for more hits
    print("  --- Matches at 10% tolerance (broader search) ---")
    matches_10 = find_ratio_matches_from_list(all_energies, KNOWN_RATIOS, tolerance=0.10)

    if matches_10:
        matches_10.sort(key=lambda m: m['rel_error'])
        seen = set()
        unique_10 = []
        for m in matches_10:
            key = (m['particle_ratio'], m['dispersion'], round(m['computed'], 4))
            if key not in seen:
                seen.add(key)
                unique_10.append(m)

        # Only show ones not already in 5% list
        new_matches = [m for m in unique_10 if m['rel_error'] > 0.05]
        if new_matches:
            print(f"  Additional matches (5-10% error): {len(new_matches)}")
            for m in new_matches[:30]:
                print(f"    {m['particle_ratio']:>25}  target={m['target']:12.4f}  computed={m['computed']:12.4f}  err={m['rel_error']*100:.1f}%  [{m['dispersion']}]")

    print()

    # ------------------------------------------------------------------- #
    # STEP 9: PHI-POWER SPECTRUM
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 9: PHI-POWER SPECTRUM — Can energies be expressed as phi^n?")
    print("=" * 80)
    print()

    print("  Eigenvalues as phi powers:")
    for mu in nonzero_mus:
        log_phi = np.log(mu) / np.log(phi)
        exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
        nearest_int = round(log_phi)
        err = abs(log_phi - nearest_int)
        marker = " <-- INTEGER!" if err < 0.01 else ""
        print(f"    mu = {exact:>14} = {mu:.8f}  log_phi(mu) = {log_phi:8.4f}  nearest int = {nearest_int}  err = {err:.4f}{marker}")

    print()

    # Check specific phi-power ratios
    print("  Checking if eigenvalue RATIOS are phi powers:")
    for i, mu_i in enumerate(nonzero_mus):
        for j, mu_j in enumerate(nonzero_mus):
            if i < j:
                r = mu_i / mu_j
                log_phi_r = np.log(r) / np.log(phi)
                nearest = round(log_phi_r)
                err = abs(log_phi_r - nearest)
                if err < 0.02:
                    name_i = exact_names.get(round(mu_i, 10), f"{mu_i:.4f}")
                    name_j = exact_names.get(round(mu_j, 10), f"{mu_j:.4f}")
                    print(f"    ({name_i})/({name_j}) = {r:.8f} = phi^{nearest} (err = {err:.6f})")

    print()

    # ------------------------------------------------------------------- #
    # STEP 10: CYCLE ENERGY RATIOS — The FULL spectrum
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 10: CYCLE ENERGY SPECTRUM — Distinct energy levels")
    print("=" * 80)
    print()

    # Use Laplacian energy (most physical: Rayleigh quotient)
    print("  Laplacian (Rayleigh) energy levels:")
    print()
    E_lap_min = min(E for E in distinct_E_lap if E > 1e-12)

    print(f"  {'#':>3}  {'E_lap':>12}  {'E/E_min':>10}  {'log_phi(E/E_min)':>18}  {'Count':>6}  {'Cycle lengths':>20}")
    print(f"  {'-'*3}  {'-'*12}  {'-'*10}  {'-'*18}  {'-'*6}  {'-'*20}")

    for idx, E in enumerate(distinct_E_lap):
        if E < 1e-12:
            continue
        ratio = E / E_lap_min
        log_phi_r = np.log(ratio) / np.log(phi) if ratio > 0 else 0
        group = energy_groups_lap[E]
        lengths = sorted(set(cd['length'] for cd in group))
        print(f"  {idx:3d}  {E:12.6f}  {ratio:10.6f}  {log_phi_r:18.6f}  {len(group):6d}  {str(lengths):>20}")

    print()

    # ------------------------------------------------------------------- #
    # STEP 11: FACE vs NON-FACE CYCLE ANALYSIS
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 11: FACE vs NON-FACE CYCLE ANALYSIS")
    print("=" * 80)
    print()

    face_cycles = [cd for cd in cycle_data if cd['symmetry']['is_face']]
    ham_cycles = [cd for cd in cycle_data if cd['symmetry']['is_hamiltonian']]

    print(f"  Face cycles (pentagonal faces): {len(face_cycles)}")
    if face_cycles:
        face_E = face_cycles[0]['E_laplacian']
        face_E_proj = face_cycles[0]['E_projection']
        print(f"    E_laplacian  = {face_E:.8f}")
        print(f"    E_projection = {face_E_proj:.8f}")
        print(f"    Projections onto eigenspaces:")
        for mu in sorted_mus:
            exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
            proj = face_cycles[0]['projections'][mu]
            print(f"      mu={exact:>14}: |c|^2 = {proj:.8f}")

    print()
    print(f"  Hamiltonian cycles: {len(ham_cycles)}")
    if ham_cycles:
        ham_E = ham_cycles[0]['E_laplacian']
        ham_E_proj = ham_cycles[0]['E_projection']
        print(f"    E_laplacian  = {ham_E:.8f}")
        print(f"    E_projection = {ham_E_proj:.8f}")
        if face_cycles:
            print(f"    E_hamiltonian / E_face = {ham_E / face_E:.8f}")
            print(f"    E_proj_ham / E_proj_face = {ham_E_proj / face_E_proj:.8f}")

    print()

    # ------------------------------------------------------------------- #
    # STEP 12: COMBINATORIAL ENERGY — Products and powers of eigenvalues
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 12: COMBINATORIAL ENERGY — Products and powers")
    print("=" * 80)
    print()

    # Build composite energies from eigenvalue products
    composite_energies = []

    # Single eigenvalues
    for mu in nonzero_mus:
        exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
        composite_energies.append((f"mu={exact}", mu))

    # Products of 2 eigenvalues
    for i, mu_i in enumerate(nonzero_mus):
        for j, mu_j in enumerate(nonzero_mus):
            if i <= j:
                name_i = exact_names.get(round(mu_i, 10), f"{mu_i:.4f}")
                name_j = exact_names.get(round(mu_j, 10), f"{mu_j:.4f}")
                composite_energies.append((f"mu({name_i})*mu({name_j})", mu_i * mu_j))

    # Products of 3 eigenvalues
    for i, mu_i in enumerate(nonzero_mus):
        for j, mu_j in enumerate(nonzero_mus):
            for k, mu_k in enumerate(nonzero_mus):
                if i <= j <= k:
                    composite_energies.append((f"3-prod", mu_i * mu_j * mu_k))

    # Powers of eigenvalues
    for mu in nonzero_mus:
        exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
        for p in range(2, 8):
            composite_energies.append((f"mu({exact})^{p}", mu**p))

    # Phi-weighted eigenvalues
    for mu in nonzero_mus:
        exact = exact_names.get(round(mu, 10), f"{mu:.4f}")
        for k in range(-4, 10):
            composite_energies.append((f"mu({exact})*phi^{k}", mu * phi**k))

    print(f"  Total composite energies: {len(composite_energies)}")

    # Check ratios
    print("  Checking ratios of composite energies against particle masses...")
    comp_matches = find_ratio_matches_from_list(composite_energies, KNOWN_RATIOS, tolerance=0.02)

    if comp_matches:
        comp_matches.sort(key=lambda m: m['rel_error'])

        seen = set()
        unique_comp = []
        for m in comp_matches:
            key = (m['particle_ratio'], round(m['computed'], 3))
            if key not in seen:
                seen.add(key)
                unique_comp.append(m)

        print(f"\n  Found {len(unique_comp)} matches at 2% tolerance:")
        print()
        print(f"  {'Particle Ratio':>25}  {'Target':>12}  {'Computed':>12}  {'Error%':>8}  {'Type':>8}  {'Num':>25}  {'Den':>25}")
        print(f"  {'-'*25}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*25}  {'-'*25}")

        for m in unique_comp[:60]:
            print(f"  {m['particle_ratio']:>25}  {m['target']:12.4f}  {m['computed']:12.4f}  {m['rel_error']*100:7.3f}%  {m['dispersion']:>8}  {m['E_numerator'][0]:>25}  {m['E_denominator'][0]:>25}")
    else:
        print("  No composite matches at 2% tolerance.")

    print()

    # ------------------------------------------------------------------- #
    # STEP 13: THE CRITICAL QUESTION — Can phi^N build mass ratios?
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("STEP 13: PURE PHI-POWER MASS RATIOS")
    print("=" * 80)
    print()

    print("  Can known mass ratios be expressed as phi^n * (lattice factor)?")
    print()

    lattice_factors = {
        "1": 1,
        "V": V,
        "E": E_edges,
        "F": F,
        "V/F": V/F,
        "E/V": E_edges/V,
        "E/F": E_edges/F,
        "V*F": V*F,
        "V*E": V*E_edges,
        "pi": np.pi,
        "2*pi": 2*np.pi,
        "4*pi": 4*np.pi,
        "pi^2": np.pi**2,
        "6*pi^5": 6*np.pi**5,
        "d": d,
        "chi": 2,
        "V/d": V/d,
    }

    for ratio_name, ratio_val in KNOWN_RATIOS.items():
        best_match = None
        best_err = 1e10

        for lf_name, lf_val in lattice_factors.items():
            for n in range(-20, 30):
                candidate = phi**n * lf_val
                if candidate > 0:
                    err = abs(candidate - ratio_val) / ratio_val
                    if err < best_err:
                        best_err = err
                        best_match = (lf_name, n, candidate, err)

        if best_match and best_match[3] < 0.10:
            lf_name, n, candidate, err = best_match
            print(f"  {ratio_name:>25} = {ratio_val:12.4f}  ~= phi^{n:3d} * {lf_name:>8} = {candidate:12.4f}  err = {err*100:.2f}%")
        else:
            # Try two phi powers
            for lf_name, lf_val in lattice_factors.items():
                for n1 in range(-10, 20):
                    for n2 in range(n1, 20):
                        candidate = (phi**n1 + phi**n2) * lf_val
                        if candidate > 0:
                            err = abs(candidate - ratio_val) / ratio_val
                            if err < best_err:
                                best_err = err
                                best_match = (f"(phi^{n1}+phi^{n2})*{lf_name}", candidate, err)

            if best_match and isinstance(best_match, tuple) and len(best_match) == 3:
                formula, candidate, err = best_match
                if err < 0.10:
                    print(f"  {ratio_name:>25} = {ratio_val:12.4f}  ~= {formula} = {candidate:12.4f}  err = {err*100:.2f}%")
                else:
                    print(f"  {ratio_name:>25} = {ratio_val:12.4f}  BEST: err = {err*100:.2f}% (> 10%)")
            else:
                if best_match:
                    print(f"  {ratio_name:>25} = {ratio_val:12.4f}  ~= phi^{best_match[1]:3d} * {best_match[0]:>8} = {best_match[2]:12.4f}  err = {best_match[3]*100:.2f}%")
                else:
                    print(f"  {ratio_name:>25} = {ratio_val:12.4f}  NO MATCH")

    print()

    # ------------------------------------------------------------------- #
    # STEP 14: SUMMARY
    # ------------------------------------------------------------------- #
    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print()

    print(f"  Dodecahedron: V={V}, E={E_edges}, F={F}, degree={d}")
    print(f"  Laplacian eigenvalues: {[round(mu, 6) for mu in sorted_mus]}")
    print(f"  Multiplicities: {[eigenspaces[mu]['multiplicity'] for mu in sorted_mus]}")
    print()
    print(f"  Total cycles found (length <= 15): {len(cycles)}")
    print(f"  Distinct Laplacian energies: {len(distinct_E_lap)}")
    print(f"  Distinct projection energies: {len(distinct_E_proj)}")
    print(f"  Face cycles: {len(face_cycles)} (all 12 pentagonal faces)")
    print(f"  Hamiltonian cycles: {len(ham_cycles)}")
    print()

    # The key eigenvalue identities
    print("  KEY EIGENVALUE IDENTITIES:")
    print(f"    3 - sqrt(5) = 2/phi^2 = {3-np.sqrt(5):.10f}")
    print(f"    3 + sqrt(5) = 2*phi^2  = {3+np.sqrt(5):.10f}")
    print(f"    ratio = (3+sqrt(5))/(3-sqrt(5)) = phi^4 = {(3+np.sqrt(5))/(3-np.sqrt(5)):.10f}")
    print(f"    phi^4 = {phi**4:.10f}")
    print(f"    Product: (3-sqrt(5))*(3+sqrt(5)) = 9-5 = 4 = 2^2")
    print(f"    Sum: (3-sqrt(5))+(3+sqrt(5)) = 6 = 2*d = 2*V/E * E")
    print()

    print("  SPECTRAL GAP STRUCTURE:")
    print(f"    Gap_1 = mu_1 - 0 = {nonzero_mus[0]:.10f} (= 3-sqrt(5) = 2/phi^2)")
    print(f"    Gap_2 = mu_2 - mu_1 = {nonzero_mus[1] - nonzero_mus[0]:.10f}")
    print(f"    Gap ratio = Gap_2/Gap_1 = {(nonzero_mus[1] - nonzero_mus[0])/nonzero_mus[0]:.10f}")
    print(f"    phi^2 - 1 = phi = {phi:.10f}")
    print(f"    Spectral gap ratio vs phi: {(nonzero_mus[1] - nonzero_mus[0])/nonzero_mus[0] / phi:.10f}")
    print()

    # Final ratio check summary
    print("  PARTICLE MASS RATIOS vs DODECAHEDRAL SPECTRUM:")
    print()

    # Compute specific ratio checks
    ratio_checks = [
        ("m_proton/m_electron", 1836.15267343),
        ("m_muon/m_electron", 206.7682830),
        ("m_tau/m_muon", 16.8170),
        ("m_W/m_Z", 0.88145),
        ("m_pion/m_muon", 1.3210),
    ]

    for name, target in ratio_checks:
        log_phi_target = np.log(target) / np.log(phi)
        # Nearest integer phi power
        n_near = round(log_phi_target)
        phi_approx = phi**n_near
        err_phi = abs(phi_approx - target) / target * 100

        # Check lattice combinations
        best_lattice = None
        best_err_l = 1e10
        for lf_name, lf_val in lattice_factors.items():
            residual = target / lf_val
            if residual > 0:
                log_phi_res = np.log(residual) / np.log(phi)
                n_res = round(log_phi_res)
                candidate = phi**n_res * lf_val
                err_l = abs(candidate - target) / target
                if err_l < best_err_l:
                    best_err_l = err_l
                    best_lattice = (lf_name, n_res, candidate)

        print(f"    {name:>25} = {target:.6f}")
        print(f"      log_phi = {log_phi_target:.6f}")
        print(f"      nearest phi^{n_near} = {phi_approx:.6f} (err = {err_phi:.2f}%)")
        if best_lattice:
            print(f"      best: phi^{best_lattice[1]} * {best_lattice[0]} = {best_lattice[2]:.6f} (err = {best_err_l*100:.2f}%)")
        print()

    print("=" * 80)
    print("COMPUTATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
