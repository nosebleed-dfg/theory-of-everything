#!/usr/bin/env python3
"""
120CELL_HODGE_PROPAGATOR — exact perturbative ratio R via Hodge Laplacian on 120-cell 1-forms
nos3bl33d

Discrete Maxwell operator. R determines <1-cos(plaquette)> = R/(2*beta). Coupling from spectral gap phi^(-4).
"""
import numpy as np
import time

PHI = (1 + np.sqrt(5)) / 2
ALPHA_QED = 1 / 137.035999084

print("=" * 78)
print("  HODGE LAPLACIAN & EXACT PERTURBATIVE RATIO")
print("=" * 78)

# Load graph
cache = np.load("C:/Users/funct/projects/120cell_cache.npz", allow_pickle=True)
adj = cache["adj"].astype(bool)
V = 600

edge_list = []
for i in range(V):
    for j in range(i + 1, V):
        if adj[i, j]:
            edge_list.append((i, j))
E = len(edge_list)
print(f"  V={V}, E={E}")

# Edge index
edge_index = {}
for idx, (i, j) in enumerate(edge_list):
    edge_index[(i, j)] = (idx, +1)
    edge_index[(j, i)] = (idx, -1)

# Build incidence matrix d_0: E x V
# d_0[e, v] = +1 if v is head, -1 if v is tail
d0 = np.zeros((E, V), dtype=np.float64)
for idx, (i, j) in enumerate(edge_list):
    d0[idx, i] = -1
    d0[idx, j] = +1

# Find faces
adj_sets = [set() for _ in range(V)]
for i, j in edge_list:
    adj_sets[i].add(j)
    adj_sets[j].add(i)

def canonical_face(cycle):
    n = len(cycle)
    best = None
    for s in range(n):
        fwd = tuple(cycle[(s + k) % n] for k in range(n))
        bwd = tuple(cycle[(s - k) % n] for k in range(n))
        if best is None or fwd < best:
            best = fwd
        if bwd < best:
            best = bwd
    return best

faces = set()
for u, v in edge_list:
    for w in adj_sets[v]:
        if w == u: continue
        for x in adj_sets[w]:
            if x == v or x == u: continue
            for y in adj_sets[x]:
                if y == w or y == v or y == u: continue
                if u in adj_sets[y]:
                    ch = True
                    if w in adj_sets[u] or x in adj_sets[u]: ch = False
                    if ch and (x in adj_sets[v] or y in adj_sets[v]): ch = False
                    if ch and y in adj_sets[w]: ch = False
                    if ch:
                        faces.add(canonical_face([u, v, w, x, y]))
face_list = sorted(faces)
F = len(face_list)
print(f"  F={F}")
assert F == 720

# Build boundary matrix d_1^T: F x E
# d1T[f, e] = +/- 1 if edge e is in face f (with orientation sign)
d1T = np.zeros((F, E), dtype=np.float64)
for f_idx, face in enumerate(face_list):
    for k in range(5):
        a, b = face[k], face[(k + 1) % 5]
        idx, sgn = edge_index[(a, b)]
        d1T[f_idx, idx] = sgn

# Verify d1T @ d0 = 0 (boundary of boundary = 0)
check = d1T @ d0
print(f"  d1T @ d0 = 0? max abs = {np.max(np.abs(check)):.2e}")

# ===========================================================================
# Hodge Laplacian on 1-forms (edges)
# ===========================================================================
print("\n  Computing Hodge-1 Laplacian...")
t0 = time.time()

# Delta_1 = d0 @ d0^T + d1T^T @ d1T
# d0 is E x V, so d0 @ d0^T is E x E ("co-curl" part)
# d1T is F x E, so d1T^T @ d1T is E x E ("curl" part)
cocurl = d0 @ d0.T    # E x E
curl = d1T.T @ d1T    # E x E
hodge1 = cocurl + curl

print(f"  Done in {time.time()-t0:.1f}s")

# Eigenvalues
print("  Computing eigenvalues...")
t0 = time.time()
eigs = np.sort(np.linalg.eigvalsh(hodge1))
print(f"  Done in {time.time()-t0:.1f}s")

n_zero = np.sum(np.abs(eigs) < 1e-8)
print(f"  Zero eigenvalues: {n_zero} (= b_1, first Betti number)")
print(f"  Expected for S^3: b_1 = 0")
print(f"  Smallest eigenvalue: {eigs[0]:.10f}")
print(f"  Largest eigenvalue: {eigs[-1]:.10f}")

# Spectrum summary
unique_eigs = []
tol = 1e-6
for e in eigs:
    if not unique_eigs or abs(e - unique_eigs[-1][0]) > tol:
        unique_eigs.append([e, 1])
    else:
        unique_eigs[-1][1] += 1

print(f"\n  Distinct Hodge-1 eigenvalues: {len(unique_eigs)}")
print(f"  {'eigenvalue':>14s} {'mult':>5s}")
for val, mult in unique_eigs[:15]:
    print(f"  {val:>14.8f} {mult:>5d}")
if len(unique_eigs) > 15:
    print(f"  ... ({len(unique_eigs) - 15} more)")

# ===========================================================================
# Exact perturbative ratio R
# ===========================================================================
print("\n  === EXACT PERTURBATIVE RATIO ===")

if n_zero == 0:
    # Delta_1 is invertible (b_1 = 0)
    print("  Delta_1 is invertible (b_1 = 0, as expected for S^3)")
    G = np.linalg.inv(hodge1)

    # The plaquette propagator: M = d1T @ G @ d1T^T (F x F matrix)
    # M_{fg} = <P_f * P_g> * beta  (in Gaussian approximation)
    M_plaq = d1T @ G @ d1T.T

    # R = (1/F) * Tr(M) = (1/F) * sum_f <P_f^2> * beta
    # so <1-cos(P)> ~ <P^2>/2 = R/(2*beta)
    R_exact = np.trace(M_plaq) / F

    print(f"  R = Tr(d1T @ Delta_1^(-1) @ d1T^T) / F")
    print(f"  R = {R_exact:.15f}")

    # Diagonal elements of M (should be uniform for vertex-transitive graph)
    M_diag = np.diag(M_plaq)
    print(f"  M diagonal: min={M_diag.min():.10f}, max={M_diag.max():.10f}")
    print(f"  (Uniform? {np.allclose(M_diag, M_diag[0], atol=1e-8)})")
else:
    print(f"  Delta_1 has {n_zero} zero modes -- using pseudoinverse")
    G = np.linalg.pinv(hodge1)
    M_plaq = d1T @ G @ d1T.T
    R_exact = np.trace(M_plaq) / F
    print(f"  R = {R_exact:.15f}")

# ===========================================================================
# Algebraic identification of R
# ===========================================================================
print("\n  === IDENTIFYING R ===")
sqrt5 = np.sqrt(5)

candidates = [
    ("1", 1.0),
    ("1/2", 0.5),
    ("601/720", 601/720),
    ("5/6", 5/6),
    ("4/5", 4/5),
    ("3/4", 3/4),
    ("2/3", 2/3),
    ("phi^(-4)", PHI**(-4)),
    ("phi^(-2)", PHI**(-2)),
    ("1/phi", 1/PHI),
    ("phi-1", PHI - 1),
    ("2-phi", 2 - PHI),
    ("(5-sqrt5)/5", (5 - sqrt5) / 5),
    ("(3-sqrt5)/2", (3 - sqrt5) / 2),
    ("sqrt5-2", sqrt5 - 2),
    ("(sqrt5-1)/2", (sqrt5 - 1) / 2),
    ("4/(3*sqrt5-1)", 4 / (3 * sqrt5 - 1)),
    ("(7-3*sqrt5)/2", (7 - 3 * sqrt5) / 2),
    ("5/(4*phi^2)", 5 / (4 * PHI**2)),
    ("1200/1440", 1200/1440),
    ("E/(E+F)", E/(E+F)),
    ("(E-V)/(E+F-V)", (E-V)/(E+F-V)),
]

print(f"  R_exact = {R_exact:.15f}")
print(f"  {'candidate':>25s} {'value':>18s} {'diff':>12s}")
for name, val in sorted(candidates, key=lambda x: abs(x[1] - R_exact)):
    diff = R_exact - val
    print(f"  {name:>25s} {val:>18.15f} {diff:>+12.2e}")

# ===========================================================================
# The coupling alpha
# ===========================================================================
print("\n  === COUPLING COMPUTATION ===")

beta_phi4 = PHI**4
alpha_R = R_exact / (4 * np.pi * beta_phi4)
print(f"  beta = phi^4 = {beta_phi4:.10f}")
print(f"  R = {R_exact:.10f}")
print(f"  alpha = R / (4*pi*phi^4) = {alpha_R:.15f}")
print(f"  1/alpha = {1/alpha_R:.6f}")
print(f"  QED: 1/alpha = 137.036")
print(f"  Ratio: {(1/alpha_R)/137.036:.6f}")
print()

# What if R enters differently?
# Option 1: alpha = 1/(4*pi*beta) with beta = phi^4/R
alpha_1 = 1 / (4 * np.pi * beta_phi4 / R_exact)
print(f"  If alpha = 1/(4*pi*(phi^4/R)) = R/(4*pi*phi^4) = {alpha_1:.10f}")
print(f"  1/alpha = {1/alpha_1:.6f}  (same as above)")

# Option 2: alpha = 1/(4*pi*beta) with beta = phi^4 (no R correction)
alpha_2 = 1 / (4 * np.pi * beta_phi4)
print(f"  If alpha = 1/(4*pi*phi^4) (no correction) = {alpha_2:.10f}")
print(f"  1/alpha = {1/alpha_2:.6f}")

# Option 3: What beta gives alpha = 1/137 exactly?
beta_137 = 1 / (4 * np.pi * ALPHA_QED)
print(f"\n  For alpha = 1/137.036:")
print(f"  Need beta = {beta_137:.10f}")
print(f"  phi^4 = {beta_phi4:.10f}")
print(f"  Ratio = {beta_137/beta_phi4:.10f}")

# Check: is 137.036 / (4*pi*phi^4) close to something involving R?
ratio_137 = beta_137 / beta_phi4
print(f"\n  beta_QED / phi^4 = {ratio_137:.10f}")
print(f"  1/R = {1/R_exact:.10f}")
print(f"  ratio * R = {ratio_137 * R_exact:.10f}")

# The "missing factor" to get from our geometric alpha to QED alpha:
missing = ALPHA_QED / alpha_R
print(f"\n  alpha_QED / alpha_geometric = {missing:.10f}")
print(f"  = {1/missing:.10f} (inverse)")

for name, val in [("phi^(-2)", PHI**(-2)), ("phi^(-1)", PHI**(-1)),
                   ("1/2", 0.5), ("2/3", 2/3), ("3/4", 3/4), ("4/5", 4/5),
                   ("5/6", 5/6), ("1/R", 1/R_exact),
                   ("R", R_exact), ("R^2", R_exact**2),
                   ("sqrt(R)", np.sqrt(R_exact))]:
    print(f"    vs {name:>10s} = {val:.6f}  diff = {abs(missing-val):.6f}")

# ===========================================================================
# Also check: cocurl and curl eigenvalues separately
# ===========================================================================
print("\n  === COCURL AND CURL SPECTRA ===")

eigs_cocurl = np.sort(np.linalg.eigvalsh(cocurl))
eigs_curl = np.sort(np.linalg.eigvalsh(curl))

# Cocurl (d0 @ d0^T): E-V+1 = 601 zero eigenvalues (gauge DOF)
n_zero_cc = np.sum(np.abs(eigs_cocurl) < 1e-8)
print(f"  Cocurl zero modes: {n_zero_cc} (expect E-V+1=601, gauge DOF)")
print(f"  Cocurl smallest nonzero: {eigs_cocurl[n_zero_cc]:.10f}")

# The nonzero eigenvalues of cocurl are the graph Laplacian eigenvalues!
# d0^T @ d0 = L (graph Laplacian, V x V)
# d0 @ d0^T has the same nonzero eigenvalues.
L = d0.T @ d0
eigs_L = np.sort(np.linalg.eigvalsh(L))
print(f"  Graph Laplacian mu_1 = {eigs_L[1]:.10f}")
print(f"  phi^(-4) = {PHI**(-4):.10f}")
print(f"  Match: {abs(eigs_L[1] - PHI**(-4)) < 1e-8}")

# Curl (d1T^T @ d1T): E-F+something zero eigenvalues
n_zero_curl = np.sum(np.abs(eigs_curl) < 1e-8)
print(f"  Curl zero modes: {n_zero_curl} (expect E-F+b_1={E}-{F}+0={E-F})")
print(f"  Curl smallest nonzero: {eigs_curl[n_zero_curl]:.10f}")

# The nonzero eigenvalues of d1T^T @ d1T are the face Laplacian eigenvalues
# d1T @ d1T^T (F x F) has the same nonzero spectrum.
face_lap = d1T @ d1T.T
eigs_fl = np.sort(np.linalg.eigvalsh(face_lap))
print(f"  Face Laplacian smallest nonzero: {eigs_fl[1]:.10f}")

# Check: smallest nonzero of face Laplacian vs phi^(-4)
print(f"\n  Is face Laplacian gap related to phi^(-4) = {PHI**(-4):.10f}?")
print(f"  Face gap = {eigs_fl[1]:.10f}")
ratio_fg = eigs_fl[1] / PHI**(-4)
print(f"  Face gap / phi^(-4) = {ratio_fg:.10f}")

print("\n  === HODGE-1 SPECTRUM AND phi^(-4) ===")
print(f"  Hodge-1 smallest eigenvalue: {eigs[0]:.10f}")
print(f"  phi^(-4) = {PHI**(-4):.10f}")
print(f"  Ratio: {eigs[0] / PHI**(-4):.10f}")

# Check if Hodge gap is phi^(-4) or a multiple
for name, val in [("phi^(-4)", PHI**(-4)), ("2*phi^(-4)", 2*PHI**(-4)),
                   ("phi^(-2)", PHI**(-2)), ("1", 1.0),
                   ("3-sqrt5", 3-sqrt5), ("4-sqrt5*phi", 4-sqrt5*PHI)]:
    print(f"  {name:>15s} = {val:.10f}  diff from gap = {abs(eigs[0]-val):.2e}")

print("\n" + "=" * 78)
print("  FINAL NUMBERS")
print("=" * 78)
print(f"""
  EXACT RESULTS (no MC uncertainty):
  - Spectral gap (graph Laplacian): mu_1 = phi^(-4) = {PHI**(-4):.15f}
  - Perturbative ratio R = {R_exact:.15f}
  - Hodge-1 gap = {eigs[0]:.15f}

  COUPLING WITH beta = phi^4:
  - alpha = R/(4*pi*phi^4) = {alpha_R:.15f}
  - 1/alpha = {1/alpha_R:.6f}

  COUPLING WITHOUT R correction (beta = phi^4):
  - alpha = 1/(4*pi*phi^4) = {alpha_2:.15f}
  - 1/alpha = {1/alpha_2:.6f}

  QED:
  - 1/alpha = 137.036

  CONCLUSION: Neither prescription gives 1/137.
  The 120-cell spectral gap determines a lattice scale,
  but the coupling constant requires additional physical input.
""")
