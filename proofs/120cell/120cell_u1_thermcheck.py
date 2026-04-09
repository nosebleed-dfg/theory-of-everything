#!/usr/bin/env python3
"""
120CELL_U1_THERMCHECK — cold vs hot start thermalization check at high beta; higher stats at beta 1-5
nos3bl33d
"""
import numpy as np
from scipy.special import i0, i1
import time

# Load precomputed lattice data
cache = np.load("C:/Users/funct/projects/120cell_cache.npz", allow_pickle=True)
adj_matrix = cache["adj"].astype(bool)
N_VERTS = 600

edge_list = []
for i in range(N_VERTS):
    for j in range(i + 1, N_VERTS):
        if adj_matrix[i, j]:
            edge_list.append((i, j))
N_EDGES = len(edge_list)

edge_index = {}
for idx, (i, j) in enumerate(edge_list):
    edge_index[(i, j)] = (idx, +1)
    edge_index[(j, i)] = (idx, -1)

adj_sets = [set() for _ in range(N_VERTS)]
for i, j in edge_list:
    adj_sets[i].add(j)
    adj_sets[j].add(i)

# Find faces
def canonical_face(cycle):
    n = len(cycle)
    best = None
    for s in range(n):
        fwd = tuple(cycle[(s+k)%n] for k in range(n))
        bwd = tuple(cycle[(s-k)%n] for k in range(n))
        if best is None or fwd < best: best = fwd
        if bwd < best: best = bwd
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
                    if ch: faces.add(canonical_face([u, v, w, x, y]))

face_list = sorted(faces)
N_FACES = len(face_list)
assert N_FACES == 720

face_edge_idx = np.zeros((N_FACES, 5), dtype=np.int32)
face_edge_sgn = np.zeros((N_FACES, 5), dtype=np.float64)
for f_idx, face in enumerate(face_list):
    for k in range(5):
        a, b = face[k], face[(k+1)%5]
        idx, sgn = edge_index[(a, b)]
        face_edge_idx[f_idx, k] = idx
        face_edge_sgn[f_idx, k] = sgn

edge_faces_eidx = np.zeros((N_EDGES, 3, 5), dtype=np.int32)
edge_faces_esgn = np.zeros((N_EDGES, 3, 5), dtype=np.float64)
efl = [[] for _ in range(N_EDGES)]
for f_idx in range(N_FACES):
    for k in range(5):
        efl[face_edge_idx[f_idx, k]].append(f_idx)
for e in range(N_EDGES):
    for k, f in enumerate(efl[e]):
        edge_faces_eidx[e, k, :] = face_edge_idx[f, :]
        edge_faces_esgn[e, k, :] = face_edge_sgn[f, :]

def compute_avg_cos(theta):
    phases = np.sum(face_edge_sgn * theta[face_edge_idx], axis=1)
    return np.mean(np.cos(phases))

def sweep(theta, beta, rng):
    proposals = rng.uniform(-np.pi, np.pi, N_EDGES)
    randoms = rng.random(N_EDGES)
    acc = 0
    for e in range(N_EDGES):
        old_phases = np.sum(edge_faces_esgn[e] * theta[edge_faces_eidx[e]], axis=1)
        old_cs = np.sum(np.cos(old_phases))
        old_val = theta[e]
        theta[e] = old_val + proposals[e]
        new_phases = np.sum(edge_faces_esgn[e] * theta[edge_faces_eidx[e]], axis=1)
        new_cs = np.sum(np.cos(new_phases))
        dS = -beta * (new_cs - old_cs)
        if dS <= 0 or randoms[e] < np.exp(-dS):
            acc += 1
        else:
            theta[e] = old_val
    return acc / N_EDGES

print("=" * 78)
print("  THERMALIZATION CHECK: COLD vs HOT START")
print("=" * 78)
print()

PHI = (1 + np.sqrt(5)) / 2

# Test at a few high-beta values
for beta in [5.0, 10.0, 20.0, 50.0]:
    t0 = time.time()
    rng_c = np.random.default_rng(42)
    rng_h = np.random.default_rng(42)

    theta_cold = np.zeros(N_EDGES)
    theta_hot = rng_h.uniform(-np.pi, np.pi, N_EDGES)

    cold_trace = []
    hot_trace = []

    n_sweeps = 500

    for s in range(n_sweeps):
        sweep(theta_cold, beta, rng_c)
        sweep(theta_hot, beta, rng_h)

        if s % 10 == 0 or s < 20:
            cold_trace.append((s, compute_avg_cos(theta_cold)))
            hot_trace.append((s, compute_avg_cos(theta_hot)))

    print(f"  beta = {beta:.1f}:")
    print(f"    Cold start trace: {cold_trace[0][1]:.6f} -> {cold_trace[-1][1]:.6f}")
    print(f"    Hot  start trace: {hot_trace[0][1]:.6f} -> {hot_trace[-1][1]:.6f}")
    print(f"    Difference at end: {abs(cold_trace[-1][1] - hot_trace[-1][1]):.6f}")

    converged = abs(cold_trace[-1][1] - hot_trace[-1][1]) < 0.005
    print(f"    Converged: {'YES' if converged else 'NO'}")
    print(f"    Time: {time.time()-t0:.0f}s")

    # Show some trajectory
    print(f"    Cold: ", end="")
    for s, c in cold_trace[:5]:
        print(f"[{s}]{c:.4f} ", end="")
    print(f"... [{cold_trace[-1][0]}]{cold_trace[-1][1]:.4f}")
    print(f"    Hot:  ", end="")
    for s, c in hot_trace[:5]:
        print(f"[{s}]{c:.4f} ", end="")
    print(f"... [{hot_trace[-1][0]}]{hot_trace[-1][1]:.4f}")
    print()

# Now do higher-stats run at the CROSSOVER REGION
print("=" * 78)
print("  HIGH-STATS CROSSOVER REGION (beta 1-5)")
print("=" * 78)

beta_fine = [1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
n_therm_hi = 500
n_meas_hi = 2000

for beta in beta_fine:
    t0 = time.time()
    rng = np.random.default_rng(42)
    theta = np.zeros(N_EDGES)

    for _ in range(n_therm_hi):
        sweep(theta, beta, rng)

    measurements = []
    for m in range(n_meas_hi):
        sweep(theta, beta, rng)
        measurements.append(compute_avg_cos(theta))

    measurements = np.array(measurements)
    mean_c = np.mean(measurements)
    # Jackknife
    n_blocks = 20
    bs = len(measurements) // n_blocks
    full = np.mean(measurements)
    jacks = np.array([np.mean(np.delete(measurements, slice(b*bs, (b+1)*bs)))
                      for b in range(n_blocks)])
    err = np.sqrt((n_blocks-1)/n_blocks * np.sum((jacks - full)**2))

    sc = i1(beta) / i0(beta)
    deficit = 1 - mean_c
    pert_ratio = deficit / (1/(2*beta)) if deficit > 0 else 0

    print(f"  beta={beta:>4.1f}: <cos>={mean_c:.7f}+/-{err:.7f}  "
          f"SC={sc:.7f}  diff={mean_c-sc:+.7f}  "
          f"R={pert_ratio:.4f}  [{time.time()-t0:.0f}s]",
          flush=True)

print()
print("  R = perturbative ratio = (1-<cos>) / (1/(2*beta))")
print("  R = 1 means free U(1); deviations = lattice corrections")
print()

# What does R converge to at weak coupling?
print("  At weak coupling (large beta), R should converge to a")
print("  lattice-specific constant determined by the graph topology.")
print("  This constant encodes how many effective DOF fluctuate")
print("  around each plaquette.")
print()
print(f"  For the 120-cell:")
print(f"    V=600, E=1200, F=720")
print(f"    Gauge DOF = E - V + 1 = 601 (after gauge fixing)")
print(f"    Plaquettes = 720")
print(f"    DOF per plaquette = 601/720 = {601/720:.6f}")
print(f"    For a square lattice: DOF/plaq = 1 (in 4D)")
print(f"    Expected R ~ DOF_per_plaq = {601/720:.4f}")
