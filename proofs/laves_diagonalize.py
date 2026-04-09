"""
LAVES_DIAGONALIZE — builds the 8-vertex Laves graph unit cell, diagonalizes adjacency, checks for phi
nos3bl33d
"""
import numpy as np
from numpy.linalg import eigvalsh
import math
from collections import deque

phi = (1 + math.sqrt(5)) / 2

print("LAVES GRAPH: BUILD, DIAGONALIZE, CHECK FOR PHI")
print("=" * 60)

# The srs (Laves) net on 8 vertices per cubic unit cell.
# Bipartite: sublattice A = {0,1,2,3}, sublattice B = {4,5,6,7}
# Each vertex in A connects to exactly 3 in B and vice versa.
# 12 edges total = F (dodecahedral faces).
#
# Known adjacency (from Wells, Three-Dimensional Nets and Polyhedra):
edges = [(0,4),(0,5),(0,6), (1,4),(1,5),(1,7),
         (2,4),(2,6),(2,7), (3,5),(3,6),(3,7)]

A = np.zeros((8,8))
for i,j in edges:
    A[i,j] = 1; A[j,i] = 1

degrees = A.sum(axis=1).astype(int)
print(f"\nDegrees: {list(degrees)}")
print(f"All degree 3: {all(d==3 for d in degrees)}")
print(f"Edges: {len(edges)} = F = {12}")
print(f"Vertices: 8 = 2^d = {2**3}")
print()

# BIPARTITE CHECK
print("BIPARTITE CHECK:")
evals_A = np.sort(eigvalsh(A))
print(f"Adjacency eigenvalues: {[f'{v:.6f}' for v in evals_A]}")

# Bipartite iff spectrum symmetric around 0
evals_rounded = sorted([round(v, 10) for v in evals_A])
neg = sorted([-v for v in evals_rounded if v < -0.001])
pos = sorted([v for v in evals_rounded if v > 0.001])
symmetric = len(neg) == len(pos) and all(abs(n-p) < 1e-8 for n,p in zip(neg,pos))
print(f"Negative eigenvalues: {[f'{v:.4f}' for v in evals_A if v < -0.001]}")
print(f"Positive eigenvalues: {[f'{v:.4f}' for v in evals_A if v > 0.001]}")
print(f"Symmetric around 0: {symmetric}")
print(f"BIPARTITE: {symmetric}")

# Verify with 2-coloring
colors = [-1]*8
colors[0] = 0
queue = deque([0])
bipartite_ok = True
while queue:
    u = queue.popleft()
    for v in range(8):
        if A[u,v]:
            if colors[v] == -1:
                colors[v] = 1 - colors[u]
                queue.append(v)
            elif colors[v] == colors[u]:
                bipartite_ok = False
print(f"2-coloring check: {bipartite_ok}")
print(f"Coloring: {colors}")
print(f"Sublattice A: {[i for i in range(8) if colors[i]==0]}")
print(f"Sublattice B: {[i for i in range(8) if colors[i]==1]}")
print()

# GRAPH LAPLACIAN
L = np.diag(degrees.astype(float)) - A
evals_L = np.sort(eigvalsh(L))
print("LAPLACIAN EIGENVALUES:")
for i, ev in enumerate(evals_L):
    checks = []
    if abs(ev) < 1e-10: checks.append("= 0")
    for val, name in [(3-math.sqrt(5), "3-sqrt5=phi^-2"),
                       (2, "2=chi"), (3, "3=d"), (4, "4=d+1"),
                       (5, "5=p"), (6, "6=2d"),
                       (3+math.sqrt(5), "3+sqrt5"),
                       (1, "1"), (phi, "phi"), (1/phi, "1/phi"),
                       (phi**2, "phi^2"), (2*phi, "2phi"),
                       (2/phi, "2/phi")]:
        if abs(ev - val) < 0.001:
            checks.append(f"~ {name}")
    label = ", ".join(checks) if checks else ""
    print(f"  lambda_{i} = {ev:.6f}  {label}")
print()

# Spectral gap
gap = min(ev for ev in evals_L if ev > 0.001)
print(f"Spectral gap: {gap:.6f}")
print(f"3 - sqrt(5) = {3-math.sqrt(5):.6f} (dodecahedron)")
print(f"phi^(-2) = {phi**(-2):.6f}")
print(f"1 = {1}")
print(f"Match dodecahedron? {abs(gap - (3-math.sqrt(5))) < 0.01}")
print()

# Tr(L+)
nonzero_evals = [ev for ev in evals_L if ev > 0.001]
tr_Lp = sum(1/ev for ev in nonzero_evals)
print(f"Tr(L+) = {tr_Lp:.6f}")
print(f"Compare dodecahedron Tr(L+) = 137/15 = {137/15:.6f}")
print()

# Does phi appear anywhere?
print("DOES PHI APPEAR IN THE SPECTRUM?")
for ev in evals_L:
    if ev > 0.001:
        # Check if ev or 1/ev relates to phi
        ratio_phi = ev / phi
        ratio_phi2 = ev / phi**2
        if abs(ratio_phi - round(ratio_phi)) < 0.01:
            print(f"  {ev:.6f} = {round(ratio_phi)}*phi")
        if abs(ratio_phi2 - round(ratio_phi2)) < 0.01:
            print(f"  {ev:.6f} = {round(ratio_phi2)}*phi^2")
        if abs(ev - round(ev)) < 0.001:
            print(f"  {ev:.6f} = {round(ev)} (integer)")

print()

# Shortest cycle
print("SHORTEST CYCLE:")
def shortest_cycle(adj, n):
    min_c = float('inf')
    for start in range(n):
        dist = [-1]*n
        parent = [-1]*n
        dist[start] = 0
        q = deque([start])
        while q:
            u = q.popleft()
            for v in range(n):
                if adj[u,v] > 0:
                    if dist[v] == -1:
                        dist[v] = dist[u] + 1
                        parent[v] = u
                        q.append(v)
                    elif parent[u] != v:
                        cycle = dist[u] + dist[v] + 1
                        min_c = min(min_c, cycle)
    return min_c

girth = shortest_cycle(A, 8)
print(f"  Girth: {girth}")
print(f"  Note: 8-vertex unit cell may create short cycles from PBC")
print(f"  True srs girth = 10 = 2p (requires larger cell)")
print()

# THE VERDICT
print("=" * 60)
print("VERDICT")
print("=" * 60)
print()
print(f"  Bipartite: YES (verified by 2-coloring and spectrum)")
print(f"  Degree: 3 = d (all vertices)")
print(f"  Edges per cell: 12 = F")
print(f"  Vertices per cell: 8 = 2^d")
print(f"  Spectral gap: {gap:.6f}")
print(f"  Laplacian eigenvalues: {[f'{v:.4f}' for v in evals_L]}")
print()

if abs(gap - (3 - math.sqrt(5))) < 0.01:
    print("  PHI APPEARS: spectral gap = 3-sqrt(5) = phi^(-2)")
    print("  The dodecahedron's spectral gap IS the Laves graph's gap!")
elif abs(gap - 1) < 0.01:
    print("  Gap = 1 (integer, not golden)")
    print("  Phi does NOT appear directly in the 8-vertex cell spectrum.")
    print("  May appear in the INFINITE Laves graph (band structure).")
else:
    print(f"  Gap = {gap:.4f}")
    print(f"  Need to check if this relates to phi.")
