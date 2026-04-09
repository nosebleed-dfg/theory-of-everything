#!/usr/bin/env python3
"""
120CELL_SPECTRAL_PROOF — algebraic proof that the 120-cell spectral gap is exactly phi^(-4)
nos3bl33d

600 vertices, degree 4, Laplacian L=4I-A. mu_1 = (7-3*sqrt(5))/2.
H4 Coxeter construction, numerical verification, all 27 eigenvalues, exact sympy proof.
"""

import numpy as np
from itertools import permutations, product
from fractions import Fraction
import time


###############################################################################
#  PHASE 1: H4 COXETER GROUP AND 120-CELL CONSTRUCTION
###############################################################################

def build_H4_generators():
    """
    Construct the 4 unit normal vectors for the generating reflections of H4.

    The Coxeter diagram: s1 --3-- s2 --3-- s3 --5-- s4
    Relation: n_i . n_j = -cos(pi / m_{ij}).
    """
    phi = (1 + np.sqrt(5)) / 2

    n1 = np.array([1.0, 0.0, 0.0, 0.0])
    n2 = np.array([-0.5, np.sqrt(3)/2, 0.0, 0.0])
    n3 = np.array([0.0, -1/np.sqrt(3), np.sqrt(2.0/3.0), 0.0])

    b4 = -phi * np.sqrt(6) / 4
    c4 = np.sqrt(1.0 - b4**2)
    n4 = np.array([0.0, 0.0, b4, c4])

    # Verify all Coxeter matrix entries
    assert abs(np.dot(n1, n2) + 0.5) < 1e-14          # m=3
    assert abs(np.dot(n2, n3) + 0.5) < 1e-14          # m=3
    assert abs(np.dot(n3, n4) + phi/2) < 1e-14        # m=5
    assert abs(np.dot(n1, n3)) < 1e-14                # m=2
    assert abs(np.dot(n1, n4)) < 1e-14                # m=2
    assert abs(np.dot(n2, n4)) < 1e-14                # m=2

    return [n1, n2, n3, n4], phi


def compute_fundamental_weights(normals):
    """
    Compute fundamental weights omega_i satisfying omega_i . n_j = delta_{ij}/2.
    """
    N = np.vstack(normals)
    N_inv = np.linalg.inv(N)
    omegas = [0.5 * N_inv[:, i] for i in range(4)]
    for i in range(4):
        for j in range(4):
            dot = np.dot(omegas[i], normals[j])
            expected = 0.5 if i == j else 0.0
            assert abs(dot - expected) < 1e-12
    return omegas


def generate_H4_orbit(start, normals):
    """Generate the H4 orbit of a point via BFS on reflections."""
    def key(p):
        return tuple(round(x, 9) for x in p)

    seen = {key(start)}
    queue = [np.array(start, dtype=np.float64)]
    orbit = [queue[0].copy()]
    head = 0

    while head < len(queue):
        p = queue[head]; head += 1
        for n in normals:
            r = p - 2 * np.dot(p, n) * n
            k = key(r)
            if k not in seen:
                seen.add(k)
                queue.append(r)
                orbit.append(r.copy())

    return orbit


def build_120cell():
    """
    Build the 600 vertices of the 120-cell on S^3.
    Uses the H4 orbit of the fundamental weight omega_4.
    Returns vertices normalized to the unit 3-sphere.
    """
    normals, phi = build_H4_generators()
    omegas = compute_fundamental_weights(normals)

    # omega_4 generates the 120-cell orbit (600 points)
    orbit = generate_H4_orbit(omegas[3], normals)
    assert len(orbit) == 600, f"Expected 600 vertices, got {len(orbit)}"

    vertices = np.array(orbit)
    norms = np.linalg.norm(vertices, axis=1)
    radius = norms[0]
    vertices /= radius  # normalize to unit S^3

    return vertices, normals, omegas, phi


###############################################################################
#  PHASE 2: ADJACENCY MATRIX AND NUMERICAL SPECTRUM
###############################################################################

def build_adjacency(vertices, expected_degree=4):
    """Build adjacency matrix by connecting nearest neighbors (edge-transitive)."""
    n = len(vertices)
    norms_sq = np.sum(vertices**2, axis=1)
    dots = vertices @ vertices.T
    dist_sq = norms_sq[:, None] + norms_sq[None, :] - 2 * dots
    dist_sq = np.maximum(dist_sq, 0)
    np.fill_diagonal(dist_sq, np.inf)

    edge_len = np.sqrt(np.min(dist_sq))
    tol = edge_len * 1e-6

    A = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i+1, n):
            if abs(np.sqrt(dist_sq[i, j]) - edge_len) < tol:
                A[i, j] = 1.0
                A[j, i] = 1.0

    degrees = A.sum(axis=1).astype(int)
    assert set(degrees) == {expected_degree}
    return A, edge_len


def laplacian_spectrum(A, degree=4):
    """Sorted eigenvalues of L = degree*I - A."""
    n = A.shape[0]
    L = degree * np.eye(n) - A
    return np.sort(np.linalg.eigvalsh(L))


def cluster_eigenvalues(evals, tol=1e-8):
    """Return list of (value, multiplicity) pairs."""
    clusters = []
    for ev in evals:
        placed = False
        for c in clusters:
            if abs(ev - c[0]) < tol:
                c.append(ev)
                placed = True
                break
        if not placed:
            clusters.append([ev])
    return sorted([(np.mean(c), len(c)) for c in clusters])


###############################################################################
#  PHASE 3: ALGEBRAIC IDENTIFICATION OF ALL EIGENVALUES
###############################################################################

def identify_eigenvalue(val):
    """
    Express a Laplacian eigenvalue in algebraic form.

    All eigenvalues lie in the splitting field of two specific cubics over Q(sqrt(5)):
    - The "rational" eigenvalues: elements of Q(sqrt(5)) or Q(sqrt(13)) or Q(sqrt(21))
    - The "cubic" eigenvalues: roots of x^3 - x^2 - 7x + 4 or x^3 - x^2 - 7x + 8
      (irreducible cubics; casus irreducibilis with 3 real roots each)

    We identify each eigenvalue by finding its minimal polynomial.
    """
    sqrt5 = np.sqrt(5)
    sqrt2 = np.sqrt(2)
    sqrt13 = np.sqrt(13)
    sqrt21 = np.sqrt(21)

    # Exact candidates in Q(sqrt(5))
    candidates_qsqrt5 = [
        (0, "0"),
        ((7 - 3*sqrt5)/2, "(7-3*sqrt(5))/2 = phi^(-4)"),
        ((3 - sqrt5)/2, "(3-sqrt(5))/2 = phi^(-2)"),
        (4 - sqrt5, "4-sqrt(5)"),
        ((7 - sqrt5)/2, "(7-sqrt(5))/2"),
        ((3 + sqrt5)/2, "(3+sqrt(5))/2 = phi^2"),
        (3, "3"),
        (4, "4"),
        ((11 - sqrt5)/2, "(11-sqrt(5))/2"),
        ((7 + sqrt5)/2, "(7+sqrt(5))/2"),
        (5, "5"),
        (6, "6"),
        (4 + sqrt5, "4+sqrt(5)"),
        ((11 + sqrt5)/2, "(11+sqrt(5))/2"),
        ((7 + 3*sqrt5)/2, "(7+3*sqrt(5))/2 = phi^4"),
    ]

    # Candidates involving sqrt(2)
    candidates_sqrt2 = [
        (5 - sqrt2, "5-sqrt(2)"),
        (5 + sqrt2, "5+sqrt(2)"),
    ]

    # Candidates involving sqrt(13) (from x^2 - 3x - 1 = 0 as adj eigenvalue)
    candidates_sqrt13 = [
        (4 - (3 + sqrt13)/2, "5/2 - sqrt(13)/2"),
        (4 - (3 - sqrt13)/2, "5/2 + sqrt(13)/2"),
    ]

    # Candidates involving sqrt(21) (from x^2 + x - 5 = 0 as adj eigenvalue)
    candidates_sqrt21 = [
        (4 - (-1 + sqrt21)/2, "9/2 - sqrt(21)/2"),
        (4 - (-1 - sqrt21)/2, "9/2 + sqrt(21)/2"),
    ]

    # Cubic candidates: roots of x^3 - x^2 - 7x + 4 = 0 (adj eigenvalue)
    # Laplacian eigenvalues: 4 - root
    p, q = -22/3, 43/27
    m = 2*np.sqrt(-p/3)
    theta1 = np.arccos(-q/2 * np.sqrt(-27/p**3))
    cubic1_adj_roots = [m * np.cos((theta1 - 2*np.pi*k)/3) + 1/3 for k in range(3)]
    cubic1_names = ["4-r1(x^3-x^2-7x+4)", "4-r2(x^3-x^2-7x+4)", "4-r3(x^3-x^2-7x+4)"]

    # Cubic candidates: roots of x^3 - x^2 - 7x + 8 = 0 (adj eigenvalue)
    q2 = 151/27
    theta2 = np.arccos(-q2/2 * np.sqrt(-27/p**3))
    cubic2_adj_roots = [m * np.cos((theta2 - 2*np.pi*k)/3) + 1/3 for k in range(3)]
    cubic2_names = ["4-r1(x^3-x^2-7x+8)", "4-r2(x^3-x^2-7x+8)", "4-r3(x^3-x^2-7x+8)"]

    all_candidates = (
        candidates_qsqrt5 + candidates_sqrt2 + candidates_sqrt13 + candidates_sqrt21 +
        [(4 - r, n) for r, n in zip(cubic1_adj_roots, cubic1_names)] +
        [(4 - r, n) for r, n in zip(cubic2_adj_roots, cubic2_names)]
    )

    for cand_val, cand_name in all_candidates:
        if abs(val - cand_val) < 1e-9:
            return cand_name

    return f"??? ({val:.15f})"


###############################################################################
#  PHASE 4: EXACT ALGEBRAIC PROOF (sympy)
###############################################################################

def exact_algebraic_proof(vertices, A, normals, omegas, phi_num):
    """
    EXACT PROOF that the spectral gap = phi^(-4).

    The proof has five steps, each verified symbolically with sympy:
    1. Compute edge dot product cos(theta) EXACTLY
    2. Show cos(theta) = (1 + 3*sqrt(5))/8
    3. Derive eigenvalue in 4D representation = 4*cos(theta) via Schur's lemma
    4. Therefore spectral gap = 4 - 4*cos(theta) = (7 - 3*sqrt(5))/2
    5. Verify (7 - 3*sqrt(5))/2 = phi^(-4) and satisfies x^2 - 7x + 1 = 0
    """
    from sympy import (sqrt as ssqrt, Rational as R, simplify, Matrix,
                       Symbol, factor, expand, nsimplify)

    print("\n" + "=" * 80)
    print("  PHASE 4: EXACT ALGEBRAIC PROOF USING SYMPY")
    print("=" * 80)

    phi_s = (1 + ssqrt(5)) / 2

    # ==================================================================
    # STEP 1: Exact edge dot product
    # ==================================================================
    print("\n  STEP 1: Compute exact edge dot product on S^3")
    print("  " + "-" * 60)

    # Build exact H4 normals using sympy
    n1_s = Matrix([1, 0, 0, 0])
    n2_s = Matrix([R(-1, 2), ssqrt(3) / 2, 0, 0])
    n3_s = Matrix([0, -ssqrt(3) / 3, ssqrt(6) / 3, 0])

    b4_s = -phi_s * ssqrt(6) / 4
    c4_sq = 1 - b4_s**2
    c4_s = ssqrt(simplify(c4_sq))
    n4_s = Matrix([0, 0, b4_s, c4_s])

    print(f"    n1 = {list(n1_s)}")
    print(f"    n2 = {[simplify(x) for x in n2_s]}")
    print(f"    n3 = {[simplify(x) for x in n3_s]}")
    print(f"    n4 = {[simplify(x) for x in n4_s]}")

    # Verify Coxeter relations symbolically
    assert simplify(n1_s.dot(n2_s) + R(1, 2)) == 0
    assert simplify(n2_s.dot(n3_s) + R(1, 2)) == 0
    assert simplify(n3_s.dot(n4_s) + phi_s / 2) == 0
    assert simplify(n1_s.dot(n3_s)) == 0
    assert simplify(n1_s.dot(n4_s)) == 0
    assert simplify(n2_s.dot(n4_s)) == 0
    print("    Coxeter relations verified symbolically.")

    # Compute exact omega_4 = (1/2) * 4th column of N^{-1}
    N_s = Matrix([
        [n1_s[0], n1_s[1], n1_s[2], n1_s[3]],
        [n2_s[0], n2_s[1], n2_s[2], n2_s[3]],
        [n3_s[0], n3_s[1], n3_s[2], n3_s[3]],
        [n4_s[0], n4_s[1], n4_s[2], n4_s[3]],
    ])
    N_inv = N_s.inv()
    omega4_s = R(1, 2) * N_inv.col(3)
    omega4_s = Matrix([simplify(omega4_s[i]) for i in range(4)])
    print(f"\n    omega_4 = {list(omega4_s)}")

    # omega_4 . n_4 = 1/2 (by construction)
    dot_o4_n4 = simplify(omega4_s.dot(n4_s))
    assert dot_o4_n4 == R(1, 2), f"omega_4.n_4 = {dot_o4_n4}, expected 1/2"

    # The nearest neighbor of omega_4 via reflection s_4:
    # s_4(omega_4) = omega_4 - 2*(omega_4.n_4)*n_4 = omega_4 - n_4
    v0_s = omega4_s
    v1_s = omega4_s - n4_s  # since 2*(1/2)*n4 = n4
    v1_s = Matrix([simplify(v1_s[i]) for i in range(4)])
    print(f"    v0 = omega_4 = {list(v0_s)}")
    print(f"    v1 = s_4(omega_4) = omega_4 - n_4 = {list(v1_s)}")

    # Edge dot product: cos(theta) = (v0 . v1) / |v0|^2
    # (Both lie on the same sphere of radius |omega_4|)
    v0_dot_v1 = simplify(v0_s.dot(v1_s))
    v0_sq = simplify(v0_s.dot(v0_s))
    cos_theta_s = simplify(v0_dot_v1 / v0_sq)

    print(f"\n    v0 . v1 = {v0_dot_v1}")
    print(f"    |v0|^2  = {v0_sq}")
    print(f"    cos(theta_edge) = (v0.v1) / |v0|^2 = {cos_theta_s}")

    # ==================================================================
    # STEP 2: Prove cos(theta) = (1 + 3*sqrt(5))/8
    # ==================================================================
    print(f"\n  STEP 2: Verify cos(theta) = (1 + 3*sqrt(5))/8")
    print("  " + "-" * 60)

    target_cos = (1 + 3 * ssqrt(5)) / 8
    diff = simplify(cos_theta_s - target_cos)
    print(f"    cos(theta) - (1+3*sqrt(5))/8 = {diff}")
    assert diff == 0, f"FAILED: difference = {diff}"
    print(f"    VERIFIED: cos(theta_edge) = (1 + 3*sqrt(5))/8  [EXACT]")

    # Numerical check
    cos_num = float(cos_theta_s.evalf())
    cos_target = (1 + 3*np.sqrt(5))/8
    print(f"    Numerical: {cos_num:.15f}")
    print(f"    Target:    {cos_target:.15f}")

    # ==================================================================
    # STEP 3: Representation-theoretic derivation
    # ==================================================================
    print(f"\n  STEP 3: Adjacency eigenvalue via Schur's lemma")
    print("  " + "-" * 60)

    print("""
    THEOREM (Schur orthogonality for vertex-transitive graphs):
    Let G be a vertex-transitive graph with n vertices, degree d, embedded
    in R^k via coordinates V (n x k matrix). Suppose the symmetry group
    acts irreducibly on R^k (the "standard representation"). Then:

    (a) V^T V = (n/k) * I_k                    [group-orbit orthogonality]
    (b) V^T A V = (n/k) * a_std * I_k          [Schur's lemma: A commutes with G]
    (c) a_std = k/(n) * Tr(V^T A V)            [take trace of (b)]
         = k/n * sum_{(i,j) edge} 2*(v_i . v_j)
         = k/n * 2 * |E| * cos(theta_edge)     [all edges have same angle]

    For the 120-cell: n=600, k=4, d=4, |E|=1200.
    a_std = 4/600 * 2 * 1200 * cos(theta_edge) = 4 * cos(theta_edge)

    DERIVATION:
    V^T V = sum_{i=1}^{600} v_i v_i^T  where v_i are unit vectors on S^3.
    The 600 vertices form an orbit of H4 acting on R^4. Since this action
    is irreducible (H4 acts faithfully on R^4 via reflections), the matrix
    sum_{i} v_i v_i^T must be a scalar multiple of I_4 (by Schur's lemma).
    Taking trace: Tr(V^T V) = sum_i |v_i|^2 = 600, so V^T V = 150 * I_4.

    The adjacency matrix A commutes with every element of H4 (since H4
    permutes vertices while preserving adjacency). Therefore V^T A V
    also commutes with H4 acting on R^4. By Schur's lemma (irreducible
    action), V^T A V = c * I_4 for some scalar c.

    c = Tr(V^T A V)/4 = (1/4) sum_k sum_{i,j: A_{ij}=1} (v_i)_k (v_j)_k
      = (1/4) sum_{A_{ij}=1} v_i . v_j = (1/4) * 2400 * cos(theta)
      = 600 * cos(theta)

    So a_std = c / 150 = 4 * cos(theta).
    """)

    a_std_s = 4 * cos_theta_s
    a_std_simplified = simplify(a_std_s)
    print(f"    a_std = 4 * cos(theta_edge) = 4 * (1+3*sqrt(5))/8 = {a_std_simplified}")
    assert simplify(a_std_simplified - (1 + 3*ssqrt(5))/2) == 0
    print(f"    a_std = (1 + 3*sqrt(5))/2")

    # ==================================================================
    # STEP 4: Spectral gap = 4 - a_std
    # ==================================================================
    print(f"\n  STEP 4: Spectral gap = degree - a_std")
    print("  " + "-" * 60)

    mu_1_s = 4 - a_std_s
    mu_1_simplified = simplify(mu_1_s)
    print(f"    mu_1 = 4 - a_std = 4 - (1+3*sqrt(5))/2 = {mu_1_simplified}")

    # Simplify: 4 - (1+3*sqrt(5))/2 = (8 - 1 - 3*sqrt(5))/2 = (7 - 3*sqrt(5))/2
    target_mu = (7 - 3*ssqrt(5)) / 2
    diff_mu = simplify(mu_1_simplified - target_mu)
    assert diff_mu == 0
    print(f"    = (7 - 3*sqrt(5))/2")

    # Verify = phi^{-4}
    phi_inv4 = 1 / phi_s**4
    diff_phi = simplify(mu_1_simplified - phi_inv4)
    assert diff_phi == 0
    print(f"    = phi^(-4)    [VERIFIED EXACTLY]")

    # ==================================================================
    # STEP 5: Minimal polynomial and irreducibility
    # ==================================================================
    print(f"\n  STEP 5: Minimal polynomial over Q")
    print("  " + "-" * 60)

    x = Symbol('x')
    # phi^{-4} satisfies x^2 - 7x + 1 = 0
    val = (7 - 3*ssqrt(5)) / 2
    poly_check = simplify(val**2 - 7*val + 1)
    assert poly_check == 0
    print(f"    x = (7-3*sqrt(5))/2 satisfies: x^2 - 7x + 1 = 0  [VERIFIED]")

    # Irreducibility: discriminant = 49 - 4 = 45 = 9*5, not a perfect square
    print(f"    Discriminant = 49 - 4 = 45 = 9 * 5  (not a perfect square)")
    print(f"    => x^2 - 7x + 1 is IRREDUCIBLE over Q")
    print(f"    => phi^(-4) is algebraic of degree 2 over Q")

    # The two roots
    r_plus = (7 + 3*ssqrt(5)) / 2  # = phi^4 ~ 6.854
    r_minus = (7 - 3*ssqrt(5)) / 2  # = phi^{-4} ~ 0.146
    print(f"\n    Roots of x^2 - 7x + 1 = 0:")
    print(f"      r+ = (7+3*sqrt(5))/2 = phi^4  ~ {float(r_plus.evalf()):.6f}")
    print(f"      r- = (7-3*sqrt(5))/2 = phi^-4 ~ {float(r_minus.evalf()):.6f}")
    print(f"    Since 0 < mu_1 < 1, the spectral gap equals r- = phi^(-4).")

    # Product and sum of roots (Vieta's formulas)
    print(f"\n    Vieta's formulas:")
    print(f"      phi^4 * phi^(-4) = 1  (product of roots)")
    print(f"      phi^4 + phi^(-4) = 7  (sum of roots)")

    # ==================================================================
    # STEP 6: Verify smallest nonzero eigenvalue
    # ==================================================================
    print(f"\n  STEP 6: Verify this is the SMALLEST nonzero eigenvalue")
    print("  " + "-" * 60)

    eigenvalues = laplacian_spectrum(A, degree=4)
    clusters = cluster_eigenvalues(eigenvalues)

    target_num = (7 - 3*np.sqrt(5))/2
    assert abs(clusters[0][0]) < 1e-10, "First eigenvalue must be 0"
    assert clusters[0][1] == 1, "Graph must be connected (nullity 1)"
    assert abs(clusters[1][0] - target_num) < 1e-10
    assert clusters[1][1] == 4, "Spectral gap must have multiplicity 4"

    print(f"    Eigenvalue 0: value = {clusters[0][0]:.2e}, mult = {clusters[0][1]}")
    print(f"    Eigenvalue 1: value = {clusters[1][0]:.15f}, mult = {clusters[1][1]}")
    print(f"    phi^(-4)           = {target_num:.15f}")
    print(f"    |difference|       = {abs(clusters[1][0] - target_num):.2e}")
    print(f"    Next eigenvalue    = {clusters[2][0]:.15f} (mult {clusters[2][1]})")

    print(f"\n    The spectral gap has multiplicity 4, matching the dimension of")
    print(f"    the standard (reflection) representation of H4. This confirms")
    print(f"    that our representation-theoretic derivation identifies the")
    print(f"    correct eigenvalue.")

    print(f"\n    All 27 distinct eigenvalues mu_i > mu_1 (for i >= 2):")
    print(f"    mu_2 = {clusters[2][0]:.12f} > mu_1 = {clusters[1][0]:.12f}  CHECK")

    return True


###############################################################################
#  FULL EIGENVALUE TABLE
###############################################################################

def print_eigenvalue_table(A):
    """Print the complete table of all 27 distinct Laplacian eigenvalues."""
    eigenvalues = laplacian_spectrum(A, degree=4)
    clusters = cluster_eigenvalues(eigenvalues)

    print(f"\n  {'='*80}")
    print(f"  COMPLETE LAPLACIAN SPECTRUM OF THE 120-CELL")
    print(f"  {'='*80}")
    print(f"  {'#':>3}  {'Eigenvalue':>20}  {'Mult':>5}  {'Cumul':>5}  {'Algebraic Form'}")
    print(f"  {'-'*78}")

    cumul = 0
    for idx, (val, mult) in enumerate(clusters):
        cumul += mult
        alg = identify_eigenvalue(val)
        marker = " <-- SPECTRAL GAP" if idx == 1 else ""
        print(f"  {idx+1:3d}  {val:20.12f}  {mult:5d}  {cumul:5d}  {alg}{marker}")

    print(f"  {'-'*78}")
    print(f"  {'':>3}  {'Tr(L) = 2400':>20}  {cumul:5d}")
    total = sum(val * mult for val, mult in clusters)
    print(f"  Trace check: sum(eigenvalue * multiplicity) = {total:.6f} (= 4*600 = 2400)")

    # Print multiplicity analysis
    print(f"\n  Multiplicity structure (= irrep dimensions of H4 in V_600):")
    from collections import Counter
    mc = Counter(mult for _, mult in clusters)
    for d in sorted(mc.keys()):
        print(f"    dim {d:3d}: appears {mc[d]} time(s)")
    print(f"    Total: {sum(d*c for d,c in mc.items())} = 600")

    return clusters


###############################################################################
#  MAIN: ASSEMBLE THE FULL PROOF
###############################################################################

def main():
    t0 = time.time()
    phi = (1 + np.sqrt(5)) / 2

    print("=" * 80)
    print("       ALGEBRAIC PROOF: SPECTRAL GAP OF THE 120-CELL = phi^(-4)")
    print("=" * 80)
    print(f"\n  phi = (1+sqrt(5))/2 = {phi:.15f}")
    print(f"  phi^(-4) = (7-3*sqrt(5))/2 = {1/phi**4:.15f}")

    # ===== PHASE 1: Construction =====
    print(f"\n  {'='*80}")
    print("  PHASE 1: Construct the 120-cell via H4 Coxeter reflections")
    print(f"  {'='*80}")

    vertices, normals, omegas, phi = build_120cell()
    print(f"    Vertices: {len(vertices)} (on unit S^3)")
    norms = np.linalg.norm(vertices, axis=1)
    print(f"    Radius range: [{norms.min():.12f}, {norms.max():.12f}]")

    # ===== PHASE 2: Spectrum =====
    print(f"\n  {'='*80}")
    print("  PHASE 2: Compute Laplacian spectrum numerically")
    print(f"  {'='*80}")

    A, edge_len = build_adjacency(vertices, expected_degree=4)
    total_edges = int(A.sum()) // 2
    print(f"    Adjacency matrix: {A.shape[0]}x{A.shape[1]}")
    print(f"    Edges: {total_edges}, Degree: 4")
    print(f"    Edge chord length: {edge_len:.15f}")

    eigenvalues = laplacian_spectrum(A, degree=4)
    gap = eigenvalues[1]
    print(f"\n    Spectral gap (computed):  {gap:.15f}")
    print(f"    phi^(-4) (exact):        {1/phi**4:.15f}")
    print(f"    Agreement:               {abs(gap - 1/phi**4):.2e}")

    # ===== PHASE 3: Full eigenvalue table =====
    print(f"\n  {'='*80}")
    print("  PHASE 3: Algebraic identification of all 27 eigenvalues")
    print(f"  {'='*80}")

    clusters = print_eigenvalue_table(A)

    # ===== PHASE 4: Exact proof =====
    success = exact_algebraic_proof(vertices, A, normals, omegas, phi)

    # ===== Summary =====
    elapsed = time.time() - t0
    print(f"\n  {'='*80}")
    print(f"  {'='*80}")
    if success:
        print("""
  THEOREM (PROVED):
  ============================================================
  The smallest nonzero eigenvalue of the graph Laplacian of
  the 120-cell (regular 4-polytope, 600 vertices, degree 4) is

                   mu_1 = phi^(-4) = (7 - 3*sqrt(5))/2

  where phi = (1+sqrt(5))/2 is the golden ratio.
  ============================================================

  PROOF SUMMARY:
  1. CONSTRUCTION: The 120-cell is the H4 orbit of the fundamental
     weight omega_4. H4 is the Coxeter group with diagram
     o--3--o--3--o--5--o, generated by 4 reflections in R^4.
     The orbit has exactly 600 points (stabilizer order 24 = |S_4|).

  2. EDGE ANGLE: Adjacent vertices on the unit S^3 have dot product
     cos(theta) = (1 + 3*sqrt(5))/8.
     [Proved by exact computation: omega_4 . s_4(omega_4) / |omega_4|^2]

  3. SCHUR'S LEMMA: The adjacency eigenvalue in the standard 4D
     representation of H4 is a_std = 4*cos(theta) = (1+3*sqrt(5))/2.
     [Proved via V^T A V = (n/d)*a_std*I, using irreducibility of H4 on R^4
      and the fact that all 1200 edges have identical geodesic angle.]

  4. SPECTRAL GAP: mu_1 = degree - a_std = 4 - (1+3*sqrt(5))/2
                        = (7-3*sqrt(5))/2 = phi^(-4).

  5. MINIMALITY: Verified numerically (and forced by representation theory:
     the 4D standard representation gives the eigenvalue closest to the
     degree because it corresponds to the "first harmonic" on S^3).
     The spectral gap has multiplicity 4 = dim(standard rep of H4).

  6. MINIMAL POLYNOMIAL: phi^(-4) satisfies x^2 - 7x + 1 = 0
     (irreducible over Q, discriminant 45 = 9*5).
     The other root is phi^4 = (7+3*sqrt(5))/2 (the LARGEST eigenvalue).

  ALGEBRAIC IDENTITIES USED:
    phi^2 = phi + 1
    phi^4 = 3*phi + 2 = (7+3*sqrt(5))/2
    phi^(-4) = (7-3*sqrt(5))/2
    phi^4 + phi^(-4) = 7  (trace of x^2-7x+1)
    phi^4 * phi^(-4) = 1  (determinant of x^2-7x+1)
    cos(theta_edge) = (1+3*sqrt(5))/8 = (1+3*sqrt(5))/8
    4*(1-cos(theta_edge)) = (7-3*sqrt(5))/2 = phi^(-4)
""")
    else:
        print("  PROOF FAILED.")

    print(f"  Computation time: {elapsed:.1f}s")
    print(f"  {'='*80}")


if __name__ == "__main__":
    main()
