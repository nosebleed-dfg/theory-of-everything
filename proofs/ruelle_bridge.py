#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RUELLE_BRIDGE — rho-twisted Ruelle zeta on S^3/2I bridging Selberg trace formula to icosahedral Artin L-function
nos3bl33d

Twisted spectrum via Frobenius, R_rho(s) evaluation, golden/order-3/involution decomposition,
Artin Gamma factor comparison, special value checks.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace',
                              line_buffering=True)

import numpy as np
from mpmath import (mp, mpf, mpc, pi as mpi, sin as mpsin, cos as mpcos,
                    exp as mpexp, log as mplog, gamma as mpgamma, sqrt as mpsqrt,
                    nstr, fac, power as mppow, diff as mpdiff, polylog, zeta as mpzeta,
                    fsum, inf as mpinf, re as mpre, im as mpim, arg as mparg,
                    atan2 as mpatan2, conj as mpconj, fabs, nprint)
import time

mp.dps = 50  # high precision throughout

# ============================================================================
# SECTION 0: Binary Icosahedral Group 2I -- Conjugacy Classes in SU(2)
# ============================================================================

# Each element of 2I in SU(2) is a rotation by angle 2*theta about some axis.
# theta is the "half-angle" parameter. Eigenvalues of the SU(2) matrix are
# exp(+/- i*theta).

# Conjugacy classes of 2I (order 120):
# (name, size, half_angle theta)
CONJ_CLASSES = [
    ("I",   1,  mpf(0)),          # identity
    ("-I",  1,  mpi),             # central element
    ("C10a", 12, mpi / 5),        # order 10, angle pi/5
    ("C10b", 12, 2 * mpi / 5),    # order 10, angle 2pi/5
    ("C10c", 12, 3 * mpi / 5),    # order 10, angle 3pi/5
    ("C10d", 12, 4 * mpi / 5),    # order 10, angle 4pi/5
    ("C6a",  20, mpi / 3),        # order 6, angle pi/3
    ("C6b",  20, 2 * mpi / 3),    # order 6, angle 2pi/3
    ("C4",   30, mpi / 2),        # order 4, angle pi/2
]

assert sum(c[1] for c in CONJ_CLASSES) == 120, "Group order must be 120"

# The golden ratio
PHI = (1 + mpsqrt(5)) / 2  # (1+sqrt(5))/2


# ============================================================================
# Character of (n+1)-dim SU(2) irrep
# ============================================================================

def chi_n(n, theta):
    """Character of the (n+1)-dim irrep of SU(2) at element with half-angle theta.
    chi_n(g) = sin((n+1)*theta) / sin(theta).
    Special cases: theta=0 -> n+1, theta=pi -> (-1)^n * (n+1).
    """
    n1 = n + 1
    if theta == 0:
        return mpf(n1)
    if fabs(theta - mpi) < mpf('1e-30'):
        return mpf((-1)**n * n1)
    return mpsin(n1 * theta) / mpsin(theta)


# ============================================================================
# SECTION 1: Character of the 2D representation rho
# ============================================================================

# rho: 2I -> SL(2,C) is the standard 2D representation.
# Since elements of 2I in SU(2) have eigenvalues exp(+/- i*theta),
# chi_rho(g) = 2*cos(theta).
# For the identity representation (n=1 irrep of SU(2)): chi_1(g) = 2*cos(theta).
# This IS the 2D rep: the defining representation of SU(2) restricted to 2I.

def chi_rho(theta):
    """Character of the 2D representation rho at element with half-angle theta."""
    return 2 * mpcos(theta)


def det_rho(theta):
    """Determinant of rho(g) for g with half-angle theta.
    Since rho: 2I -> SU(2) subset SL(2,C), det = 1 always.
    But for -I: rho(-I) = -I_2, so det(-I_2) = 1 (it's 2x2).
    """
    return mpf(1)


# Verify character values match the specification
print("=" * 80)
print("  RUELLE ZETA AND rho-TWISTED SPECTRUM FOR S^3/2I")
print("  Poincare Dodecahedral Space")
print("=" * 80)

print("\n--- Character values of rho (2D rep) ---")
print(f"{'Class':>8} {'Size':>5} {'theta':>12} {'chi_rho':>16} {'Expected':>16}")
print("-" * 60)

expected = {
    "I": mpf(2), "-I": mpf(-2),
    "C10a": PHI, "C10b": 1/PHI, "C10c": -1/PHI, "C10d": -PHI,
    "C6a": mpf(1), "C6b": mpf(-1), "C4": mpf(0),
}

for name, size, theta in CONJ_CLASSES:
    cv = chi_rho(theta)
    ev = expected[name]
    ok = "OK" if fabs(cv - ev) < mpf('1e-40') else "FAIL"
    print(f"{name:>8} {size:>5} {nstr(theta, 8):>12} {nstr(cv, 12):>16} "
          f"{nstr(ev, 12):>16}  [{ok}]")


# ============================================================================
# PART A: rho-TWISTED SPECTRUM
# ============================================================================

print("\n" + "=" * 80)
print("  PART A: rho-TWISTED SPECTRUM")
print("  mult_rho(n) = (2/120) * sum_g chi_rho(g)* x chi_n(g) x (n+1)")
print("=" * 80)

# The factor of 2 comes from: mult = dim(rho) * (1/|G|) * sum chi_rho* chi_n * dim(V_n)
# Wait, let's be precise. For the rho-twisted Laplacian on the associated vector bundle,
# the multiplicity of eigenvalue n(n+2) is:
#   mult_rho(n) = (1/|G|) * sum_{g in G} chi_rho(g)* * chi_n(g) * (n+1)
# where chi_n(g) = sin((n+1)*theta)/sin(theta) is the character of V_n,
# and (n+1) = dim(V_n). Wait no -- the Frobenius formula gives:
#   mult_rho(n) = dim(Hom_G(V_n, rho)) = (1/|G|) * sum_{g} chi_n(g)* chi_rho(g)
# But on S^3, eigenvalue n(n+2) has eigenspace V_n tensor V_n (left x right).
# The quotient M = S^3/G uses G acting on left, so:
#   mult_rho(n) = dim((V_n tensor rho)^G) * dim(V_n)
# Hmm. Actually the correct formula for the twisted spectrum is:
#   mult_rho(n) = (1/|G|) sum_{g in G} overline{chi_rho(g)} * |chi_n(g)|^2
# No wait. Let me think carefully.
#
# On S^3 = SU(2), the eigenspace for eigenvalue n(n+2) is the space of
# matrix coefficients of the (n+1)-dim irrep, which as a G x G representation
# is V_n tensor V_n* (left and right regular representation).
# The quotient S^3/G_left means we take G_left-invariants. On the rho-twisted
# bundle, we have V_rho tensor L^2(S^3) and we need:
#   mult_rho(n) = dim(Hom_G(V_rho*, V_n)) * dim(V_n)
# where the first factor is how many times rho appears in V_n, and the second
# is the right-action multiplicity.
#
# By Frobenius:
#   dim(Hom_G(rho, V_n)) = (1/|G|) sum_{g in G} chi_rho(g)* * chi_n(g)
#
# Then mult_rho(n) = dim(Hom_G(rho, V_n)) * (n+1) = (1/|G|) sum_g chi_rho(g)* chi_n(g) * (n+1)
#
# Note: chi_rho is real-valued here (all characters of 2I over reals), so chi_rho* = chi_rho.

N_MAX = 200

def mult_rho(n):
    """rho-twisted multiplicity at eigenvalue n(n+2)."""
    n1 = n + 1
    total = mpf(0)
    for name, size, theta in CONJ_CLASSES:
        total += size * chi_rho(theta) * chi_n(n, theta)
    # Frobenius inner product, then multiply by dim(V_n) = n+1
    inner = total / 120
    result = inner * n1
    return int(round(float(result)))


def mult_untwisted(n):
    """Untwisted multiplicity (trivial rep) on M = S^3/2I."""
    n1 = n + 1
    total = mpf(0)
    for name, size, theta in CONJ_CLASSES:
        total += size * chi_n(n, theta) * n1
    return int(round(float(total / 120)))


# Compute twisted spectrum
twisted_spectrum = []
for n in range(N_MAX + 1):
    lam = n * (n + 2)
    m = mult_rho(n)
    twisted_spectrum.append((n, lam, m))

# Sanity checks
print("\nSanity checks:")
total_twisted_dim = sum(m for _, _, m in twisted_spectrum)
total_untwisted_dim = sum(mult_untwisted(n) for n in range(N_MAX + 1))

# For the trivial rep, sum of mult for n=0..N should approximate |G|^{-1} * sum (n+1)^2
# For rho of dim 2, we expect approximately 2x the untwisted count
print(f"  Total rho-twisted states (n=0..{N_MAX}): {total_twisted_dim}")
print(f"  Total untwisted states  (n=0..{N_MAX}): {total_untwisted_dim}")
print(f"  Ratio twisted/untwisted: {total_twisted_dim/total_untwisted_dim:.4f}")
print(f"  (Expected ratio ~ dim(rho) = 2 if rho is faithful)")

# Print first 30 nonzero twisted multiplicities
print("\n--- First 30 nonzero rho-twisted multiplicities ---")
print(f"{'n':>5} {'lambda_n':>10} {'mult_rho(n)':>12} {'mult_M(n)':>10}")
print("-" * 42)

nonzero_count = 0
nonzero_list = []
for n, lam, m in twisted_spectrum:
    if m != 0:
        mu = mult_untwisted(n)
        print(f"{n:>5} {lam:>10} {m:>12} {mu:>10}")
        nonzero_list.append((n, lam, m))
        nonzero_count += 1
        if nonzero_count >= 30:
            break

# Full listing for analysis
all_nonzero = [(n, lam, m) for n, lam, m in twisted_spectrum if m != 0]
print(f"\n  Total nonzero twisted mults in [0,{N_MAX}]: {len(all_nonzero)}")

# Pattern analysis
print("\n--- Pattern analysis of nonzero n values ---")
nonzero_ns = [n for n, _, m in twisted_spectrum if m != 0]
if len(nonzero_ns) > 1:
    diffs = [nonzero_ns[i+1] - nonzero_ns[i] for i in range(min(30, len(nonzero_ns)-1))]
    print(f"  First 30 nonzero n values: {nonzero_ns[:30]}")
    print(f"  Consecutive differences:   {diffs}")

# Check: n mod 60 pattern (since |2I| = 120 and SU(2) double cover)
if nonzero_ns:
    mods = [n % 60 for n in nonzero_ns[:30]]
    print(f"  n mod 60:                  {mods}")


# ============================================================================
# PART B: RUELLE ZETA R_rho(s) TWISTED BY rho
# ============================================================================

print("\n" + "=" * 80)
print("  PART B: RUELLE ZETA R_rho(s)")
print("  R_rho(s) = prod_{[gamma] prim} det(I - rho(gamma) e^{-s l(gamma)})")
print("=" * 80)

# Primitive closed geodesics on S^3/2I:
# Each conjugacy class C of 2I (excluding identity) gives closed geodesics.
# The geodesic length is l = 2*theta_C (the rotation angle, i.e., 2*half_angle).
#
# For primitiveness: an element g of order d generates a cyclic subgroup <g>.
# The primitive geodesic in this free homotopy class has length l = 2*theta_min
# where theta_min is the smallest half-angle in <g>.
# Higher powers g^k give longer geodesics of length 2*k*theta_min.
#
# In 2I, the cyclic subgroups and their primitive generators:
# - Order 10 subgroups: primitive angle pi/5, powers give 2pi/5, 3pi/5, 4pi/5, pi
# - Order 6 subgroups:  primitive angle pi/3, powers give 2pi/3, pi
# - Order 4 subgroups:  primitive angle pi/2, powers give pi
#
# But for the Ruelle product, we need to identify which conjugacy classes
# are primitive (not proper powers of other group elements).
#
# A conjugacy class C is PRIMITIVE if its elements are not proper powers
# of elements in another class.
#
# For 2I:
# - C10a (theta=pi/5): primitive. These are generators of order-10 cyclic subgroups.
# - C10b (theta=2pi/5): g^2 where g in C10a. NOT primitive.
# - C10c (theta=3pi/5): g^3 where g in C10a. Could be primitive in its own right...
#   Actually g^3 generates a different cyclic subgroup than g (since gcd(3,10)=1,
#   g^3 is also a generator of <g>). So C10c elements generate the SAME cyclic
#   subgroups as C10a. The primitive geodesic for this subgroup is still C10a.
# - C10d (theta=4pi/5): g^4. gcd(4,10)=2, so g^4 generates a subgroup of order 5.
#   But g^4 = (g^2)^2, so C10d is the square of C10b.
#   Actually wait: g^2 is in C10b (angle 2pi/5). (g^2)^2 = g^4 is in C10d (angle 4pi/5).
#   So C10d is the square of C10b. But C10b is itself g^2, so C10d = g^4.
#   The element g^4 has order 10/gcd(4,10) = 5. It generates a cyclic subgroup of
#   order 5. Its primitive representative is g^4 itself (as an order-5 element).
#   Hmm but g^2 has order 5, and (g^2)^2 = g^4. So g^4 is the SQUARE of the
#   order-5 primitive g^2.
#
# Let me think about this differently using the actual Ruelle zeta formulation.
#
# The Ruelle zeta on a compact hyperbolic manifold (or spherical space form) is:
# R_rho(s) = prod over primitive closed geodesics gamma_p of
#            det(I - rho(h_gamma) * e^{-s * l(gamma_p)})
#
# where h_gamma is the holonomy (the group element corresponding to gamma_p).
#
# For S^3/Gamma, the primitive closed geodesics correspond to PRIMITIVE conjugacy
# classes: those classes [g] where g is not a proper power of any other element.
#
# In 2I, the primitive classes are determined by which elements generate
# maximal cyclic subgroups:
#
# Order structure of 2I:
#   - 1 element of order 1 (identity)
#   - 1 element of order 2 (-I)
#   - 30 elements of order 4 (the C4 class)
#   - 20 elements of order 6 (split into C6a and C6b)
#   - 24 elements of order 10 (split into C10a, C10b, C10c, C10d)
#   Wait, that's 1+1+30+20+24 = 76. Missing 44. Let me recount.
#
# Actually the orders in 2I:
#   Element in C10a (theta=pi/5): eigenvalues exp(+/- i*pi/5). Order = 2pi/(2*pi/5)... no.
#   In SU(2), element with half-angle theta has order = 2pi/(2*theta) if 2*theta divides 2*pi.
#   Wait: the element is exp(i*theta*sigma_z) = diag(exp(i*theta), exp(-i*theta)).
#   Its order is the smallest k such that k*theta is a multiple of pi (not 2*pi,
#   because exp(i*k*theta) = 1 requires k*theta = 2*pi*m, but in SU(2), -I is a
#   separate element, so we need k*theta = pi*m for some integer m, and the element
#   equals +/-I when m is even/odd).
#
#   For theta=pi/5: smallest k with k*(pi/5) = pi*m -> k=5, m=1 -> g^5 = -I.
#   So g^10 = I. Order = 10.
#
#   For theta=2pi/5: k*(2pi/5) = pi*m -> k=5, m=2 -> g^5 = (-1)^2 I = I.
#   Wait that's wrong. Let me compute directly.
#   g = diag(exp(2pi*i/5), exp(-2pi*i/5)). g^5 = diag(exp(2pi*i), exp(-2pi*i)) = I.
#   So order = 5. And g is NOT a primitive root for order 10.
#
# OK let me just systematically determine orders and primitiveness.

def element_order(theta):
    """Order of SU(2) element with half-angle theta in the group 2I.
    Element is diag(exp(i*theta), exp(-i*theta)).
    g^k = diag(exp(i*k*theta), exp(-i*k*theta)).
    g^k = I iff k*theta = 2*pi*m for some integer m.
    g^k = -I iff k*theta = pi*(2m+1).
    Order = smallest positive k with g^k = +/- I... no, order = smallest k with g^k = I.
    """
    if theta == 0:
        return 1
    # Find smallest k>0 such that k*theta is a multiple of 2*pi
    # theta = p*pi/q in lowest terms
    # We need k*p/q to be even integer, i.e., k*p = 2*q*m
    # Smallest k = 2*q / gcd(p, 2*q)
    # But we should just compute numerically for our specific angles.
    for k in range(1, 121):
        val = k * theta / mpi
        if fabs(val - round(float(val))) < mpf('1e-30'):
            rounded = int(round(float(val)))
            if rounded % 2 == 0:  # k*theta = even*pi means g^k = I
                return k
    return None


print("\n--- Conjugacy class orders and primitiveness ---")
print(f"{'Class':>8} {'Size':>5} {'theta/pi':>10} {'Order':>6} {'Primitive?':>10}")
print("-" * 45)

primitive_classes = []
non_identity_classes = []

for name, size, theta in CONJ_CLASSES:
    if name in ("I", "-I"):
        ord_val = 1 if name == "I" else 2
    else:
        ord_val = element_order(theta)

    # An element is primitive if it is NOT a proper power of another element in 2I.
    # More precisely, for the Ruelle zeta, a conjugacy class is primitive if
    # the corresponding closed geodesic is not a multiple traversal of a shorter one.
    #
    # The geodesic for class C has length l_C = 2*theta_C.
    # It's primitive if there's no class C' with theta_C' = theta_C / k for some k >= 2.
    is_prim = True
    if name not in ("I", "-I"):
        for k in range(2, 20):
            reduced_theta = theta / k
            # Check if reduced_theta matches any other class
            for name2, size2, theta2 in CONJ_CLASSES:
                if name2 not in ("I", "-I") and fabs(reduced_theta - theta2) < mpf('1e-30'):
                    is_prim = False
                    break
            if not is_prim:
                break

    prim_str = "YES" if is_prim else "no"
    theta_over_pi = nstr(theta / mpi, 6) if theta != 0 else "0"
    print(f"{name:>8} {size:>5} {theta_over_pi:>10} {ord_val:>6} {prim_str:>10}")

    if name not in ("I",):
        non_identity_classes.append((name, size, theta, ord_val, is_prim))
    if is_prim and name not in ("I", "-I"):
        primitive_classes.append((name, size, theta))


print(f"\nPrimitive classes: {[c[0] for c in primitive_classes]}")

# For the Ruelle zeta, we actually want to use ALL non-identity classes,
# but account for which are primitive and which are powers.
# The standard Ruelle zeta is:
#   R_rho(s) = prod_{gamma primitive} det(I - rho(h_gamma) e^{-s l_gamma})
#
# Since we're on a SPHERICAL space form (positive curvature), not hyperbolic,
# the product is FINITE -- it's just over the finitely many conjugacy classes.
#
# But let's also compute the full product over ALL non-identity non-central classes,
# treating each as contributing one factor per conjugacy class member.
# Actually, each conjugacy class C gives |C| closed geodesics of the same length,
# all with conjugate holonomies. Since det is a class function, they all contribute
# the same factor.
#
# More precisely: the Euler product is
#   R_rho(s) = prod_{C primitive} [det(I - rho(g_C) e^{-s l_C})]^{n_C}
# where n_C is the number of distinct primitive geodesics in class C.
# For S^3/Gamma acting freely, n_C = |C| (each element in C gives a distinct
# primitive geodesic).

# Actually for the Ruelle product on quotient spaces, each conjugacy class [g]
# (modulo the centralizer) contributes. But let me use the LOGARITHMIC form
# which is cleaner and avoids these subtleties.
#
# The key formula connecting Ruelle to spectral data (Selberg trace formula):
#   log R_rho(s) = - sum_{C != I} sum_{k=1}^{K(C)} (contribution from C^k)
# But on a spherical space, the sum over k is FINITE.
#
# Let's just directly compute both:
# (1) The "naive" product over all non-identity classes
# (2) The product over primitive classes only

def ruelle_local_factor(theta_C, s, chi_val=None):
    """Local factor det(I - rho(g_C) e^{-s l_C}) for one geodesic.
    l_C = 2*theta_C (geodesic length on unit S^3).
    det(I - rho(g) u) = 1 - chi_rho(g) u + det(rho(g)) u^2
                      = 1 - chi_rho(g) u + u^2  (since det(rho(g)) = 1 for SL2)
    """
    u = mpexp(-s * 2 * theta_C)
    if chi_val is None:
        chi_val = chi_rho(theta_C)
    return 1 - chi_val * u + u**2


def log_ruelle_naive(s):
    """log R_rho(s) using ALL non-identity conjugacy classes.
    Each class C contributes |C| * log det(I - rho(g_C) e^{-s*2*theta_C}).
    """
    result = mpf(0)
    for name, size, theta, order, is_prim in non_identity_classes:
        if name == "-I":
            # For -I: theta=pi, l=2pi, rho(-I) = -I_2
            # det(I - (-I_2) u) = det(I + u I_2) = (1+u)^2
            u = mpexp(-s * 2 * mpi)
            factor = (1 + u)**2
        else:
            factor = ruelle_local_factor(theta, s)
        if fabs(factor) > mpf('1e-100'):
            result += size * mplog(factor)
        else:
            result += size * mpf('-1e50')  # handle zeros carefully
    return result


def log_ruelle_primitive(s):
    """log R_rho(s) using only PRIMITIVE conjugacy classes.
    Each primitive class C contributes |C| * log det(I - rho(g_C) e^{-s*2*theta_C}).
    """
    result = mpf(0)
    for name, size, theta in primitive_classes:
        factor = ruelle_local_factor(theta, s)
        if fabs(factor) > mpf('1e-100'):
            result += size * mplog(factor)
        else:
            result += size * mpf('-1e50')
    return result


def ruelle_logderiv(s, use_primitive=False):
    """Compute -R_rho'/R_rho(s) = sum_C n_C sum_{k>=1} tr(rho(g_C)^k) l_C e^{-k s l_C}.
    For spherical space, we don't need the k-sum -- each class contributes directly.

    Actually, the log derivative of det(I - A u) w.r.t. u is:
    d/du log det(I - Au) = -tr(A (I - Au)^{-1})
    And d/ds [u = e^{-s l}] = -l e^{-sl}
    So d/ds log det(I - A e^{-sl}) = tr(A (I - A e^{-sl})^{-1}) l e^{-sl}

    For 2x2 with det(A)=1: (I - Au)^{-1} = (1 - chi u + u^2)^{-1} (I - Au)^adj
    tr(A (I-Au)^{-1}) = (chi - 2u) / (1 - chi u + u^2)
    """
    result = mpf(0)
    classes = primitive_classes if use_primitive else [(n, s2, t) for n, s2, t, o, p in non_identity_classes if n != "-I"]

    for name, size, theta in classes:
        l_C = 2 * theta
        u = mpexp(-s * l_C)
        chi_val = chi_rho(theta)
        denom = 1 - chi_val * u + u**2
        if fabs(denom) > mpf('1e-100'):
            numer = (chi_val - 2 * u) * l_C * u
            result += size * numer / denom

    # Handle -I if not primitive only
    if not use_primitive:
        l_neg = 2 * mpi
        u = mpexp(-s * l_neg)
        # rho(-I) = -I_2, so det(I + u I_2) = (1+u)^2
        # d/ds log(1+u)^2 = 2 * (-l) * u / (1+u) = -2*l*u/(1+u)
        # So -R'/R contribution = 2*l*u/(1+u)
        result += 1 * 2 * l_neg * u / (1 + u)

    return result


# Evaluate R_rho(s) on the real line
print("\n--- Ruelle zeta R_rho(s) on real axis (s = sigma, t = 0) ---")
print(f"{'sigma':>8} {'log|R_naive|':>18} {'log|R_prim|':>18} {'arg R_naive':>14}")
print("-" * 62)

sigma_values = np.arange(0.1, 3.01, 0.1)
ruelle_real_data = []

for sigma in sigma_values:
    s = mpf(sigma)
    lr_naive = log_ruelle_naive(s)
    lr_prim = log_ruelle_primitive(s)
    ruelle_real_data.append((sigma, float(mpre(lr_naive)), float(mpre(lr_prim))))
    print(f"{sigma:>8.2f} {float(mpre(lr_naive)):>18.10f} {float(mpre(lr_prim)):>18.10f} "
          f"{float(mpim(lr_naive)):>14.8f}")


# Find real zeros (where log|R| changes sign or is very small)
print("\n--- Searching for real zeros of R_rho(s) ---")
fine_sigmas = np.arange(0.01, 3.001, 0.01)
log_r_real = []
for sigma in fine_sigmas:
    s = mpf(sigma)
    lr = log_ruelle_naive(s)
    log_r_real.append((sigma, float(mpre(lr))))

# Find sign changes or near-zeros in exp(Re(log R))
# R(s) = 0 when log R -> -inf, which means a factor = 0
print("  Checking for zeros (local minima of Re(log R))...")
min_val = float('inf')
min_sigma = 0
for i in range(1, len(log_r_real) - 1):
    s0, v0 = log_r_real[i-1]
    s1, v1 = log_r_real[i]
    s2, v2 = log_r_real[i+1]
    if v1 < v0 and v1 < v2 and v1 < min_val:
        min_val = v1
        min_sigma = s1
    # Also check if any factor vanishes
    # Factor vanishes when 1 - chi u + u^2 = 0, i.e., u = (chi +/- sqrt(chi^2-4))/2

print(f"  Deepest minimum: sigma = {min_sigma:.4f}, log|R| = {min_val:.10f}")

# Analyze zeros of individual factors
print("\n--- Zeros of individual local factors ---")
print("  Factor det(I - rho(g) e^{-s*2*theta}) = 1 - chi*u + u^2 = 0")
print("  Solutions: u = (chi +/- sqrt(chi^2 - 4)) / 2")
print("  Since u = e^{-2*sigma*theta} > 0, need real positive u.")
print(f"{'Class':>8} {'chi':>10} {'disc':>10} {'u solutions':>30} {'sigma zeros':>20}")
print("-" * 82)

for name, size, theta in primitive_classes:
    chi_val = float(chi_rho(theta))
    disc = chi_val**2 - 4
    theta_f = float(theta)
    l_f = 2 * theta_f
    if disc >= 0:
        u1 = (chi_val + np.sqrt(disc)) / 2
        u2 = (chi_val - np.sqrt(disc)) / 2
        u_str = f"({u1:.6f}, {u2:.6f})"
        s_str = ""
        for u in [u1, u2]:
            if u > 0:
                sig = -np.log(u) / l_f
                s_str += f"{sig:.6f} "
    else:
        u_str = "complex"
        s_str = "none (complex)"
    print(f"{name:>8} {chi_val:>10.6f} {disc:>10.6f} {u_str:>30} {s_str:>20}")


# Critical line zeros
print("\n--- Ruelle zeta on critical line s = 1/2 + it ---")
print(f"{'t':>8} {'Re(log R)':>18} {'Im(log R)':>18} {'|R|':>14}")
print("-" * 62)

t_values = np.arange(0, 50.1, 1.0)
crit_line_data = []

for t in t_values:
    s = mpc(mpf('0.5'), mpf(t))
    lr = log_ruelle_naive(s)
    mod_r = mpexp(mpre(lr))
    crit_line_data.append((t, float(mpre(lr)), float(mpim(lr)), float(mod_r)))
    if t <= 20 or t % 5 == 0:
        print(f"{t:>8.1f} {float(mpre(lr)):>18.10f} {float(mpim(lr)):>18.10f} "
              f"{float(mod_r):>14.10f}")

# Find approximate zeros on critical line (|R| small)
print("\n  Approximate zeros on critical line (|R| < 0.1):")
for t, rlr, ilr, mr in crit_line_data:
    if mr < 0.1:
        print(f"    t = {t:.1f}, |R| = {mr:.8f}")


# ============================================================================
# PART C: THREE-CLASS DECOMPOSITION
# ============================================================================

print("\n" + "=" * 80)
print("  PART C: THREE-CLASS DECOMPOSITION OF -log R_rho")
print("  Golden / Order-3 / Involution")
print("=" * 80)

# Evaluate at s = 1
s_eval = mpf(1)

# Golden contribution: classes C10a, C10b, C10c, C10d
golden_contrib = mpf(0)
golden_classes = [("C10a", 12, mpi/5), ("C10b", 12, 2*mpi/5),
                  ("C10c", 12, 3*mpi/5), ("C10d", 12, 4*mpi/5)]
for name, size, theta in golden_classes:
    factor = ruelle_local_factor(theta, s_eval)
    golden_contrib += size * (-mplog(factor))

# Order-3 contribution: classes C6a, C6b
order3_contrib = mpf(0)
order3_classes = [("C6a", 20, mpi/3), ("C6b", 20, 2*mpi/3)]
for name, size, theta in order3_classes:
    factor = ruelle_local_factor(theta, s_eval)
    order3_contrib += size * (-mplog(factor))

# Involution contribution: class C4
invol_contrib = mpf(0)
invol_classes = [("C4", 30, mpi/2)]
for name, size, theta in invol_classes:
    factor = ruelle_local_factor(theta, s_eval)
    invol_contrib += size * (-mplog(factor))

# Central element -I
central_contrib = mpf(0)
u_neg = mpexp(-s_eval * 2 * mpi)
central_contrib = -mplog((1 + u_neg)**2)

total_log_r = -(golden_contrib + order3_contrib + invol_contrib + central_contrib)

print(f"\nAt s = 1:")
print(f"  -log R_rho(1) decomposition:")
print(f"    Golden contribution  (48 elements, theta=k*pi/5): {nstr(golden_contrib, 20)}")
print(f"    Order-3 contribution (40 elements, theta=pi/3,2pi/3): {nstr(order3_contrib, 20)}")
print(f"    Involution contrib   (30 elements, theta=pi/2): {nstr(invol_contrib, 20)}")
print(f"    Central (-I) contrib (1 element, theta=pi):  {nstr(central_contrib, 20)}")
print(f"    Total -log R_rho(1) = {nstr(golden_contrib + order3_contrib + invol_contrib + central_contrib, 20)}")
print(f"    R_rho(1) = {nstr(mpexp(total_log_r), 20)}")

# Weights = contribution / class size
print(f"\n  Per-element weights:")
golden_weight = golden_contrib / 48
order3_weight = order3_contrib / 40
invol_weight = invol_contrib / 30
print(f"    Golden (per element):     {nstr(golden_weight, 15)}")
print(f"    Order-3 (per element):    {nstr(order3_weight, 15)}")
print(f"    Involution (per element): {nstr(invol_weight, 15)}")

# Compare with the target numbers
print(f"\n  Comparison with A5 structure constants:")
print(f"    Golden class total 48    vs  2d = 2*dim(icosahedral) = 2*3 = 6?  (48/6 = {48/6})")
print(f"    Order-3 class total 40   vs  V = 20 (vertices)?                  (40/20 = {40/20})")
print(f"    Involution class total 30 vs dp = 15 (edge midpoints)?            (30/15 = {30/15})")
print(f"    Note: 2I has DOUBLE the elements of A5 (120 vs 60).")
print(f"    A5 class sizes: 12+12=24 (order 5), 20 (order 3), 15 (order 2), 1 (identity)")
print(f"    2I sizes: 12*4=48 (golden), 20*2=40 (order 3/6), 30 (order 4), 1+1=2 (central)")

# Ratios between contributions
print(f"\n  Weight ratios:")
print(f"    Golden / Order-3:     {nstr(golden_contrib / order3_contrib, 15)}")
print(f"    Golden / Involution:  {nstr(golden_contrib / invol_contrib, 15)}")
print(f"    Order-3 / Involution: {nstr(order3_contrib / invol_contrib, 15)}")

# At multiple s values
print(f"\n  Class contributions vs s:")
print(f"{'s':>6} {'Golden':>18} {'Order-3':>18} {'Involution':>18} {'Central':>18}")
print("-" * 82)

for s_val in [0.5, 1.0, 1.5, 2.0, 3.0]:
    s = mpf(s_val)
    gc = mpf(0)
    for _, size, theta in golden_classes:
        gc += size * (-mplog(ruelle_local_factor(theta, s)))
    oc = mpf(0)
    for _, size, theta in order3_classes:
        oc += size * (-mplog(ruelle_local_factor(theta, s)))
    ic = mpf(0)
    for _, size, theta in invol_classes:
        ic += size * (-mplog(ruelle_local_factor(theta, s)))
    cc = -mplog((1 + mpexp(-s * 2 * mpi))**2)

    print(f"{s_val:>6.1f} {nstr(gc, 14):>18} {nstr(oc, 14):>18} "
          f"{nstr(ic, 14):>18} {nstr(cc, 14):>18}")


# ============================================================================
# PART D: COMPARISON WITH ARTIN L-FUNCTION
# ============================================================================

print("\n" + "=" * 80)
print("  PART D: BRIDGE TO ARTIN L-FUNCTION")
print("  Comparing R_rho(s) structure with L(s, rho)")
print("=" * 80)

# The Artin L-function for the 2D icosahedral representation rho is:
#   L(s, rho) = prod_p det(I - rho(Frob_p) p^{-s})^{-1}
# This is an INVERSE product (Euler product), while Ruelle is a DIRECT product.
#
# The key bridge: on S^3/Gamma, the Selberg zeta Z(s) satisfies
#   Z(s) = prod_{n=0}^infty R(s+n)^{(-1)^n ...}
# The relationship involves shifted products.
#
# For our spherical space form, the spectral zeta function
#   zeta_M(s) = sum_{n: mult>0} mult(n) / (n(n+2))^s
# is related to the Ruelle/Selberg zeta via a Mellin transform.
#
# The Gamma factor for a 2D Galois representation is:
#   Gamma_rho(s) = pi^{-s} Gamma(s/2) Gamma((s+1)/2) = pi^{-s} Gamma(s) / 2^{s-1}
# (using the duplication formula).
# More precisely, for the icosahedral rep (odd, since det = sgn),
# with local factors at infinity determined by rho(complex conjugation):
#   Gamma_rho(s) = Gamma_R(s) * Gamma_R(s+1) = pi^{-s} Gamma(s/2) Gamma((s+1)/2)

print("\n--- Gamma factor analysis ---")
print("  For rho: 2D, odd (rho(-1) has eigenvalues +1, -1)")
print("  Gamma_rho(s) = pi^{-s} * Gamma(s/2) * Gamma((s+1)/2)")

for s_val in [0.5, 1.0, 1.5, 2.0]:
    s = mpf(s_val)
    gamma_factor = mpi**(-s) * mpgamma(s/2) * mpgamma((s+1)/2)
    print(f"  Gamma_rho({s_val}) = {nstr(gamma_factor, 15)}")

# Compare log derivatives
print("\n--- Log derivative comparison ---")
print("  d/ds log R_rho(s) vs d/ds log Gamma_rho(s)")
header_r = "(R_rho'/R_rho)(s)"
header_g = "(Gamma'/Gamma)(s)"
print(f"{'s':>6} {header_r:>22} {header_g:>22} {'Difference':>22}")
print("-" * 76)

def gamma_rho_logderiv(s):
    """d/ds log Gamma_rho(s) = -log(pi) + psi(s/2)/2 + psi((s+1)/2)/2."""
    from mpmath import digamma
    return -mplog(mpi) + digamma(s/2)/2 + digamma((s+1)/2)/2

for s_val in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    s = mpf(s_val)
    # Numerical log derivative of R_rho
    ds = mpf('1e-8')
    r_logderiv = (log_ruelle_naive(s + ds) - log_ruelle_naive(s - ds)) / (2 * ds)
    g_logderiv = gamma_rho_logderiv(s)

    print(f"{s_val:>6.1f} {nstr(mpre(r_logderiv), 15):>22} {nstr(g_logderiv, 15):>22} "
          f"{nstr(mpre(r_logderiv) - g_logderiv, 15):>22}")


# Spectral interpretation: compute the spectral zeta function from twisted spectrum
print("\n--- Spectral zeta from rho-twisted spectrum ---")
print("  Z_rho(s) = sum_{n: mult_rho(n)>0} mult_rho(n) / (n(n+2))^s")

for s_val in [1, 2, 3]:
    zeta_rho = mpf(0)
    for n, lam, m in twisted_spectrum:
        if lam > 0 and m > 0:
            zeta_rho += mpf(m) / mpf(lam)**s_val
    print(f"  Z_rho({s_val}) = {nstr(zeta_rho, 25)}")

# Compare spectral zeta with the untwisted one
print("\n  Comparison with untwisted spectral zeta Z_M(s):")
for s_val in [1, 2]:
    zeta_rho = mpf(0)
    zeta_M = mpf(0)
    for n in range(1, N_MAX + 1):
        lam = n * (n + 2)
        mr = mult_rho(n)
        mu = mult_untwisted(n)
        if mr > 0:
            zeta_rho += mpf(mr) / mpf(lam)**s_val
        if mu > 0:
            zeta_M += mpf(mu) / mpf(lam)**s_val
    ratio = zeta_rho / zeta_M if zeta_M != 0 else mpf(0)
    print(f"  Z_rho({s_val})/Z_M({s_val}) = {nstr(ratio, 15)}")

# The bridge equation: R_rho should encode the same information as L(s, rho)
# up to Gamma factors. Specifically:
#   Lambda(s, rho) = N^{s/2} Gamma_rho(s) L(s, rho)
# satisfies a functional equation Lambda(s) = epsilon * Lambda(1-s).
#
# For the icosahedral rep, the conductor N = 5^5 * ... depends on ramification.
# The KNOWN icosahedral Artin L-function (Buhler's computation) has conductor 800.
# Let's see if our Ruelle data is consistent.

print("\n--- Bridge equation test ---")
print("  If R_rho(s) ~ Gamma_rho(s) * L(s, rho), then")
print("  R_rho(s) / Gamma_rho(s) should be meromorphic and look like L(s, rho).")

for s_val in [1.0, 1.5, 2.0]:
    s = mpf(s_val)
    lr = log_ruelle_naive(s)
    R_val = mpexp(lr)
    gamma_val = mpi**(-s) * mpgamma(s/2) * mpgamma((s+1)/2)
    ratio = R_val / gamma_val
    print(f"  s = {s_val}: R_rho = {nstr(R_val, 12)}, Gamma = {nstr(gamma_val, 12)}, "
          f"R/Gamma = {nstr(ratio, 12)}")


# ============================================================================
# PART E: CRITICAL CHECKS -- EFFECTIVE ALPHA
# ============================================================================

print("\n" + "=" * 80)
print("  PART E: CRITICAL CHECKS -- EFFECTIVE ALPHA")
print("  Searching for alpha = 1/137.035999... in Ruelle data")
print("=" * 80)

alpha_target = mpf('0.0072973525693')  # 1/137.035999...
alpha_inv_target = mpf('137.035999084')

print(f"\n  Target: alpha = {nstr(alpha_target, 15)}")
print(f"  Target: 1/alpha = {nstr(alpha_inv_target, 15)}")

# Test various combinations
s1 = mpf(1)
lr1 = log_ruelle_naive(s1)
R1 = mpexp(lr1)

s_half = mpf('0.5')
lr_half = log_ruelle_naive(s_half)
R_half = mpexp(lr_half)

# Log derivative at s=1
ds = mpf('1e-8')
lr1p = (log_ruelle_naive(s1 + ds) - log_ruelle_naive(s1 - ds)) / (2 * ds)

print(f"\n  R_rho(1)           = {nstr(R1, 20)}")
print(f"  log|R_rho(1)|      = {nstr(mpre(lr1), 20)}")
print(f"  R_rho(1/2)         = {nstr(R_half, 20)}")
print(f"  log|R_rho(1/2)|    = {nstr(mpre(lr_half), 20)}")
print(f"  -R_rho'/R_rho(1)   = {nstr(-mpre(lr1p), 20)}")

# Try various combinations
candidates = {
    "R_rho(1)": R1,
    "1/R_rho(1)": 1/R1,
    "log|R_rho(1)|": mpre(lr1),
    "-log|R_rho(1)|": -mpre(lr1),
    "R_rho(1/2)": R_half,
    "log|R_rho(1/2)|": mpre(lr_half),
    "-R'/R(1)": -mpre(lr1p),
    "R'/R(1)": mpre(lr1p),
    "golden/total": golden_contrib / (golden_contrib + order3_contrib + invol_contrib),
    "invol/total": invol_contrib / (golden_contrib + order3_contrib + invol_contrib),
    "order3/total": order3_contrib / (golden_contrib + order3_contrib + invol_contrib),
    "golden/invol": golden_contrib / invol_contrib,
    "golden/(order3+invol)": golden_contrib / (order3_contrib + invol_contrib),
    "(golden-invol)/order3": (golden_contrib - invol_contrib) / order3_contrib,
    "golden*invol/order3^2": golden_contrib * invol_contrib / order3_contrib**2,
}

# Also try spectral quantities
zeta_rho_1 = mpf(0)
for n, lam, m in twisted_spectrum:
    if lam > 0 and m > 0:
        zeta_rho_1 += mpf(m) / mpf(lam)

candidates["Z_rho(1)"] = zeta_rho_1
candidates["1/Z_rho(1)"] = 1/zeta_rho_1

# The Green's trace is 137/15 for untwisted. What about twisted?
green_twisted = mpf(0)
for n in range(1, N_MAX + 1):
    lam = n * (n + 2)
    m = mult_rho(n)
    if m > 0:
        green_twisted += mpf(m) / mpf(lam)

candidates["G_rho (Green trace)"] = green_twisted
candidates["G_rho / (137/15)"] = green_twisted / (mpf(137)/15)
candidates["G_rho * 15/137"] = green_twisted * 15 / 137

# The untwisted trace is 137/15. What's the ratio?
green_untwisted = mpf(0)
for n in range(1, N_MAX + 1):
    lam = n * (n + 2)
    m = mult_untwisted(n)
    if m > 0:
        green_untwisted += mpf(m) / mpf(lam)

candidates["G_M (untwisted)"] = green_untwisted
candidates["G_rho / G_M"] = green_twisted / green_untwisted
candidates["(G_rho - G_M) / G_M"] = (green_twisted - green_untwisted) / green_untwisted

# Weight ratios
candidates["golden_weight * 120"] = golden_weight * 120
candidates["invol_weight * 120"] = invol_weight * 120
candidates["order3_weight * 120"] = order3_weight * 120

print(f"\n--- Candidate alpha-related quantities ---")
print(f"{'Quantity':>35} {'Value':>25} {'Ratio to alpha':>18} {'Ratio to 1/alpha':>18}")
print("-" * 100)

for label, val in candidates.items():
    val_f = float(mpre(val)) if isinstance(val, mpc) else float(val)
    if abs(val_f) > 1e-50:
        ratio_alpha = val_f / float(alpha_target)
        ratio_inv = val_f / float(alpha_inv_target)
        # Flag if close to simple rational multiple
        flag = ""
        for num in range(1, 20):
            for den in range(1, 20):
                if abs(ratio_alpha - num/den) < 0.01:
                    flag = f"  <-- ~{num}/{den} * alpha"
                if abs(ratio_inv - num/den) < 0.01:
                    flag = f"  <-- ~{num}/{den} * (1/alpha)"
                if abs(val_f * float(alpha_inv_target) - num/den) < 0.01:
                    flag = f"  <-- alpha * {num}/{den}"
        print(f"{label:>35} {val_f:>25.15f} {ratio_alpha:>18.8f} {ratio_inv:>18.8f}{flag}")


# ============================================================================
# PART F: FULL SPECTRAL COMPARISON
# ============================================================================

print("\n" + "=" * 80)
print("  PART F: SPECTRAL DATA SUMMARY AND STRUCTURAL ANALYSIS")
print("=" * 80)

# Decomposition of rho under restriction to cyclic subgroups
print("\n--- Representation theory check ---")
print("  rho is the DEFINING 2D representation of SU(2) restricted to 2I.")
print("  It corresponds to the n=1 irrep of SU(2) (the spin-1/2 rep).")
print(f"  chi_rho = chi_1 (the character of the fundamental 2D rep).")
print(f"  dim(rho) = {int(float(chi_rho(mpf(0))))}")

# Verify: mult_rho(n) should equal the multiplicity of the trivial rep
# in V_n tensor rho*, times (n+1). Since rho is the n=1 irrep,
# V_n tensor V_1 = V_{n+1} + V_{n-1} (Clebsch-Gordan).
# So mult_rho(n) = (1/120) sum chi_1* chi_n * (n+1).
# This equals (n+1) * (number of times trivial rep appears in V_n tensor V_1*).

print("\n--- Clebsch-Gordan verification ---")
print("  V_n tensor V_1 = V_{n+1} + V_{n-1} (for n >= 1)")
print("  So rho appears in V_n iff the trivial rep appears in V_n x V_1*")
print("  iff V_n contains V_1 upon restriction to 2I")

# Check first few
print(f"\n{'n':>4} {'mult_rho(n)':>12} {'mult_M(n)':>10} {'mult_M(n-1)':>12} {'mult_M(n+1)':>12}")
print("-" * 55)
for n in range(0, min(30, N_MAX + 1)):
    mr = mult_rho(n)
    mu = mult_untwisted(n)
    mu_prev = mult_untwisted(n - 1) if n > 0 else 0
    mu_next = mult_untwisted(n + 1) if n < N_MAX else 0
    if mr != 0 or mu != 0:
        print(f"{n:>4} {mr:>12} {mu:>10} {mu_prev:>12} {mu_next:>12}")

# Heat kernel comparison
print("\n--- Heat kernel for rho-twisted Laplacian ---")
print("  K_rho(t) = sum_n mult_rho(n) * exp(-n(n+2)*t)")

for t_val in [0.01, 0.1, 0.5, 1.0]:
    K_rho = mpf(0)
    K_M = mpf(0)
    for n in range(N_MAX + 1):
        lam = n * (n + 2)
        mr = mult_rho(n)
        mu = mult_untwisted(n)
        exp_val = mpexp(-mpf(lam) * mpf(t_val))
        K_rho += mr * exp_val
        K_M += mu * exp_val
    ratio = K_rho / K_M if K_M != 0 else mpf(0)
    print(f"  K_rho({t_val}) = {nstr(K_rho, 18)},  K_M({t_val}) = {nstr(K_M, 18)},  "
          f"ratio = {nstr(ratio, 12)}")

# ============================================================================
# FINAL STRUCTURAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("  STRUCTURAL SUMMARY")
print("=" * 80)

print(f"""
  Poincare Dodecahedral Space M = S^3/2I
  Binary Icosahedral Group |2I| = 120

  rho = defining 2D representation (spin-1/2, the n=1 SU(2) irrep)

  Spectrum summary (n=0..{N_MAX}):
    - Total nonzero twisted multiplicities: {len(all_nonzero)} values of n
    - First nonzero n values: {nonzero_ns[:10]}
    - Green's trace (twisted):   G_rho = {nstr(green_twisted, 18)}
    - Green's trace (untwisted): G_M   = {nstr(green_untwisted, 18)}
    - Ratio G_rho/G_M = {nstr(green_twisted/green_untwisted, 15)}
    - 137/15 = {nstr(mpf(137)/15, 18)}
    - G_M - 137/15 = {nstr(green_untwisted - mpf(137)/15, 15)} (convergence tail)

  Ruelle zeta at key points:
    - R_rho(1)   = {nstr(R1, 15)}
    - R_rho(1/2) = {nstr(R_half, 15)}

  Three-class decomposition of -log R_rho(1):
    - Golden:     {nstr(golden_contrib, 15)} ({nstr(100*golden_contrib/(golden_contrib+order3_contrib+invol_contrib+central_contrib), 6)}%)
    - Order-3:    {nstr(order3_contrib, 15)} ({nstr(100*order3_contrib/(golden_contrib+order3_contrib+invol_contrib+central_contrib), 6)}%)
    - Involution: {nstr(invol_contrib, 15)} ({nstr(100*invol_contrib/(golden_contrib+order3_contrib+invol_contrib+central_contrib), 6)}%)
    - Central:    {nstr(central_contrib, 15)} ({nstr(100*central_contrib/(golden_contrib+order3_contrib+invol_contrib+central_contrib), 6)}%)
""")

print("=" * 80)
print("  COMPUTATION COMPLETE")
print("=" * 80)
