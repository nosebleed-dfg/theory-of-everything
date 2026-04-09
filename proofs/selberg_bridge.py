"""
SELBERG_BRIDGE — Laplacian spectrum on S^3/2I; Selberg trace formula by A5 conjugacy classes; Tr(L0^-1)=137/15
nos3bl33d
"""

import math
from mpmath import mp, mpf, mpc, pi as mpi, sin as mpsin, exp as mpexp, log as mplog, fac
from mpmath import zeta as mpzeta, gamma as mpgamma, nstr, power as mppow
import numpy as np

mp.dps = 50  # 50 decimal digits precision

# ============================================================
# SECTION 1: Conjugacy classes of 2I in SU(2)
# ============================================================

# Each class: (name, size, half_angle theta)
CONJ_CLASSES = [
    ("I",   1,  mpf(0)),
    ("-I",  1,  mpi),
    ("C1", 12,  mpi / 5),
    ("C2", 12,  2 * mpi / 5),
    ("C3", 12,  3 * mpi / 5),
    ("C4", 12,  4 * mpi / 5),
    ("C5", 20,  mpi / 3),
    ("C6", 20,  2 * mpi / 3),
    ("C7", 30,  mpi / 2),
]

assert sum(c[1] for c in CONJ_CLASSES) == 120, "Group order check failed"


def chi_n(n, theta):
    """Character of the (n+1)-dim SU(2) irrep at element with half-angle theta.
    chi_n(g) = sin((n+1)*theta) / sin(theta) for theta != 0, pi.
    Special cases: theta=0 -> n+1, theta=pi -> (-1)^n * (n+1).
    """
    n1 = n + 1
    if theta == 0:
        return mpf(n1)
    if abs(theta - mpi) < mpf('1e-30'):
        return mpf((-1)**n * n1)
    return mpsin(n1 * theta) / mpsin(theta)


def multiplicity_on_M(n):
    """Multiplicity of eigenvalue lambda_n = n(n+2) on M = S^3/2I.
    mult_M(n) = (1/120) * sum_C |C| * chi_n(g_C) * (n+1)
    """
    n1 = n + 1
    total = mpf(0)
    for name, size, theta in CONJ_CLASSES:
        total += size * chi_n(n, theta) * n1
    result = total / 120
    return int(round(float(result)))


# ============================================================
# SECTION 2: Compute spectrum for n = 0..200
# ============================================================

N_MAX = 200

print("=" * 80)
print("SPECTRUM OF THE LAPLACIAN ON THE POINCARE DODECAHEDRAL SPACE")
print("M = S^3 / 2I  (Binary Icosahedral Quotient)")
print("=" * 80)

# Precompute all multiplicities
eigenvalues = []   # (n, lambda_n, mult_M)
for n in range(N_MAX + 1):
    lam = n * (n + 2)
    m = multiplicity_on_M(n)
    eigenvalues.append((n, lam, m))

# Sanity: mult should be non-negative integer
for n, lam, m in eigenvalues:
    assert m >= 0, f"Negative multiplicity at n={n}: {m}"

# ============================================================
# SECTION 3a: First 50 eigenvalues with multiplicities
# ============================================================

print("\n" + "-" * 60)
print("FIRST 50 EIGENVALUES (n, lambda_n = n(n+2), mult_M)")
print("-" * 60)
print(f"{'n':>4} {'lambda':>8} {'mult_S3':>8} {'mult_M':>8}")
print("-" * 36)
count = 0
for n, lam, m in eigenvalues:
    if count >= 50:
        break
    mult_s3 = (n + 1) ** 2
    print(f"{n:>4} {lam:>8} {mult_s3:>8} {m:>8}")
    count += 1

# Show nonzero ones
print("\n" + "-" * 60)
print("FIRST 50 NONZERO MULTIPLICITIES")
print("-" * 60)
print(f"{'n':>4} {'lambda':>8} {'mult_M':>8}")
print("-" * 28)
nonzero_count = 0
for n, lam, m in eigenvalues:
    if m > 0 and nonzero_count < 50:
        print(f"{n:>4} {lam:>8} {m:>8}")
        nonzero_count += 1

# ============================================================
# SECTION 3b: Spectral zeta function Z_M(s)
# ============================================================

print("\n" + "-" * 60)
print("SPECTRAL ZETA FUNCTION Z_M(s) = sum_{lambda>0} mult_M(n) / lambda_n^s")
print("-" * 60)

for s_val in [1, 2, 3]:
    z = mpf(0)
    for n, lam, m in eigenvalues:
        if lam > 0 and m > 0:
            z += mpf(m) / mpf(lam) ** s_val
    print(f"  Z_M({s_val}) = {nstr(z, 20)}")

# Higher N for convergence check
z1_extended = mpf(0)
for n in range(1, 2001):
    lam = n * (n + 2)
    m = multiplicity_on_M(n)
    if m > 0:
        z1_extended += mpf(m) / mpf(lam)
print(f"\n  Z_M(1) [N=2000] = {nstr(z1_extended, 20)}")

# ============================================================
# SECTION 3c: Heat kernel trace K_M(t)
# ============================================================

print("\n" + "-" * 60)
print("HEAT KERNEL TRACE K_M(t) = sum_n mult_M(n) * exp(-n(n+2)*t)")
print("-" * 60)

for t_val in [0.01, 0.1, 1.0]:
    K = mpf(0)
    for n, lam, m in eigenvalues:
        if m > 0:
            K += mpf(m) * mpexp(-mpf(lam) * mpf(t_val))
    # Compare with Weyl
    vol_M = mpi**2 / 60
    weyl_leading = vol_M / (4 * mpi * mpf(t_val))**mpf('1.5')
    print(f"  K_M({t_val}) = {nstr(K, 20)}")
    print(f"    Weyl leading term = {nstr(weyl_leading, 20)}")
    print(f"    Ratio K/Weyl = {nstr(K / weyl_leading, 15)}")

# ============================================================
# SECTION 3d: Green's function trace and 137/15 check
# ============================================================

print("\n" + "-" * 60)
print("GREEN'S FUNCTION TRACE: sum_{n>=1} mult_M(n) / (n(n+2))")
print("Expected relation to 137/15 = " + nstr(mpf(137)/15, 20))
print("-" * 60)

green_trace = mpf(0)
for n in range(1, N_MAX + 1):
    lam = n * (n + 2)
    m = multiplicity_on_M(n)
    if m > 0:
        green_trace += mpf(m) / mpf(lam)

print(f"  G_M (N={N_MAX}) = {nstr(green_trace, 25)}")
print(f"  137/15         = {nstr(mpf(137)/15, 25)}")
print(f"  Difference     = {nstr(green_trace - mpf(137)/15, 15)}")
print(f"  Ratio G_M / (137/15) = {nstr(green_trace / (mpf(137)/15), 20)}")

# Extend to N=2000 for better convergence
green_trace_ext = mpf(0)
for n in range(1, 2001):
    lam = n * (n + 2)
    m = multiplicity_on_M(n)
    if m > 0:
        green_trace_ext += mpf(m) / mpf(lam)

print(f"\n  G_M (N=2000) = {nstr(green_trace_ext, 25)}")
print(f"  Difference from 137/15 = {nstr(green_trace_ext - mpf(137)/15, 15)}")

# Try partial fraction: 1/(n(n+2)) = (1/2)(1/n - 1/(n+2))
print("\n  Alternative: sum mult_M(n) * (1/2)(1/n - 1/(n+2)):")
green_alt = mpf(0)
for n in range(1, 2001):
    m = multiplicity_on_M(n)
    if m > 0:
        green_alt += mpf(m) * (mpf(1)/n - mpf(1)/(n+2)) / 2
print(f"  = {nstr(green_alt, 25)}")

# ============================================================
# SECTION 3e: Spectral determinant via zeta regularization
# ============================================================

print("\n" + "-" * 60)
print("SPECTRAL DETERMINANT det'(Delta_M) via zeta regularization")
print("  log det' = -Z_M'(0) = -d/ds Z_M(s)|_{s=0}")
print("-" * 60)

# Numerical derivative of Z_M(s) near s=0
def zeta_M(s):
    """Spectral zeta function, using n up to 2000 for convergence at small s."""
    z = mpf(0)
    for n in range(1, 501):
        lam = n * (n + 2)
        m = multiplicity_on_M(n)
        if m > 0:
            z += mpf(m) / mppow(mpf(lam), s)
    return z

# Use finite difference for derivative
eps = mpf('1e-8')
z_prime_0 = (zeta_M(eps) - zeta_M(-eps)) / (2 * eps)
print(f"  Z_M'(0) approx = {nstr(z_prime_0, 15)}")
print(f"  -Z_M'(0)       = {nstr(-z_prime_0, 15)}")
print(f"  det'(Delta_M) = exp(-Z_M'(0)) = {nstr(mpexp(-z_prime_0), 15)}")

# Also compute Z_M(0) = number of zero modes - 1 (analytically)
z_0 = zeta_M(mpf('1e-15'))
print(f"  Z_M(~0)         = {nstr(z_0, 15)}")

# ============================================================
# SECTION 4: Trace formula decomposition by conjugacy class
# ============================================================

print("\n" + "=" * 80)
print("SELBERG TRACE FORMULA DECOMPOSITION")
print("Test function: h(n) = 1/(n(n+2)) (Green's function)")
print("=" * 80)

# Spectral side
spectral_side = mpf(0)
for n in range(1, 2001):
    m = multiplicity_on_M(n)
    if m > 0:
        spectral_side += mpf(m) / mpf(n * (n + 2))

print(f"\nSPECTRAL SIDE = {nstr(spectral_side, 25)}")

# Geometric side: (1/120) * sum_C |C| * S_C
# where S_C = sum_{n>=1} chi_n(gamma_C) * (n+1) / (n(n+2))

print("\n" + "-" * 60)
print("GEOMETRIC SIDE: Per-class contributions S_C")
print("-" * 60)

N_GEOM = 2000  # Use high N for convergence

class_contributions = {}
for name, size, theta in CONJ_CLASSES:
    S_C = mpf(0)
    for n in range(1, N_GEOM + 1):
        chi = chi_n(n, theta)
        n1 = n + 1
        S_C += chi * n1 / mpf(n * (n + 2))
    class_contributions[name] = S_C
    weighted = size * S_C / 120
    print(f"  {name:>4}: S_C = {nstr(S_C, 20):>30}   |C|*S_C/120 = {nstr(weighted, 20)}")

# Geometric total
geometric_total = mpf(0)
for name, size, theta in CONJ_CLASSES:
    geometric_total += size * class_contributions[name]
geometric_total /= 120

print(f"\nGEOMETRIC SIDE TOTAL = {nstr(geometric_total, 25)}")
print(f"SPECTRAL SIDE        = {nstr(spectral_side, 25)}")
print(f"DIFFERENCE           = {nstr(geometric_total - spectral_side, 15)}")

# ============================================================
# SECTION 5: Group by A5 classes
# ============================================================

print("\n" + "-" * 60)
print("GROUPING BY A5 CLASSES (combining +/- pairs)")
print("-" * 60)

# Identity class: {I} + {-I}
# For -I, chi_n(-I) = (-1)^n * (n+1), so combined contribution picks up only even n
S_identity = mpf(0)
for n in range(1, N_GEOM + 1):
    n1 = n + 1
    # I: chi = n+1, -I: chi = (-1)^n * (n+1)
    combined = (n1 + (-1)**n * n1) * n1 / mpf(n * (n + 2))
    S_identity += combined
S_identity_weighted = S_identity / 120

print(f"\n  IDENTITY ({'{'}I{'}'} + {'{'}-I{'}'}), 2 elements:")
print(f"    Combined S = {nstr(S_identity, 20)}")
print(f"    Weighted (2 * S)/120 = {nstr(S_identity_weighted, 20)}")

# Golden class (5-cycles): C1 + C2 + C3 + C4, 48 elements -> 24 in A5
S_golden = mpf(0)
for n in range(1, N_GEOM + 1):
    n1 = n + 1
    contrib = mpf(0)
    for theta in [mpi/5, 2*mpi/5, 3*mpi/5, 4*mpi/5]:
        contrib += 12 * chi_n(n, theta) * n1
    S_golden += contrib / mpf(n * (n + 2))
S_golden_weighted = S_golden / 120

print(f"\n  GOLDEN (C1+C2+C3+C4), 48 elements -> 24 in A5:")
print(f"    Combined S = {nstr(S_golden, 20)}")
print(f"    Weighted contribution = {nstr(S_golden_weighted, 20)}")

# Order-3 class: C5 + C6, 40 elements -> 20 in A5
S_ord3 = mpf(0)
for n in range(1, N_GEOM + 1):
    n1 = n + 1
    contrib = mpf(0)
    for theta in [mpi/3, 2*mpi/3]:
        contrib += 20 * chi_n(n, theta) * n1
    S_ord3 += contrib / mpf(n * (n + 2))
S_ord3_weighted = S_ord3 / 120

print(f"\n  ORDER-3 (C5+C6), 40 elements -> 20 in A5:")
print(f"    Combined S = {nstr(S_ord3, 20)}")
print(f"    Weighted contribution = {nstr(S_ord3_weighted, 20)}")

# Involution class: C7, 30 elements -> 15 in A5
S_invol = mpf(0)
for n in range(1, N_GEOM + 1):
    n1 = n + 1
    S_invol += 30 * chi_n(n, mpi/2) * n1 / mpf(n * (n + 2))
S_invol_weighted = S_invol / 120

print(f"\n  INVOLUTION (C7), 30 elements -> 15 in A5:")
print(f"    Combined S = {nstr(S_invol, 20)}")
print(f"    Weighted contribution = {nstr(S_invol_weighted, 20)}")

total_check = S_identity_weighted + S_golden_weighted + S_ord3_weighted + S_invol_weighted
print(f"\n  SUM OF ALL A5 CLASSES = {nstr(total_check, 25)}")
print(f"  SPECTRAL SIDE         = {nstr(spectral_side, 25)}")
print(f"  Check difference      = {nstr(total_check - spectral_side, 15)}")

# Percentages
print(f"\n  Percentage contributions:")
print(f"    Identity: {nstr(100 * S_identity_weighted / total_check, 8)}%")
print(f"    Golden:   {nstr(100 * S_golden_weighted / total_check, 8)}%")
print(f"    Order-3:  {nstr(100 * S_ord3_weighted / total_check, 8)}%")
print(f"    Invol:    {nstr(100 * S_invol_weighted / total_check, 8)}%")

# ============================================================
# SECTION 6: Check relationship to d=3, 2d=6, dp=15, V=20
# ============================================================

print("\n" + "-" * 60)
print("CHECKING RELATIONSHIPS TO ICOSAHEDRAL COMBINATORICS")
print("Edges=30, Faces=12(pentagons), Vertices=20")
print("d=3, 2d=6, dp(dual pentagon)=15, V=20")
print("-" * 60)

# The raw per-class S values (unweighted)
# "Golden" per-element: each of the 48 elements contributes S_golden_raw
S_golden_raw_per_class = {}
for name, size, theta in CONJ_CLASSES:
    if name in ["C1", "C2", "C3", "C4"]:
        S_golden_raw_per_class[name] = class_contributions[name]

print(f"\n  Golden class raw S_C values:")
for name in ["C1", "C2", "C3", "C4"]:
    print(f"    {name} (theta={name[-1]}*pi/5): S = {nstr(class_contributions[name], 20)}")

print(f"\n  Order-3 class raw S_C values:")
for name in ["C5", "C6"]:
    theta_label = "pi/3" if name == "C5" else "2pi/3"
    print(f"    {name} (theta={theta_label}): S = {nstr(class_contributions[name], 20)}")

print(f"\n  Involution class:")
print(f"    C7 (theta=pi/2): S = {nstr(class_contributions['C7'], 20)}")

# Compute ratios
print(f"\n  RATIO ANALYSIS:")
# Non-identity weighted contributions
print(f"    S_golden_weighted / S_invol_weighted = {nstr(S_golden_weighted / S_invol_weighted, 15)}")
print(f"    S_ord3_weighted / S_invol_weighted   = {nstr(S_ord3_weighted / S_invol_weighted, 15)}")
print(f"    S_golden_weighted / S_ord3_weighted   = {nstr(S_golden_weighted / S_ord3_weighted, 15)}")

# Check if weights are proportional to class sizes
print(f"\n    Class sizes: Golden=48, Ord3=40, Invol=30")
print(f"    Ratios 48:40:30 = {48/30:.4f} : {40/30:.4f} : 1")
print(f"    Actual ratios:    {nstr(S_golden_weighted / S_invol_weighted, 8)} : {nstr(S_ord3_weighted / S_invol_weighted, 8)} : 1")

# Check vs 6, 15, 20
phi = (1 + mpf(5)**mpf('0.5')) / 2
print(f"\n  Golden ratio phi = {nstr(phi, 20)}")
print(f"  phi^6 = {nstr(phi**6, 20)}")
print(f"  S_golden_weighted = {nstr(S_golden_weighted, 20)}")
print(f"  S_golden_weighted / phi^6 = {nstr(S_golden_weighted / phi**6, 15)}")

print(f"\n  15/(2*pi)^3 = {nstr(mpf(15) / (2*mpi)**3, 20)}")
print(f"  S_invol_weighted = {nstr(S_invol_weighted, 20)}")
print(f"  S_invol_weighted / (15/(2pi)^3) = {nstr(S_invol_weighted / (mpf(15)/(2*mpi)**3), 15)}")

print(f"\n  V = 20")
print(f"  S_ord3_weighted = {nstr(S_ord3_weighted, 20)}")
print(f"  S_ord3_weighted * 120 = {nstr(S_ord3_weighted * 120, 20)}")

# ============================================================
# SECTION 7: Three-class partition function Z_C(s)
# ============================================================

print("\n" + "=" * 80)
print("THREE-CLASS PARTITION FUNCTION")
print("Z_C(s) = sum_{n>=1} chi_n(gamma_C) * (n+1) / (n(n+2))^s")
print("=" * 80)

def Z_class(theta_list, sizes, s_val, N_terms=2000):
    """Compute Z_C(s) for a combined class."""
    total = mpf(0)
    for n in range(1, N_terms + 1):
        n1 = n + 1
        lam = mpf(n * (n + 2))
        chi_sum = mpf(0)
        for theta, sz in zip(theta_list, sizes):
            chi_sum += sz * chi_n(n, theta) * n1
        total += chi_sum / mppow(lam, s_val)
    return total

# Golden class (48 elements total)
thetas_golden = [mpi/5, 2*mpi/5, 3*mpi/5, 4*mpi/5]
sizes_golden = [12, 12, 12, 12]

# Order-3 class (40 elements)
thetas_ord3 = [mpi/3, 2*mpi/3]
sizes_ord3 = [20, 20]

# Involution class (30 elements)
thetas_invol = [mpi/2]
sizes_invol = [30]

for s_val in [mpf(1), mpf(2), mpf(3)]:
    Z_g = Z_class(thetas_golden, sizes_golden, s_val) / 120
    Z_o = Z_class(thetas_ord3, sizes_ord3, s_val) / 120
    Z_i = Z_class(thetas_invol, sizes_invol, s_val) / 120

    print(f"\n  s = {int(s_val)}:")
    print(f"    Z_golden  / 120 = {nstr(Z_g, 20)}")
    print(f"    Z_ord3    / 120 = {nstr(Z_o, 20)}")
    print(f"    Z_invol   / 120 = {nstr(Z_i, 20)}")
    print(f"    Z_identity/ 120 = {nstr(S_identity_weighted if s_val == 1 else 'N/A', 20)}")

    if s_val == 1:
        print(f"\n    CHECK vs icosahedral numbers:")
        print(f"    Z_golden(1) = {nstr(Z_g, 20)}")
        print(f"    phi^6 = {nstr(phi**6, 20)}")
        print(f"    Z_golden(1) / phi = {nstr(Z_g / phi, 15)}")
        print(f"    Z_golden(1) * 120 = {nstr(Z_g * 120, 15)}")

        print(f"\n    Z_invol(1) = {nstr(Z_i, 20)}")
        print(f"    15/(2pi)^3 = {nstr(mpf(15)/(2*mpi)**3, 20)}")
        print(f"    Z_invol(1) * 120 = {nstr(Z_i * 120, 15)}")

        print(f"\n    Z_ord3(1) = {nstr(Z_o, 20)}")
        print(f"    20/120 = {nstr(mpf(20)/120, 15)}")
        print(f"    Z_ord3(1) * 120 = {nstr(Z_o * 120, 15)}")

# ============================================================
# SECTION 8: Selberg zeta function
# ============================================================

print("\n" + "=" * 80)
print("SELBERG ZETA FUNCTION")
print("Z_Selberg(s) = prod_gamma prod_{k>=0} (1 - exp(-(s+k)*l(gamma)))")
print("=" * 80)

# Geodesic lengths: l = 2*theta for each non-identity conjugacy class
# Each primitive geodesic contributes. For 2I, primitive closed geodesics
# correspond to conjugacy classes (excluding identity).
geodesic_data = [
    ("C1", 12, 2*mpi/5),
    ("C2", 12, 4*mpi/5),
    ("C3", 12, 6*mpi/5),
    ("C4", 12, 8*mpi/5),
    ("C5", 20, 2*mpi/3),
    ("C6", 20, 4*mpi/3),
    ("C7", 30, mpi),
]

print("\nGeodesic lengths:")
for name, mult, length in geodesic_data:
    print(f"  {name}: l = {nstr(length, 15)} ({nstr(length/mpi, 10)} * pi), multiplicity {mult}")

# Compute |Z_Selberg(s)| along critical line s = 1/2 + it
print("\n|Z_Selberg(1/2 + it)| for t = 0 to 50:")
print(f"{'t':>8} {'|Z_Selberg|':>20}")
print("-" * 30)

K_MAX = 30  # truncation for product over k

for t_idx in range(51):
    t_float = float(t_idx)
    s_complex = mpc('0.5', str(t_float))

    log_Z = mpc(0, 0)

    for name, mult, length in geodesic_data:
        for k in range(K_MAX + 1):
            arg = -(s_complex + k) * length
            val = 1 - mpexp(arg)
            if abs(val) > mpf('1e-40'):
                log_Z += mult * mplog(val)

    z_abs = abs(mpexp(log_Z))

    if t_idx <= 10 or t_idx % 5 == 0:
        print(f"{t_float:>8.1f} {nstr(z_abs, 15):>20}")

# ============================================================
# SECTION 9: Weyl counting function
# ============================================================

print("\n" + "=" * 80)
print("WEYL COUNTING FUNCTION")
print("N_M(lambda) = #{eigenvalues <= lambda}")
print("Weyl law: N_M(lambda) ~ Vol(M)/(6*pi^2) * lambda^{3/2}")
print("Vol(M) = pi^2/60")
print("=" * 80)

vol_M = mpi**2 / 60
weyl_coeff = vol_M / (6 * mpi**2)
print(f"\nVol(M) = pi^2/60 = {nstr(vol_M, 20)}")
print(f"Weyl coefficient = Vol(M)/(6*pi^2) = 1/360 = {nstr(weyl_coeff, 20)}")
print(f"  (check: 1/360 = {nstr(mpf(1)/360, 20)})")

# Cumulative eigenvalue count
cum_count = 0
print(f"\n{'lambda':>8} {'N_M(lambda)':>12} {'Weyl':>14} {'R(lambda)':>14} {'R/Weyl%':>10}")
print("-" * 65)

checkpoint_lambdas = [3, 8, 15, 24, 35, 48, 63, 80, 99, 120, 143, 168, 195, 224,
                      300, 400, 500, 1000, 2000, 5000, 10000, 20000, 40000]

lam_idx = 0
for n, lam, m in eigenvalues:
    cum_count += m
    while lam_idx < len(checkpoint_lambdas) and lam >= checkpoint_lambdas[lam_idx]:
        check_lam = checkpoint_lambdas[lam_idx]
        # Recount up to this lambda
        actual_count = sum(mm for nn, ll, mm in eigenvalues if ll <= check_lam)
        weyl_est = weyl_coeff * mpf(check_lam) ** mpf('1.5')
        remainder = actual_count - weyl_est
        pct = 100 * remainder / weyl_est if weyl_est > 0 else mpf(0)
        print(f"{check_lam:>8} {actual_count:>12} {nstr(weyl_est, 10):>14} {nstr(remainder, 10):>14} {nstr(pct, 6):>10}")
        lam_idx += 1

# Remainder oscillation analysis
print("\n" + "-" * 60)
print("REMAINDER TERM OSCILLATION ANALYSIS")
print("Looking for periods related to geodesic lengths:")
print(f"  2*pi/l for l=2pi/5: period = 5")
print(f"  2*pi/l for l=pi:    period = 2")
print(f"  2*pi/l for l=2pi/3: period = 3")
print("-" * 60)

remainders = []
for n, lam, m in eigenvalues:
    if lam > 0:
        actual = sum(mm for nn, ll, mm in eigenvalues if ll <= lam)
        weyl_est = weyl_coeff * mpf(lam) ** mpf('1.5')
        remainders.append((n, lam, float(actual - weyl_est)))

# Check periodicity by looking at remainder at n = 0, 5, 10, 15... (period 5)
print("\n  Remainder at n = 0 mod 10 (checking period-5 component):")
for n, lam, r in remainders:
    if n <= 100 and n % 10 == 0:
        print(f"    n={n:>3}, lambda={lam:>6}, R = {r:>10.4f}")

# ============================================================
# SECTION 10: Multiplicity pattern analysis
# ============================================================

print("\n" + "=" * 80)
print("DETAILED MULTIPLICITY PATTERN")
print("=" * 80)

zero_ns = [n for n, lam, m in eigenvalues if m == 0 and n > 0]
nonzero_ns = [n for n, lam, m in eigenvalues if m > 0]

print(f"\n  Zero-multiplicity n values (n=1..{N_MAX}): {len(zero_ns)} values")
print(f"  Non-zero multiplicity values: {len(nonzero_ns)} values")
print(f"  First zero-mult n: {zero_ns[:30]}")
print(f"  Pattern of non-zero n (first 40): {[n for n, _, m in eigenvalues if m > 0][:40]}")

print("\n  Multiplicity sequence mod 60 pattern (period of 2I):")
for start in [0, 12, 20, 24, 30, 32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60]:
    if start <= N_MAX:
        m_val = eigenvalues[start][2]
        print(f"    n={start:>3}: mult = {m_val:>4}   (n+1={start+1})")

# Which n+1 values have trivial 2I representation?
# The representation ring of 2I: trivial appears in V_n iff n+1 appears as
# a dimension in the decomposition of the (n+1)-dim SU(2) rep restricted to 2I
print("\n  n values where mult_M(n) = n+1 (single trivial copy):")
single_copy = [(n, m) for n, _, m in eigenvalues if m == n+1 and n > 0]
print(f"    {[n for n, m in single_copy[:25]]}")

print("\n  n values where mult_M(n) > n+1 (multiple trivial copies):")
multi_copy = [(n, m) for n, _, m in eigenvalues if m > n+1 and n > 0]
print(f"    {[(n, m, m//(n+1)) for n, m in multi_copy[:20]]}")

# ============================================================
# SECTION 11: The 137/15 deep dive
# ============================================================

print("\n" + "=" * 80)
print("DEEP DIVE: THE 137/15 CONNECTION")
print("=" * 80)

# The raw Green's trace diverges because mult_M(n) ~ (n+1)^2/120 for large n.
# So sum mult_M(n)/(n(n+2)) ~ sum (n+1)^2/(120*n(n+2)) ~ sum 1/120 -> diverges.
#
# The 137/15 from the Pythagorean framework is Tr(L0^{-1}) where L0 is the
# FINITE DIMENSIONAL combinatorial Laplacian on the icosahedron/dodecahedron.
#
# REGULARIZATION APPROACH 1: Subtract the S^3 (identity) contribution.
# The "geometric correction" = sum of all non-identity class contributions.
# This CONVERGES because sin((n+1)*theta)/sin(theta) oscillates.

print("\n  APPROACH 1: Non-identity geometric contribution (CONVERGES)")
print("  G_geom = (1/120) * sum_{C != I,-I} |C| * S_C")

# Already computed above as golden + ord3 + invol
G_geom = S_golden_weighted + S_ord3_weighted + S_invol_weighted
print(f"  G_geom = {nstr(G_geom, 25)}")
print(f"  G_geom = Golden + Ord3 + Invol")
print(f"         = {nstr(S_golden_weighted, 15)} + ({nstr(S_ord3_weighted, 15)}) + ({nstr(S_invol_weighted, 15)})")
print(f"  G_geom * 120 = {nstr(G_geom * 120, 20)}")

# APPROACH 2: Regularized spectral trace via zeta regularization
# Z_M(s) = sum mult_M(n) / (n(n+2))^s, analytically continue to s=1
# For s>1 it converges; at s=1 it diverges. But we can subtract the pole.
print("\n  APPROACH 2: Zeta-regularized trace")
print("  Z_M(s) - (1/120)*Z_{S^3}(s) as s -> 1")

# Z_{S^3}(s) = sum_{n>=1} (n+1)^2 / (n(n+2))^s
# (n+1)^2/(n(n+2))^s = (n+1)^2 / ((n+1)^2-1)^s
# ~ 1/(n+1)^{2s-2} for large n
# So Z_{S^3}(s) has a pole at s=1 with residue related to zeta(2s-2)

# Compute regularized value: subtract identity-class divergence
# G_reg = lim_{s->1} [Z_M(s) - (1/120)*Z_{S^3}(s)]
# = (1/120) * sum_{C != I} |C| * S_C(s=1)
# = G_geom (since -I also diverges separately, let's separate it)

# Actually {I}: S_I = sum (n+1)^2/(n(n+2)) DIVERGES
# {-I}: S_{-I} = sum (-1)^n*(n+1)^2/(n(n+2))
#   = sum even: 2*(n+1)^2/(n(n+2)) from {I}+{-I} -> still diverges (n even terms)
# NO WAIT: S_{-I} by itself = sum_{n>=1} (-1)^n*(n+1)^2/(n(n+2))
# = sum_{n>=1} (-1)^n * (n+1)^2/((n+1)^2-1)
# = sum_{n>=1} (-1)^n * [1 + 1/((n+1)^2-1)]
# = sum (-1)^n + sum (-1)^n/(n(n+2))
# First sum = -1+1-1+1-... = Cesaro sum -1/2
# Second sum converges by alternating series test

# So {-I} alone: the CONVERGENT part = sum_{n>=1} (-1)^n / (n(n+2))
# = (1/2) sum (-1)^n * (1/n - 1/(n+2))
S_negI_convergent = mpf(0)
for n in range(1, 10001):
    S_negI_convergent += mpf((-1)**n) / mpf(n * (n + 2))
print(f"\n  S_{{-I}} convergent part = sum (-1)^n/(n(n+2)) = {nstr(S_negI_convergent, 25)}")
print(f"  Expected: (1/2)(sum(-1)^n/n - sum(-1)^n/(n+2))")
# sum (-1)^n/n for n>=1 = -ln(2)
# sum (-1)^n/(n+2) for n>=1 = sum (-1)^{m-2}/m for m>=3 = sum (-1)^m/m for m>=3
#   = -ln(2) - (-1 + 1/2) = -ln(2) + 1/2
# So: (1/2)(-ln(2) - (-ln(2)+1/2)) = (1/2)(-1/2) = -1/4
from mpmath import ln
print(f"  Exact value: -1/4 = {nstr(mpf(-1)/4, 25)}")
print(f"  Computed:          {nstr(S_negI_convergent, 25)}")
print(f"  Match: {abs(float(S_negI_convergent + mpf(1)/4)) < 1e-10}")

# APPROACH 3: The non-identity trace IS the geometric content.
# This equals the sum of all non-trivial conjugacy class contributions.
# From Selberg trace formula, this is an orbital integral / length spectrum sum.
print("\n  APPROACH 3: Geometric side = non-identity orbital integrals")
print("  These encode the geometry (geodesic lengths) of M")

# Each non-identity S_C with h(n) = 1/(n(n+2)):
# S_C = sum_{n>=1} chi_n(g_C) * (n+1) / (n(n+2))
# = sum_{n>=1} sin((n+1)*theta)/sin(theta) * (n+1) / (n(n+2))
#
# Using partial fractions: (n+1)/(n(n+2)) = (1/2)(1/n + 1/(n+2))
# So S_C = (1/(2sin(theta))) * sum_{n>=1} sin((n+1)*theta) * (1/n + 1/(n+2))

# Known closed form: sum_{n>=1} sin(n*alpha)/n = (pi-alpha)/2 for 0 < alpha < 2pi
# So sum_{n>=1} sin((n+1)*alpha)/n = sum_{m>=2} sin(m*alpha)/(m-1)
# This relates to Clausen function and digamma values.

# Let's compute these to high precision and try rational identification
mp.dps = 80
print("\n  HIGH PRECISION non-identity S_C values:")
for name, size, theta in CONJ_CLASSES:
    if name in ["I", "-I"]:
        continue
    S_C = mpf(0)
    for n in range(1, 20001):
        chi = chi_n(n, theta)
        n1 = n + 1
        S_C += chi * n1 / mpf(n * (n + 2))
    print(f"    {name:>4} (theta/pi={nstr(theta/mpi, 8)}, |C|={size:>2}): S_C = {nstr(S_C, 40)}")
    # Try to identify
    # Check S_C * sin(theta) -- this might be cleaner
    st = mpsin(theta)
    print(f"          S_C * sin(theta) = {nstr(S_C * st, 40)}")

mp.dps = 50

# ============================================================
# SECTION 12: The C7 (involution) class -- exact value
# ============================================================

print("\n" + "=" * 80)
print("EXACT ANALYSIS: INVOLUTION CLASS C7 (theta = pi/2)")
print("=" * 80)

# For theta = pi/2:
# chi_n(pi/2) = sin((n+1)*pi/2)/sin(pi/2) = sin((n+1)*pi/2)
# This is 0 for n even, and for n odd:
#   n=1: sin(pi) = 0
#   n=3: sin(2pi) = 0
# Wait: sin((n+1)*pi/2):
#   n=0: sin(pi/2) = 1
#   n=1: sin(pi) = 0
#   n=2: sin(3pi/2) = -1
#   n=3: sin(2pi) = 0
#   n=4: sin(5pi/2) = 1
# Pattern: chi_n = 0 for n odd, chi_n = (-1)^{n/2} for n even

# So S_{C7} = sum_{k>=1} (-1)^k * (2k+1) / (2k(2k+2))
# = sum_{k>=1} (-1)^k * (2k+1) / (4k(k+1))
# = (1/4) * sum_{k>=1} (-1)^k * (2k+1)/(k(k+1))
# = (1/4) * sum (-1)^k * (2/1 + 1/(k(k+1)))  no...
# (2k+1)/(k(k+1)) = 2k/(k(k+1)) + 1/(k(k+1)) = 2/(k+1) + 1/k - 1/(k+1)
#                  = 1/k + 1/(k+1)

# S_{C7} = (1/4) * sum_{k>=1} (-1)^k * (1/k + 1/(k+1))
# = (1/4) * [sum (-1)^k/k + sum (-1)^k/(k+1)]
# = (1/4) * [-ln(2) + sum_{m>=2} (-1)^{m-1}/m]
# = (1/4) * [-ln(2) + (-ln(2) + 1)]  since sum_{m>=2}(-1)^{m-1}/m = ln(2) - 1
# Wait: sum_{k>=1} (-1)^k/k = -ln(2)
# sum_{k>=1} (-1)^k/(k+1) = sum_{m>=2} (-1)^{m-1}/m = -(sum_{m>=2} (-1)^m/m)
#   = -(sum_{m>=1}(-1)^m/m - (-1)) = -(-ln(2)+1) = ln(2)-1

# So S_{C7} = (1/4)(-ln(2) + ln(2) - 1) = -1/4

print(f"  Analytic: S_{{C7}} = -1/4 (EXACT)")
print(f"  Computed: S_{{C7}} = {nstr(class_contributions['C7'], 30)}")
print(f"  Match: {abs(float(class_contributions['C7'] + mpf(1)/4)) < 1e-10}")

print(f"\n  Weighted: 30 * (-1/4) / 120 = -30/480 = -1/16")
print(f"  Computed weighted: {nstr(S_invol_weighted, 30)}")
# Hmm that's -0.0624... which is -1/16.01... Let me check
print(f"  -1/16 = {nstr(mpf(-1)/16, 20)}")
# -1/16 = -0.0625, computed = -0.06244 -- not quite
# Let me recheck: 30 * S_C7 / 120 where S_C7 = sum_{n>=1} chi_n * (n+1)/(n(n+2))
# chi_n(pi/2) * (n+1) = sin((n+1)pi/2)/sin(pi/2) * (n+1) = sin((n+1)pi/2) * (n+1)
# For n even = 2k: chi = (-1)^k, so chi*(n+1) = (-1)^k * (2k+1)
# S_{C7} = sum_{k>=1} (-1)^k * (2k+1) / (2k*(2k+2))

# Numerically:
S_C7_check = mpf(0)
for k in range(1, 10001):
    S_C7_check += mpf((-1)**k) * (2*k+1) / mpf(2*k*(2*k+2))
print(f"\n  Direct sum S_C7 (k=1..10000): {nstr(S_C7_check, 30)}")
print(f"  class_contributions['C7']:     {nstr(class_contributions['C7'], 30)}")

# (2k+1)/(2k(2k+2)) = (2k+1)/(4k(k+1))
# = (1/4)*(1/k + 1/(k+1))
# Nope, let me redo: (2k+1)/(4k(k+1))
# 2k+1 = k + (k+1), so (2k+1)/(4k(k+1)) = 1/(4(k+1)) + 1/(4k)
# Yes: S_{C7} = (1/4) * sum_{k>=1} (-1)^k * (1/k + 1/(k+1))
s1 = mpf(0)
s2 = mpf(0)
for k in range(1, 100001):
    s1 += mpf((-1)**k) / k
    s2 += mpf((-1)**k) / (k+1)
print(f"\n  sum (-1)^k/k = {nstr(s1, 25)} (should be -ln(2) = {nstr(-mplog(2), 25)})")
print(f"  sum (-1)^k/(k+1) = {nstr(s2, 25)} (should be ln(2)-1 = {nstr(mplog(2)-1, 25)})")
print(f"  (1/4)(s1+s2) = {nstr((s1+s2)/4, 25)}")
print(f"  (1/4)(-ln2 + ln2-1) = -1/4 = {nstr(mpf(-1)/4, 20)}")

# OK so the EXACT answer for S_{C7} IS -1/4.
# Weighted: 30*(-1/4)/120 = -30/480 = -1/16
# But computed value is -0.249750... which = -1/4 * (1 - 1/(N+1)(N+2))... truncation error
# With N=2000: -0.24975... yes, the tail decays as 1/N

# Weighted with exact S: 30 * (-1/4) / 120 = -1/16 = -0.0625 EXACTLY!
print(f"\n  EXACT weighted C7: -1/16 = {nstr(mpf(-1)/16, 20)}")

# ============================================================
# SECTION 13: Exact values for all classes
# ============================================================

print("\n" + "=" * 80)
print("EXACT / HIGH-PRECISION VALUES FOR ALL NON-IDENTITY CLASSES")
print("=" * 80)

# For theta = pi/3:
# chi_n(pi/3) = sin((n+1)*pi/3)/sin(pi/3)
# Period 6 in n: chi values cycle through 0,1,1,0,-1,-1 (times sqrt(3)/2 factor)
# Actually sin((n+1)*pi/3)/sin(pi/3):
#   n=0: sin(pi/3)/sin(pi/3) = 1
#   n=1: sin(2pi/3)/sin(pi/3) = 1
#   n=2: sin(pi)/sin(pi/3) = 0
#   n=3: sin(4pi/3)/sin(pi/3) = -1
#   n=4: sin(5pi/3)/sin(pi/3) = -1
#   n=5: sin(2pi)/sin(pi/3) = 0
# So chi_n has period 6: 1, 1, 0, -1, -1, 0

# For theta = 2pi/3:
# chi_n(2pi/3) = sin((n+1)*2pi/3)/sin(2pi/3)
#   n=0: 1, n=1: -1, n=2: 0, n=3: 1, n=4: -1, n=5: 0
# Period 3 effectively: 1, -1, 0

# For theta = pi/5:
# chi_n(pi/5) = sin((n+1)*pi/5)/sin(pi/5)
# Period 10: uses Chebyshev-like values involving phi

# Let me compute C5 and C6 exactly using the periodicity
# C5 (theta=pi/3): chi pattern (period 6): [1, 1, 0, -1, -1, 0] for n=0,1,2,3,4,5
# S_{C5} = sum_{n>=1} chi_n * (n+1) / (n(n+2))
# n=1: chi=1, contrib = 2/(1*3) = 2/3
# n=2: chi=0
# n=3: chi=-1, contrib = -4/(3*5) = -4/15
# n=4: chi=-1, contrib = -5/(4*6) = -5/24
# n=5: chi=0
# n=6k+1: chi=1, n=6k+3: chi=-1, n=6k+4: chi=-1

# Use digamma for exact evaluation
from mpmath import digamma, psi

# Generic formula: sum_{n=0}^{infty} sin((n+a)*theta) / (n+b)
# involves digamma at rational arguments

# For C7 we proved S = -1/4 exactly
# For the others, let's compute to extreme precision

mp.dps = 100
print("\n  Computing to 100 digits of precision (N=50000):")

exact_S = {}
for name, size, theta in CONJ_CLASSES:
    if name in ["I", "-I"]:
        continue
    S_C = mpf(0)
    for n in range(1, 50001):
        chi = chi_n(n, theta)
        n1 = n + 1
        S_C += chi * n1 / mpf(n * (n + 2))
    exact_S[name] = S_C
    print(f"    {name:>4}: S_C = {nstr(S_C, 50)}")

# Compute the full non-identity geometric contribution with exact C7
G_geom_exact = mpf(0)
for name, size, theta in CONJ_CLASSES:
    if name in ["I", "-I"]:
        continue
    G_geom_exact += size * exact_S[name] / 120
print(f"\n  G_geom (non-identity, N=50000) = {nstr(G_geom_exact, 40)}")

# Also compute: what is S_golden + S_ord3 + S_invol (UNweighted)?
S_nonid_unweighted = mpf(0)
for name, size, theta in CONJ_CLASSES:
    if name in ["I", "-I"]:
        continue
    S_nonid_unweighted += exact_S[name]
print(f"  Sum of raw S_C values = {nstr(S_nonid_unweighted, 40)}")

mp.dps = 50

# ============================================================
# SECTION 14: The 137/15 bridge -- finite operator trace
# ============================================================

print("\n" + "=" * 80)
print("THE 137/15 BRIDGE")
print("=" * 80)

# Tr(L0^{-1}) = 137/15 was proven for the 10x10 Gram matrix Laplacian
# of the Pythagorean framework. That is a FINITE combinatorial object.
#
# The connection to the Poincare dodecahedral space spectrum:
# The icosahedron has the same symmetry group A5 = 2I/{+/-I}
#
# Can we find 137/15 in the spectral data?
#
# Idea 1: Look at the first few eigenvalues of M only (finite truncation)
# Idea 2: Weighted sum over A5 conjugacy class contributions
# Idea 3: Specific linear combination of the class sums

print("\n  Test 1: Partial Green's trace up to n=12 (first nonzero eigenvalue)")
partial_12 = mpf(0)
for n, lam, m in eigenvalues:
    if 0 < n <= 12 and m > 0:
        partial_12 += mpf(m) / mpf(lam)
print(f"    Sum up to n=12: {nstr(partial_12, 20)}")
print(f"    This is 13/168 = {nstr(mpf(13)/168, 20)}")

print("\n  Test 2: Non-identity geometric side * 120")
print(f"    G_geom * 120 = {nstr(G_geom * 120, 20)}")
print(f"    137/15 * 120 = {nstr(mpf(137)*120/15, 10)}")

print("\n  Test 3: Ratios involving 137/15")
print(f"    G_geom / (137/15) = {nstr(G_geom / (mpf(137)/15), 15)}")
print(f"    120 * G_geom / (137/15) = {nstr(120 * G_geom / (mpf(137)/15), 15)}")

# What ARE the non-identity sums as exact fractions?
# C7: S = -1/4 exactly -> weighted = -30/(4*120) = -1/16
# C5: has period 6 pattern -> can compute exactly with Euler digamma
# C1-C4: period 10 pattern

# Test: is the sum of weighted non-identity = -137/(15*120)?
print(f"\n  Test 4: G_geom vs -137/(15*120) = {nstr(mpf(-137)/(15*120), 20)}")
print(f"          G_geom                  = {nstr(G_geom, 20)}")
print(f"    No match.")

print(f"\n  Test 5: G_geom vs simple fractions")
# Try to identify G_geom
g_val = float(G_geom)
print(f"    G_geom = {g_val}")
# Search
from fractions import Fraction
best_frac = None
best_err = 1e10
for d in range(1, 10000):
    n_approx = round(g_val * d)
    err = abs(g_val - n_approx/d)
    if err < best_err:
        best_err = err
        best_frac = Fraction(n_approx, d)
    if err < 1e-12:
        print(f"    Close fraction: {n_approx}/{d} = {n_approx/d:.15f}, err = {err:.2e}")
        break
print(f"    Best fraction found: {best_frac} = {float(best_frac):.15f}")

# ============================================================
# SECTION 15: Final summary
# ============================================================

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"\n  Eigenvalue 0: multiplicity {eigenvalues[0][2]}")
print(f"  (Connected components = 1)")

first_nonzero = [(n, lam, m) for n, lam, m in eigenvalues if m > 0 and lam > 0]
if first_nonzero:
    print(f"\n  First nonzero eigenvalue: lambda = {first_nonzero[0][1]} (n={first_nonzero[0][0]})")
    print(f"    with multiplicity {first_nonzero[0][2]}")
    print(f"    (On S^3: lambda=3 with mult=4)")

print(f"\n  SPECTRAL GAP: {first_nonzero[0][1]}")
print(f"  (On S^3 the gap is 3; on M it is {first_nonzero[0][1]}, MUCH larger)")
print(f"  Ratio gap(M)/gap(S^3) = {first_nonzero[0][1]/3:.1f}")

total_with_mult = sum(m for _, _, m in eigenvalues)
total_s3 = sum((n+1)**2 for n in range(N_MAX + 1))
print(f"\n  Total eigenvalues (with mult) up to n={N_MAX}:")
print(f"    On M:  {total_with_mult}")
print(f"    On S^3: {total_s3}")
print(f"    Ratio: {total_with_mult/total_s3:.6f} (expected 1/120 = {1/120:.6f})")

count_nonzero = sum(1 for n in range(1, N_MAX+1) if multiplicity_on_M(n) > 0)
print(f"\n  Non-zero multiplicities in n=1..{N_MAX}: {count_nonzero}/{N_MAX}")
print(f"  Fraction with eigenvalues: {count_nonzero/N_MAX:.4f}")

print(f"\n  EXACT RESULTS:")
print(f"    S_{{C7}} (involution, theta=pi/2) = -1/4 EXACTLY")
print(f"    Weighted C7 contribution = -1/16 = -0.0625 EXACTLY")
print(f"    S_{{-I}} convergent part = -1/4 EXACTLY")

print(f"\n  Geometric decomposition (non-identity contributions):")
print(f"    Golden (5-fold):  {nstr(S_golden_weighted, 20)}")
print(f"    Order-3:          {nstr(S_ord3_weighted, 20)}")
print(f"    Involution:       -1/16 = -0.0625")
print(f"    Total non-identity: {nstr(G_geom, 20)}")

# The Selberg trace formula relates spectral data to geometric data.
# The geometric side is controlled by geodesic lengths = 2*theta.
# The identity divergence = the volume term in the trace formula.

print(f"\n  SPECTRAL ZETA VALUES:")
for s_val in [1, 2, 3]:
    z = mpf(0)
    for n, lam, m in eigenvalues:
        if lam > 0 and m > 0:
            z += mpf(m) / mpf(lam) ** s_val
    print(f"    Z_M({s_val}) [N={N_MAX}] = {nstr(z, 20)}")

print(f"\n  HEAT KERNEL at t=0.01: K_M(0.01) = {nstr(mpf(0), 15)}")
K_001 = mpf(0)
for n, lam, m in eigenvalues:
    if m > 0:
        K_001 += mpf(m) * mpexp(-mpf(lam) * mpf('0.01'))
print(f"    (recomputed) = {nstr(K_001, 15)}")

print("\n" + "=" * 80)
print("COMPUTATION COMPLETE")
print("=" * 80)
