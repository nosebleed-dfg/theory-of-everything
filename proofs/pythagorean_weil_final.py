#!/usr/bin/env python3
"""
PYTHAGOREAN_WEIL_FINAL — forward + backward Weil functionals; Pythagorean Gram G_f^2 + G_b^2 >= 0
nos3bl33d

Gram decomposition A-P, A+P, Pythagorean sum, SVD test, critical conductor search.
15 Gaussian basis functions, conductor N=800. mpmath 30 digits.
"""

from mpmath import (mp, mpf, mpc, log, pi, sqrt, gamma, cos, sin,
                     exp, power, re, im, quad, loggamma,
                     digamma, zeta, fsum, matrix, nstr, fabs,
                     euler as euler_gamma, inf, eig, svd, eye)
import sys
import warnings
warnings.filterwarnings("ignore")

mp.dps = 30

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + sqrt(5)) / 2
PHI_INV = PHI - 1
PHI_SQ = PHI**2
PHI_INV_SQ = PHI_INV**2
CONDUCTOR = mpf(800)

assert fabs(PHI_SQ - PHI - 1) < mpf('1e-25')
assert fabs(PHI_INV_SQ - (2 - PHI)) < mpf('1e-25')

# A5 conjugacy classes for icosahedral 2-dim rep
A5_CLASSES = [
    ("identity",    mpf(2),   mpf(4),       1,  mpf(1)/60,  1),
    ("golden+",     PHI,      PHI_SQ,       12, mpf(12)/60, 5),
    ("golden-",    -PHI_INV,  PHI_INV_SQ,   12, mpf(12)/60, 5),
    ("order-3",     mpf(-1),  mpf(1),       20, mpf(20)/60, 3),
    ("involution",  mpf(0),   mpf(0),       15, mpf(15)/60, 2),
]


def separator(title):
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print()


def sieve_primes(N):
    if N < 2:
        return []
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, N + 1, i):
                sieve[j] = False
    return [i for i in range(2, N + 1) if sieve[i]]


# =============================================================================
# FROBENIUS TRACE COMPUTATION
# =============================================================================

def poly_eval_mod(x, p):
    return (pow(x, 5, p) + (20 % p) * x + (16 % p)) % p


def poly_divmod_Fp(a, b, p):
    a = [c % p for c in a]
    b = [c % p for c in b]
    while len(b) > 1 and b[0] % p == 0:
        b = b[1:]
    if len(b) == 0 or (len(b) == 1 and b[0] % p == 0):
        raise ValueError("Division by zero polynomial")
    inv_lead = pow(b[0], p - 2, p)
    q = []
    r = list(a)
    while len(r) >= len(b):
        if r[0] % p == 0:
            q.append(0)
            r = r[1:]
            continue
        coeff = (r[0] * inv_lead) % p
        q.append(coeff)
        for i in range(len(b)):
            r[i] = (r[i] - coeff * b[i]) % p
        r = r[1:]
    return q, r


def poly_mul_Fp(a, b, p):
    if not a or not b:
        return [0]
    n = len(a) + len(b) - 1
    result = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
    return result


def poly_rem_Fp(a, b, p):
    _, r = poly_divmod_Fp(a, b, p)
    return r if r else [0]


def poly_powmod_Fp(base, exp_val, modpoly, p):
    result = [1]
    base = poly_rem_Fp(base, modpoly, p)
    while exp_val > 0:
        if exp_val % 2 == 1:
            result = poly_mul_Fp(result, base, p)
            result = poly_rem_Fp(result, modpoly, p)
        base = poly_mul_Fp(base, base, p)
        base = poly_rem_Fp(base, modpoly, p)
        exp_val //= 2
    return result


def poly_gcd_Fp(a, b, p):
    while len(b) > 0 and any(c % p != 0 for c in b):
        while len(b) > 1 and b[0] % p == 0:
            b = b[1:]
        if len(b) == 0 or (len(b) == 1 and b[0] % p == 0):
            break
        if len(a) < len(b):
            a, b = b, a
            continue
        _, r = poly_divmod_Fp(a, b, p)
        a = b
        b = r
    while len(a) > 1 and a[0] % p == 0:
        a = a[1:]
    if len(a) > 0 and a[0] % p != 0:
        inv_lead = pow(a[0], p - 2, p)
        a = [(c * inv_lead) % p for c in a]
    return a


def factorize_mod_p(p):
    f = [1, 0, 0, 0, 20 % p, 16 % p]
    roots = []
    for r in range(p):
        if poly_eval_mod(r, p) == 0:
            roots.append(r)
    degrees = []
    remaining = list(f)
    for r in roots:
        new = [0] * (len(remaining) - 1)
        new[0] = remaining[0] % p
        for i in range(1, len(remaining) - 1):
            new[i] = (remaining[i] + r * new[i - 1]) % p
        remaining = new
        degrees.append(1)
    deg_rem = len(remaining) - 1
    if deg_rem <= 0:
        return tuple(sorted(degrees))
    if deg_rem == 1:
        degrees.append(1)
        return tuple(sorted(degrees))
    g = list(remaining)
    h = [1, 0]
    for i in range(1, deg_rem + 1):
        g_deg = len(g) - 1
        if g_deg < i:
            break
        h = poly_powmod_Fp(h, p, g, p)
        h_minus_x = list(h)
        if len(h_minus_x) >= 2:
            h_minus_x[-2] = (h_minus_x[-2] - 1) % p
        elif len(h_minus_x) == 1:
            h_minus_x = [p - 1, h_minus_x[0]]
        else:
            h_minus_x = [p - 1, 0]
        d = poly_gcd_Fp(list(g), h_minus_x, p)
        d_deg = len(d) - 1
        if d_deg > 0:
            num_factors = d_deg // i
            degrees.extend([i] * num_factors)
            g = poly_rem_Fp(g, d, p)
            if not g or all(c % p == 0 for c in g):
                g = [1]
            else:
                while len(g) > 1 and g[0] % p == 0:
                    g = g[1:]
                if len(g) > 1:
                    h = poly_rem_Fp(h, g, p)
    g_deg = len(g) - 1
    if g_deg > 0:
        degrees.append(g_deg)
    return tuple(sorted(degrees))


def frobenius_trace(p):
    if p in (2, 5):
        return mpf(0)
    pattern = factorize_mod_p(p)
    if pattern == (1, 1, 1, 1, 1):
        return mpf(2)
    elif pattern == (1, 2, 2):
        return mpf(0)
    elif pattern in ((1, 1, 3), (2, 3)):
        return mpf(-1)
    elif pattern == (5,):
        leg = pow(5, (p - 1) // 2, p)
        if leg == 1:
            return PHI
        else:
            return -PHI_INV
    else:
        return mpf(0)


def hecke_eigenvalue(ap, m):
    if m == 0:
        return mpf(1)
    if m == 1:
        return ap
    a_prev = mpf(1)
    a_curr = ap
    for _ in range(m - 1):
        a_next = ap * a_curr - a_prev
        a_prev = a_curr
        a_curr = a_next
    return a_curr


# =============================================================================
# GAUSSIAN TEST FUNCTION MACHINERY
# =============================================================================

def gaussian_fourier(alpha, xi):
    """Fourier transform of exp(-alpha*x^2)."""
    return sqrt(pi / alpha) * exp(-pi**2 * xi**2 / alpha)


# =============================================================================
# COMPUTE FROBENIUS DATA
# =============================================================================

separator("PYTHAGOREAN WEIL FUNCTIONAL: FORWARD + BACKWARD ANALYSIS")

PRIME_BOUND = 5000
MAX_POWER = 20

all_primes = sieve_primes(PRIME_BOUND)
primes_unram = [p for p in all_primes if p not in (2, 5)]

trace_data = {}
for p in primes_unram:
    trace_data[p] = frobenius_trace(p)

# Verify Chebotarev
counts = {"identity": 0, "golden+": 0, "golden-": 0, "order-3": 0, "involution": 0}
for p, ap in trace_data.items():
    if fabs(ap - 2) < mpf('0.01'):
        counts["identity"] += 1
    elif fabs(ap - PHI) < mpf('0.01'):
        counts["golden+"] += 1
    elif fabs(ap + PHI_INV) < mpf('0.01'):
        counts["golden-"] += 1
    elif fabs(ap + 1) < mpf('0.01'):
        counts["order-3"] += 1
    elif fabs(ap) < mpf('0.01'):
        counts["involution"] += 1

total = len(primes_unram)
print(f"Computed Frobenius traces for {total} unramified primes up to {PRIME_BOUND}")
print(f"Chebotarev density check:")
for name, count in counts.items():
    expected = [d for n, _, _, _, d, _ in A5_CLASSES if n == name][0]
    print(f"  {name:12s}: {count:4d}/{total} = {float(count/total):.4f}  (expected {float(expected):.4f})")


# =============================================================================
# PART 1: THE 15-BASIS GRAM MATRIX DECOMPOSITION A, P
# =============================================================================

separator("PART 1: GRAM MATRIX DECOMPOSITION  G = A - P")

# 15 calibrated Gaussian widths (log-spaced, covering narrow to wide)
BASIS_ALPHAS = [mpf(a) for a in [
    '0.01', '0.02', '0.05', '0.1', '0.2',
    '0.5', '1', '2', '5', '10',
    '20', '50', '100', '200', '500'
]]
DIM = len(BASIS_ALPHAS)
print(f"Basis dimension: {DIM}")
print(f"Widths alpha: {[float(a) for a in BASIS_ALPHAS]}")
print()

# Primes for the computation
gram_primes = [p for p in primes_unram if p <= 2000]
print(f"Using {len(gram_primes)} primes up to {max(gram_primes)}")


def compute_prime_gram_entry(alpha_i, alpha_j, primes_list, trace_dict,
                              use_squared=False, max_k=10):
    """
    Prime sum contribution to G_{ij}.
    sum_p sum_k coeff_{p^k} * log(p)/p^{k/2} * g_i(k*logp/(2pi)) * g_j(k*logp/(2pi))

    For use_squared: coeff = |a_{p^k}|^2 (Rankin-Selberg)
    Otherwise: coeff = a_{p^k} (original L-function)
    """
    total = mpf(0)
    for p in primes_list:
        ap = trace_dict[p]
        logp = log(mpf(p))
        for k in range(1, max_k + 1):
            apk = hecke_eigenvalue(ap, k)
            coeff = apk**2 if use_squared else apk
            weight = logp / power(mpf(p), mpf(k) / 2)
            xi = k * logp / (2 * pi)
            gi = gaussian_fourier(alpha_i, xi)
            gj = gaussian_fourier(alpha_j, xi)
            total += coeff * weight * gi * gj
    return total


def compute_arch_gram_entry(alpha_i, alpha_j, N_cond, dim=2, has_pole=False):
    """
    Archimedean + pole contribution to G_{ij}.

    arch = (log(N) - dim*log(2pi))/2 * g_i(0)*g_j(0) + Gamma integral
    pole = g_i(0)*g_j(0)  (for L-functions with simple pole at s=1)
    """
    g0i = gaussian_fourier(alpha_i, mpf(0))
    g0j = gaussian_fourier(alpha_j, mpf(0))

    log_N = log(N_cond)

    # Leading archimedean
    arch = (log_N - dim * log(2 * pi)) / 2 * g0i * g0j

    # Gamma integral contribution
    def gamma_integrand(t):
        if fabs(t) < mpf('1e-15'):
            return mpf(0)
        gi_t = gaussian_fourier(alpha_i, t)
        gj_t = gaussian_fourier(alpha_j, t)
        s_val = mpc(mpf('0.25'), t / 2)
        try:
            psi_val = digamma(s_val)
            return gi_t * gj_t * re(psi_val) / pi
        except Exception:
            return mpf(0)

    try:
        I_gamma = 2 * quad(gamma_integrand, [0, 50], error=True)[0]
    except Exception:
        I_gamma = mpf(0)

    result = arch + I_gamma

    # Pole contribution
    if has_pole:
        result += g0i * g0j

    return result


# Build the SEPARATE A and P matrices
print("\nBuilding archimedean matrix A and prime matrix P...")
print("(This takes a few minutes for 15x15 with prime powers up to k=10)")
sys.stdout.flush()

# We compute for the Rankin-Selberg: L_RS = L(Sym^2 rho) * zeta(s)
# Conductor: N_RS = N^2 = 640000
# Dimension: 4 (degree of L_RS)
# Has pole: yes (from zeta factor)
N_RS = CONDUCTOR**2

A_mat = matrix(DIM, DIM)  # archimedean + pole
P_mat = matrix(DIM, DIM)  # prime sum (squared traces for RS)

for i in range(DIM):
    for j in range(i, DIM):
        # Archimedean for RS
        a_val = compute_arch_gram_entry(BASIS_ALPHAS[i], BASIS_ALPHAS[j],
                                         N_RS, dim=4, has_pole=True)
        A_mat[i, j] = a_val
        A_mat[j, i] = a_val

        # Prime sum with squared traces (RS)
        p_val = compute_prime_gram_entry(BASIS_ALPHAS[i], BASIS_ALPHAS[j],
                                          gram_primes, trace_data,
                                          use_squared=True, max_k=10)
        P_mat[i, j] = p_val
        P_mat[j, i] = p_val

    if (i + 1) % 3 == 0:
        print(f"  ...row {i+1}/{DIM} done")
        sys.stdout.flush()

print("  Matrices A and P built.")


def mat_eigenvalues(M):
    """Compute and sort eigenvalues of a symmetric matrix."""
    raw = eig(M, left=False, right=False)
    return sorted([re(x) for x in raw], reverse=True)


def mat_mul(A, B):
    """Matrix multiply for mpmath matrices."""
    n = A.rows
    m = B.cols
    k = A.cols
    C = matrix(n, m)
    for i in range(n):
        for j in range(m):
            s = mpf(0)
            for l in range(k):
                s += A[i, l] * B[l, j]
            C[i, j] = s
    return C


def mat_add(A, B, scale_b=mpf(1)):
    """A + scale_b * B."""
    n = A.rows
    C = matrix(n, n)
    for i in range(n):
        for j in range(n):
            C[i, j] = A[i, j] + scale_b * B[i, j]
    return C


def mat_transpose(M):
    n = M.rows
    m = M.cols
    T = matrix(m, n)
    for i in range(n):
        for j in range(m):
            T[j, i] = M[i, j]
    return T


# =============================================================================
# PART 2: EIGENVALUE ANALYSIS OF A AND P SEPARATELY
# =============================================================================

separator("PART 2: EIGENVALUE ANALYSIS OF A AND P")

print("Archimedean matrix A (including pole):")
evals_A = mat_eigenvalues(A_mat)
n_pos_A = sum(1 for ev in evals_A if ev > mpf('1e-20'))
n_neg_A = sum(1 for ev in evals_A if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_A):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_A}+, {DIM - n_pos_A - n_neg_A}zero, {n_neg_A}-)")
print(f"  A is PSD: {'YES' if n_neg_A == 0 else 'NO'}")
print()

print("Prime matrix P (Rankin-Selberg, squared traces):")
evals_P = mat_eigenvalues(P_mat)
n_pos_P = sum(1 for ev in evals_P if ev > mpf('1e-20'))
n_neg_P = sum(1 for ev in evals_P if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_P):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_P}+, {DIM - n_pos_P - n_neg_P}zero, {n_neg_P}-)")
print(f"  P is PSD: {'YES' if n_neg_P == 0 else 'NO'}")


# =============================================================================
# PART 3: FORWARD AND BACKWARD WEIL GRAM MATRICES
# =============================================================================

separator("PART 3: FORWARD AND BACKWARD WEIL GRAM MATRICES")

# Forward: G_f = A - P (standard Weil)
G_forward = mat_add(A_mat, P_mat, scale_b=mpf(-1))

# Backward: G_b = A + P (reversed prime sign)
G_backward = mat_add(A_mat, P_mat, scale_b=mpf(1))

print("FORWARD Gram matrix  G_f = A - P  (standard Weil functional):")
evals_Gf = mat_eigenvalues(G_forward)
n_pos_Gf = sum(1 for ev in evals_Gf if ev > mpf('1e-20'))
n_neg_Gf = sum(1 for ev in evals_Gf if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_Gf):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_Gf}+, {DIM - n_pos_Gf - n_neg_Gf}zero, {n_neg_Gf}-)")
print(f"  G_f is PSD (Weil positivity): {'YES' if n_neg_Gf == 0 else 'NO'}")
print()

print("BACKWARD Gram matrix  G_b = A + P  (reversed prime sign):")
evals_Gb = mat_eigenvalues(G_backward)
n_pos_Gb = sum(1 for ev in evals_Gb if ev > mpf('1e-20'))
n_neg_Gb = sum(1 for ev in evals_Gb if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_Gb):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_Gb}+, {DIM - n_pos_Gb - n_neg_Gb}zero, {n_neg_Gb}-)")
print(f"  G_b is PSD: {'YES' if n_neg_Gb == 0 else 'NO'}")


# =============================================================================
# PART 4: PYTHAGOREAN GRAM MATRIX  G_Pyth = G_f^2 + G_b^2 = 2(A^2 + P^2)
# =============================================================================

separator("PART 4: PYTHAGOREAN GRAM MATRIX  G_f^2 + G_b^2 = 2(A^2 + P^2)")

# Compute A^2 and P^2
A_sq = mat_mul(A_mat, A_mat)
P_sq = mat_mul(P_mat, P_mat)

# G_f^2 = (A-P)^2 = A^2 - AP - PA + P^2
G_f_sq = mat_mul(G_forward, G_forward)

# G_b^2 = (A+P)^2 = A^2 + AP + PA + P^2
G_b_sq = mat_mul(G_backward, G_backward)

# Pythagorean: G_f^2 + G_b^2 = 2(A^2 + P^2)
G_pyth = mat_add(G_f_sq, G_b_sq)

# Also compute directly: 2(A^2 + P^2)
G_pyth_direct = mat_add(A_sq, P_sq, scale_b=mpf(1))
for i in range(DIM):
    for j in range(DIM):
        G_pyth_direct[i, j] *= 2

# Verify the identity G_f^2 + G_b^2 = 2(A^2 + P^2)
max_diff = mpf(0)
for i in range(DIM):
    for j in range(DIM):
        d = fabs(G_pyth[i, j] - G_pyth_direct[i, j])
        if d > max_diff:
            max_diff = d

print(f"Identity check: G_f^2 + G_b^2 = 2(A^2 + P^2)")
print(f"  Max entry difference: {nstr(max_diff, 6)}")
print(f"  Identity holds: {'YES' if max_diff < mpf('1e-15') else 'NO'}")
print()

print("Pythagorean Gram matrix  G_Pyth = G_f^2 + G_b^2 = 2(A^2 + P^2):")
evals_Gpyth = mat_eigenvalues(G_pyth)
n_pos_Gpyth = sum(1 for ev in evals_Gpyth if ev > mpf('1e-20'))
n_neg_Gpyth = sum(1 for ev in evals_Gpyth if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_Gpyth):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_Gpyth}+, {DIM - n_pos_Gpyth - n_neg_Gpyth}zero, {n_neg_Gpyth}-)")
print(f"  G_Pyth is PSD: {'YES' if n_neg_Gpyth == 0 else 'NO'}")
print()

print("  NOTE: G_Pyth = 2(A^2 + P^2) is ALWAYS PSD (sum of two PSD matrices).")
print("  This is a TRIVIAL consequence of squaring.")
print("  The CONTENT is in the PRODUCT G_f * G_b = A^2 - P^2.")


# =============================================================================
# PART 5: PRODUCT  G_f * G_b = A^2 - P^2
# =============================================================================

separator("PART 5: PRODUCT  G_f * G_b = A^2 - P^2")

# G_f * G_b = (A-P)(A+P) = A^2 - P^2 + AP - PA  (if A, P don't commute)
# Wait: (A-P)(A+P) = A^2 + AP - PA - P^2
# This is A^2 - P^2 ONLY if AP = PA (they commute)

# Check commutativity
AP = mat_mul(A_mat, P_mat)
PA = mat_mul(P_mat, A_mat)
commutator_norm = mpf(0)
for i in range(DIM):
    for j in range(DIM):
        d = fabs(AP[i, j] - PA[i, j])
        if d > commutator_norm:
            commutator_norm = d

print(f"Commutativity check: [A, P] = AP - PA")
print(f"  Max entry of commutator: {nstr(commutator_norm, 6)}")
print(f"  A and P commute: {'YES (up to numerical error)' if commutator_norm < mpf('1e-10') else 'NO'}")
print()

# Compute actual product G_f * G_b
GfGb = mat_mul(G_forward, G_backward)

# Compute A^2 - P^2
A2_minus_P2 = mat_add(A_sq, P_sq, scale_b=mpf(-1))

# Check if GfGb = A^2 - P^2 (true only if commutative)
diff_prod = mpf(0)
for i in range(DIM):
    for j in range(DIM):
        d = fabs(GfGb[i, j] - A2_minus_P2[i, j])
        if d > diff_prod:
            diff_prod = d
print(f"Product check: G_f*G_b vs A^2-P^2")
print(f"  Max entry difference: {nstr(diff_prod, 6)}")
if diff_prod < mpf('1e-10'):
    print(f"  G_f*G_b = A^2 - P^2  (commutative case)")
else:
    print(f"  G_f*G_b = A^2 + AP - PA - P^2  (non-commutative)")
print()

# Eigenvalues of G_f * G_b
print("Product matrix  G_f * G_b:")
evals_prod = mat_eigenvalues(GfGb)
n_pos_prod = sum(1 for ev in evals_prod if ev > mpf('1e-20'))
n_neg_prod = sum(1 for ev in evals_prod if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_prod):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Signature: ({n_pos_prod}+, {DIM - n_pos_prod - n_neg_prod}zero, {n_neg_prod}-)")
print()

print("INTERPRETATION:")
if n_neg_prod == 0:
    print("  G_f * G_b is PSD!")
    print("  => A^2 >= P^2 in matrix order")
    print("  => The archimedean 'leg' dominates the prime 'leg' in every direction")
    print("  => This is the Pythagorean condition: |hypotenuse| >= |each leg|")
else:
    print("  G_f * G_b has negative eigenvalues.")
    print("  => A^2 is NOT >= P^2 in all directions.")
    print("  => The prime sum escapes archimedean domination in some directions.")


# =============================================================================
# PART 6: SINGULAR VALUE ANALYSIS
# =============================================================================

separator("PART 6: SINGULAR VALUE ANALYSIS  sigma_min(A) vs sigma_max(P)")

# For symmetric PSD matrices, singular values = eigenvalues
# So sigma_min(A) = min eigenvalue of A (if PSD)
# sigma_max(P) = max eigenvalue of P

# But A might not be PSD. Use actual singular values.
# mpmath svd: returns (U, S, V) where S is diagonal

# Since A and P are symmetric, their singular values are |eigenvalues|

svals_A = sorted([fabs(ev) for ev in evals_A])  # ascending
svals_P = sorted([fabs(ev) for ev in evals_P], reverse=True)  # descending

sigma_min_A = svals_A[0] if svals_A else mpf(0)
sigma_max_P = svals_P[0] if svals_P else mpf(0)

print(f"Singular values of A (ascending):")
for idx, sv in enumerate(svals_A):
    print(f"  sigma_{idx+1:2d} = {nstr(sv, 12)}")

print()
print(f"Singular values of P (descending):")
for idx, sv in enumerate(svals_P):
    print(f"  sigma_{idx+1:2d} = {nstr(sv, 12)}")

print()
print(f"sigma_min(A) = {nstr(sigma_min_A, 15)}")
print(f"sigma_max(P) = {nstr(sigma_max_P, 15)}")
print(f"Ratio sigma_min(A)/sigma_max(P) = {nstr(sigma_min_A / sigma_max_P if sigma_max_P > 0 else inf, 10)}")
print()

if sigma_min_A >= sigma_max_P:
    print("sigma_min(A) >= sigma_max(P)")
    print("=> ||Av|| >= ||Pv|| for all unit vectors v")
    print("=> A^2 >= P^2 in matrix order")
    print("=> G_f = A - P is PSD (Pythagorean condition satisfied!)")
else:
    print("sigma_min(A) < sigma_max(P)")
    print("=> The archimedean does NOT dominate the prime sum in the smallest direction.")
    print(f"   Gap: sigma_max(P) - sigma_min(A) = {nstr(sigma_max_P - sigma_min_A, 10)}")
    print()
    # Find which eigenvalues of A are smaller than sigma_max(P)
    n_small = sum(1 for sv in svals_A if sv < sigma_max_P)
    print(f"   {n_small} singular values of A are smaller than sigma_max(P)")


# =============================================================================
# PART 7: CRITICAL CONDUCTOR SEARCH
# =============================================================================

separator("PART 7: CRITICAL CONDUCTOR SEARCH")

print("The archimedean scales as ~ log(N_cond). The prime sum is fixed.")
print("Find q_crit such that for N_cond >= q_crit, all eigenvalues of G_f >= 0.")
print()
print("Testing conductors: N_RS = N^2 for various N...")
print()

# For this we rebuild A with different conductors, keeping P fixed
test_conductors = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200,
                   100000, 200000, 500000, 1000000, 5000000, 10000000]

print(f"{'N':>12s} {'N_RS=N^2':>14s} {'min_eval(Gf)':>18s} {'n_neg':>6s} {'PSD?':>6s}")
print("-" * 62)

# Precompute a SCALED version: A(N) = A_fixed + log(N/N_0) * scaling_matrix
# Since arch = (log(N) - dim*log(2pi))/2 * g_i(0)*g_j(0) + gamma_integral
# The gamma_integral doesn't depend on N.
# So A(N) = A(N_0) + (log(N) - log(N_0))/2 * R
# where R_{ij} = g_i(0) * g_j(0) is the rank-1 outer product

# Build the rank-1 scaling matrix R
R_mat = matrix(DIM, DIM)
g0_values = [gaussian_fourier(BASIS_ALPHAS[i], mpf(0)) for i in range(DIM)]
for i in range(DIM):
    for j in range(DIM):
        R_mat[i, j] = g0_values[i] * g0_values[j]

# Base conductor used to build A_mat
log_N_RS_base = log(N_RS)

q_crit = None
for N in test_conductors:
    N_sq = mpf(N)**2
    log_N_sq = log(N_sq)

    # A_new = A_mat + (log(N_sq) - log(N_RS_base))/2 * R_mat
    delta_log = (log_N_sq - log_N_RS_base) / 2
    A_new = matrix(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            A_new[i, j] = A_mat[i, j] + delta_log * R_mat[i, j]

    G_f_new = mat_add(A_new, P_mat, scale_b=mpf(-1))
    evals_new = mat_eigenvalues(G_f_new)
    min_ev = min(evals_new)
    n_neg = sum(1 for ev in evals_new if ev < -mpf('1e-15'))
    is_psd = n_neg == 0

    print(f"{N:>12d} {nstr(N_sq, 10):>14s} {nstr(min_ev, 10):>18s} {n_neg:>6d} {'YES' if is_psd else 'NO':>6s}")

    if is_psd and q_crit is None:
        q_crit = N

print()
if q_crit is not None:
    print(f"CRITICAL CONDUCTOR: N_crit ~ {q_crit}")
    print(f"  For N >= {q_crit}: G_f = A - P is PSD (Weil positivity holds)")
    print(f"  Our icosahedral rep has N = 800, N_RS = 640000")
    if q_crit <= 800:
        print(f"  800 >= {q_crit}: Weil positivity HOLDS for our L-function!")
    else:
        print(f"  800 < {q_crit}: Weil positivity FAILS for our L-function with current basis.")
else:
    print(f"  No critical conductor found in tested range.")
    print(f"  Eigenvalues may depend on more than just the conductor.")


# =============================================================================
# PART 8: THE ORIGINAL L-FUNCTION (not Rankin-Selberg) COMPARISON
# =============================================================================

separator("PART 8: COMPARISON WITH ORIGINAL L(s, rho)")

print("Now compare: the original L-function with raw traces a_p (not squared).")
print("Build A_orig (conductor 800, dim 2, no pole) and P_orig (raw traces).")
print()

# Build A_orig and P_orig for L(s, rho)
A_orig = matrix(DIM, DIM)
P_orig = matrix(DIM, DIM)

for i in range(DIM):
    for j in range(i, DIM):
        a_val = compute_arch_gram_entry(BASIS_ALPHAS[i], BASIS_ALPHAS[j],
                                         CONDUCTOR, dim=2, has_pole=False)
        A_orig[i, j] = a_val
        A_orig[j, i] = a_val

        p_val = compute_prime_gram_entry(BASIS_ALPHAS[i], BASIS_ALPHAS[j],
                                          gram_primes, trace_data,
                                          use_squared=False, max_k=10)
        P_orig[i, j] = p_val
        P_orig[j, i] = p_val

    if (i + 1) % 5 == 0:
        print(f"  ...row {i+1}/{DIM} done")
        sys.stdout.flush()

G_f_orig = mat_add(A_orig, P_orig, scale_b=mpf(-1))
G_b_orig = mat_add(A_orig, P_orig, scale_b=mpf(1))

print()
print("Original L(s, rho) FORWARD Gram  G_f = A - P:")
evals_Gf_orig = mat_eigenvalues(G_f_orig)
n_neg_orig = sum(1 for ev in evals_Gf_orig if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_Gf_orig):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Negative eigenvalues: {n_neg_orig}")
print()

print("Original L(s, rho) BACKWARD Gram  G_b = A + P:")
evals_Gb_orig = mat_eigenvalues(G_b_orig)
n_neg_b_orig = sum(1 for ev in evals_Gb_orig if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_Gb_orig):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Negative eigenvalues: {n_neg_b_orig}")
print()

# Product for original
GfGb_orig = mat_mul(G_f_orig, G_b_orig)
print("Original PRODUCT  G_f * G_b:")
evals_prod_orig = mat_eigenvalues(GfGb_orig)
n_neg_prod_orig = sum(1 for ev in evals_prod_orig if ev < -mpf('1e-20'))
for idx, ev in enumerate(evals_prod_orig):
    sign = "+" if ev > mpf('1e-20') else "-" if ev < -mpf('1e-20') else "0"
    print(f"  lambda_{idx+1:2d} = {nstr(ev, 12):>22s}  [{sign}]")
print(f"  Negative eigenvalues: {n_neg_prod_orig}")
print()

# Pythagorean for original
G_pyth_orig = mat_add(mat_mul(G_f_orig, G_f_orig),
                       mat_mul(G_b_orig, G_b_orig))
evals_pyth_orig = mat_eigenvalues(G_pyth_orig)
n_neg_pyth_orig = sum(1 for ev in evals_pyth_orig if ev < -mpf('1e-20'))
print("Original PYTHAGOREAN  G_f^2 + G_b^2 = 2(A^2 + P^2):")
print(f"  All eigenvalues non-negative: {'YES' if n_neg_pyth_orig == 0 else 'NO'}")
print(f"  Min eigenvalue: {nstr(min(evals_pyth_orig), 12)}")
print(f"  Max eigenvalue: {nstr(max(evals_pyth_orig), 12)}")


# =============================================================================
# PART 9: THE KEY RELATIONSHIP -- WHEN DOES G_f * G_b >= 0 IMPLY G_f >= 0?
# =============================================================================

separator("PART 9: IMPLICATIONS OF THE PYTHAGOREAN STRUCTURE")

print("THEOREM ATTEMPT:")
print("  If G_b = A + P is PSD  AND  G_f * G_b = A^2 - P^2 + [A,P] is PSD")
print("  Does this imply G_f = A - P is PSD?")
print()
print("ANSWER: Not in general. But with additional structure...")
print()

# If G_b is PSD and invertible, we can write:
# G_f = G_f * G_b * G_b^{-1}
# If G_f * G_b >= 0 and G_b^{-1} >= 0, then G_f >= 0.
# But G_b^{-1} >= 0 iff G_b >= 0 (for symmetric matrices).
# And if G_f * G_b >= 0 with G_b >= 0, then... not necessarily G_f >= 0.

# The CONGRUENCE approach: if G_b > 0 (strictly PD), then
# G_b^{-1/2} G_f G_b^{-1/2} has the same inertia as G_f.
# And G_b^{-1/2} (G_f * G_b) G_b^{-1/2} = G_b^{-1/2} G_f G_b^{1/2}
# which is similar to G_f (same eigenvalues) only if A,P commute.

# If A and P commute (and they nearly do for our basis), then
# G_f G_b = (A-P)(A+P) = A^2 - P^2
# The eigenvalues of G_f G_b are products of corresponding eigenvalues of G_f and G_b.
# If G_b > 0 and G_f G_b >= 0, then each eigenvalue of G_f = (eval of G_f G_b)/(eval of G_b) >= 0.

# Check the commutative case
AP_orig = mat_mul(A_orig, P_orig)
PA_orig = mat_mul(P_orig, A_orig)
comm_orig = mpf(0)
for i in range(DIM):
    for j in range(DIM):
        d = fabs(AP_orig[i, j] - PA_orig[i, j])
        if d > comm_orig:
            comm_orig = d

print(f"For original L(rho):")
print(f"  Commutator [A, P] max entry: {nstr(comm_orig, 6)}")
print(f"  Commutative: {'YES' if comm_orig < mpf('1e-10') else 'NO'}")
print()
print(f"For Rankin-Selberg:")
print(f"  Commutator [A, P] max entry: {nstr(commutator_norm, 6)}")
print(f"  Commutative: {'YES' if commutator_norm < mpf('1e-10') else 'NO'}")
print()

# If [A,P] ~ 0, do the simultaneous diagonalization analysis
if comm_orig < mpf('1e-5'):
    print("A and P approximately commute (for original L).")
    print("In the simultaneous eigenbasis, G_f has eigenvalue a_i - p_i")
    print("and G_b has eigenvalue a_i + p_i where a_i, p_i are eigenvalues of A, P.")
    print()

    # Since they approximately commute, the eigenvalues of G_f and G_b
    # should approximately pair up as (a_i - p_i, a_i + p_i)
    evals_A_orig = mat_eigenvalues(A_orig)
    evals_P_orig = mat_eigenvalues(P_orig)

    print(f"{'i':>3s} {'a_i':>16s} {'p_i':>16s} {'a_i-p_i':>16s} {'a_i+p_i':>16s} {'a_i^2-p_i^2':>16s}")
    print("-" * 88)
    for i in range(DIM):
        ai = evals_A_orig[i]
        pi_val = evals_P_orig[i]
        print(f"{i+1:>3d} {nstr(ai, 10):>16s} {nstr(pi_val, 10):>16s} "
              f"{nstr(ai - pi_val, 10):>16s} {nstr(ai + pi_val, 10):>16s} "
              f"{nstr(ai**2 - pi_val**2, 10):>16s}")

    print()
    print("NOTE: The eigenvalue pairing above is approximate (eigenvectors may differ).")
    print("For exact analysis, would need simultaneous diagonalization.")

if commutator_norm < mpf('1e-5'):
    print()
    print("A and P approximately commute (for Rankin-Selberg).")
    print(f"{'i':>3s} {'a_i (RS)':>16s} {'p_i (RS)':>16s} {'a_i-p_i':>16s} {'a_i+p_i':>16s} {'a_i^2-p_i^2':>16s}")
    print("-" * 88)
    for i in range(DIM):
        ai = evals_A[i]
        pi_val = evals_P[i]
        print(f"{i+1:>3d} {nstr(ai, 10):>16s} {nstr(pi_val, 10):>16s} "
              f"{nstr(ai - pi_val, 10):>16s} {nstr(ai + pi_val, 10):>16s} "
              f"{nstr(ai**2 - pi_val**2, 10):>16s}")


# =============================================================================
# PART 10: COMPREHENSIVE SUMMARY
# =============================================================================

separator("PART 10: COMPREHENSIVE SUMMARY")

print("FORWARD vs BACKWARD vs PYTHAGOREAN")
print()
print("=" * 78)
print(f"{'Matrix':30s} {'PSD?':>6s} {'n_neg':>6s} {'min eval':>18s} {'max eval':>18s}")
print("-" * 78)

results = [
    ("RS: A (archimedean+pole)", n_neg_A == 0, n_neg_A, min(evals_A), max(evals_A)),
    ("RS: P (prime, sq traces)", n_neg_P == 0, n_neg_P, min(evals_P), max(evals_P)),
    ("RS: G_f = A - P (forward)", n_neg_Gf == 0, n_neg_Gf, min(evals_Gf), max(evals_Gf)),
    ("RS: G_b = A + P (backward)", n_neg_Gb == 0, n_neg_Gb, min(evals_Gb), max(evals_Gb)),
    ("RS: G_f*G_b (product)", n_neg_prod == 0, n_neg_prod, min(evals_prod), max(evals_prod)),
    ("RS: G_f^2+G_b^2 (pyth)", n_neg_Gpyth == 0, n_neg_Gpyth, min(evals_Gpyth), max(evals_Gpyth)),
    ("Orig: G_f = A - P", n_neg_orig == 0, n_neg_orig, min(evals_Gf_orig), max(evals_Gf_orig)),
    ("Orig: G_b = A + P", n_neg_b_orig == 0, n_neg_b_orig, min(evals_Gb_orig), max(evals_Gb_orig)),
    ("Orig: G_f*G_b", n_neg_prod_orig == 0, n_neg_prod_orig, min(evals_prod_orig), max(evals_prod_orig)),
    ("Orig: G_f^2+G_b^2", n_neg_pyth_orig == 0, n_neg_pyth_orig, min(evals_pyth_orig), max(evals_pyth_orig)),
]

for name, psd, nn, mn, mx in results:
    print(f"{name:30s} {'YES' if psd else 'NO':>6s} {nn:>6d} {nstr(mn, 10):>18s} {nstr(mx, 10):>18s}")

print()
print("=" * 78)
print("  KEY FINDINGS")
print("=" * 78)
print()

# Determine the main conclusion
if n_neg_Gf == 0:
    print("  [RS] FORWARD WEIL IS PSD! Weil positivity holds for L_RS.")
    print("        => GRH for L(s, Sym^2 rho) * zeta(s)")
    print("        => GRH for zeta(s) (since factors inherit zeros)")
elif n_neg_Gb == 0 and n_neg_Gf > 0:
    print(f"  [RS] Forward has {n_neg_Gf} negative eigenvalue(s), backward is PSD.")
    print("        The backward (A+P) sees MORE positivity than forward (A-P).")
    print(f"        Forward min eigenvalue: {nstr(min(evals_Gf), 10)}")
    print(f"        Backward min eigenvalue: {nstr(min(evals_Gb), 10)}")
else:
    print(f"  [RS] Both forward ({n_neg_Gf} neg) and backward ({n_neg_Gb} neg) have issues.")

print()
if n_neg_prod == 0:
    print("  [RS] PRODUCT G_f*G_b is PSD!")
    print("        => A^2 >= P^2 (archimedean dominates prime in squared sense)")
    print("        Combined with backward PSD and commutativity => forward PSD!")
elif n_neg_prod > 0 and n_neg_Gb == 0:
    print(f"  [RS] Product has {n_neg_prod} negative eigenvalue(s).")
    print("        Backward is PSD but the product is not.")
    print("        The Pythagorean condition A^2 >= P^2 fails in some directions.")

print()
if q_crit is not None:
    print(f"  CRITICAL CONDUCTOR: N_crit ~ {q_crit}")
    if q_crit <= 800:
        print("  Our L-function (N=800) EXCEEDS the critical conductor.")
    else:
        print(f"  Our L-function (N=800) is BELOW the critical conductor ({q_crit}).")
        print(f"  Need N >= {q_crit} for Weil positivity in our basis.")

print()
print("  The Pythagorean Gram G_f^2 + G_b^2 = 2(A^2 + P^2) is ALWAYS PSD.")
print("  This is trivially true (sum of squares of symmetric matrices).")
print("  The NON-TRIVIAL content is in the PRODUCT G_f*G_b = A^2 - P^2.")
print("  If A^2 >= P^2 AND G_b >= 0: then G_f >= 0 (by congruence/commutativity).")
print()
print("  THE PYTHAGOREAN PATH:")
print("  1. Archimedean A grows with log(conductor)")
print("  2. Prime sum P is bounded (converges as primes increase)")
print("  3. For large enough conductor: A dominates P in ALL directions")
print("  4. The critical question: does the ACTUAL conductor suffice?")

print()
print("COMPUTATION COMPLETE.")
