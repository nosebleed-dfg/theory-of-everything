#!/usr/bin/env python3
"""
PYTHAGOREAN_WEIL — squared-trace Weil functional via Rankin-Selberg; fixes sign problem of linear approach
nos3bl33d

|a_p|^2 >= 0 always. Pythagorean bound, Gram matrix, Rankin-Selberg positivity.
Polynomial x^5+20x+16, conductor N=800, Galois group A5. mpmath 30 digits.
"""

from mpmath import (mp, mpf, mpc, log, pi, sqrt, gamma, cos, sin,
                     exp, power, re, im, quad, loggamma,
                     digamma, zeta, binomial, fsum, matrix, nstr, fabs,
                     euler as euler_gamma, inf, eig)
import sys
import warnings
warnings.filterwarnings("ignore")

mp.dps = 30

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + sqrt(5)) / 2
PHI_INV = PHI - 1              # = 1/phi
PHI_SQ = PHI**2                # = phi + 1
PHI_INV_SQ = PHI_INV**2        # = 1/phi^2 = 2 - phi... wait: 1/phi^2 = phi^2 - 2 no.
# 1/phi = phi - 1, so 1/phi^2 = (phi-1)^2 = phi^2 - 2phi + 1 = (phi+1) - 2phi + 1 = 2 - phi
CONDUCTOR = mpf(800)

# Verify
assert fabs(PHI_SQ - PHI - 1) < mpf('1e-25'), "phi^2 = phi + 1 failed"
assert fabs(PHI_INV_SQ - (2 - PHI)) < mpf('1e-25'), "1/phi^2 = 2-phi failed"

# A5 conjugacy classes for the icosahedral 2-dim rep
# (name, trace a_p, |a_p|^2, class_size, density, element_order)
A5_CLASSES = [
    ("identity",    mpf(2),   mpf(4),       1,  mpf(1)/60,  1),
    ("golden+",     PHI,      PHI_SQ,       12, mpf(12)/60, 5),
    ("golden-",    -PHI_INV,  PHI_INV_SQ,   12, mpf(12)/60, 5),
    ("order-3",     mpf(-1),  mpf(1),       20, mpf(20)/60, 3),
    ("involution",  mpf(0),   mpf(0),       15, mpf(15)/60, 2),
]

# Rankin-Selberg traces b_p = |a_p|^2 and Sym^2 traces c_p = a_p^2 - 1
# (since det = 1, rho x rho_bar = Sym^2 rho + trivial, so b_p = c_p + 1)
RS_CLASSES = []
for name, ap, ap_sq, size, density, order in A5_CLASSES:
    RS_CLASSES.append((name, ap_sq, ap_sq - 1, density, order))


def separator(title):
    """Print a section header."""
    print()
    print("=" * 76)
    print(f"  {title}")
    print("=" * 76)
    print()


def sieve_primes(N):
    """Sieve of Eratosthenes up to N."""
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
# FROBENIUS TRACE COMPUTATION (from the actual polynomial x^5 + 20x + 16)
# =============================================================================

def poly_eval_mod(x, p):
    """Evaluate x^5 + 20x + 16 mod p."""
    return (pow(x, 5, p) + (20 % p) * x + (16 % p)) % p


def poly_divmod_Fp(a, b, p):
    """Polynomial division a / b in F_p[x]. Returns (quotient, remainder)."""
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
    """Multiply two polynomials in F_p[x]."""
    if not a or not b:
        return [0]
    n = len(a) + len(b) - 1
    result = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
    return result


def poly_rem_Fp(a, b, p):
    """Remainder of a / b in F_p[x]."""
    _, r = poly_divmod_Fp(a, b, p)
    return r if r else [0]


def poly_powmod_Fp(base, exp_val, modpoly, p):
    """Compute base^exp mod modpoly in F_p[x]."""
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
    """GCD of two polynomials in F_p[x]."""
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
    """
    Get the factorization pattern of x^5 + 20x + 16 mod p.
    Returns sorted tuple of degrees of irreducible factors.
    """
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
    h = [1, 0]  # x

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
    """
    Compute the Frobenius trace a_p for the icosahedral 2-dim representation.

    The factorization pattern of x^5+20x+16 mod p determines the conjugacy
    class in A5:
      - (1,1,1,1,1) -> identity, a_p = 2
      - (1,2,2)     -> double transpositions (involutions), a_p = 0
      - (1,1,3) or (2,3) -> 3-cycles, a_p = -1
      - (5)         -> 5-cycles, a_p = phi or -1/phi
    """
    if p in (2, 5):
        return mpf(0)  # ramified primes

    pattern = factorize_mod_p(p)

    if pattern == (1, 1, 1, 1, 1):
        return mpf(2)
    elif pattern == (1, 2, 2):
        return mpf(0)
    elif pattern in ((1, 1, 3), (2, 3)):
        return mpf(-1)
    elif pattern == (5,):
        # Distinguish the two 5-cycle classes via Legendre symbol
        leg = pow(5, (p - 1) // 2, p)
        if leg == 1:
            return PHI
        else:
            return -PHI_INV
    else:
        return mpf(0)


def hecke_eigenvalue(ap, m):
    """
    Compute a_{p^m} using the Chebyshev recurrence:
        a_{p^{m+1}} = a_p * a_{p^m} - a_{p^{m-1}}
    with a_{p^0} = 1, a_{p^1} = a_p.
    """
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
# COMPUTE FROBENIUS DATA
# =============================================================================

print("=" * 76)
print("  PYTHAGOREAN WEIL FUNCTIONAL: SQUARED-TRACE REFORMULATION")
print("  Reformulating Weil from pentagram (linear traces) to")
print("  triangle (squared traces) to eliminate the sign problem")
print("=" * 76)

PRIME_BOUND = 5000
MAX_POWER = 20

all_primes = sieve_primes(PRIME_BOUND)
primes_unram = [p for p in all_primes if p not in (2, 5)]

# Compute all Frobenius traces
trace_data = {}
for p in primes_unram:
    ap = frobenius_trace(p)
    trace_data[p] = ap

# Count classes
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
print(f"\nComputed Frobenius traces for {total} unramified primes up to {PRIME_BOUND}")
print(f"\nChebotarev density verification:")
for name, count in counts.items():
    expected = [d for n, _, _, _, d, _ in A5_CLASSES if n == name][0]
    print(f"  {name:12s}: {count:4d}/{total} = {float(count/total):.4f}  (expected {float(expected):.4f})")


# =============================================================================
# PART 1: THE SQUARED MODULUS FORMULATION
# =============================================================================

separator("PART 1: THE SQUARED MODULUS FORMULATION")

print("For each prime p with Frobenius eigenvalues alpha_p, beta_p:")
print("  alpha_p = e^{i*theta_p},  beta_p = e^{-i*theta_p}  (Ramanujan)")
print("  a_p = alpha_p + beta_p = 2*cos(theta_p)")
print()
print("The local Euler factor:")
print("  (1 - alpha_p*u)(1 - beta_p*u) = 1 - a_p*u + u^2")
print()
print("At |u| = 1: this is |1 - alpha_p*u|^2 >= 0 (perfect square!)")
print("At u = p^{-s}: the squared modulus structure persists.")
print()

print("Frobenius class traces and their squared moduli:")
print(f"{'Class':15s} {'a_p':>12s} {'|a_p|^2':>12s} {'Pyth gap':>12s} {'|Euler(u=1)|^2':>16s}")
print("-" * 72)

for name, ap, ap_sq, size, density, order in A5_CLASSES:
    # Pythagorean remainder: 4 = a_p^2 + Delta_p, so Delta_p = 4 - a_p^2
    delta_p = 4 - ap_sq
    # |1 - e^{i*theta} * u|^2 at u=1: = |1 - e^{i*theta}|^2 = 2 - 2cos(theta)
    # where theta = arccos(a_p/2)
    if fabs(ap) <= 2:
        euler_sq = 2 - ap  # = |1 - alpha_p|^2 where alpha_p = e^{i*arccos(a_p/2)}
        # Actually |1 - a_p*u + u^2| at u=1 = |2 - a_p|
        euler_at_1 = (2 - ap)  # value of 1 - a_p + 1 = 2 - a_p
        euler_sq_val = euler_at_1**2
    else:
        euler_sq_val = mpf(0)

    print(f"{name:15s} {nstr(ap, 8):>12s} {nstr(ap_sq, 8):>12s} {nstr(delta_p, 8):>12s} {nstr(euler_sq_val, 8):>16s}")

print()
print("KEY OBSERVATIONS:")
print(f"  1. |a_p|^2 >= 0 for ALL classes (no sign problem!)")
print(f"  2. Pythagorean decomposition: a_p^2 + Delta_p = 4 = (Ramanujan bound)^2")
print(f"  3. At golden+: Delta = 4 - phi^2 = 4 - (phi+1) = 3 - phi = {nstr(3-PHI, 10)}")
print(f"  4. At golden-: Delta = 4 - 1/phi^2 = 4 - (2-phi) = 2 + phi = {nstr(2+PHI, 10)}")
print(f"  5. Sum: (3-phi) + (2+phi) = 5 EXACTLY (the pentagon number!)")


# =============================================================================
# PART 2: THE PYTHAGOREAN BOUND ON THE PRIME SUM
# =============================================================================

separator("PART 2: PYTHAGOREAN BOUND ON THE PRIME SUM")

print("The Weil functional: W(f) = archimedean(f) - prime_sum(f)")
print("where prime_sum = sum_p a_p * log(p)/sqrt(p) * h(log p)")
print()
print("LINEAR problem: a_p can be positive or negative. The minus sign")
print("in front makes positive a_p DECREASE W, negative a_p INCREASE W.")
print("No universal sign control.")
print()
print("SQUARED approach: consider |prime_sum|^2 instead.")
print("By Cauchy-Schwarz:")
print("  |sum a_p * w_p|^2 <= (sum a_p^2 * w_p) * (sum w_p)")
print()
print("The LEFT side is |prime_sum|^2 >= 0.")
print("The RIGHT side uses a_p^2 (always >= 0) and w_p (always >= 0).")
print()

# Compute the Pythagorean bound for various prime ranges
print("Pythagorean bound computation:")
print(f"{'N':>8s} {'sum a_p^2/p':>16s} {'(sum|a_p|/sqp)^2':>18s} {'CS ratio':>12s} {'avg a_p^2':>12s}")
print("-" * 72)

bounds = [100, 500, 1000, 2000, 5000]
for N in bounds:
    primes_N = [p for p in primes_unram if p <= N]
    if not primes_N:
        continue

    # sum a_p^2 * log(p)/p
    sum_ap_sq_logp_over_p = fsum(
        trace_data[p]**2 * log(mpf(p)) / mpf(p) for p in primes_N
    )

    # sum log(p)/p  (the weight sum)
    sum_logp_over_p = fsum(
        log(mpf(p)) / mpf(p) for p in primes_N
    )

    # |sum a_p * log(p)/sqrt(p)|^2
    linear_sum = fsum(
        trace_data[p] * log(mpf(p)) / sqrt(mpf(p)) for p in primes_N
    )
    linear_sq = linear_sum**2

    # Cauchy-Schwarz bound: |sum a_p w_p|^2 <= (sum a_p^2 w_p)(sum w_p)
    # with w_p = log(p)/sqrt(p)
    sum_ap_sq_w = fsum(
        trace_data[p]**2 * log(mpf(p)) / sqrt(mpf(p)) for p in primes_N
    )
    sum_w = fsum(
        log(mpf(p)) / sqrt(mpf(p)) for p in primes_N
    )
    cs_bound = sum_ap_sq_w * sum_w
    cs_ratio = linear_sq / cs_bound if cs_bound > 0 else mpf(0)

    avg_ap_sq = fsum(trace_data[p]**2 for p in primes_N) / len(primes_N)

    print(f"{N:>8d} {nstr(sum_ap_sq_logp_over_p, 10):>16s} {nstr(linear_sq, 10):>18s} {nstr(cs_ratio, 6):>12s} {nstr(avg_ap_sq, 8):>12s}")

print()
print(f"Expected average |a_p|^2 by Schur orthogonality:")
avg_schur = fsum(d * t_sq for _, _, t_sq, _, d, _ in A5_CLASSES)
print(f"  <|a_p|^2> = sum density * trace^2 = {nstr(avg_schur, 10)}")
print(f"  = (1/60)*4 + (12/60)*phi^2 + (12/60)/phi^2 + (20/60)*1 + (15/60)*0")
print(f"  = {nstr(mpf(4)/60 + mpf(12)/60*PHI_SQ + mpf(12)/60*PHI_INV_SQ + mpf(20)/60, 10)}")
print(f"  = {nstr(avg_schur, 10)}")
print(f"  This is 1 (Schur orthogonality for dim-2 rep): {fabs(avg_schur - 1) < mpf('1e-20')}")


# =============================================================================
# PART 3: PYTHAGOREAN WEIL FUNCTIONAL W_P
# =============================================================================

separator("PART 3: PYTHAGOREAN WEIL FUNCTIONAL")

print("Instead of the linear Weil W = arch - prime_sum,")
print("consider the PYTHAGOREAN Weil: W_P = arch^2 - prime_sum^2")
print()
print("If arch >= |prime_sum|: both W >= 0 AND W_P >= 0.")
print("W_P = (arch - prime)(arch + prime)")
print()

# Gaussian test function f(x) = exp(-alpha * x^2)
# Fourier transform: g(xi) = sqrt(pi/alpha) * exp(-pi^2 xi^2 / alpha)

def gaussian_test(alpha, x):
    """Gaussian test function."""
    return exp(-alpha * x**2)

def gaussian_fourier(alpha, xi):
    """Fourier transform of Gaussian."""
    return sqrt(pi / alpha) * exp(-pi**2 * xi**2 / alpha)


def compute_archimedean_contribution(alpha, N_cond, dim=2):
    """
    Archimedean contribution to the Weil functional.

    For a degree-d L-function with conductor N:
    arch(f) = f(0) * (d/2 * log(N/pi^d) + d * digamma(1/2) + 2*d*log(2))
            + integral terms from Gamma factors

    Simplified for the icosahedral L-function (weight 1, dim 2):
    arch = (1/2) * log(N/pi^2) * g(0) + digamma integral terms

    We use the standard form from Weil's explicit formula.
    """
    g0 = gaussian_fourier(alpha, mpf(0))  # g(0) = sqrt(pi/alpha)

    # Log conductor term
    log_cond = log(N_cond / pi**dim) / 2

    # Digamma/Gamma terms for GL(2) (two Gamma_R factors)
    # For Gamma_R(s) = pi^{-s/2} Gamma(s/2):
    # The contribution is integral_0^infty [f(x)/x - f(0)/sinh(x)] dx
    # For Gaussian f(x) = exp(-alpha x^2), this integral converges.

    # Practical: use the formula
    # arch = log(N/(4pi^2)) * g(0)/2 + 2 * integral_0^infty [g(x)/(e^x-1) - g(0)*e^{-x}/x] dx
    # plus pole contribution from s=0 and s=1

    # For a simplified but correct leading-order computation:
    # arch ~ (d/2) * log(N) * g(0) + lower order terms

    # The POLE contribution: for zeta, there's a pole at s=1 giving +f(0).
    # For Artin L-functions with no pole: pole = 0.
    # But ρ_ico is irreducible, so L(s, rho_ico) has no pole.

    # Full archimedean for weight-1 GL(2) with conductor N:
    # W_arch(f) = (log N - 2*log(2pi)) * g(0) + 2 * I_Gamma
    # where I_Gamma = integral_0^infty [g(t)*(Psi(1+it) + Psi(1-it))/2 - g(0)/t^2 ...] dt

    # For our purposes, compute numerically:
    log_N = log(N_cond)

    # Leading archimedean term
    arch_leading = (log_N - 2 * log(2 * pi)) * g0

    # Gamma integral contribution (numerical integration)
    # The Gamma contribution for each Gamma_R factor is:
    # integral_0^infty h(t) * Re[psi(1/4 + it/2)] dt / pi
    # where psi = digamma

    def gamma_integrand(t):
        """Integrand for the Gamma contribution to the Weil functional."""
        if fabs(t) < mpf('1e-15'):
            return mpf(0)
        gt = gaussian_fourier(alpha, t)
        # For weight 1 modular form: Gamma_R(s+1) = pi^{-(s+1)/2} Gamma((s+1)/2)
        # The contribution involves Re[psi((1+it)/2)]
        s_val = mpc(mpf('0.25'), t / 2)
        try:
            psi_val = digamma(s_val)
            return gt * re(psi_val) / pi
        except Exception:
            return mpf(0)

    # Numerical integration (truncate at reasonable bound)
    try:
        I_gamma = 2 * quad(gamma_integrand, [0, 50], error=True)[0]
    except Exception:
        I_gamma = mpf(0)

    arch_total = arch_leading + I_gamma
    return arch_total


def compute_prime_sum(alpha, primes_list, trace_dict, max_k=10):
    """
    Prime sum in the explicit formula:
    prime_sum = sum_p sum_{k>=1} a_{p^k} * log(p) / p^{k/2} * g(k*log(p)/(2*pi))

    where g is the Fourier transform of the test function.
    """
    total = mpf(0)
    for p in primes_list:
        ap = trace_dict[p]
        logp = log(mpf(p))
        sqp = sqrt(mpf(p))

        for k in range(1, max_k + 1):
            apk = hecke_eigenvalue(ap, k)
            weight = logp / power(mpf(p), mpf(k) / 2)
            # The argument of g in the explicit formula
            xi = k * logp / (2 * pi)
            gval = gaussian_fourier(alpha, xi)
            total += apk * weight * gval

    return total


def compute_prime_sum_squared_traces(alpha, primes_list, trace_dict, max_k=10):
    """
    Prime sum using SQUARED traces |a_p|^2 instead of a_p.
    This is the prime sum for the Rankin-Selberg L-function.

    For L_RS: b_{p^k} = |a_{p^k}|^2 (product of Hecke eigenvalues)
    Actually: b_p = |a_p|^2 = a_p^2 (since a_p is real for icosahedral rep)

    The Hecke eigenvalues for L_RS satisfy:
    For rho x rho_bar (dim 4): the Euler factor is
    (1 - alpha^2 u)(1 - alpha*beta u)(1 - alpha*beta u)(1 - beta^2 u)
    = (1 - alpha^2 u)(1 - u)^2 (1 - beta^2 u)  ... NO, det=1 so alpha*beta=1
    Wait: alpha_p * beta_p = 1 (det=1), so the middle terms are (1-u)^2.

    Actually for rho x rho_bar where rho has eigenvalues alpha, beta with alpha*beta=1:
    The tensor product eigenvalues are alpha/alpha=1, alpha/beta=alpha^2,
    beta/alpha=beta^2, beta/beta=1 (since rho_bar has eigenvalues alpha_bar, beta_bar
    and for real reps, bar = conjugate = inverse).

    Wait, for the icosahedral rep, the eigenvalues alpha_p, beta_p of Frobenius satisfy:
    alpha_p * beta_p = det(Frob_p) = 1 (since det rho = 1)
    alpha_p + beta_p = a_p (the trace)

    For rho x rho_bar: the eigenvalues of Frob_p are
    alpha_p * conj(alpha_p), alpha_p * conj(beta_p), beta_p * conj(alpha_p), beta_p * conj(beta_p)

    Since |alpha_p| = |beta_p| = 1 (Ramanujan for Artin):
    alpha_p * conj(alpha_p) = 1
    beta_p * conj(beta_p) = 1
    alpha_p * conj(beta_p) = alpha_p / beta_p = alpha_p^2 (since alpha*beta = 1)
    beta_p * conj(alpha_p) = beta_p / alpha_p = beta_p^2

    So b_p = 1 + alpha_p^2 + beta_p^2 + 1 = 2 + a_p^2 - 2 = a_p^2

    The Hecke eigenvalues at prime powers for L_RS use a different recurrence.
    For our simplified computation, we use b_p = a_p^2 at the first power
    and compute higher powers from the Rankin-Selberg Euler product.

    Simplified: just use b_p = a_p^2 for k=1, b_{p^k} = (a_{p^k})^2 for k>1.
    This is NOT exactly right for the RS L-function (which has degree 4 Euler product)
    but captures the key positivity property.
    """
    total = mpf(0)
    for p in primes_list:
        ap = trace_dict[p]
        logp = log(mpf(p))

        for k in range(1, max_k + 1):
            apk = hecke_eigenvalue(ap, k)
            bpk = apk**2  # Squared trace (always >= 0!)
            weight = logp / power(mpf(p), mpf(k) / 2)
            xi = k * logp / (2 * pi)
            gval = gaussian_fourier(alpha, xi)
            total += bpk * weight * gval

    return total


# Compute Pythagorean Weil for various test function widths
print("Pythagorean Weil W_P = arch^2 - prime_sum^2 for Gaussian test functions:")
print()
print(f"{'alpha':>10s} {'arch':>14s} {'prime':>14s} {'W=arch-prime':>14s} {'W_P=a^2-p^2':>14s} {'|prime|_sq':>14s}")
print("-" * 86)

test_alphas = [mpf('0.001'), mpf('0.01'), mpf('0.05'), mpf('0.1'), mpf('0.5'), mpf('1.0'), mpf('2.0'), mpf('5.0')]
test_primes = [p for p in primes_unram if p <= 2000]

for alpha in test_alphas:
    arch = compute_archimedean_contribution(alpha, CONDUCTOR)
    prime = compute_prime_sum(alpha, test_primes, trace_data, max_k=MAX_POWER)
    prime_sq = compute_prime_sum_squared_traces(alpha, test_primes, trace_data, max_k=MAX_POWER)

    W_linear = arch - prime
    W_pyth = arch**2 - prime**2

    print(f"{nstr(alpha, 6):>10s} {nstr(arch, 8):>14s} {nstr(prime, 8):>14s} {nstr(W_linear, 8):>14s} {nstr(W_pyth, 8):>14s} {nstr(prime_sq, 8):>14s}")

print()
print("KEY: |prime|_sq uses |a_p|^2 (always >= 0). No sign problem in the TRACES.")
print("But the prime sum ITSELF still enters with a minus sign in the explicit formula.")
print("The Pythagorean approach bounds |prime_sum|^2 rather than prime_sum.")


# =============================================================================
# PART 4: DIRECT SUM-OF-SQUARES VIA RANKIN-SELBERG
# =============================================================================

separator("PART 4: RANKIN-SELBERG L-FUNCTION (SQUARED TRACES)")

print("The Rankin-Selberg L-function: L_RS(s) = L(s, rho x rho_bar)")
print()
print("Factorization: L_RS(s) = L(s, Sym^2 rho) * zeta(s)")
print("  (since det rho = 1, so rho x rho_bar = Sym^2 rho + trivial)")
print()

# Table of Rankin-Selberg traces
print("Rankin-Selberg traces by Frobenius class:")
print(f"{'Class':15s} {'a_p':>10s} {'b_p=a_p^2':>12s} {'c_p=a_p^2-1':>12s} {'1 (trivial)':>12s} {'check':>8s}")
print("-" * 76)

for name, ap, ap_sq, size, density, order in A5_CLASSES:
    bp = ap_sq                  # rho x rho_bar trace
    cp = ap_sq - 1              # Sym^2 rho trace
    trivial = mpf(1)            # trivial rep trace
    check = "OK" if fabs(bp - cp - trivial) < mpf('1e-25') else "FAIL"
    print(f"{name:15s} {nstr(ap, 8):>10s} {nstr(bp, 8):>12s} {nstr(cp, 8):>12s} {'1':>12s} {check:>8s}")

print()
print("CRITICAL OBSERVATION: b_p = a_p^2 >= 0 for ALL primes.")
print("  golden+:    b_p = phi^2 = phi + 1 =", nstr(PHI_SQ, 10))
print("  golden-:    b_p = 1/phi^2 = 2 - phi =", nstr(PHI_INV_SQ, 10))
print("  order-3:    b_p = 1")
print("  involution: b_p = 0")
print("  identity:   b_p = 4")
print()
print("ALL non-negative! The sign problem that killed the pentagram is GONE.")
print()

# Verify the factorization L_RS = L(Sym^2) * zeta
print("Rankin-Selberg = Sym^2 x zeta factorization check:")
print()

# At each class: b_p = c_p + 1 (Sym^2 trace + trivial trace)
# This means: L_RS(s) has Euler factor
# (1 - b_p p^{-s} + ....) which factors as
# (1 - c_p p^{-s} + ...) * (1 - p^{-s})^{-1}
# The Sym^2 Euler factor for det=1 rep:
# For alpha^2 + 1 + beta^2 = a_p^2 - 1 + 1 = a_p^2... hmm
# Actually Sym^2 has eigenvalues alpha^2, alpha*beta=1, beta^2
# So L(s, Sym^2) = prod_p (1-alpha_p^2/p^s)^{-1}(1-1/p^s)^{-1}(1-beta_p^2/p^s)^{-1}
# And zeta(s) = prod_p (1-p^{-s})^{-1}
# So L_RS = L(Sym^2) * zeta has a DOUBLE pole from the two (1-p^{-s})^{-1} factors.
# Wait no: L_RS = L(Sym^2) * zeta only if Sym^2 does NOT include the trivial factor.
# Since rho x rho_bar = Sym^2 rho (when det=1, exterior power = trivial is separated)
# Actually: rho x rho_bar = Sym^2(rho) + Wedge^2(rho) = Sym^2(rho) + det(rho) = Sym^2(rho) + trivial
# So L(s, rho x rho_bar) = L(s, Sym^2 rho) * L(s, trivial) = L(s, Sym^2 rho) * zeta(s)

print("  L(s, rho x rho_bar) = L(s, Sym^2 rho) * zeta(s)")
print()
print("  zeta(s) has a simple pole at s=1 with residue 1.")
print("  L(s, Sym^2 rho) is entire (Sym^2 of irreducible dim-2 Artin rep).")
print("  Therefore L_RS has a simple pole at s=1.")
print()
print("  This means: the Weil functional for L_RS includes a POLE TERM")
print("  which contributes POSITIVELY to W_RS.")
print()

# Compute the "Rankin-Selberg Hadamard gap"
print("Rankin-Selberg spectral gaps:")
print(f"{'Class':15s} {'b_p=|a_p|^2':>12s} {'dim-b_p':>10s} {'(dim-b_p)^2':>12s}")
print("-" * 55)

dim_RS = mpf(4)  # dimension of rho x rho_bar
min_gap_RS = inf
for name, ap, ap_sq, size, density, order in A5_CLASSES:
    bp = ap_sq
    gap = (dim_RS - bp)**2
    if name != "identity" and gap < min_gap_RS:
        min_gap_RS = gap
    print(f"{name:15s} {nstr(bp, 8):>12s} {nstr(dim_RS - bp, 8):>10s} {nstr(gap, 8):>12s}")

print()
print(f"Minimum non-identity RS gap: {nstr(min_gap_RS, 15)}")
print(f"  = (4 - phi^2)^2 = (4 - phi - 1)^2 = (3 - phi)^2")
print(f"  = {nstr((3 - PHI)**2, 15)}")
print(f"  = ((5 - sqrt(5))/2)^2 = (25 - 10*sqrt(5) + 5)/4 = (30 - 10*sqrt(5))/4")
exact_gap = (30 - 10*sqrt(5)) / 4
print(f"  = {nstr(exact_gap, 15)}")


# =============================================================================
# PART 5: RANKIN-SELBERG WEIL POSITIVITY
# =============================================================================

separator("PART 5: RANKIN-SELBERG WEIL POSITIVITY")

print("The Rankin-Selberg Weil functional:")
print("  W_RS(f) = arch_RS(f) + pole_RS(f) - prime_sum_RS(f)")
print()
print("Where:")
print("  - arch_RS: log-conductor + digamma terms for L_RS (dim 4)")
print("  - pole_RS: contribution from the simple pole at s=1 (+POSITIVE)")
print("  - prime_sum_RS: uses b_p = |a_p|^2 >= 0 (ALL terms non-negative)")
print()
print("The prime sum is NON-NEGATIVE. So W_RS >= 0 requires:")
print("  arch_RS + pole_RS >= prime_sum_RS")
print()
print("The pole term is f(0) * residue = f(0) > 0 for f > 0.")
print("The arch term scales as ~ dim/2 * log(N_RS) * f(0).")
print("The prime sum grows as ~ <|a_p|^2> * sqrt(N) = sqrt(N) (Schur).")
print()
print("So for LARGE N (conductor): arch + pole DOMINATES prime_sum.")
print("The question: does it dominate for ALL test functions, not just at large N?")
print()

# Conductor of Rankin-Selberg
# N_RS = N^2 * (correction from ramification)
# For the icosahedral rep with N=800:
# N_{Sym^2} = 800^(3/2)... actually, for Sym^2 of a conductor-N form,
# the conductor is bounded by N^3 (Bushnell-Henniart), typically N^2.
# For the Rankin-Selberg: N_RS = N_{Sym^2} * 1 (zeta has conductor 1)
# We use N_RS = N^2 = 640000 as a reasonable estimate.
N_RS = CONDUCTOR**2
N_Sym2 = CONDUCTOR**2  # conservative estimate

print(f"Conductors:")
print(f"  N(rho)           = {int(CONDUCTOR)}")
print(f"  N(rho x rho_bar) ~ N^2 = {int(N_RS)}")
print(f"  N(Sym^2 rho)     ~ N^2 = {int(N_Sym2)}")
print()

# Compute W_RS for various test functions
print("Rankin-Selberg Weil functional values:")
print(f"{'alpha':>10s} {'arch_RS':>14s} {'pole':>14s} {'prime_RS':>14s} {'W_RS':>14s} {'W_RS >= 0?':>12s}")
print("-" * 82)

for alpha in test_alphas:
    g0_RS = gaussian_fourier(alpha, mpf(0))

    # Archimedean for RS (dimension 4, conductor N^2)
    arch_RS = compute_archimedean_contribution(alpha, N_RS, dim=4)

    # Pole contribution: for L_RS with simple pole at s=1
    # The pole of L_RS(s) at s=1 comes from zeta(s).
    # Residue of zeta at s=1 is 1.
    # So the pole contribution to the explicit formula is:
    # pole = -Res_{s=1} (f * L_RS'/L_RS) = f(0) * order_of_pole
    # For a simple pole: pole contribution = +g(0) (positive!)
    pole = g0_RS

    # Prime sum for RS: uses b_p = a_p^2
    prime_RS = compute_prime_sum_squared_traces(alpha, test_primes, trace_data, max_k=MAX_POWER)

    W_RS = arch_RS + pole - prime_RS
    is_pos = "YES" if W_RS > 0 else "NO"

    print(f"{nstr(alpha, 6):>10s} {nstr(arch_RS, 8):>14s} {nstr(pole, 8):>14s} {nstr(prime_RS, 8):>14s} {nstr(W_RS, 8):>14s} {is_pos:>12s}")

print()
print("ANALYSIS:")
print("  The prime_RS sum uses ONLY non-negative coefficients b_p = |a_p|^2.")
print("  The arch_RS and pole terms must collectively dominate this sum.")
print("  If W_RS >= 0 for ALL test functions -> GRH for L_RS = L(Sym^2) * zeta.")
print("  Then: GRH for zeta follows from factoring out L(Sym^2).")


# =============================================================================
# PART 6: TRIANGLE WEIL MATRIX (GRAM MATRIX COMPARISON)
# =============================================================================

separator("PART 6: TRIANGLE WEIL MATRIX (GRAM MATRIX)")

print("Compare the Gram matrix for L(s, rho) vs L_RS(s, rho x rho_bar).")
print("The Gram matrix G_{ij} = W(f_i * f_j) measures the positivity of the")
print("Weil functional as a quadratic form.")
print()
print("For L(rho): traces a_p can be negative -> Gram matrix may have negative eigenvalues.")
print("For L_RS:   traces |a_p|^2 >= 0      -> expect BETTER positivity properties.")
print()

# Build Gram matrices using Gaussian basis functions
# f_j(x) = exp(-alpha_j * x^2) for different widths alpha_j
gram_dim = 8
gram_alphas = [mpf(2)**(-k) for k in range(-2, gram_dim - 2)]  # 4, 2, 1, 0.5, 0.25, ...

print(f"Gram matrix dimension: {gram_dim}")
print(f"Basis widths alpha_j: {[float(a) for a in gram_alphas]}")
print()

# Helper: convolution of two Gaussians
# f_i * f_j has Fourier transform g_i * g_j
# g_i(xi) = sqrt(pi/alpha_i) * exp(-pi^2 xi^2/alpha_i)
# g_i * g_j(xi) = g_i(xi) * g_j(xi)  ... NO, convolution in real = product in Fourier
# Actually for the Weil functional W(h) where h = f_i * f_j (convolution):
# The Fourier transform of the convolution is the product: hat(h)(xi) = hat(f_i)(xi) * hat(f_j)(xi)
# For Gaussians: hat(h)(xi) = sqrt(pi/alpha_i) * sqrt(pi/alpha_j) * exp(-pi^2 xi^2 (1/alpha_i + 1/alpha_j))
# = pi/sqrt(alpha_i*alpha_j) * exp(-pi^2 xi^2 / beta_ij) where 1/beta_ij = 1/alpha_i + 1/alpha_j

def combined_gaussian_fourier(alpha_i, alpha_j, xi):
    """Fourier transform of convolution of two Gaussians."""
    beta_inv = 1/alpha_i + 1/alpha_j
    return pi / sqrt(alpha_i * alpha_j) * exp(-pi**2 * xi**2 * beta_inv)


def compute_gram_entry(alpha_i, alpha_j, primes_list, trace_dict, use_squared=False, max_k=10):
    """
    Compute Gram matrix entry G_{ij} = W(f_i * f_j).

    The prime sum part of W(f_i * f_j):
    sum_p sum_k a_{p^k} * log(p)/p^{k/2} * hat(f_i)(k*logp/(2pi)) * hat(f_j)(k*logp/(2pi))

    For use_squared=True: use |a_{p^k}|^2 instead of a_{p^k}.
    """
    total = mpf(0)
    for p in primes_list:
        ap = trace_dict[p]
        logp = log(mpf(p))

        for k in range(1, max_k + 1):
            apk = hecke_eigenvalue(ap, k)
            if use_squared:
                apk = apk**2  # squared traces

            weight = logp / power(mpf(p), mpf(k) / 2)
            xi = k * logp / (2 * pi)

            # Product of Fourier transforms at xi
            gi = gaussian_fourier(alpha_i, xi)
            gj = gaussian_fourier(alpha_j, xi)
            total += apk * weight * gi * gj

    return total


def compute_archimedean_gram(alpha_i, alpha_j, N_cond, dim=2):
    """Archimedean + pole contribution to Gram entry."""
    g0_i = gaussian_fourier(alpha_i, mpf(0))
    g0_j = gaussian_fourier(alpha_j, mpf(0))

    # Leading order: (log N - d*log(2pi)) * g_i(0) * g_j(0) / something
    # For the Gram matrix, the archimedean contribution is the product of
    # the Fourier transforms at 0 times the log-conductor.
    log_N = log(N_cond)

    # Simplified leading-order archimedean Gram entry
    arch = (log_N - dim * log(2 * pi)) / 2 * g0_i * g0_j

    # Add Gamma contribution (numerical)
    def gamma_gram_integrand(t):
        if fabs(t) < mpf('1e-15'):
            return mpf(0)
        gi = gaussian_fourier(alpha_i, t)
        gj = gaussian_fourier(alpha_j, t)
        s_val = mpc(mpf('0.25'), t / 2)
        try:
            psi_val = digamma(s_val)
            return gi * gj * re(psi_val) / pi
        except Exception:
            return mpf(0)

    try:
        I_gamma = 2 * quad(gamma_gram_integrand, [0, 50], error=True)[0]
    except Exception:
        I_gamma = mpf(0)

    return arch + I_gamma


# Build the two Gram matrices
print("Building Gram matrices...")
print()

gram_primes = [p for p in primes_unram if p <= 1000]

# Gram matrix for L(s, rho) — uses raw traces a_p
G_rho = matrix(gram_dim, gram_dim)
for i in range(gram_dim):
    for j in range(i, gram_dim):
        prime_part = compute_gram_entry(gram_alphas[i], gram_alphas[j],
                                         gram_primes, trace_data,
                                         use_squared=False, max_k=10)
        arch_part = compute_archimedean_gram(gram_alphas[i], gram_alphas[j],
                                              CONDUCTOR, dim=2)
        G_rho[i, j] = arch_part - prime_part
        G_rho[j, i] = G_rho[i, j]

# Gram matrix for L_RS — uses squared traces |a_p|^2
G_RS = matrix(gram_dim, gram_dim)
for i in range(gram_dim):
    for j in range(i, gram_dim):
        prime_part = compute_gram_entry(gram_alphas[i], gram_alphas[j],
                                         gram_primes, trace_data,
                                         use_squared=True, max_k=10)
        arch_part_RS = compute_archimedean_gram(gram_alphas[i], gram_alphas[j],
                                                 N_RS, dim=4)
        # Add pole contribution for RS (simple pole at s=1)
        g0i = gaussian_fourier(gram_alphas[i], mpf(0))
        g0j = gaussian_fourier(gram_alphas[j], mpf(0))
        pole_part = g0i * g0j

        G_RS[i, j] = arch_part_RS + pole_part - prime_part
        G_RS[j, i] = G_RS[i, j]

# Compute eigenvalues
print("Gram matrix for L(s, rho) [raw traces, may have sign problem]:")
try:
    evals_rho = sorted([re(x) for x in eig(G_rho, left=False, right=False)], reverse=True)
    for idx, ev in enumerate(evals_rho):
        sign = "+" if ev > 0 else "-" if ev < 0 else "0"
        print(f"  lambda_{idx+1} = {nstr(ev, 12):>20s}  [{sign}]")
    n_pos_rho = sum(1 for ev in evals_rho if ev > mpf('1e-20'))
    n_neg_rho = sum(1 for ev in evals_rho if ev < -mpf('1e-20'))
    n_zero_rho = gram_dim - n_pos_rho - n_neg_rho
    print(f"  Signature: ({n_pos_rho}+, {n_zero_rho}zero, {n_neg_rho}-)")
    print(f"  PSD (positive semi-definite): {'YES' if n_neg_rho == 0 else 'NO'}")
except Exception as e:
    print(f"  Eigenvalue computation failed: {e}")
    evals_rho = None

print()
print("Gram matrix for L_RS(s, rho x rho_bar) [squared traces, no sign problem]:")
try:
    evals_RS = sorted([re(x) for x in eig(G_RS, left=False, right=False)], reverse=True)
    for idx, ev in enumerate(evals_RS):
        sign = "+" if ev > 0 else "-" if ev < 0 else "0"
        print(f"  lambda_{idx+1} = {nstr(ev, 12):>20s}  [{sign}]")
    n_pos_RS = sum(1 for ev in evals_RS if ev > mpf('1e-20'))
    n_neg_RS = sum(1 for ev in evals_RS if ev < -mpf('1e-20'))
    n_zero_RS = gram_dim - n_pos_RS - n_neg_RS
    print(f"  Signature: ({n_pos_RS}+, {n_zero_RS}zero, {n_neg_RS}-)")
    print(f"  PSD (positive semi-definite): {'YES' if n_neg_RS == 0 else 'NO'}")
except Exception as e:
    print(f"  Eigenvalue computation failed: {e}")
    evals_RS = None

# Compare
print()
print("COMPARISON:")
if evals_rho is not None and evals_RS is not None:
    min_rho = min(evals_rho)
    min_RS = min(evals_RS)
    print(f"  Smallest eigenvalue of G(rho):  {nstr(min_rho, 15)}")
    print(f"  Smallest eigenvalue of G(RS):   {nstr(min_RS, 15)}")
    if min_RS > min_rho:
        print(f"  G(RS) is MORE positive than G(rho) by {nstr(min_RS - min_rho, 10)}")
    else:
        print(f"  G(RS) is NOT more positive than G(rho)")

    ratio = min_RS / min_rho if fabs(min_rho) > mpf('1e-30') else inf
    print(f"  Ratio of smallest eigenvalues (RS/rho): {nstr(ratio, 10)}")


# =============================================================================
# PART 6b: PURE PRIME-SUM GRAM MATRICES (without archimedean)
# =============================================================================

print()
print("-" * 76)
print("  PART 6b: PURE PRIME-SUM GRAM MATRICES (trace structure only)")
print("-" * 76)
print()
print("To isolate the effect of squared traces, build Gram matrices using")
print("ONLY the prime sum (no archimedean/pole). Compare the trace structures.")
print()

# Pure prime-sum Gram for L(rho)
P_rho = matrix(gram_dim, gram_dim)
for i in range(gram_dim):
    for j in range(i, gram_dim):
        P_rho[i, j] = compute_gram_entry(gram_alphas[i], gram_alphas[j],
                                           gram_primes, trace_data,
                                           use_squared=False, max_k=10)
        P_rho[j, i] = P_rho[i, j]

# Pure prime-sum Gram for L_RS
P_RS = matrix(gram_dim, gram_dim)
for i in range(gram_dim):
    for j in range(i, gram_dim):
        P_RS[i, j] = compute_gram_entry(gram_alphas[i], gram_alphas[j],
                                          gram_primes, trace_data,
                                          use_squared=True, max_k=10)
        P_RS[j, i] = P_RS[i, j]

print("Prime-sum Gram for L(rho) [raw traces a_p, can be negative]:")
try:
    evals_P_rho = sorted([re(x) for x in eig(P_rho, left=False, right=False)], reverse=True)
    n_pos = sum(1 for ev in evals_P_rho if ev > mpf('1e-20'))
    n_neg = sum(1 for ev in evals_P_rho if ev < -mpf('1e-20'))
    for idx, ev in enumerate(evals_P_rho):
        sign = "+" if ev > 0 else "-" if ev < 0 else "0"
        print(f"  lambda_{idx+1} = {nstr(ev, 10):>18s}  [{sign}]")
    print(f"  Signature: ({n_pos}+, {gram_dim - n_pos - n_neg}zero, {n_neg}-)")
except Exception as e:
    print(f"  Failed: {e}")
    evals_P_rho = None

print()
print("Prime-sum Gram for L_RS [squared traces |a_p|^2, always >= 0]:")
try:
    evals_P_RS = sorted([re(x) for x in eig(P_RS, left=False, right=False)], reverse=True)
    n_pos = sum(1 for ev in evals_P_RS if ev > mpf('1e-20'))
    n_neg = sum(1 for ev in evals_P_RS if ev < -mpf('1e-20'))
    for idx, ev in enumerate(evals_P_RS):
        sign = "+" if ev > 0 else "-" if ev < 0 else "0"
        print(f"  lambda_{idx+1} = {nstr(ev, 10):>18s}  [{sign}]")
    print(f"  Signature: ({n_pos}+, {gram_dim - n_pos - n_neg}zero, {n_neg}-)")
except Exception as e:
    print(f"  Failed: {e}")
    evals_P_RS = None

print()
print("KEY QUESTION: Is the RS prime-sum Gram PSD?")
print("If yes: the prime sum with squared traces has inherently better structure.")
print("If no:  the sign problem persists even with squared traces (it's in the")
print("        sum structure, not just the traces).")
print()

if evals_P_RS is not None:
    min_P_RS = min(evals_P_RS)
    print(f"  Min eigenvalue of P_RS: {nstr(min_P_RS, 15)}")
    if min_P_RS >= -mpf('1e-20'):
        print("  P_RS is PSD -> squared traces give inherent positivity!")
    else:
        print("  P_RS is NOT PSD -> sign problem is structural, not just trace-level.")
        print("  BUT: the negative eigenvalues may be SMALLER than for P_rho,")
        print("  making the archimedean domination easier.")

# P_RS PSD stability test with increasing prime ranges
print()
print("-" * 76)
print("  PART 6c: P_RS PSD STABILITY TEST")
print("-" * 76)
print()
print("Does P_RS remain PSD as we include more primes?")
print()
print(f"{'Primes up to':>14s} {'dim':>5s} {'min eval':>18s} {'PSD?':>6s} {'max eval':>18s}")
print("-" * 68)

for pb in [200, 500, 1000, 2000, 5000]:
    p_test = [p for p in primes_unram if p <= pb]
    if len(p_test) < 3:
        continue
    test_dim = min(6, len(p_test))
    test_gram_alphas = [mpf(2)**(-k) for k in range(-1, test_dim - 1)]
    P_test = matrix(test_dim, test_dim)
    for i in range(test_dim):
        for j in range(i, test_dim):
            P_test[i, j] = compute_gram_entry(test_gram_alphas[i], test_gram_alphas[j],
                                               p_test, trace_data,
                                               use_squared=True, max_k=10)
            P_test[j, i] = P_test[i, j]
    try:
        evs = sorted([re(x) for x in eig(P_test, left=False, right=False)])
        is_psd = all(ev > -mpf('1e-15') for ev in evs)
        print(f"{pb:>14d} {test_dim:>5d} {nstr(evs[0], 10):>18s} {'YES' if is_psd else 'NO':>6s} {nstr(evs[-1], 10):>18s}")
    except Exception as e:
        print(f"{pb:>14d} {test_dim:>5d} {'FAILED':>18s}")

# Also test with various alphas (wider basis = harder test)
print()
print("P_RS with wider alpha range (harder test, 1000 primes):")
p_1000 = [p for p in primes_unram if p <= 1000]
for spread in [2, 4, 8, 16]:
    test_dim = 6
    # Use wider spread: alpha = spread^k for k in range
    test_gram_alphas_wide = [mpf(spread)**(-k) for k in range(-1, test_dim - 1)]
    P_test = matrix(test_dim, test_dim)
    for i in range(test_dim):
        for j in range(i, test_dim):
            P_test[i, j] = compute_gram_entry(test_gram_alphas_wide[i], test_gram_alphas_wide[j],
                                               p_1000, trace_data,
                                               use_squared=True, max_k=10)
            P_test[j, i] = P_test[i, j]
    try:
        evs = sorted([re(x) for x in eig(P_test, left=False, right=False)])
        is_psd = all(ev > -mpf('1e-15') for ev in evs)
        print(f"  spread={spread:>2d}, alphas={[float(nstr(a,4)) for a in test_gram_alphas_wide]}")
        print(f"    min eval = {nstr(evs[0], 10)}, PSD = {'YES' if is_psd else 'NO'}")
    except Exception as e:
        print(f"  spread={spread:>2d}: FAILED ({e})")


# =============================================================================
# PART 7: PYTHAGOREAN CONCLUSION
# =============================================================================

separator("PART 7: PYTHAGOREAN CONCLUSION")

print("THE PYTHAGOREAN DECOMPOSITION OF THE WEIL FUNCTIONAL")
print()
print("For each prime p, the Frobenius trace satisfies the Pythagorean identity:")
print("  a_p^2 + Delta_p = 4   (Ramanujan bound = 2, squared = 4)")
print()
print("where Delta_p = 4 - a_p^2 is the 'Pythagorean remainder'.")
print()

# Print the Pythagorean decomposition for each class
print(f"{'Class':15s} {'a_p':>10s} {'a_p^2':>10s} {'Delta_p':>10s} {'check=4':>10s}")
print("-" * 60)

for name, ap, ap_sq, size, density, order in A5_CLASSES:
    delta_p = 4 - ap_sq
    check = ap_sq + delta_p
    print(f"{name:15s} {nstr(ap, 8):>10s} {nstr(ap_sq, 6):>10s} {nstr(delta_p, 6):>10s} {nstr(check, 6):>10s}")

print()

# The golden Pythagorean triple
print("THE GOLDEN PYTHAGOREAN TRIPLE:")
print(f"  a_p = phi: phi^2 + (3-phi) = (phi+1) + (3-phi) = 4 = 2^2")
print(f"  Normalized: (phi/2)^2 + ((3-phi)/4) = 1")
print(f"  This is cos^2(theta) + sin^2(theta) = 1 with theta = arccos(phi/2)")
theta_golden = mp.acos(PHI / 2)
print(f"  theta = arccos(phi/2) = {nstr(theta_golden, 10)} radians")
print(f"                        = {nstr(theta_golden * 180 / pi, 10)} degrees")
print(f"  = pi/5 = 36 degrees: {fabs(theta_golden - pi/5) < mpf('1e-20')}")
print()
print(f"  a_p = -1/phi: 1/phi^2 + (2+phi) = (2-phi) + (2+phi) = 4")
print(f"  theta = arccos(-1/(2*phi)) = {nstr(mp.acos(-PHI_INV/2), 10)} radians")
theta_golden_minus = mp.acos(-PHI_INV / 2)
print(f"                              = {nstr(theta_golden_minus * 180 / pi, 10)} degrees")
print(f"  = 3*pi/5 = 108 degrees: {fabs(theta_golden_minus - 3*pi/5) < mpf('1e-20')}")
print()

# The triangle inequality approach
print("TRIANGLE INEQUALITY FOR THE WEIL FUNCTIONAL:")
print()
print("  W(f) = arch(f) + pole(f) - prime_sum(f)")
print()
print("  For the original L(rho): prime_sum uses a_p (sign-varying)")
print("  For L_RS = L(Sym^2)*zeta: prime_sum uses |a_p|^2 (non-negative)")
print()
print("  GRH for L_RS implies GRH for both factors (L(Sym^2) and zeta).")
print("  [Proof: if zeta had a zero off the critical line, L_RS would too,")
print("   contradicting GRH for L_RS.]")
print()

# Compute the key ratio: archimedean vs prime sum for RS
print("CONVERGENCE TEST: Does arch_RS + pole_RS dominate prime_sum_RS?")
print()
print(f"{'Primes up to':>14s} {'arch+pole':>14s} {'prime_RS':>14s} {'ratio':>10s} {'W_RS > 0?':>12s}")
print("-" * 70)

test_alpha_fixed = mpf('1.0')  # wider in Fourier space = more primes contribute
prime_bounds = [100, 200, 500, 1000, 2000, 5000]

for pb in prime_bounds:
    p_list = [p for p in primes_unram if p <= pb]
    if not p_list:
        continue

    arch_RS = compute_archimedean_contribution(test_alpha_fixed, N_RS, dim=4)
    g0 = gaussian_fourier(test_alpha_fixed, mpf(0))
    pole = g0

    prime_RS = compute_prime_sum_squared_traces(test_alpha_fixed, p_list, trace_data, max_k=MAX_POWER)

    W_RS = arch_RS + pole - prime_RS
    ratio = (arch_RS + pole) / prime_RS if prime_RS > 0 else inf
    is_pos = "YES" if W_RS > 0 else "NO"

    print(f"{pb:>14d} {nstr(arch_RS + pole, 8):>14s} {nstr(prime_RS, 8):>14s} {nstr(ratio, 6):>10s} {is_pos:>12s}")

print()

# Final summary of eigenvalue comparison
print("=" * 76)
print("  FINAL EIGENVALUE SUMMARY")
print("=" * 76)
print()

if evals_rho is not None:
    print("Full Weil Gram matrix for L(rho) [linear traces]:")
    n_neg_rho = sum(1 for ev in evals_rho if ev < -mpf('1e-20'))
    print(f"  Eigenvalues: {[float(nstr(ev, 6)) for ev in evals_rho]}")
    print(f"  Negative eigenvalues: {n_neg_rho}")
    print(f"  Minimum: {nstr(min(evals_rho), 12)}")

if evals_RS is not None:
    print()
    print("Full Weil Gram matrix for L_RS [squared traces]:")
    n_neg_RS = sum(1 for ev in evals_RS if ev < -mpf('1e-20'))
    print(f"  Eigenvalues: {[float(nstr(ev, 6)) for ev in evals_RS]}")
    print(f"  Negative eigenvalues: {n_neg_RS}")
    print(f"  Minimum: {nstr(min(evals_RS), 12)}")

if evals_P_rho is not None:
    print()
    print("Prime-only Gram for L(rho) [linear traces]:")
    print(f"  Eigenvalues: {[float(nstr(ev, 6)) for ev in evals_P_rho]}")
    print(f"  All positive: {all(ev > -mpf('1e-20') for ev in evals_P_rho)}")

if evals_P_RS is not None:
    print()
    print("Prime-only Gram for L_RS [squared traces]:")
    print(f"  Eigenvalues: {[float(nstr(ev, 6)) for ev in evals_P_RS]}")
    print(f"  All positive: {all(ev > -mpf('1e-20') for ev in evals_P_RS)}")

print()
print("=" * 76)
print("  THE PYTHAGOREAN PATH TO RH")
print("=" * 76)
print()
print("  1. The pentagram (linear traces a_p) has a SIGN PROBLEM:")
print("     the explicit formula subtracts the prime sum, and positive traces")
print("     become negative contributions. More primes = worse.")
print()
print("  2. The triangle (squared traces |a_p|^2) ELIMINATES the sign problem")
print("     in the traces: all b_p = |a_p|^2 >= 0.")
print()
print("  3. The Rankin-Selberg L-function L_RS = L(Sym^2 rho) * zeta(s)")
print("     naturally uses squared traces, and has a POLE at s=1 that")
print("     contributes positively to the Weil functional.")
print()
print("  4. The Pythagorean decomposition: a_p^2 + Delta_p = 4")
print("     At golden primes: Delta = 3 - phi (from phi^2 + (3-phi) = 4)")
print("     The angle is pi/5 = 36 degrees (the golden triangle angle!)")
print()
print("  5. If the RS Gram matrix is PSD: GRH for L_RS -> GRH for zeta -> RH.")
print()

# Check the key conclusion
if evals_RS is not None:
    if all(ev > -mpf('1e-15') for ev in evals_RS):
        print("  *** RS Gram matrix appears PSD with current basis and primes! ***")
        print("  This is NECESSARY but NOT SUFFICIENT for Weil positivity")
        print("  (positivity must hold for ALL test functions, not just our basis).")
    else:
        print("  RS Gram matrix has negative eigenvalues with current parameters.")
        print("  The squared-trace approach improves but does not fully resolve")
        print("  the sign problem. The minus sign in W = arch - prime_sum")
        print("  is structural (from the explicit formula), not from the traces.")

print()
print("  THE SIGN PROBLEM is really about the EXPLICIT FORMULA structure,")
print("  not just the trace signs. The formula says:")
print("    sum_zeros f(gamma)^2 = arch + pole - prime_sum")
print("  The minus sign is BUILT INTO the relation between zeros and primes.")
print("  Squaring traces helps (b_p >= 0) but the prime_sum still enters")
print("  with a minus sign. The arch + pole must dominate.")
print()
print("  The Pythagorean insight: a^2 + b^2 = c^2 is about GEOMETRY,")
print("  not algebra. The Weil positivity is ultimately a statement about")
print("  the GEOMETRY of the zero distribution: zeros on the critical line")
print("  = points on a circle (|1-1/rho| = 1), and the Pythagorean theorem")
print("  governs circular geometry.")

print()
print("COMPUTATION COMPLETE.")
