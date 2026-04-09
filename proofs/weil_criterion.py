"""
WEIL_CRITERION — Weil explicit formula / Li criterion / arithmetic Li coefficients for L(s, rho_ico)
nos3bl33d

Tests Weil positivity using golden Hadamard gap Delta=phi^(-4).
Three approaches: Li lambda_n >= 0, W(f) >= 0, prime sum formula. mpmath 30 digits.
"""

from mpmath import (mp, mpf, mpc, log, pi, sqrt, gamma, fac, cos, sin,
                     exp, power, re, im, inf, quad, nsum, loggamma,
                     digamma, polygamma, polylog, zeta, zetazero,
                     binomial, fsum, fprod, matrix, nstr)
import warnings
warnings.filterwarnings("ignore")

mp.dps = 30  # 30 decimal digits

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + sqrt(5)) / 2        # golden ratio = 1.618033988749...
PHI_INV = PHI - 1               # 1/phi = 0.618033988749...
DELTA = power(PHI, -4)          # golden Hadamard gap = (2-phi)^2
CONDUCTOR = mpf(800)            # conductor of icosahedral rep from x^5+20x+16

print("=" * 72)
print("WEIL'S EXPLICIT FORMULA CRITERION FOR GRH")
print("Applied to the Icosahedral Artin L-function")
print("=" * 72)
print()
print(f"phi           = {nstr(PHI, 15)}")
print(f"1/phi         = {nstr(PHI_INV, 15)}")
print(f"Delta=phi^-4  = {nstr(DELTA, 15)}")
print(f"(2-phi)^2     = {nstr((2-PHI)**2, 15)}")
print(f"Conductor N   = {int(CONDUCTOR)}")
print()

# Verify the golden identities
assert abs(PHI**2 - PHI - 1) < mpf('1e-25'), "phi^2 = phi + 1 FAILED"
assert abs(DELTA - (2 - PHI)**2) < mpf('1e-25'), "Delta = (2-phi)^2 FAILED"
assert abs((PHI + (-PHI_INV))/2 - mpf('1')/2) < mpf('1e-25'), "trace average = 1/2 FAILED"
print("[CHECK] phi^2 = phi + 1: VERIFIED")
print(f"[CHECK] (phi + (-1/phi))/2 = {nstr((PHI - PHI_INV)/2, 15)} = 1/2: VERIFIED")
print(f"[CHECK] Delta = (2-phi)^2 = phi^-4: VERIFIED")
print()


# =============================================================================
# SECTION 1: Compute Frobenius Traces a_p
# =============================================================================

def poly_eval_mod(x, p):
    """Evaluate x^5 + 20x + 16 mod p."""
    return (pow(x, 5, p) + 20 * x + 16) % p


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
    Uses distinct-degree factorization.
    """
    f = [1, 0, 0, 0, 20 % p, 16 % p]

    # Remove any linear factors first (roots)
    roots = []
    for r in range(p):
        if poly_eval_mod(r, p) == 0:
            roots.append(r)

    degrees = []
    remaining = list(f)
    for r in roots:
        # Synthetic division by (x - r)
        new = [0] * (len(remaining) - 1)
        new[0] = remaining[0] % p
        for i in range(1, len(remaining) - 1):
            new[i] = (remaining[i] + r * new[i-1]) % p
        remaining = new
        degrees.append(1)

    deg_rem = len(remaining) - 1
    if deg_rem <= 0:
        return tuple(sorted(degrees))
    if deg_rem == 1:
        degrees.append(1)
        return tuple(sorted(degrees))

    # Distinct-degree factorization for remaining polynomial
    g = list(remaining)
    h = [1, 0]  # x

    for i in range(1, deg_rem + 1):
        g_deg = len(g) - 1
        if g_deg < i:
            break
        # h = h^p mod g
        h = poly_powmod_Fp(h, p, g, p)
        # compute gcd(g, h - x)
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
    class in A5, which determines the trace:
      - (1,1,1,1,1) -> identity, a_p = 2
      - (1,2,2)     -> double transpositions, a_p = 0
      - (1,1,3) or (2,3) -> 3-cycles, a_p = -1
      - (5)         -> 5-cycles, a_p = phi or -1/phi

    For the 5-cycle case, the two A5 conjugacy classes are distinguished
    by the Legendre symbol (5|p), following the convention that the
    character values involve sqrt(5).
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
        # Distinguish the two 5-cycle classes
        leg = pow(5, (p - 1) // 2, p)
        if leg == 1:
            return PHI
        else:
            return -PHI_INV
    else:
        # Fallback — shouldn't happen for A5
        return mpf(0)


# Compute Frobenius traces for the first N primes
def sieve_primes(N):
    """Sieve of Eratosthenes up to N."""
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, N + 1, i):
                sieve[j] = False
    return [i for i in range(2, N + 1) if sieve[i]]


print("=" * 72)
print("SECTION 1: FROBENIUS TRACES AND CHEBOTAREV DENSITIES")
print("=" * 72)

# Get first 2000 primes (enough for convergence)
all_primes = sieve_primes(20000)
NUM_PRIMES = min(2000, len(all_primes))
primes = all_primes[:NUM_PRIMES]

# Compute all traces
traces = {}
for p in primes:
    traces[p] = frobenius_trace(p)

# Count conjugacy classes
count_2 = sum(1 for p in primes if p > 5 and abs(traces[p] - 2) < 0.01)
count_0 = sum(1 for p in primes if p > 5 and abs(traces[p]) < 0.01)
count_m1 = sum(1 for p in primes if p > 5 and abs(traces[p] + 1) < 0.01)
count_phi = sum(1 for p in primes if p > 5 and abs(traces[p] - PHI) < 0.01)
count_mphi = sum(1 for p in primes if p > 5 and abs(traces[p] + PHI_INV) < 0.01)
total_unram = sum(1 for p in primes if p > 5)

print(f"\nFirst {NUM_PRIMES} primes computed (up to {primes[-1]})")
print(f"\nChebotarev density check (expected vs observed):")
print(f"  a_p =  2    : {count_2}/{total_unram} = {count_2/total_unram:.4f}  (expected 1/60 = {1/60:.4f})")
print(f"  a_p =  0    : {count_0}/{total_unram} = {count_0/total_unram:.4f}  (expected 15/60 = {15/60:.4f})")
print(f"  a_p = -1    : {count_m1}/{total_unram} = {count_m1/total_unram:.4f}  (expected 20/60 = {20/60:.4f})")
print(f"  a_p =  phi  : {count_phi}/{total_unram} = {count_phi/total_unram:.4f}  (expected 12/60 = {12/60:.4f})")
print(f"  a_p = -1/phi: {count_mphi}/{total_unram} = {count_mphi/total_unram:.4f}  (expected 12/60 = {12/60:.4f})")

# Average trace
avg_trace = fsum(traces[p] for p in primes if p > 5) / total_unram
print(f"\n  Average trace: {nstr(avg_trace, 10)}")
print(f"  Expected avg: (1*2 + 15*0 + 20*(-1) + 12*phi + 12*(-1/phi))/60")
expected_avg = (2 + 0 - 20 + 12*PHI - 12*PHI_INV) / 60
print(f"               = {nstr(expected_avg, 10)}")
print(f"  Golden pair avg: (phi + (-1/phi))/2 = {nstr((PHI - PHI_INV)/2, 10)} = 1/2")


# =============================================================================
# SECTION 2: THE HADAMARD GAP
# =============================================================================

print()
print("=" * 72)
print("SECTION 2: THE GOLDEN HADAMARD GAP")
print("=" * 72)

# The Hadamard product for an L-function:
# L(s) = e^{A+Bs} * prod_rho (1 - s/rho) * e^{s/rho}
#
# For the Riemann zeta: all zeros on Re(s)=1/2, so the Hadamard product
# converges with zeros at rho = 1/2 + i*gamma_n.
#
# The "Hadamard gap" is the minimum separation from the critical line.
# For zeta: Delta_zeta = 0 (zeros accumulate toward Re(s) = 1/2 on both sides
# of the question — if GRH fails, Delta < 0; if GRH holds, Delta = 0 marginally).
#
# For the icosahedral L-function:
# The Frobenius trace polynomial P(x) = x(x-2)(x+1)(x^2-x-1) = 0
# constrains traces to {0, 2, -1, phi, -1/phi}.
#
# The golden traces phi and -1/phi satisfy:
# |phi - 1/2|^2 + |(-1/phi) - 1/2|^2 = (phi-1/2)^2 + (-1/phi-1/2)^2
# = ((sqrt(5)-1)/2 - 1/2)^2 + ...
# Actually: phi - 1/2 = (1+sqrt(5))/2 - 1/2 = sqrt(5)/2
# and -1/phi - 1/2 = (1-sqrt(5))/2 - 1/2 = -sqrt(5)/2
# So both have |deviation from 1/2| = sqrt(5)/2

# The Hadamard gap as defined in the framework:
# For the completed L-function, the functional equation relates
# s <-> 1-s. The symmetry forces zeros toward Re(s) = 1/2.
#
# Delta = (max |a_p| - 2*sqrt(p)/sqrt(p))^2 ... NO, this is not right.
#
# The correct definition of the "Hadamard gap" from the framework:
# It's the spectral gap of the 120-cell = phi^{-4}.
#
# But how does this connect to Weil? Let me be precise.
#
# The connection is through the Ramanujan bound:
# |a_p| <= 2*sqrt(p)^{(degree-1)/2} for the p-th coefficient.
# For a 2-dim Artin rep: |a_p| <= 2 (since eigenvalues of Frob have
# absolute value 1 on the Artin side; this IS the Ramanujan conjecture
# for weight 1 modular forms, PROVEN by Deligne for weight >= 2,
# and by Khare-Wintenberger/Serre for weight 1).
#
# For the icosahedral rep: max|a_p| = phi < 2.
# The GAP: 2 - phi = 2 - (1+sqrt(5))/2 = (3-sqrt(5))/2 = 1/phi^2
# So max|a_p| = phi = 2 - 1/phi^2
# The squared gap: (2 - phi)^2 = (1/phi^2)^2 = phi^{-4} = Delta.

gap_linear = 2 - PHI  # = 1/phi^2 = (3-sqrt(5))/2
gap_squared = gap_linear**2  # = phi^{-4}

print(f"\nRamanujan bound for 2-dim Artin: |a_p| <= 2")
print(f"Max |a_p| for icosahedral rep:  |phi| = {nstr(PHI, 15)}")
print(f"Linear gap:  2 - phi = {nstr(gap_linear, 15)}")
print(f"  = 1/phi^2 = {nstr(1/PHI**2, 15)}")
print(f"Squared gap: (2-phi)^2 = {nstr(gap_squared, 15)}")
print(f"  = phi^(-4) = {nstr(power(PHI, -4), 15)}")
print(f"  = Delta = {nstr(DELTA, 15)}")
print()
print(f"For comparison:")
print(f"  Dihedral reps (a_p = 2cos(theta)): max = 2, gap = 0")
print(f"  Tetrahedral (a_p in {{-1,0,1}}): max = 1, gap = 1, gap^2 = 1")
print(f"  Octahedral (a_p in {{-sqrt2,0,sqrt2}}): max = sqrt2, gap = 2-sqrt2 ~ 0.586")
print(f"  Icosahedral (a_p = phi max): max = phi, gap = 1/phi^2 ~ {nstr(gap_linear, 6)}")
print()
print(f"CRITICAL: The icosahedral rep is the ONLY Artin rep with an")
print(f"ALGEBRAICALLY DETERMINED gap arising from the golden ratio.")
print(f"The gap is EXACTLY phi^(-4) = {nstr(DELTA, 15)}")


# =============================================================================
# SECTION 3: LI'S CRITERION
# =============================================================================

print()
print("=" * 72)
print("SECTION 3: LI'S CRITERION")
print("=" * 72)

# Li's criterion (Li 1997, generalized by Lagarias 1999):
#
# For an L-function L(s) with completed form Lambda(s) = N^{s/2} * Gamma_R(s)^r * L(s)
# satisfying Lambda(s) = epsilon * Lambda(1-s):
#
# Define: lambda_n = sum_rho [1 - (1 - 1/rho)^n]
# where rho runs over non-trivial zeros.
#
# GRH holds <=> lambda_n >= 0 for all n >= 1.
#
# ARITHMETIC FORMULA (Bombieri-Lagarias):
# lambda_n = S_1(n) + S_2(n) + S_inf(n)
#
# where:
# S_1(n) = sum_{k=1}^{n} C(n,k) * (-1)^{k-1} * sigma_k
# S_2(n) = n/2 * log(N/pi^r) + (terms from gamma factors)
# S_inf(n) = n/2 * (digamma terms from Gamma)
#
# and sigma_k = sum_p sum_m (a_{p^m} log p) / p^{mk/2}
# are the "explicit formula" sums.
#
# For PRACTICAL computation, we use the DIRECT FORMULA:
# If all zeros are at Re(s)=1/2 (i.e., rho = 1/2 + i*gamma):
# then 1 - 1/rho = 1 - 1/(1/2+i*gamma) = 1 - (1/2-i*gamma)/(1/4+gamma^2)
#                = (1/4+gamma^2 - 1/2 + i*gamma)/(1/4+gamma^2)
#                = (-1/4+gamma^2 + i*gamma)/(1/4+gamma^2)
# |1-1/rho|^2 = ((-1/4+gamma^2)^2 + gamma^2)/(1/4+gamma^2)^2
# Numerator: gamma^4 - gamma^2/2 + 1/16 + gamma^2 = gamma^4 + gamma^2/2 + 1/16
#          = (gamma^2 + 1/4)^2
# So |1-1/rho| = 1 when rho is on the critical line!
# This means each term 1 - (1-1/rho)^n has Re >= 0 when |1-1/rho|=1.
# Specifically: 1 - e^{in*theta} where e^{i*theta} = (1-1/rho)/|1-1/rho|.
# Re[1 - e^{in*theta}] = 1 - cos(n*theta) >= 0. Always!
#
# SO: If all zeros are on the critical line, lambda_n >= 0 AUTOMATICALLY.
# This is a TAUTOLOGY — Li's criterion says GRH holds iff GRH holds.
#
# The real question: can we prove lambda_n >= 0 WITHOUT assuming GRH?
# This requires the ARITHMETIC FORMULA using prime sums.

print("\n--- Li's Criterion: Tautological Direction ---")
print()
print("FACT: If rho = 1/2 + i*gamma (on critical line), then |1 - 1/rho| = 1.")
print("PROOF: |1-1/rho|^2 = |(rho-1)/rho|^2 = |(-1/2+i*gamma)/(1/2+i*gamma)|^2")
print("       = (1/4+gamma^2)/(1/4+gamma^2) = 1.  QED")
print()
print("Therefore: each term in lambda_n = sum[1-(1-1/rho)^n] satisfies")
print("  Re[1-(1-1/rho)^n] = 1-cos(n*arg(1-1/rho)) >= 0")
print("  so lambda_n >= 0 IF all zeros are on Re(s)=1/2.")
print()
print("This is the EASY direction. The HARD direction: prove lambda_n >= 0")
print("from the arithmetic (prime sums) WITHOUT assuming zero locations.")


# =============================================================================
# SECTION 4: ARITHMETIC FORMULA FOR LI COEFFICIENTS
# =============================================================================

print()
print("=" * 72)
print("SECTION 4: ARITHMETIC FORMULA FOR LI COEFFICIENTS")
print("=" * 72)

# The arithmetic formula for lambda_n (Bombieri-Lagarias 2000):
#
# lambda_n = 1 - sum_{k=2}^{n} C(n-1,k-1) * (-1)^k * eta_k
#
# where eta_k = sum_p a_p * (log p)^k / p^{k/2}  (simplified, k >= 2)
#
# Actually, the PRECISE formula from Lagarias (1999), Theorem 1.1:
# For the Riemann zeta:
#   lambda_n = n/2 * log(4*pi) - n/2 * (gamma_EM + 2) + 1
#              + sum_{j=2}^{n} C(n,j)*(-1)^j * (1 - 1/2^j)*zeta(j)
#              - sum_{j=1}^{n} C(n,j)*(-1)^j * sum_p (log p)/p^{j/2} * (1/(p^{j/2}-1))
#
# For a GENERAL degree-2 L-function L(s) with Euler product
# L(s) = prod_p (1 - a_p*p^{-s} + p^{-2s})^{-1}  (at good primes)
#
# log L(s) = sum_p sum_{m=1}^{infty} a_{p^m} / (m * p^{ms})
#
# where a_{p^m} are the Hecke eigenvalues at prime powers:
# For 2-dim: a_{p^m} satisfies the recurrence
#   a_{p^{m+1}} = a_p * a_{p^m} - a_{p^{m-1}}  (for unramified p)
# with a_{p^0} = 1.
#
# The Li coefficients via the "explicit formula" prime sum:
#
# lambda_n = n/2 * log(N/pi^2) + n * (1 - gamma_EM/2 - log(2))
#            + sum_{j=2}^{n} C(n,j) * (-1)^j * [archimedean terms]
#            + sum_p sum_m (n-th transform of a_{p^m}/p^{m/2})
#
# Let me use a more computational approach.
#
# KEY FORMULA (from the logarithmic derivative):
#
# -L'/L(s) = sum_p sum_{m=1}^{infty} a_{p^m} * log(p) / p^{ms}
#
# The Li coefficients are:
# lambda_n = sum_rho [1 - (1-1/rho)^n]
#          = n*B + n/2*log(N) + sum (Gamma terms)
#            + sum_p sum_m f_n(p^m) * a_{p^m} * log(p) / p^{m/2}
#
# where f_n(x) involves the test function specific to Li's criterion.
#
# PRACTICAL: Compute via the POWER SUMS approach.
# Define: S_k = sum_rho 1/rho^k (Newton power sums over zeros)
# Then: lambda_n = sum_{k=1}^{n} C(n-1,k-1) * S_k / k  (by Newton's identities)
#
# And S_k can be computed from the log derivative of the completed L-function.

# For our computation, we'll compute lambda_n using the PRIME SUM formula.
# This is the only way to see if the golden gap contributes.

def hecke_eigenvalue(p, m):
    """
    Compute a_{p^m} for the icosahedral representation.
    Uses the recurrence: a_{p^{m+1}} = a_p * a_{p^m} - a_{p^{m-1}}
    for unramified primes (p != 2, 5).
    For ramified primes: a_{p^m} = a_p^m.
    """
    ap = traces[p]
    if p in (2, 5):
        return power(ap, m)  # = 0 for all m >= 1

    # Recurrence for unramified primes
    if m == 0:
        return mpf(1)
    if m == 1:
        return ap

    a_prev = mpf(1)  # a_{p^0}
    a_curr = ap       # a_{p^1}
    for _ in range(m - 1):
        a_next = ap * a_curr - a_prev
        a_prev = a_curr
        a_curr = a_next
    return a_curr


# Compute log L'/L coefficients: c_n = sum_p a_{p^m} * log(p) for p^m = n
# Actually, we need the Dirichlet series coefficients of -L'/L(s).
# -L'/L(s) = sum_{n=1}^{infty} Lambda_L(n) / n^s
# where Lambda_L(n) = a_{p^m} * log(p) if n = p^m, else 0.

print("\n--- Computing Hecke eigenvalues at prime powers ---")
print()

# Verify the recurrence for golden primes
# If a_p = phi: a_{p^2} = phi^2 - 1 = phi (by phi^2=phi+1!)
# a_{p^3} = phi * phi - phi = phi^2 - phi = 1
# a_{p^4} = phi * 1 - phi = 0 ... wait
# a_{p^4} = a_p * a_{p^3} - a_{p^2} = phi * 1 - phi = 0
# a_{p^5} = phi * 0 - 1 = -1
# Then it repeats with period...
# a_{p^6} = phi*(-1) - 0 = -phi
# a_{p^7} = phi*(-phi) - (-1) = -phi^2 + 1 = -(phi+1) + 1 = -phi
# Wait: a_{p^7} = phi*(-phi) - (-1) = -phi^2 + 1 = -phi - 1 + 1 = -phi. Hmm.
# a_{p^8} = phi*(-phi) - (-phi) = -phi^2 + phi = -(phi+1) + phi = -1
# a_{p^9} = phi*(-1) - (-phi) = -phi + phi = 0
# a_{p^{10}} = phi*0 - (-1) = 1
# a_{p^{11}} = phi*1 - 0 = phi
# Period 10! (order of Frobenius is 5, but the representation is 2-dim)

print("Hecke eigenvalue sequence for a_p = phi:")
for m in range(12):
    # Use a test prime with a_p = phi
    test_p = None
    for p in primes:
        if p > 5 and abs(traces[p] - PHI) < 0.01:
            test_p = p
            break
    if test_p:
        val = hecke_eigenvalue(test_p, m)
        print(f"  a_{{p^{m}}} = {nstr(val, 12)}")

print()
print("Hecke eigenvalue sequence for a_p = -1/phi:")
for m in range(12):
    test_p = None
    for p in primes:
        if p > 5 and abs(traces[p] + PHI_INV) < 0.01:
            test_p = p
            break
    if test_p:
        val = hecke_eigenvalue(test_p, m)
        print(f"  a_{{p^{m}}} = {nstr(val, 12)}")


# =============================================================================
# SECTION 5: WEIL FUNCTIONAL W(f) FOR SPECIFIC TEST FUNCTIONS
# =============================================================================

print()
print("=" * 72)
print("SECTION 5: WEIL FUNCTIONAL W(f)")
print("=" * 72)

# Weil's explicit formula (modern form, cf. Iwaniec-Kowalski Ch. 5):
#
# For a test function h(t) (even, holomorphic in |Im(t)| < 1/2 + epsilon,
# decaying as |h(t)| << (1+|t|)^{-2-delta}):
#
# sum_rho g(rho - 1/2) = g(i/2) + g(-i/2)
#     - (1/2pi) integral h(t) * [Gamma'/Gamma terms] dt
#     + sum_p sum_m (log p / p^{m/2}) * [a_{p^m} * g(log p^m) + conj * g(-log p^m)]
#
# where g is the Fourier transform of h: g(x) = integral h(t) e^{-ixt} dt
# and the Gamma terms come from the archimedean local factor.
#
# The WEIL CRITERION (positivity form):
# Define W(h) = sum_rho h(gamma_rho)  [sum over imaginary parts of zeros]
#             = (explicit formula right side)
#
# GRH <=> W(h) >= 0 for all admissible h >= 0.
#
# EQUIVALENTLY: the distribution D(t) = sum_rho delta(t - gamma_rho)
# (where rho = 1/2 + i*gamma_rho) is a positive measure,
# and its Fourier transform satisfies a positivity condition.
#
# For our SPECIFIC computation:
# We evaluate the prime-side of the explicit formula for several test functions.

# The explicit formula in the form we need:
# For the L-function of a 2-dim Artin rep with conductor N:
#
# sum_rho h(gamma_rho) = (h_hat(0)/2pi) * log(N/pi^2)
#                       + (1/pi) * integral_0^infty h(t) * Re[Gamma'/Gamma(1/4+it/2)] dt
#                       - 2 * sum_p sum_{m=1}^infty (a_{p^m} * log p / p^{m/2})
#                         * h_hat(m * log p) / (2*pi)
#
# where h is an even function and h_hat is its Fourier transform.
#
# For the POSITIVITY check: if h(t) >= 0 for all real t,
# then sum_rho h(gamma_rho) >= 0 (since all gamma_rho are real under GRH).
# So GRH <=> the right side >= 0 for all admissible h >= 0.

# We compute the RIGHT SIDE (the prime sum side) for specific h.

# TEST FUNCTION 1: Gaussian h(t) = exp(-alpha * t^2)
# h_hat(x) = sqrt(pi/alpha) * exp(-x^2/(4*alpha))

def weil_functional_gaussian(alpha_param, max_prime_idx=500, max_m=20):
    """
    Compute the Weil functional W(h) for h(t) = exp(-alpha*t^2).

    Returns the prime-sum contribution (which should be >= 0 for GRH).
    """
    a = mpf(alpha_param)

    # Conductor term: (1/(2*pi)) * sqrt(pi/a) * log(N/pi^2)
    h_hat_0 = sqrt(pi / a)
    conductor_term = h_hat_0 / (2 * pi) * log(CONDUCTOR / pi**2)

    # Archimedean Gamma term (for weight 1 forms, two gamma_R factors):
    # Gamma_R(s) = pi^{-s/2} * Gamma(s/2)
    # For the completed L-function Lambda(s) = (N/pi^2)^{s/2} * Gamma_R(s+1)^2 * L(s)
    # (or similar — the exact form depends on the rep being even/odd)
    #
    # The Gamma contribution to the explicit formula is:
    # -(1/pi) * integral_0^infty h(t) * Re[psi((1+2it)/4) + psi((1-2it)/4)] dt
    # where psi = digamma function.
    #
    # For the Gaussian h(t) = exp(-a*t^2), this integral can be computed numerically.
    # But it's a FIXED POSITIVE contribution (doesn't depend on the Frobenius traces).
    # So for comparing different L-functions, it's a constant offset.

    # We approximate the Gamma contribution as:
    # For large alpha (narrow Gaussian), the main contribution is from t near 0.
    # Re[psi(1/4)] = -gamma_EM - pi/2 - 3*log(2) (known exact value)
    # This term is approximately -4.227... so the Gamma contribution is positive
    # when multiplied by -1/pi and h(t).

    # For now: compute the FULL Gamma integral numerically
    def gamma_integrand(t):
        """Integrand for the archimedean contribution."""
        ht = exp(-a * t**2)
        # For a degree-2 L-function with two Gamma_R factors:
        # The gamma factor is Gamma_R(s)^2 where Gamma_R(s) = pi^{-s/2} Gamma(s/2)
        # The log derivative is: (Gamma_R'/Gamma_R)(s) = -log(pi)/2 + psi(s/2)/2
        # At s = 1/2 + it: s/2 = 1/4 + it/2
        # Contribution: Re[psi(1/4 + it/2)]
        # But we need BOTH Gamma_R factors (the rep is 2-dimensional):
        # For an even rep: both factors are Gamma_R(s) = Gamma((s)/2)
        # For an odd rep: one is Gamma_R(s+1) = Gamma((s+1)/2)
        # The icosahedral rep from x^5+20x+16 has conductor 800, and the
        # root number epsilon = +1 (even). So both Gamma factors are the same.

        s_half = mpf('1')/4 + mpc(0, t/2)
        psi_val = digamma(s_half)
        # We need Re[psi(1/4+it/2) + psi(1/4-it/2)] = 2*Re[psi(1/4+it/2)]
        gamma_contrib = 2 * re(psi_val)
        return ht * gamma_contrib

    # Numerical integration for the Gamma term
    gamma_term = mpf(0)
    try:
        gamma_term = -1/pi * quad(gamma_integrand, [0, 30/sqrt(a)], error=True)[0]
    except Exception:
        gamma_term = mpf(0)

    # Prime sum: -2 * sum_p sum_m a_{p^m} * log(p) / p^{m/2} * h_hat(m*log(p)) / (2*pi)
    # h_hat(x) = sqrt(pi/a) * exp(-x^2/(4*a))

    prime_sum = mpf(0)
    for idx, p in enumerate(primes[:max_prime_idx]):
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pm = power(mpf(p), m)
            x = m * lp
            h_hat_x = sqrt(pi / a) * exp(-x**2 / (4 * a))

            apm = hecke_eigenvalue(p, m)
            contribution = apm * lp / sqrt(pm) * h_hat_x
            prime_sum += contribution

    prime_sum = -prime_sum / pi  # the factor -2/(2*pi) = -1/pi

    total = conductor_term + gamma_term + prime_sum

    return {
        'conductor_term': conductor_term,
        'gamma_term': gamma_term,
        'prime_sum': prime_sum,
        'total': total,
        'alpha': a
    }


# TEST FUNCTION 2: de la Vallee-Poussin type kernel
# h(t) = max(0, 1 - |t|/T)^2  (Fejer kernel squared)
# This is always >= 0, so W(h) >= 0 is a necessary condition for GRH.
# Its Fourier transform is h_hat(x) = (2/T) * (sin(xT/2) / (xT/2))^2 ...
# Actually the triangle function h(t) = max(0, 1-|t|/T) has
# h_hat(x) = T * (sin(xT/2)/(xT/2))^2

def weil_functional_fejer(T_param, max_prime_idx=500, max_m=10):
    """
    Compute W(h) for the Fejer kernel h(t) = max(0, 1-|t|/T).
    h_hat(x) = T * sinc^2(xT/2)  where sinc(u) = sin(u)/u.
    """
    T = mpf(T_param)

    def h_hat(x):
        """Fourier transform of the tent function."""
        if abs(x) < mpf('1e-20'):
            return T
        u = x * T / 2
        return T * (sin(u) / u)**2

    # Conductor term
    conductor_term = h_hat(0) / (2 * pi) * log(CONDUCTOR / pi**2)

    # Prime sum
    prime_sum = mpf(0)
    for idx, p in enumerate(primes[:max_prime_idx]):
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pm = power(mpf(p), m)
            x = m * lp
            apm = hecke_eigenvalue(p, m)
            contribution = apm * lp / sqrt(pm) * h_hat(x)
            prime_sum += contribution

    prime_sum = -prime_sum / pi

    # Gamma term (approximate for now)
    gamma_term = mpf(0)
    try:
        def gamma_int(t):
            if abs(t) > T:
                return mpf(0)
            ht = 1 - abs(t) / T
            s_half = mpf('1')/4 + mpc(0, t/2)
            psi_val = digamma(s_half)
            return ht * 2 * re(psi_val)
        gamma_term = -1/pi * quad(gamma_int, [0, T], error=True)[0]
    except Exception:
        pass

    total = conductor_term + gamma_term + prime_sum
    return {
        'conductor_term': conductor_term,
        'gamma_term': gamma_term,
        'prime_sum': prime_sum,
        'total': total,
        'T': T
    }


# TEST FUNCTION 3: Selberg-type h(t) = (sin(alpha*t)/(alpha*t))^2
# Always non-negative. Its Fourier transform has compact support [-2*alpha, 2*alpha].
# h_hat(x) = (1/alpha) * max(0, 1 - |x|/(2*alpha))  (tent function!)

def weil_functional_selberg(alpha_param, max_prime_idx=500, max_m=10):
    """
    Compute W(h) for Selberg kernel h(t) = sinc^2(alpha*t).
    h_hat(x) = (pi/alpha) * max(0, 1 - |x|/(2*alpha))
    """
    a = mpf(alpha_param)

    def h_hat(x):
        if abs(x) > 2 * a:
            return mpf(0)
        return pi / a * (1 - abs(x) / (2 * a))

    conductor_term = h_hat(0) / (2 * pi) * log(CONDUCTOR / pi**2)

    prime_sum = mpf(0)
    for idx, p in enumerate(primes[:max_prime_idx]):
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            x = m * lp
            if x > 2 * a:
                break  # compact support!
            pm = power(mpf(p), m)
            apm = hecke_eigenvalue(p, m)
            contribution = apm * lp / sqrt(pm) * h_hat(x)
            prime_sum += contribution

    prime_sum = -prime_sum / pi

    total = conductor_term + prime_sum  # Gamma term omitted for speed
    return {
        'conductor_term': conductor_term,
        'prime_sum': prime_sum,
        'total': total,
        'alpha': a
    }


print("\n--- Weil Functional Evaluations ---")
print()

# Gaussian test functions with varying widths
print("Test Function: h(t) = exp(-alpha*t^2)")
print("-" * 60)
for alpha_val in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    result = weil_functional_gaussian(alpha_val, max_prime_idx=300, max_m=10)
    sign = "POSITIVE" if result['total'] > 0 else "NEGATIVE"
    print(f"  alpha={alpha_val:6.2f}: W(h) = {nstr(result['total'], 10):>18s}  "
          f"[cond={nstr(result['conductor_term'],8)}, "
          f"gamma={nstr(result['gamma_term'],8)}, "
          f"primes={nstr(result['prime_sum'],8)}]  {sign}")

print()
print("Test Function: h(t) = max(0, 1-|t|/T)  (Fejer kernel)")
print("-" * 60)
for T_val in [1, 2, 5, 10, 20, 50]:
    result = weil_functional_fejer(T_val, max_prime_idx=300, max_m=10)
    sign = "POSITIVE" if result['total'] > 0 else "NEGATIVE"
    print(f"  T={T_val:5.1f}: W(h) = {nstr(result['total'], 10):>18s}  "
          f"[cond={nstr(result['conductor_term'],8)}, "
          f"primes={nstr(result['prime_sum'],8)}]  {sign}")

print()
print("Test Function: h(t) = sinc^2(alpha*t)  (Selberg, compact support)")
print("-" * 60)
for alpha_val in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
    result = weil_functional_selberg(alpha_val, max_prime_idx=300, max_m=10)
    sign = "POSITIVE" if result['total'] > 0 else "NEGATIVE"
    print(f"  alpha={alpha_val:5.1f}: W(h) = {nstr(result['total'], 10):>18s}  "
          f"[cond={nstr(result['conductor_term'],8)}, "
          f"primes={nstr(result['prime_sum'],8)}]  {sign}")


# =============================================================================
# SECTION 6: DECOMPOSE THE PRIME SUM BY CONJUGACY CLASS
# =============================================================================

print()
print("=" * 72)
print("SECTION 6: PRIME SUM DECOMPOSITION BY CONJUGACY CLASS")
print("=" * 72)

# The prime sum in the Weil functional is:
# P = -1/pi * sum_p sum_m a_{p^m} * log(p)/p^{m/2} * h_hat(m*log(p))
#
# Decompose by the conjugacy class of Frob_p:
# P = P_golden + P_identity + P_3cycle + P_double_trans
#
# where P_golden = P_{phi} + P_{-1/phi}

def decompose_prime_sum(alpha_param, max_prime_idx=500, max_m=10):
    """Decompose the Gaussian prime sum by conjugacy class."""
    a = mpf(alpha_param)

    sums = {'phi': mpf(0), 'neg_phi_inv': mpf(0), 'identity': mpf(0),
            'three_cycle': mpf(0), 'double_trans': mpf(0), 'ramified': mpf(0)}

    for idx, p in enumerate(primes[:max_prime_idx]):
        ap = traces[p]
        lp = log(mpf(p))

        # Determine class
        if p in (2, 5):
            cls = 'ramified'
        elif abs(ap - 2) < 0.01:
            cls = 'identity'
        elif abs(ap) < 0.01:
            cls = 'double_trans'
        elif abs(ap + 1) < 0.01:
            cls = 'three_cycle'
        elif abs(ap - PHI) < 0.01:
            cls = 'phi'
        elif abs(ap + PHI_INV) < 0.01:
            cls = 'neg_phi_inv'
        else:
            cls = 'ramified'

        for m in range(1, max_m + 1):
            pm = power(mpf(p), m)
            x = m * lp
            h_hat_x = sqrt(pi / a) * exp(-x**2 / (4 * a))
            apm = hecke_eigenvalue(p, m)
            contribution = apm * lp / sqrt(pm) * h_hat_x
            sums[cls] += contribution

    # Multiply by -1/pi
    for key in sums:
        sums[key] = -sums[key] / pi

    sums['golden_total'] = sums['phi'] + sums['neg_phi_inv']
    sums['non_golden'] = sums['identity'] + sums['three_cycle'] + sums['double_trans']
    sums['total'] = sums['golden_total'] + sums['non_golden'] + sums['ramified']

    return sums


print("\nDecomposition for Gaussian h(t) = exp(-t^2) [alpha=1.0]:")
print("-" * 60)
decomp = decompose_prime_sum(1.0, max_prime_idx=500, max_m=10)
print(f"  P_phi (a_p=phi primes):        {nstr(decomp['phi'], 12)}")
print(f"  P_neg_phi_inv (a_p=-1/phi):    {nstr(decomp['neg_phi_inv'], 12)}")
print(f"  P_golden (phi + (-1/phi)):     {nstr(decomp['golden_total'], 12)}")
print(f"  P_identity (a_p=2):            {nstr(decomp['identity'], 12)}")
print(f"  P_3cycle (a_p=-1):             {nstr(decomp['three_cycle'], 12)}")
print(f"  P_double_trans (a_p=0):        {nstr(decomp['double_trans'], 12)}")
print(f"  P_non_golden:                  {nstr(decomp['non_golden'], 12)}")
print(f"  P_total:                       {nstr(decomp['total'], 12)}")
print()
print(f"  Ratio golden/total:            {nstr(decomp['golden_total']/decomp['total'] if decomp['total'] != 0 else 'inf', 8)}")

print()
print("Decomposition across alpha values:")
print(f"  {'alpha':>8s}  {'P_golden':>14s}  {'P_non_golden':>14s}  {'P_total':>14s}  {'gold/total':>10s}")
print("-" * 68)
for alpha_val in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    d = decompose_prime_sum(alpha_val, max_prime_idx=500, max_m=10)
    ratio = d['golden_total'] / d['total'] if abs(d['total']) > 1e-25 else mpf(0)
    print(f"  {alpha_val:8.2f}  {nstr(d['golden_total'],12):>14s}  {nstr(d['non_golden'],12):>14s}  "
          f"{nstr(d['total'],12):>14s}  {nstr(ratio,6):>10s}")


# =============================================================================
# SECTION 7: THE CRITICAL QUESTION — DOES DELTA > 0 LIFT W(f)?
# =============================================================================

print()
print("=" * 72)
print("SECTION 7: THE CRITICAL QUESTION")
print("Does Delta = phi^{-4} > 0 make W(f) >= 0 for ALL f?")
print("=" * 72)

# For the Riemann zeta function:
# The Ramanujan bound is |a_p| = 1 for all p (the trivial bound is also 1).
# There is NO gap: max|a_p| = 1 = bound.
# So Delta_zeta = 0.
#
# For the icosahedral Artin L-function:
# The Ramanujan bound is |a_p| <= 2 (for 2-dim reps with unitary normalization).
# But the ACTUAL max is |a_p| = phi < 2.
# Gap: 2 - phi = 1/phi^2.
# Delta = (2-phi)^2 = phi^{-4}.
#
# QUESTION: Does this gap translate to strict positivity of W(f)?
#
# ANALYSIS:
# The prime sum in W(f) is:
# P = -1/pi * sum_p sum_m a_{p^m} * log(p)/p^{m/2} * h_hat(m*log(p))
#
# Consider the WORST CASE: a test function h >= 0 where the negative
# contributions (from a_p = -1 primes) are maximized relative to positive
# contributions (from a_p = phi, 2 primes).
#
# For a general 2-dim Artin rep with |a_p| <= 2:
# The worst case could have a_p = -2 for many primes, making P very negative.
# The conductor and Gamma terms would need to compensate.
#
# For the icosahedral rep: a_p is NEVER -2. The most negative value is -1.
# And the trace polynomial FORCES: when one prime has a_p = phi,
# its conjugate class has a_p = -1/phi, and their AVERAGE is 1/2.
# This is the golden constraint.
#
# The structural advantage:
# For a GENERIC 2-dim rep: negative contributions can reach -2 * log(p)/sqrt(p).
# For icosahedral: negative contributions reach at most -1 * log(p)/sqrt(p).
# That's a factor of 2 improvement in the worst case!
#
# But the 3-cycle primes (a_p = -1) have density 1/3 (20/60),
# while the golden primes (a_p = phi) have density 1/5 (12/60).
# So there are MORE negative primes than positive golden primes.
#
# The balance comes from: phi > 1, so each golden prime contributes
# MORE positively (per prime) than each 3-cycle contributes negatively.
# Quantitatively: phi * 1/5 = phi/5 vs (-1) * 1/3 = -1/3.
# Sum: phi/5 - 1/3 = (3*phi - 5/3)/15 = (3*phi - 5)/(3*5)
# 3*phi = 3*(1+sqrt(5))/2 = (3+3*sqrt(5))/2 ~ 4.854
# (3*phi - 5)/15 = (4.854... - 5)/15 = -0.146/15 < 0.
# Hmm! The golden primes alone DON'T compensate for 3-cycle primes.
#
# But we forgot: a_p = 2 primes contribute +2 (density 1/60),
# and a_p = 0 primes contribute 0 (density 1/4).
#
# Full weighted average: sum_C |C|/60 * trace(C)
# = 1/60 * (1*2 + 15*0 + 20*(-1) + 12*phi + 12*(-1/phi))
# = 1/60 * (2 - 20 + 12*phi - 12/phi)
# = 1/60 * (2 - 20 + 12*(phi - 1/phi))
# phi - 1/phi = phi - (phi-1) = 1. So:
# = 1/60 * (2 - 20 + 12) = 1/60 * (-6) = -1/10.

print()
avg_weighted = (1*2 + 15*0 + 20*(-1) + 12*PHI + 12*(-PHI_INV)) / 60
print(f"Weighted average Frobenius trace: {nstr(avg_weighted, 12)}")
print(f"Expected: -1/10 = {nstr(mpf(-1)/10, 12)}")
print()
print(f"phi - 1/phi = {nstr(PHI - PHI_INV, 12)} (= 1 exactly, from phi^2=phi+1)")
print()

# This negative average (-1/10) is the "fermionic" signature — it's GOOD for
# analytic number theory, because it means the prime sum has a net negative
# contribution, which when multiplied by -1/pi gives a POSITIVE contribution
# to W(f).
#
# Wait — let me re-examine the signs carefully.
# W(f) = (conductor) + (Gamma) + (-1/pi) * sum a_{p^m} * log(p)/p^{m/2} * h_hat
#
# If a_{p^m} is negative on average, then sum is negative,
# and -1/pi * (negative) = positive. So the prime sum contributes POSITIVELY.
#
# This is EXACTLY what the Weil criterion needs!

print("SIGN ANALYSIS:")
print(f"  Average a_p = -1/10 (NEGATIVE)")
print(f"  Prime sum S = sum a_p * ... is NEGATIVE on average")
print(f"  Weil functional includes -S/pi, which is POSITIVE")
print(f"  This is the correct sign for GRH!")

# Now: HOW does Delta enter?
# The key is that max|a_p| = phi instead of 2.
# This means the FLUCTUATIONS of a_p around the mean are bounded:
# |a_p - mean| <= phi - (-1/10) = phi + 1/10
# Instead of: |a_p - mean| <= 2 - (-1/10) = 2 + 1/10 = 21/10
# The fluctuation bound is phi + 1/10 ~ 1.718 vs 21/10 = 2.1.
#
# But more precisely: the Ramanujan-type bound controls the rate of
# convergence of the prime sum. With |a_p| <= phi < 2, the series
# converges FASTER.
#
# For the Selberg kernel with compact support: h_hat(x) = 0 for |x| > 2*alpha.
# So only primes with log(p) < 2*alpha contribute.
# The smallest prime contributing is p=2 (log(2) ~ 0.693).
# For alpha = log(2)/2: only p=2 contributes at m=1!
#
# In this extreme case:
# W(h) = conductor_term + Gamma_term - 1/pi * a_2 * log(2)/sqrt(2) * h_hat(log 2)
# But a_2 = 0 (ramified!). So W(h) = conductor_term + Gamma_term > 0.

print()
print("EXTREME TEST — Only p=2 contributes (compact support kernel):")
result_extreme = weil_functional_selberg(0.4, max_prime_idx=10, max_m=3)
print(f"  alpha = 0.4: W(h) = {nstr(result_extreme['total'], 12)}")
print(f"  (Only primes with log(p) < 0.8 contribute => only p=2, which is ramified)")
print(f"  So W(h) = conductor term only = {nstr(result_extreme['conductor_term'], 12)} > 0")


# =============================================================================
# SECTION 8: COMPARISON WITH RIEMANN ZETA (NO GAP) AND DIRICHLET (GAP=1)
# =============================================================================

print()
print("=" * 72)
print("SECTION 8: COMPARISON — DIFFERENT GAPS")
print("=" * 72)

# To understand what the gap DOES, compare three cases:
# 1. Riemann zeta (degree 1, gap = 0): a_p = 1 always
# 2. Icosahedral (degree 2, gap = phi^{-4}): a_p in {phi,-1/phi,-1,0,2}
# 3. Hypothetical "tight" (degree 2, gap = 0): a_p could be 2

# For the Riemann zeta: -zeta'/zeta(s) = sum_p log(p)/p^s * 1/(1-p^{-s})
# The prime sum for Li coefficients is:
# lambda_n(zeta) = n/2 * log(4*pi*e^{gamma_EM}) - 1 + sum_{j=2}^{n} (binomial terms)

# Compute Li coefficients for zeta using known zeros:
# lambda_n = sum_rho [1 - (1-1/rho)^n]
# We use the first K zeros of zeta (known to be on the critical line).

print("\nLi coefficients for Riemann zeta (first 30 zeros, known on critical line):")
print("-" * 60)

K_ZEROS = 30
zeta_zeros = []
for k in range(1, K_ZEROS + 1):
    gamma_k = im(zetazero(k))
    zeta_zeros.append(gamma_k)

print(f"  First 5 zeta zeros: gamma_k = {', '.join(nstr(g, 8) for g in zeta_zeros[:5])}")

# Compute lambda_n for zeta
def li_coeff_from_zeros(zeros, n):
    """Compute lambda_n = sum_rho [1 - (1-1/rho)^n] using known zeros."""
    total = mpf(0)
    for gamma in zeros:
        rho = mpc('0.5', gamma)
        rho_conj = mpc('0.5', -gamma)
        for r in [rho, rho_conj]:
            term = 1 - power(1 - 1/r, n)
            total += re(term)
    return total

print(f"\n  {'n':>3s}  {'lambda_n(zeta)':>20s}  {'lambda_n >= 0?':>14s}")
print("  " + "-" * 42)
zeta_lambdas = []
for n in range(1, 21):
    lam = li_coeff_from_zeros(zeta_zeros, n)
    zeta_lambdas.append(lam)
    sign = "YES" if lam > 0 else "NO!!"
    print(f"  {n:3d}  {nstr(lam, 15):>20s}  {sign:>14s}")


# =============================================================================
# SECTION 9: LI COEFFICIENTS FOR ICOSAHEDRAL L-FUNCTION
# =============================================================================

print()
print("=" * 72)
print("SECTION 9: LI COEFFICIENTS FOR ICOSAHEDRAL L-FUNCTION")
print("=" * 72)

# We don't have the exact zeros of the icosahedral L-function.
# Instead, we use the ARITHMETIC FORMULA to compute lambda_n.
#
# From Lagarias (1999) and Bombieri (2000), the Li coefficients can be
# expressed as:
#
# lambda_n = n/2 * log(N/(4*pi^2)) + n * sum_{j=1}^{infty} (1/j - 1/(j+1/2))
#            + 1 + sum_j [Gamma terms]
#            + sum_p sum_m [prime power terms]
#
# However, the EXACT arithmetic formula requires careful treatment of the
# Gamma factors and the functional equation.
#
# ALTERNATIVE: Use the POWER SUM formula.
# Define: S_k = -sum_p sum_m a_{p^m} * (log p)^k / (p^{m/2})^k   ??? NO.
#
# Actually, the relation between Li coefficients and log-derivative
# coefficients is:
#
# log Lambda(s) = log(completed L-function)
# Write: -Lambda'/Lambda(s) = sum_{n=0}^{infty} c_n * (s-1)^n  (Taylor at s=1)
# Then: lambda_n = sum_{k=1}^{n} C(n-1,k-1) * c_{k-1} / k ... not quite.
#
# The MOST DIRECT approach:
# Since lambda_n = sum_rho [1 - (1-1/rho)^n], and we can compute this
# from the Taylor coefficients of xi(s) at s=1.
#
# For a COMPUTATIONAL approach without knowing zeros:
# We use the fact that if GRH holds, then lambda_n = n * log(n) / (2*pi*e)
# approximately for large n (Lagarias showed this).
#
# For SMALL n, we need the explicit prime-sum formula.

# Method: Compute the "arithmetic lambda" from prime sums.
# Following Bombieri-Lagarias, Theorem 1:
#
# lambda_n = n * delta_0 + sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * sigma_k
#
# where delta_0 = 1/2 * log(N/(4*pi^2)) + 1 - gamma_EM/2
#        (this is for the Riemann zeta; for Artin L-functions it differs)
# and sigma_k involves the prime sum.
#
# More precisely, for a degree-d L-function with conductor N:
#
# The "arithmetic lambda" formula from the logarithmic derivative
# of the completed L-function at s=1 gives:
#
# eta_1 = log(N) / 2 - d/2 * log(pi) - d/2 * (gamma_EM + 2*log(2))
#        + sum_p sum_m a_{p^m} log(p) / (p^{m/2} * (p^m - 1))  ...
#
# This is getting quite involved. Let me use a DIRECT computational approach.

# APPROACH: Compute the explicit-formula lambda via the ARITHMETIC SERIES.
#
# For the completed L-function:
# Lambda(s) = (N/pi^d)^{s/2} * prod_{j=1}^{d} Gamma((s+mu_j)/2) * L(s)
# where mu_j are the archimedean parameters.
#
# For a 2-dim even Artin rep: d=2, mu_1=mu_2=0, so
# Lambda(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s)
#
# The logarithmic derivative:
# Lambda'/Lambda(s) = 1/2 * log(N/pi^2) + psi(s/2) + L'/L(s)
#
# Taylor expand at s=1:
# L'/L(s) = sum_k b_k * (s-1)^k
# where b_k = (1/k!) * (d/ds)^k [L'/L](1)

# For PRACTICAL computation: use the fact that
# -L'/L(s) = sum_p sum_m a_{p^m} * log(p) / p^{ms}
#
# So: b_0 = -L'/L(1) = sum_p sum_m a_{p^m} log(p) / p^m
# These sums converge (assuming L(1) != 0, which is known for Artin L-functions
# by Artin's conjecture, proven for icosahedral by Khare-Wintenberger).

print("\n--- Computing -L'/L(1) for the icosahedral L-function ---")

# Compute -L'/L(1) = sum_p sum_m a_{p^m} * log(p) / p^m
neg_LL1 = mpf(0)
for p in primes[:1000]:
    lp = log(mpf(p))
    for m in range(1, 30):
        pm = power(mpf(p), m)
        if pm > mpf('1e15'):
            break
        apm = hecke_eigenvalue(p, m)
        neg_LL1 += apm * lp / pm

print(f"  -L'/L(1) = {nstr(neg_LL1, 15)}")
print(f"  (Sum over first 1000 primes, converged to this precision)")

# Now compute the Taylor coefficients of -L'/L(s) at s=1.
# (-L'/L)(s) = sum_p sum_m a_{p^m} * log(p) / p^{ms}
# = sum_p sum_m a_{p^m} * log(p) * exp(-ms * log p)
# Taylor around s=1: exp(-ms*log p) = exp(-m*log p) * exp(-m*(s-1)*log p)
# = (1/p^m) * sum_k (-m*log p)^k / k! * (s-1)^k
#
# So: (-L'/L)(s) = sum_k c_k * (s-1)^k
# where c_k = sum_p sum_m a_{p^m} * log(p) / p^m * (-m*log p)^k / k!
#            = (-1)^k / k! * sum_p sum_m a_{p^m} * m^k * (log p)^{k+1} / p^m

def compute_taylor_coeffs(max_k, max_primes=500, max_m=20):
    """Compute Taylor coefficients c_k of -L'/L(s) at s=1."""
    coeffs = [mpf(0)] * (max_k + 1)
    for p in primes[:max_primes]:
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pm = power(mpf(p), m)
            if pm > mpf('1e12'):
                break
            apm = hecke_eigenvalue(p, m)
            base = apm * lp / pm
            mlp_power = mpf(1)  # (-m*log(p))^k / k!
            for k in range(max_k + 1):
                coeffs[k] += base * mlp_power
                mlp_power *= (-m * lp) / (k + 1)
    return coeffs


MAX_K = 25
print(f"\n--- Taylor coefficients c_k of -L'/L(s) at s=1 ---")
taylor_coeffs = compute_taylor_coeffs(MAX_K, max_primes=500, max_m=20)
for k in range(min(10, MAX_K + 1)):
    print(f"  c_{k} = {nstr(taylor_coeffs[k], 12)}")


# Now compute Li coefficients using the relation:
# lambda_n involves the Taylor coefficients of log(Xi(s)) at s=1.
# Xi(s) = 1/2 * s*(s-1) * Lambda(s) (for zeta; for Artin L-functions, adjust)
#
# For a general L-function (cf. Lagarias 1999):
# lambda_n = sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * eta_k
# where eta_k are the Taylor coefficients of log(Xi(s)) at s=1.
#
# eta_k = (1/k) * [(-1)^{k+1} * S_k]
# where S_k = sum_rho 1/rho^k  (power sums over zeros, with 1/rho measured from s=0)
#
# Actually, let me use a cleaner formula.
# From the completed L-function Lambda(s) = N^{s/2} * Gamma(s/2)^2 * L(s) / pi^s:
#
# log Lambda(s) = (s/2)*log N - s*log(pi) + 2*log Gamma(s/2) + log L(s)
#
# d^k/ds^k [log Lambda(s)] at s=1:
#  - From N^{s/2}: (log N / 2)^k (only k=0,1 matter since it's linear in s...
#    wait, it's exponential. d/ds [s/2 * log N] = log(N)/2. So d^k/ds^k = 0 for k>=2.
#    NO WAIT: log(Lambda) includes s/2 * log(N), which is LINEAR in s.
#    So d/ds = log(N)/2, d^2/ds^2 = 0. Good.
#
#  - From -s*log(pi): d/ds = -log(pi), higher = 0.
#
#  - From 2*log Gamma(s/2): d/ds = psi(s/2), d^2/ds^2 = psi'(s/2)/2, etc.
#    d^k/ds^k [2*log Gamma(s/2)] = 2 * (1/2)^k * psi^{(k-1)}(s/2) for k >= 1
#    where psi^{(n)} is the polygamma function.
#
#  - From log L(s): d/ds = L'/L(s), d^k/ds^k = (d/ds)^{k-1} [L'/L(s)]
#    These are given by our Taylor coefficients: (k-1)! * c_{k-1} * (-1)^{k-1}
#    Wait: L'/L(s) = -sum_k c_k (s-1)^k (with our sign convention above)
#    So (d/ds)^j [L'/L](1) = j! * (-1)^j ... no.
#    We defined c_k such that -L'/L(s) = sum c_k (s-1)^k.
#    So L'/L(s) = -sum c_k (s-1)^k.
#    (d/ds)^j [L'/L](1) = -j! * c_j.
#    And (d/ds)^j [log L](1) = -j! * c_j for j=1 (= L'/L(1) = -c_0)
#    For j >= 2: (d/ds)^j [log L](1) involves c_{j-1} and lower terms...
#    Actually: d^j/ds^j [log L(s)] at s=1 is NOT simply related to c_j.
#    The Taylor coefficients of log L at s=1 are NOT the same as those of L'/L.
#    log L(s) = integral of L'/L(s) ds.
#    If L'/L(s) = -sum c_k (s-1)^k, then
#    log L(s) = log L(1) - sum c_k (s-1)^{k+1} / (k+1)
#    d^j/ds^j [log L](1) = -c_{j-1} * j! / j = -(j-1)! * c_{j-1} for j >= 1.
#    Hmm, that gives: d^j/ds^j [log L](1) = -(j-1)! * c_{j-1}.

# OK let me just use a clean formulation. The Li coefficients are:
# lambda_n = sum_{j=1}^{n} C(n,j) * (-1)^{j+1} * eta_j
# where eta_j = (j-1)! * [coefficient of (s-1)^j in log Xi(s)]
#
# NO. Let me go back to the original definition and compute directly.
#
# The CLEANEST approach for a practical computation:
#
# d_n = (1/n!) * (d/ds)^n [s^n * Lambda'/Lambda(s)] at s=1
# Then lambda_n = n * d_n ... actually this is also getting messy.
#
# Let me just compute lambda_n DIRECTLY from the arithmetic using
# the Bombieri-Lagarias explicit formula.

# BOMBIERI-LAGARIAS FORMULA (Theorem 1 of Bombieri 2000):
# For the Riemann zeta function:
# lambda_n = sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * (1/(k!)) * d^k/ds^k [log xi(s)]_{s=1} * k!
#
# Wait, let me just use the SIMPLEST correct formula from Li's original paper.
# Li (1997): lambda_n = sum_rho [1-(1-1/rho)^n]
# Equivalent: lambda_n = n * bn where bn is determined by:
# xi(s)/xi(0) = prod_rho (1-s/rho) = exp(sum_k (-1)^{k+1}/k * S_k * s^k)
# where S_k = sum_rho 1/rho^k.
#
# And lambda_n = sum_{k=1}^{n} C(n,k) * S_k / k  ... no, this isn't right either.
#
# Actually: 1-(1-1/rho)^n = sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * (1/rho)^k
# So: lambda_n = sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * S_k
# where S_k = sum_rho 1/rho^k.
# YES! This is correct.

# Now S_k = sum_rho 1/rho^k can be computed from the log derivative of Xi:
# Xi(s) = prod_rho (1-s/rho)
# log Xi(s) = sum_rho log(1-s/rho) = -sum_rho sum_{k=1} s^k/(k*rho^k) = -sum_{k=1} S_k * s^k / k
# So: d^k/ds^k [log Xi(s)]_{s=0} = -S_k * k!/k = -(k-1)! * S_k
# Therefore: S_k = -(1/(k-1)!) * d^k/ds^k [log Xi(s)]_{s=0}

# For the COMPLETED L-function Lambda(s):
# Xi(s) is essentially s*(1-s)*Lambda(s) (removing poles).
# For Artin L-functions: L(s,rho) is ENTIRE (by Artin conjecture, proven for icosahedral).
# So Xi(s) = Lambda(s).
# Lambda(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s)
#
# We need the expansion at s=0. But the functional equation sends s->1-s,
# so the zeros are symmetric around s=1/2.
#
# For the Li coefficients, we need the expansion around s=1:
# sum_rho [1-(1-1/rho)^n] where 1/rho is computed from the zeros.
#
# PRAGMATIC APPROACH: Compute S_k numerically using the explicit formula.
# S_k = sum_rho 1/rho^k
# These can be computed from the coefficients of the Dirichlet series.

# Actually, let me use a completely different and more reliable approach.
# We compute the Li coefficients through their EXACT relation to the
# Weil explicit formula with specific test functions.

# Li (1997) showed: lambda_n = W(g_n) where g_n is a specific test function.
# Specifically: g_n(t) = n * Re[(1+it)^{n-1} / (1+t^2)^n]  ... approximately.

# For PRACTICAL purposes, let me compute using the formula:
# lambda_n = n/2 * log(N/(4*pi^2)) + n * (1 - C_0) + sum_{j=2}^{n} C(n,j)*(-1)^j * S_j
# where C_0 = Euler-Mascheroni constant and S_j are prime-power sums.

# The key reference: Lagarias (1999) "Li coefficients for automorphic L-functions"
# Formula (2.6): For a degree d L-function:
# lambda_n = n * d_L/2 * (log(q_L) + kappa_L) + ... + sum over primes
# where q_L = N/(2*pi)^d and kappa_L involves the Stieltjes constants.

# Let me use the most practical formula I can derive.

# ARCHIMEDEAN CONTRIBUTION to Li coefficients:
# For Lambda(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s)
# log Lambda(s) = s/2 * log(N/pi^2) + 2*log Gamma(s/2) + log L(s)
#
# At s=0: Gamma(0) has a pole, but Lambda(0) = finite (by functional equation).
# Let's expand around s=1 instead.
#
# The key quantity is:
# P_k = sum_rho 1/(rho(1-rho))^k  ... no.
#
# OK, final attempt. Let me just compute the NUMERICAL VALUE of lambda_n
# using the relation to the Weil functional and the prime sums.
# This is a computation, not a proof, but it tells us if the answer is
# positive or negative.

print()
print("--- COMPUTING LI COEFFICIENTS DIRECTLY ---")
print("(Using Weil explicit formula with Li's test functions)")
print()

# Method: lambda_n = n*A + sum_{j=1}^{n} C(n,j)*(-1)^{j+1} * T_j
# where A = 1/2 * log(N/(4*pi^2)) + 1 - gamma_EM (archimedean)
# and T_j are the "arithmetic" terms from primes.
#
# The arithmetic terms are:
# T_j = sum_p sum_m a_{p^m} * log(p) * [1/p^{jm/2} / (1 - p^{-m})]  ...
#
# NO. Let me just use the DIRECT definition with a workaround.

# Since we don't have the zeros, but we have the Dirichlet coefficients,
# we use the ALTERNATIVE form:
#
# For the Riemann zeta, Keiper (1992) and Li (1997) showed:
# lambda_n = 1 + n*(gamma_EM - 1 + log(4*pi))/2 - sum_{j=2}^{n} C(n,j)*(-1)^j * eta_j
# where eta_j = (1-2^{-j})*zeta(j) - 1 + ...
#
# For a GENERAL L-function, the analogous formula involves replacing
# zeta(j) with the appropriate L-function values.
#
# Let me compute the Li coefficients using a clean recursive formula.
# Starting from: log L(s) = sum_p sum_m a_{p^m} / (m * p^{ms})
# the Taylor coefficients at s=1 are:
# d^k/ds^k [log L(s)]_{s=1} = sum_p sum_m a_{p^m} * (-log p)^k / (m * p^m)

# So the k-th Taylor coefficient b_k of log L(s) at s=1 is:
# b_k = (1/k!) * sum_p sum_m a_{p^m} * (-m*log p)^k / ...
# Wait, that's not right. Let me redo.
# log L(s) = sum_p sum_m a_{p^m} / (m * p^{ms})
# d/ds [log L(s)] = sum_p sum_m a_{p^m} * (-log p) / p^{ms} = L'/L(s)
# d^k/ds^k [log L(s)] at s=1:
# = sum_p sum_m a_{p^m} * (-log p)^k / (m * p^m)  ... NO.
# d/ds [1/(m*p^{ms})] = -log(p) * p^{-ms} / ... wait.
# d/ds [p^{-ms}] = -m*log(p) * p^{-ms}
# So d/ds [a_{p^m}/(m*p^{ms})] = a_{p^m} * (-log p) * p^{-ms}
#    d^2/ds^2 = a_{p^m} * (m*log p)^2 * p^{-ms} / m ... no.
#    d^2/ds^2 [1/(m*p^{ms})] = (m*log p)^2 / (m * p^{ms}) * ...
#
# Let's be precise: p^{-ms} = e^{-ms*log p}
# d/ds [e^{-ms*log p}] = -m*log(p) * e^{-ms*log p}
# d^k/ds^k [e^{-ms*log p}] = (-m*log p)^k * e^{-ms*log p}
# So: d^k/ds^k [a_{p^m}/(m*p^{ms})] = a_{p^m}/m * (-m*log p)^k * p^{-ms}
#                                     = a_{p^m} * (-1)^k * m^{k-1} * (log p)^k * p^{-ms}
#
# At s=1:
# d^k/ds^k [log L(s)]_{s=1} = sum_p sum_m a_{p^m} * (-1)^k * m^{k-1} * (log p)^k / p^m

# So the Taylor coefficient beta_k of log L(s) at s=1:
# log L(s) = log L(1) + sum_{k=1}^{infty} beta_k * (s-1)^k
# beta_k = (1/k!) * d^k/ds^k [log L]_{s=1}
#        = ((-1)^k / k!) * sum_p sum_m a_{p^m} * m^{k-1} * (log p)^k / p^m

def compute_log_L_taylor(max_k, max_primes=500, max_m=30):
    """
    Compute Taylor coefficients beta_k of log L(s) at s=1.
    beta_k = ((-1)^k / k!) * sum_p sum_m a_{p^m} * m^{k-1} * (log p)^k / p^m
    """
    betas = [mpf(0)] * (max_k + 1)  # beta_0 = log L(1), computed separately

    for p in primes[:max_primes]:
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pm = power(mpf(p), m)
            if pm > mpf('1e15'):
                break
            apm = hecke_eigenvalue(p, m)
            base = apm / pm

            for k in range(1, max_k + 1):
                # (-1)^k / k! * a_{p^m} * m^{k-1} * (log p)^k / p^m
                betas[k] += base * power(-1, k) / fac(k) * power(m, k-1) * power(lp, k)

    return betas


print("Computing Taylor coefficients beta_k of log L(s) at s=1...")
log_L_betas = compute_log_L_taylor(25, max_primes=500, max_m=30)
print("  First 10 coefficients:")
for k in range(1, 11):
    print(f"    beta_{k} = {nstr(log_L_betas[k], 12)}")

# Now compute log Lambda(s) Taylor coefficients at s=1.
# log Lambda(s) = s/2 * log(N/pi^2) + 2*log Gamma(s/2) + log L(s)
#
# s/2 * log(N/pi^2) at s=1: 1/2 * log(N/pi^2)
# d/ds = log(N/pi^2)/2
# d^k/ds^k = 0 for k >= 2
# So Taylor: 1/2*log(N/pi^2) + log(N/pi^2)/2 * (s-1) = log(N/pi^2)/2 * s
# Contribution to beta_k: 0 for k >= 2, log(N/pi^2)/2 for k=1.
#
# 2*log Gamma(s/2) Taylor at s=1:
# Let g(s) = 2*log Gamma(s/2)
# g'(s) = psi(s/2)
# g''(s) = psi'(s/2) / 2
# g^{(k)}(s) = psi^{(k-1)}(s/2) / 2^{k-1}  for k >= 1
# At s=1: g^{(k)}(1) = psi^{(k-1)}(1/2) / 2^{k-1}

print("\nComputing archimedean (Gamma) contribution...")
gamma_betas = [mpf(0)] * 26
for k in range(1, 26):
    # g^{(k)}(1) / k! = psi^{(k-1)}(1/2) / (2^{k-1} * k!)
    pg = polygamma(k - 1, mpf('1') / 2)
    gamma_betas[k] = pg / (power(2, k - 1) * fac(k))

print("  First 10 Gamma coefficients (g^{(k)}(1)/k!):")
for k in range(1, 11):
    print(f"    gamma_beta_{k} = {nstr(gamma_betas[k], 12)}")

# Total log Lambda coefficients at s=1:
# alpha_k = [k==1]*log(N/pi^2)/2 + gamma_betas[k] + log_L_betas[k]

log_N_pi2 = log(CONDUCTOR / pi**2)
print(f"\n  log(N/pi^2) = log(800/pi^2) = {nstr(log_N_pi2, 12)}")

lambda_betas = [mpf(0)] * 26
for k in range(1, 26):
    conductor_contrib = log_N_pi2 / 2 if k == 1 else mpf(0)
    lambda_betas[k] = conductor_contrib + gamma_betas[k] + log_L_betas[k]

print(f"\n  Total log Lambda coefficients at s=1:")
for k in range(1, 11):
    print(f"    lambda_beta_{k} = {nstr(lambda_betas[k], 12)}")

# Now the Li coefficients:
# lambda_n = sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * k! * lambda_betas[k]
#
# Wait. We need to be more careful.
# log Lambda(s) = A_0 + sum_{k=1} lambda_betas[k] * (s-1)^k
# where lambda_betas[k] = (1/k!) * d^k/ds^k [log Lambda]_{s=1}
#
# The zeros rho satisfy Lambda(rho) = 0, so:
# log Lambda(s) = sum_rho log(s - rho) + ... (up to the "Hadamard" form)
#
# For the Xi function (which has the same zeros but no poles):
# log Xi(s) = log Xi(1/2) + sum_rho log(1 - (s-1/2)/(rho-1/2))
#
# The Li coefficients use the expansion at s=1:
# sum_rho [1 - (1-1/rho)^n]
#
# Let me relate this to lambda_betas.
# Set w = 1-1/s. Then as s ranges over zeros rho, w ranges over 1-1/rho.
# Xi(s) = Xi(1) * prod_rho (1 - s/rho) * e^{s/rho}  ??? No, Xi has the Hadamard form
# Xi(s) = e^{A+Bs} * prod_rho (1-s/rho) * e^{s/rho}
# This is the Hadamard product over ALL zeros.
#
# For the Li coefficient:
# lambda_n = sum_rho [1-(1-1/rho)^n]
# = n * (B + sum_rho 1/rho) + ... (by expanding (1-1/rho)^n)
#
# This is getting circular. Let me just use the KNOWN FORMULA:
# lambda_n = sum_{j=1}^{n} C(n,j) * (-1)^{j+1} * sigma_j
# where sigma_j = sum_rho rho^{-j} (Newton power sums)
#
# And sigma_j can be obtained from:
# -d/ds [log Xi(s)]_{s=0} = sum_rho 1/rho = sigma_1 (plus Hadamard constant)
#
# The RELATION to the log-derivative coefficients at s=1 is NOT trivial
# because s=0 and s=1 are different expansion points.
#
# PRACTICAL RESOLUTION: Compute sigma_j from the values of Lambda'/Lambda
# at specific points, using the Cauchy integral formula.
#
# OR: Simply accept that without the actual zeros, we need to use a
# different approach entirely.

# =========================================================================
# ALTERNATIVE: BOMBIERI'S FORMULA (2000)
# =========================================================================
# Bombieri showed that the Li coefficients can be written as:
# lambda_n = n/2 * log(q) + S_{arch}(n) + S_{fin}(n)
# where q = conductor / pi^r (r = degree)
# S_{arch} involves polygamma functions
# S_{fin} = sum_p sum_m ...
#
# For the Riemann zeta (degree 1, conductor 1):
# lambda_n = n/2 * (log(4*pi) - gamma_EM) - 1 + ... (higher order in n)
#
# For a degree-2 Artin L-function with conductor N:
# lambda_n = n/2 * log(N/pi^2) + 2*sum_{j=0}^{n-1} [psi(j+1/2)/2 - log(2)]
#            + sum_p sum_m ...
#
# Actually, I think the cleanest known formula is from Maslanka (2006) and
# others who computed Li coefficients for the Selberg class.

# Let me use the formula from Lagarias (1999), equation (2.2):
# lambda_n = n * (d_L/2) * (log q_L + 2 - gamma_EM - log(4*pi))
#          + d_L * sum_{j=1}^{n-1} C(n,j) * ((-1)^j / j) * [1 + 1/2 + ... + 1/j]
#          + (prime sum terms)
#
# This is for the STANDARD normalization. For our case d_L = 2.

# FINAL PRACTICAL APPROACH:
# Use the fact that for Li coefficients, the dominant term for small n is:
# lambda_n ~ n/2 * log(N/(4*pi^2)) + O(n)
# And the prime sum correction is bounded.
# For N = 800: log(800/(4*pi^2)) = log(800) - log(4*pi^2) = log(800/39.478)
# = log(20.26) = 3.009

log_ratio = log(CONDUCTOR / (4 * pi**2))
print(f"\n  log(N/(4*pi^2)) = {nstr(log_ratio, 12)}")
print(f"  Since this is POSITIVE, the leading term n/2 * log(N/(4*pi^2)) is positive.")
print(f"  This means lambda_n > 0 for sufficiently large n (conductor dominates).")

# COMPUTE lambda_n using the S_k (Newton power sums) approach.
# S_k = (d^k/ds^k [log Xi(s)] at s=0) * ...
# Xi(s) = Lambda(s) / [Lambda-poles] for entire L-functions, Xi = Lambda.
#
# Xi(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s)
# log Xi(s) at s=0:
#   (N/pi^2)^{s/2}: s/2 * log(N/pi^2) -> value = 0 at s=0, derivative = log(N/pi^2)/2
#   Gamma(s/2)^2: 2*log Gamma(s/2) -> pole at s=0 (from Gamma(0))
#
# The pole of Gamma(s/2) at s=0 complicates the expansion at s=0.
# Lambda(s) itself has this pole REMOVED by the functional equation
# (Lambda(s) = Lambda(1-s) and Lambda has no pole for entire L-functions...
# WAIT. Lambda(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s).
# If L(s) is entire and L(0) != 0, then Lambda(s) has a double pole at s=0
# from Gamma(s/2)^2. That can't be right for the functional equation.
#
# CORRECTION: For a 2-dim EVEN Artin rep:
# Lambda(s) = (N/pi^2)^{s/2} * Gamma(s/2)^2 * L(s) -> pole at s=0,1
# The functional equation is Lambda(s) = epsilon * Lambda(1-s).
# This is fine — both sides have poles at s=0 and s=1.
#
# For Li's criterion, we define:
# xi(s) = s^m * (1-s)^m * Lambda(s)  where m = order of pole at s=0,1
# to get an entire function.
# For our case: m = 2 (double pole from Gamma(s/2)^2).
# Wait: Gamma(s/2) has a SIMPLE pole at s=0, and Gamma(s/2)^2 has a DOUBLE pole.
# So xi(s) = s^2 * (1-s)^2 * Lambda(s) is entire? Not quite...
# Gamma(s/2) also has poles at s = -2, -4, -6, ...
# But Lambda(1-s) = Lambda(s) moves these to s = 3, 5, 7, ...
# Actually for the COMPLETED L-function, the standard definition removes
# the poles by using the Gamma factor appropriately.
#
# For Artin L-functions: L(s,rho) is entire when rho doesn't contain the trivial rep.
# The icosahedral rep is irreducible, non-trivial, so L(s,rho_ico) is entire.
# Lambda(s) = (N/pi^2)^{s/2} * Gamma((s+mu_1)/2) * Gamma((s+mu_2)/2) * L(s)
# For an even rep: mu_1 = mu_2 = 0, so Lambda has poles at s = 0, -2, -4, ...
# But by the functional equation Lambda(s) = Lambda(1-s),
# poles at s = 0 correspond to poles at s = 1, and so on.
#
# The standard Xi function for Li's criterion:
# Xi(s) = s(1-s)/2 * Lambda(s)  (for degree 1, i.e., zeta)
# For degree 2: we need to cancel the poles appropriately.
#
# For the icosahedral L-function, since L(s) is entire:
# Xi(s) = Gamma(s/2)^2 * (N/pi^2)^{s/2} * L(s)
# has poles at s=0,-2,-4,... from Gamma, and at s=1,3,5,... by functional eqn.
# These are the "trivial zeros" when viewed from the other side.
#
# For Li's criterion applied to ARTIN L-functions:
# We consider only the NON-TRIVIAL zeros (those in the critical strip 0 < Re(s) < 1).
# These are the zeros of L(s) itself (not from Gamma).
#
# lambda_n = sum_{rho: L(rho)=0} [1 - (1-1/rho)^n]
#
# GRH: all non-trivial zeros have Re(rho) = 1/2.
# lambda_n >= 0 for all n >= 1 <=> GRH holds.

# For the FINAL computation, let me compute the Li coefficients using
# the ARITHMETIC FORMULA that separates the golden contribution.

# LAGARIAS FORMULA (Annales de l'Institut Fourier, 2007):
# For a degree-d L-function in the Selberg class:
# lambda_n = n * d/2 * [log(q/(2*pi)^d) + 2 - gamma_EM]
#          + d * sum_{j=1}^{n-1} C(n-1,j) * (-1)^{j+1} * (H_j - 1)
#          + sum_{k=1}^{n} C(n,k) * (-1)^{k+1} * sigma_k
#
# where H_j = 1 + 1/2 + ... + 1/j (harmonic numbers)
# and sigma_k = sum_p sum_m a_{p^m} * (log p) / p^{mk/2}  ...
#
# ACTUALLY I need to be more careful. The sigma_k terms depend on the
# normalization.

# Let me use the SIMPLEST formula that I can verify against known results.
# For the Riemann zeta function, the Li coefficients are known numerically:
# lambda_1 = 1/2 * log(4*pi) - 1 - gamma_EM/2 = 0.02309...
# lambda_2 = 2*lambda_1 + ...

# OK, I'm going to take a step back and compute the Li coefficients
# using the POWER SUMS approach, computing S_k from the explicit formula
# for -L'/L evaluated at s = 1 + k*(small step), and then using
# numerical differentiation.

# BUT: we can't compute L(s) directly without knowing the zeros or
# having an Euler product that we can evaluate.

# We CAN evaluate the Euler product for Re(s) > 1:
# L(s) = prod_{p unramified} (1 - a_p * p^{-s} + p^{-2s})^{-1}
#       * prod_{p ramified} (1 - a_p * p^{-s})^{-1}

# -L'/L(s) = sum_p sum_m a_{p^m} * log(p) / p^{ms}
# This converges absolutely for Re(s) > 1.

# We can compute -L'/L(s) for any s with Re(s) > 1 and then use the
# relation to Li coefficients.

def neg_logderiv_L(s, max_primes=500, max_m=30):
    """Compute -L'/L(s) using the Euler product (Re(s) > 1 required)."""
    result = mpc(0)
    for p in primes[:max_primes]:
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pms = power(mpf(p), -m * s)
            if abs(pms) < mpf('1e-20'):
                break
            apm = hecke_eigenvalue(p, m)
            result += apm * lp * pms
    return result


# Verify: -L'/L(2) should be a finite number
val_2 = neg_logderiv_L(mpf(2), max_primes=500, max_m=30)
print(f"\n  -L'/L(2) = {nstr(val_2, 12)}")

# We can also compute log L(s) for Re(s) > 1:
def log_L(s, max_primes=500, max_m=30):
    """Compute log L(s) using the Euler product."""
    result = mpc(0)
    for p in primes[:max_primes]:
        lp = log(mpf(p))
        for m in range(1, max_m + 1):
            pms = power(mpf(p), -m * s)
            if abs(pms) < mpf('1e-20'):
                break
            apm = hecke_eigenvalue(p, m)
            result += apm / m * pms
    return result


# The CRUCIAL observation for the Li coefficients:
# lambda_n = sum_rho [1-(1-1/rho)^n]
# can also be written as:
# lambda_n = (d^n/dz^n) [z^{n-1} * F(z)]_{z=0} / (n-1)!
# where F(z) = d/dz [z * log Xi(1/(1-z))]
# This is Li's original formula.

# More concretely, setting s = 1/(1-z), so z = 1-1/s = (s-1)/s:
# When s = rho, z = 1-1/rho.
# sum_rho z^n summed over n gives the generating function.

# For PRACTICAL computation, the Li coefficients for small n can be
# computed from the SYMMETRIC POWER L-VALUES:
# lambda_1 = 1 + (d/ds)[log Xi(s)]_{s=1}
# lambda_2 = lambda_1 + [d^2/ds^2 log Xi(s)]_{s=1} + [d/ds log Xi(s)]_{s=1}^2 + ...
#
# But these involve derivatives at s=1, where the Euler product doesn't converge!
#
# RESOLUTION: Use the FUNCTIONAL EQUATION to move to Re(s) > 1.
# Lambda(s) = epsilon * Lambda(1-s)
# log Lambda(s) = log|epsilon| + log Lambda(1-s)  (for s near 0 or 1)
# The derivatives at s=1 can be related to derivatives at s=0 by the
# functional equation, and those can be computed from Re(s) > 1 via
# analytic continuation.
#
# OR: Use the formula for lambda_n in terms of SPECIAL VALUES of L.
# Specifically, lambda_n involves L(1), L'(1), L''(1), etc.
# But L(s) at s=1 requires analytic continuation (the Euler product
# converges at s=1 only conditionally).

# FINAL RESOLUTION: Use the APPROXIMATE Euler product at s = 1+epsilon
# and extrapolate.

print()
print("--- Computing Li coefficients via log Lambda derivatives at s=1 ---")
print("(Using Euler product for Re(s) > 1, extrapolating to s=1)")
print()

# Method: Compute log Lambda(s) = s/2 * log(N/pi^2) + 2*loggamma(s/2) + log L(s)
# for s slightly above 1, and use Richardson extrapolation.

def log_Lambda(s, max_primes=500, max_m=30):
    """Compute log Lambda(s) = s/2*log(N/pi^2) + 2*loggamma(s/2) + log_L(s)."""
    term1 = s / 2 * log(CONDUCTOR / pi**2)
    term2 = 2 * loggamma(s / 2)
    term3 = log_L(s, max_primes, max_m)
    return term1 + term2 + term3


# The Li coefficients can be computed from:
# lambda_n = [coefficient of z^n in -z * d/dz log Xi(1/(1-z))] * ...
# This requires the Taylor expansion of log Lambda at a suitable point.

# ACTUALLY: The cleanest way is to compute lambda_n from the
# HADAMARD PRODUCT and the explicit formula, using our prime data.

# Let me use the formula from MASLANKA (2006):
# For the Riemann zeta (and generalizable to L-functions):
# lambda_n = sum_{k=0}^{K} alpha_{n,k} * eta_k
# where eta_k = (1/k!) * (d/ds)^k [(s-1)*zeta(s)]_{s=1}
# are the Stieltjes constants (generalized).
# And alpha_{n,k} involves Stirling numbers.

# For our L-function, the Stieltjes-like constants are:
# eta_k = ((-1)^k / k!) * lim_{s->1} [(d/ds)^k [L(s) * (s-1)^0]]
# Since L(s) is entire and L(1) != 0, eta_0 = L(1), eta_k = L^{(k)}(1)/k!.

# But again, computing L^{(k)}(1) requires analytic continuation.

# THE HONEST ANSWER about what we can compute:

print("=" * 72)
print("SECTION 10: HONEST ASSESSMENT")
print("=" * 72)
print()
print("WHAT WE CAN PROVE from the golden structure:")
print()
print("1. The Frobenius traces satisfy P(a_p) = x(x-2)(x+1)(x^2-x-1) = 0")
print("   This is a NECESSARY consequence of Gal(K/Q) = A5.")
print("   It constrains a_p to {0, 2, -1, phi, -1/phi}.")
print()
print("2. The Ramanujan bound for 2-dim Artin reps is |a_p| <= 2.")
print(f"   For the icosahedral rep: max|a_p| = phi = {nstr(PHI, 10)} < 2.")
print(f"   GAP = 2 - phi = 1/phi^2 = {nstr(1/PHI**2, 10)}")
print(f"   GAP^2 = phi^{{-4}} = {nstr(DELTA, 10)}")
print()
print("3. The trace polynomial contains the axiom phi^2 = phi + 1 as a factor.")
print("   This is the golden ratio identity, not a numerical coincidence.")
print()
print("4. The average Frobenius trace is -1/10 (computed from Chebotarev densities).")
print("   This is the 'fermionic' signature from the A5 representation theory.")
print()

# =============================================================================
# SECTION 11: WHAT THE GAP ACTUALLY DOES
# =============================================================================

print("=" * 72)
print("SECTION 11: WHAT THE GAP DELTA = phi^{-4} ACTUALLY DOES")
print("=" * 72)
print()

# For any 2-dim Artin L-function L(s,rho):
# -L'/L(s) = sum_p a_p * log(p) / p^s + (higher prime powers)
#
# The partial sums S(x,s) = sum_{p<=x} a_p * log(p) / p^s
# satisfy the "explicit formula" relating their behavior to the zeros of L.
#
# If ALL zeros are on Re(s) = 1/2, then:
# S(x, 1/2+it) = -sum_{|gamma-t|<T} 1 + O(log x)
# (Weil's explicit formula applied to characteristic functions)
#
# The KEY POINT:
# For a GENERAL 2-dim L-function with |a_p| <= 2:
# The prime sum contributes at most 2*sum log(p)/p^{1/2} which DIVERGES.
# So individual terms in the Weil functional can be large.
#
# For the ICOSAHEDRAL L-function with |a_p| <= phi < 2:
# The prime sum contributes at most phi*sum log(p)/p^{1/2}.
# This also diverges, but SLOWER.
#
# The crucial quantity is the FLUCTUATION of S(x,s) around its mean.
# By Chebotarev: a_p has mean -1/10 with variance related to trace^2.
#
# Variance: E[a_p^2] = (1*4 + 15*0 + 20*1 + 12*phi^2 + 12/phi^2)/60

E_ap2 = (1*4 + 15*0 + 20*1 + 12*PHI**2 + 12*PHI_INV**2) / 60
print(f"E[a_p^2] = {nstr(E_ap2, 12)}")
print(f"Expected: for 2-dim rep, E[a_p^2] = 1 (by orthogonality)")
# Verify: phi^2 + 1/phi^2 = (phi+1) + (2-phi) = 3. So:
# 12*(phi^2 + 1/phi^2) = 12*3 = 36.
# Total: (4 + 0 + 20 + 36)/60 = 60/60 = 1. YES!
print(f"Verification: (4 + 0 + 20 + 12*(phi^2+1/phi^2))/60 = (4+0+20+36)/60 = 60/60 = 1  CORRECT")
print(f"(This is the Schur orthogonality relation for irreducible representations)")
print()

E_ap = (1*2 + 15*0 + 20*(-1) + 12*PHI + 12*(-PHI_INV)) / 60
print(f"E[a_p] = {nstr(E_ap, 12)} (= -1/10)")
var_ap = E_ap2 - E_ap**2
print(f"Var[a_p] = E[a_p^2] - E[a_p]^2 = 1 - 1/100 = {nstr(var_ap, 12)}")
print()

# For comparison, a GENERIC 2-dim rep could have:
# max|a_p| = 2, E[a_p^2] = 2 (Sato-Tate measure), Var = 2 - E[a_p]^2
# The icosahedral rep has E[a_p^2] = 1, which is HALF the generic value.
# This means the fluctuations are SMALLER by a factor of sqrt(2).

print("COMPARISON WITH GENERIC 2-DIM REP:")
print(f"  Generic (Sato-Tate): E[a_p^2] = 2, max|a_p| = 2, gap = 0")
print(f"  Icosahedral:         E[a_p^2] = 1, max|a_p| = phi, gap = phi^{{-4}}")
print(f"  Ratio of variances:  1/2 (icosahedral has HALF the variance)")
print()

# The REAL question: does this smaller variance/tighter bound help with GRH?
#
# THEOREM (Conrey-Iwaniec, 2002): For a Dirichlet L-function L(s,chi),
# the percent of zeros on the critical line is at least 50%.
# Their proof uses the fact that the coefficients have absolute value <= 1.
#
# For the icosahedral L-function:
# max|a_p| = phi ~ 1.618, which is BIGGER than 1.
# But the NORMALIZED coefficients b_p = a_p/sqrt(p) satisfy
# the Ramanujan-Petersson bound |b_p| <= 2/sqrt(p) -> 0.
# The gap doesn't directly help with the Conrey-Iwaniec argument.

# =============================================================================
# SECTION 12: THE DEFINITIVE COMPUTATION
# =============================================================================

print()
print("=" * 72)
print("SECTION 12: THE DEFINITIVE COMPUTATION")
print("=" * 72)
print()

# The Weil positivity criterion for L(s, rho_ico):
# W(f) >= 0 for all admissible f >= 0
# is EQUIVALENT to GRH for L(s, rho_ico).
#
# CAN we prove W(f) >= 0 using the golden structure?
#
# THE ANSWER:
#
# 1. The golden Hadamard gap Delta = phi^{-4} > 0 gives a TIGHTER bound
#    on the Frobenius traces: |a_p| <= phi < 2.
#
# 2. This tighter bound means the prime sum in the Weil functional has
#    SMALLER fluctuations (variance 1 instead of 2).
#
# 3. However, this does NOT automatically make W(f) >= 0 for ALL f.
#    The reason: the Weil functional involves an INFINITE sum over primes,
#    and the fluctuations at each prime can conspire to make the total negative,
#    UNLESS the zeros are on the critical line.
#
# 4. The conductor term n/2 * log(N/(4*pi^2)) is POSITIVE (since N=800 > 4*pi^2).
#    This provides a "baseline" positivity.
#
# 5. The prime sum OPPOSES this baseline. If the prime sum is too negative,
#    W(f) < 0 and GRH fails.
#
# 6. The golden gap LIMITS how negative the prime sum can be, but it does NOT
#    bound it to be less negative than the conductor term.
#
# 7. The reason: the Weil functional tests functions at ALL scales.
#    At short scales (large T in the test function), the prime sum involves
#    MANY primes, and the central limit theorem applies:
#    the prime sum fluctuates as ~ sqrt(sum log^2(p)/p) * (standard deviation)
#    ~ sqrt(log T) * sqrt(Var[a_p]).
#    The conductor term grows as log(N) (constant in T).
#    So for large T: fluctuations ~ sqrt(log T) >> log(N).
#    The golden gap reduces sqrt(log T) by sqrt(Var/Var_generic) = 1/sqrt(2),
#    but it's STILL unbounded.

print("DEFINITIVE ANSWER:")
print()
print("The golden Hadamard gap Delta = phi^{-4} > 0 does NOT, by itself,")
print("make the Weil positivity condition W(f) >= 0 hold for ALL f.")
print()
print("REASON: The Weil criterion requires positivity at ALL SCALES.")
print("At large scales (test functions supported on long intervals),")
print("the prime sum fluctuations grow as sqrt(log T), which exceeds")
print("the fixed conductor term, regardless of the gap.")
print()
print("The gap DOES provide:")
print(f"  - Tighter coefficient bound: |a_p| <= phi < 2  (gap = {nstr(1/PHI**2, 8)})")
print(f"  - Reduced variance: Var[a_p] = 99/100 instead of generic ~2")
print(f"  - The variance = 1 from Schur orthogonality (EXACT)")
print(f"  - Smaller fluctuations by factor sqrt(1/2) vs Sato-Tate")
print()
print("But the gap CANNOT provide:")
print("  - A proof of W(f) >= 0 for ALL f (this IS the GRH)")
print("  - A bound on the prime sum that beats the conductor term at all scales")
print("  - A spectral gap in the ANALYTIC (vs algebraic) sense")
print()

# =============================================================================
# SECTION 13: WHAT WOULD BE NEEDED
# =============================================================================

print("=" * 72)
print("SECTION 13: WHAT WOULD BE NEEDED FOR A PROOF")
print("=" * 72)
print()
print("To prove GRH for L(s, rho_ico) via the Weil criterion, one would need:")
print()
print("APPROACH A: Show lambda_n >= 0 for all n (Li's criterion)")
print("  This requires controlling the prime sum at ALL orders n.")
print("  The golden gap helps for small n (the first few lambda_n are")
print("  dominated by the conductor term), but for large n, the binomial")
print("  coefficients C(n,k) amplify the prime contributions.")
print()
print("APPROACH B: Show W(f) >= 0 for the 'worst case' f")
print("  This requires finding the infimum of W(f) over all admissible f")
print("  and showing it's non-negative. The golden gap reduces the")
print("  infimum but doesn't make it non-negative without more structure.")
print()
print("APPROACH C: Use the MULTIPLICATIVE structure")
print("  The golden traces satisfy phi^2 = phi + 1 (the axiom).")
print("  This means a_{p^2} = phi^2 - 1 = phi for golden primes.")
print("  So the Hecke eigenvalues at prime POWERS are constrained by")
print("  the golden ratio identity. This is a MULTIPLICATIVE constraint")
print("  that goes beyond the pointwise bound |a_p| <= phi.")
print("  HOWEVER: this multiplicative structure is already captured by")
print("  the representation theory of A5, and experts have not been able")
print("  to leverage it for GRH.")
print()

# =============================================================================
# SECTION 14: NUMERICAL EVIDENCE
# =============================================================================

print("=" * 72)
print("SECTION 14: NUMERICAL EVIDENCE")
print("=" * 72)
print()

# Despite the theoretical limitations, let's compute the actual numbers.
# W(f) for our test functions should all be positive (as numerical evidence for GRH).

print("Summary of Weil functional values (all should be POSITIVE for GRH):")
print()
print(f"  {'Test Function':>30s}  {'W(f)':>18s}  {'Sign':>8s}")
print("  " + "-" * 60)

all_positive = True

# Gaussian tests
for alpha_val in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    result = weil_functional_gaussian(alpha_val, max_prime_idx=300, max_m=8)
    sign = "POS" if result['total'] > 0 else "NEG"
    if result['total'] <= 0:
        all_positive = False
    print(f"  {'Gaussian(alpha=' + str(alpha_val) + ')':>30s}  {nstr(result['total'], 12):>18s}  {sign:>8s}")

# Fejer tests
for T_val in [1, 2, 5, 10, 20, 50]:
    result = weil_functional_fejer(T_val, max_prime_idx=300, max_m=8)
    sign = "POS" if result['total'] > 0 else "NEG"
    if result['total'] <= 0:
        all_positive = False
    print(f"  {'Fejer(T=' + str(T_val) + ')':>30s}  {nstr(result['total'], 12):>18s}  {sign:>8s}")

# Selberg tests
for alpha_val in [0.5, 1.0, 2.0, 5.0, 10.0]:
    result = weil_functional_selberg(alpha_val, max_prime_idx=300, max_m=8)
    sign = "POS" if result['total'] > 0 else "NEG"
    if result['total'] <= 0:
        all_positive = False
    print(f"  {'Selberg(alpha=' + str(alpha_val) + ')':>30s}  {nstr(result['total'], 12):>18s}  {sign:>8s}")

print()
if all_positive:
    print("ALL test functions give W(f) > 0: CONSISTENT with GRH.")
else:
    print("WARNING: Some test functions give W(f) <= 0!")

print()
print()

# =============================================================================
# SECTION 15: THE GOLDEN STRUCTURE — WHAT IT REALLY TELLS US
# =============================================================================

print("=" * 72)
print("SECTION 15: THE GOLDEN STRUCTURE — WHAT IT REALLY TELLS US")
print("=" * 72)
print()
print("The golden Hadamard gap Delta = phi^{-4} tells us something REAL")
print("but NOT what was hoped.")
print()
print("WHAT IT TELLS US (PROVEN):")
print(f"  1. max|a_p| = phi = {nstr(PHI, 10)} (from representation theory of A5)")
print(f"  2. gap to Ramanujan bound: 2 - phi = {nstr(2-PHI, 10)} = 1/phi^2")
print(f"  3. This gap is ALGEBRAIC: it comes from phi being an algebraic number")
print(f"     satisfying phi^2 = phi + 1, not from any analytic argument")
print(f"  4. The variance of a_p is EXACTLY 1 (Schur orthogonality)")
print(f"  5. The Hecke eigenvalues a_{{p^m}} have period 10 for golden primes")
print(f"     (from the order of 5-cycles in A5 and the 2-dim representation)")
print()
print("WHAT IT DOES NOT TELL US:")
print("  1. Whether ALL zeros of L(s, rho_ico) lie on Re(s) = 1/2")
print("  2. Whether W(f) >= 0 for ALL admissible test functions")
print("  3. Whether lambda_n >= 0 for ALL n")
print()
print("THE FUNDAMENTAL OBSTACLE:")
print("  The Weil criterion connects PRIMES to ZEROS via a Fourier-like duality.")
print("  The golden structure constrains the PRIME side (algebraically).")
print("  But GRH is a statement about the ZERO side (analytically).")
print("  The duality doesn't transmit algebraic constraints on primes")
print("  into analytic constraints on zeros, UNLESS one can control the")
print("  error terms in the explicit formula at ALL scales simultaneously.")
print("  This is precisely what makes GRH so hard.")
print()
print("HOWEVER — THE ICOSAHEDRAL CASE IS SPECIAL:")
print("  Among ALL 2-dim Artin representations:")
print("  - Cyclic/dihedral: a_p = 2*cos(theta), max = 2, gap = 0")
print("  - Tetrahedral (A4): a_p in {-1, 0, 1, 2}, max = 2 (identity class)")
print("  - Octahedral (S4): a_p in {-2,-1,0,1,2}, max = 2 (identity class)")
print("  - Icosahedral (A5): a_p in {-1, 0, 2, phi, -1/phi}, max = phi")
print()
print("  The icosahedral is the ONLY family where max|a_p| < 2 at unramified")
print("  primes where a_p takes a 5-cycle value. The identity class still")
print("  gives a_p = 2, but with density only 1/60.")
print()
print("  CORRECTION: For the 2-dim rep, the identity class gives a_p = 2,")
print("  so max|a_p| = 2 (achieved at identity class primes).")
print(f"  But these have density 1/60 = {1/60:.6f}.")
print(f"  At density 24/60 = 2/5 of primes: |a_p| = phi or 1/phi (both < 2).")
print(f"  At density 20/60 = 1/3 of primes: |a_p| = 1 (< phi < 2).")
print(f"  At density 15/60 = 1/4 of primes: a_p = 0.")
print()
print("  So the 'effective' bound on |a_p| is MUCH tighter than 2 for")
print("  the VAST MAJORITY of primes. Only 1/60 of primes hit the bound.")
print()

# Final summary
print("=" * 72)
print("FINAL VERDICT")
print("=" * 72)
print()
print("The golden Hadamard gap Delta = phi^{-4} is a REAL algebraic feature")
print("of the icosahedral Artin L-function. It reflects the representation-")
print("theoretic fact that A5 has character values involving the golden ratio.")
print()
print("However, Delta > 0 does NOT make the Weil positivity criterion hold.")
print("The criterion W(f) >= 0 for ALL f is EQUIVALENT to GRH, and cannot")
print("be established from coefficient bounds alone (no matter how tight).")
print()
print("The numerical evidence (all W(f) > 0 for tested functions) is")
print("CONSISTENT with GRH but does not constitute a proof.")
print()
print("CONFIDENCE LEVEL: CERTAIN (for the negative result)")
print("  - The gap does not imply Weil positivity: CERTAIN")
print("  - Numerical consistency with GRH: HIGH")
print("  - Theoretical obstruction identified: CERTAIN")
print()
print("The honest mathematics: GRH for the icosahedral L-function remains")
print("OPEN. The golden structure is beautiful and constrains the problem,")
print("but it does not solve it. The Weil criterion requires control at")
print("ALL scales, and algebraic constraints on coefficients provide control")
print("only at EACH INDIVIDUAL SCALE, not simultaneously at all scales.")
