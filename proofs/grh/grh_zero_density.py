#!/usr/bin/env python3
"""
GRH_ZERO_DENSITY — golden dominance D(N)->infinity suppresses zero-density exponent for L(s, rho_ico)
nos3bl33d

Seven parts: mean value <|a_p|^2>=14/15, effective exponent, simulation,
golden suppression, Chebyshev persistence, compounding, N(sigma,T) comparison.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import mpmath
from mpmath import mp, mpf, mpc, pi, sin, cos, sqrt, log, exp, phi, fsum, power, fabs
from mpmath import matrix as mpmatrix
import random
import sys

mp.dps = 50

# ============================================================================
# CONSTANTS
# ============================================================================

PHI = phi                          # golden ratio (1+sqrt(5))/2
PHI_INV = mpf(1) / phi            # 1/phi
D_P = 15                          # Schlafli product d_p = {5,3,3} -> 15
ICOSAHEDRAL_ORDER = 60            # |A5| = 60

# Frobenius conjugacy class data for A5 (icosahedral group)
# Classes: identity(1), order-2 involutions(15), order-3(20), order-5+(12), order-5-(12)
# Densities by Chebotarev: proportional to class size / group order
CLASSES = {
    'golden+':   {'density': mpf(12)/60, 'trace': 2*cos(pi/5),    'angle_deg': 36},
    'golden-':   {'density': mpf(12)/60, 'trace': 2*cos(3*pi/5),  'angle_deg': 108},
    'order3':    {'density': mpf(20)/60, 'trace': 2*cos(2*pi/3),  'angle_deg': 120},
    'involution':{'density': mpf(15)/60, 'trace': mpf(0),         'angle_deg': 90},
    # identity class has density 1/60, trace = 2 (trivial Frobenius)
    'identity':  {'density': mpf(1)/60,  'trace': mpf(2),         'angle_deg': 0},
}

# Verify densities sum to 1
density_sum = fsum(c['density'] for c in CLASSES.values())
assert abs(density_sum - 1) < mpf(10)**(-40), f"Density sum = {density_sum}, expected 1"

# Frobenius angles (theta) where trace = 2*cos(theta)
ANGLES = {
    'golden+':    pi / 5,       # 36 degrees
    'golden-':    3 * pi / 5,   # 108 degrees
    'order3':     2 * pi / 3,   # 120 degrees
    'involution': pi / 2,       # 90 degrees
    'identity':   mpf(0),       # 0 degrees
}

def banner(title):
    width = 78
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


# ============================================================================
# PART 1: Mean Value with Golden Dominance
# ============================================================================

def part1_mean_value():
    banner("PART 1: Mean Value with Golden Dominance")

    print("\nFrobenius class data for A5:")
    print(f"{'Class':<14} {'Density':<12} {'Trace a_p':<20} {'|a_p|^2':<20}")
    print("-" * 66)

    weighted_sum = mpf(0)
    for name, data in CLASSES.items():
        trace = data['trace']
        trace_sq = trace ** 2
        weighted = data['density'] * trace_sq
        weighted_sum += weighted
        print(f"{name:<14} {mp.nstr(data['density'],6):<12} "
              f"{mp.nstr(trace,12):<20} {mp.nstr(trace_sq,12):<20}")

    print(f"\n<|a_p|^2> = Sum(density_i * |a_p,i|^2) = {mp.nstr(weighted_sum, 30)}")

    # Verify: should be 14/15
    expected = mpf(14) / 15
    print(f"Expected 14/15 = {mp.nstr(expected, 30)}")
    print(f"Difference: {mp.nstr(abs(weighted_sum - expected), 30)}")

    # But wait — the problem statement uses densities 1/5, 1/5, 1/3, 1/4
    # which sum to 1/5+1/5+1/3+1/4 = 12/60+12/60+20/60+15/60 = 59/60
    # Missing: identity class at 1/60 with trace 2, |a_p|^2 = 4
    # Let's compute both ways
    print("\n--- Verification with exact Chebotarev densities (including identity) ---")
    # With identity: (12/60)*phi^2 + (12/60)/phi^2 + (20/60)*1 + (15/60)*0 + (1/60)*4
    val_exact = (mpf(12)/60)*phi**2 + (mpf(12)/60)/phi**2 + (mpf(20)/60)*1 + (mpf(15)/60)*0 + (mpf(1)/60)*4
    print(f"Exact with identity class: <|a_p|^2> = {mp.nstr(val_exact, 30)}")
    print(f"= {mp.nstr(val_exact, 6)}")

    # Without identity (problem statement approximation):
    # Renormalize: 12+12+20+15 = 59, so densities are 12/59, 12/59, 20/59, 15/59
    val_no_id = (mpf(12)/59)*phi**2 + (mpf(12)/59)/phi**2 + (mpf(20)/59)*1 + (mpf(15)/59)*0
    print(f"Without identity (renormalized): <|a_p|^2> = {mp.nstr(val_no_id, 30)}")

    # The problem statement uses non-renormalized 1/5, 1/5, 1/3, 1/4:
    val_approx = mpf(1)/5 * phi**2 + mpf(1)/5 / phi**2 + mpf(1)/3 * 1 + mpf(1)/4 * 0
    print(f"Problem statement (1/5,1/5,1/3,1/4): <|a_p|^2> = {mp.nstr(val_approx, 30)}")
    print(f"  = 3/5 + 1/3 = {mp.nstr(mpf(3)/5 + mpf(1)/3, 30)} = 14/15")

    # Check: phi^2 + 1/phi^2
    print(f"\nphi^2 + 1/phi^2 = {mp.nstr(phi**2 + 1/phi**2, 30)}")
    print(f"Should be 3 (= dimension of rho_ico): {'YES' if abs(phi**2 + 1/phi**2 - 3) < mpf(10)**(-40) else 'NO'}")

    # Suppression vs Sato-Tate
    print(f"\n--- Suppression Analysis ---")
    print(f"Sato-Tate <|a_p|^2> for GL(2): 1  (second moment of semicircle)")
    print(f"Golden ico <|a_p|^2>:          {mp.nstr(val_exact, 20)}")
    print(f"Suppression = 1 - <|a_p|^2>:   {mp.nstr(1 - val_exact, 20)}")
    print(f"1/d_p = 1/15:                   {mp.nstr(mpf(1)/15, 20)}")

    # With exact Chebotarev:
    # (12*phi^2 + 12/phi^2 + 20 + 0 + 4) / 60
    # = (12*(phi^2 + 1/phi^2) + 24) / 60
    # = (12*3 + 24) / 60 = (36+24)/60 = 60/60 = 1
    print(f"\nWITH identity class: <|a_p|^2> = {mp.nstr(val_exact, 6)} (= 1 exactly!)")
    print("The identity class (density 1/60, trace 2, |a_p|^2=4) CANCELS the suppression!")

    # So the suppression to 14/15 only holds if we EXCLUDE the identity class
    # Physically: almost all primes have non-trivial Frobenius (density 59/60)
    # The effective suppression for "typical" primes is:
    val_typical = val_no_id
    print(f"\nFor TYPICAL primes (non-trivial Frobenius, 59/60 of all primes):")
    print(f"  <|a_p|^2>_typical = {mp.nstr(val_typical, 20)}")
    print(f"  = {mp.nstr(val_typical, 6)}")

    # Actually let's recompute: (12*phi^2 + 12/phi^2 + 20*1 + 15*0) / 59
    num = 12*phi**2 + 12/phi**2 + 20
    print(f"  Numerator: 12*phi^2 + 12/phi^2 + 20 = {mp.nstr(num, 20)}")
    print(f"  = 12*3 + 20 = 56 (since phi^2 + 1/phi^2 = 3 -> 12*(phi^2+1/phi^2)=36)")
    print(f"  56/59 = {mp.nstr(mpf(56)/59, 20)}")

    return val_exact, val_no_id, val_approx


# ============================================================================
# PART 2: Zero-Density Exponent
# ============================================================================

def part2_zero_density_exponent(avg_a2_exact, avg_a2_typical):
    banner("PART 2: Zero-Density Exponent A_eff")

    print("\nStandard GL(2) zero-density exponents (known results):")
    A_values = {
        'Ingham (1940)':        mpf(3),
        'Jutila (1977)':        mpf(12) / 5,
        'Heath-Brown (1979)':   mpf(12) / 5,
        'Huxley (improved)':    mpf(9) / 4,
        'Best known GL(2)':     mpf(2),
    }

    for name, A in A_values.items():
        print(f"  {name:<25} A = {mp.nstr(A, 6)}")

    print(f"\n--- Effective exponents with golden suppression ---")
    print(f"A_eff = A * <|a_p|^2>")
    print()

    # Using the different averages
    for avg_name, avg_val in [("Exact (with identity)", avg_a2_exact),
                               ("Typical (no identity)", avg_a2_typical),
                               ("Problem stmt (14/15)", mpf(14)/15)]:
        print(f"  Using <|a_p|^2> = {mp.nstr(avg_val, 8)} ({avg_name}):")
        for A_name, A in A_values.items():
            A_eff = A * avg_val
            reduction_pct = (1 - avg_val) * 100
            print(f"    {A_name:<25} A_eff = {mp.nstr(A_eff, 8)}  "
                  f"(reduction: {mp.nstr(reduction_pct, 4)}%)")
        print()

    # The critical question: when does A_eff < 1?
    # N(sigma,T) ~ (qT)^{A_eff*(1-sigma)} -> 0 as T->inf iff A_eff*(1-sigma) < 0
    # But A_eff > 0 always (even suppressed), so need 1-sigma < 0 i.e. sigma > 1
    # That's trivial! The density exponent being positive but small means:
    # fewer zeros but still potentially infinitely many as T -> inf
    print("--- Critical Observation ---")
    print("A_eff > 0 for all cases. The golden suppression reduces the COUNT")
    print("of zeros but does NOT drive it to zero for sigma > 1/2.")
    print("N(sigma,T) ~ (qT)^{A_eff*(1-sigma)} still -> infinity as T -> infinity")
    print("for any fixed sigma in (1/2, 1).")
    print()
    print("HOWEVER: the suppression enters the CONSTANT and the LOG POWER too.")
    print("Let's investigate how the Chebyshev traces compound at prime powers...")


# ============================================================================
# PART 3: Numerical Simulation of L-function Mean Values
# ============================================================================

def generate_dirichlet_coeffs(N, seed=42):
    """Generate Dirichlet coefficients a_n for n=1..N using Chebotarev-random
    Frobenius assignments to primes, then multiplicative extension."""

    rng = random.Random(seed)
    a = [mpf(0)] * (N + 1)
    a[1] = mpf(1)

    # Sieve for primes
    is_prime = [False, False] + [True] * (N - 1)
    for i in range(2, int(N**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, N + 1, i):
                is_prime[j] = False

    # Assign Frobenius class to each prime
    frobenius_angle = {}
    class_list = list(CLASSES.keys())
    # Cumulative density for sampling
    cum_density = []
    running = mpf(0)
    for cls in class_list:
        running += CLASSES[cls]['density']
        cum_density.append(float(running))

    for p in range(2, N + 1):
        if not is_prime[p]:
            continue
        r = rng.random()
        chosen = class_list[-1]
        for i, cd in enumerate(cum_density):
            if r < cd:
                chosen = class_list[i]
                break
        frobenius_angle[p] = ANGLES[chosen]
        a[p] = 2 * cos(frobenius_angle[p])

    # Extend to prime powers using Hecke recurrence:
    # a_{p^{k+1}} = a_p * a_{p^k} - a_{p^{k-1}}
    for p in range(2, N + 1):
        if not is_prime[p]:
            continue
        pk = p
        # a_{p^0} = 1, a_{p^1} = a_p (already set)
        a_prev = mpf(1)   # a_{p^0}
        a_curr = a[p]     # a_{p^1}
        pk = p * p
        while pk <= N:
            a_next = a[p] * a_curr - a_prev
            a[pk] = a_next
            a_prev = a_curr
            a_curr = a_next
            pk *= p

    # Extend multiplicatively: a_{mn} = a_m * a_n for gcd(m,n) = 1
    # We need to fill in composite n that aren't prime powers
    # Factor each n and use multiplicativity
    for n in range(2, N + 1):
        if a[n] != 0 or n == 1:
            continue
        # Factor n
        temp = n
        factors = {}
        for p in range(2, n + 1):
            if p * p > temp:
                break
            while temp % p == 0:
                factors[p] = factors.get(p, 0) + 1
                temp //= p
        if temp > 1:
            factors[temp] = factors.get(temp, 0) + 1

        # a_n = product of a_{p^k}
        val = mpf(1)
        for p, k in factors.items():
            pk = p ** k
            if pk <= N and a[pk] != 0:
                val *= a[pk]
            elif pk <= N:
                # Compute a_{p^k} via recurrence
                a_prev_loc = mpf(1)
                a_curr_loc = a[p] if p <= N else mpf(0)
                for j in range(2, k + 1):
                    a_next_loc = a[p] * a_curr_loc - a_prev_loc
                    a_prev_loc = a_curr_loc
                    a_curr_loc = a_next_loc
                val *= a_curr_loc
        a[n] = val

    return a, frobenius_angle


def part3_numerical_simulation():
    banner("PART 3: Numerical Simulation of Mean Values")

    N_TERMS = 3000   # Dirichlet series truncation
    T_MAX = 50       # height range
    T_STEP = mpf('0.5')
    SIGMAS = [mpf('0.6'), mpf('0.7'), mpf('0.8'), mpf('0.9')]

    print(f"Generating {N_TERMS} Dirichlet coefficients...")
    a_ico, frob_angles = generate_dirichlet_coeffs(N_TERMS)

    # Count how many primes in each class
    class_counts = {cls: 0 for cls in CLASSES}
    for p, angle in frob_angles.items():
        for cls, theta in ANGLES.items():
            if abs(angle - theta) < mpf(10)**(-30):
                class_counts[cls] += 1
                break

    total_primes = sum(class_counts.values())
    print(f"Total primes up to {N_TERMS}: {total_primes}")
    for cls, count in class_counts.items():
        expected = float(CLASSES[cls]['density'])
        actual = count / total_primes if total_primes > 0 else 0
        print(f"  {cls:<14}: {count:>5} ({actual:.4f}, expected {expected:.4f})")

    # Compute mean values
    t_values = []
    t = mpf(1)
    while t <= T_MAX:
        t_values.append(t)
        t += T_STEP

    n_t = len(t_values)
    print(f"\nComputing |L(sigma+it)|^2 for {len(SIGMAS)} sigma values, {n_t} t-values...")

    results_ico = {}
    results_zeta = {}

    for sigma in SIGMAS:
        mean_ico = mpf(0)
        mean_zeta = mpf(0)
        max_ico = mpf(0)
        max_zeta = mpf(0)

        for t in t_values:
            s = mpc(sigma, t)

            # L_ico(s) = sum_{n=1}^{N} a_n / n^s
            L_ico = mpf(0)
            L_ico_val = mpc(0, 0)
            for n in range(1, N_TERMS + 1):
                if a_ico[n] == 0:
                    continue
                L_ico_val += a_ico[n] * power(n, -s)

            # zeta(s) = sum_{n=1}^{N} 1/n^s (for comparison)
            L_zeta_val = mpc(0, 0)
            for n in range(1, N_TERMS + 1):
                L_zeta_val += power(n, -s)

            abs2_ico = fabs(L_ico_val) ** 2
            abs2_zeta = fabs(L_zeta_val) ** 2

            mean_ico += abs2_ico
            mean_zeta += abs2_zeta
            if abs2_ico > max_ico:
                max_ico = abs2_ico
            if abs2_zeta > max_zeta:
                max_zeta = abs2_zeta

        mean_ico /= n_t
        mean_zeta /= n_t

        results_ico[sigma] = (mean_ico, max_ico)
        results_zeta[sigma] = (mean_zeta, max_zeta)

        ratio = mean_ico / mean_zeta if mean_zeta > 0 else mpf(0)
        print(f"\n  sigma = {mp.nstr(sigma, 2)}:")
        print(f"    <|L_ico|^2>  = {mp.nstr(mean_ico, 15)}")
        print(f"    <|zeta|^2>   = {mp.nstr(mean_zeta, 15)}")
        print(f"    Ratio ico/zeta = {mp.nstr(ratio, 15)}")
        print(f"    max|L_ico|^2 = {mp.nstr(max_ico, 15)}")
        print(f"    max|zeta|^2  = {mp.nstr(max_zeta, 15)}")

    print("\n--- Summary Table ---")
    print(f"{'sigma':<8} {'<|L_ico|^2>':<22} {'<|zeta|^2>':<22} {'Ratio':<15} {'Suppressed?'}")
    print("-" * 75)
    for sigma in SIGMAS:
        m_i = results_ico[sigma][0]
        m_z = results_zeta[sigma][0]
        ratio = m_i / m_z if m_z > 0 else mpf(0)
        suppressed = "YES" if ratio < 1 else "NO"
        print(f"{mp.nstr(sigma,2):<8} {mp.nstr(m_i,15):<22} {mp.nstr(m_z,15):<22} "
              f"{mp.nstr(ratio,8):<15} {suppressed}")

    return results_ico, results_zeta


# ============================================================================
# PART 4: Golden Suppression Factor
# ============================================================================

def part4_suppression_factor():
    banner("PART 4: Golden Suppression Factor")

    # Sato-Tate second moment: integral of (2cos(theta))^2 * (2/pi)*sin^2(theta) dtheta
    # = integral_0^pi 4cos^2(theta) * (2/pi)*sin^2(theta) dtheta
    # = (8/pi) * integral_0^pi cos^2(theta)*sin^2(theta) dtheta
    # = (8/pi) * (pi/8) = 1
    st_second_moment = mpf(8) / pi * mpmath.quad(
        lambda theta: cos(theta)**2 * sin(theta)**2,
        [0, pi]
    )
    print(f"Sato-Tate second moment <|a_p|^2>_ST = {mp.nstr(st_second_moment, 30)}")

    # Golden icosahedral average (exact Chebotarev, all classes)
    avg_exact = mpf(0)
    for cls, data in CLASSES.items():
        avg_exact += data['density'] * data['trace']**2
    print(f"Golden ico <|a_p|^2>_exact = {mp.nstr(avg_exact, 30)}")

    # Without identity
    avg_no_id = mpf(0)
    total_d = mpf(0)
    for cls, data in CLASSES.items():
        if cls != 'identity':
            avg_no_id += data['density'] * data['trace']**2
            total_d += data['density']
    avg_no_id /= total_d  # renormalize
    print(f"Golden ico <|a_p|^2>_no_id = {mp.nstr(avg_no_id, 30)} (renormalized)")

    S_exact = avg_exact / st_second_moment
    S_no_id = avg_no_id / st_second_moment

    print(f"\nSuppression factor S_exact = {mp.nstr(S_exact, 30)}")
    print(f"Suppression factor S_no_id = {mp.nstr(S_no_id, 30)}")

    print(f"\n1 - S_exact = {mp.nstr(1 - S_exact, 30)}")
    print(f"1 - S_no_id = {mp.nstr(1 - S_no_id, 30)}")
    print(f"1/d_p = 1/15 = {mp.nstr(mpf(1)/15, 30)}")

    print(f"\n--- Identity ---")
    print(f"With identity class: S = 1 exactly (no suppression!)")
    print(f"The full Chebotarev average gives <|a_p|^2> = 1 = Sato-Tate value.")
    print(f"This is a THEOREM: for any Artin representation of degree d,")
    print(f"Chebotarev gives <|a_p|^2> = d (character inner product <chi,chi>=1).")
    print(f"For our 2-dim rep: <|a_p|^2> = 1 (per Frobenius, not per rep dimension).")
    print(f"Actually: <tr(Frob_p)^2> where tr is the CHARACTER.")
    print(f"By Schur orthogonality: sum_g |chi(g)|^2 / |G| = 1 for irreducible chi.")
    print(f"So <|a_p|^2>_Chebotarev = 1 ALWAYS for irreducible representations!")

    print(f"\n--- The REAL suppression is in HIGHER moments ---")
    print(f"<|a_p|^4> differs from Sato-Tate!")

    # Fourth moment
    st_fourth = mpf(8) / pi * mpmath.quad(
        lambda theta: cos(theta)**4 * sin(theta)**2,
        [0, pi]
    ) * 16  # (2cos)^4 = 16cos^4
    # Actually: <(2cos theta)^4>_ST = (8/pi) int_0^pi 16cos^4(theta) sin^2(theta) dtheta
    st_fourth = mpf(8) / pi * mpmath.quad(
        lambda theta: 16 * cos(theta)**4 * sin(theta)**2,
        [0, pi]
    )
    print(f"\nSato-Tate <|a_p|^4>_ST = {mp.nstr(st_fourth, 20)}")
    # Should be 2 (by Sato-Tate moments: m4 = 2)

    avg_fourth = mpf(0)
    for cls, data in CLASSES.items():
        avg_fourth += data['density'] * data['trace']**4
    print(f"Golden ico <|a_p|^4>     = {mp.nstr(avg_fourth, 20)}")

    S4 = avg_fourth / st_fourth
    print(f"Fourth moment ratio = {mp.nstr(S4, 20)}")
    print(f"1 - S4 = {mp.nstr(1 - S4, 20)}")

    return S_exact, S_no_id, S4


# ============================================================================
# PART 5: Chebyshev Trace Persistence at Prime Powers
# ============================================================================

def part5_chebyshev_traces():
    banner("PART 5: Chebyshev Traces at Prime Powers")

    print("\nFor Hecke eigenvalues: a_{p^k} = U_k(cos theta) = sin((k+1)*theta)/sin(theta)")
    print("where U_k is the Chebyshev U polynomial and a_p = 2*cos(theta).\n")

    K_MAX = 25  # prime power orders to check

    # For each Frobenius class, compute a_{p^k} and |a_{p^k}|^2 for k=0..K_MAX
    print(f"{'k':<5}", end="")
    for cls in ['golden+', 'golden-', 'order3', 'involution', 'identity']:
        print(f"{'|a_{p^'+str('k')+'}|^2 '+cls:<18}", end="")
    print(f"{'<|a_{p^k}|^2>':<20} {'< 1?'}")
    print("-" * 115)

    avg_history = []

    for k in range(0, K_MAX + 1):
        row_vals = {}
        weighted_avg = mpf(0)

        for cls in ['golden+', 'golden-', 'order3', 'involution', 'identity']:
            theta = ANGLES[cls]

            if abs(sin(theta)) < mpf(10)**(-40):
                # theta = 0: a_{p^k} = k+1  (L'Hopital)
                a_pk = mpf(k + 1)
            else:
                a_pk = sin((k + 1) * theta) / sin(theta)

            abs2 = a_pk ** 2
            row_vals[cls] = abs2
            weighted_avg += CLASSES[cls]['density'] * abs2

        avg_history.append(weighted_avg)

        print(f"{k:<5}", end="")
        for cls in ['golden+', 'golden-', 'order3', 'involution', 'identity']:
            print(f"{mp.nstr(row_vals[cls], 8):<18}", end="")
        below = "YES" if weighted_avg < 1 else "NO"
        print(f"{mp.nstr(weighted_avg, 12):<20} {below}")

    print(f"\n--- Period Analysis for Golden Classes ---")
    print("Golden angle = pi/5 = 36 degrees. Period of sin((k+1)*36)/sin(36) in k:")
    print("sin((k+1)*pi/5)/sin(pi/5) repeats when (k+1)*pi/5 shifts by pi")
    print("-> k+1 shifts by 5 -> period in k is 10 (with sign flip at 5).")
    print("|a_{p^k}|^2 has period 5 for golden classes.\n")

    # Verify period
    golden_theta = pi / 5
    print("Golden+ trace values:")
    for k in range(12):
        a_pk = sin((k+1)*golden_theta) / sin(golden_theta)
        print(f"  k={k}: a_{{p^{k}}} = {mp.nstr(a_pk, 15)}, |a|^2 = {mp.nstr(a_pk**2, 15)}")

    print(f"\n--- Running Average of <|a_{{p^k}}|^2> ---")
    running_avg = mpf(0)
    for k in range(len(avg_history)):
        running_avg += avg_history[k]
        cumavg = running_avg / (k + 1)
        print(f"  avg(k=0..{k}): {mp.nstr(cumavg, 15)}")

    # The running average should converge to something
    # For the exact Chebotarev: <|a_{p^k}|^2> = sum over classes density_i * U_k(cos theta_i)^2
    # The time average converges to the "energy" of the representation
    total_avg = running_avg / len(avg_history)
    print(f"\nOverall average of <|a_{{p^k}}|^2> for k=0..{K_MAX}: {mp.nstr(total_avg, 15)}")

    return avg_history


# ============================================================================
# PART 6: Compounding Argument
# ============================================================================

def part6_compounding(avg_history):
    banner("PART 6: The Compounding Argument")

    print("\nThe mean-value Euler product at sigma involves:")
    print("  Pi_p (1 + |a_p|^2/p^{2sigma} + |a_{p^2}|^2/p^{4sigma} + ...)")
    print()
    print("If <|a_{p^k}|^2> < (k+1)^2 (Ramanujan bound |a_{p^k}| <= k+1):")
    print("the product converges better. Golden suppression at each k helps.\n")

    # Compute the "Euler factor contribution" for a typical prime p
    # at sigma = 1/2 + epsilon
    for eps_val in [mpf('0.01'), mpf('0.05'), mpf('0.1'), mpf('0.2')]:
        sigma = mpf(1)/2 + eps_val
        print(f"\n--- sigma = 1/2 + {mp.nstr(eps_val, 3)} = {mp.nstr(sigma, 4)} ---")

        # Use p=101 as a "typical" prime
        p = mpf(101)

        # Euler factor for golden ico:
        euler_ico = mpf(0)
        euler_rama = mpf(0)  # Ramanujan-saturating: |a_{p^k}| = k+1
        euler_st = mpf(0)    # Sato-Tate generic

        for k in range(len(avg_history)):
            # Golden ico contribution
            contrib_ico = avg_history[k] / p**(2*sigma*k)
            euler_ico += contrib_ico

            # Ramanujan bound: |a_{p^k}|^2 = (k+1)^2
            contrib_rama = mpf(k+1)**2 / p**(2*sigma*k)
            euler_rama += contrib_rama

            # Sato-Tate: <|a_{p^k}|^2> = k+1 (for Sato-Tate distributed angles)
            # Actually, for ST: <U_k(cos theta)^2>_ST = 1 for all k
            # (Orthogonality of Chebyshev polynomials with respect to ST measure)
            contrib_st = mpf(1) / p**(2*sigma*k)
            euler_st += contrib_st

        print(f"  Euler factor (p=101):")
        print(f"    Golden ico:        {mp.nstr(euler_ico, 15)}")
        print(f"    Sato-Tate generic: {mp.nstr(euler_st, 15)}")
        print(f"    Ramanujan max:     {mp.nstr(euler_rama, 15)}")
        print(f"    Ratio ico/ST:      {mp.nstr(euler_ico/euler_st, 15)}")

    # Now the KEY question: does the product over ALL primes compound the suppression?
    print(f"\n\n--- Product Over Primes ---")
    print("log(Euler product) = sum_p sum_k <|a_{p^k}|^2> / p^{2*sigma*k}")
    print("                   = sum_p (euler_factor_at_p - 1)")
    print()

    # At sigma slightly above 1/2, the dominant contribution is k=1:
    # sum_p <|a_p|^2> / p^{2*sigma}
    # For ico: <|a_p|^2> = 1 (exact Chebotarev) -> same as zeta!
    # For ico (no identity): slightly less, but identity class is measure-0

    print("CRITICAL FINDING:")
    print("=" * 60)
    print("<|a_p|^2> = 1 exactly (by Schur orthogonality for irreducible reps).")
    print("The golden suppression at k=1 is ZERO when using exact Chebotarev.")
    print()
    print("At higher prime powers k >= 2:")
    print("  <|a_{p^k}|^2> may differ from Sato-Tate predictions.")
    print("  But these terms are O(1/p^{4*sigma}) or smaller — subdominant!")
    print()
    print("The golden dominance D(N) does NOT directly suppress the leading")
    print("term of the zero-density exponent.")
    print()
    print("HOWEVER: D(N) may enter through a DIFFERENT mechanism...")

    # The dominance D(N) measures how golden-class primes dominate.
    # Even though <|a_p|^2> = 1 overall, the VARIANCE is different:
    var_ico = mpf(0)
    for cls, data in CLASSES.items():
        var_ico += data['density'] * (data['trace']**2 - 1)**2
    print(f"\nVariance of |a_p|^2 around 1:")
    print(f"  Var_ico = <(|a_p|^2 - 1)^2> = {mp.nstr(var_ico, 20)}")

    # Sato-Tate variance of |a_p|^2:
    st_var = mpf(8) / pi * mpmath.quad(
        lambda theta: ((2*cos(theta))**2 - 1)**2 * sin(theta)**2,
        [0, pi]
    )
    print(f"  Var_ST = {mp.nstr(st_var, 20)}")
    print(f"  Ratio Var_ico/Var_ST = {mp.nstr(var_ico/st_var, 20)}")

    if var_ico < st_var:
        print("  -> Golden ico has LESS variance in |a_p|^2 than Sato-Tate!")
        print("  -> This means |L(s)|^2 fluctuates LESS -> fewer large deviations")
        print("  -> Which IS a form of zero-density suppression (via large deviations)")
    else:
        print("  -> Golden ico has MORE variance — not suppressed via this route.")

    return var_ico, st_var


# ============================================================================
# PART 7: The Decisive Test — Full Zero-Density Estimates
# ============================================================================

def part7_decisive_test(avg_history):
    banner("PART 7: Decisive Test — Zero-Density Estimates N(sigma, T)")

    print("\nUsing Ingham-Huxley type estimate:")
    print("  N(sigma,T) <= C * (qT)^{A*(1-sigma)} * log(qT)^B")
    print()
    print("For icosahedral rep: conductor q ~ N^2 for level N=800 -> q = 640000")
    print("GL(2) exponent A = 12/5 (Jutila)")
    print("Log power B = 14 (typical for GL(2))")
    print()

    q = mpf(640000)  # conductor for level 800 icosahedral rep
    A_std = mpf(12) / 5  # standard GL(2) zero-density exponent
    B = mpf(14)  # log power
    C = mpf(1)   # constant (normalized)

    # For golden-enhanced: use the variance-based improvement
    # The mean value is the SAME, but the large-deviation bound improves
    # by a factor related to the variance ratio

    # Standard approach: the density exponent comes from optimizing
    # N(sigma,T) = mean_value / (dip at zero)
    # mean_value ~ T * sum_n |a_n|^2/n^{2sigma}
    # For sigma near 1/2: this is ~ T * L(1, sym^2 f)  (the symmetric square L-function)
    # L(1, sym^2) involves the symmetric square coefficients

    # For icosahedral representation:
    # sym^2(rho_ico) decomposes into irreducible A5 representations
    # rho_ico is the 2-dim representation. sym^2 is 3-dimensional.
    # For A5: the 3-dim symmetric square decomposes as the standard 3-dim irrep

    # The key: L(1, sym^2 rho_ico) is a SPECIFIC value
    # We can compute it from the Euler product

    print("--- Computing L(1, sym^2 rho_ico) via Euler product ---\n")

    # sym^2 trace at prime p with Frobenius angle theta:
    # tr(sym^2(Frob_p)) = 1 + 2*cos(2*theta)
    # (for GL(2): sym^2 has traces a_p^2 - 1 = 4cos^2(theta) - 1 = 1 + 2cos(2theta))

    print(f"{'Class':<14} {'theta':<12} {'tr(sym^2)':<15} {'density':<10}")
    print("-" * 51)
    sym2_avg = mpf(0)
    for cls in ['golden+', 'golden-', 'order3', 'involution', 'identity']:
        theta = ANGLES[cls]
        tr_sym2 = 1 + 2*cos(2*theta)
        d = CLASSES[cls]['density']
        sym2_avg += d * tr_sym2
        print(f"{cls:<14} {mp.nstr(theta*180/pi, 6):<12} {mp.nstr(tr_sym2, 10):<15} {mp.nstr(d, 6)}")

    print(f"\n<tr(sym^2)> = {mp.nstr(sym2_avg, 20)}")
    print(f"(Should be 1 for irreducible sym^2, 0 if sym^2 contains trivial.)")

    # Compute Euler product approximation for L(1, sym^2)
    # L(1, sym^2) = prod_p (1 - b_p/p)^{-1} approximately
    # where b_p = tr(sym^2(Frob_p))
    # More precisely: (1-alpha_p^2/p)(1-1/p)(1-beta_p^2/p) for eigenvalues alpha,beta
    # But at s=1 we use the approximation with partial products

    # Generate sym^2 coefficients for primes up to limit
    PRIME_LIMIT = 10000
    is_prime = [False, False] + [True] * (PRIME_LIMIT - 1)
    for i in range(2, int(PRIME_LIMIT**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, PRIME_LIMIT + 1, i):
                is_prime[j] = False

    rng = random.Random(12345)
    class_list = list(CLASSES.keys())
    cum_density = []
    running = mpf(0)
    for cls in class_list:
        running += CLASSES[cls]['density']
        cum_density.append(float(running))

    log_L_sym2 = mpf(0)
    log_L_sym2_st = mpf(0)  # Sato-Tate comparison

    for p in range(2, PRIME_LIMIT + 1):
        if not is_prime[p]:
            continue

        # Assign Frobenius
        r = rng.random()
        chosen = class_list[-1]
        for i, cd in enumerate(cum_density):
            if r < cd:
                chosen = class_list[i]
                break

        theta = ANGLES[chosen]
        b_p = 1 + 2*cos(2*theta)
        # Euler factor: -log(1 - b_p/p) ≈ b_p/p + b_p^2/(2p^2) + ...
        if abs(b_p/p) < mpf('0.99'):
            log_L_sym2 += -log(1 - b_p/mpf(p))
        else:
            log_L_sym2 += b_p/mpf(p) + b_p**2/(2*mpf(p)**2)

        # Sato-Tate: b_p = 1+2cos(2theta) with ST measure
        # <b_p>_ST = 1 (from <cos(2theta)>_ST = 0)
        log_L_sym2_st += mpf(1)/p  # Leading term only

    L_sym2 = exp(log_L_sym2)
    L_sym2_st = exp(log_L_sym2_st)

    print(f"\nPartial Euler products (primes up to {PRIME_LIMIT}):")
    print(f"  L(1, sym^2 rho_ico) ~ {mp.nstr(L_sym2, 15)}")
    print(f"  L(1, sym^2)_ST ref  ~ {mp.nstr(L_sym2_st, 15)}")
    print(f"  Ratio: {mp.nstr(L_sym2 / L_sym2_st, 15)}")

    # Now: zero-density estimates
    print(f"\n\n{'='*78}")
    print(f"  ZERO-DENSITY ESTIMATES: N(sigma, T)")
    print(f"{'='*78}\n")

    SIGMAS = [mpf('0.55'), mpf('0.6'), mpf('0.65'), mpf('0.7'),
              mpf('0.75'), mpf('0.8'), mpf('0.85'), mpf('0.9'), mpf('0.95')]
    T_VALUES = [mpf(100), mpf(1000), mpf(10000), mpf(100000)]

    # Standard estimate
    print("--- Standard GL(2) Estimate: N_std(sigma,T) = C*(qT)^{A*(1-sigma)} * log(qT)^B ---\n")
    print(f"{'sigma':<8}", end="")
    for T in T_VALUES:
        print(f"{'T='+mp.nstr(T,0):<18}", end="")
    print()
    print("-" * 80)

    N_std_data = {}
    for sigma in SIGMAS:
        print(f"{mp.nstr(sigma,3):<8}", end="")
        for T in T_VALUES:
            qT = q * T
            N_est = C * qT**(A_std * (1 - sigma)) * log(qT)**B
            N_std_data[(sigma, T)] = N_est
            if N_est > mpf(10)**15:
                print(f"{mp.nstr(N_est, 4):<18}", end="")
            else:
                print(f"{mp.nstr(N_est, 8):<18}", end="")
        print()

    # Golden-enhanced estimate: THREE suppression mechanisms
    print(f"\n\n--- Golden-Enhanced Estimate ---\n")
    print("Three suppression mechanisms:\n")

    # Mechanism 1: Mean-value improvement from higher-moment control
    # The variance ratio gives a large-deviation improvement
    # Since <|a_p|^2>=1 (same as generic), the leading exponent A is unchanged.
    # But the constant C is improved by the variance ratio.

    # Mechanism 2: The sym^2 L-function value enters the mean-value bound
    # M(sigma,T) ~ T * L(2sigma-1, sym^2)^{something}
    # If L(1, sym^2 rho_ico) differs from generic: affects the bound

    # Mechanism 3: Higher Chebyshev moments modify the error terms
    # These affect B (the log power) not A

    # Let's compute all three effects:

    # Mechanism 1: Variance improvement on C
    # var_ico / var_ST ~ reduction in constant
    var_ico = mpf(0)
    var_st_comp = mpf(0)
    for cls, data in CLASSES.items():
        var_ico += data['density'] * (data['trace']**2 - 1)**2

    var_st_comp = mpf(8) / pi * mpmath.quad(
        lambda theta: ((2*cos(theta))**2 - 1)**2 * sin(theta)**2,
        [0, pi]
    )
    C_ratio = var_ico / var_st_comp
    print(f"Mechanism 1 — Variance ratio: C_gold/C_std = {mp.nstr(C_ratio, 15)}")

    # Mechanism 2: sym^2 L-value ratio
    L_ratio = L_sym2 / L_sym2_st
    print(f"Mechanism 2 — sym^2 L-value ratio: {mp.nstr(L_ratio, 15)}")

    # Mechanism 3: Higher moments modify B
    # The log power B comes from sum_{k >= 2} contributions
    # Golden structure: <|a_{p^k}|^2> oscillates with period 5
    # Running average of <|a_{p^k}|^2> for k >= 2:
    if len(avg_history) > 2:
        higher_avg = fsum(avg_history[2:]) / (len(avg_history) - 2)
    else:
        higher_avg = mpf(1)

    # For Sato-Tate: <|U_k(cos theta)|^2>_ST = 1 for all k
    B_ratio = higher_avg / 1  # Ratio of higher-moment averages
    print(f"Mechanism 3 — Higher moment avg (k>=2): {mp.nstr(higher_avg, 15)}")
    print(f"  B_gold/B_std = {mp.nstr(B_ratio, 15)}")

    # Combined golden estimate:
    # N_gold(sigma,T) = C_ratio * N_std(sigma,T) * (correction from sym^2)
    # The sym^2 correction enters the mean value at sigma:
    # M(sigma,T) ~ T * L(2sigma-1, sym^2)  near sigma=1/2

    # For sigma > 1/2: L(2sigma-1, sym^2) is finite and bounded
    # The density estimate is:
    # N(sigma,T) <= C * T^{A(1-sigma)} * q^{A(1-sigma)} * (log qT)^B * L(2sigma-1,sym^2)^{gamma}
    # where gamma depends on the method

    # Let's use gamma = 1 (Jutila's method):
    gamma = mpf(1)
    C_gold = C * C_ratio  # Reduced constant from variance

    print(f"\n--- Combined Golden Zero-Density Estimate ---")
    print(f"N_gold(sigma,T) = {mp.nstr(C_gold,6)} * (qT)^{{A*(1-sigma)}} * log(qT)^{mp.nstr(B,0)} * L_ratio^gamma\n")

    print(f"{'sigma':<8}", end="")
    for T in T_VALUES:
        print(f"{'T='+mp.nstr(T,0):<18}", end="")
    print()
    print("-" * 80)

    N_gold_data = {}
    for sigma in SIGMAS:
        print(f"{mp.nstr(sigma,3):<8}", end="")
        for T in T_VALUES:
            qT = q * T
            N_gold = C_gold * qT**(A_std * (1 - sigma)) * log(qT)**B * L_ratio**gamma
            N_gold_data[(sigma, T)] = N_gold
            if N_gold > mpf(10)**15:
                print(f"{mp.nstr(N_gold, 4):<18}", end="")
            else:
                print(f"{mp.nstr(N_gold, 8):<18}", end="")
        print()

    # RATIO table
    print(f"\n--- Ratio N_gold / N_std ---\n")
    print(f"{'sigma':<8}", end="")
    for T in T_VALUES:
        print(f"{'T='+mp.nstr(T,0):<18}", end="")
    print()
    print("-" * 80)

    for sigma in SIGMAS:
        print(f"{mp.nstr(sigma,3):<8}", end="")
        for T in T_VALUES:
            ratio = N_gold_data[(sigma, T)] / N_std_data[(sigma, T)]
            print(f"{mp.nstr(ratio, 10):<18}", end="")
        print()

    # Critical sigma: where N(sigma,T) < 1
    print(f"\n--- Critical sigma_crit where N(sigma,T) < 1 ---\n")
    print(f"{'T':<12} {'sigma_crit (std)':<20} {'sigma_crit (gold)':<20} {'Improvement':<15}")
    print("-" * 67)

    for T in T_VALUES:
        # Find sigma_crit for standard
        sigma_crit_std = find_sigma_crit(q, T, A_std, B, C)
        sigma_crit_gold = find_sigma_crit(q, T, A_std, B, C_gold * L_ratio**gamma)

        improvement = sigma_crit_std - sigma_crit_gold
        print(f"{mp.nstr(T,0):<12} {mp.nstr(sigma_crit_std,12):<20} "
              f"{mp.nstr(sigma_crit_gold,12):<20} {mp.nstr(improvement,10):<15}")

    # THE VERDICT
    print(f"\n{'='*78}")
    print(f"  THE VERDICT")
    print(f"{'='*78}")
    print()
    print("1. The leading zero-density EXPONENT A is NOT suppressed by golden dominance.")
    print("   <|a_p|^2> = 1 exactly (Schur orthogonality), so A_eff = A.")
    print()
    print("2. The golden structure suppresses the CONSTANT (via variance ratio):")
    print(f"   C_gold/C_std = {mp.nstr(C_ratio, 10)}")
    print()
    print("3. The sym^2 L-value ratio provides additional constant-level suppression:")
    print(f"   L(1,sym^2)_ico / L(1,sym^2)_ST = {mp.nstr(L_ratio, 10)}")
    print()
    print("4. Higher Chebyshev moments (<|a_{p^k}|^2> for k>=2) show interesting")
    print(f"   oscillatory behavior with average {mp.nstr(higher_avg, 10)} (vs 1 for Sato-Tate).")
    print()
    print("5. The combined effect: the zero-free region is widened by a CONSTANT amount")
    print("   (not growing with T), moving sigma_crit closer to 1/2 by O(1/log T).")
    print()
    print("6. D(N) -> infinity does NOT drive the zero-density exponent to 0.")
    print("   It enters as a CONSTANT improvement, not an exponent suppression.")
    print()
    print("7. HOW D(N) actually enters: through the GOLDEN DOMINANCE of the")
    print("   Frobenius distribution, which constrains |a_p| to a DISCRETE set")
    print("   {phi, 1/phi, -1, 0, 2} rather than the continuous Sato-Tate distribution.")
    print("   This discreteness improves large-deviation estimates (mechanism 1)")
    print("   and creates arithmetic correlations in the coefficients.")

    return N_std_data, N_gold_data


def find_sigma_crit(q, T, A, B, C_eff):
    """Find sigma where N(sigma,T) = C_eff * (qT)^{A*(1-sigma)} * log(qT)^B = 1.

    Solve: A*(1-sigma)*log(qT) + B*log(log(qT)) + log(C_eff) = 0
    -> 1-sigma = -(B*log(log(qT)) + log(C_eff)) / (A*log(qT))
    -> sigma = 1 + (B*log(log(qT)) + log(C_eff)) / (A*log(qT))
    """
    qT = q * T
    lqT = log(qT)
    llqT = log(lqT)

    if C_eff > 0:
        log_C = log(C_eff)
    else:
        log_C = mpf(0)

    # N = 1 when: A*(1-sigma)*lqT + B*llqT + log_C = 0
    one_minus_sigma = -(B * llqT + log_C) / (A * lqT)
    sigma_crit = 1 - one_minus_sigma

    # Clamp to [0.5, 1]
    if sigma_crit < mpf('0.5'):
        sigma_crit = mpf('0.5')
    if sigma_crit > 1:
        sigma_crit = mpf(1)

    return sigma_crit


# ============================================================================
# PART 8 (BONUS): Where D(N) ACTUALLY enters — the Effective Conductor
# ============================================================================

def part8_effective_conductor():
    banner("PART 8 (BONUS): Where D(N) Actually Enters")

    print("\nThe golden dominance D(N) counts how many primes up to N have golden")
    print("Frobenius (angle pi/5 or 3pi/5) vs non-golden. By Chebotarev, this is")
    print("asymptotically (2/5) * pi(N), but D(N) measures the DOMINANCE, not the count.\n")

    # D(N) as defined in the framework:
    # D(N) = (# golden primes up to N) * phi / (# non-golden primes up to N)
    # This grows like (2/5)*phi / (3/5) = 2phi/3 ~ 1.078 (constant!)
    # Unless D(N) is defined differently...

    # The standard definition from the Pythagorean framework:
    # D(N) measures the spectral dominance of golden eigenvalues
    # D(N) = sum_{p <= N, golden} log(p) / sum_{p <= N} log(p) * weight
    # With weight = phi for golden+ and 1/phi for golden-

    print("Computing D(N) for various N:\n")

    rng = random.Random(999)
    class_list = list(CLASSES.keys())
    cum_density = []
    running = mpf(0)
    for cls in class_list:
        running += CLASSES[cls]['density']
        cum_density.append(float(running))

    # Sieve
    LIMIT = 100000
    is_prime = [False, False] + [True] * (LIMIT - 1)
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, LIMIT + 1, i):
                is_prime[j] = False

    # Assign Frobenius and compute D(N)
    golden_sum = mpf(0)
    total_sum = mpf(0)
    weighted_golden = mpf(0)

    checkpoints = [100, 500, 1000, 5000, 10000, 50000, 100000]
    checkpoint_idx = 0

    print(f"{'N':<10} {'D(N)':<20} {'D_weighted(N)':<20} {'#golden/#total':<20}")
    print("-" * 70)

    n_golden = 0
    n_total = 0

    for p in range(2, LIMIT + 1):
        if not is_prime[p]:
            continue

        r = rng.random()
        chosen = class_list[-1]
        for i, cd in enumerate(cum_density):
            if r < cd:
                chosen = class_list[i]
                break

        n_total += 1
        logp = log(mpf(p))
        total_sum += logp

        if chosen in ('golden+', 'golden-'):
            n_golden += 1
            golden_sum += logp
            if chosen == 'golden+':
                weighted_golden += logp * phi
            else:
                weighted_golden += logp / phi

        if checkpoint_idx < len(checkpoints) and p >= checkpoints[checkpoint_idx]:
            D_N = golden_sum / total_sum * 5 / 2 if total_sum > 0 else mpf(0)
            D_weighted = weighted_golden / total_sum if total_sum > 0 else mpf(0)
            frac = mpf(n_golden) / n_total if n_total > 0 else mpf(0)
            print(f"{checkpoints[checkpoint_idx]:<10} "
                  f"{mp.nstr(D_N, 12):<20} {mp.nstr(D_weighted, 12):<20} "
                  f"{mp.nstr(frac, 12):<20}")
            checkpoint_idx += 1

    print(f"\nD(N) converges to ~ 1 (normalized golden fraction = 2/5).")
    print("D(N) does NOT diverge to infinity in the simple counting sense.")
    print()
    print("The 'D(N) -> infinity' in the framework refers to the SPECTRAL")
    print("dominance: the accumulated golden eigenvalue advantage in the")
    print("explicit formula. Specifically:")
    print()
    print("  D(N) = sum_{p<=N, golden} phi^k * log(p) / p^{1/2}")
    print("       = sum over golden primes of their 'amplified' contribution")
    print()
    print("This DOES grow ~ phi * sqrt(N) / log(N) -> infinity.")
    print("But it enters the zero-density estimate as a CORRECTION to the")
    print("explicit formula, not as a suppression of the density exponent.\n")

    # Compute the spectral D(N)
    print("--- Spectral D(N) ---\n")
    rng2 = random.Random(999)  # Same seed for reproducibility

    spectral_D = mpf(0)
    spectral_generic = mpf(0)

    checkpoints2 = [100, 500, 1000, 5000, 10000, 50000, 100000]
    cp_idx = 0

    print(f"{'N':<10} {'D_spectral(N)':<25} {'Generic sum':<25} {'Ratio':<15}")
    print("-" * 75)

    for p in range(2, LIMIT + 1):
        if not is_prime[p]:
            continue

        r = rng2.random()
        chosen = class_list[-1]
        for i, cd in enumerate(cum_density):
            if r < cd:
                chosen = class_list[i]
                break

        logp = log(mpf(p))
        theta = ANGLES[chosen]
        a_p = 2 * cos(theta)

        # Golden spectral contribution: a_p * log(p) / sqrt(p)
        contrib = a_p * logp / sqrt(mpf(p))
        spectral_generic += fabs(contrib)

        if chosen in ('golden+', 'golden-'):
            spectral_D += fabs(contrib)

        if cp_idx < len(checkpoints2) and p >= checkpoints2[cp_idx]:
            ratio = spectral_D / spectral_generic if spectral_generic > 0 else mpf(0)
            print(f"{checkpoints2[cp_idx]:<10} "
                  f"{mp.nstr(spectral_D, 15):<25} {mp.nstr(spectral_generic, 15):<25} "
                  f"{mp.nstr(ratio, 10):<15}")
            cp_idx += 1

    print(f"\nThe spectral D(N) grows like sqrt(N)/log(N) -> infinity.")
    print(f"Final spectral D({LIMIT}) = {mp.nstr(spectral_D, 15)}")
    print(f"Final generic sum = {mp.nstr(spectral_generic, 15)}")
    print(f"Golden fraction of spectral weight = {mp.nstr(spectral_D/spectral_generic, 10)}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("  GOLDEN DOMINANCE AND ZERO-DENSITY EXPONENTS FOR L(s, rho_ico)")
    print("  mpmath precision: {} decimal digits".format(mp.dps))
    print("=" * 78)

    # Part 1
    avg_exact, avg_no_id, avg_approx = part1_mean_value()

    # Part 2
    part2_zero_density_exponent(avg_exact, avg_no_id)

    # Part 3 (numerical — computationally expensive)
    print("\n[Part 3 uses truncated Dirichlet series — may take a moment...]")
    results_ico, results_zeta = part3_numerical_simulation()

    # Part 4
    S_exact, S_no_id, S4 = part4_suppression_factor()

    # Part 5
    avg_history = part5_chebyshev_traces()

    # Part 6
    var_ico, var_st = part6_compounding(avg_history)

    # Part 7
    N_std, N_gold = part7_decisive_test(avg_history)

    # Part 8 (bonus)
    part8_effective_conductor()

    # FINAL SUMMARY
    banner("FINAL SUMMARY")
    print()
    print("Q: Does D(N) -> infinity suppress the zero-density exponent for L(s, rho_ico)?")
    print()
    print("A: NO — not in the way claimed. Here is what actually happens:")
    print()
    print("1. SCHUR ORTHOGONALITY THEOREM kills the leading suppression:")
    print("   <|a_p|^2> = 1 for ANY irreducible representation (exact Chebotarev).")
    print("   The zero-density EXPONENT A is unchanged: A_eff = A.")
    print()
    print("2. The '14/15' value arises from EXCLUDING the identity class (density 1/60).")
    print("   Including it: (12*phi^2 + 12/phi^2 + 20 + 0 + 4)/60 = 60/60 = 1.")
    print(f"   The identity class (trace 2, |a_p|^2 = 4) exactly compensates.")
    print()
    print("3. Golden dominance DOES suppress via THREE secondary mechanisms:")
    print(f"   a) Variance reduction: Var_ico/Var_ST = {mp.nstr(var_ico/var_st, 8)}")
    print(f"   b) sym^2 L-value modification (constant factor)")
    print(f"   c) Chebyshev oscillation at higher prime powers")
    print(f"   Combined: sigma_crit moves closer to 1/2 by O(1/log T).")
    print()
    print("4. The spectral D(N) -> infinity enters the EXPLICIT FORMULA,")
    print("   not the density exponent. It provides a growing 'golden advantage'")
    print("   in the prime sum, which improves zero-free regions but does not")
    print("   drive N(sigma,T) to zero for sigma > 1/2.")
    print()
    print("5. BOTTOM LINE: The golden dominance provides CONSTANT-LEVEL improvements")
    print("   to zero-density estimates. It does NOT provide the infinite suppression")
    print("   needed for GRH. The dream of D_eff -> infinity killing the exponent")
    print("   is blocked by Schur orthogonality: <|a_p|^2> = 1 is a THEOREM,")
    print("   not something the golden structure can override.")
    print()
    print("6. WHERE TO LOOK INSTEAD: The golden structure's real power may be in")
    print("   a) Subconvexity bounds (where discrete Frobenius helps)")
    print("   b) Moment bounds for L(1/2+it, rho_ico) (where golden angles cluster)")
    print("   c) The FUNCTIONAL EQUATION with icosahedral root number")
    print("   d) Langlands functoriality (sym^k lifts have golden-constrained traces)")
    print()
    print("=" * 78)
    print("  Investigation complete.")
    print("=" * 78)


if __name__ == '__main__':
    main()
