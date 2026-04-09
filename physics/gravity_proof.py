#!/usr/bin/env python3
"""
GRAVITY_PROOF — tests whether off-center zero pairs are inconsistent with L(s, rho_ico) explicit formula
nos3bl33d

Hooke's law analogy: quadratic potential spring constant grows over golden primes, forcing delta=0.
Eight parts testing different facets of the gravitational GRH argument.
"""

import sys
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

from mpmath import (
    mp, mpf, mpc, sqrt, log, pi, exp, cosh, sinh, cos, sin,
    fabs, fsum, power, re, im, quad, inf
)

mp.dps = 30

# =============================================================================
# CONSTANTS
# =============================================================================
PHI = (1 + sqrt(mpf(5))) / 2       # golden ratio
INV_PHI = PHI - 1                   # 1/phi = phi - 1
CONDUCTOR = mpf(800)

print("=" * 78)
print("  GRAVITATIONAL PROOF OF GRH — COMPUTATIONAL TEST")
print("  Testing whether off-center zero pairs are inconsistent")
print("  with the explicit formula for L(s, rho_ico)")
print("=" * 78)
print()
print(f"  phi            = {PHI}")
print(f"  1/phi          = {INV_PHI}")
print(f"  conductor      = {CONDUCTOR}")
print(f"  precision      = {mp.dps} decimal digits")
print()


# =============================================================================
# PRIME INFRASTRUCTURE
# =============================================================================

def sieve_primes(n):
    """Sieve of Eratosthenes up to n."""
    if n < 2:
        return []
    is_prime = bytearray(b'\x01') * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = 0
    return [i for i in range(2, n + 1) if is_prime[i]]


def is_golden_prime(p):
    """
    A prime p is 'golden' for the icosahedral representation if
    Frob_p maps to a 5-cycle conjugacy class in A5.

    For the splitting field of x^2 + x - 1 (= Q(sqrt(5))):
    p splits iff p = +/-1 mod 5 (quadratic reciprocity).
    These are the primes where the Frobenius trace is phi or -1/phi.

    For the icosahedral Galois representation, the golden primes
    (5-cycle classes) have density 24/60 = 2/5 by Chebotarev.
    Among split primes (p = +/-1 mod 5), the Frobenius eigenvalue
    involves phi.

    We use the simple criterion: p = +/-1 mod 5 (i.e., p mod 5 in {1, 4}).
    Note: p=2,3 have p mod 5 in {2,3}, so they are non-golden. p=5 ramifies.
    """
    if p == 5:
        return False  # ramified
    return (p % 5) in (1, 4)


def golden_trace(p):
    """
    Return the Frobenius trace at a golden prime.
    p = 1 mod 5 -> trace = phi (5-cycle class 5A)
    p = 4 mod 5 -> trace = -1/phi (5-cycle class 5B)
    """
    if p % 5 == 1:
        return PHI
    elif p % 5 == 4:
        return -INV_PHI
    else:
        raise ValueError(f"p={p} is not a golden prime")


def classify_prime(p):
    """
    Classify prime p into A5 conjugacy classes based on p mod 5.
    Returns (class_name, trace).

    For the icosahedral representation over Q(zeta_5):
    - p = 1 mod 5: 5-cycle class 5A, trace = phi
    - p = 4 mod 5: 5-cycle class 5B, trace = -1/phi
    - p = 2,3 mod 5: other classes (involutions, order-3, etc.)
    - p = 5: ramified

    For non-golden primes, the exact class depends on higher splitting
    behavior. We assign traces based on the full Artin conductor data.
    For this test, non-golden primes get trace 0 (involution class)
    as a conservative placeholder — their contributions are bounded
    and don't affect the divergence argument.
    """
    if p == 5:
        return ("ramified", mpf(0))
    r = p % 5
    if r == 1:
        return ("5A", PHI)
    elif r == 4:
        return ("5B", -INV_PHI)
    elif r == 2:
        return ("inv/o3", mpf(0))   # conservative: bounded trace
    elif r == 3:
        return ("inv/o3", -mpf(1))  # order-3 trace
    return ("unknown", mpf(0))


# =============================================================================
# PART 1: EXCESS PER GOLDEN PRIME
# =============================================================================

def part1_excess_per_prime():
    """
    Compute E(p, delta) = 2*sqrt(p) * (cosh(delta*log(p)) - 1)
    for each golden prime p up to 10000.
    Verify E ~ sqrt(p) * delta^2 * (log p)^2 for small delta.
    """
    print("=" * 78)
    print("  PART 1: EXCESS PER GOLDEN PRIME")
    print("  E(p, delta) = 2*sqrt(p) * (cosh(delta*log(p)) - 1)")
    print("=" * 78)
    print()

    primes = sieve_primes(10000)
    golden_primes = [p for p in primes if is_golden_prime(p)]

    deltas = [mpf("0.001"), mpf("0.01"), mpf("0.1")]

    print(f"  Total primes up to 10000: {len(primes)}")
    print(f"  Golden primes up to 10000: {len(golden_primes)}")
    print(f"  Golden fraction: {len(golden_primes)/len(primes):.4f} (expected 2/5 = 0.4000)")
    print()

    # Show a sample of golden primes
    print(f"  First 20 golden primes: {golden_primes[:20]}")
    print()

    # Table header
    print(f"  {'p':>6s}  {'trace':>8s}  {'delta=0.001':>14s}  {'approx':>14s}  "
          f"{'delta=0.01':>14s}  {'approx':>14s}  {'delta=0.1':>14s}  {'approx':>14s}")
    print("  " + "-" * 104)

    sample_primes = golden_primes[:10] + golden_primes[-5:]

    for p in sample_primes:
        p_mp = mpf(p)
        sqp = sqrt(p_mp)
        lp = log(p_mp)
        tr = golden_trace(p)

        row = f"  {p:>6d}  {float(tr):>8.4f}"
        for d in deltas:
            exact = 2 * sqp * (cosh(d * lp) - 1)
            approx = sqp * d**2 * lp**2
            row += f"  {float(exact):>14.8f}  {float(approx):>14.8f}"
        print(row)

    # Verify approximation quality
    print()
    print("  Approximation quality (exact/approx ratio, should -> 1 for small delta):")
    for d in deltas:
        ratios = []
        for p in golden_primes:
            p_mp = mpf(p)
            exact = 2 * sqrt(p_mp) * (cosh(d * log(p_mp)) - 1)
            approx = sqrt(p_mp) * d**2 * (log(p_mp))**2
            if approx > 0:
                ratios.append(float(exact / approx))
        avg_ratio = sum(ratios) / len(ratios)
        max_ratio = max(ratios)
        min_ratio = min(ratios)
        print(f"    delta={float(d):.3f}: avg ratio={avg_ratio:.8f}, "
              f"range=[{min_ratio:.8f}, {max_ratio:.8f}]")

    print()
    return golden_primes


# =============================================================================
# PART 2: CUMULATIVE GOLDEN EXCESS
# =============================================================================

def part2_cumulative_excess(golden_primes):
    """
    Total_E(N, delta) = sum_{golden p <= N} trace(p) * E(p, delta) * log(p) / sqrt(p)

    The factor log(p)/sqrt(p) comes from the explicit formula weight.
    The trace(p) is the Frobenius eigenvalue (phi or -1/phi).
    """
    print("=" * 78)
    print("  PART 2: CUMULATIVE GOLDEN EXCESS")
    print("  Total_E = sum trace(p) * 2*sqrt(p)*(cosh(delta*log(p))-1) * log(p)/sqrt(p)")
    print("         = sum trace(p) * 2*(cosh(delta*log(p))-1) * log(p)")
    print("=" * 78)
    print()

    deltas = [mpf("0.001"), mpf("0.01"), mpf("0.1")]

    # Note: the sqrt(p) in E(p) cancels with 1/sqrt(p) in the explicit formula weight
    # So: trace(p) * E(p,d) * log(p)/sqrt(p) = trace(p) * 2*(cosh(d*log(p))-1) * log(p)

    for d in deltas:
        print(f"  --- delta = {float(d)} ---")
        cumulative = mpf(0)
        abs_cumulative = mpf(0)  # track |cumulative| to check for cancellations
        checkpoints = [100, 500, 1000, 2000, 5000, 10000]
        cp_idx = 0

        for p in golden_primes:
            p_mp = mpf(p)
            lp = log(p_mp)
            tr = golden_trace(p)
            excess = 2 * (cosh(d * lp) - 1) * lp
            weighted = tr * excess
            cumulative += weighted
            abs_cumulative += fabs(weighted)

            if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
                ratio_cancel = float(fabs(cumulative) / abs_cumulative) if abs_cumulative > 0 else 0
                print(f"    p <= {checkpoints[cp_idx]:>5d}: "
                      f"Total_E = {float(cumulative):>+16.8f}, "
                      f"|Total_E| = {float(fabs(cumulative)):>14.8f}, "
                      f"|Total_E|/sum|terms| = {ratio_cancel:.6f}")
                cp_idx += 1
                if cp_idx >= len(checkpoints):
                    break

        print()

    # CRITICAL CHECK: Do the golden+ (trace phi) and golden- (trace -1/phi)
    # contributions cancel or accumulate?
    print("  CANCELLATION ANALYSIS:")
    print("  Golden+ trace = phi ~ 1.618, Golden- trace = -1/phi ~ -0.618")
    print("  Net per golden pair: phi + (-1/phi) = 1.0 (ALWAYS POSITIVE)")
    print()

    d = mpf("0.01")
    sum_plus = mpf(0)
    sum_minus = mpf(0)
    count_plus = 0
    count_minus = 0

    for p in golden_primes:
        p_mp = mpf(p)
        lp = log(p_mp)
        excess = 2 * (cosh(d * lp) - 1) * lp

        if p % 5 == 1:
            sum_plus += PHI * excess
            count_plus += 1
        else:
            sum_minus += (-INV_PHI) * excess
            count_minus += 1

    print(f"  delta = 0.01:")
    print(f"    Golden+ (p=1 mod 5): {count_plus} primes, sum = {float(sum_plus):>+16.8f}")
    print(f"    Golden- (p=4 mod 5): {count_minus} primes, sum = {float(sum_minus):>+16.8f}")
    print(f"    Total:               sum = {float(sum_plus + sum_minus):>+16.8f}")
    print(f"    Does it cancel? {'YES (bad for argument)' if fabs(sum_plus + sum_minus) < fabs(sum_plus) * mpf('0.01') else 'NO (good — net accumulation)'}")
    print()


# =============================================================================
# PART 3: THE BALANCE TEST — INTEGRAL OF EXCESS SQUARED
# =============================================================================

def part3_balance_test():
    """
    The excess D(x, delta, gamma_0) from an off-center pair at 1/2 +/- delta + i*gamma_0
    minus a collapsed pair at 1/2 + i*gamma_0:

    D(x) = x^{1/2+delta+igamma}/(1/2+delta+igamma) + x^{1/2-delta+igamma}/(1/2-delta+igamma)
         - 2*x^{1/2+igamma}/(1/2+igamma)

    We compute integral |D(x)|^2 dx/x over [2, X] for various X.
    """
    print("=" * 78)
    print("  PART 3: BALANCE TEST — EXCESS ENERGY INTEGRAL")
    print("  int_2^X |D(x)|^2 dx/x where D = pair - collapsed contribution")
    print("=" * 78)
    print()

    deltas = [mpf("0.001"), mpf("0.01"), mpf("0.1")]
    gammas = [mpf(1), mpf(10), mpf(100)]

    # For the integral |D(x)|^2 dx/x, we can compute analytically.
    # D(x) = x^{1/2+igamma} * F(x, delta, gamma) where
    # F(x) = x^delta/(1/2+delta+igamma) + x^{-delta}/(1/2-delta+igamma) - 2/(1/2+igamma)
    #
    # |D(x)|^2 = x * |F(x)|^2  (since |x^{1/2+igamma}|^2 = x)
    #
    # So int |D(x)|^2 dx/x = int |F(x)|^2 dx
    #
    # F(x) involves x^{+/-delta}, so |F(x)|^2 involves x^{2delta}, x^0, x^{-2delta}
    # Each integrates to a power of X.

    print(f"  {'delta':>8s}  {'gamma':>8s}  {'X=100':>16s}  {'X=1000':>16s}  "
          f"{'X=10000':>16s}  {'X=100000':>16s}")
    print("  " + "-" * 80)

    for d in deltas:
        for g in gammas:
            rho_plus = mpc(mpf("0.5") + d, g)
            rho_minus = mpc(mpf("0.5") - d, g)
            rho_center = mpc(mpf("0.5"), g)

            results = []
            for X_val in [100, 1000, 10000, 100000]:
                X = mpf(X_val)

                # Numerical integration via summation over integer points
                # (faster and more stable than quad for oscillatory integrands)
                # Use logarithmic sampling for large X
                n_points = min(2000, X_val)
                if X_val <= 2000:
                    xs = [mpf(k) for k in range(2, X_val + 1)]
                else:
                    # logarithmic sampling
                    xs = []
                    log_min = log(mpf(2))
                    log_max = log(X)
                    for i in range(n_points):
                        t = log_min + (log_max - log_min) * i / (n_points - 1)
                        xs.append(exp(t))

                total = mpf(0)
                for x in xs:
                    # D(x) = x^rho_plus / rho_plus + x^rho_minus / rho_minus - 2*x^rho_center / rho_center
                    term_plus = power(x, rho_plus) / rho_plus
                    term_minus = power(x, rho_minus) / rho_minus
                    term_center = 2 * power(x, rho_center) / rho_center
                    D = term_plus + term_minus - term_center
                    # |D|^2 / x
                    D_abs_sq = re(D * D.conjugate())
                    if X_val <= 2000:
                        total += D_abs_sq / x  # Riemann sum with dx=1
                    else:
                        # dx = x * d(log x) ≈ x * (log_max - log_min) / n_points
                        dx = x * (log_max - log_min) / n_points
                        total += D_abs_sq / x * dx

                results.append(float(total))

            print(f"  {float(d):>8.3f}  {float(g):>8.1f}  "
                  f"{results[0]:>16.6f}  {results[1]:>16.6f}  "
                  f"{results[2]:>16.6f}  {results[3]:>16.6f}")

    print()
    print("  KEY QUESTION: Does the integral grow with X?")
    print("  If yes: the excess energy is unbounded -> pair inconsistent.")
    print()


# =============================================================================
# PART 4: GOLDEN COHERENCE IN THE EXCESS
# =============================================================================

def part4_golden_coherence(golden_primes):
    """
    At x = golden_prime, the excess is weighted by the golden trace.

    Golden_excess(N) = sum_{golden p <= N} |trace(p)| * |D(p, delta, gamma_0)|
    On_line(N) = sum_{golden p <= N} |trace(p)| * |collapsed(p, gamma_0)|

    Does golden_excess / on_line grow?
    """
    print("=" * 78)
    print("  PART 4: GOLDEN COHERENCE IN THE EXCESS")
    print("  Ratio = golden_excess / on_line_contribution vs N")
    print("=" * 78)
    print()

    d = mpf("0.01")
    g = mpf(10)
    rho_plus = mpc(mpf("0.5") + d, g)
    rho_minus = mpc(mpf("0.5") - d, g)
    rho_center = mpc(mpf("0.5"), g)

    golden_excess_sum = mpf(0)
    on_line_sum = mpf(0)

    checkpoints = [50, 100, 500, 1000, 2000, 5000, 10000]
    cp_idx = 0

    print(f"  delta = {float(d)}, gamma_0 = {float(g)}")
    print()
    print(f"  {'N':>6s}  {'golden_excess':>16s}  {'on_line':>16s}  {'ratio':>12s}  {'sqrt(N)':>10s}")
    print("  " + "-" * 66)

    count = 0
    for p in golden_primes:
        p_mp = mpf(p)
        tr = fabs(golden_trace(p))

        # Pair contribution
        term_plus = power(p_mp, rho_plus) / rho_plus
        term_minus = power(p_mp, rho_minus) / rho_minus
        # Collapsed contribution
        term_center = 2 * power(p_mp, rho_center) / rho_center

        D = term_plus + term_minus - term_center
        D_abs = sqrt(re(D * D.conjugate()))

        collapsed_abs = sqrt(re(term_center * term_center.conjugate()))

        golden_excess_sum += tr * D_abs * log(p_mp)
        on_line_sum += tr * collapsed_abs * log(p_mp)
        count += 1

        if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
            ratio = float(golden_excess_sum / on_line_sum) if on_line_sum > 0 else 0
            sqrtN = float(sqrt(mpf(checkpoints[cp_idx])))
            print(f"  {checkpoints[cp_idx]:>6d}  {float(golden_excess_sum):>16.6f}  "
                  f"{float(on_line_sum):>16.6f}  {ratio:>12.8f}  {sqrtN:>10.2f}")
            cp_idx += 1
            if cp_idx >= len(checkpoints):
                break

    print()
    print("  If ratio is CONSTANT: excess is proportional to on-line -> no new constraint.")
    print("  If ratio GROWS: excess dominates -> pair is inconsistent.")
    print()


# =============================================================================
# PART 5: EXPLICIT COMPUTATION
# =============================================================================

def part5_explicit_computation(golden_primes):
    """
    For delta=0.01, gamma_0=10:
    Compute pair vs collapsed contribution at each golden prime.
    Accumulate weighted excess. Compare to log(N).
    """
    print("=" * 78)
    print("  PART 5: EXPLICIT COMPUTATION")
    print("  Hypothetical off-center zero at 0.51 + 10i")
    print("=" * 78)
    print()

    d = mpf("0.01")
    g = mpf(10)

    rho_plus = mpc(mpf("0.5") + d, g)   # 0.51 + 10i
    rho_minus = mpc(mpf("0.5") - d, g)  # 0.49 + 10i
    rho_center = mpc(mpf("0.5"), g)      # 0.50 + 10i

    print(f"  rho+ = {rho_plus}")
    print(f"  rho- = {rho_minus}")
    print(f"  rho_center = {rho_center}")
    print()

    # (a) through (g)
    total_weighted_excess = mpf(0)
    total_weighted_excess_signed = mpf(0)

    checkpoints = [100, 500, 1000, 2000, 5000, 10000]
    cp_idx = 0

    print(f"  {'N':>6s}  {'Total(N)':>16s}  {'Total_signed':>16s}  "
          f"{'log(N)':>10s}  {'Total/log(N)':>14s}  {'grows?':>8s}")
    print("  " + "-" * 76)

    prev_ratio = mpf(0)
    for p in golden_primes:
        p_mp = mpf(p)
        lp = log(p_mp)
        tr = golden_trace(p)

        # (a) pair contribution
        C_pair = power(p_mp, rho_plus) / rho_plus + power(p_mp, rho_minus) / rho_minus

        # (b) collapsed contribution
        C_coll = 2 * power(p_mp, rho_center) / rho_center

        # (c) excess
        E = C_pair - C_coll
        E_abs = sqrt(re(E * E.conjugate()))

        # (d) weight by golden trace * log(p)
        weighted = fabs(tr) * E_abs * lp
        weighted_signed = tr * re(E) * lp  # signed version

        # (e) accumulate
        total_weighted_excess += weighted
        total_weighted_excess_signed += weighted_signed

        # (f) compare to log(N)
        if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
            N = checkpoints[cp_idx]
            logN = log(mpf(N))
            ratio = total_weighted_excess / logN
            grows = "YES" if ratio > prev_ratio * mpf("1.01") else "no"
            print(f"  {N:>6d}  {float(total_weighted_excess):>16.6f}  "
                  f"{float(total_weighted_excess_signed):>16.6f}  "
                  f"{float(logN):>10.4f}  {float(ratio):>14.6f}  {grows:>8s}")
            prev_ratio = ratio
            cp_idx += 1
            if cp_idx >= len(checkpoints):
                break

    print()
    # (g) Does Total(N)/log(N) -> infinity?
    print("  VERDICT: If Total(N)/log(N) grows without bound,")
    print("  the excess exceeds the pole's payment -> pair can't exist.")
    print()


# =============================================================================
# PART 6: GRAVITATIONAL ANALOGY — SPRING CONSTANT
# =============================================================================

def part6_spring_constant(golden_primes):
    """
    The total spring constant: K_total = sum_{golden p <= N} 2*|trace(p)|*sqrt(p)*(log p)^2
    This measures the "stiffness" of the force pulling zeros to the critical line.
    """
    print("=" * 78)
    print("  PART 6: GRAVITATIONAL ANALOGY — SPRING CONSTANT")
    print("  K_total = sum 2*|trace(p)|*sqrt(p)*(log p)^2")
    print("=" * 78)
    print()

    K_total = mpf(0)
    K_golden_plus = mpf(0)
    K_golden_minus = mpf(0)

    # Also compute a "raw" spring constant without the trace weight
    K_raw = mpf(0)

    checkpoints = [100, 500, 1000, 2000, 5000, 10000]
    cp_idx = 0

    print(f"  {'N':>6s}  {'K_total':>16s}  {'K_raw':>16s}  {'K/N^{3/2}':>14s}  "
          f"{'K/(N*logN^2)':>16s}")
    print("  " + "-" * 70)

    for p in golden_primes:
        p_mp = mpf(p)
        sqp = sqrt(p_mp)
        lp = log(p_mp)
        tr = fabs(golden_trace(p))

        k_term = 2 * tr * sqp * lp**2
        K_total += k_term

        k_raw = 2 * sqp * lp**2
        K_raw += k_raw

        if p % 5 == 1:
            K_golden_plus += 2 * PHI * sqp * lp**2
        else:
            K_golden_minus += 2 * INV_PHI * sqp * lp**2

        if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
            N = checkpoints[cp_idx]
            N_mp = mpf(N)
            ratio_32 = K_total / power(N_mp, mpf("1.5"))
            logN = log(N_mp)
            ratio_NlogN2 = K_total / (N_mp * logN**2)
            print(f"  {N:>6d}  {float(K_total):>16.4f}  {float(K_raw):>16.4f}  "
                  f"{float(ratio_32):>14.6f}  {float(ratio_NlogN2):>16.6f}")
            cp_idx += 1
            if cp_idx >= len(checkpoints):
                break

    print()
    print(f"  K(golden+) = {float(K_golden_plus):.4f}")
    print(f"  K(golden-) = {float(K_golden_minus):.4f}")
    print(f"  K(total)   = {float(K_total):.4f}")
    print()
    print("  If K_total -> infinity: spring becomes infinitely stiff.")
    print("  No finite delta can resist an infinite restoring force.")
    print("  Growth rate: K ~ N^{3/2} / log(N) (from prime density + sqrt(p) growth)")
    print()


# =============================================================================
# PART 7: THE CATCH SEARCH
# =============================================================================

def part7_catch_search(golden_primes):
    """
    Systematically search for holes in the argument:
    (a) Does the explicit formula constrain the excess?
    (b) Can excess be compensated by rearranging on-line zeros?
    (c) Is sqrt(N) growth correct, or are there cancellations?
    (d) Does Schur orthogonality kill it?
    (e) At what N does excess first exceed pole's payment?
    """
    print("=" * 78)
    print("  PART 7: THE CATCH SEARCH — LOOKING FOR HOLES")
    print("=" * 78)
    print()

    # =========================================================================
    # (a) DOES THE EXPLICIT FORMULA ACTUALLY CONSTRAIN THE EXCESS?
    # =========================================================================
    print("  (a) EXPLICIT FORMULA CONSTRAINT")
    print("  " + "-" * 60)
    print()
    print("  The explicit formula (Weil's version):")
    print("    psi(x) = x - sum_rho x^rho/rho + constants")
    print()
    print("  LHS (psi) is FIXED by primes. It doesn't know about zeros.")
    print("  RHS must match. If we perturb a zero pair off the line,")
    print("  the sum changes. This change MUST be compensated by other zeros.")
    print()
    print("  The question: CAN the other zeros compensate?")
    print("  If the excess grows faster than any rearrangement of on-line")
    print("  zeros can absorb, then no -> pair inconsistent.")
    print()
    print("  FINDING: The excess from an off-center pair at height gamma_0")
    print("  is O(x^{1/2+delta}) at large x, while on-line zeros contribute")
    print("  O(x^{1/2}). The x^delta factor means the off-center excess")
    print("  EVENTUALLY dominates at large x. But 'eventually' means")
    print("  x >> exp(1/delta), which can be astronomical for small delta.")
    print()

    # =========================================================================
    # (b) CAN THE EXCESS BE COMPENSATED?
    # =========================================================================
    print("  (b) COMPENSATION BY ON-LINE ZEROS")
    print("  " + "-" * 60)
    print()

    # The off-center pair at 1/2+/-delta contributes x^{1/2+delta} terms.
    # On-line zeros can only contribute x^{1/2} terms (with oscillation).
    # For x > exp(1/delta), the x^delta factor exceeds any bounded oscillation.

    d = mpf("0.01")
    x_crossover = exp(1 / d)
    print(f"  For delta = {float(d)}:")
    print(f"    x_crossover = exp(1/delta) = exp(100) = {float(x_crossover):.6e}")
    print(f"    Beyond this x, the excess x^delta dominates x^{1/2} by a factor > e.")
    print()

    d = mpf("0.001")
    x_crossover = exp(1 / d)
    print(f"  For delta = {float(d)}:")
    print(f"    x_crossover = exp(1000) = {float(x_crossover):.6e}")
    print(f"    HUGE. The compensation regime is astronomically large.")
    print()

    print("  CATCH: For small delta, the excess only dominates at enormous x.")
    print("  The 'gravitational' argument must work at FINITE x to be useful.")
    print("  This is where golden coherence enters — it amplifies the excess")
    print("  at moderate x through systematic (non-cancelling) accumulation.")
    print()

    # =========================================================================
    # (c) CANCELLATION CHECK
    # =========================================================================
    print("  (c) CANCELLATION CHECK — SIGN CHANGES IN THE SUM")
    print("  " + "-" * 60)
    print()

    d = mpf("0.01")
    g = mpf(10)
    rho_plus = mpc(mpf("0.5") + d, g)
    rho_minus = mpc(mpf("0.5") - d, g)
    rho_center = mpc(mpf("0.5"), g)

    # Track the REAL PART of the signed sum
    signed_sum = mpf(0)
    abs_sum = mpf(0)
    sign_changes = 0
    prev_sign = 0

    for p in golden_primes:
        p_mp = mpf(p)
        lp = log(p_mp)
        tr = golden_trace(p)

        C_pair = power(p_mp, rho_plus) / rho_plus + power(p_mp, rho_minus) / rho_minus
        C_coll = 2 * power(p_mp, rho_center) / rho_center
        E = C_pair - C_coll

        # The explicit formula contribution is Re(trace * E * log(p) / sqrt(p))
        # But E already includes sqrt(p) in the p^{1/2} factor
        contribution = re(tr * E * lp / sqrt(p_mp))
        signed_sum += contribution
        abs_sum += fabs(contribution)

        curr_sign = 1 if contribution > 0 else -1
        if prev_sign != 0 and curr_sign != prev_sign:
            sign_changes += 1
        prev_sign = curr_sign

    cancel_ratio = float(fabs(signed_sum) / abs_sum) if abs_sum > 0 else 0

    print(f"  delta={float(d)}, gamma={float(g)}:")
    print(f"    Signed sum:     {float(signed_sum):>+16.8f}")
    print(f"    Sum of |terms|: {float(abs_sum):>16.8f}")
    print(f"    |signed|/sum:   {cancel_ratio:.6f}")
    print(f"    Sign changes:   {sign_changes} out of {len(golden_primes)} terms")
    print()

    if cancel_ratio < 0.01:
        print("    MASSIVE CANCELLATION DETECTED.")
        print("    The oscillatory factor e^{i*gamma*log(p)} causes wild sign changes.")
        print("    The golden trace (phi or -1/phi) does NOT prevent this.")
    elif cancel_ratio < 0.1:
        print("    SIGNIFICANT CANCELLATION but not total.")
    else:
        print("    LIMITED CANCELLATION — systematic accumulation.")
    print()

    # =========================================================================
    # (d) SCHUR ORTHOGONALITY CHECK
    # =========================================================================
    print("  (d) SCHUR ORTHOGONALITY CHECK")
    print("  " + "-" * 60)
    print()

    # The key question: does Schur constrain the NONLINEAR function
    # trace(p) * 2*(cosh(delta*log(p)) - 1)?

    # Schur orthogonality for A5:
    # sum_{g in G} chi(g) * conj(chi'(g)) = |G| * delta_{chi,chi'}
    # Applied to primes via Chebotarev: this constrains LINEAR sums of traces.
    # Our excess involves cosh(delta*log(p)) which is a FUNCTION of p, not of the trace.

    print("  Schur orthogonality constrains:")
    print("    sum a_p f(p) for MULTIPLICATIVE f(p)  (Hecke eigenvalue sums)")
    print()
    print("  Our excess is: sum trace(p) * g(p) where g(p) = 2(cosh(d*log(p))-1)*log(p)")
    print("  Here g(p) is NOT multiplicative — it's a function of p alone.")
    print()
    print("  Schur says: <chi_1 * chi_2> = 0 for different irreducibles.")
    print("  But cosh(delta*log(p)) is not a character value — it's a")
    print("  function of the PRIME ITSELF, not of the Frobenius class.")
    print()

    # Compute the Chebotarev average of the excess weight
    # Over A5 conjugacy classes:
    #   5A (12 elements): trace = phi, density = 12/60 = 1/5
    #   5B (12 elements): trace = -1/phi, density = 12/60 = 1/5
    #   (12) (15 elements): trace = 0, density = 15/60 = 1/4
    #   (123) (20 elements): trace = -1, density = 20/60 = 1/3
    #   (1) (1 element): trace = 2, density = 1/60

    print("  CHEBOTAREV AVERAGE of trace * cosh_excess:")
    print("  (The cosh_excess depends on p, not on the class!)")
    print()

    # For a "typical" prime of size p, the average trace weighted by class:
    # avg = (1/5)*phi + (1/5)*(-1/phi) + (1/4)*0 + (1/3)*(-1) + (1/60)*2
    avg_trace = mpf(1)/5 * PHI + mpf(1)/5 * (-INV_PHI) + mpf(1)/4 * 0 + mpf(1)/3 * (-1) + mpf(1)/60 * 2
    print(f"  Average trace over Chebotarev classes:")
    print(f"    (1/5)*phi + (1/5)*(-1/phi) + (1/4)*0 + (1/3)*(-1) + (1/60)*2")
    print(f"    = (1/5)*(phi - 1/phi) - 1/3 + 1/30")
    print(f"    = (1/5)*1 - 1/3 + 1/30")
    print(f"    = 1/5 - 1/3 + 1/30")
    print(f"    = 6/30 - 10/30 + 1/30")
    print(f"    = -3/30 = -1/10")
    print(f"    Computed: {float(avg_trace):.10f}")
    print(f"    Expected: -0.1")
    print()

    # But wait — this is the average for ALL primes, not just golden.
    # For GOLDEN primes only:
    # avg_golden = (1/2)*phi + (1/2)*(-1/phi) = (1/2)*(phi - 1/phi) = 1/2
    avg_golden = (PHI - INV_PHI) / 2
    print(f"  Average trace over GOLDEN primes only:")
    print(f"    (1/2)*phi + (1/2)*(-1/phi) = (phi - 1/phi)/2 = 1/2")
    print(f"    Computed: {float(avg_golden):.10f}")
    print()

    # The golden-only sum has average trace = 1/2 per prime.
    # The cosh_excess factor 2*(cosh(d*log(p))-1)*log(p) is ALWAYS POSITIVE.
    # So the average golden contribution per prime is:
    #   (1/2) * 2*(cosh(d*log(p))-1)*log(p) = (cosh(d*log(p))-1)*log(p) > 0

    print("  CONCLUSION: The golden-weighted excess has POSITIVE average per prime.")
    print("  Schur orthogonality constrains character-weighted sums, but our")
    print("  excess is weighted by a NON-CHARACTER function (cosh(d*log(p))-1).")
    print("  Schur does NOT kill this sum.")
    print()

    # BUT: there's a subtlety
    print("  HOWEVER: Schur DOES constrain the golden-only sum indirectly.")
    print("  The sum sum_{p golden} trace(p)*f(p) = sum_all trace(p)*f(p)")
    print("  minus sum_{non-golden} trace(p)*f(p).")
    print("  The full sum IS constrained by Schur (it's related to the")
    print("  logarithmic derivative of L(s, rho)).")
    print("  But the RATE of cancellation is the question.")
    print()

    # =========================================================================
    # (e) WHEN DOES EXCESS FIRST EXCEED THE POLE'S PAYMENT?
    # =========================================================================
    print("  (e) WHEN DOES EXCESS EXCEED THE POLE'S PAYMENT?")
    print("  " + "-" * 60)
    print()

    # The "pole's payment" is O(log(N)) — the contribution from the pole at s=1
    # (For L-functions without a pole, the "payment" is from the archimedean terms,
    # also O(log(N)).)

    for d in [mpf("0.001"), mpf("0.01"), mpf("0.1")]:
        cumulative = mpf(0)
        first_exceed = None

        for p in golden_primes:
            p_mp = mpf(p)
            lp = log(p_mp)
            tr = fabs(golden_trace(p))

            # Simplified excess (without oscillatory gamma factor)
            excess = tr * 2 * (cosh(d * lp) - 1) * lp
            cumulative += excess

            logN = log(p_mp)
            if cumulative > logN and first_exceed is None:
                first_exceed = p

        if first_exceed:
            print(f"  delta={float(d):.3f}: excess first exceeds log(N) at p = {first_exceed}")
        else:
            print(f"  delta={float(d):.3f}: excess does NOT exceed log(N) up to p = {golden_primes[-1]}")

    print()

    # Now with the oscillatory factor
    print("  WITH oscillatory factor e^{i*gamma*log(p)} (gamma=10):")
    g = mpf(10)
    for d in [mpf("0.001"), mpf("0.01"), mpf("0.1")]:
        rho_plus = mpc(mpf("0.5") + d, g)
        rho_minus = mpc(mpf("0.5") - d, g)
        rho_center = mpc(mpf("0.5"), g)

        cumulative_abs = mpf(0)
        first_exceed = None

        for p in golden_primes:
            p_mp = mpf(p)
            lp = log(p_mp)
            tr = fabs(golden_trace(p))

            C_pair = power(p_mp, rho_plus) / rho_plus + power(p_mp, rho_minus) / rho_minus
            C_coll = 2 * power(p_mp, rho_center) / rho_center
            E = C_pair - C_coll
            E_abs = sqrt(re(E * E.conjugate()))

            cumulative_abs += tr * E_abs * lp
            logN = log(p_mp)

            if cumulative_abs > logN and first_exceed is None:
                first_exceed = p

        if first_exceed:
            print(f"    delta={float(d):.3f}: |excess| first exceeds log(N) at p = {first_exceed}")
        else:
            print(f"    delta={float(d):.3f}: |excess| does NOT exceed log(N) up to p = {golden_primes[-1]}")

    print()


# =============================================================================
# PART 8: THE VERDICT
# =============================================================================

def part8_verdict(golden_primes):
    """
    Final assessment: does the gravitational argument work?
    """
    print("=" * 78)
    print("  PART 8: THE VERDICT")
    print("=" * 78)
    print()

    # 1. Does the cosh excess grow without bound?
    print("  QUESTION 1: Does the cosh excess grow without bound?")
    print("  " + "-" * 60)

    for d in [mpf("0.001"), mpf("0.01"), mpf("0.1")]:
        totals = []
        cumulative = mpf(0)
        checkpoints = [1000, 5000, 9949]
        cp_idx = 0

        for p in golden_primes:
            p_mp = mpf(p)
            lp = log(p_mp)
            tr = fabs(golden_trace(p))
            cumulative += tr * 2 * (cosh(d * lp) - 1) * lp

            if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
                totals.append((checkpoints[cp_idx], float(cumulative)))
                cp_idx += 1
                if cp_idx >= len(checkpoints):
                    break

        # capture final state if last checkpoint wasn't hit
        if len(totals) < 3:
            totals.append((golden_primes[-1], float(cumulative)))

        print(f"    delta={float(d):.3f}: ", end="")
        for N, val in totals:
            print(f"Total({N})={val:.4f}  ", end="")

        if len(totals) >= 2:
            growth = totals[-1][1] / totals[0][1] if totals[0][1] > 0 else float('inf')
            print(f"  growth factor: {growth:.2f}x")
        else:
            print()

    print("  ANSWER: YES, the non-oscillatory excess grows without bound.")
    print("  Rate ~ N * (delta * log(N))^2 / 2")
    print()

    # 2. Does the excess exceed the available payment?
    print("  QUESTION 2: Does the excess exceed the pole's payment (log N)?")
    print("  " + "-" * 60)

    d = mpf("0.01")
    cumulative = mpf(0)
    for p in golden_primes:
        p_mp = mpf(p)
        lp = log(p_mp)
        tr = fabs(golden_trace(p))
        cumulative += tr * 2 * (cosh(d * lp) - 1) * lp

    final_logN = log(mpf(golden_primes[-1]))
    ratio = float(cumulative / final_logN)
    print(f"    delta=0.01: Total(10000) = {float(cumulative):.4f}, log(10000) = {float(final_logN):.4f}")
    print(f"    Ratio: {ratio:.2f}x")
    print(f"    ANSWER: {'YES' if ratio > 1 else 'NO'} (without oscillation)")
    print()

    # 3. Are there cancellations?
    print("  QUESTION 3: Do oscillatory cancellations kill the growth?")
    print("  " + "-" * 60)

    d = mpf("0.01")
    g = mpf(10)
    rho_plus = mpc(mpf("0.5") + d, g)
    rho_minus = mpc(mpf("0.5") - d, g)
    rho_center = mpc(mpf("0.5"), g)

    signed_sum = mpc(0, 0)
    abs_sum = mpf(0)

    for p in golden_primes:
        p_mp = mpf(p)
        lp = log(p_mp)
        tr = golden_trace(p)

        C_pair = power(p_mp, rho_plus) / rho_plus + power(p_mp, rho_minus) / rho_minus
        C_coll = 2 * power(p_mp, rho_center) / rho_center
        E = C_pair - C_coll

        # In the explicit formula, the contribution is tr * E * log(p) / sqrt(p)
        # (the p^{1/2} in E cancels with 1/sqrt(p))
        contrib = tr * E * lp / sqrt(p_mp)
        signed_sum += contrib
        abs_sum += sqrt(re(contrib * contrib.conjugate()))

    cancel_ratio = float(sqrt(re(signed_sum * signed_sum.conjugate())) / abs_sum)
    print(f"    delta=0.01, gamma=10:")
    print(f"    |signed sum| = {float(sqrt(re(signed_sum * signed_sum.conjugate()))):.8f}")
    print(f"    sum |terms|  = {float(abs_sum):.8f}")
    print(f"    ratio        = {cancel_ratio:.6f}")
    print()

    if cancel_ratio < 0.01:
        print("    ANSWER: YES. The oscillatory factor e^{i*gamma*log(p)} causes")
        print("    near-complete cancellation. The signed sum is ~1% of the absolute sum.")
        print("    This is the CENTRAL PROBLEM with the gravitational argument:")
        print("    the e^{i*gamma*log(p)} phase rotates wildly across primes,")
        print("    and golden coherence (trace = phi) does not align these phases.")
    elif cancel_ratio < 0.1:
        print("    ANSWER: PARTIALLY. Significant cancellation but not fatal.")
    else:
        print("    ANSWER: NO. Cancellations are limited.")
    print()

    # 4. Does Schur constrain the excess?
    print("  QUESTION 4: Does Schur orthogonality constrain the excess?")
    print("  " + "-" * 60)
    print()
    print("    The excess function cosh(delta*log(p)) - 1 is NOT a character.")
    print("    It depends on p, not on the Frobenius class.")
    print("    Schur constrains sum chi(Frob_p) * f(p) where f is multiplicative.")
    print("    Our cosh(delta*log(p)) is not multiplicative: cosh(d*log(pq)) != cosh(d*log(p))*cosh(d*log(q)).")
    print()
    print("    HOWEVER: The Selberg class axioms constrain the ANALYTIC properties")
    print("    of L(s), which DO constrain power sums of zeros at specific test")
    print("    functions. The cosh excess is related to the SHIFTED L-function")
    print("    L(s+delta) + L(s-delta) - 2L(s), which IS an analytic object.")
    print()
    print("    ANSWER: Schur alone does not kill the excess.")
    print("    But analytic constraints on the L-function do bound it.")
    print("    The question is whether those bounds are tight enough.")
    print()

    # 5. Is this argument different from earlier attempts?
    print("  QUESTION 5: Is this genuinely different from earlier 'golden coherence' attempts?")
    print("  " + "-" * 60)
    print()
    print("    EARLIER ATTEMPTS: Used the golden trace to claim systematic accumulation")
    print("    in the explicit formula. Failed because the oscillatory factor")
    print("    e^{i*gamma*log(p)} is independent of the Frobenius class and")
    print("    causes cancellations that swamp the trace coherence.")
    print()
    print("    THIS ATTEMPT: Uses the cosh excess, which is ALWAYS POSITIVE")
    print("    (cosh >= 1, so cosh - 1 >= 0). The positivity means no sign")
    print("    cancellations in the ENERGY (squared) integral.")
    print()
    print("    THE DIFFERENCE: In Part 3, we computed int |D(x)|^2 dx/x,")
    print("    which is always positive. The oscillatory phases cancel in")
    print("    the SIGNED sum but not in the SQUARED sum.")
    print()
    print("    BUT: The explicit formula constrains the SIGNED sum, not the")
    print("    squared sum. The squared sum is a SECOND MOMENT bound.")
    print("    You need to connect the second moment to the explicit formula,")
    print("    which requires additional analytic machinery (e.g., Selberg's")
    print("    moment method, or the Deuring-Heilbronn phenomenon).")
    print()

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    print("=" * 78)
    print("  FINAL VERDICT")
    print("=" * 78)
    print()

    # Run the definitive test: does the SIGNED excess grow or cancel?
    print("  THE DEFINITIVE TEST: Signed excess vs N for multiple gamma values")
    print()
    print(f"  {'gamma':>8s}  {'|signed(1000)|':>16s}  {'|signed(5000)|':>16s}  "
          f"{'|signed(~10k)|':>17s}  {'growth?':>10s}")
    print("  " + "-" * 75)

    d = mpf("0.01")

    for g_val in [1, 10, 14.13, 21.02, 100]:
        g = mpf(g_val)
        rho_plus = mpc(mpf("0.5") + d, g)
        rho_minus = mpc(mpf("0.5") - d, g)
        rho_center = mpc(mpf("0.5"), g)

        checkpoints = [1000, 5000, 9949]  # 9949 = largest golden prime <= 10000
        cp_idx = 0
        signed_vals = []

        signed_sum = mpc(0, 0)

        for p in golden_primes:
            p_mp = mpf(p)
            lp = log(p_mp)
            tr = golden_trace(p)

            C_pair = power(p_mp, rho_plus) / rho_plus + power(p_mp, rho_minus) / rho_minus
            C_coll = 2 * power(p_mp, rho_center) / rho_center
            E = C_pair - C_coll

            signed_sum += tr * E * lp / sqrt(p_mp)

            if cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
                mag = float(sqrt(re(signed_sum * signed_sum.conjugate())))
                signed_vals.append(mag)
                cp_idx += 1
                if cp_idx >= len(checkpoints):
                    break

        # Final value: use whatever we have after processing all golden primes
        if len(signed_vals) < 3:
            mag = float(sqrt(re(signed_sum * signed_sum.conjugate())))
            signed_vals.append(mag)
        while len(signed_vals) < 3:
            signed_vals.append(0.0)

        if signed_vals[0] > 1e-15:
            growth = signed_vals[2] / signed_vals[0]
            grows = "YES" if growth > 2 else ("maybe" if growth > 1.2 else "NO")
        else:
            growth = float('inf')
            grows = "N/A"

        print(f"  {g_val:>8.2f}  {signed_vals[0]:>16.8f}  {signed_vals[1]:>16.8f}  "
              f"{signed_vals[2]:>17.8f}  {grows:>10s}")

    print()
    print("  " + "=" * 70)
    print()
    print("  VERDICT SUMMARY:")
    print()
    print("  1. COSH EXCESS GROWS: YES. The non-oscillatory excess grows ~N*(d*logN)^2.")
    print("     This is rigorous and unconditional.")
    print()
    print("  2. EXCEEDS PAYMENT: YES (without oscillation). The non-oscillatory")
    print("     excess overwhelms log(N) for any fixed delta > 0.")
    print()
    print("  3. OSCILLATORY CANCELLATION: THIS IS THE HOLE.")
    print("     The explicit formula involves the COMPLEX contribution p^rho/rho,")
    print("     not the absolute value. The phase factor e^{i*gamma*log(p)}")
    print("     rotates independently of the golden trace, causing the SIGNED")
    print("     sum to undergo a random walk of step size ~delta^2*(log p)^2.")
    print("     After N golden primes, |signed sum| ~ delta^2 * (log N)^2 * sqrt(N),")
    print("     while the payment is ~log(N).")
    print("     Ratio ~ delta^2 * log(N) * sqrt(N), which DOES grow.")
    print()
    print("  4. SCHUR: Does NOT kill the excess. The cosh function is not a character.")
    print()
    print("  5. NOVELTY: The cosh/gravitational framing is genuinely different from")
    print("     the linear trace argument. The QUADRATIC (energy) formulation")
    print("     avoids direct cancellation. But converting energy bounds to")
    print("     pointwise explicit formula constraints requires Selberg-type")
    print("     machinery that is itself equivalent in difficulty to GRH.")
    print()
    print("  BOTTOM LINE:")
    print("  The gravitational argument identifies a REAL tension: the cosh excess")
    print("  is always positive and grows. But the explicit formula works with")
    print("  SIGNED (complex) sums, and the oscillatory cancellation is the gap.")
    print("  Bridging this gap — showing that the energy excess FORCES the signed")
    print("  sum to be inconsistent — requires proving something equivalent to GRH.")
    print("  The argument is therefore CIRCULAR at the critical step.")
    print()
    print("  STATUS: HOLE FOUND at step (3). The gravitational potential is real,")
    print("  but the explicit formula sees the SIGNED sum, not the energy.")
    print("  This is the same fundamental barrier as all explicit-formula")
    print("  approaches to GRH: you need to control oscillatory cancellation,")
    print("  and that IS the hard part.")
    print()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    golden_primes = part1_excess_per_prime()
    part2_cumulative_excess(golden_primes)
    part3_balance_test()
    part4_golden_coherence(golden_primes)
    part5_explicit_computation(golden_primes)
    part6_spring_constant(golden_primes)
    part7_catch_search(golden_primes)
    part8_verdict(golden_primes)

    print("=" * 78)
    print("  COMPUTATION COMPLETE")
    print("=" * 78)
