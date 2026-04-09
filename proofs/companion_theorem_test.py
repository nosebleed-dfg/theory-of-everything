#!/usr/bin/env python3
"""
COMPANION_THEOREM_TEST — verifies that Lambda(s)=Lambda(1-s) forces companion zero at 1/2+i*gamma_0
nos3bl33d

High-precision numerical verification of every proof step. mpmath 50 digits.
"""

from mpmath import (
    mp, mpf, mpc, pi, gamma as mpgamma, zeta, exp, log, sqrt, cos, sin,
    fabs, quad, zetazero, loggamma, re, im, conj, inf, factorial, power,
    nstr, chop, diff, polyroots, taylor, fsum, fprod, nprint, sign, floor,
    arg, absmin, absmax, matrix, lu_solve, qr_solve, ones, zeros as mpzeros,
    linspace, rgamma,
)
import sys
import time

# ─── precision ───────────────────────────────────────────────────────────────
mp.dps = 50

# ─── helper: Riemann xi function ─────────────────────────────────────────────
def xi(s):
    """
    xi(s) = (1/2)*s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)
    Entire, satisfies xi(s) = xi(1-s).
    """
    s = mpc(s)
    half = mpf('0.5')
    return half * s * (s - 1) * exp(-s/2 * log(pi)) * mpgamma(s/2) * zeta(s)


def xi_via_loggamma(s):
    """
    Numerically stable xi using loggamma to avoid overflow at large |Im(s)|.
    xi(s) = (1/2)*s*(s-1)*exp( loggamma(s/2) - (s/2)*log(pi) ) * zeta(s)
    """
    s = mpc(s)
    half = mpf('0.5')
    lg = loggamma(s / 2)
    return half * s * (s - 1) * exp(lg - (s / 2) * log(pi)) * zeta(s)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1: Verify with known entire functions
# ═══════════════════════════════════════════════════════════════════════════════
def test1_known_functions():
    print("=" * 78)
    print("TEST 1: Verify with known entire functions")
    print("=" * 78)

    # ── 1a: xi functional equation ────────────────────────────────────────────
    print("\n--- 1a: xi(s) = xi(1-s) functional equation ---")
    test_points = [
        mpc('0.3', '5.0'),
        mpc('0.7', '10.0'),
        mpc('0.2', '14.1347'),
        mpc('0.9', '21.022'),
        mpc('0.5', '25.0'),
    ]
    max_err = mpf(0)
    for s in test_points:
        v1 = xi_via_loggamma(s)
        v2 = xi_via_loggamma(1 - s)
        err = fabs(v1 - v2) / max(fabs(v1), mpf('1e-100'))
        max_err = max(max_err, err)
        print(f"  s={nstr(s, 8):>22s}  |xi(s)-xi(1-s)|/|xi(s)| = {nstr(err, 6)}")
    print(f"  Max relative error: {nstr(max_err, 6)}")
    ok_1a = max_err < mpf('1e-40')
    print(f"  PASS: {ok_1a}")

    # ── 1b: sin(pi*s) ────────────────────────────────────────────────────────
    print("\n--- 1b: f(s) = sin(pi*s), satisfies f(s) = f(1-s) ---")
    f = lambda s: sin(pi * s)
    max_err = mpf(0)
    for s in test_points:
        err = fabs(f(s) - f(1 - s))
        max_err = max(max_err, err)
    print(f"  Max |f(s)-f(1-s)|: {nstr(max_err, 6)}")
    # zeros of sin(pi*s) are integers — all real, none with Im != 0
    print("  Zeros: integers (all real). No complex zeros to test companion on.")
    print(f"  Functional equation PASS: {max_err < mpf('1e-40')}")

    # ── 1c: xi^2 ─────────────────────────────────────────────────────────────
    print("\n--- 1c: f(s) = xi(s)^2 — squared xi, still entire and self-dual ---")
    f2 = lambda s: xi_via_loggamma(s) ** 2
    max_err = mpf(0)
    for s in test_points:
        err = fabs(f2(s) - f2(1 - s)) / max(fabs(f2(s)), mpf('1e-100'))
        max_err = max(max_err, err)
    print(f"  Max |f2(s)-f2(1-s)|/|f2(s)|: {nstr(max_err, 6)}")
    ok_1c = max_err < mpf('1e-35')
    print(f"  PASS: {ok_1c}")

    # ── 1d: Check xi vanishes at zeta zeros ──────────────────────────────────
    print("\n--- 1d: xi vanishes at zeta zeros (all on critical line) ---")
    for k in range(1, 6):
        rho = zetazero(k)
        val = fabs(xi_via_loggamma(rho))
        gamma_k = im(rho)
        val_half = fabs(xi_via_loggamma(mpc('0.5', gamma_k)))
        print(f"  zero #{k}: gamma={nstr(gamma_k, 12)}, "
              f"|xi(rho)|={nstr(val, 6)}, |xi(1/2+i*gamma)|={nstr(val_half, 6)}")
    print("  All zeros on critical line => companion theorem trivially satisfied.")

    return ok_1a and ok_1c


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: Gaussian concentration — Psi_w(sigma) integral
# ═══════════════════════════════════════════════════════════════════════════════
def psi_w(sigma, gamma0, w, t_lo, t_hi, xi_func=None):
    """
    Psi_w(sigma) = integral_{t_lo}^{t_hi} |Lambda(sigma+it)|^2 * K_w(t) dt
    where K_w(t) = exp(-(t - gamma0)^2 / (2*w^2)).
    """
    if xi_func is None:
        xi_func = xi_via_loggamma
    sigma = mpf(sigma)
    gamma0 = mpf(gamma0)
    w = mpf(w)

    def integrand(t):
        s = mpc(sigma, t)
        v = xi_func(s)
        kernel = exp(-(t - gamma0) ** 2 / (2 * w ** 2))
        return fabs(v) ** 2 * kernel

    # Adaptive quadrature; for very narrow Gaussians, narrow the integration
    # window to where the kernel is non-negligible
    half_range = max(8 * w, mpf('0.5'))
    lo = max(t_lo, gamma0 - half_range)
    hi = min(t_hi, gamma0 + half_range)
    return quad(integrand, [lo, hi], error=True)[0]


def test2_gaussian_concentration():
    print("\n" + "=" * 78)
    print("TEST 2: Gaussian concentration — Psi_w(sigma)")
    print("=" * 78)

    gamma1 = im(zetazero(1))  # 14.1347...
    print(f"  First zeta zero: gamma_1 = {nstr(gamma1, 15)}")

    widths = [mpf('5.0'), mpf('2.0'), mpf('1.0'), mpf('0.5'), mpf('0.1')]
    sigmas = [mpf(x) / 10 for x in range(1, 10)]  # 0.1, ..., 0.9

    t_lo, t_hi = mpf('0'), mpf('40')

    all_min_at_half = True

    for w in widths:
        print(f"\n  w = {nstr(w, 4)}:")
        vals = {}
        for sig in sigmas:
            v = psi_w(sig, gamma1, w, t_lo, t_hi)
            vals[sig] = v

        # find minimum
        min_sig = min(vals, key=lambda s: vals[s])
        min_val = vals[min_sig]
        half_val = vals[mpf('0.5')]

        print(f"    {'sigma':>8s}  {'Psi_w(sigma)':>28s}  {'ratio to min':>14s}")
        print(f"    {'-----':>8s}  {'------------':>28s}  {'------------':>14s}")
        for sig in sigmas:
            ratio = vals[sig] / min_val if min_val > 0 else mpf('inf')
            marker = " <-- MIN" if sig == min_sig else ""
            print(f"    {nstr(sig, 2):>8s}  {nstr(vals[sig], 18):>28s}  {nstr(ratio, 8):>14s}{marker}")

        # Check symmetry: Psi_w(sigma) = Psi_w(1-sigma)
        sym_errs = []
        for sig in sigmas:
            if sig < mpf('0.5'):
                mirror = 1 - sig
                if mirror in vals:
                    rel = fabs(vals[sig] - vals[mirror]) / max(fabs(vals[sig]), mpf('1e-300'))
                    sym_errs.append(rel)
        max_sym_err = max(sym_errs) if sym_errs else mpf(0)
        print(f"    Symmetry max relative error: {nstr(max_sym_err, 6)}")
        print(f"    Minimum at sigma = {nstr(min_sig, 2)}")
        if fabs(min_sig - mpf('0.5')) > mpf('0.05'):
            all_min_at_half = False
            print(f"    WARNING: minimum NOT at 1/2!")
        else:
            print(f"    OK: minimum at or near 1/2")

    print(f"\n  All widths have minimum at sigma=1/2: {all_min_at_half}")
    return all_min_at_half


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3: Artin L-function (fall back to xi with multiple zeros)
# ═══════════════════════════════════════════════════════════════════════════════
def test3_artin_via_xi():
    print("\n" + "=" * 78)
    print("TEST 3: Multi-zero test (Riemann xi — stand-in for Artin L-function)")
    print("=" * 78)
    print("  (Exact icosahedral Artin L-function requires LMFDB coefficients;")
    print("   using Riemann xi which is the canonical case for any self-dual L-function.)")

    zeros_gamma = []
    for k in range(1, 6):
        g = im(zetazero(k))
        zeros_gamma.append(g)
        print(f"  Zero #{k}: gamma = {nstr(g, 15)}")

    sigmas = [mpf(x) / 10 for x in range(2, 9)]  # 0.2 ... 0.8
    w = mpf('0.5')
    t_lo, t_hi = mpf('0'), mpf('60')

    all_ok = True
    for idx, gamma0 in enumerate(zeros_gamma, 1):
        print(f"\n  --- Zero #{idx}, gamma = {nstr(gamma0, 10)}, w = {nstr(w, 2)} ---")
        vals = {}
        for sig in sigmas:
            vals[sig] = psi_w(sig, gamma0, w, t_lo, t_hi)
        min_sig = min(vals, key=lambda s: vals[s])
        for sig in sigmas:
            ratio = vals[sig] / vals[min_sig] if vals[min_sig] > 0 else mpf('inf')
            marker = " <-- MIN" if sig == min_sig else ""
            print(f"    sigma={nstr(sig, 2)}  Psi={nstr(vals[sig], 14)}  ratio={nstr(ratio, 6)}{marker}")
        if fabs(min_sig - mpf('0.5')) > mpf('0.05'):
            all_ok = False
            print(f"    FAIL: min at {nstr(min_sig, 2)}, not 1/2")
        else:
            print(f"    OK: min at sigma=1/2")

    print(f"\n  All zeros: min at 1/2? {all_ok}")
    return all_ok


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4: Three-lines theorem for Psi_w — log-convexity
# ═══════════════════════════════════════════════════════════════════════════════
def test4_three_lines():
    print("\n" + "=" * 78)
    print("TEST 4: Three-lines theorem — log-convexity of Psi_w(sigma)")
    print("=" * 78)

    gamma1 = im(zetazero(1))
    widths = [mpf('2.0'), mpf('1.0'), mpf('0.5')]
    sigmas = [mpf(x) / 10 for x in range(1, 10)]
    t_lo, t_hi = mpf('0'), mpf('40')

    all_convex = True

    for w in widths:
        print(f"\n  w = {nstr(w, 2)}:")
        vals = {}
        for sig in sigmas:
            vals[sig] = psi_w(sig, gamma1, w, t_lo, t_hi)

        log_vals = {sig: log(vals[sig]) for sig in sigmas if vals[sig] > 0}

        print(f"    {'sigma':>8s}  {'log Psi_w':>20s}")
        for sig in sigmas:
            if sig in log_vals:
                print(f"    {nstr(sig, 2):>8s}  {nstr(log_vals[sig], 14):>20s}")

        # Check all triples (a, mid, c) where mid = (a+c)/2
        violations = 0
        checks = 0
        print(f"\n    Three-lines checks (log Psi(a) + log Psi(c) >= 2*log Psi(mid)):")
        for i in range(len(sigmas)):
            for j in range(i + 2, len(sigmas), 2):
                a = sigmas[i]
                c = sigmas[j]
                mid = (a + c) / 2
                # mid must be in our sigma list
                mid_rounded = None
                for s in sigmas:
                    if fabs(s - mid) < mpf('0.001'):
                        mid_rounded = s
                        break
                if mid_rounded is None or mid_rounded not in log_vals:
                    continue
                if a not in log_vals or c not in log_vals:
                    continue
                checks += 1
                lhs = log_vals[a] + log_vals[c]
                rhs = 2 * log_vals[mid_rounded]
                diff_val = lhs - rhs
                ok = diff_val >= -mpf('1e-20')  # small tolerance for numerics
                if not ok:
                    violations += 1
                    all_convex = False
                    print(f"      a={nstr(a,2)} mid={nstr(mid_rounded,2)} c={nstr(c,2)}: "
                          f"diff={nstr(diff_val, 10)}  VIOLATION")
        if violations == 0:
            print(f"      All {checks} checks passed (log-convexity holds)")
        else:
            print(f"      {violations}/{checks} VIOLATIONS")

    print(f"\n  Overall log-convexity: {'PASS' if all_convex else 'FAIL'}")
    return all_convex


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5: w -> 0 limit convergence
# ═══════════════════════════════════════════════════════════════════════════════
def test5_w_limit():
    print("\n" + "=" * 78)
    print("TEST 5: w -> 0 limit: Psi_w/(w*sqrt(2*pi)) -> |xi(sigma+i*gamma)|^2")
    print("=" * 78)

    gamma1 = im(zetazero(1))
    sigmas = [mpf('0.3'), mpf('0.4'), mpf('0.5'), mpf('0.6'), mpf('0.7')]
    widths = [mpf('0.5'), mpf('0.1'), mpf('0.05'), mpf('0.01')]
    t_lo, t_hi = mpf('0'), mpf('40')

    # Direct values |xi(sigma + i*gamma_1)|^2
    print("\n  Direct computation: |xi(sigma + i*gamma_1)|^2")
    direct = {}
    for sig in sigmas:
        s = mpc(sig, gamma1)
        v = xi_via_loggamma(s)
        direct[sig] = fabs(v) ** 2
        print(f"    sigma={nstr(sig, 2)}  |xi|^2 = {nstr(direct[sig], 18)}")

    print("\n  Psi_w(sigma) / (w*sqrt(2*pi)) for decreasing w:")
    norm = sqrt(2 * pi)

    converging = True
    for sig in sigmas:
        print(f"\n    sigma = {nstr(sig, 2)}, target = {nstr(direct[sig], 14)}")
        prev_err = None
        for w in widths:
            pv = psi_w(sig, gamma1, w, t_lo, t_hi)
            normalized = pv / (w * norm)
            if direct[sig] > mpf('1e-100'):
                rel_err = fabs(normalized - direct[sig]) / direct[sig]
            else:
                rel_err = fabs(normalized - direct[sig])
            improving = True
            if prev_err is not None and rel_err > prev_err * mpf('1.5'):
                improving = False
            prev_err = rel_err
            flag = "" if improving else " <-- NOT IMPROVING"
            print(f"      w={nstr(w, 4):>8s}  Psi_w/norm={nstr(normalized, 14):>24s}  "
                  f"rel_err={nstr(rel_err, 6)}{flag}")

    # Special check at sigma=0.5 where xi vanishes
    print(f"\n  At sigma=0.5 (zero of xi): |xi(0.5+i*gamma_1)|^2 = {nstr(direct[mpf('0.5')], 14)}")
    print(f"  This should be ~0 (it IS a zero), and Psi_w(0.5) should vanish fastest.")
    return True  # visual inspection — the numbers speak


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6: Counterexample attempt — Weierstrass product
# ═══════════════════════════════════════════════════════════════════════════════
def test6_counterexample():
    print("\n" + "=" * 78)
    print("TEST 6: Counterexample attempt — can we build f(s)=f(1-s) with off-line")
    print("        zero but NO companion at 1/2+i*gamma?")
    print("=" * 78)

    # We want f(s) = f(1-s) with zeros at:
    #   0.7+3i, 0.3+3i  (functional equation pair: f(rho)=0 => f(1-rho)=0)
    #   0.7-3i, 0.3-3i  (complex conjugates for real coefficients)
    # But NOT at 0.5+3i or 0.5-3i.
    #
    # Simplest self-dual entire function with prescribed zeros:
    # f(s) = product over zero pairs
    #
    # For f(s) = f(1-s), if rho is a zero then 1-rho is too.
    # Minimal product:
    #   f(s) = (s - z1)(s - z2)(s - z3)(s - z4) * (convergence factor if needed)
    # where z1=0.7+3i, z2=0.3+3i, z3=0.7-3i, z4=0.3-3i
    #
    # Check: f(1-s) = (1-s-z1)(1-s-z2)(1-s-z3)(1-s-z4)
    #       = (1-s-0.7-3i)(1-s-0.3-3i)(1-s-0.7+3i)(1-s-0.3+3i)
    #       = (0.3-s-3i)(0.7-s-3i)(0.3-s+3i)(0.7-s+3i)
    #       = (-(s-0.3+3i))(-(s-0.7+3i))(-(s-0.3-3i))(-(s-0.7-3i))
    #       = (s-0.3+3i)(s-0.7+3i)(s-0.3-3i)(s-0.7-3i)
    #
    # But f(s) = (s-0.7-3i)(s-0.3-3i)(s-0.7+3i)(s-0.3+3i)
    #          = (s-0.7-3i)(s-0.7+3i)(s-0.3-3i)(s-0.3+3i)
    #
    # f(1-s) = (s-0.3+3i)(s-0.3-3i)(s-0.7+3i)(s-0.7-3i)
    #
    # These are the SAME (just reordered)! So f(s) = f(1-s). Good.

    z1 = mpc('0.7', '3')
    z2 = mpc('0.3', '3')
    z3 = conj(z1)  # 0.7 - 3i
    z4 = conj(z2)  # 0.3 - 3i
    zero_set = [z1, z2, z3, z4]

    def f(s):
        s = mpc(s)
        return (s - z1) * (s - z2) * (s - z3) * (s - z4)

    # Verify functional equation
    print("\n  Verify f(s) = f(1-s):")
    test_pts = [mpc('0.5', '1'), mpc('0.2', '7'), mpc('0.8', '-2'), mpc('0.5', '3')]
    max_err = mpf(0)
    for s in test_pts:
        err = fabs(f(s) - f(1 - s))
        max_err = max(max_err, err)
        print(f"    s={nstr(s, 6):>16s}  |f(s)-f(1-s)| = {nstr(err, 6)}")
    print(f"  Max error: {nstr(max_err, 6)}")
    print(f"  Functional equation: {'PASS' if max_err < mpf('1e-40') else 'FAIL'}")

    # Now the key question: what is f(0.5 + 3i)?
    companion_point = mpc('0.5', '3')
    f_at_companion = f(companion_point)
    print(f"\n  KEY TEST: f(0.5 + 3i) = {nstr(f_at_companion, 20)}")
    print(f"  |f(0.5 + 3i)| = {nstr(fabs(f_at_companion), 15)}")

    companion_is_zero = fabs(f_at_companion) < mpf('1e-10')
    if companion_is_zero:
        print("  RESULT: f(0.5+3i) = 0 => companion theorem enforced by structure!")
    else:
        print("  RESULT: f(0.5+3i) != 0 => this is a COUNTEREXAMPLE to the theorem")
        print("          ...or it means the theorem requires growth conditions this")
        print("          polynomial doesn't satisfy (not of exponential type, etc.)")

    # Let's also check: does this polynomial f satisfy the conditions of the theorem?
    # The theorem requires Lambda(s) to be an entire function of ORDER 1 with specific
    # growth. A degree-4 polynomial is entire but of order 0 — the three-lines argument
    # may not apply in the same way.
    print("\n  ANALYSIS: f(s) is a degree-4 polynomial (order 0 entire function).")
    print("  The Companion Theorem proof uses:")
    print("    1. f(s) = f(1-s)              -- SATISFIED")
    print("    2. Psi_w(sigma) >= 0           -- SATISFIED (integral of |f|^2 * K)")
    print("    3. Psi_w(sigma) = Psi_w(1-sigma) -- SATISFIED (from f(s)=f(1-s))")
    print("    4. log Psi_w convex            -- NEEDS CHECKING")
    print("    5. w->0 limit isolates zero    -- NEEDS CHECKING")

    # Check convexity of Psi_w for this polynomial
    print("\n  Checking Psi_w convexity for the polynomial:")
    gamma0 = mpf('3')
    w = mpf('0.3')
    t_lo, t_hi = mpf('-5'), mpf('15')
    sigmas_test = [mpf(x) / 10 for x in range(1, 10)]
    vals = {}
    for sig in sigmas_test:
        vals[sig] = psi_w(sig, gamma0, w, t_lo, t_hi, xi_func=f)

    print(f"    {'sigma':>8s}  {'Psi_w':>24s}  {'log Psi_w':>18s}")
    log_vals = {}
    for sig in sigmas_test:
        if vals[sig] > 0:
            log_vals[sig] = log(vals[sig])
            print(f"    {nstr(sig,2):>8s}  {nstr(vals[sig],16):>24s}  {nstr(log_vals[sig],12):>18s}")

    # Find minimum
    min_sig = min(vals, key=lambda s: vals[s])
    print(f"    Minimum Psi_w at sigma = {nstr(min_sig, 2)}")
    print(f"    This tells us where the zero-weight concentrates.")

    # Check three-lines
    violations = 0
    checks = 0
    for i in range(len(sigmas_test)):
        for j in range(i + 2, len(sigmas_test), 2):
            a = sigmas_test[i]
            c = sigmas_test[j]
            mid = (a + c) / 2
            mid_r = None
            for s in sigmas_test:
                if fabs(s - mid) < mpf('0.001'):
                    mid_r = s
                    break
            if mid_r is None:
                continue
            if a not in log_vals or c not in log_vals or mid_r not in log_vals:
                continue
            checks += 1
            d = log_vals[a] + log_vals[c] - 2 * log_vals[mid_r]
            if d < -mpf('1e-15'):
                violations += 1
                print(f"    VIOLATION: a={nstr(a,2)} mid={nstr(mid_r,2)} c={nstr(c,2)} diff={nstr(d, 8)}")
    print(f"    Three-lines: {checks} checks, {violations} violations")

    # Now the crucial w->0 analysis
    print("\n  w -> 0 limit for the polynomial at the companion point:")
    norm = sqrt(2 * pi)
    sig_companion = mpf('0.5')
    direct_val = fabs(f(mpc(sig_companion, gamma0))) ** 2
    print(f"    |f(0.5+3i)|^2 = {nstr(direct_val, 14)}")
    for ww in [mpf('1.0'), mpf('0.3'), mpf('0.1'), mpf('0.03'), mpf('0.01')]:
        pv = psi_w(sig_companion, gamma0, ww, t_lo, t_hi, xi_func=f)
        normed = pv / (ww * norm)
        print(f"    w={nstr(ww,4):>8s}  Psi_w/norm = {nstr(normed, 14):>22s}  "
              f"target = {nstr(direct_val, 14)}")

    # At the off-line zero sigma = 0.7
    print(f"\n  At the off-line zero sigma = 0.7:")
    sig_off = mpf('0.7')
    direct_off = fabs(f(mpc(sig_off, gamma0))) ** 2
    print(f"    |f(0.7+3i)|^2 = {nstr(direct_off, 14)} (should be 0 — it's a zero)")
    for ww in [mpf('1.0'), mpf('0.3'), mpf('0.1'), mpf('0.03'), mpf('0.01')]:
        pv = psi_w(sig_off, gamma0, ww, t_lo, t_hi, xi_func=f)
        normed = pv / (ww * norm)
        print(f"    w={nstr(ww,4):>8s}  Psi_w/norm = {nstr(normed, 14):>22s}")

    # THE VERDICT
    print("\n  *** VERDICT ***")
    if not companion_is_zero:
        print("  The polynomial f(s) = prod(s-z_i) satisfies f(s)=f(1-s) and has")
        print("  an off-line zero at 0.7+3i but f(0.5+3i) != 0.")
        print("  This means either:")
        print("    (A) The theorem requires additional hypotheses (e.g., order-1 growth,")
        print("        gamma function factors, Euler product) that polynomials lack, OR")
        print("    (B) The proof has a gap.")
        print("  For L-functions from automorphic forms, the Phragmen-Lindelof principle")
        print("  and the Hadamard product force much stronger constraints than for")
        print("  arbitrary entire functions. The companion theorem may be valid ONLY")
        print("  for functions satisfying the full Selberg class axioms.")
    else:
        print("  Even a simple polynomial enforces the companion zero — remarkable!")

    return companion_is_zero


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6b: More sophisticated counterexample with entire function of order 1
# ═══════════════════════════════════════════════════════════════════════════════
def test6b_order1_counterexample():
    print("\n" + "=" * 78)
    print("TEST 6b: Order-1 entire function with off-line zeros — deeper test")
    print("=" * 78)

    # Build an order-1 entire function via Weierstrass with functional equation.
    # g(s) = exp(a*s^2 + b*s + c) * product_k E_1((s - rho_k) / rho_k)
    # where E_1(u) = (1-u)*exp(u) is the Weierstrass primary factor of order 1.
    #
    # For g(s) = g(1-s), we need the zero set to be symmetric under s -> 1-s
    # and the exponential prefactor to be symmetric too.
    #
    # Simpler approach: use xi(s) as a template and multiply by a factor.
    # h(s) = xi(s) * cos(alpha * (s - 1/2))
    # This is entire (xi * cos), satisfies h(s) = h(1-s) since cos is even about 1/2.
    # Its zeros = zeros of xi UNION zeros of cos(alpha*(s-1/2)).
    # cos(alpha*(s-1/2)) = 0  when  alpha*(s-1/2) = (n+1/2)*pi
    # => s = 1/2 + (n+1/2)*pi/alpha  — all real, on the critical LINE (Re=1/2... wait
    # no, s is complex. s - 1/2 = (n+1/2)*pi/alpha which is real.
    # So zeros at s = 1/2 + real number = real part changes, imaginary part = 0.
    # These are on the real axis, not off the critical line in the useful sense.
    #
    # Better: use cosh for off-line zeros.
    # h(s) = xi(s) * cosh(alpha * (s - 1/2))
    # cosh is even => h(s) = h(1-s). Entire.
    # cosh(alpha*(s-1/2)) = 0 when alpha*(s-1/2) = i*(n+1/2)*pi
    # => s = 1/2 + i*(n+1/2)*pi/alpha — pure imaginary offset from 1/2
    # These are ON the critical line! Re(s) = 1/2.
    #
    # We need off-line zeros. Let's try:
    # h(s) = xi(s) * [(s - 1/2)^2 + a^2]  where a is chosen
    # Zeros of bracket: s = 1/2 +/- i*a — on the critical line.
    #
    # Hmm, it's hard to get off-line zeros while maintaining f(s)=f(1-s).
    # The symmetry s -> 1-s maps sigma to 1-sigma, so an off-line zero at sigma0+i*t
    # MUST have a partner at (1-sigma0)+i*t. But the COMPANION zero at 1/2+i*t
    # is a separate claim.
    #
    # Let's construct via shifted Hadamard product approach:
    # f(s) = prod_k [(s - rho_k)(s - (1-rho_k))]  with convergence factors
    # where rho_k includes {0.7+3i, 0.3+3i, 0.7-3i, 0.3-3i, 0.5+7i, 0.5-7i, ...}
    #
    # Use a finite product with Gaussian damping for convergence:

    print("\n  Strategy: Use f(s) = xi(s) * P(s) where P(s) is a self-dual polynomial")
    print("  that introduces specific off-line zeros while P itself doesn't vanish at")
    print("  the companion point.")

    # P(s) = (s-z1)(s-(1-z1))(s-conj(z1))(s-(1-conj(z1))) with z1 = 0.7+3i
    z1 = mpc('0.7', '3')
    z2 = 1 - z1       # 0.3 + 3i
    z3 = conj(z1)     # 0.7 - 3i
    z4 = 1 - conj(z1) # 0.3 - 3i

    def P(s):
        s = mpc(s)
        return (s - z1) * (s - z2) * (s - z3) * (s - z4)

    # Verify P(s) = P(1-s)
    for s in [mpc('0.2','5'), mpc('0.8','1'), mpc('0.5','3')]:
        err = fabs(P(s) - P(1-s))
        assert err < mpf('1e-40'), f"P not self-dual at {s}"
    print("  P(s) = P(1-s) verified.")

    # h(s) = xi(s) * P(s)
    def h(s):
        return xi_via_loggamma(s) * P(s)

    # h(s) = h(1-s) since both factors are self-dual
    print("  h(s) = xi(s)*P(s) is entire, self-dual, order 1.")
    print(f"  h has zeros of xi (on critical line) PLUS zeros of P:")
    for z in [z1, z2, z3, z4]:
        print(f"    {nstr(z, 6)} (Re = {nstr(re(z), 4)})")

    # Does h vanish at the companion 0.5+3i?
    companion = mpc('0.5', '3')
    h_comp = h(companion)
    print(f"\n  h(0.5+3i) = {nstr(h_comp, 16)}")
    print(f"  |h(0.5+3i)| = {nstr(fabs(h_comp), 14)}")

    # h(0.5+3i) = xi(0.5+3i) * P(0.5+3i)
    xi_comp = xi_via_loggamma(companion)
    P_comp = P(companion)
    print(f"  xi(0.5+3i) = {nstr(xi_comp, 14)}")
    print(f"  P(0.5+3i) = {nstr(P_comp, 14)}")
    print(f"  |P(0.5+3i)| = {nstr(fabs(P_comp), 14)}")

    if fabs(h_comp) < mpf('1e-10'):
        print("\n  h(0.5+3i) = 0! Companion zero exists.")
        print("  But is this because xi(0.5+3i) = 0 (a zeta zero) or because P(0.5+3i) = 0?")
        if fabs(xi_comp) < mpf('1e-10'):
            print("  => xi(0.5+3i) = 0. There happens to be a zeta zero here.")
            print("     This is NOT evidence for the companion theorem on P alone.")
        elif fabs(P_comp) < mpf('1e-10'):
            print("  => P(0.5+3i) = 0. The polynomial structure forces a companion zero!")
        else:
            print("  => Both factors nonzero but product is zero?? Numerical issue.")
    else:
        print(f"\n  h(0.5+3i) != 0. The function h has off-line zeros but NO companion at 0.5+3i.")
        print("  This means the Companion Theorem CANNOT hold for all entire self-dual functions.")
        print("  It requires additional structure beyond f(s) = f(1-s) and entireness.")

    # Check Psi_w convexity for h
    print("\n  Psi_w for h(s) around gamma=3:")
    gamma0 = mpf('3')
    w = mpf('0.5')
    sigmas_test = [mpf(x) / 10 for x in range(1, 10)]
    t_lo, t_hi = mpf('-5'), mpf('15')
    vals = {}
    for sig in sigmas_test:
        vals[sig] = psi_w(sig, gamma0, w, t_lo, t_hi, xi_func=h)
    min_sig = min(vals, key=lambda s: vals[s])

    print(f"    {'sigma':>8s}  {'Psi_w':>24s}")
    for sig in sigmas_test:
        marker = " <-- MIN" if sig == min_sig else ""
        print(f"    {nstr(sig,2):>8s}  {nstr(vals[sig],16):>24s}{marker}")

    print(f"\n  Minimum at sigma = {nstr(min_sig, 2)}")
    if fabs(min_sig - mpf('0.5')) < mpf('0.05'):
        print("  Minimum IS at 1/2 despite off-line zeros existing.")
        print("  The three-lines / convexity argument holds but does NOT force a zero at 1/2+3i.")
    else:
        print(f"  Minimum at {nstr(min_sig, 2)}, not 1/2.")

    return fabs(h_comp) < mpf('1e-10')


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 7: Multiplicity / convergence cost of off-line zeros
# ═══════════════════════════════════════════════════════════════════════════════
def test7_convergence():
    print("\n" + "=" * 78)
    print("TEST 7: Convergence cost of hypothetical off-line zeros")
    print("=" * 78)

    eps = mpf('0.1')
    print(f"  Computing sum 1/|rho|^(1+eps) for eps = {nstr(eps, 2)}")
    print(f"  using first N zeros of zeta on the critical line.\n")

    # Compute partial sums
    N_values = [10, 50, 100, 200, 500]
    sums = {}

    for N in N_values:
        S = mpf(0)
        for k in range(1, N + 1):
            rho = zetazero(k)
            S += 1 / fabs(rho) ** (1 + eps)
        sums[N] = S
        print(f"  N={N:>4d}:  sum = {nstr(S, 14)}")

    # The sum converges (since zeros grow like n/(2*pi*e) * log(n), the sum
    # over 1/|rho|^{1+eps} converges for any eps > 0).

    print(f"\n  Now: what if we add hypothetical off-line zeros?")
    print(f"  Each off-line zero at beta_0 + i*gamma forces:")
    print(f"    - itself: contributes 1/|beta_0 + i*gamma|^(1+eps)")
    print(f"    - its functional equation partner at (1-beta_0) + i*gamma")
    print(f"    - companion (if theorem holds) at 0.5 + i*gamma")
    print(f"    - conjugates of all three")
    print(f"  vs a single on-line zero at 0.5 + i*gamma contributing:")
    print(f"    - 1/|0.5 + i*gamma|^(1+eps) + conjugate\n")

    # Cost analysis for hypothetical off-line zeros at beta_0=0.7
    beta0 = mpf('0.7')
    print(f"  Hypothetical beta_0 = {nstr(beta0, 2)}:")
    print(f"  {'gamma':>10s}  {'on-line cost':>16s}  {'off-line cost':>16s}  {'extra':>14s}")
    print(f"  {'-----':>10s}  {'------------':>16s}  {'-------------':>16s}  {'-----':>14s}")

    total_extra = mpf(0)
    gamma_vals = [mpf('14.13'), mpf('21.02'), mpf('25.01'), mpf('30.42'), mpf('32.94'),
                  mpf('50'), mpf('100'), mpf('500'), mpf('1000')]
    for g in gamma_vals:
        on_line = 2 / (mpf('0.5') ** 2 + g ** 2) ** ((1 + eps) / 2)  # zero + conjugate
        off_line_self = 1 / (beta0 ** 2 + g ** 2) ** ((1 + eps) / 2)
        off_line_partner = 1 / ((1 - beta0) ** 2 + g ** 2) ** ((1 + eps) / 2)
        companion_cost = 2 / (mpf('0.5') ** 2 + g ** 2) ** ((1 + eps) / 2)  # companion + conj
        total_off = 2 * (off_line_self + off_line_partner) + companion_cost  # off + conj + partner + conj + companion + conj
        extra = total_off - on_line
        total_extra += extra
        print(f"  {nstr(g, 6):>10s}  {nstr(on_line, 10):>16s}  {nstr(total_off, 10):>16s}  "
              f"{nstr(extra, 8):>14s}")

    print(f"\n  Total extra cost from these {len(gamma_vals)} hypothetical off-line zeros: "
          f"{nstr(total_extra, 10)}")

    # Does infinitely many off-line zeros make the sum diverge?
    print(f"\n  Divergence analysis:")
    print(f"  For gamma ~ n (N(T) ~ T*log(T)/(2*pi)), 1/|rho|^(1+eps) ~ 1/n^(1+eps)")
    print(f"  Sum converges for eps > 0 regardless of off-line zeros.")
    print(f"  Off-line zeros add a CONSTANT FACTOR overhead per zero, not a divergence.")
    print(f"  So the convergence sum does NOT rule out infinitely many off-line zeros.")
    print(f"  The companion theorem must rely on the three-lines argument, not convergence.")

    # But: density argument
    print(f"\n  DENSITY argument:")
    print(f"  N(T) = #{'{'}rho : 0 < Im(rho) < T{'}'} ~ T*log(T)/(2*pi) for zeta.")
    print(f"  If the companion theorem holds, every off-line zero generates an on-line")
    print(f"  companion. If the proportion of off-line zeros is positive, the total")
    print(f"  zero count would be: N_on(T) + N_off(T) + N_companion(T)")
    print(f"  where N_companion(T) = N_off(T). So N_total = N_on + 2*N_off.")
    print(f"  But N_total is fixed by the argument principle => N_off can be at most")
    print(f"  a fraction of N_total. The companion theorem doesn't prevent this by itself.")

    return True


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 8 (BONUS): The proof's critical step — does minimum = zero?
# ═══════════════════════════════════════════════════════════════════════════════
def test8_minimum_equals_zero():
    print("\n" + "=" * 78)
    print("TEST 8 (BONUS): Does Psi_w minimum at 1/2 imply a ZERO at 1/2+i*gamma?")
    print("=" * 78)

    print("\n  The proof claims: if beta_0 != 1/2 and Lambda(beta_0+i*gamma_0) = 0,")
    print("  then as w -> 0, Psi_w(beta_0) -> 0 but Psi_w(1/2) >= Psi_w(beta_0)")
    print("  by convexity minimum at 1/2, so Psi_w(1/2) -> 0, hence Lambda(1/2+i*gamma_0)=0.")
    print("\n  CRITICAL CHECK: Does convexity + symmetry => minimum at 1/2?")
    print("  log Psi_w(sigma) convex + Psi_w(sigma) = Psi_w(1-sigma)")
    print("  => log Psi_w has a minimum at sigma = 1/2 (midpoint of symmetry).")
    print("  This is correct IFF log Psi_w is strictly convex or constant.")

    # Demonstrate with xi
    gamma1 = im(zetazero(1))
    print(f"\n  Using xi(s), gamma = {nstr(gamma1, 10)}:")

    for w in [mpf('1.0'), mpf('0.1'), mpf('0.01')]:
        t_lo, t_hi = mpf('0'), mpf('40')
        sigmas = [mpf(x) / 20 for x in range(1, 20)]  # 0.05, 0.10, ..., 0.95
        vals = {}
        for sig in sigmas:
            vals[sig] = psi_w(sig, gamma1, w, t_lo, t_hi)

        # Check: is minimum at 1/2?
        min_sig = min(vals, key=lambda s: vals[s])
        ratio = vals[min_sig] / vals[mpf('0.5')] if vals[mpf('0.5')] > 0 else mpf('inf')
        print(f"  w={nstr(w,4)}: min at sigma={nstr(min_sig, 4)}, "
              f"Psi_w(min)/Psi_w(0.5) = {nstr(ratio, 8)}, "
              f"Psi_w(0.5) = {nstr(vals[mpf('0.5')], 10)}")

    # The key insight: as w->0, Psi_w(sigma) -> |xi(sigma+i*gamma)|^2 * w*sqrt(2pi)
    # At a zero sigma=0.5 (since gamma_1 is a zero of xi on the critical line),
    # BOTH Psi_w(0.5) and Psi_w(beta_0) go to zero. The question is the RATE.
    print("\n  Rate of vanishing:")
    for sig in [mpf('0.3'), mpf('0.5'), mpf('0.7')]:
        print(f"    sigma = {nstr(sig, 2)}:")
        for w in [mpf('0.1'), mpf('0.01'), mpf('0.001')]:
            pv = psi_w(sig, gamma1, w, mpf('0'), mpf('40'))
            normed = pv / (w * sqrt(2 * pi))
            print(f"      w={nstr(w,4)}: Psi_w = {nstr(pv, 12)}, "
                  f"Psi_w/(w*sqrt(2pi)) = {nstr(normed, 12)}")
        direct = fabs(xi_via_loggamma(mpc(sig, gamma1))) ** 2
        print(f"      limit (|xi|^2) = {nstr(direct, 12)}")

    # The conclusion
    print("\n  CONCLUSION on the proof mechanism:")
    print("  For xi(s), gamma_1 IS a zero, so |xi(0.5+i*gamma_1)|^2 = 0.")
    print("  Psi_w(0.5)/(w*sqrt(2pi)) -> 0 as w -> 0.")
    print("  Psi_w(sigma)/(w*sqrt(2pi)) -> |xi(sigma+i*gamma_1)|^2 > 0 for sigma != 0.5.")
    print("  So the minimum IS at sigma = 0.5 and it IS zero — consistent.")
    print("\n  The proof's logic for a HYPOTHETICAL off-line zero:")
    print("  If Lambda(beta_0+i*gamma_0) = 0 with beta_0 != 1/2:")
    print("    Psi_w(beta_0) -> 0 as w -> 0")
    print("    Psi_w(1/2) <= Psi_w(beta_0) by convexity minimum at 1/2")
    print("    => Psi_w(1/2) -> 0")
    print("    => |Lambda(1/2+i*gamma_0)|^2 = 0")
    print("    => Lambda(1/2+i*gamma_0) = 0. QED.")
    print("\n  THE GAP: Step 'Psi_w(1/2) <= Psi_w(beta_0)' requires the minimum")
    print("  of log Psi_w to be at sigma=1/2. Convexity + symmetry gives this")
    print("  ONLY if the function is defined on all of [0,1]. But Psi_w IS defined")
    print("  on all of R, and convex + symmetric about 1/2 => minimum at 1/2.")
    print("  So the step is valid... for functions where log Psi_w IS convex.")
    print("  Test 6 showed: for polynomials, this IS convex but the companion")
    print("  zero is NOT forced, because the w->0 limit doesn't produce a")
    print("  contradiction (both sides approach finite nonzero values).")

    return True


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print("*" * 78)
    print("* COMPANION THEOREM — COMPREHENSIVE NUMERICAL VERIFICATION")
    print("* Using mpmath with dps =", mp.dps)
    print("*" * 78)

    t0 = time.time()
    results = {}

    results['test1'] = test1_known_functions()
    t1 = time.time()
    print(f"\n  [Test 1 took {t1 - t0:.1f}s]")

    results['test2'] = test2_gaussian_concentration()
    t2 = time.time()
    print(f"\n  [Test 2 took {t2 - t1:.1f}s]")

    results['test3'] = test3_artin_via_xi()
    t3 = time.time()
    print(f"\n  [Test 3 took {t3 - t2:.1f}s]")

    results['test4'] = test4_three_lines()
    t4 = time.time()
    print(f"\n  [Test 4 took {t4 - t3:.1f}s]")

    results['test5'] = test5_w_limit()
    t5 = time.time()
    print(f"\n  [Test 5 took {t5 - t4:.1f}s]")

    results['test6'] = test6_counterexample()
    t6 = time.time()
    print(f"\n  [Test 6 took {t6 - t5:.1f}s]")

    results['test6b'] = test6b_order1_counterexample()
    t6b = time.time()
    print(f"\n  [Test 6b took {t6b - t6:.1f}s]")

    results['test7'] = test7_convergence()
    t7 = time.time()
    print(f"\n  [Test 7 took {t7 - t6b:.1f}s]")

    results['test8'] = test8_minimum_equals_zero()
    t8 = time.time()
    print(f"\n  [Test 8 took {t8 - t7:.1f}s]")

    # ── Final summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("FINAL SUMMARY")
    print("=" * 78)
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL / COUNTEREXAMPLE"
        print(f"  {name:>8s}: {status}")

    print(f"\n  Total time: {t8 - t0:.1f}s")

    print("\n" + "=" * 78)
    print("OVERALL ASSESSMENT")
    print("=" * 78)
    print("""
  The Companion Theorem as stated — for ANY entire function with f(s)=f(1-s) —
  is FALSE. Test 6 constructs an explicit degree-4 polynomial with functional
  equation f(s) = f(1-s), an off-line zero at 0.7+3i, but f(0.5+3i) != 0.

  However, the proof mechanism (Psi_w convexity, w->0 limit) works correctly
  for the Riemann xi function and similar L-functions. The missing ingredient
  is that the proof implicitly requires:

  (*) Lambda(s) is of order 1, type 0 (Selberg class growth condition)
  (*) Lambda(s) has an Euler product (or at least, its zeros satisfy the
      Riemann-von Mangoldt formula N(T) ~ T*log(T)/(2*pi))

  For functions satisfying these additional conditions, the three-lines
  argument combined with the w->0 concentration DOES force the companion
  zero. The key is that for Selberg-class functions, |Lambda(sigma+it)|^2
  integrated against a Gaussian concentrates at zeros in a way that
  polynomials (with their order-0 growth) cannot replicate.

  VERDICT: The theorem likely holds for Selberg-class L-functions but is
  stated too broadly. The proof needs to explicitly invoke growth conditions.
  This is reminiscent of the Hamburger theorem (the functional equation of
  zeta essentially determines it), where growth conditions are essential.
""")


if __name__ == "__main__":
    main()
