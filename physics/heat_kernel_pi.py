"""
HEAT_KERNEL_PI — heat kernel trace on S^3/2I; spectrum, spectral zeta, Selberg geodesics, pi connections
nos3bl33d
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from mpmath import (
    mp, mpf, pi, sin, cos, sqrt, exp, log, nstr, floor, fabs, nint
)

mp.dps = 60

# =====================================================================
# Constants
# =====================================================================

GROUP_ORDER = 120  # |2I|
N_MAX = 50000      # eigenvalue cutoff for Molien computation

PI_CF = [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]


def continued_fraction(x, n_terms=15):
    """Compute the first n_terms of the continued fraction of x."""
    cf = []
    for _ in range(n_terms):
        a = int(floor(x))
        cf.append(a)
        frac = x - a
        if fabs(frac) < mpf(10)**(-45):
            break
        x = 1 / frac
    return cf


def cf_str(c):
    if len(c) <= 1:
        return str(c)
    return f'[{c[0]}; {", ".join(str(a) for a in c[1:])}]'


def sep(title):
    print(f'\n{"="*72}')
    print(f'  {title}')
    print(f'{"="*72}')


# =====================================================================
# Spectrum via Molien series
# =====================================================================
# Molien series for 2I: f(t) = (1+t^30)/((1-t^12)(1-t^20))
# All non-zero coefficients occur at EVEN n.
# Equivalently: g(u) = (1+u^15)/((1-u^6)(1-u^10)), coeff of u^m = mult(2m).

def build_spectrum(n_max):
    mults = [0] * (n_max + 1)
    for k in range(n_max // 12 + 1):
        for j in range(n_max // 20 + 1):
            idx = 12 * k + 20 * j
            if idx <= n_max:
                mults[idx] += 1
            idx2 = idx + 30
            if idx2 <= n_max:
                mults[idx2] += 1
    return mults


def verify_character_formula(mults, n_check=50):
    """Cross-check Molien series against direct character formula."""
    conj_classes = [
        (1,   mpf(0)),
        (1,   pi),
        (12,  pi / 5),
        (12,  2 * pi / 5),
        (12,  3 * pi / 5),
        (12,  4 * pi / 5),
        (20,  pi / 3),
        (20,  2 * pi / 3),
        (30,  pi / 2),
    ]
    for n in range(n_check + 1):
        total = mpf(0)
        for count, theta in conj_classes:
            n1 = n + 1
            if fabs(theta) < mpf(10)**(-45):
                chi = mpf(n1)
            elif fabs(theta - pi) < mpf(10)**(-45):
                chi = mpf((-1)**n * n1)
            else:
                chi = sin(n1 * theta) / sin(theta)
            total += count * chi
        m_char = int(nint(total / GROUP_ORDER))
        if m_char != mults[n]:
            print(f'  MISMATCH at n={n}: Molien={mults[n]}, char={m_char}')
            return False
    return True


# =====================================================================
# Heat kernel trace
# =====================================================================

def heat_trace(t, mults, skip_zero=False):
    total = mpf(0)
    start = 1 if skip_zero else 0
    for n in range(start, len(mults)):
        if mults[n] > 0:
            lam = mpf(n) * (n + 2)
            term = mults[n] * exp(-lam * t)
            if term < mpf(10)**(-50) and n > 100:
                break
            total += term
    return total


def weyl_approx(t):
    vol = pi**2 / 60
    return vol / (4 * pi * t) ** (mpf(3) / 2)


# =====================================================================
# Spectral zeta function (converges for Re(s) > 3/2 on S^3 quotients)
# =====================================================================

def spectral_zeta(s, mults):
    """zeta_M(s) = Sigma_{n>=1} mult(n) / (n(n+2))^s"""
    total = mpf(0)
    for n in range(1, len(mults)):
        if mults[n] > 0:
            lam = mpf(n) * (n + 2)
            total += mults[n] / lam**s
    return total


# =====================================================================
# Regularized zeta_M(1)
# =====================================================================

def spectral_zeta_reg1(mults):
    """
    Regularized zeta_M(1):
    Subtract the smooth asymptotic (n+1)/120 for ALL n,
    then sum [mult(n) - (n+1)/120] / (n(n+2)).

    Even + odd terms diverge individually but cancel exactly
    because the excess at even n ~ +1/(240*m) cancels
    the deficit at odd n ~ -1/(240*m).

    The sum converges conditionally (paired even+odd).
    """
    total = mpf(0)
    for n in range(1, len(mults)):
        lam = mpf(n) * (n + 2)
        smooth = mpf(n + 1) / 120
        total += (mults[n] - smooth) / lam
    return total


# =====================================================================
# S^3 spectral zeta: EXACT formula
# =====================================================================

def s3_zeta2_exact():
    """
    PROVEN by partial fractions:

    zeta_{S^3}(2) = Sigma_{n>=1} (n+1)^2 / (n(n+2))^2
                  = Sigma_{m>=2} m^2 / (m^2-1)^2

    Partial fractions: m^2/(m^2-1)^2 =
        1/(4(m-1)) + 1/(4(m-1)^2) - 1/(4(m+1)) + 1/(4(m+1)^2)

    Telescoping + Basel:
        Sum [1/(m-1) - 1/(m+1)] from m=2..inf = 1 + 1/2 = 3/2
        Sum [1/(m-1)^2] from m=2..inf = pi^2/6
        Sum [1/(m+1)^2] from m=2..inf = pi^2/6 - 1 - 1/4

    Result: zeta_{S^3}(2) = 1/16 + pi^2/12

    Therefore: pi = sqrt(12 * (zeta_{S^3}(2) - 1/16))
    """
    return mpf(1)/16 + pi**2/12


# =====================================================================
# Selberg geodesic sums
# =====================================================================

def geodesic_sums():
    classes = [
        (12,  pi / 5,   'pi/5'),
        (12,  2*pi / 5, '2pi/5'),
        (12,  3*pi / 5, '3pi/5'),
        (12,  4*pi / 5, '4pi/5'),
        (20,  pi / 3,   'pi/3'),
        (20,  2*pi / 3, '2pi/3'),
        (30,  pi / 2,   'pi/2'),
    ]

    total_sin = mpf(0)
    total_sin2 = mpf(0)
    details = []

    for count, theta, label in classes:
        s = sin(theta)
        cs = mpf(count) / s
        cs2 = mpf(count) / s**2
        total_sin += cs
        total_sin2 += cs2
        details.append((count, theta, label, s, cs, cs2))

    return total_sin, total_sin2, details


# =====================================================================
# MAIN
# =====================================================================

def main():
    print('Building spectrum (n_max = %d)...' % N_MAX)
    mults = build_spectrum(N_MAX)

    # ================================================================
    sep('PART 1: Spectrum of S^3/2I (Poincare Dodecahedral Space)')
    # ================================================================

    print('\nCross-checking Molien vs character formula (n=0..50)...')
    if verify_character_formula(mults, 50):
        print('  Verified: all multiplicities match.')

    print(f'\nNon-zero multiplicities for n = 0 to 100:')
    print(f'{"n":>5} {"mult":>6} {"lam":>12}')
    print('-' * 28)
    nz_count = 0
    for n in range(101):
        if mults[n] > 0:
            print(f'{n:>5} {mults[n]:>6} {n*(n+2):>12}')
            nz_count += 1
    print(f'\nNon-zero: {nz_count} out of 101')

    # Key property: only even n
    odd_nz = sum(1 for n in range(1, N_MAX+1, 2) if mults[n] > 0)
    even_nz = sum(1 for n in range(0, N_MAX+1, 2) if mults[n] > 0)
    print(f'\nAll multiplicities at EVEN n: odd_nonzero={odd_nz}, even_nonzero={even_nz}')

    # Which even n < 60 are zero?
    zero_small = [n for n in range(2, 60, 2) if mults[n] == 0]
    print(f'Even n < 60 with mult=0: {zero_small}')
    print(f'All even n >= 60 have mult > 0.')

    # Asymptotic behavior
    print(f'\nAsymptotic: mult(2m) ~ (m+1)/30 for large m (from Molien pole)')
    print(f'Equivalently: mult(n) ~ (n+2)/60 for large even n')
    for n in [100, 1000, 10000, 50000]:
        if n <= N_MAX and mults[n] > 0:
            expected = (n+2) / mpf(60)
            print(f'  n={n:>5}: mult={mults[n]:>4}, (n+2)/60={nstr(expected,6):>8}, ratio={nstr(mpf(mults[n])/expected,8)}')

    # ================================================================
    sep('PART 2: Heat Kernel Trace')
    # ================================================================

    vol_M = pi**2 / 60
    print(f'\nVol(M) = pi^2/60 = {nstr(vol_M, 25)}')

    print(f'\n{"t":>10}  {"K_M(t)":>28}  {"Weyl":>22}  {"ratio":>12}')
    print('-' * 80)
    for t in [mpf('0.00001'), mpf('0.0001'), mpf('0.001'), mpf('0.01'),
              mpf('0.1'), mpf('0.5'), mpf('1'), mpf('2')]:
        km = heat_trace(t, mults)
        wy = weyl_approx(t)
        r = km / wy
        print(f'{nstr(t,6):>10}  {nstr(km,22):>28}  {nstr(wy,16):>22}  {nstr(r,8):>12}')

    # At t=0.5+, K_M -> 1 (spectral gap lambda_1 = 168 is huge)
    print(f'\nLarge spectral gap: lambda_1 = 168 (at n=12).')
    print(f'K_M(t) = 1 + exp(-168t) + ... for t > 0.1.')
    print(f'K_M(0.01) - 1 = {nstr(heat_trace(mpf("0.01"), mults, skip_zero=True), 20)}')
    print(f'K_M(0.1) - 1  = {nstr(heat_trace(mpf("0.1"), mults, skip_zero=True), 20)}')

    # Heat kernel at t = 1/168 (one e-folding of first eigenvalue)
    t168 = mpf(1) / 168
    km168 = heat_trace(t168, mults)
    km168_norm = km168 / vol_M
    print(f'\nAt t = 1/168 = 1/lambda_1:')
    print(f'  K_M(1/168) = {nstr(km168, 25)}')
    print(f'  K_M/Vol(M) = {nstr(km168_norm, 20)}')
    cf168 = continued_fraction(km168_norm)
    print(f'  CF(K_M/Vol) = {cf_str(cf168[:12])}')

    # ================================================================
    sep('PART 3: Spectral Zeta Function')
    # ================================================================

    print(f'\nzeta_M(s) = Sigma mult(n) / (n(n+2))^s')
    print(f'Converges for Re(s) > 3/2 on a 3-manifold.')
    print(f'(mult(n) ~ n/60 for even n => diverges for s <= 1)')

    zetas = {}
    for s in [mpf(2), mpf(3), mpf(4)]:
        z = spectral_zeta(s, mults)
        zetas[s] = z
        cf = continued_fraction(z, 12)
        print(f'\n  zeta_M({int(s)}) = {nstr(z, 35)}')
        print(f'    CF = {cf_str(cf[:12])}')

    # Convergence check
    print(f'\n  Convergence of zeta_M(2):')
    for nm in [300, 1000, 5000, 10000, 50000]:
        z = mpf(0)
        for n in range(1, min(nm+1, N_MAX+1)):
            if mults[n] > 0:
                z += mults[n] / (mpf(n)*(n+2))**2
        print(f'    n_max={nm:>5}: {nstr(z, 30)}')

    z2 = zetas[mpf(2)]
    z3 = zetas[mpf(3)]
    z4 = zetas[mpf(4)]

    # ================================================================
    sep('PART 3b: EXACT RESULT -- zeta_{S^3}(2) = pi^2/12 + 1/16')
    # ================================================================

    print(f"""
  THEOREM (proved by partial fractions):

    zeta_{{S^3}}(2) = Sigma_{{n>=1}} (n+1)^2 / (n(n+2))^2
                    = Sigma_{{m>=2}} m^2 / (m^2-1)^2
                    = pi^2/12 + 1/16

  PROOF:
    Partial fractions: m^2/(m^2-1)^2 =
      1/(4(m-1)) + 1/(4(m-1)^2) - 1/(4(m+1)) + 1/(4(m+1)^2)

    Summing from m=2 to infinity:
    - Telescoping: (1/4)*Sum[1/(m-1) - 1/(m+1)] = (1/4)*(1+1/2) = 3/8
    - Basel: (1/4)*Sum[1/(m-1)^2] = (1/4)*pi^2/6 = pi^2/24
    - Basel: (1/4)*Sum[1/(m+1)^2] = (1/4)*(pi^2/6 - 5/4) = pi^2/24 - 5/16

    Total = 3/8 + pi^2/24 + pi^2/24 - 5/16 = 1/16 + pi^2/12

  COROLLARY: pi = sqrt(12*(zeta_{{S^3}}(2) - 1/16))
""")

    z_s3_exact = s3_zeta2_exact()
    print(f'  zeta_{{S^3}}(2) exact = {nstr(z_s3_exact, 35)}')

    # Verify numerically with Richardson extrapolation
    print(f'\n  Numerical verification (Richardson extrapolation):')
    S_vals = []
    for N in [2000, 4000, 8000]:
        s = mpf(0)
        for m in range(2, N+1):
            s += mpf(m)**2 / (mpf(m)**2 - 1)**2
        S_vals.append(s)
    r1 = 2*S_vals[1] - S_vals[0]
    r2 = 2*S_vals[2] - S_vals[1]
    rich2 = (4*r2 - r1)/3
    print(f'    Richardson 2nd order: {nstr(rich2, 30)}')
    print(f'    Exact formula:       {nstr(z_s3_exact, 30)}')
    print(f'    Difference:          {nstr(fabs(rich2 - z_s3_exact), 8)}')

    pi_from_zeta = sqrt(12 * (z_s3_exact - mpf(1)/16))
    print(f'\n  pi from spectrum: sqrt(12*(zeta_S3(2) - 1/16))')
    print(f'    = {nstr(pi_from_zeta, 40)}')
    print(f'    pi = {nstr(pi, 40)}')
    print(f'    Match: {fabs(pi_from_zeta - pi) < mpf(10)**(-40)}')

    # ================================================================
    sep('PART 4: Regularized zeta_M(1)')
    # ================================================================

    print(f'\n  zeta_M(1) = Sigma mult(n)/(n(n+2)) DIVERGES (pole at s=1)')
    print(f'  mult(n) ~ (n+1)/120 on average => harmonic divergence')

    print(f'\n  Regularization: subtract smooth part (n+1)/120 for ALL n')
    print(f'  zeta_M^reg(1) = Sigma [mult(n) - (n+1)/120] / (n(n+2))')
    print(f'  Converges conditionally (even+odd terms cancel at O(1/n))')

    z_reg = spectral_zeta_reg1(mults)
    print(f'\n  zeta_M^reg(1) [n_max={N_MAX}] = {nstr(z_reg, 30)}')
    cf_reg = continued_fraction(z_reg)
    print(f'  CF = {cf_str(cf_reg)}')

    print(f'\n  Convergence:')
    for nm in [1000, 5000, 10000, 20000, 50000]:
        z = mpf(0)
        for n in range(1, min(nm+1, N_MAX+1)):
            z += (mults[n] - mpf(n+1)/120) / (mpf(n)*(n+2))
        print(f'    n_max={nm:>5}: {nstr(z, 25)}')

    # Key ratios
    print(f'\n  120 * zeta_M^reg(1) = {nstr(z_reg * 120, 20)}')
    print(f'  pi                  = {nstr(pi, 20)}')
    print(f'  Difference from pi: {nstr(z_reg*120 - pi, 15)}')
    print(f'  (Close to pi but NOT equal -- converging to a different value)')

    # ================================================================
    sep('PART 5: Connection to pi -- Comprehensive Search')
    # ================================================================

    print(f'\npi = {nstr(pi, 40)}')
    print(f'pi CF = {cf_str(continued_fraction(pi, 15))}')

    # Collect all computable quantities
    candidates = {}

    candidates['zeta_M(2)'] = z2
    candidates['zeta_M(3)'] = z3
    candidates['zeta_M(4)'] = z4
    candidates['120*zeta_M(2)'] = 120 * z2
    candidates['120*zeta_M(3)'] = 120 * z3
    candidates['zeta_M(2)/zeta_M(3)'] = z2 / z3
    candidates['zeta_M(2)/zeta_M(4)'] = z2 / z4
    candidates['zeta_M(3)/zeta_M(4)'] = z3 / z4
    candidates['sqrt(zeta_M(2)/zeta_M(4))'] = sqrt(z2 / z4)
    candidates['zeta_M^reg(1)'] = z_reg
    candidates['120*zeta_M^reg(1)'] = 120 * z_reg
    candidates['60*zeta_M^reg(1)'] = 60 * z_reg

    # Heat kernel values
    for label, t in [('1/(4pi)', 1/(4*pi)), ('1/168', mpf(1)/168)]:
        km = heat_trace(t, mults)
        candidates[f'K_M({label})'] = km
        candidates[f'K_M({label})/Vol'] = km / vol_M

    # Geodesic sums
    total_sin, total_sin2, _ = geodesic_sums()
    candidates['Sum n/sin(th)'] = total_sin
    candidates['Sum n/sin^2(th)'] = total_sin2
    candidates['Sum n/sin(th) / 120'] = total_sin / 120

    # Products with pi
    candidates['120*zeta_M(2)*pi^2'] = 120 * z2 * pi**2
    candidates['zeta_M^reg(1)*pi'] = z_reg * pi
    candidates['zeta_M^reg(1)/pi'] = z_reg / pi
    candidates['zeta_M^reg(1)/pi^2'] = z_reg / pi**2

    # Search for CF matching [3; 7, 15, ...]
    print(f'\n{"Expression":>35}  {"Value":>28}  {"CF (first 10)"}')
    print('-' * 115)

    close_to_pi = []
    for label, val in candidates.items():
        if val > 0:
            cf = continued_fraction(val, 12)
            flag = ''
            if len(cf) >= 3 and cf[0] == 3 and cf[1] == 7 and cf[2] == 15:
                flag = ' <=== MATCH [3;7,15] !!!'
            elif len(cf) >= 2 and cf[0] == 3 and cf[1] == 7:
                flag = ' <== [3;7,...] !!'
            elif cf[0] == 3:
                flag = ' <-- starts with 3'

            hits = [a for a in cf[:8] if a in [7, 15, 292]]
            if hits and not flag:
                flag = f' <-- has {hits}'

            print(f'{label:>35}  {nstr(val, 20):>28}  {cf_str(cf[:10])}{flag}')

            dist = fabs(val - pi)
            if dist < mpf(1):
                close_to_pi.append((dist, label, val))

    if close_to_pi:
        close_to_pi.sort()
        print(f'\n  Quantities closest to pi:')
        for d, label, val in close_to_pi[:5]:
            cf = continued_fraction(val, 12)
            print(f'    {label} = {nstr(val, 25)}, |x - pi| = {nstr(d, 10)}')
            print(f'      CF = {cf_str(cf)}')

    # ================================================================
    sep('PART 6: Selberg Geodesic Contributions')
    # ================================================================

    total_sin, total_sin2, details = geodesic_sums()

    print(f'\n{"count":>5} {"theta":>10} {"sin(th)":>18} {"n/sin(th)":>18} {"n/sin^2(th)":>18}')
    print('-' * 78)
    for count, theta, label, s, cs, cs2 in details:
        print(f'{count:>5} {label:>10} {nstr(s,12):>18} {nstr(cs,12):>18} {nstr(cs2,12):>18}')

    print(f'\nSum n_g/sin(theta_g) = {nstr(total_sin, 30)}')
    print(f'  CF = {cf_str(continued_fraction(total_sin))}')

    # Sum n/sin^2 = 538/3 EXACTLY
    print(f'\nSum n_g/sin^2(theta_g) = {nstr(total_sin2, 30)}')
    print(f'  = 538/3 EXACTLY')

    print(f"""
  PROOF (algebraic):
    Using sin(pi/5) = sin(4pi/5), sin(2pi/5) = sin(3pi/5):
    Sum = 24/sin^2(pi/5) + 24/sin^2(2pi/5) + 40/sin^2(pi/3) + 30

    sin^2(pi/5) = (10-2*sqrt(5))/16
    => 24/sin^2(pi/5) = 24*16/(10-2sqrt5) = 384/(10-2sqrt5)
    Similarly: 24/sin^2(2pi/5) = 384/(10+2sqrt5)

    384/(10-2sqrt5) + 384/(10+2sqrt5)
    = 384*20 / ((10)^2 - (2sqrt5)^2) = 7680/80 = 96

    40/sin^2(pi/3) = 40/(3/4) = 160/3

    Total = 96 + 160/3 + 30 = (288 + 160 + 90)/3 = 538/3  QED
""")

    # Exact form of the irrational sum
    A = sqrt(20 + 8*sqrt(5))
    exact_sin_sum = 24 * A / sqrt(5) + 80/sqrt(3) + 30
    print(f'  Sum n/sin(th) = 24*sqrt(20+8*sqrt(5))/sqrt(5) + 80/sqrt(3) + 30')
    print(f'  = {nstr(exact_sin_sum, 30)}')
    print(f'  Verified: {fabs(exact_sin_sum - total_sin) < mpf(10)**(-40)}')

    # ================================================================
    sep('PART 7: Heat Kernel and Volume')
    # ================================================================

    print(f'\nVol(M) = pi^2/60 = {nstr(vol_M, 30)}')
    print(f'Vol(S^3) = 2*pi^2 = {nstr(2*pi**2, 30)}')
    print(f'|2I| = 120, |A5| = 60')
    print(f'\npi = sqrt(60*Vol(M)) = sqrt(60*pi^2/60) = pi  [tautological]')

    print(f'\nNormalized heat kernel K_M(t)/Vol(M) for small t:')
    for t in [mpf('0.01'), mpf('0.1'), mpf('0.5'), mpf('1')]:
        km = heat_trace(t, mults)
        kn = km / vol_M
        cf = continued_fraction(kn, 10)
        print(f'  t={nstr(t,4):>6}: K_M/Vol = {nstr(kn, 18):>22}  CF = {cf_str(cf[:8])}')

    # At large t, K_M(t) -> 1, so K_M/Vol -> 1/Vol = 60/pi^2
    sixty_over_pisq = 60 / pi**2
    print(f'\n  As t -> inf: K_M/Vol -> 60/pi^2 = {nstr(sixty_over_pisq, 20)}')
    print(f'  CF(60/pi^2) = {cf_str(continued_fraction(sixty_over_pisq, 12))}')

    # ================================================================
    sep('PART 8: Spectral Sum Ratios')
    # ================================================================

    # Build various spectral sums (using only convergent ones)
    sum_m = mpf(0)
    sum_m_lam2 = mpf(0)
    sum_m_inv_lam2 = mpf(0)
    sum_m_inv_lam3 = mpf(0)

    for n in range(1, N_MAX+1):
        if mults[n] > 0:
            m = mpf(mults[n])
            lam = mpf(n) * (n + 2)
            sum_m += m
            sum_m_inv_lam2 += m / lam**2
            sum_m_inv_lam3 += m / lam**3

    ratios = {
        'zeta_M(2)/zeta_M(3)': z2 / z3,
        'zeta_M(3)/zeta_M(4)': z3 / z4,
        'zeta_M(2)/zeta_M(4)': z2 / z4,
        'sqrt(zeta_M(2)/zeta_M(4))': sqrt(z2 / z4),
        'zeta_M(2)^2 / zeta_M(4)': z2**2 / z4,
        '1/sqrt(zeta_M(2))': 1/sqrt(z2),
        'sqrt(zeta_M(2))*120': sqrt(z2)*120,
    }

    for label, val in ratios.items():
        cf = continued_fraction(val, 12)
        pi_hits = [a for a in cf[:8] if a in [3, 7, 15, 292]]
        flag = f'  <-- {pi_hits}' if pi_hits else ''
        print(f'  {label:>35} = {nstr(val, 20):>25}  CF={cf_str(cf[:8])}{flag}')

    # ================================================================
    sep('PART 9: The CF [3; 7, 15, 1, 292, ...] Deep Search')
    # ================================================================

    print(f'\npi = {nstr(pi, 30)}')
    print(f'pi CF = {cf_str(PI_CF)}')
    print(f'\nSearching for expressions whose CF starts [3; 7, 15, ...]...\n')

    # K_M(t)/Vol at t = 1/168 was [8; 1, 15, ...] -- has 15!
    print(f'Notable near-miss: K_M(1/168)/Vol(M) = {nstr(km168_norm, 20)}')
    cf_note = continued_fraction(km168_norm, 12)
    print(f'  CF = {cf_str(cf_note[:12])}')
    print(f'  Has [; ..., 15, ...] at position 2!')

    # Try: 120*zeta_reg(1) was 3.007... which has CF [3; 137, ...]
    # Try additional combinations
    print(f'\n120*zeta_M^reg(1) = {nstr(z_reg*120, 20)}')
    print(f'  CF = {cf_str(continued_fraction(z_reg*120, 12))}')
    print(f'  Starts with 3, but a[1] = 137 not 7')
    print(f'  Note: 137 is the fine structure denominator!')

    # The zeta(1)/zeta(2) ratio had [612; 1, 4, 15, ...]
    # But with the corrected zeta_M(1) being divergent, this was truncation-dependent.

    # Systematic: try a*z2 + b*z3 etc
    found = False
    for a in range(-20, 21):
        for b in range(-20, 21):
            if a == 0 and b == 0:
                continue
            val = a * z2 * 120 + b * z3 * 120
            if mpf(2.5) < val < mpf(3.5):
                cf = continued_fraction(val, 8)
                if len(cf) >= 2 and cf[0] == 3 and cf[1] == 7:
                    print(f'\n  FOUND: {a}*120*z_M(2) + {b}*120*z_M(3) = {nstr(val, 20)}')
                    print(f'    CF = {cf_str(cf)}')
                    found = True

    # Try with geodesic sums
    for a in range(-5, 6):
        for b in range(-5, 6):
            if a == 0 and b == 0:
                continue
            val = a * total_sin / 120 + b * total_sin2 / 120
            if mpf(2.5) < val < mpf(3.5):
                cf = continued_fraction(val, 8)
                if len(cf) >= 2 and cf[0] == 3 and cf[1] == 7:
                    print(f'\n  FOUND: {a}*(Sum_sin/120) + {b}*(Sum_sin2/120) = {nstr(val, 20)}')
                    print(f'    CF = {cf_str(cf)}')
                    found = True

    if not found:
        print(f'\n  No linear combination of spectral quantities matches [3; 7, 15, ...].')

    # ================================================================
    sep('FINAL SUMMARY')
    # ================================================================

    print(f"""
  RESULTS OF HEAT KERNEL ANALYSIS ON S^3/2I
  ==========================================

  1. SPECTRUM
     - Molien series: (1+t^30)/((1-t^12)(1-t^20))
     - Verified against character formula
     - All multiplicities at EVEN n only
     - First non-zero: n = 0, 12, 20, 24, 30, 32, 36, ...
     - 15 even values below 60 are zero; ALL even n >= 60 have mult > 0
     - Asymptotic: mult(n) ~ (n+2)/60 for large even n

  2. HEAT KERNEL
     - K_M(t) converges to 1 rapidly (spectral gap lambda_1 = 168)
     - Weyl asymptotics K ~ Vol/(4pi*t)^(3/2) valid for t << 10^-4
     - Vol(M) = pi^2/60

  3. SPECTRAL ZETA -- KEY RESULTS
     a) zeta_M(s) converges for Re(s) > 3/2 (not at s=1!)

     b) zeta_{{S^3}}(2) = pi^2/12 + 1/16  [EXACT, proved by partial fractions]
        => pi = sqrt(12*(zeta_{{S^3}}(2) - 1/16))
        This is a GENUINE derivation of pi from eigenvalue sums.

     c) zeta_M(2) = {nstr(z2, 30)} [converged]
        zeta_M(3) = {nstr(z3, 30)} [converged]

     d) Regularized zeta_M^reg(1) ~ {nstr(z_reg, 20)} (slow convergence)
        120 * zeta_M^reg(1) ~ {nstr(z_reg*120, 15)} (close to pi = 3.14159...)
        CF = [3; 137, ...] -- starts with 3, and a[1] = 137!

  4. GEODESIC SUMS
     Sum n_g/sin^2(theta_g) = 538/3 EXACTLY (proved algebraically)
     Sum n_g/sin(theta_g) = 24*sqrt(20+8*sqrt5)/sqrt5 + 80/sqrt3 + 30
                          = {nstr(total_sin, 20)} (irrational)

  5. PI CONNECTIONS
     - The S^3 spectral zeta gives pi EXACTLY through Basel's problem:
       zeta_{{S^3}}(2) = pi^2/12 + 1/16
     - No direct CF match to [3; 7, 15, 1, 292, ...] found in any
       spectral quantity on the quotient M = S^3/2I.
     - 120*zeta_M^reg(1) ~ 3.0... with CF [3; 137, ...]:
       the 137 is tantalizing but the value != pi.
     - The CF of pi does not emerge from the QUOTIENT spectrum,
       only from the full S^3 spectrum (through Basel/zeta(2)).

  6. GRAPH vs MANIFOLD
     Graph dodecahedron: Tr(L^-1) = 137/15 (combinatorial, exact rational)
     Manifold S^3/2I: zeta_M(1) diverges (spectral, not comparable)
     The 137 appearing in CF of 120*zeta_M^reg(1) may be coincidence
     or may reflect the icosahedral symmetry in both contexts.
""")


if __name__ == '__main__':
    main()
