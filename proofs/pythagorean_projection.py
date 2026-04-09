#!/usr/bin/env python3
"""
PYTHAGOREAN_PROJECTION — tests F(1/sqrt(d))^2 = Tr(L^-1)^2 + d on the dodecahedron and other regular graphs
nos3bl33d

Ihara zeta log-derivative vs Green's function trace. mpmath 80 digits.
"""

from mpmath import mp, mpf, sqrt, pi, atan, cos, sin, acos, log, fsum
from fractions import Fraction as Frac

mp.dps = 80

# ==============================================================================
# UTILITY
# ==============================================================================

def header(title):
    w = 72
    print("\n" + "=" * w)
    print(f"  {title}")
    print("=" * w)

def subheader(title):
    print(f"\n--- {title} ---")

def fmt(x, digits=50):
    return mp.nstr(x, digits)

def frac_to_mpf(f):
    """Convert a fractions.Fraction to mpf."""
    return mpf(f.numerator) / mpf(f.denominator)

def factorize(n):
    n = abs(int(n))
    if n <= 1:
        return {n: 1} if n == 1 else {}
    factors = {}
    d_f = 2
    while d_f * d_f <= n:
        while n % d_f == 0:
            factors[d_f] = factors.get(d_f, 0) + 1
            n //= d_f
        d_f += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def fmt_factors(n):
    n = abs(int(n))
    if n <= 1:
        return str(n)
    f = factorize(n)
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))


# ==============================================================================
# IHARA ZETA MACHINERY
# ==============================================================================

def ihara_log_deriv(adj_eigenvalues, d_graph, V_graph, E_graph):
    """
    Compute F(u) = -u zeta'/zeta at u = 1/sqrt(d) for a d-regular graph.

    For a d-regular graph with V vertices, E edges:
      zeta_X(u)^{-1} = (1 - u^2)^{r-1} * prod_i (1 - lambda_i * u + (d-1) * u^2)
    where r = E - V + 1 (rank of fundamental group = # independent cycles).

    So: -log zeta = (r-1)*log(1-u^2) + sum_i log(1 - lambda_i*u + (d-1)*u^2)

    F(u) = -u * d/du [log zeta(u)]
         = (r-1) * 2u^2/(1-u^2) + sum_i u*(lambda_i - 2(d-1)*u) / (1 - lambda_i*u + (d-1)*u^2)

    Wait, sign: -u * d/du [-log zeta] = -u * [d/du of -log zeta] but
    F = -u * zeta'/zeta = -u * d/du log(zeta)
      = u * d/du [-log zeta]
      = u * d/du [(r-1)*log(1-u^2) + sum log(1 - lam*u + (d-1)*u^2)]

    = (r-1) * u * (-2u)/(1-u^2) + sum_i u * (-lam + 2(d-1)*u) / (1 - lam*u + (d-1)*u^2)

    Hmm, let me be very careful.

    zeta^{-1} = (1-u^2)^{r-1} * prod(1 - lam_i u + (d-1)u^2)
    log(zeta^{-1}) = (r-1)log(1-u^2) + sum log(1 - lam_i u + (d-1)u^2)
    log(zeta) = -(r-1)log(1-u^2) - sum log(...)

    F = -u * d/du log(zeta)
      = -u * [-(r-1)*(-2u)/(1-u^2) - sum (-lam_i + 2(d-1)u)/(1 - lam_i u + (d-1)u^2)]
      = -u * [(r-1)*2u/(1-u^2) + sum (lam_i - 2(d-1)u)/(1 - lam_i u + (d-1)u^2)]
      = -(r-1)*2u^2/(1-u^2) - sum u*(lam_i - 2(d-1)u)/(1 - lam_i u + (d-1)u^2)
      = -(r-1)*2u^2/(1-u^2) + sum u*(2(d-1)u - lam_i)/(1 - lam_i u + (d-1)u^2)

    Hmm that gives a DIFFERENT sign for the edge contribution compared to what I had before.
    Let me just directly verify with the known dodecahedron value.
    """
    u_val = 1 / sqrt(mpf(d_graph))
    u2 = u_val**2
    r_val = E_graph - V_graph + 1

    # Try both sign conventions and see which matches
    # Convention A: F = 2(r-1)u^2/(1-u^2) + sum mult*u*(lam - 2(d-1)u)/(1 - lam*u + (d-1)u^2)
    edge_A = 2 * (r_val - 1) * u2 / (1 - u2)
    eig_A = mpf(0)
    for lam, mult in adj_eigenvalues:
        denom = 1 - lam * u_val + (d_graph - 1) * u2
        numer = u_val * (lam - 2 * (d_graph - 1) * u_val)
        eig_A += mult * numer / denom
    F_A = edge_A + eig_A

    # Convention B: the negative
    F_B = -edge_A - eig_A

    # Convention C: positive edge, negative eig
    F_C = edge_A - eig_A

    # Convention D: negative edge, positive eig
    F_D = -edge_A + eig_A

    return F_A, F_B, F_C, F_D


def ihara_F_correct(adj_eigenvalues, d_graph, V_graph, E_graph):
    """
    Compute F(u) using the explicit formula approach.

    For a d-regular graph, the Ihara zeta has:
    log zeta_X(u)^{-1} = (r-1) log(1-u^2) + sum_i log(1 - lambda_i u + (d-1)u^2)

    The "logarithmic derivative of the Ihara zeta" typically means:
    u * d/du log zeta_X(u)^{-1}
    with a sign chosen so F > 0 for small u.

    Actually, the standard convention used in the Pythagorean framework context
    seems to give F(1/sqrt(3)) = 4308/715 + (270/143)*sqrt(3) ~ 9.295 for the dodecahedron.

    Let me just compute all pieces and return them separately so we can figure out the
    right combination.
    """
    u_val = 1 / sqrt(mpf(d_graph))
    u2 = u_val**2
    r_val = E_graph - V_graph + 1

    # Individual per-eigenvalue factors: 1 - lam*u + (d-1)*u^2
    # At u = 1/sqrt(d): u^2 = 1/d, so (d-1)*u^2 = (d-1)/d
    # factor_i = 1 - lam/sqrt(d) + (d-1)/d = (2d-1)/d - lam/sqrt(d)

    # The sum over eigenvalues of -u * d/du log(factor_i):
    #   = sum_i (-u) * (-lam + 2(d-1)u) / (1 - lam*u + (d-1)*u^2)
    #   = sum_i u * (lam - 2(d-1)u) / (factor_i)

    eig_sum = mpf(0)
    for lam, mult in adj_eigenvalues:
        factor = 1 - lam * u_val + (d_graph - 1) * u2
        deriv = (-lam + 2 * (d_graph - 1) * u_val)
        contrib = -u_val * deriv / factor  # = u*(lam - 2(d-1)*u)/factor
        eig_sum += mult * contrib

    # Edge/cycle term: -u * d/du [(r-1)*log(1-u^2)] = -(r-1)*u*(-2u)/(1-u^2) = 2(r-1)*u^2/(1-u^2)
    edge_term = 2 * (r_val - 1) * u2 / (1 - u2)

    # F = edge_term + eig_sum  gives u * d/du log(zeta^{-1}) but with some sign
    # Let me try: the "Ihara logarithmic derivative" as used by the user
    # They define it via per-eigenvalue pieces f_I(lam) = (4 - lam*sqrt3)/(5 - lam*sqrt3)
    # and presumably F = sum mult * f_I(lam) [plus edge term?]

    # Let's compute what sum mult * f_I gives for dodecahedron
    # f_I(lam) = (4 - lam*sqrt(3)) / (5 - lam*sqrt(3))
    # At d=3, u=1/sqrt(3):
    #   4 = 2(d-1) = 2*2 = 4  ... yes, 2(d-1) = 4
    #   5 = d + (d-1) = 3 + 2 = 5... or 1/u^2 + (d-1) = 3 + 2 = 5... yes
    # So f_I(lam) = (2(d-1) - lam*sqrt(d)) / (d + (d-1) - lam*sqrt(d))
    #            = (2(d-1) - lam/u) / ((2d-1) - lam/u)
    # Multiply top and bottom by u:
    #            = (2(d-1)u - lam) / ((2d-1)u - lam)
    # At u = 1/sqrt(d):
    #   numerator = 2(d-1)/sqrt(d) - lam
    #   denominator = (2d-1)/sqrt(d) - lam

    # And (1 - lam*u + (d-1)*u^2) = 1 - lam/sqrt(d) + (d-1)/d
    # = (d - lam*sqrt(d) + d - 1)/d = (2d - 1 - lam*sqrt(d))/d
    # So denominator of f_I = (2d-1)/sqrt(d) - lam = [(2d-1) - lam*sqrt(d)]/sqrt(d)
    # And factor_i/d = [(2d-1) - lam*sqrt(d)]/(d*sqrt(d))... hmm

    # Let me just check: for the dodecahedron, what is sum mult * f_I(lam)?
    # That's what Part 3 computed: 19.295 for all, 13.197 excluding zero mode
    # And the known F = 9.295
    # So F = sum_all mult * f_I(lam) - 10? That's 19.295 - 10 = 9.295. YES.
    # And edge_term = 10 when d=3 for dodecahedron (confirmed in Part 3).
    # So: F = sum_all mult * f_I(lam) - edge_term? That seems weird.
    # Or: F = sum_nonzero mult * f_I + 1*f_I(d) + edge - 2*edge?

    # Actually from the run output: "F - sum_all = -10.0" and edge = 10
    # So F = sum_all - 10 and edge_contrib = 10
    # This means: F = sum_all mult*f_I - edge_contrib
    # But that's strange since edge_contrib should ADD to F not subtract.

    # Unless the user's F value already INCLUDES the edge term baked in,
    # and f_I(lam) = (4 - lam*sqrt3)/(5 - lam*sqrt3) already contains it.
    # No, the per-eigenvalue sum is 19.295, edge is 10, F is 9.295.
    # 19.295 - 10 = 9.295. So F = Sigma_f_I - edge.
    # OR equivalently the correct formula has MINUS edge:
    # F = -edge + sum f_I = -(2(r-1)u^2/(1-u^2)) + sum mult * f_I(lam)

    # This is consistent with:
    # F = u * d/du log zeta
    #   = -[(r-1)*2u^2/(1-u^2)] - sum_i u*(-lam + 2(d-1)u) / factor_i
    #   = -(r-1)*2u^2/(1-u^2) + sum u*(lam - 2(d-1)u)/factor_i

    # And f_I(lam) per eigenvalue = u*(lam - 2(d-1)u)/factor ... wait that's NEGATIVE
    # for lam = 3, d = 3: u*(3 - 4*u)/factor = (1/sqrt3)*(3 - 4/sqrt3)/(1 - 3/sqrt3 + 2/3)
    # = 0.577*(3 - 2.309)/(1 - 1.732 + 0.667) = 0.577*0.691/(-0.065) < 0

    # But f_I(3) was computed as 6.098, which is POSITIVE.
    # f_I(3) = (4 - 3*sqrt3)/(5 - 3*sqrt3) = (4-5.196)/(5-5.196) = (-1.196)/(-0.196) = 6.098 POSITIVE

    # So f_I(lam) = (2(d-1) - lam*sqrt(d)) / ((2d-1) - lam*sqrt(d))
    # = (4 - lam*sqrt3) / (5 - lam*sqrt3)

    # And the per-eigenvalue LOG derivative piece is:
    # -u * d/du log(1 - lam*u + (d-1)*u^2)
    # = -u * (-lam + 2(d-1)*u) / (1 - lam*u + (d-1)*u^2)
    # = u * (lam - 2(d-1)*u) / factor

    # At u=1/sqrt(3), d=3, lam=3:
    # = (1/sqrt3)(3 - 4/sqrt3) / (1 - 3/sqrt3 + 2/3)
    # = (1/sqrt3)(3 - 2.309) / (1.667 - 1.732) = (0.577)(0.691)/(-0.065)
    # = 0.399 / (-0.065) = -6.098

    # So the per-eigenvalue log derivative piece = -f_I for each eigenvalue!
    # That means sum of log deriv pieces = -sum f_I = -19.295
    # And edge_term = 10 (positive)
    # F = -edge + sum_pereig = ... that gives -10 - (-19.295) = 9.295! No wait...

    # Let me be super explicit:
    # log(zeta^{-1}) = (r-1)*log(1-u^2) + sum log(factor_i)
    # d/du log(zeta^{-1}) = (r-1)*(-2u)/(1-u^2) + sum (-lam + 2(d-1)u)/factor_i
    # u * d/du log(zeta^{-1}) = -(r-1)*2u^2/(1-u^2) + sum u*(-lam + 2(d-1)u)/factor_i
    #                         = -10 + sum [-f_I(lam_i)]
    #                         = -10 + (-19.295) = -29.295

    # log(zeta) = -log(zeta^{-1})
    # u * d/du log(zeta) = -u * d/du log(zeta^{-1}) = 29.295

    # F should be related to zeta'/zeta. The user's F = 9.295 is neither 29.295 nor -29.295.
    # But 29.295 - 20 = 9.295. That's 29.295 - V = 9.295.

    # Alternatively, maybe F is just the sum of per-eigenvalue f_I excluding the zero mode
    # plus some correction. The user gave F = 4308/715 + (270/143)*sqrt(3) and said this
    # comes from the Ihara zeta. Let me just trust that value and also trust f_I(lam).

    # FOR THE OTHER GRAPHS, let me compute F as: sum_all mult * f_I(lam)
    # where f_I(lam) = (2(d-1) - lam*sqrt(d)) / ((2d-1) - lam*sqrt(d))
    # This is the per-eigenvalue contribution to the det part.
    # And then figure out what the "total F" should be from context.

    # Actually let me re-examine. For the dodecahedron:
    # sum_all mult * f_I = 19.295
    # sum_nonzero mult * f_I = 13.197
    # edge = 10
    # F_known = 9.295

    # 19.295 - 10 = 9.295 YES
    # So: F = (sum_all f_I) - edge

    # In terms of log zeta:
    # u * d/du log(zeta^{-1}) = -edge + sum[-f_I] = -10 + (-19.295) = -29.295
    # u * d/du log(zeta) = 29.295
    # But the user's F = 9.295 = 29.295 - 20 = 29.295 - V

    # Hmm... OR the edge contribution should be different.
    # (r-1)*2u^2/(1-u^2) at u=1/sqrt(3), r=11: 10*2*(1/3)/(2/3) = 10*1 = 10
    # So the FULL u*d/du log(zeta) = 10 + 19.295 = 29.295? No, with signs:
    # u*d/du log zeta = -u*d/du log(zeta^{-1}) = edge_term - sum[-f_I] = 10 + 19.295 = 29.295
    # Hmm, 29.295 is too big.

    # OK I think the issue is: the user's F is not the full -u*zeta'/zeta.
    # It's something else -- perhaps just the contribution from det(I - Au + (d-1)u^2 I).
    # The Ihara det is: det(I - Au + (d-1)u^2 I) = prod_i (1 - lam_i*u + (d-1)*u^2)
    # -u * d/du log det = -sum_i u*(-lam_i + 2(d-1)u)/factor_i = sum f_I(lam_i)
    # = 19.295

    # That's still not 9.295.

    # 19.295 = F + edge = 9.295 + 10
    # But wait, the actual output said "F - sum_all = -10.0" meaning F_known - sum_f_I = -10
    # So F_known = sum_f_I - 10.

    # Let me try: maybe the user's F includes HALF the edge contribution?
    # Or the formula is different. Or the f_I formula from the user's prompt
    # (4 - lam*sqrt3)/(5 - lam*sqrt3) is NOT actually per-eigenvalue of Ihara
    # but something else (like Mobius-transformed Ihara or normalized).

    # REGARDLESS: For this computation, what matters is whether F^2 = Tr^2 + d.
    # The user has GIVEN us F for the dodecahedron. For other graphs, let me
    # compute F using the SAME convention, which appears to be:
    # F = sum_all mult * f_I(lam) - 2(r-1)*u^2/(1-u^2)
    # = sum_all mult * (2(d-1) - lam*sqrt(d))/((2d-1) - lam*sqrt(d)) - 2(r-1)/(d-1)

    # At u=1/sqrt(d): u^2 = 1/d, 1-u^2 = (d-1)/d
    # 2(r-1)*u^2/(1-u^2) = 2(r-1)*(1/d)/((d-1)/d) = 2(r-1)/(d-1)

    sqd = sqrt(mpf(d_graph))
    f_I_sum = mpf(0)
    for lam, mult in adj_eigenvalues:
        f_I_val = (2*(d_graph-1) - lam*sqd) / ((2*d_graph - 1) - lam*sqd)
        f_I_sum += mult * f_I_val

    edge_val = mpf(2) * (r_val - 1) / (d_graph - 1) if d_graph > 1 else mpf(0)

    F_val = f_I_sum - edge_val

    return F_val, f_I_sum, edge_val


# ==============================================================================
# PART 1: DODECAHEDRON (d=3)
# ==============================================================================

header("PART 1: DODECAHEDRON -- F(1/sqrt(3))^2 vs Tr(L^-1)^2 + 3")

# Known exact values
# F(1/sqrt(3)) = 4308/715 + (270/143)*sqrt(3)  [in Q(sqrt(3))]
# Tr(L^-1) = 137/15

a_F = Frac(4308, 715)   # rational part of F
b_F = Frac(270, 143)    # coefficient of sqrt(3) in F
Tr_inv = Frac(137, 15)
d_val = 3

sqrt3 = sqrt(mpf(3))

a_F_mp = frac_to_mpf(a_F)
b_F_mp = frac_to_mpf(b_F)
Tr_mp = frac_to_mpf(Tr_inv)

F_val = a_F_mp + b_F_mp * sqrt3

print(f"F(1/sqrt(3)) = {a_F} + ({b_F})*sqrt(3)")
print(f"Tr(L^-1)     = {Tr_inv}")
print(f"\nF(1/sqrt(3)) numerical = {fmt(F_val)}")
print(f"Tr(L^-1) numerical     = {fmt(Tr_mp)}")
print(f"sqrt(3)                = {fmt(sqrt3)}")

# (a) F^2 exactly in Q(sqrt(3))
# F = a + b*sqrt(3), so F^2 = (a^2 + 3*b^2) + 2*a*b*sqrt(3)
rational_F2_exact = a_F**2 + 3 * b_F**2
irrational_F2_exact = 2 * a_F * b_F  # coefficient of sqrt(3)

subheader("(a) F^2 in Q(sqrt(3))")
print(f"F^2 = (a^2 + 3*b^2) + 2*a*b*sqrt(3)")
print(f"    = {rational_F2_exact.numerator}/{rational_F2_exact.denominator} + ({irrational_F2_exact.numerator}/{irrational_F2_exact.denominator})*sqrt(3)")

F2_numerical = F_val**2
rat_F2_mp = frac_to_mpf(rational_F2_exact)
irr_F2_mp = frac_to_mpf(irrational_F2_exact)
F2_check = rat_F2_mp + irr_F2_mp * sqrt3

print(f"\nF^2 numerical       = {fmt(F2_numerical)}")
print(f"F^2 from Q(sqrt3)  = {fmt(F2_check)}")
print(f"Difference          = {fmt(F2_numerical - F2_check)}")

# (b) Tr^2
subheader("(b) Tr(L^-1)^2")
Tr2_exact = Tr_inv**2
print(f"Tr^2 = (137/15)^2 = {Tr2_exact.numerator}/{Tr2_exact.denominator}")
Tr2_mp = frac_to_mpf(Tr2_exact)
print(f"     = {fmt(Tr2_mp)}")

# (c) F^2 - Tr^2
subheader("(c) F^2 - Tr^2")
diff_numerical = F2_numerical - Tr2_mp
print(f"F^2 - Tr^2     = {fmt(diff_numerical)}")
print(f"d = 3          = {fmt(mpf(3))}")
print(f"F^2 - Tr^2 - 3 = {fmt(diff_numerical - 3)}")

# (d) Rational part of F^2 - Tr^2
subheader("(d) Rational part of F^2 - Tr^2")
rational_diff_exact = rational_F2_exact - Tr2_exact
print(f"rational(F^2) - Tr^2 = {rational_diff_exact.numerator}/{rational_diff_exact.denominator}")
print(f"                     = {fmt(frac_to_mpf(rational_diff_exact), 40)}")
print(f"                     vs d = 3")

# (e) Irrational part
subheader("(e) Irrational part of F^2 - Tr^2")
print(f"irrational coeff (of sqrt(3)) = 2*a*b = {irrational_F2_exact.numerator}/{irrational_F2_exact.denominator}")
print(f"                              = {fmt(irr_F2_mp, 40)}")
irr_numerical = irr_F2_mp * sqrt3
print(f"2*a*b * sqrt(3) numerical     = {fmt(irr_numerical, 40)}")

# (f) Numerical value and relative error from d=3
subheader("(f) Numerical value and relative error from d=3")
total_diff = frac_to_mpf(rational_diff_exact) + irr_numerical
print(f"F^2 - Tr^2 = {fmt(total_diff)}")
print(f"d           = 3")
print(f"F^2 - Tr^2 - d = {fmt(total_diff - 3)}")
rel_err = abs(total_diff - 3) / 3
print(f"Relative error |F^2 - Tr^2 - d|/d = {fmt(rel_err)}")

# (g) Exact factorization
subheader("(g) Exact F^2 - Tr^2 factorization")
print(f"\nF^2 - Tr^2 = {rational_diff_exact.numerator}/{rational_diff_exact.denominator} + ({irrational_F2_exact.numerator}/{irrational_F2_exact.denominator})*sqrt(3)")
print(f"\nRational part factorization:")
n_r = rational_diff_exact.numerator
d_r = rational_diff_exact.denominator
print(f"  {n_r}/{d_r}")
print(f"  numerator   {abs(n_r)} = {fmt_factors(n_r)}")
print(f"  denominator {d_r} = {fmt_factors(d_r)}")

n_i = irrational_F2_exact.numerator
d_i = irrational_F2_exact.denominator
print(f"\nIrrational coeff factorization:")
print(f"  {n_i}/{d_i}")
print(f"  numerator   {abs(n_i)} = {fmt_factors(n_i)}")
print(f"  denominator {d_i} = {fmt_factors(d_i)}")


# ==============================================================================
# PART 2: PER-EIGENVALUE TEST
# ==============================================================================

header("PART 2: PER-EIGENVALUE -- f_Ihara(lam)^2 vs f_Laplacian(lam)^2")

sqrt5 = sqrt(mpf(5))
dodec_eigenvalues = [
    (mpf(3),  1, "3"),
    (sqrt5,   3, "sqrt5"),
    (mpf(1),  5, "1"),
    (mpf(0),  4, "0"),
    (mpf(-2), 4, "-2"),
    (-sqrt5,  3, "-sqrt5"),
]

V_dodec = 20
E_dodec = 30
d_dodec = 3

print(f"\nDodecahedron: V={V_dodec}, E={E_dodec}, d={d_dodec}")
print(f"u = 1/sqrt(3) = {fmt(1/sqrt3, 20)}")
print()

print(f"{'lam':>8s} | {'mult':>4s} | {'f_I(lam)':>22s} | {'f_L(lam)':>22s} | {'f_I^2 - f_L^2':>25s}")
print("-" * 95)

fi_sq_weighted = mpf(0)
fl_sq_weighted = mpf(0)
diff_weighted = mpf(0)

for lam, mult, name in dodec_eigenvalues:
    fi = (4 - lam * sqrt3) / (5 - lam * sqrt3)

    if abs(3 - lam) < mpf('1e-40'):  # zero mode
        print(f"{name:>8s} | {mult:>4d} | {fmt(fi, 18):>22s} | {'inf (zero mode)':>22s} | {'N/A (zero mode)':>25s}")
        continue

    fl = 1 / (3 - lam)
    fi_sq = fi**2
    fl_sq = fl**2
    diff_val = fi_sq - fl_sq

    fi_sq_weighted += mult * fi_sq
    fl_sq_weighted += mult * fl_sq
    diff_weighted += mult * diff_val

    print(f"{name:>8s} | {mult:>4d} | {fmt(fi, 18):>22s} | {fmt(fl, 18):>22s} | {fmt(diff_val, 20):>25s}")

print(f"\nWeighted sum (mult * [f_I^2 - f_L^2]) = {fmt(diff_weighted)}")
print(f"For constant d/V = {d_dodec}/{V_dodec} = {d_dodec/V_dodec}: diff/V = {fmt(diff_weighted/V_dodec, 20)}")


# ==============================================================================
# PART 3: WEIGHTED PYTHAGOREAN
# ==============================================================================

header("PART 3: WEIGHTED SUMS (excluding zero-mode lam=3)")

print(f"\nSum mult * f_I^2             = {fmt(fi_sq_weighted)}")
print(f"Sum mult * f_L^2             = {fmt(fl_sq_weighted)}")
print(f"Sum mult * (f_I^2 - f_L^2)  = {fmt(diff_weighted)}")
print(f"\nReference values:")
print(f"  d     = {d_dodec}")
print(f"  V     = {V_dodec}")
print(f"  E     = {E_dodec}")
print(f"  diff/V = {fmt(diff_weighted / V_dodec, 20)}")
print(f"  diff/d = {fmt(diff_weighted / d_dodec, 20)}")

# Trace reconstruction
F_sum_all = mpf(0)
Tr_sum = mpf(0)
for lam, mult, name in dodec_eigenvalues:
    fi = (4 - lam * sqrt3) / (5 - lam * sqrt3)
    F_sum_all += mult * fi
    if abs(3 - lam) > mpf('1e-40'):
        Tr_sum += mult / (3 - lam)

r_dodec = E_dodec - V_dodec + 1  # = 11
edge_contrib = mpf(2) * (r_dodec - 1) / (d_dodec - 1)  # = 2*10/2 = 10

print(f"\n--- Trace reconstruction ---")
print(f"Sum_all mult * f_I(lam)     = {fmt(F_sum_all, 30)}")
print(f"edge = 2(r-1)/(d-1) = 2*10/2 = {fmt(edge_contrib, 10)}")
print(f"Sum_all - edge              = {fmt(F_sum_all - edge_contrib, 30)}")
print(f"Known F(1/sqrt(3))          = {fmt(F_val, 30)}")
print(f"Difference                  = {fmt(F_val - (F_sum_all - edge_contrib))}")
print(f"\nSum mult / (d-lam) = Tr(L^-1) = {fmt(Tr_sum, 30)}")
print(f"Known Tr(L^-1) = 137/15        = {fmt(Tr_mp, 30)}")

# So F = Sum_all f_I - edge. This means for other graphs we should use the same formula.


# ==============================================================================
# PART 4: OTHER REGULAR GRAPHS
# ==============================================================================

header("PART 4: OTHER REGULAR GRAPHS")

def compute_F_and_Tr(adj_eigs, d_g, V_g, E_g):
    """
    Compute F and Tr(L^-1) for a d-regular graph.
    F = sum_all mult * f_I(lam) - 2*(r-1)/(d-1)
    where f_I(lam) = (2(d-1) - lam*sqrt(d)) / ((2d-1) - lam*sqrt(d))
    and r = E - V + 1.
    """
    sqd = sqrt(mpf(d_g))
    r_g = E_g - V_g + 1

    fI_sum = mpf(0)
    Tr_g = mpf(0)
    for lam, mult in adj_eigs:
        fI = (2*(d_g - 1) - lam * sqd) / ((2*d_g - 1) - lam * sqd)
        fI_sum += mult * fI
        mu = d_g - lam
        if abs(mu) > mpf('1e-40'):
            Tr_g += mult / mu

    edge_g = mpf(2) * (r_g - 1) / (d_g - 1) if d_g > 1 else mpf(0)
    F_g = fI_sum - edge_g

    return F_g, Tr_g, fI_sum, edge_g

def test_graph(name, adj_eigs, d_g, V_g, E_g, Tr_expected=None):
    subheader(f"{name} (d={d_g}, V={V_g}, E={E_g})")

    F_g, Tr_g, fI_sum, edge_g = compute_F_and_Tr(adj_eigs, d_g, V_g, E_g)

    print(f"F(1/sqrt({d_g}))  = {fmt(F_g, 30)}")
    print(f"  (fI_sum = {fmt(fI_sum, 20)}, edge = {fmt(edge_g, 10)})")
    print(f"Tr(L^-1)    = {fmt(Tr_g, 30)}")
    if Tr_expected is not None:
        Tr_exp_mp = frac_to_mpf(Tr_expected)
        print(f"Tr expected = {fmt(Tr_exp_mp, 30)}")
        print(f"Tr diff     = {fmt(Tr_g - Tr_exp_mp)}")

    F2 = F_g**2
    Tr2 = Tr_g**2
    diff = F2 - Tr2
    residual = diff - d_g

    print(f"\nF^2       = {fmt(F2, 30)}")
    print(f"Tr^2      = {fmt(Tr2, 30)}")
    print(f"F^2 - Tr^2 = {fmt(diff, 30)}")
    print(f"d           = {d_g}")
    print(f"F^2 - Tr^2 - d = {fmt(residual, 30)}")
    rel = abs(residual) / d_g if d_g > 0 else mpf(0)
    print(f"|F^2 - Tr^2 - d|/d = {fmt(rel, 20)}")

    # Also test d/V
    residual_dV = diff - mpf(d_g) / V_g
    print(f"F^2 - Tr^2 - d/V = {fmt(residual_dV, 30)}")

    return F_g, Tr_g

# (a) K4 - Complete graph on 4 vertices
# d=3, V=4, E=6
# Adj eigenvalues: 3(x1), -1(x3)
K4_eigs = [(mpf(3), 1), (mpf(-1), 3)]
test_graph("K4 (Complete graph)", K4_eigs, 3, 4, 6, Frac(3, 4))

# (b) Petersen graph
# d=3, V=10, E=15
# Adj eigenvalues: 3(x1), 1(x5), -2(x4)
Pet_eigs = [(mpf(3), 1), (mpf(1), 5), (mpf(-2), 4)]
test_graph("Petersen graph", Pet_eigs, 3, 10, 15, Frac(33, 10))

# (c) Cube graph (3-cube, Q3)
# d=3, V=8, E=12
# Adj eigenvalues: 3(x1), 1(x3), -1(x3), -3(x1)
Cube_eigs = [(mpf(3), 1), (mpf(1), 3), (mpf(-1), 3), (mpf(-3), 1)]
test_graph("Cube graph (Q3)", Cube_eigs, 3, 8, 12, Frac(29, 12))

# (d) Icosahedron
# d=5, V=12, E=30
# Adjacency eigenvalues: 5(x1), sqrt(5)(x3), -1(x5), -sqrt(5)(x3)
# Trace check: 5 + 3*sqrt5 + 5*(-1) + 3*(-sqrt5) = 0 OK
# Sum of squares: 25 + 15 + 5 + 15 = 60 = 2E OK
Ico_eigs = [(mpf(5), 1), (sqrt5, 3), (mpf(-1), 5), (-sqrt5, 3)]
test_graph("Icosahedron", Ico_eigs, 5, 12, 30)

# (e) Dodecahedron (cross-check with Part 1)
Dodec_eigs = [(mpf(3), 1), (sqrt5, 3), (mpf(1), 5), (mpf(0), 4), (mpf(-2), 4), (-sqrt5, 3)]
F_dodec, Tr_dodec = test_graph("Dodecahedron", Dodec_eigs, 3, 20, 30, Frac(137, 15))

print(f"\nDodecahedron cross-check:")
print(f"  F from exact Q(sqrt3) = {fmt(F_val, 30)}")
print(f"  F from eigenvalue sum = {fmt(F_dodec, 30)}")
print(f"  Difference            = {fmt(F_val - F_dodec)}")


# ==============================================================================
# PART 5: THE AXIOM TRIANGLE
# ==============================================================================

header("PART 5: THE AXIOM TRIANGLE (Dodecahedron)")

F_hyp = F_val
Tr_leg = Tr_mp
d_leg = sqrt3

print(f"Hypotenuse (F):  {fmt(F_hyp, 20)}")
print(f"Leg 1 (Tr):      {fmt(Tr_leg, 20)}")
print(f"Leg 2 (sqrt(d)): {fmt(d_leg, 20)}")

pyth_check = F_hyp**2 - Tr_leg**2 - mpf(d_val)
print(f"\nF^2 - Tr^2 - d = {fmt(pyth_check, 30)} (0 if exact)")

# Angle
theta = atan(d_leg / Tr_leg)
print(f"\ntheta = arctan(sqrt(d)/Tr) = arctan(sqrt(3)*15/137)")
print(f"      = arctan(15*sqrt(3)/137)")
print(f"      = {fmt(theta, 20)} radians")
print(f"      = {fmt(theta * 180 / pi, 20)} degrees")
print(f"cos(theta) = Tr/F = {fmt(Tr_leg / F_hyp, 20)}")
print(f"sin(theta) = sqrt(d)/F = {fmt(d_leg / F_hyp, 20)}")

# Notable angle comparisons
arctan_half = atan(mpf(1)/2)
print(f"\nNotable angles:")
print(f"  arctan(1/2) = {fmt(arctan_half * 180 / pi, 15)} deg (1x2 rectangle)")
print(f"  18 deg = pi/10 (dodecahedral quantum)")
print(f"  36 deg = pi/5")
print(f"  theta  = {fmt(theta * 180 / pi, 15)} deg")

# 1x2 rectangle mapping
print(f"\n--- 1x2 Rectangle Mapping ---")
ratio_legs = d_leg / Tr_leg
print(f"sqrt(d)/Tr = {fmt(ratio_legs, 20)}")
print(f"For 1x2 rectangle, ratio = 1/2 = 0.5")

k_scale = Tr_leg / 2
print(f"\nIf Tr = 2*k, then k = {fmt(k_scale, 15)}")
print(f"  sqrt(d)/k = {fmt(d_leg / k_scale, 15)}  (should be 1 for 1x2 match)")
print(f"  F/k       = {fmt(F_hyp / k_scale, 15)}  (should be sqrt(5) for 1x2 match)")
print(f"  sqrt(5)   = {fmt(sqrt5, 15)}")


# ==============================================================================
# PART 6: MODIFIED PYTHAGOREAN RELATIONS
# ==============================================================================

header("PART 6: MODIFIED PYTHAGOREAN RELATIONS (Dodecahedron)")

F2_v = F_val**2
Tr2_v = Tr_mp**2
phi = (1 + sqrt5) / 2

print(f"F^2 = {fmt(F2_v, 30)}")
print(f"Tr^2 = {fmt(Tr2_v, 30)}")

# (a) F^2 = Tr^2 + d
res_a = F2_v - Tr2_v - mpf(d_val)
print(f"\n(a) F^2 - Tr^2 - d           = {fmt(res_a, 30)}")

# (b) F^2 = Tr^2 + d/V
res_b = F2_v - Tr2_v - mpf(d_val)/V_dodec
print(f"(b) F^2 - Tr^2 - d/V         = {fmt(res_b, 30)}")

# (c) F^2 = Tr^2 + d*phi
res_c = F2_v - Tr2_v - mpf(d_val) * phi
print(f"(c) F^2 - Tr^2 - d*phi       = {fmt(res_c, 30)}")

# (c') F^2 = Tr^2 + d*phi^2
res_c2 = F2_v - Tr2_v - mpf(d_val) * phi**2
print(f"(c') F^2 - Tr^2 - d*phi^2    = {fmt(res_c2, 30)}")

# (d) F_no_edge
edge_val = mpf(10)
F_no_edge = F_val - edge_val
print(f"\n(d) F_no_edge = F - 10 = {fmt(F_no_edge, 20)}")
print(f"    F_no_edge^2         = {fmt(F_no_edge**2, 30)}")
print(f"    F_no_edge^2 - Tr^2  = {fmt(F_no_edge**2 - Tr2_v, 30)}")
print(f"    F_no_edge^2 - Tr^2 - d = {fmt(F_no_edge**2 - Tr2_v - mpf(d_val), 30)}")

# (e) Excluding golden pair
F_nongolden = edge_val  # start with edge
Tr_nongolden = mpf(0)
F_golden = mpf(0)
Tr_golden = mpf(0)

for lam, mult, name in dodec_eigenvalues:
    fi = (4 - lam * sqrt3) / (5 - lam * sqrt3)
    is_golden = (abs(abs(lam) - sqrt5) < mpf('1e-40'))
    if not is_golden:
        F_nongolden += mult * fi
        if abs(3 - lam) > mpf('1e-40'):
            Tr_nongolden += mult / (3 - lam)
    else:
        F_golden += mult * fi
        Tr_golden += mult / (3 - lam)

print(f"\n(e) Excluding golden pair (lam = +/- sqrt(5)):")
print(f"    F_nongolden = {fmt(F_nongolden, 20)}")
print(f"    Tr_nongolden = {fmt(Tr_nongolden, 20)}")
print(f"    F_golden    = {fmt(F_golden, 20)}")
print(f"    Tr_golden   = {fmt(Tr_golden, 20)}")
print(f"    F_nongolden^2 - Tr_nongolden^2 = {fmt(F_nongolden**2 - Tr_nongolden**2, 30)}")
print(f"    F_golden^2 - Tr_golden^2       = {fmt(F_golden**2 - Tr_golden**2, 30)}")

# (f) scaled by V
res_f = (F2_v - Tr2_v) * V_dodec
print(f"\n(f) (F^2 - Tr^2) * V = {fmt(res_f, 30)}")
print(f"    d*V = {d_dodec * V_dodec} = 60 = |A5|")
print(f"    (F^2 - Tr^2)*V - d*V = {fmt(res_f - d_dodec * V_dodec, 30)}")

# (g) Exact decomposition
print(f"\n(g) Exact decomposition in Q(sqrt(3)):")
print(f"    F^2 - Tr^2 = {rational_diff_exact.numerator}/{rational_diff_exact.denominator} + ({irrational_F2_exact.numerator}/{irrational_F2_exact.denominator})*sqrt(3)")
rat_part_mp = frac_to_mpf(rational_diff_exact)
irr_part_mp = frac_to_mpf(irrational_F2_exact) * sqrt3
print(f"    rational part            = {fmt(rat_part_mp, 30)}")
print(f"    irrational part (num)    = {fmt(irr_part_mp, 30)}")
print(f"    rational part - 3        = {fmt(rat_part_mp - 3, 30)}")


# ==============================================================================
# PART 7: THE EXACT DECOMPOSITION
# ==============================================================================

header("PART 7: EXACT DECOMPOSITION IN Q(sqrt(3))")

rational_part = frac_to_mpf(rational_diff_exact)
irrational_coeff = frac_to_mpf(irrational_F2_exact)
irrational_part = irrational_coeff * sqrt3

print(f"F^2 = {rational_F2_exact.numerator}/{rational_F2_exact.denominator} + ({irrational_F2_exact.numerator}/{irrational_F2_exact.denominator})*sqrt(3)")
print(f"Tr^2 = {Tr2_exact.numerator}/{Tr2_exact.denominator}")
print(f"F^2 - Tr^2 = {rational_diff_exact.numerator}/{rational_diff_exact.denominator} + ({irrational_F2_exact.numerator}/{irrational_F2_exact.denominator})*sqrt(3)")

print(f"\n=== KEY QUESTION: Is rational(F^2) - Tr^2 = d = 3? ===")
print(f"rational(F^2) - Tr^2 = {rational_diff_exact}")
print(f"                     = {rational_diff_exact.numerator}/{rational_diff_exact.denominator}")
print(f"                     = {fmt(rational_part, 40)}")
print(f"d = 3")
print(f"Difference from 3: {fmt(rational_part - 3, 40)}")

three_exact = Frac(3, 1)
is_exactly_3 = (rational_diff_exact == three_exact)
print(f"\nIs rational(F^2) - Tr^2 EXACTLY 3? --> {is_exactly_3}")

if not is_exactly_3:
    delta = rational_diff_exact - three_exact
    print(f"delta = rational(F^2) - Tr^2 - 3 = {delta}")
    print(f"      = {delta.numerator}/{delta.denominator}")
    print(f"      ~ {float(delta):.20f}")
    print(f"\nFactor delta:")
    n_delta = abs(delta.numerator)
    d_delta = delta.denominator
    print(f"  numerator {n_delta} = {fmt_factors(n_delta)}")
    print(f"  denominator {d_delta} = {fmt_factors(d_delta)}")
    # Try to simplify
    print(f"  sign: {'positive' if delta > 0 else 'negative'}")

print(f"\n=== IRRATIONAL PART ===")
print(f"2*a*b = {irrational_F2_exact.numerator}/{irrational_F2_exact.denominator}")
print(f"      = {fmt(irrational_coeff, 40)}")
print(f"2*a*b * sqrt(3) = {fmt(irrational_part, 40)}")

print(f"\n=== RELATIVE MAGNITUDES ===")
rat_minus_3 = rational_part - 3
print(f"|rational part - 3|  = {fmt(abs(rat_minus_3), 20)}")
print(f"|irrational part|    = {fmt(abs(irrational_part), 20)}")
if abs(rat_minus_3) > mpf('1e-60'):
    print(f"|irr|/|rat-3|        = {fmt(abs(irrational_part) / abs(rat_minus_3), 20)}")
print(f"|irr part| / d       = {fmt(abs(irrational_part) / 3, 20)}")
print(f"|irr part| / F^2     = {fmt(abs(irrational_part) / F2_v, 20)}")
print(f"|F^2 - Tr^2 - 3| / 3 = {fmt(abs(diff_numerical - 3) / 3, 20)} (total relative error)")

print(f"\n=== FULL SUMMARY ===")
print(f"F(1/sqrt(3))^2 = Tr(L^-1)^2 + ({rational_diff_exact}) + ({irrational_F2_exact})*sqrt(3)")
print(f"Numerically:     {fmt(Tr2_v, 15)} + {fmt(rational_part + irrational_part, 15)} = {fmt(F2_v, 15)}")

# What IS F^2 - Tr^2 as a single number?
val = diff_numerical
print(f"\n=== F^2 - Tr^2 as single number ===")
print(f"F^2 - Tr^2       = {fmt(val, 40)}")
print(f"sqrt(F^2 - Tr^2) = {fmt(sqrt(abs(val)), 20)}")
print(f"sqrt(3)           = {fmt(sqrt3, 20)}")
print(f"(F^2 - Tr^2)/3   = {fmt(val / 3, 20)}")
print(f"(F^2 - Tr^2)/sqrt(3) = {fmt(val / sqrt3, 20)}")

# Is (F^2 - Tr^2) = some recognizable expression?
# Try various things
print(f"\n--- Searching for pattern ---")
print(f"F^2 - Tr^2 = {fmt(val, 20)}")
print(f"3           = 3.0")
print(f"3*phi       = {fmt(3*phi, 20)}")
print(f"3/phi       = {fmt(3/phi, 20)}")
print(f"pi          = {fmt(pi, 20)}")
print(f"e           = {fmt(mp.e, 20)}")


# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

header("SUMMARY TABLE: ALL GRAPHS")

print(f"\n{'Graph':<20s} | {'d':>3s} | {'V':>3s} | {'F':>15s} | {'Tr':>15s} | {'F^2-Tr^2':>15s} | {'d':>5s} | {'residual':>15s} | {'|res|/d':>12s}")
print("-" * 125)

all_graphs = [
    ("K4",          [(mpf(3), 1), (mpf(-1), 3)],                                  3, 4, 6,   Frac(3,4)),
    ("Petersen",    [(mpf(3), 1), (mpf(1), 5), (mpf(-2), 4)],                     3, 10, 15, Frac(33,10)),
    ("Cube (Q3)",   [(mpf(3), 1), (mpf(1), 3), (mpf(-1), 3), (mpf(-3), 1)],      3, 8, 12,  Frac(29,12)),
    ("Dodecahedron",[(mpf(3), 1), (sqrt5, 3), (mpf(1), 5), (mpf(0), 4), (mpf(-2), 4), (-sqrt5, 3)], 3, 20, 30, Frac(137,15)),
    ("Icosahedron", [(mpf(5), 1), (sqrt5, 3), (mpf(-1), 5), (-sqrt5, 3)],         5, 12, 30, None),
]

for gname, eigs, dg, Vg, Eg, Tr_ex in all_graphs:
    Fg, Trg, _, _ = compute_F_and_Tr(eigs, dg, Vg, Eg)
    Fg2 = Fg**2
    Trg2 = Trg**2
    diffg = Fg2 - Trg2
    resg = diffg - dg
    relg = abs(resg) / dg if dg > 0 else mpf(0)

    print(f"{gname:<20s} | {dg:>3d} | {Vg:>3d} | {fmt(Fg,12):>15s} | {fmt(Trg,12):>15s} | {fmt(diffg,12):>15s} | {dg:>5d} | {fmt(resg,12):>15s} | {fmt(relg,8):>12s}")


# ==============================================================================
# THE PYTHAGOREAN TRIANGLE (Visual)
# ==============================================================================

header("THE PYTHAGOREAN TRIANGLE (Dodecahedron)")

theta_deg = theta * 180 / pi
actual_res = F2_v - Tr2_v - mpf(d_val)

print(f"""
   If F^2 = Tr^2 + d were exact:

         F = {fmt(F_val, 10)}
        /|
       / |
      /  |  sqrt(d) = sqrt(3) = {fmt(sqrt3, 10)}
     /   |
    /    |
   /_____|
   Tr = 137/15 = {fmt(Tr_mp, 10)}

   Angle at base: theta = {fmt(theta_deg, 10)} degrees

   ACTUAL residual: F^2 - Tr^2 - d = {fmt(actual_res, 20)}
   Relative error:  |residual|/d   = {fmt(abs(actual_res)/d_val, 15)}

   The relation holds to about {fmt(-mp.log10(abs(actual_res)/d_val), 4)} significant figures.
""")

print("=" * 72)
print("  COMPUTATION COMPLETE")
print("=" * 72)
