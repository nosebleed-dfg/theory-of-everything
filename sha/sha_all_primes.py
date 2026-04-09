#!/usr/bin/env python3
"""
SHA_ALL_PRIMES — tests all 8 SHA-256 prime polynomials against the linearized round char_poly
nos3bl33d

Offset-shift hypothesis, golden connections, GF(p) factorizations, axiom residual M^2-M-I.
"""

from sympy import (
    Symbol, Matrix, Poly, sqrt, floor, Rational, factor, gcd,
    GF, ZZ, QQ, roots, nroots, eye, pprint, simplify, expand,
    resultant, prod, div, rem, lcm, minimal_polynomial, nsimplify,
    AlgebraicNumber
)
from sympy.polys.polytools import Poly
from sympy.polys.domains import FiniteField
from sympy import factorint, isprime
import sys

x = Symbol('x')

# ============================================================
# SHA-256 Linearized Round Matrix
# ============================================================
M = Matrix([
    [1, 0, 0, 0, 1, 0, 0, 1],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 1],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0]
])

# SHA-256 primes
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]

# Golden ratio
from sympy import GoldenRatio as phi_sym
phi_val = (1 + 5**0.5) / 2

def separator(title):
    print(f"\n{'='*72}")
    print(f"  {title}")
    print(f"{'='*72}\n")


# ============================================================
# PART 0: Verify the characteristic polynomial
# ============================================================
separator("PART 0: CHARACTERISTIC POLYNOMIAL VERIFICATION")

char_poly_matrix = M.charpoly(x)
char_poly = Poly(char_poly_matrix.as_expr(), x, domain='ZZ')
print(f"Char poly of M:  {char_poly.as_expr()}")

# Expected: x^8 - 2x^7 + x^6 - x^4 - 1
expected = Poly(x**8 - 2*x**7 + x**6 - x**4 - 1, x, domain='ZZ')
print(f"Expected:        {expected.as_expr()}")
print(f"Match: {char_poly == expected}")

det_M = M.det()
print(f"\ndet(M) = {det_M}")
print(f"M is invertible: {det_M != 0}")

# Try to factor over Q
print(f"\nFactor over Q: {factor(char_poly.as_expr(), x)}")

# Check irreducibility over Q
try:
    is_irred = char_poly.is_irreducible
    print(f"Irreducible over Q: {is_irred}")
except:
    print("Irreducibility check: skipped")


# ============================================================
# PART 1: Check all 8 prime polynomials
# ============================================================
separator("PART 1: ALL 8 PRIME POLYNOMIALS vs CHAR POLY")

print("For each SHA-256 prime p, frac(sqrt(p)) satisfies:")
print("  q_p(x) = x^2 + 2k*x + (k^2 - p)  where k = floor(sqrt(p))")
print()

prime_polys = {}
remainders = {}

for p in PRIMES:
    k = int(p**0.5)  # floor(sqrt(p))
    frac_val = p**0.5 - k

    # q_p(x) = x^2 + 2kx + (k^2 - p)
    q_p = Poly(x**2 + 2*k*x + (k**2 - p), x, domain='ZZ')
    prime_polys[p] = q_p

    # Verify the root
    root_check = frac_val**2 + 2*k*frac_val + (k**2 - p)

    # Compute remainder of char_poly divided by q_p
    quotient, remainder = div(char_poly, q_p, x, domain='ZZ')
    remainders[p] = remainder

    print(f"p={p:2d}: k={k}, q_p = {q_p.as_expr()}")
    print(f"       frac(sqrt({p})) = {frac_val:.10f}")
    print(f"       root check: q_p(frac) = {root_check:.2e}")
    print(f"       char_poly mod q_p: remainder = {remainder.as_expr()}")

    if remainder == Poly(0, x, domain='ZZ') or remainder.is_zero:
        print(f"       *** DIVIDES! q_{p} is a factor of char_poly! ***")
    else:
        # Evaluate remainder at the actual root
        rem_expr = remainder.as_expr()
        rem_at_root = float(rem_expr.subs(x, frac_val))
        print(f"       remainder evaluated at root: {rem_at_root:.10f}")
    print()


# ============================================================
# PART 2: The offset shift — find delta such that
#         x^2 - x - 1 divides char_poly(x + delta)
# ============================================================
separator("PART 2: OFFSET-SHIFT HYPOTHESIS")

# Compute all 8 roots numerically
char_roots = nroots(char_poly.as_expr(), n=15)
print("All 8 roots of char_poly:")
for i, r in enumerate(char_roots):
    r_re = complex(r).real
    r_im = complex(r).imag
    is_real = abs(r_im) < 1e-12
    marker = " <-- REAL" if is_real else ""
    print(f"  r_{i}: {r_re:+.12f} {r_im:+.12f}i{marker}")

print(f"\nphi = {phi_val:.12f}")
print()

# For each real root, compute delta = root - phi
print("Distance from each REAL root to phi:")
real_roots = []
for i, r in enumerate(char_roots):
    r_c = complex(r)
    if abs(r_c.imag) < 1e-12:
        real_roots.append(r_c.real)
        delta = r_c.real - phi_val
        print(f"  r_{i} = {r_c.real:.12f},  delta = r - phi = {delta:+.15f}")

# The dominant eigenvalue
dominant = max(real_roots, key=abs)
delta_dominant = dominant - phi_val
print(f"\nDominant eigenvalue: {dominant:.12f}")
print(f"phi:                {phi_val:.12f}")
print(f"Delta (dominant - phi): {delta_dominant:.15f}")
print(f"Delta as fraction of phi: {delta_dominant/phi_val:.6f}")

# Check: does x^2 - x - 1 divide char_poly(x + delta) for specific deltas?
golden_poly = Poly(x**2 - x - 1, x, domain='ZZ')
print("\n--- Testing specific shifts ---")

test_deltas = [
    ("dominant - phi", delta_dominant),
    ("0.029", 0.029),
    ("0.032", 0.032),
    ("1/phi^3", 1/phi_val**3),
    ("sqrt(5)-2", 5**0.5 - 2),
    ("2-phi", 2 - phi_val),
]

for name, delta in test_deltas:
    shifted = char_poly.as_expr().subs(x, x + delta)
    shifted_at_phi = float(shifted.subs(x, phi_val))
    shifted_at_phi2 = float(shifted.subs(x, 1 - phi_val))  # other root of x^2-x-1
    print(f"  delta = {name} ({delta:.10f})")
    print(f"    char_poly(phi + delta) = {shifted_at_phi:.2e}")
    print(f"    char_poly((1-phi) + delta) = {shifted_at_phi2:.2e}")

# Systematic: find delta such that char_poly(phi + delta) = 0
# This means phi + delta is a root, so delta = root - phi
print("\n--- Exact shift values (root - phi) ---")
for i, r in enumerate(char_roots):
    r_c = complex(r)
    delta_c = r_c.real - phi_val + 1j*r_c.imag
    print(f"  delta_{i} = {delta_c.real:+.12f} {delta_c.imag:+.12f}i  (|delta| = {abs(delta_c):.10f})")


# ============================================================
# PART 3: ALL QUADRATIC FACTORS — remainder analysis
# ============================================================
separator("PART 3: QUADRATIC REMAINDER ANALYSIS")

print("Remainders char_poly mod q_p(x) = a*x + b:")
print(f"{'p':>4s} {'a':>20s} {'b':>20s} {'|remainder|':>15s}")
print("-" * 65)

for p in PRIMES:
    q_p = prime_polys[p]
    _, rem = div(char_poly, q_p, x, domain='ZZ')
    coeffs = rem.all_coeffs() if not rem.is_zero else [0, 0]
    if len(coeffs) == 1:
        a, b = 0, coeffs[0]
    elif len(coeffs) == 2:
        a, b = coeffs
    else:
        a, b = 0, 0

    k = int(p**0.5)
    frac_val = p**0.5 - k
    rem_val = float(a * frac_val + b)
    print(f"{p:4d} {str(a):>20s} {str(b):>20s} {abs(rem_val):15.6f}")

# Also check golden poly x^2 - x - 1
print("\nAlso checking x^2 - x - 1 (golden ratio polynomial):")
_, rem_golden = div(char_poly, golden_poly, x, domain='ZZ')
print(f"  char_poly mod (x^2-x-1) = {rem_golden.as_expr()}")
print(f"  Evaluated at phi: {float(rem_golden.as_expr().subs(x, phi_val)):.10f}")


# ============================================================
# PART 4: PRODUCT of all 8 prime polynomials
# ============================================================
separator("PART 4: PRODUCT OF ALL 8 PRIME POLYNOMIALS")

product_poly = Poly(1, x, domain='ZZ')
for p in PRIMES:
    product_poly = product_poly * prime_polys[p]

print(f"P(x) = product of all 8 prime polys")
print(f"Degree of P(x): {product_poly.degree()}")
print(f"P(x) = {product_poly.as_expr()}")

# GCD of char_poly and product
g = gcd(char_poly.as_expr(), product_poly.as_expr(), x)
print(f"\ngcd(char_poly, P(x)) = {g}")

if g != 1:
    print("*** SHARED FACTOR FOUND! ***")
else:
    print("No shared factor over Q.")

# Reduction of char_poly modulo P(x)
_, char_mod_product = div(char_poly, product_poly, x, domain='ZZ')
print(f"\nchar_poly mod P(x) = {char_mod_product.as_expr()}")

# Also: resultant
res = resultant(char_poly.as_expr(), product_poly.as_expr(), x)
print(f"\nResultant(char_poly, P(x)) = {res}")
print(f"  = {factorint(abs(int(res)))}")


# ============================================================
# PART 5: RECIPROCAL POLYNOMIAL (char poly of M^{-1})
# ============================================================
separator("PART 5: RECIPROCAL POLYNOMIAL (MIRROR)")

# reciprocal: x^8 * char_poly(1/x)
recip_expr = x**8 * char_poly.as_expr().subs(x, 1/x)
recip_expanded = expand(recip_expr)
recip_poly = Poly(recip_expanded, x, domain='ZZ')
print(f"Reciprocal poly: {recip_poly.as_expr()}")

# Check if x^2 - x - 1 divides reciprocal
_, rem_recip_golden = div(recip_poly, golden_poly, x, domain='ZZ')
print(f"\nreciprocal mod (x^2-x-1) = {rem_recip_golden.as_expr()}")
print(f"  Evaluated at phi: {float(rem_recip_golden.as_expr().subs(x, phi_val)):.10f}")

# Check all 8 prime polys against reciprocal
print("\nReciprocal poly mod each prime poly:")
for p in PRIMES:
    _, rem_r = div(recip_poly, prime_polys[p], x, domain='ZZ')
    print(f"  p={p:2d}: remainder = {rem_r.as_expr()}")

# Roots of reciprocal
recip_roots = nroots(recip_poly.as_expr(), n=15)
print("\nRoots of reciprocal poly (= 1/roots of char_poly):")
for i, r in enumerate(recip_roots):
    r_c = complex(r)
    is_real = abs(r_c.imag) < 1e-12
    marker = " <-- REAL" if is_real else ""
    print(f"  {r_c.real:+.12f} {r_c.imag:+.12f}i{marker}")


# ============================================================
# PART 6: AXIOM RESIDUAL M^2 - M - I
# ============================================================
separator("PART 6: AXIOM RESIDUAL M^2 - M - I")

I8 = eye(8)
residual = M**2 - M - I8
print("M^2 - M - I =")
pprint(residual)

print(f"\nFrobenius norm of residual: {float(sum(residual[i,j]**2 for i in range(8) for j in range(8))**0.5):.6f}")
print(f"Max absolute entry: {max(abs(residual[i,j]) for i in range(8) for j in range(8))}")
print(f"Trace of residual: {residual.trace()}")
print(f"Det of residual: {residual.det()}")
print(f"Rank of residual: {residual.rank()}")

# Char poly of the residual
res_char = Poly(residual.charpoly(x).as_expr(), x, domain='ZZ')
print(f"\nChar poly of (M^2-M-I): {res_char.as_expr()}")
print(f"Factored: {factor(res_char.as_expr(), x)}")

# Check if residual is nilpotent
print(f"\n(M^2-M-I)^2 zero? {(residual**2).is_zero_matrix}")
print(f"(M^2-M-I)^3 zero? {(residual**3).is_zero_matrix}")
print(f"(M^2-M-I)^4 zero? {(residual**4).is_zero_matrix}")

# Eigenvalues of residual
res_eigenvals = residual.eigenvals()
print(f"\nEigenvalues of (M^2-M-I):")
for ev, mult in res_eigenvals.items():
    print(f"  {ev} (multiplicity {mult})")


# ============================================================
# PART 7: FACTORIZATION OVER FINITE FIELDS
# ============================================================
separator("PART 7: FINITE FIELD FACTORIZATIONS")

test_moduli = [2, 3, 5, 7, 11, 13, 17, 19, 137]

golden_poly_expr = x**2 - x - 1

for p in test_moduli:
    print(f"\n--- mod {p} ---")

    try:
        # Factor char_poly mod p
        cp_modp = Poly(char_poly.as_expr(), x, modulus=p)
        factors_modp = cp_modp.factor_list()

        print(f"  char_poly mod {p} factors:")
        for fac, mult in factors_modp[1]:
            print(f"    ({fac.as_expr()})^{mult}")

        # Factor golden poly mod p
        gp_modp = Poly(golden_poly_expr, x, modulus=p)
        gp_factors = gp_modp.factor_list()
        print(f"  golden poly mod {p} factors:")
        for fac, mult in gp_factors[1]:
            print(f"    ({fac.as_expr()})^{mult}")

        # Check if golden poly divides char_poly mod p
        _, rem_modp = div(cp_modp, gp_modp, x)
        divides = rem_modp.is_zero if hasattr(rem_modp, 'is_zero') else (rem_modp == 0)

        # More robust check: see if any factor of char_poly matches golden poly factors
        char_factor_exprs = [fac.as_expr() for fac, mult in factors_modp[1]]
        golden_factor_exprs = [fac.as_expr() for fac, mult in gp_factors[1]]

        shared = []
        for cf in factors_modp[1]:
            for gf in gp_factors[1]:
                if cf[0] == gf[0]:
                    shared.append(cf[0].as_expr())

        if shared:
            print(f"  *** SHARED FACTOR(S) mod {p}: {shared} ***")
        elif divides:
            print(f"  *** x^2-x-1 DIVIDES char_poly mod {p}! ***")
        else:
            print(f"  No shared factors mod {p}")

    except Exception as e:
        print(f"  Error: {e}")

# Special: mod 5 analysis (golden prime)
separator("PART 7b: DEEP ANALYSIS mod 5")

print("mod 5:")
print(f"  x^2 - x - 1 mod 5 = x^2 + 4x + 4 = (x+2)^2 mod 5")

# Verify
gp_mod5 = Poly(x**2 - x - 1, x, modulus=5)
print(f"  Sympy factors of x^2-x-1 mod 5: {gp_mod5.factor_list()}")

cp_mod5 = Poly(char_poly.as_expr(), x, modulus=5)
print(f"  Char poly mod 5 factors: {cp_mod5.factor_list()}")

# Check if (x+2)^2 = (x-3)^2 divides char_poly mod 5
xp2_sq = Poly((x+2)**2, x, modulus=5)
_, rem5 = div(cp_mod5, xp2_sq, x)
print(f"  char_poly mod (x+2)^2 mod 5: remainder = {rem5}")

# Check (x+2) alone
xp2 = Poly(x+2, x, modulus=5)
_, rem5_lin = div(cp_mod5, xp2, x)
print(f"  char_poly mod (x+2) mod 5: remainder = {rem5_lin}")

# Evaluate char_poly at x = -2 mod 5 = 3 mod 5
val_at_3 = sum(c * pow(3, i, 5) for i, c in enumerate(reversed(char_poly.all_coeffs()))) % 5
print(f"  char_poly(3) mod 5 = {val_at_3}")
val_at_neg2 = sum(c * pow(-2, i, 5) for i, c in enumerate(reversed(char_poly.all_coeffs()))) % 5
print(f"  char_poly(-2) mod 5 = {val_at_neg2 % 5}")


# ============================================================
# PART 8: Cross-connections and summary
# ============================================================
separator("PART 8: CROSS-CONNECTIONS")

print("--- Resultants of char_poly with each prime poly ---")
for p in PRIMES:
    res_p = resultant(char_poly.as_expr(), prime_polys[p].as_expr(), x)
    print(f"  Res(char_poly, q_{p}) = {res_p} = {factorint(abs(int(res_p)))}")

print(f"\n--- Resultant with golden poly ---")
res_golden = resultant(char_poly.as_expr(), golden_poly_expr, x)
print(f"  Res(char_poly, x^2-x-1) = {res_golden} = {factorint(abs(int(res_golden)))}")

# Discriminant of char_poly
disc = char_poly.discriminant()
print(f"\n--- Discriminant of char_poly ---")
print(f"  disc = {disc}")
print(f"  factored: {factorint(abs(int(disc)))}")

# Check: is the discriminant divisible by 5? By primes in SHA?
for p in PRIMES:
    if disc % p == 0:
        print(f"  disc divisible by {p}: YES (multiplicity {0})")
        d = abs(int(disc))
        count = 0
        while d % p == 0:
            d //= p
            count += 1
        print(f"    {p}^{count} divides disc")
    else:
        print(f"  disc divisible by {p}: NO")

# Final check: evaluate char_poly at sqrt(p) for each p
print("\n--- char_poly evaluated at sqrt(p) ---")
for p in PRIMES:
    val = float(char_poly.as_expr().subs(x, p**0.5))
    print(f"  char_poly(sqrt({p:2d})) = {val:+.6f}")

print("\n--- char_poly evaluated at frac(sqrt(p)) ---")
for p in PRIMES:
    k = int(p**0.5)
    frac = p**0.5 - k
    val = float(char_poly.as_expr().subs(x, frac))
    print(f"  char_poly(frac(sqrt({p:2d}))) = {val:+.10f}")


# ============================================================
# PART 9: Galois group structure hints
# ============================================================
separator("PART 9: GALOIS GROUP STRUCTURE")

# The splitting field: examine how char_poly splits over extensions
# If it factors over Q(sqrt(d)) for some d, the Galois group is smaller

from sympy import sqrtdenest

print("Testing factorization over Q(sqrt(d)) for d = 2,3,5,7,11,13,17,19:")
for d in PRIMES:
    try:
        # Try to factor char_poly over Q(sqrt(d))
        # Approach: substitute x = a + b*sqrt(d), solve for rational a,b
        # Simpler: check if disc is a perfect square mod the extension

        # Method: compute char_poly mod (t^2 - d) treating coefficients in Q[t]/(t^2-d)
        # Actually let's just check: does char_poly have a root in Q(sqrt(d))?
        # A root in Q(sqrt(d)) means char_poly(a + b*sqrt(d)) = 0 for rational a,b
        # This expands to: P(a,b,d) + Q(a,b,d)*sqrt(d) = 0
        # So both P = 0 and Q = 0.

        # Instead: resultant approach. char_poly(x) and x = a + b*sqrt(d)
        # means (x-a)^2 = b^2*d, so x^2 - 2ax + a^2 - b^2*d = 0
        # For char_poly to have a root satisfying this quadratic:
        # gcd(char_poly, x^2 - 2ax + a^2 - b^2*d) != 1
        # This is hard to solve in general. Let's just check numerically.

        # Numeric check: do any roots of char_poly lie in Q(sqrt(d))?
        found = False
        for r in char_roots:
            r_c = complex(r)
            if abs(r_c.imag) < 1e-10:
                # Real root. Check if (r - a)/b = sqrt(d) for some rational a,b
                # i.e., is (r_c.real)^2 or (r_c.real - k)^2 = rational * d for small integers?
                r_val = r_c.real
                for a_num in range(-5, 6):
                    for a_den in range(1, 6):
                        a = a_num / a_den
                        diff = r_val - a
                        b_sq_d = diff**2
                        if abs(b_sq_d) < 1e-15:
                            continue
                        b_sq = b_sq_d / d
                        b = b_sq**0.5
                        # Check if b is rational (small fraction)
                        for b_den in range(1, 10):
                            b_num = b * b_den
                            if abs(b_num - round(b_num)) < 1e-8 and abs(round(b_num)) > 0:
                                found = True
                                a_exact = f"{a_num}/{a_den}" if a_den > 1 else str(a_num)
                                b_exact = f"{int(round(b_num))}/{b_den}" if b_den > 1 else str(int(round(b_num)))
                                # Verify
                                recon = a + (round(b_num)/b_den) * (d**0.5)
                                if abs(recon - r_val) < 1e-8:
                                    print(f"  d={d:2d}: root {r_val:.10f} = {a_exact} + ({b_exact})*sqrt({d})  *** HIT ***")
                                    found = True
                                    break
                            if found:
                                break
                    if found:
                        break
                if found:
                    break
        if not found:
            print(f"  d={d:2d}: no roots found in Q(sqrt({d}))")
    except Exception as e:
        print(f"  d={d:2d}: error - {e}")


# ============================================================
# PART 10: Summary statistics
# ============================================================
separator("SUMMARY")

print("SHA-256 Char Poly:  x^8 - 2x^7 + x^6 - x^4 - 1")
print(f"Determinant:        {det_M}")
print(f"Irreducible over Q: {char_poly.is_irreducible}")
print()

print("Divisibility by prime polynomials: NONE (all remainders nonzero)")
print()

print("Closest root to phi:")
min_delta = min(abs(complex(r).real - phi_val) for r in char_roots if abs(complex(r).imag) < 1e-10)
closest = [r for r in char_roots if abs(complex(r).imag) < 1e-10 and abs(complex(r).real - phi_val - min_delta) < 1e-12 or abs(complex(r).real - phi_val + min_delta) < 1e-12]
print(f"  delta = {min_delta:.15f}")
print(f"  as fraction of phi: {min_delta/phi_val:.10f}")
print()

print("Golden poly divides char_poly mod p:")
for p in test_moduli:
    try:
        cp_modp = Poly(char_poly.as_expr(), x, modulus=p)
        gp_modp = Poly(golden_poly_expr, x, modulus=p)
        _, rem_modp = div(cp_modp, gp_modp, x)
        divides = rem_modp.is_zero if hasattr(rem_modp, 'is_zero') else False
        print(f"  mod {p:3d}: {'YES ***' if divides else 'no'}")
    except:
        print(f"  mod {p:3d}: error")

print()
print("Axiom residual M^2 - M - I:")
print(f"  Rank: {residual.rank()}")
print(f"  Trace: {residual.trace()}")
print(f"  Frobenius norm: {float(sum(residual[i,j]**2 for i in range(8) for j in range(8))**0.5):.6f}")

print("\n" + "="*72)
print("  DONE")
print("="*72)
