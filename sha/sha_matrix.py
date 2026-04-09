"""
SHA_MATRIX — constructs the linearized SHA-256 round matrix; tests whether char_poly contains x^2-x-1
nos3bl33d

Multiple linearization strategies. If golden poly divides char_poly: Fibonacci reduction applies.
"""

import sympy
from sympy import (
    Matrix, Symbol, Rational, factor, rem, Poly, sqrt,
    eye, zeros, pprint, simplify, GF, ZZ, QQ, latex
)
from sympy.polys.polytools import gcd as poly_gcd
import numpy as np
from fractions import Fraction

x = Symbol('x')
phi = (1 + sqrt(5)) / 2
golden_poly = x**2 - x - 1

# ==========================================================================
# PART 1: SHA-256 Round Function Matrices (Multiple Linearizations)
# ==========================================================================

def banner(title):
    sep = "=" * 72
    print(f"\n{sep}")
    print(f"  {title}")
    print(sep)

def analyze_matrix(name, M):
    """Full characteristic polynomial analysis for a matrix."""
    banner(f"Analysis: {name}")
    print(f"\nMatrix ({M.rows}x{M.cols}):")
    pprint(M)

    # Characteristic polynomial
    cp = M.charpoly(x)
    cp_expr = cp.as_expr()
    print(f"\nCharacteristic polynomial:")
    print(f"  {cp_expr}")

    # Factor over Q
    factored = factor(cp_expr)
    print(f"\nFactored over Q:")
    print(f"  {factored}")

    # Factor as polynomial object for clean factor list
    cp_poly = Poly(cp_expr, x, domain='ZZ')
    try:
        factor_list = cp_poly.factor_list()
        print(f"\nFactor list:")
        for fac, mult in factor_list[1]:
            print(f"  ({fac.as_expr()})^{mult}")
    except Exception as e:
        print(f"  (factor_list failed: {e})")

    # Check divisibility by x^2 - x - 1
    gp = Poly(golden_poly, x, domain='ZZ')
    remainder = rem(cp_expr, golden_poly, x)
    print(f"\nDivisibility by x^2 - x - 1 (golden polynomial):")
    print(f"  Remainder: {remainder}")
    print(f"  DIVIDES: {simplify(remainder) == 0}")

    # Also check x^2 + x - 1 (the "negative" golden poly)
    neg_golden = x**2 + x - 1
    remainder_neg = rem(cp_expr, neg_golden, x)
    print(f"\nDivisibility by x^2 + x - 1:")
    print(f"  Remainder: {remainder_neg}")
    print(f"  DIVIDES: {simplify(remainder_neg) == 0}")

    # GCD with golden poly
    g = poly_gcd(Poly(cp_expr, x), gp)
    print(f"\nGCD(char_poly, x^2-x-1) = {g.as_expr()}")

    # Eigenvalues (symbolic)
    print(f"\nEigenvalues (symbolic):")
    try:
        eigs = M.eigenvals()
        for eig, mult in eigs.items():
            eig_simplified = simplify(eig)
            print(f"  {eig_simplified}  (multiplicity {mult})")
            # Check if eigenvalue equals phi or -1/phi
            diff_phi = simplify(eig - phi)
            diff_neg_phi = simplify(eig + 1/phi)
            if diff_phi == 0:
                print(f"    *** THIS IS PHI! ***")
            if diff_neg_phi == 0:
                print(f"    *** THIS IS -1/phi! ***")
    except Exception as e:
        print(f"  (symbolic eigenvalues failed: {e})")

    # Numerical eigenvalues
    print(f"\nEigenvalues (numerical):")
    M_float = np.array(M.tolist(), dtype=float)
    eigvals = np.linalg.eigvals(M_float)
    for ev in sorted(eigvals, key=lambda z: (-abs(z), z.real)):
        if abs(ev.imag) < 1e-12:
            print(f"  {ev.real:.10f}  |eigenvalue| = {abs(ev.real):.10f}")
        else:
            print(f"  {ev.real:.10f} + {ev.imag:.10f}i  |eigenvalue| = {abs(ev):.10f}")

    # Check if any eigenvalue modulus equals sqrt(5)
    sqrt5 = np.sqrt(5)
    phi_val = (1 + sqrt5) / 2
    print(f"\nEigenvalue moduli vs golden constants:")
    print(f"  phi = {phi_val:.10f}")
    print(f"  1/phi = {1/phi_val:.10f}")
    print(f"  sqrt(5) = {sqrt5:.10f}")
    for i, ev in enumerate(eigvals):
        mod = abs(ev)
        if abs(mod - phi_val) < 1e-6:
            print(f"  eigenvalue {i}: |lambda| = {mod:.10f} ~ phi!")
        if abs(mod - 1/phi_val) < 1e-6:
            print(f"  eigenvalue {i}: |lambda| = {mod:.10f} ~ 1/phi!")
        if abs(mod - sqrt5) < 1e-6:
            print(f"  eigenvalue {i}: |lambda| = {mod:.10f} ~ sqrt(5)!")

    # Determinant and trace
    det_val = M.det()
    tr_val = M.trace()
    print(f"\ndet(M) = {det_val}")
    print(f"tr(M) = {tr_val}")
    print(f"det relates to 5? det = {det_val}, det mod 5 = {det_val % 5 if isinstance(det_val, (int, sympy.Integer)) else 'N/A'}")

    return cp_expr, eigvals


# --------------------------------------------------------------------------
# VERSION A: M_simple (Ch=0, Maj=0, rotations=identity)
# --------------------------------------------------------------------------
# a_new = a + e + h
# b_new = a
# c_new = b
# d_new = c
# e_new = d + e + h
# f_new = e
# g_new = f
# h_new = g

M_simple = Matrix([
    [1, 0, 0, 0, 1, 0, 0, 1],   # a = a + e + h
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 1, 0, 0, 1],   # e = d + e + h
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])

# --------------------------------------------------------------------------
# VERSION B: Ch(e,f,g) ~ g, Maj(a,b,c) ~ a (first-argument linear approx)
# --------------------------------------------------------------------------
# T1 = h + e + g  (Sigma1(e)~e, Ch(e,f,g)~g)
# T2 = a + a = 2a  (Sigma0(a)~a, Maj(a,b,c)~a) ... but in GF(2), 2a=0. Over Z, 2a.
# For Z/2^32: keep 2a.
# a_new = T1 + T2 = 2a + e + g + h
# e_new = d + T1 = d + e + g + h

M_version_B = Matrix([
    [2, 0, 0, 0, 1, 0, 1, 1],   # a = 2a + e + g + h
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 1, 0, 1, 1],   # e = d + e + g + h
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])

# --------------------------------------------------------------------------
# VERSION C: Pure permutation/shift only (no additions at all)
# --------------------------------------------------------------------------
# a_new = 0 (no input -- but that's degenerate)
# Actually: the pure shift is a=a, e=d, rest shifts.
# If we strip ALL additions, just the permutation structure:
# a->b->c->d, e->f->g->h, with a and e getting "input"
# Pure permutation: b=a, c=b, d=c, f=e, g=f, h=g, a=a, e=d
# This gives a permutation matrix

M_perm = Matrix([
    [1, 0, 0, 0, 0, 0, 0, 0],   # a = a (identity, no additions)
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 0, 0, 0, 0],   # e = d
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])

# --------------------------------------------------------------------------
# VERSION D: Full linear approximation with Ch(e,f,g)=e+f+g, Maj(a,b,c)=a+b+c
# --------------------------------------------------------------------------
# T1 = h + e + (e + f + g) = h + 2e + f + g  (over Z)
# T2 = a + (a + b + c) = 2a + b + c  (over Z)
# a_new = T1 + T2 = 2a + b + c + 2e + f + g + h
# e_new = d + T1 = d + 2e + f + g + h

M_full_linear = Matrix([
    [2, 1, 1, 0, 2, 1, 1, 1],   # a = 2a + b + c + 2e + f + g + h
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 2, 1, 1, 1],   # e = d + 2e + f + g + h
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])

# --------------------------------------------------------------------------
# VERSION E: GF(2) linearization (Ch=e+f+g, Maj=a+b+c, all mod 2)
# --------------------------------------------------------------------------
# Over GF(2): 2x = 0 for any x
# T1 = h + e + e + f + g = h + f + g  (2e vanishes)
# T2 = a + a + b + c = b + c  (2a vanishes)
# a_new = T1 + T2 = b + c + f + g + h
# e_new = d + T1 = d + f + g + h

M_gf2 = Matrix([
    [0, 1, 1, 0, 0, 1, 1, 1],   # a = b + c + f + g + h (mod 2)
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 0, 1, 1, 1],   # e = d + f + g + h (mod 2)
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])

# --------------------------------------------------------------------------
# VERSION F: Minimal "feedback" matrix -- just the shift + h-to-a feedback
# --------------------------------------------------------------------------
# Stripping everything except the shift structure and the h->a,e feedback:
# a = a + h, e = d + h, rest = shifts. Simplest nontrivial version.

M_feedback = Matrix([
    [1, 0, 0, 0, 0, 0, 0, 1],   # a = a + h
    [1, 0, 0, 0, 0, 0, 0, 0],   # b = a
    [0, 1, 0, 0, 0, 0, 0, 0],   # c = b
    [0, 0, 1, 0, 0, 0, 0, 0],   # d = c
    [0, 0, 0, 1, 0, 0, 0, 1],   # e = d + h
    [0, 0, 0, 0, 1, 0, 0, 0],   # f = e
    [0, 0, 0, 0, 0, 1, 0, 0],   # g = f
    [0, 0, 0, 0, 0, 0, 1, 0]    # h = g
])


# ==========================================================================
# PART 2: Run All Analyses
# ==========================================================================

matrices = [
    ("Version A: M_simple (Ch=0, Maj=0, rotations=id)", M_simple),
    ("Version B: Ch~g, Maj~a", M_version_B),
    ("Version C: Pure permutation/shift", M_perm),
    ("Version D: Full linear (Ch=e+f+g, Maj=a+b+c) over Z", M_full_linear),
    ("Version E: GF(2) linearization", M_gf2),
    ("Version F: Minimal feedback (shift + h->a,e)", M_feedback),
]

all_results = {}
for name, mat in matrices:
    cp, eigs = analyze_matrix(name, mat)
    all_results[name] = (cp, eigs)


# ==========================================================================
# PART 3: The Golden Constant Connection -- h2 = sqrt(5)
# ==========================================================================

banner("PART 3: SHA-256 Initial Values & Golden Ratio")

# SHA-256 initial hash values (fractional parts of sqrt of first 8 primes)
h_init = [
    0x6A09E667,  # sqrt(2)
    0xBB67AE85,  # sqrt(3)
    0x3C6EF372,  # sqrt(5) <-- PHI CONNECTION
    0xA54FF53A,  # sqrt(7)
    0x510E527F,  # sqrt(11)
    0x9B05688C,  # sqrt(13)
    0x1F83D9AB,  # sqrt(17)
    0x5BE0CD19,  # sqrt(19)
]

print("\nSHA-256 initial hash values:")
for i, h in enumerate(h_init):
    frac_val = h / (2**32)
    prime_idx = [2, 3, 5, 7, 11, 13, 17, 19][i]
    actual_frac = float(sympy.sqrt(prime_idx)) - int(float(sympy.sqrt(prime_idx)))
    print(f"  h{i} = 0x{h:08X} = {h:>10d}  frac(sqrt({prime_idx:2d})) = {frac_val:.10f}  (actual: {actual_frac:.10f})")

h2_val = h_init[2]
h2_frac = h2_val / (2**32)
print(f"\nKey value h2 = 0x{h2_val:08X}")
print(f"  h2/2^32 = {h2_frac:.15f}")
print(f"  sqrt(5) - 2 = {float(sympy.sqrt(5)) - 2:.15f}")
print(f"  1/phi^3 = {float(1/phi**3):.15f}")
print(f"  sqrt(5) = 2*phi - 1, so frac(sqrt(5)) = 2*phi - 3 = 2*(1+sqrt(5))/2 - 3 = sqrt(5) - 2")
print(f"  And sqrt(5) - 2 = 1/(sqrt(5)+2) = (sqrt(5)-2) ... = phi^(-3) + phi^(-4)")

# The actual check: is sqrt(5)-2 = 1/phi^3?
# phi^3 = phi^2 * phi = (phi+1)*phi = phi^2 + phi = (phi+1) + phi = 2*phi + 1 = sqrt(5) + 2
# So 1/phi^3 = 1/(sqrt(5)+2) = (sqrt(5)-2)/((sqrt(5))^2 - 4) = (sqrt(5)-2)/1 = sqrt(5)-2  YES!
print(f"\n  CONFIRMED: h2/2^32 = sqrt(5) - 2 = 1/phi^3")
print(f"  phi^3 = 2*phi + 1 = sqrt(5) + 2 = {float(phi**3):.10f}")
print(f"  1/phi^3 = sqrt(5) - 2 = {float(1/phi**3):.10f}")

# Golden hash constant check
golden_const = 0x9E3779B9
print(f"\n  MurmurHash golden constant: 0x{golden_const:08X} = {golden_const}")
print(f"  2^32/phi = {2**32/float(phi):.2f}")
print(f"  floor(2^32/phi) = {int(2**32/float(phi))}")
print(f"  0x9E3779B9 = {golden_const} vs floor(2^32/phi) = {int(2**32/float(phi))}")
print(f"  Match: {golden_const == int(2**32/float(phi))}")

# Also: 2^32 * (sqrt(5)-1)/2 = 2^32 / phi = the golden constant
# And: 2^32 * (sqrt(5)-2) = h2. So h2 = 2^32/phi^3 = golden_const / phi^2
print(f"\n  h2 = 2^32/phi^3")
print(f"  golden_const = 2^32/phi")
print(f"  Ratio: golden_const/h2 = phi^2 = {golden_const/h2_val:.10f}")
print(f"  phi^2 = {float(phi**2):.10f}")
print(f"  Match: {abs(golden_const/h2_val - float(phi**2)) < 1e-6}")


# ==========================================================================
# PART 4: Pythagorean Matrix Connection
# ==========================================================================

banner("PART 4: Pythagorean Matrix Connection")

M_pyth = Matrix([[1, 2], [-2, 1]])
print("\nPythagorean matrix M_pyth:")
pprint(M_pyth)
cp_pyth = M_pyth.charpoly(x).as_expr()
print(f"char poly: {cp_pyth}")
print(f"det = {M_pyth.det()}")
print(f"eigenvalues: {M_pyth.eigenvals()}")

# Check: does M_simple decompose as tensor product involving M_pyth?
# M_pyth is 2x2 with det=5. If M_simple = M_pyth tensor M_4x4, then det(M_simple) = det(M_pyth)^4 * det(M_4x4)^2
# Or: M_simple = M_4x4 tensor M_pyth, then det(M_simple) = det(M_4x4)^2 * det(M_pyth)^4

det_simple = M_simple.det()
print(f"\ndet(M_simple) = {det_simple}")
print(f"det(M_pyth)^4 = {M_pyth.det()**4}")
print(f"Is det(M_simple) divisible by 5? {det_simple % 5 == 0 if isinstance(det_simple, (int, sympy.Integer)) else 'check needed'}")

# Block structure analysis
# M_simple has a 2x2 "connection" structure:
# Upper-left 4x4 (a,b,c,d) and lower-right 4x4 (e,f,g,h)
# connected by off-diagonal blocks
print("\nBlock decomposition of M_simple:")
A11 = M_simple[:4, :4]
A12 = M_simple[:4, 4:]
A21 = M_simple[4:, :4]
A22 = M_simple[4:, 4:]
print("Upper-left 4x4 (a,b,c,d self-interaction):")
pprint(A11)
print("Upper-right 4x4 (e,f,g,h -> a,b,c,d):")
pprint(A12)
print("Lower-left 4x4 (a,b,c,d -> e,f,g,h):")
pprint(A21)
print("Lower-right 4x4 (e,f,g,h self-interaction):")
pprint(A22)

# The 2x2 "super-structure": treat each 4x4 block as a scalar
# The structure is: [[A11, A12], [A21, A22]]
# where A11 and A22 are shift matrices with a self-loop on position 0,
# and A12, A21 are coupling matrices
print("\nA11 eigenvalues:", [complex(e) for e in A11.eigenvals().keys()])
print("A22 eigenvalues:", [complex(e) for e in A22.eigenvals().keys()])


# ==========================================================================
# PART 5: Fibonacci Structure in Powers of M_simple
# ==========================================================================

banner("PART 5: Fibonacci Structure in M_simple^n")

# If M satisfies x^2 = x + 1 (Cayley-Hamilton for a 2x2 with char poly x^2-x-1),
# then M^n = F_n * M + F_{n-1} * I
# For 8x8, the char poly is degree 8, so the recurrence is 8-term, not 2-term.
# But: if x^2-x-1 DIVIDES the char poly, then the RESTRICTION of M to the
# corresponding eigenspace satisfies the Fibonacci recurrence.

# Let's compute M^n and check for Fibonacci patterns
def fib(n):
    """Fibonacci: F(0)=0, F(1)=1, F(2)=1, ..."""
    if n <= 0:
        return 0
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b

print("\nFibonacci numbers F(1) through F(20):")
fibs = [fib(n) for n in range(21)]
print(f"  {fibs}")

# Compute M_simple^n for n=1..20 and check if M^n = F_n*A + F_{n-1}*B
# for some fixed A, B
print("\nComputing M_simple^n for n=1..20...")
M_pow = [None] * 21
M_pow[0] = eye(8)
M_pow[1] = M_simple

for n in range(2, 21):
    M_pow[n] = M_simple * M_pow[n-1]

# If M^n = F_n*A + F_{n-1}*B, then:
# M^1 = F_1*A + F_0*B = A  => A = M
# M^2 = F_2*A + F_1*B = A + B => B = M^2 - M
# Check: M^3 = F_3*A + F_2*B = 2M + (M^2 - M) = M + M^2 = ?

A_fib = M_pow[1]
B_fib = M_pow[2] - M_pow[1]
print(f"\n  If M^n = F_n*M + F_{'{'}n-1{'}'} * B where B = M^2 - M:")
print(f"  B = M^2 - M =")
pprint(B_fib)

fib_check_passes = True
for n in range(1, 21):
    predicted = fib(n) * A_fib + fib(n-1) * B_fib
    actual = M_pow[n]
    diff = actual - predicted
    if diff != zeros(8):
        fib_check_passes = False
        if n <= 5:  # Only show first few failures
            print(f"\n  M^{n}: Fibonacci prediction FAILS")
            print(f"    Max diff entry: {max(abs(e) for e in diff)}")

if fib_check_passes:
    print(f"\n  *** FIBONACCI REDUCTION WORKS! M^n = F_n*M + F_(n-1)*(M^2-M) for all n=1..20 ***")
else:
    print(f"\n  Fibonacci 2-term reduction does NOT hold (expected for 8x8 matrix).")
    print(f"  The minimal polynomial is likely degree > 2.")

    # Find the actual minimal polynomial
    print(f"\n  Computing minimal polynomial of M_simple...")
    try:
        min_poly = M_simple.minpoly(x)
        print(f"  Minimal polynomial: {min_poly.as_expr()}")
        print(f"  Degree: {Poly(min_poly.as_expr(), x).degree()}")

        # Factor minimal polynomial
        min_factored = factor(min_poly.as_expr())
        print(f"  Factored: {min_factored}")

        # Check golden poly divisibility of minimal poly
        min_rem = rem(min_poly.as_expr(), golden_poly, x)
        print(f"  Remainder of min_poly / (x^2-x-1): {min_rem}")
        print(f"  x^2-x-1 divides minimal poly: {simplify(min_rem) == 0}")
    except Exception as e:
        print(f"  (minimal polynomial computation failed: {e})")


# ==========================================================================
# PART 6: The Cayley-Hamilton Recurrence
# ==========================================================================

banner("PART 6: Cayley-Hamilton Recurrence for M_simple")

# The characteristic polynomial gives a recurrence for M^n.
# If char_poly = x^8 + c7*x^7 + ... + c0, then
# M^8 + c7*M^7 + ... + c0*I = 0
# So M^8 = -c7*M^7 - ... - c0*I

cp_simple = M_simple.charpoly(x).as_expr()
cp_poly_simple = Poly(cp_simple, x)
coeffs = cp_poly_simple.all_coeffs()  # [leading, ..., constant]
print(f"\nChar poly coefficients (x^8 down to x^0):")
for i, c in enumerate(coeffs):
    print(f"  x^{8-i}: {c}")

# Verify Cayley-Hamilton
print(f"\nVerifying Cayley-Hamilton M^8 + c7*M^7 + ... + c0*I = 0:")
CH_sum = zeros(8)
for i, c in enumerate(coeffs):
    power = 8 - i
    CH_sum += c * M_pow[power] if power <= 20 else c * (M_simple**power)
print(f"  Sum = zero matrix: {CH_sum == zeros(8)}")


# ==========================================================================
# PART 7: GF(2) Analysis -- The TRUE Linear Structure
# ==========================================================================

banner("PART 7: GF(2) Characteristic Polynomial of M_gf2")

# Over GF(2), the characteristic polynomial factors differently.
# This is the mathematically "correct" linearization for XOR-based operations.

print("\nM_gf2 (mod 2):")
M_gf2_mod2 = M_gf2.applyfunc(lambda e: e % 2)
pprint(M_gf2_mod2)

# Compute char poly over GF(2) using sympy
print("\nCharacteristic polynomial of M_gf2 over GF(2):")
# Work with the integer matrix, compute char poly, then reduce mod 2
cp_gf2_Z = M_gf2_mod2.charpoly(x).as_expr()
print(f"  Over Z: {cp_gf2_Z}")

# Reduce coefficients mod 2
cp_gf2_poly = Poly(cp_gf2_Z, x, domain='ZZ')
coeffs_gf2 = cp_gf2_poly.all_coeffs()
coeffs_mod2 = [c % 2 for c in coeffs_gf2]
cp_mod2_expr = sum(c * x**(len(coeffs_mod2)-1-i) for i, c in enumerate(coeffs_mod2))
print(f"  Over GF(2): {cp_mod2_expr}")

# Factor over GF(2) using sympy's finite field support
print(f"\n  Factoring over GF(2):")
try:
    cp_gf2_ff = Poly(cp_mod2_expr, x, modulus=2)
    gf2_factors = cp_gf2_ff.factor_list()
    print(f"  Factor list:")
    for fac, mult in gf2_factors[1]:
        print(f"    ({fac.as_expr()})^{mult}")
except Exception as e:
    print(f"  (GF(2) factorization failed: {e})")

# Check if x^2 + x + 1 divides (which is x^2-x-1 over GF(2), since -1=1 mod 2)
golden_gf2 = x**2 + x + 1  # x^2 - x - 1 = x^2 + x + 1 (mod 2)
print(f"\n  Over GF(2): x^2-x-1 = x^2+x+1")
try:
    rem_gf2 = Poly(cp_mod2_expr, x, modulus=2).rem(Poly(golden_gf2, x, modulus=2))
    print(f"  Remainder of char_poly / (x^2+x+1) over GF(2): {rem_gf2.as_expr()}")
    print(f"  x^2+x+1 divides char_poly over GF(2): {rem_gf2.is_zero}")
except Exception as e:
    print(f"  (GF(2) division failed: {e})")

# Also: order of M_gf2 over GF(2) -- the smallest n such that M^n = I mod 2
print(f"\n  Computing order of M_gf2 mod 2 (smallest n with M^n = I mod 2):")
I8 = eye(8)
M_curr = M_gf2_mod2.copy()
order = None
for n in range(1, 300):
    M_curr_mod2 = M_curr.applyfunc(lambda e: e % 2)
    if M_curr_mod2 == I8:
        order = n
        break
    M_curr = (M_gf2_mod2 * M_curr).applyfunc(lambda e: e % 2)

if order:
    print(f"  Order = {order}")
    # Check if order is a Pisano period
    # Pisano period pi(2) = 3 (Fibonacci mod 2 repeats every 3)
    print(f"  Pisano period pi(2) = 3")
    print(f"  Order mod Pisano period: {order % 3}")
    print(f"  Order is multiple of pi(2): {order % 3 == 0}")
    # Pisano period pi(2^32) = 3 * 2^31
    pisano_2_32 = 3 * (2**31)
    print(f"  Pisano period pi(2^32) = 3 * 2^31 = {pisano_2_32}")
else:
    print(f"  Order > 299 (or does not exist)")


# ==========================================================================
# PART 8: Parametric Analysis -- When DOES x^2-x-1 divide char poly?
# ==========================================================================

banner("PART 8: Parametric Search for Golden Structure")

# Parametric matrix: vary the feedback coefficients
# M(alpha, beta, gamma) where:
# a_new = alpha*a + beta*e + gamma*h
# e_new = delta*d + epsilon*e + gamma*h
# rest = shifts

alpha, beta, gamma, delta, epsilon = sympy.symbols('alpha beta gamma delta epsilon')

M_param = Matrix([
    [alpha, 0, 0, 0, beta, 0, 0, gamma],
    [1,     0, 0, 0, 0,    0, 0, 0],
    [0,     1, 0, 0, 0,    0, 0, 0],
    [0,     0, 1, 0, 0,    0, 0, 0],
    [0,     0, 0, 1, epsilon, 0, 0, gamma],
    [0,     0, 0, 0, 1,    0, 0, 0],
    [0,     0, 0, 0, 0,    1, 0, 0],
    [0,     0, 0, 0, 0,    0, 1, 0]
])

print("\nParametric matrix M(alpha, beta, gamma, epsilon):")
print("  a_new = alpha*a + beta*e + gamma*h")
print("  e_new = d + epsilon*e + gamma*h")
print("  b=a, c=b, d=c, f=e, g=f, h=g")

# For M_simple: alpha=1, beta=1, gamma=1, epsilon=1
# For x^2-x-1 to divide the char poly, we need phi to be an eigenvalue.
# det(M - phi*I) = 0

print("\nComputing det(M_param - phi*I) symbolically...")
print("(This finds the constraint on alpha,beta,gamma,epsilon for phi to be an eigenvalue)")

# This is expensive symbolically. Let's do a numerical search instead.
print("\nNumerical search: scanning (alpha, beta, gamma, epsilon) in {-2..2}^4")
print("for matrices whose char poly is divisible by x^2-x-1...\n")

phi_val = (1 + np.sqrt(5)) / 2
golden_roots = [phi_val, -1/phi_val]  # roots of x^2-x-1

hits = []
for a_val in range(-2, 3):
    for b_val in range(-2, 3):
        for g_val in range(-2, 3):
            for e_val in range(-2, 3):
                M_test = np.array([
                    [a_val, 0, 0, 0, b_val, 0, 0, g_val],
                    [1,     0, 0, 0, 0,     0, 0, 0],
                    [0,     1, 0, 0, 0,     0, 0, 0],
                    [0,     0, 1, 0, 0,     0, 0, 0],
                    [0,     0, 0, 1, e_val,  0, 0, g_val],
                    [0,     0, 0, 0, 1,     0, 0, 0],
                    [0,     0, 0, 0, 0,     1, 0, 0],
                    [0,     0, 0, 0, 0,     0, 1, 0]
                ], dtype=float)
                eigvals = np.linalg.eigvals(M_test)
                # Check if phi is an eigenvalue (within tolerance)
                for gr in golden_roots:
                    diffs = np.abs(eigvals - gr)
                    if np.min(diffs) < 1e-8:
                        hits.append((a_val, b_val, g_val, e_val, eigvals))
                        break

if hits:
    print(f"Found {len(hits)} parameter sets with golden eigenvalue!\n")
    for h in hits[:10]:  # Show first 10
        a_v, b_v, g_v, e_v, evs = h
        print(f"  alpha={a_v}, beta={b_v}, gamma={g_v}, epsilon={e_v}")
        real_evs = sorted([ev.real if abs(ev.imag) < 1e-10 else complex(ev) for ev in evs],
                         key=lambda z: -abs(z) if isinstance(z, float) else -abs(z))
        print(f"    eigenvalues: {[f'{e:.6f}' if isinstance(e, float) else f'{e:.6f}' for e in evs]}")
        # For these, compute char poly symbolically and verify
        M_sym = Matrix([
            [a_v, 0, 0, 0, b_v, 0, 0, g_v],
            [1,   0, 0, 0, 0,   0, 0, 0],
            [0,   1, 0, 0, 0,   0, 0, 0],
            [0,   0, 1, 0, 0,   0, 0, 0],
            [0,   0, 0, 1, e_v, 0, 0, g_v],
            [0,   0, 0, 0, 1,   0, 0, 0],
            [0,   0, 0, 0, 0,   1, 0, 0],
            [0,   0, 0, 0, 0,   0, 1, 0]
        ])
        cp_h = M_sym.charpoly(x).as_expr()
        rem_h = rem(cp_h, golden_poly, x)
        fac_h = factor(cp_h)
        print(f"    char poly: {cp_h}")
        print(f"    factored: {fac_h}")
        print(f"    remainder mod (x^2-x-1): {rem_h}")
        print(f"    GOLDEN DIVISOR: {simplify(rem_h) == 0}")
        print()
    if len(hits) > 10:
        print(f"  ... and {len(hits) - 10} more.")
else:
    print("No parameter sets found with golden eigenvalue in {-2..2}^4.")
    print("The SHA-256 round structure does NOT naturally produce golden eigenvalues")
    print("in this parametric family.")


# ==========================================================================
# PART 9: What DOES the SHA round matrix characteristic poly look like?
# ==========================================================================

banner("PART 9: Detailed Structure of M_simple Characteristic Polynomial")

cp_expr_simple = all_results[list(all_results.keys())[0]][0]
print(f"\nChar poly of M_simple: {cp_expr_simple}")

# Check for OTHER notable polynomial divisors
notable_polys = {
    "x^2 - x - 1 (golden)": x**2 - x - 1,
    "x^2 + x - 1 (neg golden)": x**2 + x - 1,
    "x^2 + 1 (Gaussian)": x**2 + 1,
    "x^2 - 2 (sqrt 2)": x**2 - 2,
    "x^2 - 3 (sqrt 3)": x**2 - 3,
    "x^2 - 5 (sqrt 5)": x**2 - 5,
    "x^2 - x + 1 (6th root of unity)": x**2 - x + 1,
    "x^2 + x + 1 (cube root of unity)": x**2 + x + 1,
    "x^2 - 2x - 1 (silver ratio)": x**2 - 2*x - 1,
    "x^2 - 3x + 1": x**2 - 3*x + 1,
    "x^4 - x^3 - 1": x**4 - x**3 - 1,
    "x^4 + 1": x**4 + 1,
}

print(f"\nDivisibility check against notable polynomials:")
for name, poly in notable_polys.items():
    r = rem(cp_expr_simple, poly, x)
    divides = simplify(r) == 0
    if divides:
        print(f"  {name}: DIVIDES!! remainder = {r}")
    else:
        print(f"  {name}: no (remainder = {r})")

# Roots of char poly (all of them, precisely)
print(f"\nAll roots of char poly (via sympy solve):")
roots = sympy.solve(cp_expr_simple, x)
for i, r in enumerate(roots):
    r_simplified = simplify(r)
    r_float = complex(r_simplified)
    print(f"  root {i}: {r_simplified}")
    print(f"          = {r_float}")
    # Check: is this root related to phi?
    ratio_to_phi = simplify(r / phi)
    if ratio_to_phi.is_rational:
        print(f"          ratio to phi: {ratio_to_phi} (RATIONAL!)")


# ==========================================================================
# PART 10: Comprehensive Summary
# ==========================================================================

banner("COMPREHENSIVE SUMMARY")

print("""
Question: Does x^2 - x - 1 (the golden polynomial) divide the characteristic
polynomial of the linearized SHA-256 round matrix?

Results for each linearization:
""")

for name, (cp, eigs) in all_results.items():
    r = rem(cp, golden_poly, x)
    divides = simplify(r) == 0
    det_val = None
    for nm, mat in matrices:
        if nm == name:
            det_val = mat.det()
            break
    print(f"  {name}:")
    print(f"    char poly = {cp}")
    print(f"    det = {det_val}")
    print(f"    x^2-x-1 divides: {'YES!!' if divides else 'NO'}")
    print(f"    remainder: {r}")
    print()

print("""
Connection to SHA-256 initial values:
  h2 = floor(frac(sqrt(5)) * 2^32) = 0x3C6EF372
  sqrt(5) - 2 = 1/phi^3
  So h2 encodes phi^(-3) in the initial state.
  The MurmurHash constant 0x9E3779B9 = floor(2^32/phi) is NOT used in SHA-256.
  But phi appears implicitly through sqrt(5) in h2.
""")

print("="*72)
print("  ANALYSIS COMPLETE")
print("="*72)
