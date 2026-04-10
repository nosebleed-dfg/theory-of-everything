"""
Exact Tr(L^+) computation for 4D regular polytopes.
Uses sympy for exact rational arithmetic where eigenvalues are known exactly.
For irrational eigenvalues, identifies them and uses Galois conjugate cancellation.
"""

from fractions import Fraction
from sympy import isprime, factorint, sqrt, Rational, simplify
import numpy as np

phi_f = (1 + np.sqrt(5)) / 2

print("=" * 70)
print("EXACT TRACE COMPUTATIONS FOR 4D REGULAR POLYTOPES")
print("=" * 70)
print()

# ============================================================
# 1. 5-CELL
# ============================================================
print("1. 5-CELL (4-simplex)")
print("   Eigenvalues: 0(1), 5(4)")
tr = Fraction(4, 5)
print(f"   Tr(L^+) = 4/5")
print(f"   Numerator = {tr.numerator}")
print(f"   Prime? {isprime(tr.numerator)}")
if not isprime(tr.numerator):
    print(f"   Factorization: {factorint(tr.numerator)}")
print()

# ============================================================
# 2. 8-CELL (TESSERACT)
# ============================================================
print("2. 8-CELL (Tesseract)")
print("   Eigenvalues: 0(1), 2(4), 4(6), 6(4), 8(1)")
tr = Fraction(4, 2) + Fraction(6, 4) + Fraction(4, 6) + Fraction(1, 8)
print(f"   Tr(L^+) = 4/2 + 6/4 + 4/6 + 1/8 = {tr}")
print(f"   Numerator = {tr.numerator}")
print(f"   Prime? {isprime(tr.numerator)}")
print()

# ============================================================
# 3. 16-CELL
# ============================================================
print("3. 16-CELL")
print("   Eigenvalues: 0(1), 6(4), 8(3)")
tr = Fraction(4, 6) + Fraction(3, 8)
print(f"   Tr(L^+) = 4/6 + 3/8 = {tr}")
print(f"   Numerator = {tr.numerator}")
print(f"   Prime? {isprime(tr.numerator)}")
if not isprime(tr.numerator):
    print(f"   Factorization: {factorint(tr.numerator)}")
print()

# ============================================================
# 4. 24-CELL
# ============================================================
print("4. 24-CELL")
print("   Eigenvalues: 0(1), 4(4), 8(9), 10(8), 12(2)")
tr = Fraction(4, 4) + Fraction(9, 8) + Fraction(8, 10) + Fraction(2, 12)
print(f"   Tr(L^+) = 4/4 + 9/8 + 8/10 + 2/12")
print(f"           = 1 + 9/8 + 4/5 + 1/6")
print(f"           = {tr}")
print(f"   Numerator = {tr.numerator}")
print(f"   Prime? {isprime(tr.numerator)}")
if not isprime(tr.numerator):
    print(f"   Factorization: {factorint(tr.numerator)}")
print()

# ============================================================
# 5. 600-CELL
# ============================================================
print("5. 600-CELL")
print("   Numerical eigenvalue groups:")
# From the computation:
# -0.0000000000 x 1
# 2.2917960675 x 4    -> need to identify
# 5.5278640450 x 9    -> need to identify
# 9.0000000000 x 16
# 12.0000000000 x 25
# 14.0000000000 x 36
# 14.4721359550 x 9   -> need to identify
# 15.0000000000 x 16
# 15.7082039325 x 4   -> need to identify

# The irrational eigenvalues involve phi = (1+sqrt(5))/2
# Let me identify them:
# 2.2917960675: let's check if this is 3 - phi = 3 - 1.618 = 1.382, no.
# Try: 6 - 2*sqrt(5) = 6 - 4.472 = 1.528, no.
# Try: phi^2 = phi + 1 = 2.618, no.
# Try: 3*phi - 2*phi = phi = 1.618, no.
# Try: 9 - 3*sqrt(5) = 9 - 6.708 = 2.292... YES! 9 - 3*sqrt(5) = 2.2917960675

val1 = 9 - 3 * np.sqrt(5)
val2 = 9 + 3 * np.sqrt(5)  # = 15.708...
print(f"   2.2918 ~= 9 - 3*sqrt(5)? {abs(val1 - 2.2917960675) < 1e-6}")
print(f"   15.7082 ~= 9 + 3*sqrt(5)? {abs(val2 - 15.7082039325) < 1e-6}")

# 5.5278640450: try 3*(3-sqrt(5)) = 9-3*sqrt(5) = 2.292, no.
# Try: 3*phi = 4.854, no.
# Try: 2*phi^2 = 5.236, no.
# Try: 3*phi + 1 = 5.854, no.
# Try: 12 - 3*sqrt(5) = 12 - 6.708 = 5.292, no.
# Try: phi^3 = phi^2*phi = 2.618*1.618 = 4.236, no.
# Try: 3*(3-phi) = 3*1.382 = 4.146, no.
# Try: 6 - sqrt(5)/2 = 6 - 1.118 = 4.882, no.
# Try: 2+sqrt(5) = 4.236, no.
# Try: 15/2 - 3*sqrt(5)/2 = 7.5 - 3.354 = 4.146, no.
# Try: 6 - phi = 6 - 1.618 = 4.382, no.
# Try: 6 - phi^(-1) = 6 - 0.618 = 5.382, no.
# Try: 4*phi = 6.472, no.
# Try: 3*phi + phi^(-1) = 4.854 + 0.618 = 5.472, no.
# Hmm. sqrt(5) = 2.2360679...
# 5.5278640450 / sqrt(5) = 2.472..., not clean
# 5.5278640450 - 3 = 2.5278... not clean
# Actually: 5.527864 = 6 - 0.472136 = 6 - (sqrt(5)-2) = 8 - sqrt(5)? 8-2.236 = 5.764, no.
# Try: 15/2 - sqrt(5) = 7.5 - 2.236 = 5.264, no.
# Try: 5 + sqrt(5)/2 = 5 + 1.118 = 6.118, no.
# Try: 3 + sqrt(5) = 5.236, no.
# Try: 18/2 - sqrt(5) = 9 - 2.236 = 6.764, no.
# Let me compute (15 - sqrt(45))/2: sqrt(45) = 3*sqrt(5) = 6.708; (15-6.708)/2 = 4.146, no.
# Try: (15-sqrt(5))/2 = (15-2.236)/2 = 6.382, no.
# Actually let me just try: 6 - sqrt(5)/...
# 5.527864045 * 2 = 11.05572809
# 5.527864045 * 3 = 16.583592135
# 14.472135955 * 3 / ... hmm
# The conjugate of 5.527864 should be 14.472136 (they sum to 20, diff is...)
# 5.527864 + 14.472136 = 20. So the sum of conjugates = 20.
# 14.472136 - 5.527864 = 8.944... = 4*sqrt(5)? 4*2.236 = 8.944 YES!
# So they are (20 +- 4*sqrt(5))/2 = 10 +- 2*sqrt(5)
val3 = 10 - 2 * np.sqrt(5)
val4 = 10 + 2 * np.sqrt(5)
print(f"   5.5279 ~= 10 - 2*sqrt(5)? {abs(val3 - 5.5278640450) < 1e-6}")
print(f"   14.4721 ~= 10 + 2*sqrt(5)? {abs(val4 - 14.4721359550) < 1e-6}")

print()
print("   EXACT Eigenvalues of 600-cell Laplacian:")
print("      0             x 1")
print("      9 - 3*sqrt(5) x 4")
print("      10 - 2*sqrt(5) x 9")
print("      9             x 16")
print("      12            x 25")
print("      14            x 36")
print("      10 + 2*sqrt(5) x 9")
print("      15            x 16")
print("      9 + 3*sqrt(5) x 4")
print()

# Check total multiplicity: 1+4+9+16+25+36+9+16+4 = 120 YES!
total_mult = 1 + 4 + 9 + 16 + 25 + 36 + 9 + 16 + 4
print(f"   Total multiplicity: {total_mult} (should be 120)")

# Compute Tr(L^+):
# Rational parts:
rational_part = Fraction(16, 9) + Fraction(25, 12) + Fraction(36, 14) + Fraction(16, 15)
print(f"   Rational part: 16/9 + 25/12 + 36/14 + 16/15")
print(f"                = {Fraction(16,9)} + {Fraction(25,12)} + {Fraction(18,7)} + {Fraction(16,15)}")
print(f"                = {rational_part}")

# Irrational part: conjugate pairs cancel!
# Pair 1: 4/(9-3*sqrt(5)) + 4/(9+3*sqrt(5))
# = 4 * [(9+3*sqrt(5)) + (9-3*sqrt(5))] / [(9-3*sqrt(5))(9+3*sqrt(5))]
# = 4 * 18 / (81 - 45)
# = 4 * 18 / 36
# = 72 / 36 = 2
pair1 = Fraction(72, 36)
print(f"   Pair 1: 4/(9-3sqrt5) + 4/(9+3sqrt5) = 72/36 = {pair1}")

# Pair 2: 9/(10-2*sqrt(5)) + 9/(10+2*sqrt(5))
# = 9 * [(10+2*sqrt(5)) + (10-2*sqrt(5))] / [(10-2*sqrt(5))(10+2*sqrt(5))]
# = 9 * 20 / (100 - 20)
# = 9 * 20 / 80
# = 180 / 80 = 9/4
pair2 = Fraction(180, 80)
print(f"   Pair 2: 9/(10-2sqrt5) + 9/(10+2sqrt5) = 180/80 = {pair2}")

total = rational_part + pair1 + pair2
print(f"   Tr(L^+) = {rational_part} + {pair1} + {pair2} = {total}")
print(f"   Tr(L^+) = {total}")
print(f"   Decimal: {float(total):.10f}")
print(f"   Numerical check: 11.7492063492")
print(f"   Numerator = {total.numerator}")
print(f"   Denominator = {total.denominator}")
print(f"   Prime? {isprime(total.numerator)}")
if not isprime(total.numerator):
    print(f"   Factorization: {factorint(total.numerator)}")
print()

# ============================================================
# 6. 120-CELL
# ============================================================
print("6. 120-CELL")
print("   This has 27 distinct eigenvalues. Let me identify them all.")
print()

# From the numerical computation:
# -0.0000000000 x 1
# 0.1458980338 x 4
# 0.3819660113 x 9
# 0.6972243623 x 16
# 1.0745770822 x 25
# 1.4818013007 x 36
# 1.7639320225 x 24
# 2.2087121525 x 16
# 2.3819660113 x 24
# 2.6180339887 x 9
# 2.8218059200 x 36
# 3.0000000000 x 40
# 3.4480705698 x 25
# 3.5857864376 x 48
# 4.0000000000 x 18
# 4.3027756377 x 16
# 4.3819660113 x 30
# 4.6180339887 x 24
# 5.0000000000 x 8
# 6.0000000000 x 8
# 6.2360679775 x 24
# 6.4142135624 x 48
# 6.4773523480 x 25
# 6.6180339887 x 30
# 6.6963927793 x 36
# 6.7912878475 x 16
# 6.8541019662 x 4

# Total multiplicity check:
mults = [1, 4, 9, 16, 25, 36, 24, 16, 24, 9, 36, 40, 25, 48, 18, 16, 30, 24, 8, 8, 24, 48, 25, 30, 36, 16, 4]
print(f"   Total multiplicity: {sum(mults)} (should be 600)")

# The eigenvalues involve sqrt(2), sqrt(5), and combinations.
# Let me identify each:
eig_vals_num = [
    0.0, 0.1458980338, 0.3819660113, 0.6972243623, 1.0745770822,
    1.4818013007, 1.7639320225, 2.2087121525, 2.3819660113, 2.6180339887,
    2.8218059200, 3.0, 3.4480705698, 3.5857864376, 4.0, 4.3027756377,
    4.3819660113, 4.6180339887, 5.0, 6.0, 6.2360679775, 6.4142135624,
    6.4773523480, 6.6180339887, 6.6963927793, 6.7912878475, 6.8541019662
]

# Key constants:
s5 = np.sqrt(5)  # 2.2360679...
s2 = np.sqrt(2)  # 1.4142135...
phi_v = (1 + s5) / 2  # 1.6180339...
iphi = (s5 - 1) / 2  # 0.6180339...

print("   Attempting eigenvalue identification:")
print()

# Let me try to identify each eigenvalue
candidates = {
    "0": 0,
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5": 5,
    "6": 6,
    "7": 7,
    "8": 8,
    "phi": phi_v,
    "1/phi": iphi,
    "phi^2": phi_v**2,
    "2*phi": 2*phi_v,
    "3*phi": 3*phi_v,
    "sqrt(5)": s5,
    "2*sqrt(5)": 2*s5,
    "3*sqrt(5)": 3*s5,
    "sqrt(2)": s2,
    "2*sqrt(2)": 2*s2,
    "3*sqrt(2)": 3*s2,
    "phi*sqrt(2)": phi_v * s2,
    "sqrt(2)/phi": s2 / phi_v,
}

# Build more candidates from combinations
extended = {}
for name1, v1 in candidates.items():
    for name2, v2 in candidates.items():
        if v1 + v2 > 0 and v1 + v2 < 8:
            extended[f"{name1}+{name2}"] = v1 + v2
        if v1 - v2 > 0 and v1 - v2 < 8:
            extended[f"{name1}-{name2}"] = v1 - v2

# Also try half-integer combinations
for n in range(1, 15):
    for d in [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30]:
        v = n / d
        if 0 < v < 8:
            extended[f"{n}/{d}"] = v
        # With sqrt(5)
        v2 = n / d + s5
        if 0 < v2 < 8:
            extended[f"{n}/{d}+sqrt5"] = v2
        v2 = n / d - s5
        if 0 < v2 < 8:
            extended[f"{n}/{d}-sqrt5"] = v2
        # With sqrt(2)
        v2 = n / d + s2
        if 0 < v2 < 8:
            extended[f"{n}/{d}+sqrt2"] = v2
        v2 = n / d - s2
        if 0 < v2 < 8:
            extended[f"{n}/{d}-sqrt2"] = v2

# Also: a + b*sqrt(5) for various a, b
for a_num in range(-6, 14):
    for a_den in [1, 2, 3, 4, 5, 6, 10, 12]:
        a = a_num / a_den
        for b_num in range(-6, 6):
            for b_den in [1, 2, 3, 4, 5, 6, 10, 12]:
                b = b_num / b_den
                v = a + b * s5
                if 0 < v < 7.5 and abs(v) > 1e-6:
                    extended[f"{a_num}/{a_den}+{b_num}/{b_den}*s5"] = v
                v = a + b * s2
                if 0 < v < 7.5 and abs(v) > 1e-6:
                    extended[f"{a_num}/{a_den}+{b_num}/{b_den}*s2"] = v

all_candidates = {**candidates, **extended}

for i, ev in enumerate(eig_vals_num):
    if abs(ev) < 1e-6:
        print(f"   lambda = 0 (x {mults[i]})")
        continue
    best_name = None
    best_err = 1e-4
    for name, val in all_candidates.items():
        err = abs(ev - val)
        if err < best_err:
            best_err = err
            best_name = name
    if best_name:
        print(f"   lambda = {ev:.10f} ~= {best_name} = {all_candidates[best_name]:.10f} (err={best_err:.2e}) x {mults[i]}")
    else:
        print(f"   lambda = {ev:.10f} -- UNIDENTIFIED x {mults[i]}")

print()
print("   === Attempting direct numerical Tr(L^+) as fraction ===")
# Using the numerical trace value: 253.8130530514
tr_num = 253.8130530514
print(f"   Numerical Tr(L^+) = {tr_num}")

# Try to identify as p/q for small q
for q in range(1, 10000):
    p = round(tr_num * q)
    if abs(p / q - tr_num) < 1e-6:
        f = Fraction(p, q)
        print(f"   Candidate: {f} = {float(f):.10f} (err={abs(float(f) - tr_num):.2e})")
        print(f"   Numerator = {f.numerator}, Denominator = {f.denominator}")
        print(f"   Prime? {isprime(f.numerator)}")
        if not isprime(f.numerator):
            print(f"   Factorization: {factorint(f.numerator)}")
        break

print()
print("=" * 70)
print()
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'Polytope':<25} {'Tr(L^+)':<20} {'Numerator':<12} {'Prime?'}")
print("-" * 70)

results = [
    ("5-cell (simplex)", "4/5", 4, False),
    ("8-cell (tesseract)", "103/24", 103, True),
    ("16-cell", "25/24", 25, False),
    ("24-cell", "?/?", "?", "?"),
    ("600-cell", "?/?", "?", "?"),
    ("120-cell", "?/?", "?", "?"),
]

# Actually compute 24-cell exactly
tr_24 = Fraction(4, 4) + Fraction(9, 8) + Fraction(8, 10) + Fraction(2, 12)
results[3] = ("24-cell", str(tr_24), tr_24.numerator, bool(isprime(tr_24.numerator)))

for name, trace, num, prime in results:
    print(f"   {name:<25} {str(trace):<20} {str(num):<12} {prime}")
