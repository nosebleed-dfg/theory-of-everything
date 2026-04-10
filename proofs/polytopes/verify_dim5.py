"""
CRITICAL: Dimension 5 shows all 3 regular polytopes having prime numerators.
Need to verify this carefully and check wider range.
"""

from fractions import Fraction
from sympy import isprime, factorint, binomial

print("DIMENSION 5 VERIFICATION:")
print()

# 5-simplex (K6)
tr_s = Fraction(5, 6)
print(f"5-simplex (K6): Tr(L+) = {tr_s}, numerator = {tr_s.numerator}, prime? {isprime(tr_s.numerator)}")

# 5-cube (Q5)
tr_h = Fraction(0)
for k in range(1, 6):
    c = int(binomial(5, k))
    term = Fraction(c, 2*k)
    tr_h += term
    print(f"   k={k}: C(5,{k})={c}, term={c}/(2*{k})={term}")
print(f"5-cube (Q5): Tr(L+) = {tr_h}, numerator = {tr_h.numerator}, prime? {isprime(tr_h.numerator)}")

# 5-cross-polytope (CP5)
raw_num = 2*25 - 10 + 1  # = 41
raw_den = 2*5*4  # = 40
tr_c = Fraction(raw_num, raw_den)
print(f"5-cross-polytope (CP5): Tr(L+) = {tr_c}, numerator = {tr_c.numerator}, prime? {isprime(tr_c.numerator)}")

print()
print("All three are prime: 5, 887, 41")
print()

# But wait: in 2D, there are infinitely many regular polygons.
# In 3D, there are exactly 5 Platonic solids.
# In 4D, there are exactly 6.
# In d >= 5, there are exactly 3.
#
# The "all prime" property in d=5 is for 3 polytopes.
# The "all prime" property in d=3 is for 5 polytopes.
# So d=3 is still special (more polytopes, all prime).
# But d=5 also achieves 100%.

# Let me check ALL dimensions from 2 to 50 systematically
print("=" * 70)
print("COMPLETE SURVEY d=2 to d=50")
print("=" * 70)
print()

results = {}

for d in range(2, 51):
    if d == 2:
        # Check first 3 regular polygons (triangle, square, pentagon)
        # But really there are infinitely many, so we can't say "all prime"
        # Skip this dimension for the "all prime" comparison
        pass
    elif d == 3:
        # 5 Platonic solids - known all prime
        results[d] = (5, 5)
    elif d == 4:
        # 6 regular polytopes - known 2 prime
        results[d] = (6, 2)
    else:
        # d >= 5: exactly 3 regular polytopes
        n = d + 1
        tr_s = Fraction(n-1, n)
        prime_s = bool(isprime(tr_s.numerator))

        tr_h = Fraction(0)
        for k in range(1, d+1):
            tr_h += Fraction(int(binomial(d, k)), 2*k)
        prime_h = bool(isprime(tr_h.numerator))

        raw_num = 2*d*d - 2*d + 1
        tr_c = Fraction(raw_num, 2*d*(d-1))
        prime_c = bool(isprime(tr_c.numerator))

        count = sum([prime_s, prime_h, prime_c])
        results[d] = (3, count)

print(f"{'Dim':<5} {'# Polytopes':<15} {'# Prime':<10} {'All Prime?'}")
print("-" * 45)

all_prime_dims = []
for d in sorted(results.keys()):
    total, primes = results[d]
    all_prime = primes == total
    if all_prime:
        all_prime_dims.append(d)
    marker = " ***" if all_prime else ""
    print(f"   {d:<5} {total:<15} {primes:<10} {'YES' + marker if all_prime else 'no'}")

print()
print(f"Dimensions where ALL regular polytopes have prime Tr(L+) numerator:")
print(f"   {all_prime_dims}")
print()

# Check which specific polytopes give primes in d=5
print("Detail for all-prime dimensions:")
for d in all_prime_dims:
    print(f"\n   Dimension {d}:")
    if d == 3:
        print("      Tetrahedron: 3, Cube: 29, Octahedron: 13, Icosahedron: 7, Dodecahedron: 137")
    else:
        n = d + 1
        tr_s = Fraction(n-1, n)
        tr_h = Fraction(0)
        for k in range(1, d+1):
            tr_h += Fraction(int(binomial(d, k)), 2*k)
        tr_c = Fraction(2*d*d-2*d+1, 2*d*(d-1))
        print(f"      Simplex K_{n}: num = {tr_s.numerator}")
        print(f"      Hypercube Q_{d}: num = {tr_h.numerator}")
        print(f"      Cross-polytope CP_{d}: num = {tr_c.numerator}")
