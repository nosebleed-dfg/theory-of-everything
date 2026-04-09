"""
The cross-polytope Tr(L+) has numerator (d-1)/(2d) + d/(2(d-1))
= [(d-1)^2 + d^2] / [2d(d-1)]
= [2d^2 - 2d + 1] / [2d(d-1)]

So the numerator (after reduction) comes from 2d^2 - 2d + 1.
Check: when does gcd(2d^2-2d+1, 2d(d-1)) > 1?
"""

from fractions import Fraction
from sympy import isprime, factorint
from math import gcd

print("Cross-polytope: Tr(L+) = (2d^2 - 2d + 1) / (2d(d-1))")
print()
print(f"{'d':<4} {'2d^2-2d+1':<12} {'2d(d-1)':<12} {'Tr(L+)':<20} {'Num after red':<15} {'Prime?'}")
print("-" * 75)

for d in range(2, 30):
    raw_num = 2*d*d - 2*d + 1
    raw_den = 2*d*(d-1)
    tr = Fraction(raw_num, raw_den)
    num = tr.numerator
    prime = bool(isprime(num))
    g = gcd(raw_num, raw_den)
    prime_str = "YES" if prime else f"no  ({factorint(num)})"
    print(f"   {d:<4} {raw_num:<12} {raw_den:<12} {str(tr):<20} {num:<15} {prime_str}")

print()
print("The sequence 2d^2-2d+1 for d=2,3,...:")
seq = [2*d*d - 2*d + 1 for d in range(2, 30)]
print(seq)
print()
print("These are the 'centered square numbers' or 'Delannoy number d(1,n-1)'")
print()

# Check: does gcd ever reduce?
print("GCD reductions:")
for d in range(2, 100):
    raw_num = 2*d*d - 2*d + 1
    raw_den = 2*d*(d-1)
    g = gcd(raw_num, raw_den)
    if g > 1:
        print(f"   d={d}: gcd({raw_num}, {raw_den}) = {g}")

print()
print("If no output above, then gcd is always 1 and numerator is always 2d^2-2d+1.")
print()

# Interesting: d=19 gives 2*361-38+1 = 685 = 5*137
print("d=19: 2*19^2-2*19+1 = 685 = 5*137")
print(f"   137 appears as a factor!")
print()

# Also: 2d^2-2d+1 is prime for which d?
primes_at = [d for d in range(2, 200) if isprime(2*d*d - 2*d + 1)]
print(f"d where 2d^2-2d+1 is prime: {primes_at[:20]}...")
# These are d where the cross-polytope numerator is prime
