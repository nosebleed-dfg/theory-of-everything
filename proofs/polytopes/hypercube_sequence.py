"""
Check: is the numerator of Tr(L+) for the d-dimensional hypercube Q_d always prime?

For Q_d, the Laplacian eigenvalues are 2k for k=0,...,d with multiplicity C(d,k).

Tr(L+) = sum_{k=1}^{d} C(d,k) / (2k)
"""

from fractions import Fraction
from sympy import isprime, factorint, binomial

print("Hypercube Tr(L+) sequence:")
print()
print(f"{'d':<4} {'Vertices':<10} {'Tr(L+)':<25} {'Numerator':<15} {'Prime?'}")
print("-" * 70)

for d in range(2, 20):
    tr = Fraction(0)
    for k in range(1, d + 1):
        tr += Fraction(int(binomial(d, k)), 2 * k)

    n_vertices = 2**d
    num = tr.numerator
    prime = bool(isprime(num))
    prime_str = "YES" if prime else f"no  ({factorint(num)})"

    print(f"   {d:<4} {n_vertices:<10} {str(tr):<25} {num:<15} {prime_str}")

print()
print("Just the numerators:")
nums = []
for d in range(2, 20):
    tr = Fraction(0)
    for k in range(1, d + 1):
        tr += Fraction(int(binomial(d, k)), 2 * k)
    nums.append(tr.numerator)

print(nums)
print()
print("Prime mask:", [bool(isprime(n)) for n in nums])
