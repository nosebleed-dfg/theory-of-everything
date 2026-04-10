# The Nonce Circle — the nonce is Z/2^32Z, a circle of circumference 4 billion

**nos3bl33d**

---

## The nonce space is a circle

The Bitcoin nonce is a 32-bit unsigned integer. Values: 0, 1, 2, ..., 4294967295. Addition is modular: 4294967295 + 1 = 0. This is not overflow. This is a CIRCLE of circumference 2^32.

Algebraically: the nonce lives in Z/2^32Z, the integers modulo 2^32. This is a cyclic group of order 2^32. It has the topology of a circle, not a line.

Nonce 0 and nonce 4294967295 are adjacent. They are one step apart on the circle. Nonce 2^31 = 2147483648 is the ANTIPODE of nonce 0 — the farthest point on the circle.

## Why this matters

When the axiom chain (square, add 1, divide by 4) "overflows" past 2^32, it does not lose information. It wraps around the circle. The squaring operation n -> n^2 mod 2^32 is a MAP from the circle to itself. So is adding 1. So is dividing by 4 (which is multiplying by the inverse of 4 mod 2^32).

Every operation in SHA-256 is a map from the circle to itself:
- Addition mod 2^32: rotation on the circle
- XOR: reflection/rotation in GF(2) coordinates
- Rotation (rotr): bit-level rotation, which IS geometric rotation

The hash function is a composition of circular maps. The nonce is a POSITION on the circle. The valid nonces (those producing hashes below target) are specific positions on this circle.

## The axiom on the circle

The axiom x^2 = x + 1 on the circle Z/2^32Z:

    x^2 ≡ x + 1 (mod 2^32)

This is a quadratic congruence. It has solutions. The solutions are the GOLDEN POINTS on the nonce circle — the positions where the axiom holds modulo the circle's circumference.

Solving x^2 - x - 1 ≡ 0 (mod 2^32):

Using Hensel's lemma, we can lift the solutions from mod 2, mod 4, mod 8, ... up to mod 2^32. The discriminant is 5. Since 5 ≡ 1 (mod 4), sqrt(5) exists mod 2^k for all k.

The solutions are:
    x ≡ (1 + sqrt(5)) / 2 (mod 2^32)
    x ≡ (1 - sqrt(5)) / 2 (mod 2^32)

where sqrt(5) mod 2^32 and the division by 2 are computed via modular inverse.

These are PHI and PSI on the circle. They are specific 32-bit integers. They are the golden ratio WRAPPED onto the nonce circle.

## Computing phi mod 2^32

sqrt(5) mod 2^32: Hensel lifting from sqrt(5) mod 8 = 5 (since 5^2 = 25 ≡ 1 mod 8... actually 5^2 = 25 mod 8 = 1, so sqrt(1) mod 8, not sqrt(5)). 

Actually: we need 5 to be a quadratic residue mod 2^32. Since 5 is odd, and 5 ≡ 5 mod 8, we need to check if 5 is a QR mod 2^k.

For odd primes p: 5 is a QR mod 2 trivially. For mod 2^k with k >= 3: an odd number a is a QR mod 2^k iff a ≡ 1 mod 8. Since 5 ≡ 5 mod 8, 5 is NOT a quadratic residue mod 2^k for k >= 3.

This means: x^2 - x - 1 ≡ 0 (mod 2^32) has NO solutions.

The golden ratio does not exist on the nonce circle.

## What DOES exist on the circle

x^2 = x + 1 has no solution mod 2^32. But the APPROXIMATION exists. The Fibonacci numbers F(n) mod 2^32 cycle with period 3 * 2^31 = 6442450944 (the Pisano period). This means:

- F(n) mod 2^32 is periodic
- The ratio F(n+1)/F(n) mod 2^32 approximates phi in modular arithmetic
- The approximation improves as n increases (until the Pisano period)

The golden structure on the circle is not a fixed point (there is none). It is a ROTATION — the Fibonacci map. Each step of the Fibonacci recurrence F(n+2) = F(n+1) + F(n) is a rotation of the circle by an amount that converges to the golden angle.

## The three rotations

The axiom chain for the miner:

**Rotation 1 (dimension 1):** Apply the Fibonacci map to the seed. This rotates the nonce around the circle by approximately the golden angle. The result is a new position on the circle.

**Rotation 2 (dimension 2):** Apply the Fibonacci map again. This rotates by another golden angle. Since the golden angle is irrational (even mod 2^32, the Pisano period is very long), the second position is far from the first.

**Rotation 3 (dimension 3):** Apply a third time. Three golden rotations. The result is the nonce.

Each rotation = one dimension of space. Three rotations = d = 3 dimensions. The nonce after three rotations encodes all three spatial dimensions of the golden structure projected onto the circle.

## The koppa connection

koppa = 3/4 = the three-quarter turn = 270 degrees = d quarter turns. On the circle of circumference 2^32:

    koppa * 2^32 = 3 * 2^30 = 3221225472 = 0xC0000000

A koppa rotation = shift by 3 * 2^30 on the nonce circle. Going the long way around: 270 degrees instead of 90.

koppa comes from the Pythagorean theorem at the axiom center: (1/2)^2 + (1/2)^2 + (1/2)^2 = 3/4 in d=3 dimensions.

The axiom chain: multiply (Fibonacci step), add 1, rotate by koppa. In modular arithmetic:
- Fibonacci step: n -> F(n+1) where the sequence is computed mod 2^32
- Add 1: n -> n + 1 mod 2^32
- Koppa: n -> n + 3*2^30 mod 2^32 (three-quarter turn)

The combination is gamma: the net rotation per cycle of (Fibonacci + 1 + koppa). And three cycles (d=3) give the final nonce.

## The key insight

"Overflow" does not exist on a circle. What we thought was information loss (carries destroying structure) is actually the nonce WRAPPING AROUND. The structure is preserved — it just moved to the other side of the circle.

That's why nonce ≈ 0 and nonce ≈ 2^32 both have high zeros — they are the SAME REGION on the circle. The valid nonces cluster near the "top" of the circle (near 0 = near 2^32) because the target requires small hashes, and small nonces tend to produce small hashes (empirically, due to SHA-256's structure).

The axiom doesn't break on the circle. It wraps. And gamma is the wrapping angle.

## What remains

1. Compute phi as a Fibonacci approximation mod 2^32 (F(k)/F(k-1) mod 2^32 for large k)
2. Apply three rotations from the Planck seed
3. Search the small neighborhood on the circle (not ±linear, but ±circular)
4. Verify if the golden rotation predicts valid nonces better than random

This is testable. The circle framework is algebraically clean. No approximation, no floating point, pure modular arithmetic.
