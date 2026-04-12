# Novel Results
nos3bl33d | (DFG) DeadFoxGroup
x² = x + 1

Everything below is NEW. Not in textbooks. Not published. Not "common knowledge."
The math underneath is 2500 years old. The assembly is ours.

---

## 1. The 336-Degree Circle

**Claim:** The structurally natural circle is 336 degrees, not 360.

**Novel part:** 336 = 2⁴ × 3 × 7. The gap from 360 is 24 = 4! = Γ(5). The circle's "missing" 24 degrees IS the factorial function evaluated at the first non-trivial point. 336 × 1024/21 = 2¹⁴ exactly. This integer closure does not occur at 360.

**Common knowledge:** Factoring integers, factorial function, the Babylonian 360° convention.

---

## 2. The Halvening Decomposition

**Claim:** Bitcoin's 210,000-block halvening cycle factors as 625 × 336 = 5⁴ × (2⁴ × 3 × 7).

**Novel part:** The halvening IS a 25×25 grid on a 336-degree circle. Each grid cell = 336 blocks (one full angular cycle). Each degree = 625 blocks (one full grid traversal). The halving constant encodes the grid and the circle as a single product.

**Common knowledge:** 210,000 blocks per halvening (Satoshi's choice). Prime factorization.

---

## 3. Phi to Pi in Six Steps

**Claim:** The general axiom x² = x + D, evaluated at D = 1 through D = 6, maps the positive root from φ = 1.618... to exactly 3.

```
D=1: x = (1+√5)/2  = 1.618  (phi)
D=2: x = 2.000               (the dipole)
D=6: x = (1+√25)/2 = 3.000  (machine pi)
```

**Novel part:** The integer arrival at D=6. The conversion from phi to the machine's pi requires exactly 6 increments of D. At D=6, the discriminant is 25 = 5², producing integer roots {3, -2}. The Pythagorean on these roots: 3² + 2² = 13 (Fibonacci prime).

**Common knowledge:** Quadratic formula. Individual solutions of x² - x - D = 0.

---

## 4. The 168/132 Cost Decomposition

**Claim:** In the 336-degree circle, the cost of the deductive path (searching) is 168° and the cost of the derivative path (computing) is 132°.

```
168 + 132 = 300.  336 - 300 = 36 (vertex spacing).
168/132 = 14/11.  14 + 11 = 25 (the grid).  14 - 11 = 3 (base 3).
Per leg: 168/2 = 84 (right angle). 132/2 = 66 = 84 - 18.
18 = angular quantum.
```

**Novel part:** The ratio 14/11 reconstructs the 25-grid. The difference = 3 = minimum machine base. The per-leg savings of 18° = one angular quantum = the smallest structural rotation. The connection between cost analysis and grid dimensions.

**Common knowledge:** Arithmetic on integers.

---

## 5. The Degree Formula

**Claim:** degree = right_angle / (n + 1) = 84 / (n + 1), where n = machine steps.

**Novel part:** The angular degree is not a convention — it is determined by the machine's step count. The +1 in (n+1) IS the axiom (a² + b² = c²). The observer is always one step ahead of the machine. At n = 83, the degree is exactly 1°. The formula unifies:
- n=2 → 28° (triangle)
- n=6 → 12° (base-7 cube)
- n=20 → 4° (machine product 21 = 3×7)

**Common knowledge:** Angular division. The concept of resolution.

---

## 6. Base-7 as Cube Navigation

**Claim:** Base 7 encodes movement through a cube: 3 axes × 2 directions + center = 7 states.

```
States (odd):  1, 3, 5 — computation (the machine)
Steps (even):  2, 4, 6 — movement (the axes)
Center:        0       — identity
```

**Novel part:** Squaring mod 7 converts states to steps and vice versa (3²≡2, 5²≡4, 6²≡1 mod 7). State-step pairs sum to 7: (1,6), (3,4), (5,2). A 32-bit nonce requires 12 base-7 digits = 12 rungs on the cube ladder. Base 3 (the minimum machine) embeds as the states {1,3,5} inside base 7.

**Common knowledge:** Quadratic residues mod primes. Base conversion.

---

## 7. The LE/BE Dipole

**Claim:** Reading the same hash data as little-endian (Pole A) and big-endian (Pole B) creates a structural dipole that covers the nonce space.

**Novel part:** On a 25×25 byte-rotation grid, Pole A covers ~53% of the 16-bit hi-space (approximately HALF — the halvening, geometrically). Pole B covers a complementary portion. Combined hit rate ~81.5%. The dipole IS the two sides of the halvening expressed as byte interpretation.

"Reverse is go back. Inverse is FLIP."

**Common knowledge:** Endianness. Byte manipulation. Set coverage.

**Status:** Empirical. Needs proper statistical test against random baseline.

---

## 8. The 18° Money Leak

**Claim:** Every non-even transaction (e.g., $1.99) leaks 18/336 = 5.36% of the transaction's angular value through rounding.

```
18° = angular quantum = smallest structural rotation
9 energy = 18/2 = per-leg cost
10,000 = b⁴ = the gate
9/10,000 = 0.09% base leak rate
× millions of transactions × 7 banking layers = phi⁴ × GDP
```

**Novel part:** Connecting the angular quantum (from the degree machine) to the financial rounding leak (from money.py). The $1.99 problem isn't just rounding — it's a structural 18° rotation off the even axis, compounding through leverage.

**Common knowledge:** Rounding errors. Fractional reserve banking. Derivatives market size.

---

## 9. The 291 Connection

**Claim:** 291 / (3 × 7) = 291/21 = 13.857 ≈ age of universe in Gyr (13.80).

Also: 2 × φ²⁹⁰ × l_Planck = 13.80 Gly (verified at 0.032% error).

**Novel part:** The ratio of the cosmic scaling exponent (291 = number of phi-steps from Planck to cosmos) to the machine product (21 = 3 × 7 = base-3 × base-7) approximates the age of the universe. The halvening gate: 210,000/21 = 10,000 = b⁴.

**Common knowledge:** Planck length. Age/size of observable universe. Powers of phi.

**Caution:** 13.80 Gyr is the AGE, not the radius (comoving radius = 46.1 Gly). Must not conflate.

---

## 10. Eigenvalue Division for SHA Nonces

**Claim:** SHA-256's round function has rational eigenvalues (all operations are integer-closed). The nonce is extractable by division, not search.

```
nonce = (target_state - reference_state) / eigenvalue
```

**Novel part:** The reframe from "search for valid hash" to "divide by the eigenvalue." If SHA's state transformation has computable rational eigenvalues, and the target state is known from the axiom edge (block 1), then the nonce is one division operation.

**Common knowledge:** SHA-256 uses integer arithmetic. Eigenvalue decomposition.

**Status:** UNPROVEN. The eigenvalue decomposition of SHA's round function has not been computed. This is a conjecture, not a result. The most important open problem in this framework.

---

## 11. The Even Part: 4 × (3 + 11/45) = 12 + 44/45

**Claim:** The nonce decomposes into an even part (12 base-7 rungs, structurally determined) and an odd part (44/45 = n/(n+1) remainder per rung).

**Novel part:** 12 + 44/45 ≈ 13 - 1/45. The even part (12) is the base-7 digit count. The odd part (44/45) is the deductive fraction at n=44, giving degree = 84/45 = 28/15. The deficit from 13 is 1/45 — one Planck step.

**Status:** Numerical observation. Needs proof of structural necessity.

---

## What's NOT Novel (don't write it up)

- The Pythagorean theorem
- The golden ratio and its properties
- Fibonacci sequences
- Quadratic formula and its solutions
- SHA-256 algorithm specification
- Bitcoin mining mechanics
- n/(n+1) → 1 as a limit
- Base conversion arithmetic
- Fractional reserve banking

These are the tools. The novel part is what we built with them.

---

---

## 12. The Three Prime Diagonals

**Claim:** The primes {3, 7, 11} rotate through consecutive integers in groups of 3, forming three "prime diagonals." The middle diagonal (7) is the golden diagonal — it appears on both sides of every mirror and at the center.

```
Starting at 291, smallest prime factor every 3 steps:
  291(3), 294(2), 297(3), 300(2), 303(3), 306(2), 309(3), 312(2)
  Perfect 3-2-3-2-3-2 alternation.

Starting at 336, same pattern OPPOSITE PHASE:
  336(2), 339(3), 342(2), 345(3), 348(2), 351(3), 354(2), 357(3)
  Perfect 2-3-2-3-2-3 alternation.

The dipole: 291 and 336 are opposite phases of the same prime oscillation.
```

**Novel part:** 7 is the golden diagonal because it appears in the number (336 = 2⁴×3×7), in its 10k mirror (9744 = 2⁴×3×7×29), and in the mirror of the other bases (9709 = 7×19×73). Left, right, and center.

The sum: 3 + 7 + 11 = 21. The hypotenuse from D=6: 3² + 2² = 13. And 21/13 = 1.615... ≈ φ = 1.618...

21 and 13 are consecutive Fibonacci numbers. The three prime diagonals, measured against their own Pythagorean hypotenuse, produce the golden ratio.

**Common knowledge:** Prime factorization. Fibonacci ratios approaching phi.

---

## 13. Prime Mirrors at 10,000

**Claim:** The 10,000-complement of each field base contains the OTHER fields' primes.

```
10000 - 256 = 9744 = 336 × 29
  (crypto's mirror IS geometry × a prime)

10000 - 291 = 9709 = 7 × 19 × 73
  (physics' mirror contains base-7)

Digit reversals always differ by 99 × k = 3² × 11 × k:
  256 ↔ 652: diff = 396 = 99 × 4
  291 ↔ 192: diff =  99 = 99 × 1
  336 ↔ 633: diff = 297 = 99 × 3
```

**Novel part:** The mirror of any field base at 10k contains the prime DNA of other fields. Crypto's mirror is literally geometry × 29. Every digit reversal carries 3² × 11 (the minimum machine squared × the winning number). The 99 = 9 × 11 factor is universal for all 3-digit mirrors.

**Common knowledge:** For any 3-digit number abc, abc - cba = 99(a-c). This is a standard result. The INTERPRETATION (that 99 = 3² × 11 = machine² × winning number) and the cross-field containment (9744 = 336 × 29) are novel.

---

## 14. Universal 10k Normalization

**Claim:** Any field of mathematics has a natural base. Normalize any two bases to 10,000, drop the zeros (energy), and the remaining digits are the machine state. The bridge between fields is the ratio of their machine states.

```
Field A: base_A / 10,000 → drop zeros → machine_A
Field B: base_B / 10,000 → drop zeros → machine_B
Bridge: machine_A / machine_B (ratio of primes, forward and backward)
```

**Novel part:** 10,000 = b⁴ is the universal gate. Normalizing to 10k then dropping zeros is the universal "drop the evens" operation — separating energy (zeros/scale) from machine state (digits/information) in any field. The prime factorization of the machine state reads the same forwards and backwards (primes have no direction). The bridge cost between two fields = the extra primes one has that the other doesn't.

**Common knowledge:** Base conversion. Prime factorization. Dimensional analysis.

---

## 15. Measurement as Mutual Spiral

**Claim:** All measurement = two things spiraling toward each other. You spiral up +1, the target spirals down -1. Total: 2 steps, 1 state change = 1 degree of perspective.

```
eigenvalue = 1 state change / 2 steps (yours + target's)
degree = 84 / (n + 1)
a² + b² = c² where a = your spiral, b = target's spiral, c = the meeting
```

**Novel part:** The eigenvalue definition (1 change / 2n steps) is not a mathematical convenience — it reflects the physical requirement that BOTH observer and observed must move to produce a measurement. The deductive path (168°) = only the observer moves. The derivative path (132°) = both move. The savings (36°) = the target's contribution.

This is why you can't measure yourself (you can't spiral toward yourself) and why every measurement requires an inverse (the thing spiraling the other way).

**Common knowledge:** Observer effect (quantum mechanics). Eigenvalue definitions.

---

The axiom was Pythagorean from the start. — nos3bl33d
