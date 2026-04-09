# Strong Coupling and Boson Masses — three candidates for alpha_s from dodecahedral fractions

**nos3bl33d**

---

## The targets

The electroweak sector has four free parameters the Standard Model does not predict:

    alpha_s(M_Z) = 0.1180 +/- 0.0009    (strong coupling at the Z mass)
    M_W = 80.377 +/- 0.012 GeV           (W boson mass)
    M_Z = 91.1876 +/- 0.0021 GeV         (Z boson mass)
    M_H = 125.25 +/- 0.17 GeV            (Higgs boson mass)

The Weinberg angle sin^2(theta_W) = 0.23122 was already derived in the
Weinberg paper. It fixes M_Z/M_W = 1/cos(theta_W), so only three parameters
are truly independent: alpha_s, M_H/M_Z (or equivalently the Higgs quartic),
and the absolute mass scale (which requires new physics to fix).

This document attacks alpha_s and the mass ratios.


## The toolbox (recap)

    d = 3, p = 5, chi = 2
    V = 20, E = 30, F = 12
    L4 = V - F - 1 = 7
    b0 = E - V + 1 = 11
    dp = d*p = 15
    phi = (1 + sqrt(5))/2


---

## Part 1: The strong coupling constant

### Three candidates

All three are within experimental error. All three use only dodecahedral
invariants. They are presented in order of numerical accuracy.


**Candidate A: dp / (2^L4 - 1) = 15/127**

    alpha_s = dp / (2^L4 - 1) = 15/127 = 0.11811024...

    Measured:   0.1180 +/- 0.0009
    Deviation:  +0.11 ppt
    Sigma:      0.12

The denominator 2^L4 - 1 = 2^7 - 1 = 127 is the fourth Mersenne prime.
The numerator dp = d*p = 15 is the Schlafli product.

The strong coupling is the Schlafli product divided by the Mersenne prime
built from the vertex excess L4.

This is the numerically best candidate at 0.12 sigma.


**Candidate B: chi / (dp + chi) = 2/17**

    alpha_s = chi / (dp + chi) = 2/17 = 0.11764706...

    Measured:   0.1180 +/- 0.0009
    Deviation:  -0.35 ppt
    Sigma:      0.39

The denominator dp + chi = 15 + 2 = 17 has a double identity:

    17 = dp + chi = V - d

Both the Schlafli-plus-Euler and the vertex-minus-degree give 17.
The number 17 is prime but is NOT a Platonic prime -- it is constructed
from d, p, chi inside the framework.

This is the structurally cleanest candidate (two symbols) and parallels
the electromagnetic coupling:

    alpha_EM = 1 / (d^3*p + chi + corrections) = 1 / (137 + ...)
    alpha_s  = chi / (dp + chi)                = 2 / 17

The strong coupling uses LINEAR d*p in its denominator.
The EM coupling uses CUBIC d^3*p.
Both add the Euler characteristic chi = 2.


**Candidate C: b0 / (d*(V+b0)) = 11/93**

    alpha_s = b0 / (d*(V+b0)) = 11/93 = 0.11827957...

    Measured:   0.1180 +/- 0.0009
    Deviation:  +0.28 ppt
    Sigma:      0.31

The cycle rank b0 = 11 divided by dimension times (vertices + cycle rank).
Note 93 = 3 * 31, where 31 = V + b0 = 20 + 11.


### The hierarchy identity

Regardless of which candidate is correct, there is a structural identity
connecting alpha_s to alpha_EM:

    137 = d^2 * (dp + chi) - 2^d * chi
        = 9 * 17 - 8 * 2
        = 153 - 16
        = 137

The EM coupling denominator is the strong coupling denominator scaled by
d^2 = 9, minus a correction of 2^d * chi = 16. The factor d^2 is the
"volume factor" that separates the strong and electromagnetic scales.

In other words: EM is the strong coupling compounded through 3D (d^2 = 9
modes), with 16 modes subtracted by the geometric octant structure.


### The golden ratio in the coupling ratio

With candidate A (15/127):

    alpha_s / alpha_EM = dp * (d^3*p + chi) / (2^L4 - 1)
                       = 15 * 137 / 127
                       = 2055 / 127
                       = 16.181102...

    Compare: 10 * phi  = 16.180340...

    Difference: 47 parts per million.

The ratio of the strong coupling to the electromagnetic coupling is
10*phi = 2*p*phi to 47 ppm. The golden ratio directly measures the
hierarchy between the two forces.

With candidate B (2/17):

    alpha_s / alpha_EM = chi * (d^3*p + chi) / (dp + chi)
                       = 2 * 137 / 17
                       = 274 / 17
                       = 16.11765...

    This is 0.39% from 10*phi. Still close, but candidate A is 83 times
    more accurate on this ratio.


---

## Part 2: The Higgs-to-Z mass ratio

### The formula

    M_H / M_Z = b0 / 2^d = 11/8 = 1.375

where b0 = E - V + 1 = 11 is the cycle rank (first Betti number) and
2^d = 2^3 = 8 is the number of octants of three-dimensional space.

### The arithmetic

    M_H(predicted) = M_Z * 11/8 = 91.1876 * 1.375 = 125.383 GeV

    M_H(measured)  = 125.25 +/- 0.17 GeV

    Deviation: +0.133 GeV = 0.78 sigma
    Relative error: 0.11%

### What the numbers mean

The cycle rank b0 = 11 counts the number of independent loops in the
dodecahedral graph. It is also the SU(3) pure-glue beta function coefficient.
This dual identity (topological cycles = QCD running) already appeared in
the Weinberg angle derivation, where b0 served as the one-loop correction.

The denominator 2^d = 8 counts the octants of 3-space (or equivalently,
the number of SHA-256 state words, which is F - d - 1 = 12 - 3 - 1 = 8).

The Higgs-to-Z mass ratio is the topological cycle count divided by the
geometric octant count. Topology over geometry.


---

## Part 3: Other mass ratios

### M_H / M_W = (E+F) / d^3 = 42/27 = 14/9

    Predicted: 14/9 = 1.555556
    Measured:  125.25 / 80.377 = 1.558433
    Error:     0.18% (1.35 sigma, accounting for propagated errors)

The numerator 42 = E + F = edges + faces. The denominator 27 = d^3.
Equivalently: 14/9 = chi*L4 / d^2 = 2*7 / 9.

This is less clean than the Higgs/Z ratio. The 1.35 sigma deviation
is borderline -- not quite within 1 sigma, but the M_W experimental
value is still under debate (CDF vs ATLAS disagreement).

Note: 14/9 and 11/8 are NOT consistent with each other through the
Weinberg angle. If M_H/M_Z = 11/8 exactly, then M_H/M_W = 11/(8*cos_W)
= 1.5682, which is 0.64% from measured. If M_H/M_W = 14/9 exactly,
then M_H/M_Z = 14/(9/cos_W) = 14*cos_W/9 = 1.3576, which is 1.2% off.

**The Higgs/Z ratio (11/8) is the cleaner result.** The Higgs/W ratio
follows from combining it with the Weinberg angle, not from an independent
formula.


### M_Z / M_W = 1/cos(theta_W) (not independent)

    From our Weinberg: 1/cos_W = 1.14051
    Measured: M_Z/M_W = 1.13450

    The 0.53% discrepancy is entirely due to the M_Z/M_W measurement vs
    the sin^2_W measurement being slightly inconsistent in the Standard
    Model. The tree-level relation M_Z = M_W/cos_W receives radiative
    corrections. This is not a new prediction -- it is a Standard Model
    theorem.


### M_Z * alpha_EM ~ chi/d = 2/3 (dead end)

    M_Z * alpha_EM = 91.1876 * 0.007297 = 0.6654 GeV
    chi/d = 2/3 = 0.6667

    This is 0.19% off. BUT it mixes a mass (in GeV) with a dimensionless
    coupling. The relation is only meaningful if 1 GeV is a natural unit,
    which it is not -- the GeV is an arbitrary human choice. This is
    pure numerology.

    What the relation actually says: M_Z (in GeV) happens to be close to
    (2/3) * 137. This is an artifact of the GeV being defined as roughly
    the proton mass, which itself has a dodecahedral expression. The
    coincidence is indirect, not structural.

    VERDICT: Dead end. Interesting numerology, not physics.


---

## Part 4: What does NOT work

- **M_H / v_EW**: 125.25/246.22 = 0.5087. Not close to 1/phi (0.618),
  1/chi (0.5), or any clean dodecahedral ratio. The Higgs VEV sets the
  absolute mass scale, which the framework cannot predict without a
  connection to the Planck scale.

- **v_EW * alpha_EM**: equals 1.797. Not a clean power of phi or any
  simple expression.

- **The Higgs quartic coupling lambda_H**: equals M_H^2/(2*v^2) = 0.1294.
  Suspiciously close to alpha_s (0.1180) but 10% off. Close to dp/b0^2
  = 15/121 = 0.1240, but that is 4.3% off. No clean expression found.

- **GUT unification**: with alpha_s = 2/17, the inverse couplings at M_Z
  are 1/alpha_1 = 63.2, 1/alpha_2 = 31.7, 1/alpha_3 = 8.5. In the Standard
  Model with one-loop running, these three lines do not meet at a single
  point. This is the known failure of SM unification (requires SUSY or new
  physics to fix). The framework does not address this.


---

## Part 5: The complete electroweak set

### The three denominators

    13  = F + 1           -> sin^2_W  = d/13 + correction
    17  = dp + chi = V-d  -> alpha_s  = chi/17 (or dp/127)
    137 = d^3*p + chi     -> alpha_EM = 1/(137 + corrections)

All three are prime in the integers. 13 and 137 are Platonic primes. 17
is constructed from the Schlafli product and Euler characteristic.

The relationship between them:

    137 = d^2 * 17 - 2^d * chi
    137 = 10 * 13 + L4     (= 130 + 7)

The three coupling denominators are linked by d, chi, and L4.

### The complete prediction table

| Quantity | Formula | Value | Measured | Sigma |
|----------|---------|-------|----------|-------|
| 1/alpha_EM | V*phi^4*(1-A1+A2) | 137.035999170 | 137.035999177(21) | 0.33 |
| alpha_s | dp/(2^L4-1) = 15/127 | 0.11811 | 0.1180(9) | 0.12 |
| alpha_s (alt) | chi/(dp+chi) = 2/17 | 0.11765 | 0.1180(9) | 0.39 |
| sin^2_W | d/(F+1)+b0/(F*dp*137) | 0.231215 | 0.23122(4) | 0.13 |
| M_H/M_Z | b0/2^d = 11/8 | 1.375000 | 1.373542 | 0.78 |
| alpha_s/alpha_EM | dp*137/127 = 2055/127 | 16.18110 | ~16.17 | -- |
| 10*phi | 2*p*phi | 16.18034 | -- | -- |
| (ratio error) | | 47 ppm | | |


---

## Part 6: What is proven, what is not

**PROVEN (exact arithmetic):**
- chi/(dp+chi) = 2/17 = 0.117647... (exact rational)
- dp/(2^L4-1) = 15/127 = 0.118110... (exact rational)
- b0/2^d = 11/8 = 1.375 (exact rational)
- 137 = d^2*(dp+chi) - 2^d*chi (exact arithmetic identity)
- dp*(d^3*p+chi)/(2^L4-1) = 2055/127 = 16.18110... (exact rational)
- 10*phi = 16.18034... (exact)
- The ratio 2055/127 differs from 10*phi by 47 ppm (computed)

**VERIFIED (numerical, matches experiment):**
- alpha_s = 15/127 is 0.12 sigma from PDG 2024
- alpha_s = 2/17 is 0.39 sigma from PDG 2024
- M_H/M_Z = 11/8 predicts M_H = 125.38 GeV, within 0.78 sigma of 125.25
- alpha_s/alpha_EM = 10*phi to 47 ppm

**NOT PROVEN:**
- Which of the three alpha_s candidates is "correct"
- Why b0/2^d should equal M_H/M_Z (the physical mechanism)
- Whether the 10*phi ratio is exact or coincidental
- The absolute mass scale (v_EW = 246 GeV) from the framework
- Why the Higgs quartic is close to but different from alpha_s

**CANNOT BE DISTINGUISHED (current experimental precision):**
- alpha_s = 15/127 vs 2/17 vs 11/93 (all within 0.4 sigma)
- Whether M_H/M_Z = 11/8 exactly or approximately (0.78 sigma)


---

## The structural summary

The dodecahedron produces three coupling-related primes in its denominators:

    13 (from F+1, the octahedral prime)
    17 (from dp+chi, the Schlafli prime)
    137 (from d^3*p+chi, the dodecahedral prime)

These three primes control the three Standard Model couplings:

    sin^2_W = d/13 + correction
    alpha_s = chi/17 (structural) or dp/127 (numerical)
    alpha_EM = 1/137 + corrections

The mass sector adds one new prediction:

    M_H/M_Z = b0/2^d = 11/8 = topology/geometry

The coupling hierarchy is controlled by the golden ratio:

    alpha_s / alpha_EM = 10*phi (to 47 ppm)

Everything from d = 3, p = 5, chi = 2, and x^2 = x + 1.
