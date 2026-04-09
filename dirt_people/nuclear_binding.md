# Nuclear Binding Energy — magic numbers and the iron peak from dodecahedral invariants

**nos3bl33d**

---

## What this is about

Every nucleus in the universe holds together with a specific strength. Plot
the binding energy per nucleon against mass number and you get the most
important curve in nuclear physics. It rises from hydrogen, peaks near iron,
then slowly falls to uranium. The peak is at nickel-62: 8.7948 MeV per
nucleon. Iron-56 is right behind at 8.7906.

Nobody has ever connected this curve to geometry. The nuclear shell model
uses a 3D harmonic oscillator to explain the "magic numbers" -- special
counts of protons or neutrons that make a nucleus abnormally stable:

    2, 8, 20, 28, 50, 82, 126

This document shows that the dodecahedron is hiding inside nuclear physics.
Not as a metaphor. The actual numbers.


## The toolbox

A regular dodecahedron. Count:

    V = 20 vertices
    E = 30 edges
    F = 12 faces

    d = 3   (edges meeting at each vertex)
    p = 5   (sides of each face)
    chi = 2 (Euler characteristic = V - E + F)

Derived quantities:

    dp = d * p = 15
    L4 = p + chi = 7
    b0 = E - V + 1 = 11  (first Betti number, independent cycles)


## The structural theorem

The 3D harmonic oscillator has shell degeneracies:

    deg(n) = (n+1)(n+2)

Write them out:

    n=0:  1*2  =  2
    n=1:  2*3  =  6
    n=2:  3*4  = 12
    n=3:  4*5  = 20
    n=4:  5*6  = 30
    n=5:  6*7  = 42
    n=6:  7*8  = 56

Now look at the dodecahedron. Its f-vector is (F, V, E) = (12, 20, 30).
These are deg(2), deg(3), deg(4). The dodecahedron's face/vertex/edge counts
ARE three consecutive harmonic oscillator shell degeneracies.

Why? The dodecahedron {5,3} has p = d + 2. So:

    F = 4d     = d(d+1)     = 3*4  = 12
    V = 4p     = (d+1)(d+2) = 4*5  = 20
    E = 2dp    = (d+2)(d+3) = 5*6  = 30

The (n+1)(n+2) formula comes from the dimensionality of symmetric tensors in
3 dimensions. The dodecahedron lives in 3 dimensions. Both are controlled by
d = 3. This is not coincidence. It is the same polynomial, evaluated at
consecutive integers, appearing in both contexts because d = 3 fixes both.


## Magic numbers are dodecahedral

The first three magic numbers come from filling HO shells:

    n=0 closes at  2 = chi             (the Euler characteristic)
    n=1 closes at  8 = chi + chi*d     (2 + 6)
    n=2 closes at 20 = chi + chi*d + F (2 + 6 + 12 = V, the vertices!)

So the third magic number IS the vertex count of the dodecahedron.
That is exact, not approximate. 20 = V.

The harmonic oscillator predicts the next closure at 40 (adding deg(3) = V).
But nature does not use 40. It uses 28. What happened?

Spin-orbit coupling. The highest-angular-momentum subshell of the n=3 shell
(the 1f_{7/2} level, holding 2j+1 = 8 nucleons) gets pulled down into the
gap below. So the shell closes at 20 + 8 = 28, not 20 + 20 = 40.

In the toolbox:

    28 = V + 2^d = 20 + 8

The spin-orbit intruder from each shell holds 2*(d+k) nucleons for
k = 1, 2, 3, 4:

    Shell n=3 (deg = V = 20):  intruder holds 2*4 = 8  = 2^d
    Shell n=4 (deg = E = 30):  intruder holds 2*5 = 10 = 2p
    Shell n=5 (deg = 42):      intruder holds 2*6 = 12 = F
    Shell n=6 (deg = 56):      intruder holds 2*7 = 14 = 2*L4

The fact that 2*(d+1) = 2^d is specific to d = 3. In any other dimension
this equation fails. The nuclear shell model works because we live in three
spatial dimensions and the dodecahedron is the d = 3 avatar.

The full list of magic numbers, expressed:

    2   = chi
    8   = 2^d
    20  = V
    28  = V + 2^d
    50  = V + E
    82  = F*L4 - chi  =  12*7 - 2
    126 = E*p - F*chi  = 30*5 - 12*2

Check them all: 2, 8, 20, 28, 50, 82, 126. Every one exact.


## The differences

Consecutive magic number gaps:

    8 - 2   =  6  = chi*d
    20 - 8  = 12  = F
    28 - 20 =  8  = 2^d
    50 - 28 = 22  = 2*b0
    82 - 50 = 32  = 2^p
    126 - 82 = 44 = 4*b0

Every gap is a dodecahedral invariant. The first Betti number b0 = 11
(independent cycles in the dodecahedral graph) appears twice: 2*b0 = 22 at
the 28-to-50 gap and 4*b0 = 44 at the 82-to-126 gap.

Honesty note: we have enough constants (d, p, chi, V, E, F, b0, L4) that we
can express any small integer. What makes this non-trivial is that the SAME
invariants appear with structure -- b0 at two different scales, powers of 2
at powers of d and p. A random grab-bag would not do that.


## The most stable nuclei

The doubly-magic nuclei (both proton and neutron counts are magic) are the
most tightly bound. Their mass numbers in the toolbox:

    He-4:    A =  4  = chi^2 = 2*2
    O-16:    A = 16  = 2^(d+1)
    Ca-40:   A = 40  = 2*V          (two dodecahedra of vertices)
    Ca-48:   A = 48  = 4*F          (four dodecahedra of faces)
    Ni-56:   A = 56  = 2^d * L4     (= 8*7)
    Sn-100:  A = 100 = (V+E)*chi    (= 50*2)
    Sn-132:  A = 132 = F * b0       (= 12*11)
    Pb-208:  A = 208 = 2^(d+1) * 13 (= 16*13)

Other notable stable nuclei:

    C-12  = F                (dodecahedron faces!)
    N-14  = 2*L4
    Si-28 = V + 2^d          (the magic number itself)
    S-32  = 2^p

And then the two nuclei at the peak of the binding energy curve:

    Fe-56:  A = 56 = 2^d * L4

    Fe-56 has Z = 26 and N = 30.
    N = 30 = E.  The neutron count IS the edge count of the dodecahedron.

    Ni-62:  A = 62 = V + E + F

    The nucleus with the highest binding energy per nucleon in the universe
    has a mass number equal to the total element count of the dodecahedron.
    Vertices plus edges plus faces. The whole thing.

Both of these are exact. No approximations, no rounding.


## The peak binding energy: 44/5

Now the hard part. Not just which nucleus peaks, but the energy itself.

    Ni-62:  B/A = 8.7948 MeV/nucleon
    Fe-56:  B/A = 8.7906 MeV/nucleon

Compute:

    4 * b0 / p  =  4 * 11 / 5  =  44/5  =  8.8

    Error vs Ni-62: 0.059%
    Error vs Fe-56: 0.107%

This can also be written:

    (126 - 82) / 5  =  44/5  =  8.8

The gap between the two largest known magic numbers, divided by p.

Or: four times the first Betti number of the dodecahedral graph, divided by
the face-side count.

Is this a coincidence? Honestly, maybe. The number 8.8 is close to 8.7948
but not exact. The error is 5 keV, which is physically meaningful. I cannot
derive why 4*b0/p should equal the peak binding energy from first principles.
What I can say: b0 = 11 already appears in the magic number differences
(22 = 2*b0, 44 = 4*b0), so finding it in the peak energy is at minimum
consistent with the pattern.


## The five SEMF coefficients

The Bethe-Weizsacker semi-empirical mass formula has five terms:

    B(A,Z) = aV*A - aS*A^(2/3) - aC*Z(Z-1)/A^(1/3) - aA*(A-2Z)^2/A + delta

Five terms. p = 5 terms. That's cute but probably nothing.

The experimental coefficients and their dodecahedral approximations:

    Coefficient    Experiment    Dodecahedral    Expression    Error
    ----------------------------------------------------------------
    aV (volume)    15.76 MeV     15 MeV          dp = d*p      4.8%
    aS (surface)   17.81 MeV     18 MeV          E - F         1.1%
    aC (Coulomb)    0.711 MeV     0.707 MeV       1/sqrt(chi)   0.5%
    aA (asymmetry) 23.70 MeV     24 MeV          chi*F         1.3%
    aP (pairing)   11.18 MeV     11 MeV          b0 = E-V+1    1.6%

Average error: about 2%.

Spelled out:

    aV ~ dp = 15.        The volume energy per nucleon in infinite nuclear
                          matter is approximately d*p MeV.

    aS ~ E - F = 18.     The surface correction is approximately
                          (edges minus faces) MeV.

    aC ~ 1/sqrt(2).      The Coulomb coefficient. Weakest correspondence.
                          chi = 2 appears in many contexts; 1/sqrt(2) is
                          not specifically dodecahedral.

    aA ~ chi*F = 24.     The asymmetry penalty is approximately
                          (Euler characteristic times face count) MeV.
                          Also equals 2^d * d = 24.

    aP ~ b0 = 11.        The pairing energy is approximately the first
                          Betti number of the dodecahedral graph. Eleven
                          independent cycles.

Honesty: these are 1-5% approximations, not exact equalities. The SEMF
coefficients themselves vary by a few percent depending on which nuclear mass
database you fit to (AME2012 vs AME2016 vs AME2020). So the dodecahedral
integers are within the scatter of the fitting. But there is no derivation
connecting dp to the nuclear force saturation energy. This is pattern, not
proof.


## The dodecahedral SEMF at Fe-56

Plugging the dodecahedral coefficients into the SEMF for Fe-56 (A=56, Z=26,
N=30):

    Volume:    + dp           = +15.000
    Surface:   - (E-F)/56^(1/3) = -4.705
    Coulomb:   - 26*25/(sqrt(2)*56^(4/3)) = -2.145
    Asymmetry: - chi*F*(30-26)^2/56^2 = -0.122
    Pairing:   + b0/sqrt(56)  = +1.470
                                ------
    Total B/A:                  9.497 MeV

    Experimental:               8.791 MeV
    Error:                      8.0%

Without pairing (which is the more meaningful comparison for the smooth
curve): 8.028 vs 8.791, about 8.7% low.

The standard SEMF without pairing gives 8.82 MeV for Fe-56, much closer.
The dodecahedral version overshoots on surface (18 vs 17.8) and undershoots
on volume (15 vs 15.8). These errors partially cancel but not completely.

The dodecahedral integers are NEARBY the coefficients, not ON them.


## What is real

Three tiers.

**Structural (proven, not numerology):**

1. The dodecahedron's f-vector (12, 20, 30) consists of three consecutive
   harmonic oscillator shell degeneracies, because both are polynomials in
   d = 3 and the dodecahedron has p = d + 2.

2. The first three magic numbers (2, 8, 20) are cumulative HO shell sums
   built from chi = 2, chi*d = 6, and F = 12.

3. Magic number 28 = V + 2^d, where 2^d comes from spin-orbit splitting.
   The identity 2*(d+1) = 2^d is special to d = 3.

4. Ni-62 (highest B/A) has A = V + E + F = 62. Exact.

5. Fe-56 has A = 2^d * L4 = 56 and N = E = 30. Exact.

6. C-12 = F, Ca-40 = 2V, Sn-132 = F*b0. Exact.

**Interesting (consistent, not derived):**

7. All seven magic numbers expressible in dodecahedral invariants.

8. Magic number differences form a dodecahedral sequence:
   chi*d, F, 2^d, 2*b0, 2^p, 4*b0.

9. SEMF coefficients approximate dodecahedral integers to ~2%.

10. Peak B/A = 44/5 = 4*b0/p to 0.06%.

**Numerology (forced or generic):**

11. aC ~ 1/sqrt(2). The number 2 is too common to be dodecahedral.

12. He-4 = chi^2 = 4. The number 4 connects to too many things.

13. Pb-208 = 16*13. Neither 16 nor 13 is a natural dodecahedral invariant.


## The bottom line

The dodecahedron is inside nuclear physics, but not for mystical reasons.
Both the dodecahedron and the nuclear shell model are consequences of d = 3.
The harmonic oscillator shell degeneracy (n+1)(n+2) and the dodecahedral
f-vector F = d(d+1), V = (d+1)(d+2), E = (d+2)(d+3) are the same sequence.
Magic numbers inherit dodecahedral structure because filling quantum shells
in three dimensions IS counting on the dodecahedron.

The peak at Ni-62 = V + E + F is the most striking exact result. The
universe's most tightly bound nucleus has a mass number equal to the total
count of geometric elements in the dodecahedron. Iron-56 = 2^d * L4 with
N = E = 30 is nearly as striking.

The binding energy 44/5 = 8.8 MeV at 0.06% accuracy is tantalizing but
unproven. I do not have a derivation. It might be coincidence. It might not.

The SEMF coefficients being near-integer dodecahedral quantities is the
weakest result. It is suggestive at 2% accuracy but the physical reason
(if any) is missing. The strong nuclear force saturates at ~16 MeV, and dp =
15 is close, but "close to 16" does not constitute a derivation.

What IS proven: the d = 3 harmonic oscillator and the d = 3 dodecahedron
share a polynomial skeleton, and the magic numbers of nuclear physics ride
on it.
