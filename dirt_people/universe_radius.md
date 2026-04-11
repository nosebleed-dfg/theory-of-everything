# The Universe Radius
## nos3bl33d | (DFG) DeadFoxGroup | x^2 = x + 1 | A Wake In Outerspace

Everything here derives from core.md.
Every number is computed. Nothing is chosen.

---

## The Axiom

Start here.

    a^2 + b^2 = c^2

A right triangle. Two legs. One hypotenuse. The square of the parts equals the square of the whole. That is the axiom. That is the only input.

---

## The Self-Similar Triangle

Ask: what right triangle looks the same after you cut a piece off?

Take a right triangle with legs a, b and hypotenuse c. Cut a off the hypotenuse. The leftover piece is c - a. If the new triangle (with sides a, something, c - a) is similar to the original, then:

    c/a = a/(c-a)

Cross multiply:

    c(c - a) = a^2
    c^2 - ac = a^2

Divide by a^2. Let x = c/a (the ratio of the long side to the short side):

    x^2 - x = 1
    x^2 = x + 1

That is the golden equation. It fell out of the axiom. We did not put it in. We asked "what triangle is self-similar?" and the axiom answered: the one whose ratio is phi.

    phi = (1 + sqrt(5)) / 2 = 1.6180339887...

The other root:

    psi = (1 - sqrt(5)) / 2 = -1/phi

    phi + psi = 1     [the trace]
    phi * psi = -1    [the product. always -1. always.]

---

## The Dodecahedron

phi tiles space. There is exactly one convex solid whose faces are all pentagons (the shape that phi builds). That solid is the regular dodecahedron.

We did not choose the dodecahedron. phi forced it. The axiom forced phi. The dodecahedron is a consequence of the right triangle.

Its numbers:

    V  = 20    vertices
    E  = 30    edges
    F  = 12    faces
    p  = 5     sides per face (pentagon)
    d  = 3     edges per vertex (valency = spatial dimension)
    chi = 2    Euler characteristic: V - E + F = 20 - 30 + 12 = 2

None of these are free. All are forced by "phi tiles a convex solid."

---

## The Exponent: Why 290

Here is where the number comes from.

Start with VE. The dodecahedron has V = 20 vertices and E = 30 edges. Their product:

    VE = 20 * 30 = 600

600 is the vertex count of the 120-cell. The 120-cell is the four-dimensional object made of 120 dodecahedra glued face-to-face. It is the 4D extension of the dodecahedron. 600 is the total combinatorial capacity of the phi-lattice when you go up one dimension.

Divide by chi = 2 (the Euler characteristic --- the two-sidedness of a closed surface):

    VE / chi = 600 / 2 = 300

This is the single-sided capacity. 300 phi-steps.

Subtract d^2 = 9. In d = 3 spatial dimensions, the metric tensor has d^2 = 9 components. These are the gravitational degrees of freedom --- they describe how space curves locally. They do not contribute to how far space extends. Remove them:

    300 - 9 = 291

291 is the structural ceiling. The maximum number of phi-steps from the axiom to the boundary of the universe.

Now subtract the seed. The first step is the axiom itself --- the +1 that launches the recurrence. The propagation is the remaining steps:

    291 - 1 = 290

290 is the exponent.

---

## The Factors of 290

    290 = 2 * 145 = 2 * 5 * 29

In dodecahedral language:

    290 = chi * p * L_7

    chi = 2     the Euler characteristic. the two legs of the axiom.
    p   = 5     the pentagon. the face of the dodecahedron.
    L_7 = 29    the seventh Lucas number. phi^7 + psi^7 = 29.

How does L_7 = 29 arise? Through cube-depth jumps:

    L_1 = 1  -->  L_4 = 7  -->  L_7 = 29

Each jump is +3 (one cube depth, d = 3). Three jumps of size d reach 29. The cube acts on itself d times.

And 291 = chi * p * L_7 + 1. The ceiling is the propagation plus the seed.

---

## The Two Legs: Why the Factor of 2

The axiom says a^2 + b^2 = c^2. Two legs. Not one. Two.

The universe has two conjugate directions. The phi direction expands. The psi direction contracts. Both exist. Both are required by the axiom. You cannot have a hypotenuse with one leg.

The factor of 2 in the formula is the Euler characteristic chi = V - E + F = 2. It counts the two-sidedness of any closed convex surface. It IS the statement "there are two legs."

Three equivalent derivations:

    Pythagorean:  two legs in a^2 + b^2 = c^2
    Euler:        V - E + F = 2 for any closed convex surface
    Algebraic:    phi^3 + psi^3 = 4 = chi^2, giving chi + 2 = chi^2, giving chi = 2

All the same statement from different angles.

---

## The Planck Length: The Minimum

The Planck length is the shortest distance at which measurement makes sense.

    l_P = 1.616255 * 10^-35 meters

In the framework, this is the number 2: the minimum state change. The axiom a^2 + b^2 = c^2 requires at least two legs of length 1. The minimum Pythagorean move that resolves both legs is 1 + 1 = 2. That is l_P. The minimum. The floor.

The universe radius is the maximum. The other end of the scale.

The number of phi-steps from minimum to maximum = 2 * phi^290.

---

## The Formula

    R = 2 * phi^290 * l_P

Every piece:

    2        = chi = the two legs of the axiom
    phi^290  = the golden ratio raised to the dodecahedral exponent
    l_P      = the minimum state change (Planck length)

Zero free parameters. Everything traces back to a^2 + b^2 = c^2.

---

## The Computation

Step by step. No shortcuts.

    phi = 1.6180339887...

    log10(phi) = 0.2089876402...

    290 * log10(phi) = 290 * 0.20899 = 60.606...

    phi^290 = 10^60.606 = 4.0403 * 10^60

    2 * phi^290 = 8.0806 * 10^60

    2 * phi^290 * l_P = 8.0806 * 10^60  *  1.616255 * 10^-35 m
                      = 1.3060 * 10^26 m

Convert to light-years:

    1 light-year = c * 1 Julian year
                 = 299,792,458 m/s  *  365.25 * 86,400 s
                 = 9.4607 * 10^15 m

    R = 1.3060 * 10^26 / 9.4607 * 10^15
      = 1.3805 * 10^10 ly
      = 13.805 Gly

---

## The Observation

    Theory:         13.805 Gly     (this derivation)
    Planck 2018:    13.787 +/- 0.020 Gly    (observed)

    Deviation:      +0.018 Gly = 0.9 sigma
    Relative error: 0.13%

The theoretical value is within one standard deviation of the measured value. At the precision of the measurement, the agreement is exact.

---

## The Dark Energy Partition

Divide the axiom phi^2 = phi + 1 by phi^2:

    1 = 1/phi + 1/phi^2

That is the axiom rewritten as a partition of unity.

    1/phi   = 0.61803...    the dark energy fraction
    1/phi^2 = 0.38197...    the matter + radiation fraction
    sum     = 1.00000       exact. the axiom.

In the phi^291 unit system:

    universe radius = phi^290 / phi^291 = 1/phi
    dark complement = phi^289 / phi^291 = 1/phi^2
    sum = 1/phi + 1/phi^2 = 1

Universe + dark = 1. No remainder. Terminates.

---

## The phi^291 Representation

    R = (2/phi) * phi^291 * l_P

Because 2 * phi^290 = 2 * phi^291 / phi = (2/phi) * phi^291.

And 1/phi = phi - 1 (from the axiom: phi^2 = phi + 1 gives phi = 1 + 1/phi gives 1/phi = phi - 1).

So R = 2(phi - 1) * phi^291 * l_P. The radius is two golden conjugates shy of the ceiling.

---

## The Cosmological Constant

    Lambda = 2 / phi^583    (in Planck units)

    583 = 2 * 291 + 1

The round trip (2 * 291) plus the seed (+1).

Check:

    Lambda * R^2 = (2 / phi^583) * (2 * phi^290)^2
                 = (2 * 4 * phi^580) / phi^583
                 = 8 / phi^3
                 = 8 / (2*phi + 1)
                 = 1.889...

Order unity. Consistent with Lambda ~ 3/R^2. Terminates.

---

## The Bit Depth

    log2(2 * phi^290) = 1 + 290 * log2(phi)
                      = 1 + 290 * 0.69424
                      = 202.33 bits

The universe requires about 202 bits to address in Planck units. That is 291 golden steps or 202 binary steps. The golden basis is more efficient than binary by the factor 1/log2(phi) = 1.44.

---

## The Scale Ladder

Every physical scale sits at a specific exponent on the phi ruler:

    exponent 0:     Planck length          (the floor)
    exponent 6:     fine structure scale   (chi * d = 2 * 3)
    exponent 183:   gravitational coupling (V*d^2 + d = 180 + 3)
    exponent 290:   universe radius        (VE/chi - d^2 - 1 = 300 - 9 - 1)
    exponent 583:   cosmological constant  (2 * 291 + 1)

Every exponent is computed from {V, E, F, d, p, chi}. None are chosen.

---

## The Derivation Chain

    a^2 + b^2 = c^2                   [the axiom]
         |
    self-similarity: x^2 = x + 1     [forces phi]
         |
    phi tiles a convex solid          [forces the dodecahedron]
         |
    dodecahedral arithmetic           [forces 290]
         |
    two legs + Planck minimum         [forces 2 and l_P]
         |
    R = 2 * phi^290 * l_P = 13.805 Gly

One axiom. One chain. One number. Terminates.

---

    x^2 = x + 1.
