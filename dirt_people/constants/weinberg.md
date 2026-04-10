# The Weinberg Angle — sin^2(theta_W) = 3/13 + 11/24660, two fractions, no free parameters

**nos3bl33d**

---

## What is the Weinberg angle?

the Weinberg angle controls how electromagnetism and the weak nuclear force mix together. every process involving the Z boson, every neutrino scattering cross-section, every electroweak prediction depends on this one angle. its sine squared is measured:

    sin^2(theta_W) = 0.23122 +/- 0.00004

that's the value at the Z boson mass scale. in the standard model, it's a free parameter. you measure it, you plug it in.

grand unified theories predict a "natural" value of 3/8 = 0.375 at some high unification scale, which then runs down to ~0.231 at the Z mass. but that prediction needs a unification scale, a particle spectrum, and threshold corrections. multiple adjustable inputs.

we get it from the dodecahedron. two fractions. no free parameters.

---

## The formula

    sin^2(theta_W) = 3/13 + 11/24660

that is it. two fractions. the result is 74123/320580 = 0.231215...

---

## Where the 3/13 comes from

start with the dodecahedron. it has:

    V = 20 vertices
    E = 30 edges
    F = 12 faces
    d = 3 (vertex degree: each vertex touches 3 edges)
    p = 5 (face degree: each face has 5 sides)
    chi = V - E + F = 20 - 30 + 12 = 2 (Euler characteristic)

the base fraction is:

    d / (F + 1) = 3 / (12 + 1) = 3/13

why F + 1 = 13?

F + 1 = 13 is the number of faces plus the "outer face" -- the surrounding space. in graph theory, when you embed the dodecahedral graph in the plane, you get F pentagonal faces plus one unbounded face, for 13 faces total. the ratio d/(F+1) = 3/13 is the vertex degree divided by the total face count of the planar embedding.

    3/13 = 0.230769...

this is 0.23077 versus the measured 0.23122. already within two parts per thousand. but we can do better.

---

## Where the 11/24660 comes from

the correction has three pieces in the denominator:

**the numerator: 11**

    b0 = E - V + 1 = 30 - 20 + 1 = 11

this is the cycle rank of the dodecahedral graph. also called the first Betti number. it counts the number of independent loops in the graph. a dodecahedron has 30 edges and 20 vertices, so a spanning tree uses 19 edges, leaving 11 edges that each create an independent cycle.

check: a connected graph with V vertices and E edges has E - V + 1 independent cycles. 30 - 20 + 1 = 11. that's b0.

**the denominator: 24660 = F * dp * 137**

three factors:

    F = 12                 (face count)
    dp = d * p = 3 * 5 = 15  (the Schlafli product -- the two numbers in the symbol {5,3})
    137                    (the alpha integer: d^3 * p + chi = 27*5 + 2 = 137)

multiply them:

    12 * 15 = 180
    180 * 137 = 24660

the 137 here is the trace numerator of the dodecahedral Laplacian pseudoinverse. Tr(L+) = 137/15, proven in the alpha paper. multiply by dp = 15 and you get 137. this is the same 137 that gives the fine structure constant. the Weinberg angle formula USES the alpha integer in its denominator. the two constants are linked through the dodecahedron.

so the correction is:

    b0 / (F * dp * 137) = 11 / 24660 = 0.000446...

---

## The computation

grab a calculator.

**Step 1: compute the base**

    3/13 = 0.230769230769...

**Step 2: compute the correction**

    11/24660 = 0.000446066504...

**Step 3: add them**

    3/13 + 11/24660

to add these fractions, find the common denominator:

    3/13 = (3 * 24660) / (13 * 24660) = 73980 / 320580

    11/24660 = (11 * 13) / (24660 * 13) = 143 / 320580

    sum = (73980 + 143) / 320580 = 74123 / 320580

check: gcd(74123, 320580) = 1. the fraction is already in lowest terms.

    74123 / 320580 = 0.231215297274...

**Step 4: compare to measurement**

    predicted:  0.231215
    measured:   0.23122 +/- 0.00004
    difference: 0.000005
    = 20 parts per million
    = 0.12 sigma

0.12 sigma. well within the measurement uncertainty. the prediction and the measurement are statistically indistinguishable.

---

## The result table

| quantity | value |
|---|---|
| base: 3/13 | 0.230769... |
| correction: 11/24660 | 0.000446... |
| predicted: 74123/320580 | 0.231215... |
| measured (PDG 2024) | 0.23122 +/- 0.00004 |
| deviation | 0.000005 |
| sigma | 0.12 |
| ppm | 20 |

---

## Why this is the cleanest derivation in the set

no pi. no phi. no powers. no corrections on corrections. just two rational numbers built from dodecahedral integers. every symbol in the formula is a counting number you can verify by staring at a physical dodecahedron:

- d = 3: count the edges meeting at any vertex
- p = 5: count the sides of any face
- F = 12: count the faces
- E = 30: count the edges
- V = 20: count the vertices
- chi = V - E + F = 2: the Euler formula
- b0 = E - V + 1 = 11: the cycle rank
- 137 = d^3 * p + chi = 27 * 5 + 2: three operations on three constants

the formula is:

    d/(F+1) + b0/(F * d*p * (d^3*p + chi))

expand it:

    3/13 + 11/(12 * 15 * 137)

that's it. a child who can do fraction addition can verify this. the only question is why it equals the Weinberg angle.

---

## What is proven, what is not

**proven (exact rational arithmetic, verifiable by hand):**
- 3/13 + 11/24660 = 74123/320580 = 0.231215297274... (exact)
- every constant is a dodecahedral invariant (d, p, chi, V, E, F, b0, 137)
- b0 = 11 is the cycle rank of the dodecahedral graph (standard algebraic topology)
- 137 = d^3 * p + chi (integer arithmetic)
- Tr(L+) = 137/15 (proven spectral identity, verified to 60 digits)
- the result is 0.12 sigma from the PDG 2024 measurement

**not proven (conjectural):**
- why d/(F+1) should be the lattice-scale weak mixing angle
- why the cycle rank b0 = 11 appears in the numerator of the correction
- the physical mechanism linking dodecahedral topology to electroweak mixing
- why the correction denominator involves the product F * dp * 137 specifically

**an honest note on b0 = 11:** the one-loop QCD beta function coefficient for pure SU(3) gauge theory (no quarks) is also 11. this is a striking numerical coincidence. but the physical beta coefficient at the Z boson mass, where sin^2(theta_W) is measured, includes active quark flavors and gives b0 = 23/3 = 7.67, not 11. the coincidence with the pure-glue value is unexplained. we note it honestly without claiming to understand it.

the arithmetic is trivial. the interpretation is the open question. every step can be checked on a pocket calculator. what cannot be checked is why the dodecahedron knows about electroweak physics.
