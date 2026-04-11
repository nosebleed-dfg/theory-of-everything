# The Alpha Derivation
## nos3bl33d | (DFG) DeadFoxGroup | x^2 = x + 1 | A Wake In Outerspace

Everything here derives from core.md.
One axiom. One shape. One number. That number is 1/137.036.

---

## Start With the Axiom

    x^2 = x + 1

Two roots. phi = (1+sqrt(5))/2. psi = (1-sqrt(5))/2 = -1/phi.

phi + psi = 1. phi * psi = -1. phi^2 + psi^2 = 3.

That is the whole seed. Everything below falls from these three lines.

---

## The Shape

phi gives you the pentagon. The diagonal-to-side ratio of a regular pentagon is exactly phi. No other regular polygon has this property. phi selects the pentagon and nothing else.

Three pentagons fit at a vertex. That is the only way to tile a sphere with pentagons without overlap or gap. Three pentagons at a vertex, twelve faces total. That is the dodecahedron.

    faces     F = 12     [pentagonal faces]
    edges     E = 30     [each edge shared by two faces]
    vertices  V = 20     [each vertex shared by three faces]
    V - E + F = 20 - 30 + 12 = 2 = chi    [Euler characteristic]

    face degree   p = 5    [sides per face: pentagon]
    vertex degree d = 3    [edges per vertex: trivalent]

phi forced the pentagon. The pentagon forced the dodecahedron. The dodecahedron forced all of {V, E, F, d, p, chi}. One axiom. One shape. Six integers.

Also: pi = 5 * arccos(phi/2). The pentagon gives you pi for free.

---

## The Wiring Diagram

Take the 20 vertices of the dodecahedron. Connect the ones that share an edge. You get a graph. 20 nodes, 30 edges, every node has exactly 3 neighbors (because d = 3).

This graph has a Laplacian. The Laplacian is a 20-by-20 matrix L = 3I - A, where A is the adjacency matrix (1 if connected, 0 if not) and 3I puts the degree on the diagonal.

The Laplacian asks: at each vertex, what is the difference between this vertex's value and the average of its neighbors? It measures how information flows across the dodecahedron. It is the heat equation, the diffusion equation, the spectral geometry of the shape --- all encoded in one matrix.

---

## The Spectrum

The Laplacian is a 20-by-20 symmetric matrix. It has 20 eigenvalues. One of them is always 0 (the constant vector: if every vertex has the same value, there is no flow).

The dodecahedron has the symmetry group A5 --- the group of 60 rotations that map it to itself. This symmetry forces the eigenvalues to come in packets. Eigenvalues that are related by the symmetry must be equal. So instead of 20 independent eigenvalues, you get a few distinct values with multiplicities.

Here is how you find them.

The dodecahedron is a distance-regular graph. From any vertex, the distance partition is (1, 3, 6, 6, 3, 1) --- palindromic because the dodecahedron has an antipodal map. The intersection array is {3, 2, 1, 1, 1; 1, 1, 1, 2, 3}.

From this array you build a 6-by-6 tridiagonal matrix (one row per distance class). The eigenvalues of this small matrix ARE the distinct eigenvalues of the full 20-by-20 adjacency matrix. You are reducing a 20-dimensional problem to a 6-dimensional one, using symmetry alone.

The characteristic polynomial of this tridiagonal matrix factors as:

    mu * (mu - 3) * (mu + 2) * (mu - 1) * (mu^2 - 5) = 0

Six roots:

    mu = 0, 3, -2, 1, sqrt(5), -sqrt(5)

These are adjacency eigenvalues. The Laplacian eigenvalues are lambda = 3 - mu:

    lambda = 3, 0, 5, 2, 3-sqrt(5), 3+sqrt(5)

Now look at the golden ones.

    3 - sqrt(5) = 2/phi^2     [because phi^2 = phi+1, so 2/(phi+1) = 2/phi^2]
    3 + sqrt(5) = 2*phi^2     [= 2*(phi+1) = 2phi + 2]

The axiom is inside the spectrum. sqrt(5) = 2*phi - 1 = phi - psi. The golden ratio is literally an eigenvalue of this shape.

---

## The Multiplicities

The 20 eigenvalues are 6 distinct values with multiplicities that must sum to 20. The constraints are:

    sum of multiplicities = 20
    trace of A = 0  (no self-loops)
    trace of A^2 = 2*E = 60  (counts edge pairs)

The trace condition splits into rational and irrational parts. The irrational part forces the sqrt(5) eigenvalues to have equal multiplicity. Working through:

    mu = 3:       multiplicity 1    (connected graph)
    mu = sqrt(5): multiplicity 3    
    mu = 1:       multiplicity 5    
    mu = 0:       multiplicity 4    
    mu = -sqrt(5):multiplicity 3    
    mu = -2:      multiplicity 4    

Check: 1 + 3 + 5 + 4 + 3 + 4 = 20. Good.
Trace: 3 + 3*sqrt(5) + 5 + 0 - 3*sqrt(5) - 8 = 0. Good.
Trace of A^2: 9 + 15 + 5 + 0 + 15 + 16 = 60. Good.

In Laplacian language:

    lambda = 0:            mult 1    (the constant mode)
    lambda = 3-sqrt(5):    mult 3    (the golden low mode)
    lambda = 2:            mult 5    
    lambda = 3:            mult 4    
    lambda = 3+sqrt(5):    mult 3    (the golden high mode)
    lambda = 5:            mult 4    

Five non-zero distinct eigenvalues. 19 non-zero eigenvalues total. One zero eigenvalue.

---

## The Trace

The Laplacian pseudoinverse L0^{-1} is the matrix that inverts L on everything except the zero eigenspace. Its trace is the sum of reciprocal non-zero eigenvalues, weighted by multiplicity:

    Tr(L0^{-1}) = 3/(3-sqrt(5)) + 5/2 + 4/3 + 3/(3+sqrt(5)) + 4/5

Start with the golden pair. They rationalize:

    3/(3-sqrt(5)) + 3/(3+sqrt(5))
    = 3 * [(3+sqrt(5)) + (3-sqrt(5))] / [(3)^2 - (sqrt(5))^2]
    = 3 * 6 / (9 - 5)
    = 3 * 6 / 4
    = 18/4
    = 9/2

The denominator 9 - 5 = 4 = chi^2. The numerator 6 = chi * d. Multiplied by d = 3 gives d^2/chi = 9/2. Every piece is framework data.

Now sum everything:

    9/2 + 5/2 + 4/3 + 4/5

Common denominator = lcm(2, 3, 5) = 30:

    135/30 + 75/30 + 40/30 + 24/30 = 274/30 = 137/15

There it is.

    Tr(L0^{-1}) = 137/15

The numerator is 137. The denominator is 15 = d * p = 3 * 5.

---

## Why 137

The numerator 274 = 2 * 137, and the denominator 30 = 2 * 15, so the factor of chi = 2 cancels and you get 137/15.

But where does 137 come from inside the sum?

The golden pair contributed 135/30 to the numerator. 135 = d^3 * p = 27 * 5. The cube of the dimension times the face degree. That is the big piece.

The remaining rational terms contribute 75 + 40 + 24 = 139 to the numerator. Total: 135 + 139 = 274 = 2 * 137.

After the chi = 2 cancels: 137 = the combined spectral weight of ALL eigenspaces of the dodecahedral Laplacian, measured in units of dp = 15. The golden eigenvalues supply the d^3 * p = 135 piece. The topology (chi = 2) supplies the extra 2. Together: 137 = d^3 * p + chi = 135 + 2.

This is not numerology. The eigenvalues of the Laplacian are DETERMINED by the combinatorics of the dodecahedron. The combinatorics are DETERMINED by phi. phi is DETERMINED by the axiom. The number 137 falls out of the spectrum because there is nowhere else for it to go.

---

## From 137 to 137.036

The trace gives you the integer skeleton: 137. But alpha is not an integer. It is 1/137.035999177...

The tree-level coupling is:

    V * phi^4 = 20 * phi^4

Compute phi^4 from the axiom:

    phi^2 = phi + 1
    phi^3 = phi * phi^2 = phi(phi+1) = phi^2 + phi = 2*phi + 1
    phi^4 = phi * phi^3 = phi(2*phi+1) = 2*phi^2 + phi = 2(phi+1) + phi = 3*phi + 2

    phi^4 = 3*phi + 2 = 3*(1+sqrt(5))/2 + 2 = (7 + 3*sqrt(5))/2 = 6.85410...

    V * phi^4 = 20 * 6.85410... = 137.0820393...

The integer part is 137. Same integer. The fractional part 0.082... comes from the irrational part 30*sqrt(5) = 67.082..., which is the golden contribution. The integer 70 + fractional 67.082 = 137.082.

Now the tree level overshoots. The measured value is 137.036, not 137.082. You need two corrections to bring it down.

---

## The Corrections

Two corrections. Both built from dodecahedral constants. No free parameters.

**One-loop (A1):**

    A1 = d / (chi * phi^(chi*d) * (2*pi)^d)
       = 3 / (2 * phi^6 * (2*pi)^3)

The exponents:
- chi*d = 6: the one-loop exponent. Also the number of faces of the cube.
- (2*pi)^d = (2*pi)^3: the phase-space volume of 3D momentum space.

Computing:
    phi^5 = 5*phi + 3
    phi^6 = 8*phi + 5 = 17.94427...
    (2*pi)^3 = 248.05021...
    denominator = 2 * 17.94427 * 248.05021 = 8901.44...
    A1 = 3/8901.44 = 3.36997 * 10^{-4}

This is small. It shifts the tree level down by V*phi^4 * A1 = 137.082 * 0.000337 = 0.0462. That takes 137.082 down to 137.036. Already in the right neighborhood.

**Two-loop (A2):**

    A2 = 1 / (chi * phi^(d^3))
       = 1 / (2 * phi^27)

The exponent: d^3 = 27. The cube of the dimension. phi^27 is large:

    phi^n = F_n * phi + F_{n-1}    [Fibonacci representation]
    phi^27 = 196418 * phi + 121393 = 439204.01...
    A2 = 1 / (2 * 439204.01) = 1.13842 * 10^{-6}

This is tiny. Sub-ppm. It adds back a small correction:
    V*phi^4 * A2 = 137.082 * 0.00000114 = 0.000156.

---

## Assembly

    1/alpha = V * phi^4 * (1 - A1 + A2)
            = 137.0820393 * (1 - 0.000336997 + 0.00000113842)
            = 137.0820393 * 0.999664142
            = 137.035 999 170

CODATA 2022 measurement: 137.035 999 177(21).

Deviation: -0.007. That is -0.051 ppb. That is 0.33 sigma.

Within experimental uncertainty.

---

## The Skeleton and the Flesh

The Laplacian trace gives you the skeleton: 137/15. Multiply by dp = 15 and you get 137. The integer. The frame.

The tree-level formula V*phi^4 gives you the flesh: 137.082..., including the fractional golden overshoot.

The two corrections A1 and A2 give you the skin: shaving 0.046 off and adding 0.000156 back, landing at 137.035999170.

All three layers use the same five constants: {V, d, chi, phi, pi}. All five come from the dodecahedron. The dodecahedron comes from phi. phi comes from x^2 = x + 1.

One axiom. One shape. One number.

---

## The Parameter Count

    V   = 20              from the dodecahedron
    d   = 3               vertex degree
    chi = 2               Euler characteristic
    phi = (1+sqrt(5))/2   the axiom
    pi  = 5*arccos(phi/2) pentagon geometry

All five are determined by x^2 = x + 1. No external input. No fitting. No tuning. The formula has zero free parameters.

The fine structure constant is not a free parameter of nature. It is the spectral fingerprint of the dodecahedron. The dodecahedron is the geometric shadow of the axiom.

The number 137 was never mysterious. It was always sitting in the trace of the inverse Laplacian, waiting for someone to invert the right matrix.

---

## The Two Routes

Two independent ways to get 137 from the same shape:

**Route 1 (spectral):** Sum the reciprocal eigenvalues of the dodecahedral Laplacian. Get 137/15. Multiply by dp. Get 137.

**Route 2 (combinatorial):** Multiply the vertex count by phi^4. Get 137.082. Floor it. Get 137.

Same integer. Different paths. Both forced by the dodecahedron. The spectral route sees 137 through eigenvalues. The combinatorial route sees 137 through vertex counts and golden powers.

They agree because they are reading the same structure from two sides. One side is the spectrum (how the shape vibrates). The other side is the algebra (how the shape counts). Both sides say 137.

    137 = d^3 * p + chi = 135 + 2.
    137/15 = Tr(L0^{-1}).
    137.082... = V * phi^4.
    137.036... = V * phi^4 * (1 - A1 + A2) = 1/alpha.

One axiom. One shape. One number. Done.

    x^2 = x + 1.
