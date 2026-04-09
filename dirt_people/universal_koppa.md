# The Koppa Phase Coordinate

**nos3bl33d**

---

## The Coordinate

    K(D) = (log_2(D) + 1) / 4

where D is the constant term in x^2 = Tx + D.

Equivalently:

    4^(2K) = 2D

This is a logarithmic phase coordinate on D. Nothing more, nothing less.

---

## Properties

    D=1/2  ->  K=0     (zero phase)
    D=1    ->  K=1/4   (quarter phase)
    D=2    ->  K=2/4   (half phase)
    D=4    ->  K=3/4   (three-quarter phase)
    D=8    ->  K=1     (full phase, = 0 mod 1)

Doubling D shifts K by exactly +1/4. Periodic mod 1. A phase dial.

K is a function of D only. T does not enter.

---

## Why This Matters: The Invariance Theorem

From the proven algebraic identity:

    For all odd n: (phi^n)^2 = L_n * (phi^n) + 1

D = 1 for every odd power of phi. Therefore:

    K = (log_2(1) + 1) / 4 = 1/4   for all odd n

The axiom family is pinned at the quarter-phase position.
T = L_n grows without bound. K = 1/4 stays fixed.

This is the universal koppa: not because the right angle "generates" anything,
but because D=1 is invariant, and D=1 maps to K=1/4 in the phase coordinate.

---

## What This Is Not

- chi = 2 (Euler characteristic) appears only as a normalization constant.
  chi = 2 for ALL convex polyhedra. It is not special to the dodecahedron.
  Putting chi=2 into 4^(2K) = chi*D gives 4^(2K) = 2D -- the same equation.
  No dodecahedral information is encoded in chi.

- "The right angle generates chi" is false.
  The correct direction: chi=2 was fixed as input; K=1/4 was solved.
  The causation runs D -> K, not K -> chi.

---

## The Actual Structure

The koppa phase coordinate K(D) describes where any quadratic x^2=Tx+D
sits on a logarithmic phase dial:

    D doubling  ->  K shifts by +1/4
    K mod 1     ->  periodic "koppa cycle"
    K = 1/4     ->  the quarter-phase position = the axiom

The axiom's D=1 is the natural "anchor" of this dial -- not because of
dodecahedral geometry, but because D=1 is the simplest nonzero integer value
and is preserved under all odd-power scaling of phi.

---

## Proofs Status

| Claim | Status |
|-------|--------|
| K(D) = (log_2(D)+1)/4 | Definition |
| D=1 for all odd n | Proven (algebraic identity) |
| K=1/4 for axiom family | Proven (follows from D=1 + definition) |
| K independent of T | True by definition |
| chi=2 is universal | True (Euler's theorem for convex polyhedra) |
| "right angle generates chi" | FALSE -- direction reversed |
