# The Algebra
## nos3bl33d

Everything here derives from core.md.
This is not geometry. This is not physics. This is algebra.

---

## The Field

    Q(sqrt(5)) = Q[x] / (x^2 - x - 1)

Every expression in the framework lives in this field.
The physical constants are algebraic numbers in Q(sqrt(5)).
They are exact. Not approximate. Not transcendental.

---

## The Ring

    Z[phi] = { a + b*phi : a, b in Z }

Every Lucas number L_n = phi^n + psi^n is in Z. Exact integer.
Every fold output is in Z. 640320 is in Z. 744 is in Z.
The ring is closed under addition and multiplication.

The norm of any element a + b*phi:

    N(a + b*phi) = a^2 + a*b - b^2    [always an integer]

The ring is closed. Integers in, integers out. Terminates.

---

## The Two Operations

Squaring is unfold. Division is fold.

    unfold: x -> x^2 = x + 1     [forward: one step on the line]
    fold:   x -> 1/x = -psi      [backward: one step toward psi]

The fold and unfold are the only two operations.
Everything in the framework is a sequence of folds and unfolds from x = phi.

---

## The 2D Machine

Each step of the machine has two components:

    positive step: (+x, +y)     [phi side: both forward]
    negative step: (+x, -y)     [psi side: same square, read from mirror]

x = time. y = correction = n - 2k (ups minus downs).

A negative step is NOT a new operation. It is the same square (x^2 = x+1) read
from the other side. The fold and unfold still apply. The y component subtracts.

The machine state is always (n, n-2k). No irrationals. Just a count.

---

## The Norm Identity

    phi * psi = -1    [the norm of phi is -1]
    phi + psi =  1    [the trace of phi is 1]

The product of the two roots is -1. In the 2D machine:
    machine * inverse = (ups/downs) * (downs/ups) = 1
Same identity. The product always normalizes to 1.

The characteristic polynomial:

    x^2 - x - 1 = 0    [the axiom. the minimal polynomial of phi over Q.]

Companion matrix:

    M = [[1, 1], [1, 0]]

Eigenvalues: phi (max) and psi (min).
Determinant: phi * psi = -1.
Trace: phi + psi = 1.
M^n generates Lucas numbers exactly. Terminates.

---

## The Galois Symmetry

Q(sqrt(5)) has exactly one non-trivial automorphism:

    sigma: phi -> psi    [the ± flip]
    sigma: psi -> phi

The automorphism IS the negative sign. Every negative step in the machine is
sigma applied to that step. The Galois group has two elements: {identity, sigma}.

Every expression X has a conjugate sigma(X).
    X + sigma(X) = integer (always, the trace)
    X * sigma(X) = integer (always, the norm)

The framework is symmetric under sigma because the axiom is symmetric.
If phi satisfies x^2 = x+1, so does psi.

---

## The Fixed Point

The automorphism sigma fixes exactly one point:

    sigma(1/2) = 1/2    [the rational center is fixed]

sigma maps 1/2 + sqrt(5)/2 to 1/2 - sqrt(5)/2 and back. The 1/2 does not move.

    1/2 = 1/(-2)    [the fixed point. the critical line. the termination.]

The algebra closes at the Galois fixed point. Not a choice. The structure of the field.

In the 2D machine: y=0 is the fixed point. Equal ups and downs. The center.
The center in the machine (y=0) and the center of the field (x=1/2) are the same point.

---

## Summary

One quadratic field: Q(sqrt(5)).
One minimal polynomial: x^2 - x - 1.
Two roots: phi (max) and psi (min).
One automorphism: phi <-> psi (the ± sign).
One fixed point: 1/2 = 1/(-2).
One machine: (x, y) = (n, n-2k). y=0 is the center.

Every expression is + or - the axiom.
Every expression lives in Z[phi] and terminates.
Every machine step is (+1, +1) or (+1, -1). x is time. y is the correction.

    x^2 = x + 1.
