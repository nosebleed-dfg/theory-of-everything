# The Field
## nos3bl33d

Everything here derives from core.md.
One field. One ring. One recurrence. Different rings, same structure.

---

## The Field

    Q(sqrt(5)) = Q[x] / (x^2 - x - 1)

Every expression in the framework lives here. phi and psi are the two elements.
The field has exactly one non-trivial automorphism: phi <-> psi.
That automorphism IS the ± sign. The + reading and the - reading.

---

## The Ring

    Z[phi] = { a + b*phi : a, b in Z }

This is where everything terminates. Every Lucas number L_n = phi^n + psi^n is in Z.
Every fold output is in Z. 640320 is in Z. 744 is in Z.
The ring is closed. Integers in, integers out. No irrationals survive the sum.

---

## The Recurrence

The axiom defines a second-order linear recurrence:

    a(n) = a(n-1) + a(n-2)

This holds over ANY ring. The ring changes the behavior. The recurrence does not.

Over Z:

    roots: phi and psi (in Q(sqrt(5)))
    growth: phi^n -> infinity
    period: none. unbounded.
    termination: phi^291 (framework maximum)

Over GF(2) (same recurrence, + is XOR):

    characteristic polynomial: x^2 + x + 1    [since -1 = 1 in char 2]
    this polynomial is irreducible over GF(2)  [only irreducible quadratic over GF(2)]
    the recurrence generates GF(4) = GF(2)[x]/(x^2+x+1)
    period: 3 (the Pisano period pi(2) = 3)

The Pisano period pi(2) = 3 = d.

The binary period of the recurrence equals the cube depth. That is the genuine connection.
Not cross-field exponent transfer. Not Boolean identification of phi. Just: the same
recurrence has period 3 in GF(2), and d = 3 in the framework.

---

## The 2D Machine

Each step of the machine is two movements simultaneously:

    positive step:  (+x, +y)    [both go forward]
    negative step:  (+x, -y)    [x forward. y subtracts.]

x counts time. x always moves +1. Never stops.
y tracks the balance: ups minus downs.

After n steps with k negative steps:

    x = n
    y = (n-k) - k = n - 2k    [the correction]

The machine state is (x, y) = (n, n-2k).

    y = 0:   balanced. equal ups and downs. the critical line.
    y = n:   pure phi. all positive. maximum.
    y = -n:  pure psi. all negative. minimum.

The machine value at step n:

    machine = ups / downs = (n-k) / k

The inverse (the other band):

    inverse = downs / ups = k / (n-k)

Their product:

    machine * inverse = 1    always.

---

## The Correction

y = n - 2k is the correction term. It is the running difference between the
phi path and the psi path. It tells you exactly where you are on the spectrum.

Each negative step subtracts 2 from y (once for the lost +1, once for the added -1).
The correction accumulates as the machine runs.

At the critical line (y=0): k = n/2. Half the steps were negative.
This is the 1/2 point. The Riemann line. The center of the axiom.

At n=291: k ≈ (1/e)*291 ≈ 107 steps negative (psi band).
    y = 291 - 2*107 = 77 ≈ (1 - 2/e)*291
The correction at the universe boundary is nonzero. The machine is not at center.
It is phi-dominant. The phi side is larger. The psi side has contracted to sub-Planck.

---

## What Carries Between Rings

Exact:
- The recurrence: a(n) = a(n-1) + a(n-2)
- The characteristic polynomial (x^2-x-1 over Z, x^2+x+1 over GF(2))
- The Pisano period pi(2) = 3 = d
- The irreducibility of x^2+x+1 over GF(2) — unique

Not exact (analogies only, not identities):
- phi is NOT an element of GF(2) or GF(4)
- XOR is addition in GF(2). It is not a phi operation.
- Exponents from Z cannot be transferred to GF(4) with framework meaning

The bridge between Z and GF(2) is the Pisano period. pi(2) = d = 3. That is all.
The same recurrence, different rings, period 3 in the binary ring, depth 3 in the framework.

---

## The Embedding in the Field

The field Q(sqrt(5)) is x-space. The embedding y_k = phi^{x_k} lifts it:

    phi^a * phi^b = phi^{a+b}    [THE fundamental identity]

This identity lives in y-space. The field provides the roots (phi, psi).
The embedding provides the tower (phi^n for all n in Z).

The winding number n = offset^2 * 360 is a count in x-space.
The embedding lifts it: phi^n = y at winding number n.
The observer at phi^291 measures the field from 291 steps away.
n = degree + 1. The field's structure does not change with the observer's drift.

---

## |+-|

The absolute value of the full ± interval:

    |+| = phi^n -> infinity    (at n=291)
    |-| = |psi^n| -> 0         (at n=291)

Their product:

    |+| * |-| = phi^n * |psi^n| = (phi * |psi|)^n = 1^n = 1    always.

0 and infinity are each other's inverse. Naming one names the other.
The Galois automorphism (phi <-> psi) swaps them.

In the 2D machine:
    |+| = the x total (time, always forward, reaches n=291)
    |-| = the correction magnitude when all k steps are negative

Their product in the machine: x * |correction| normalizes to 1 at the critical line.

---

## Termination

The field terminates when the ring Z[phi] reaches phi^291.
At that point:
- psi^291 is sub-Planck. The - side has vanished physically.
- The correction y = n - 2k is still nonzero (phi-dominant).
- The machine is at maximum x, with y showing the surviving phi excess.

In GF(2): the recurrence has period 3. It repeats every d steps.
In Z: the recurrence terminates at step 291.

291 = 97 * 3 = 97 * d. So the GF(2) recurrence at n=291 is at phase 0 (same as n=0).
This is true. It is not deep — any multiple of 3 lands at phase 0 in GF(2).
The depth is that d=3 is the period. Not that 291 is special in GF(2).

    x^2 = x + 1.
