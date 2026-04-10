# The Backwards
## nos3bl33d

Everything here derives from core.md.
The machine runs forward and backward with one formula.
Running backward derives all irrationals from ^291.

---

## The Direction

Forward: 1 -> 6 -> 11 -> ... -> 291    [add p=5 each step]
Backward: 291 -> 286 -> 281 -> ... -> 1  [subtract p=5 each step]

Same formula. s(n) = p + n. s^{-1}(n) = p + n. No new machinery.

The machine is self-inverse. Running it backwards from phi^291 generates every
phi^n for all n, in decreasing order. The irrationals appear as intermediate values.

---

## Irrationals as Fractions of ^291

In ^291 units, phi^291 = 1 (the unit. the max. the reference).

Every phi^n expressed in ^291 units:

    phi^n = phi^291 / phi^(291-n)    [exact. no approximation.]

    n=290:  phi^290 = phi^291 / phi^1   = phi^291 / phi      [one step below unit]
    n=289:  phi^289 = phi^291 / phi^2   = phi^291 / (phi+1)  [two steps below]
    n=1:    phi^1   = phi^291 / phi^290                       [291 steps below unit]

Every irrational phi^n is phi^291 (known, the unit) divided by another phi power.
The ratio of two irrationals is itself a ratio within the same tower.

In ^291 units: there are no irrationals. Everything is a ratio of powers of the unit.
The irrationality lives in the base (phi itself), not in the relationships.

---

## sqrt(5) from ^291

    phi + psi = 1           [rational: the trace]
    phi - psi = sqrt(5)     [irrational: the difference]

    phi^n + psi^n = L_n     [exact integer: Lucas number]
    phi^n - psi^n = sqrt(5) * F_n    [sqrt(5) times Fibonacci]

So:

    sqrt(5) = (phi^n - psi^n) / F_n    [for every n. exact.]

At n=291:

    sqrt(5) = (phi^291 - psi^291) / F_291

phi^291 is the unit. psi^291 is sub-Planck (≈ 0 at this scale). F_291 is an exact integer.

    sqrt(5) ≈ phi^291 / F_291    [at n=291. psi^291 negligible.]

The irrational sqrt(5) is the unit divided by the 291st Fibonacci number.
Exact in the limit. Terminates.

---

## The Backward Derivation

Start at phi^291 = 1 (the unit).
Run s^{-1}(n) = n - p backwards:

    n=291:  phi^291 = 1                        [the unit]
    n=286:  phi^286 = 1 / phi^5 = psi^5 + ...  [five steps back]
    n=281:  phi^281 = 1 / phi^10              [ten steps back]
    ...
    n=1:    phi^1 = phi = 1 / phi^290          [290 steps back, 58 p-steps]

Each backward p-step divides by phi^p = phi^5. Exact.

The chain of divisions from phi^291 back to phi^1 generates every value.
No new formula. No irrationals introduced. Just division by phi^5 at each step.

The irrational phi IS the result of running the backward chain all the way to n=1.
phi = phi^291 / phi^290. One step above the bottom.

---

## Psi from the Backward Chain

psi = -1/phi. In ^291 units:

    psi^n = (-1)^n / phi^n = (-1)^n * phi^{-n}

Running the backward chain gives phi^{-n} = 1/phi^n. The psi values are just the
sign-flipped inverses of the phi values. The negative sign (from phi*psi = -1) toggles
at each step.

    backward chain of phi:   phi^291, phi^286, ..., phi^1
    corresponding psi chain: psi^291, psi^286, ..., psi^1   [same magnitudes, sign alternates]

Both chains come from the same formula s(n) = p + n run in reverse.
Both start at the unit phi^291 = 1.
The psi chain is the phi chain with sign alternation.

---

## The Tie

Every irrational in the framework:
    = phi^n for some n
    = phi^291 / phi^(291-n)    [ratio involving the unit]
    = the backward chain value at step n

Every irrational IS a position in the backward chain starting from ^291.
The backward chain uses one formula (subtract p each step) and one starting point (phi^291).
No irrationals are imported from outside. They all fall from the unit by the self-inverse step.

    sqrt(5) = (phi^291 - psi^291) / F_291    [unit minus conjugate, over integer]
    phi     = phi^291 / phi^290              [unit over one-step-back]
    psi     = -1/phi = -phi^290/phi^291      [negated inverse of phi]

All derived. All exact. All from ^291 running backwards.

    x^2 = x + 1.
