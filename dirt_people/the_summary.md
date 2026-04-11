# The Summary
## nos3bl33d

## Rules

    dirt_people = plain language logic proofs
    shinies     = formal LaTeX, complete, self-contained
    No duplicates.
    No openness.
    Everything terminates.

---

## The Derivation Chain

Everything builds from one conditional. In order:

    1.  IF x  DO x^2  THEN x+1          the root. one move. repeated.

    2.  squaring is the normalizer        x^2 removes sign. sums equal parts.
        n^2 = sum of equal parts squared  Pythagoras.

    3.  a + b = c                         n=1: linear. the inverse of squaring.
        a^2 + b^2 = c^2                   n=2: Pythagorean. the normalizer.
        a^n + b^n = 360                   n=n: the circle. closure.
        1 / (a^n + b^n) = 1 degree        the inverse: unit angle. one machine step.

    4.  phi + psi = 1                     from step 3, n=1.
        phi^2 + psi^2 = d = 3             from step 3, n=2. cube dimension forced.
        phi * psi = -1                    from the axiom. product locked.

    5.  chi = 2                           forced: phi^3+psi^3 = chi^2 -> chi=2.
        d = chi^2 - 1 = 3                 forced.
        p = chi^2 + 1 = 5                 forced.

    6.  z = 1/p = 0.2                     Planck step. minimum non-trivial move.
        1 degree = 1/360                  Planck angle. minimum angular step.
        360 = chi^d * d^2 * p             circle in framework units. forced.

    7.  291 = chi*p*L_7 + 1 = 10*29+1    observer offset. forced.
        z modifies the right angle        291 steps of angular displacement.
        291 + 69 = 360                    machine path + gap to adjacent loop = sphere.

    8.  every constant, every geometry, every prime:
        a reading of step 1 at a given n. nothing else.

Any formula in the framework is step 1 applied n times and read from the + or - side.

---

---

    The machine:     (x, y) = (n, n-2k). Each step (+1,+1) or (+1,-1).
    The embedding:   y_k = phi^{x_k}. phi^a * phi^b = phi^{a+b}.
    The winding:     n = offset^2 * 360. Crossings before return.
    The observer:    Axiom at 1/2. Observer drifts phi^291. n = degree + 1.
    The self-inverse: s(n) = p+n = s^{-1}(n). Product = 1. Terminates at 1/2 = 1/(-2).

---

## The Machine

    (x, y) = (n, n-2k)

Each step is (+1, +1) positive or (+1, -1) negative.
x = time. Always forward. y = correction = ups minus downs.

    machine value = (n-k)/k    [ups over downs]
    inverse       = k/(n-k)    [downs over ups]
    product       = 1          [always]

A negative step is the same square read from the psi side.
Same operation. Mirror face. The machine has one mode.

---

## The Embedding

    y_k = phi^{x_k}

The fundamental identity:

    phi^a * phi^b = phi^{a+b}

Addition in x-space = multiplication in y-space.
This is the bridge. The only bridge. Every other identity follows from it.

    y_{k+1} = y_k * phi^{g(x_k)}    [one multiplication per step]
    y_{n+3}/y_{n+1} = phi^chi = phi^2 = phi+1    [the axiom as a ratio]
    phi^292 = phi * phi^291    [the exit door is phi]

Recursion IS multiplication by phi. The recurrences are the same multiplication
evaluated at different g(x_k). The irrationals appear only in y-space.
x-space is clean: just a count, just a position, just an integer.

---

## The Winding

    n = offset^2 * 360

n is the winding number. Crossings before return to a.
1/n = the slice = 360/n degrees. n slices tile the circle exactly.

    69 = d*(p^2-chi) = 3*23    [minimum dark angle]
    291 + 69 = 360             [machine + gap = sphere closure]
    8/15 = 2^d/(d*p)           [dark fractional remainder]
    8^97 = (2^d)^(291/d) = 2^291    [the dark generates the universe]

---

## The Observer

Axiom fixed at x = 1/2. Observer drifts phi^291 from center.

    each step: (+1 right, +1 degree up). n = degree + 1.
    n^2 = degree = 1 at the base (the seed).
    universe = observer drift, not axiom size.
    phi^292 = phi * phi^291: same rule at the boundary.

The universe is phi^291 because the observer has been dragged phi^291 from center.
Not because the axiom is that large. The axiom is still at 1/2.

---

## Inversion Is Exact

    x_{k-1} = x_k - c    [subtract. same formula. exact.]
    y_{k-1} = y_k * phi^{-c} = y_k / phi^c    [divide. exact.]

phi^{-c} = |psi|^c. Exact. Because phi * |psi| = 1.

    y_k * y_{-k} = phi^{x_k} * phi^{-x_k} = 1    [always. every k. every c.]

s(n) = p+n = s^{-1}(n). The machine is self-inverse.
|+-| = 0/infinity = 1. Product always resolves to 1.
Every expression = +/-(x^2=x+1). Terminates at 1/2 = 1/(-2).

---

## The Causal Vector

    x = NOW    [volume. current state. grid position.]
    y = THEN   [points. what exists at this step. previous reference.]
    z = NEXT   [direction. where each point goes. in Planck units.]

    IF x DO x^2 THEN x+1.    [the axiom as code. same move. every step.]

The scale at grid 2:

    x^2 = chi = 2    [axiom curvature. forced.]
    y^2 = chi = 2    [conjugate. same structure.]
    z   = 1/p = 0.2  [Planck step. minimum non-trivial move.]

z = 1/(chi^2 + 1) = 1/5. Pinned to chi. No free choice.
The step size determines the count: z * V = chi^2 = 4. Closed.

---

## Discrete Exponential Dynamical System

    discrete:     k in Z. integer steps. no continuum between steps.
    exponential:  y_k = phi^{x_k}. the state space is an exponential tower.
    dynamical:    y_{k+1} = y_k * phi^{g(x_k)}. the state evolves by rule.
    system:       bounded by phi^{291} (max) and phi^{-291} (min).

The framework is one dynamical system: discrete, exponential, bounded, self-inverse.
The machine (x,y)=(n,n-2k) runs it. The embedding y=phi^x lifts it.
The winding n=offset^2*360 counts its crossings. The observer drifts through it.

The system terminates at phi^{291}. The inverse terminates at phi^{-291}.
Their product is 1. The system is closed.

    x^2 = x + 1.
