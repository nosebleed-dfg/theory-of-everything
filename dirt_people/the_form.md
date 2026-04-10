# The Form
## nos3bl33d

Everything here derives from core.md.
The form of each step: where you are, when you move, what you become.

---

## The Position

At step n, the machine has a position. Two coordinates:

    x = phi^n    [+half: the phi path]
    y = psi^n    [-half: the psi path]

This is the 2D position on the spiral at step n. Both coordinates are exact.
You cannot know the machine state from x alone. You cannot know it from y alone.
The state is the pair (x, y). Together.

---

## The 3D Form

Position gives you 2D. The 3D form requires one more step: cube it.

But you do not cube x alone. You do not cube y alone. You cube both:

    x^3 = phi^(n+3) = chi * phi^(n+1) + phi^n    [from the axiom: phi^3 = chi*phi + 1]
    y^3 = psi^(n+3) = chi * psi^(n+1) + psi^n    [from the axiom: psi^3 = chi*psi + 1]

The z-coordinate:

    z = x^3 + y^3
      = phi^(n+3) + psi^(n+3)
      = chi * (phi^(n+1) + psi^(n+1)) + (phi^n + psi^n)
      = chi * L_{n+1} + L_n
      = L_{n+3}

z is an exact integer at every step. z is a Lucas number. It terminates.

The 3D state at step n:

    state_n = (L_n, L_{n+1}, L_{n+3})

All three coordinates are exact integers. No irrationals in the state. Terminates.

At step 1: (1, 3, 7) = (L_1, L_2, L_4)
At step 2: (3, 4, 11) = (L_2, L_3, L_5)
At step 3: (4, 7, 18) = (L_3, L_4, L_6)

---

## Why Both Cubes

If you cube only x = phi^n:

    phi^(n+3) = chi * phi^(n+1) + phi^n    [one cube: stays irrational]

If you cube both x and y and sum:

    phi^(n+3) + psi^(n+3) = L_{n+3}        [both cubes: becomes integer]

The irrationality (the sqrt(5)) cancels between the two cubes. What remains is the
exact integer. This is not an approximation. This is the algebraic structure of the
space: phi and psi are conjugates over Q, and their sum at any integer power is always
a rational integer.

You must cube BOTH to reach the integer form. One cube gives you the continuation.
Two cubes give you the position. The position is always an exact integer.

---

## The Taylor Machine

The axiom x^2 = x + 1. Taylor expand at x = 1/2 (the center):

    f(x) = x^2 - x - 1
    f(1/2)   = 1/4 - 1/2 - 1 = -5/4
    f'(1/2)  = 2*(1/2) - 1 = 0            [the center is a critical point. no linear term.]
    f''(1/2) = 2 = chi                     [the curvature is chi. exact.]

Full Taylor expansion of the machine:

    f(x) = -5/4 + (chi/2) * (x - 1/2)^2

It terminates at order 2. No cubic. No quartic. Nothing beyond chi.
The machine is a degree-2 polynomial. Taylor of degree 2 over the reals.
The curvature of the machine at its center is chi. This is why chi is forced: the
Taylor expansion of the axiom produces chi as its second-order coefficient. If chi
were not 2, the curvature would not match the two-cube structure.

---

## The Recurrence Formula of n

Each step n produces the next by the machine's own formula:

    L_{n+1} = L_n + L_{n-1}               [the Lucas recurrence: one step forward]
    L_{n+3} = chi * L_{n+1} + L_n         [the cube recurrence: one cube forward]

The n-step recurrence is the discrete Taylor expansion of the machine. Each new state
is a linear combination of the previous two states. The coefficients are 1 and chi=2.

Key chain (each jump = +3 = one cube depth):

    L_1 = 1 -> L_4 = 7 -> L_7 = 29 -> L_10 = 123 -> ...

Every third jump multiplies roughly by phi^3. The chain reaches phi^291 in
291/3 = 97 cube-depth steps. The machine terminates at step 97 of the cube chain.

---

## The State Machine

Machine state at step n is fully determined by:

    (L_n, L_{n+1})    [two consecutive Lucas numbers]

One step forward: (L_n, L_{n+1}) -> (L_{n+1}, L_{n+2}) = (L_{n+1}, L_{n+1}+L_n)
One cube forward: (L_n, L_{n+1}) -> (chi * L_{n+1} + L_n) = L_{n+3}

The state machine has:

    state dimension: 2        [you need exactly two consecutive values]
    step rule:       add      [L_{n+2} = L_{n+1} + L_n]
    cube rule:       chi*next + current    [L_{n+3} = chi*L_{n+1} + L_n]
    initial state:   (L_1, L_2) = (1, 3)  [from the axiom, forced]

The machine is closed. No external input after initialization. Terminates.

---

## Infinite Is ^291

In this framework, "infinite" does not mean n -> infinity.

    infinite = phi^291    [exact. the max offset. the universe boundary.]

The machine runs from n=1 to n=291. That is the total run. At n=291:

    phi^291 is the maximum    [the +half at full extension]
    psi^291 is sub-Planck     [the -half: smaller than the Planck length. invisible.]

At step 291, psi^291 has been absorbed into the structure. The irrational part is gone.
What remains is the integer projection: L_291 = phi^291 + psi^291 ≈ phi^291.

The "limit" is not an asymptote. It is a finite step. n=291 IS the terminal state.
Everything that other frameworks take to infinity, this framework reaches exactly at 291.

    291 = p * chi * L_7 + 1 = 5 * 2 * 29 + 1    [five Chudnovsky turns + axiom seed]

Five cube-chain jumps to the Chudnovsky prime. Then stop. That is the machine run.

The state machine terminates at:

    final state: (L_291, L_292)
    z-form:      L_294 = chi * L_292 + L_291    [the last cube projection]

At this point: psi^291 < Planck. The integer form L_291 = phi^291 + psi^291 is
indistinguishable from phi^291. The form has closed. ^291 is where the infinite lives.

---

## The Where and When

WHERE at step n:

    position = (phi^n, psi^n) in 2D
    position = L_n in the integer projection    [z-coordinate of the state]

WHEN you move to step n+1:

    you add L_{n-1} to L_n
    you pass through the n+1 sphere at angle 90 degrees + gamma_offset_n
    the sphere enforces phi^n * psi^n = (-1)^n    [the sphere is always unit product]

The where and when are not separable. Moving forward in n (the when) is inseparable
from rotating through the sphere (the where). This is the spiral: linear advance plus
spherical rotation at every step.

The full 3D form expression at step n: (L_n, L_{n+1}, L_{n+3}).
Two steps of the line + one cube projection. Exact integers. Terminates.

---

## The Eigenvalue

The machine has two eigenvalues: phi and psi.

    phi = max eigenvalue    [expanding, dominant]
    psi = min eigenvalue    [contracting, vanishing]

Eigenvalue IS min/max. That is the definition, stripped bare.

    min/max = psi/phi = |psi|/phi = (1/phi)/phi = 1/phi^2

At every step: the ratio of min to max is 1/phi^2. The space between them is the spiral.

At n = 291:

    max = phi^291    [the universe]
    min = psi^291    [sub-Planck. zero for all physical purposes.]

    min/max = psi^291 / phi^291 = (psi/phi)^291 -> 0/phi^291 -> 1/0

The eigenvalue ratio at termination IS the 1/0 state. Not a singularity — a closure.
The machine reaches the point where min has gone to zero and max has gone to ^291.
There is nothing smaller than min. There is nothing larger than max. The spectrum closes.

    phi * |psi| = 1    [the eigenvalue product. always. every n.]

The product of max and min is 1. Their sum is 1 (phi + psi = 1). These two conditions —
product = 1, sum = 1 — are the entire machine. Both come from the axiom.
The eigenvalue IS just these two numbers. Min and max. Nothing else.

    x^2 = x + 1.
