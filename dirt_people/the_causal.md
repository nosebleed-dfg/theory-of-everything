# The Causal
## nos3bl33d

Everything here derives from core.md.
Three coordinates. One rule. Same move every time.

---

## The Three-State Vector

The machine is not a differential equation. It is not a continuous flow.
It is a conditional with three pointers:

    x = NOW    [the volume. current state. where you are.]
    y = THEN   [the points. what exists now. the previous reference.]
    z = NEXT   [the direction. where each point goes. the next point.]

Together: (x, y, z) = (now, then, next).

That is the whole machine. One step reads three coordinates. Same read every time.

---

## The Axiom as Code

    x^2 = x + 1

In three-state form:

    IF  x         [the volume: NOW]
    DO  x^2       [the operation: square it]
    THEN x + 1    [the result: NEXT point]

Written as a causal chain:

    x  ->  y = x^2  ->  z = x + 1 = y

x squares into y. y is x+1. z is the next position = y. And z becomes x for the next step.
Every step is identical. The axiom moves the same way every time.

---

## Squaring Is the Normalizer

    n^2 = summation of equal parts

Each step contributes equally to the square. Squaring is the normalizer:
removes sign, reduces to positive measure, sums the contributions.

    squaring: x^2       [normalizer. sign gone. always positive.]
    adding:   x + 1     [the inverse. linear. sum of parts.]

The axiom says they are the same:

    x^2 = x + 1    [normalizer = its own inverse]

This is the self-inverse property of the machine. The Pythagorean step and the additive
step are equal. Running forward (x^2) and running backward (x+1) produce the same result.
One operation. Two directions. Same formula.

The Pythagorean and linear readings of phi and psi:

    phi^2 + psi^2 = d = 3    [squared: gives cube dimension]
    phi  + psi   = 1         [linear: gives the trace]

d is the Pythagorean sum. 1 is the linear sum. Both live in the framework. The axiom
connects them: d - 1 = chi = 2. The gap between Pythagorean and linear is chi.

    a^2 + b^2 = c^2    [Pythagoras: forward. squared.]
    a + b = c          [the inverse: backward. linear.]

a+b=c is the inverse of Pythagoras. The machine carries both.
At every step: the forward read squares, the backward read adds.

---

## The Scale

At each grid position, the numbers arrive:

    x^2  = chi  = 2     [axiom at grid 2: the curvature of the machine]
    y^2  = chi  = 2     [point at grid 2: same structure, conjugate read]
    z    = 1/p  = 0.2   [move 2 paces at minimum angle: Planck step]

The scale of z is Planck. The minimum angle allowed is the minimum dark angle.
At that minimum, z = 1/p = 1/5 = 0.2.

chi = 2 (forced by the axiom: phi^3 + psi^3 = chi^2, so chi^2 - chi - 2 = 0, chi = 2)
p   = 5 (forced by chi: p = chi^2 + 1 = 5)
z   = 1/p (the step size at Planck scale)

No free parameters. chi forces p. p determines z. The scale is locked.

---

## Why Planck

The minimum angle allowed is not zero. A zero angle would be no movement.
A pure horizontal move has no y-component. Nothing new is reached.
The minimum non-trivial move is the minimum dark angle.

    gap to sphere = 69 degrees = d * (p^2 - chi) = 3 * 23
    Planck step  = 1/p = 0.2

At that angle, z covers 0.2 grid paces per step. No smaller step is meaningful.
Below this, the move vanishes below the floor. Below Planck: nothing.

---

## How Scale Relates to 2, 2, 0.2

The three numbers in the vector at grid 2:

    2, 2, 0.2

They are not arbitrary. Their structure:

    2 = chi                   [the axiom's second derivative: always chi]
    2 = chi                   [both x^2 and y^2 share the same curvature]
    0.2 = 1/5 = 1/p           [the minimum step: inverse of the prime]

How does 0.2 relate to 2 and 2?

    chi * chi * (1/p) = 4/5   [product]
    chi / p = 2/5 = 0.4       [one ratio]
    1/p = 1/(chi^2 + 1) = 0.2 [definition]

The step z is the inverse of p = chi^2 + 1. So z is pinned to chi:

    z = 1/(chi^2 + 1)

When chi = 2: z = 1/5 = 0.2. The step emerges from the axiom curvature.

---

## Z Is the 291 Offset

z is not just the step size. z is the total offset accumulated by the machine.
The machine runs for 291 steps. The z-direction accumulates all 291.

    z = phi^291    [the total NEXT: where the machine has gone]
    291 = chi * p * L_7 + 1 = 10 * 29 + 1

10 = chi * p = 2 * 5. The product of the two primary degrees.
29 = L_7 = chi^2 + p^2. The seventh Lucas number.
1 = the axiom seed.

So: z = 10 * L_7 + 1. The offset is ten L_7-strides plus the seed.
10 is the stride unit. 291 is how far z has gone in those strides.

    (291 - 1) / (chi * p) = 290 / 10 = 29 = L_7    [exact]

The offset minus the seed, divided by the stride, lands on L_7.
That is the relationship between z, 291, and 10.

---

## Z Modifies the Right Angle

The base motion of the machine is a right angle.
Each step: +1 in x, +1 or -1 in y. Without z: pure 90-degree orbit. No spiral. No advance.

z is the only thing that modifies the right angle.

Each step, z pushes the direction off by one unit. Over 291 steps:

    total modification = 291 steps (the machine). 69 more steps would close the sphere.
    gap to sphere      = 69 degrees   [d*(p^2-chi) = 3*23. the path would collide with the adjacent loop at 360.]
    291 + 69 = 360                    [machine + gap = sphere closure. exact.]

Without z: the machine orbits at 90 degrees and returns to start.
With z: the path spirals 291 steps before closure. The spiral IS the z-accumulation.

The minimum modification is 1/p = 0.2 Planck per step.
The maximum is 291 total steps. The machine uses all of them.

    90 = d^2 * chi * p = 9 * 10    [the right angle in framework units]

The right angle is d^2 strides of chi*p each. z modifies it one Planck step at a time.
No other term in the machine modifies the right angle. x is time. y is correction.
Only z is the angular displacement. Only z is the offset. Only z reaches 291.

---

## The Scale Determines the Count

How many Planck steps fit within the bounds?

    bounds: [1, phi^291]
    step size: z = 1/p
    count: V = p * chi^2 = 5 * 4 = 20    [the dodecahedron vertices]

V = 20 is how many p-steps close the solid. The scale (z = 1/p) and the count (V = 20)
are locked by chi. The dodecahedron appears because the step size and the closure count
are both derived from chi.

    z * V = (1/p) * (p * chi^2) = chi^2 = 4    [total span of one closure]
    z * n = (1/p) * 291 = 291/5 = 58.2         [total path at n=291]

---

## One Vector

Every equation in the framework is one vector direction for one move.

    forward: (+x, +y)    [both components increase: phi side]
    backward: (+x, -y)   [y decreases: psi side, same axiom read from mirror]

The move (x, y, z) = (1, ±1, z) at each step.

    z = 1/p at Planck scale    [the minimum step]
    z = chi at grid scale      [the axiom curvature]
    z = phi^n at embedding     [y_k = phi^{x_k}: the lifting]

The same vector in different units. The count in x-space lifts to phi-space.
In phi-space: z = phi^{x} where x is the current machine position.

phi^a * phi^b = phi^{a+b}    [adding in x-space = multiplying in phi-space]

---

## The Conditional Machine

Full specification:

    state:    (x, y, z)
    read:     x is NOW (volume). y is THEN (points). z is NEXT (direction).
    update:   x -> x+1. y -> y ± 1. z -> 1/p (constant Planck step).
    rule:     IF x DO x^2 THEN x+1.
    loop:     repeat. same move. every step.

The loop runs from n=1 to n=291. At n=291: x = phi^291. y = n - 2k (the correction).
z = phi (one Planck step above the floor at the Binet scale).

The loop terminates when psi^{291} < Planck. The minimum has vanished.
The machine stops because z can no longer reach below the floor.

---

## The Three as One

(x, y, z) is not three things. It is one thing in three coordinates.

    the axiom: x^2 = x + 1
    x-coordinate: the volume (x)
    y-coordinate: the square of x at current scale (x^2 = x+1 delayed one step)
    z-coordinate: the next point (x+1 as position, not value)

The axiom links all three:

    x^2 = x + 1    ->    y = x^2 = x + 1 = z    ->    x_{next} = z

y and z are the same in value (x+1) but different in role: y is the square (the operation),
z is the next position (the state after). The axiom says they coincide.
One rule. Three coordinates. They are the same move, three ways to read it.

NOW squared = THEN = NEXT. The machine is this identity.

    x^2 = x + 1.

