# The Machine
## nos3bl33d

Everything here derives from core.md.
One machine. One recursion. The negative sign is an instruction, not a new system.

---

## The Machine

Each move is one x and one y movement simultaneously.

    positive move:  (+x, +y)    [both go forward]
    negative move:  (+x, -y)    [x still goes forward. y subtracts.]

x always moves +1. always. every step.
y moves +1 (positive) or -1 (negative, reading from the other face).

After n steps with k negative moves:

    x total = n              [x never stops. never goes back.]
    y total = (n-k) - k = n - 2k    [ups minus downs]

The negative moves get subtracted from the system. Not stored separately.
Not sent to another machine. Subtracted. The y total accumulates the difference.

    machine state = (n, n-2k) = (x total, y total)

When k=0:     y = n.     pure phi. maximum.
When k=n/2:   y = 0.     balanced. the critical line. 1/2.
When k=n:     y = -n.    pure psi. minimum.

The machine value:

    value   = ups / downs = (n-k) / k
    inverse = downs / ups = k / (n-k)
    product = value * inverse = 1    [always]

The machine runs forward in x. The y tracks the running difference between
positive and negative steps. The y-value IS the machine's position on the phi/psi spectrum.

The embedding lifts the machine into phi-space:

    y_k = phi^{x_k}
    phi^a * phi^b = phi^{a+b}    [addition in x = multiplication in y]

Each machine step (+1 in x) becomes multiplication by phi in y-space.
The machine counts. The embedding lives.

---

## The Negative Sign

A negative is just a square read on the other side.

The axiom:

    x^2 = x + 1    [reading from the + side: phi]
    x^2 = x + 1    [reading from the - side: psi]

Same operation. Same square. Different side.

    phi^2 = phi + 1    [+ reading]
    psi^2 = psi + 1    [- reading: identical rule, mirror position]

When you see a negative in the machine output, you are not in a new system.
You are reading the square from the psi side instead of the phi side.
The operation is the same. The side flipped.

    phi * psi = -1

The -1 says: phi is being read from the psi side (or psi from the phi side).
It is a square product across the mirror. Same square. Opposite sign. That is all.

---

## Chi Is the Mirror Recursion Count

chi = 2 is not a free constant. It is the number of recursions you run on the
mirror side when you hit a negative.

    one forward step:    +1    (the normal move)
    mirror step:         chi   (recursions on the reflected side)
    return:              +1    (come back to original side)
    total depth:         d = chi + 1 = 3    [NO. see below.]

Actually:

    one forward step:    n -> n+1    (cross the mirror)
    mirror recursions:   chi = 2 steps on the y side
    total:               n -> n+3 = d steps

The cube depth d = 3 IS the cost of one negative crossing:
    1 step to the mirror + 2 steps on the mirror = 3 total = d.

---

## The n+3 Formula Is the Mirror Rule

    L_{n+3} = chi * L_{n+1} + L_n

Read this as the machine instruction for hitting a negative:

    L_{n+1}: you are at n+1. the sign flips here (odd step: phi^(n+1)*psi^(n+1) = -1).
    chi * L_{n+1}: run chi=2 recursions at the mirror point.
    + L_n: carry the previous step back.
    = L_{n+3}: you are now 3 steps ahead. the crossing is complete.

The formula is not a coincidence. It is the exact algebraic form of the rule:
    go forward one step, hit the negative, run chi mirror steps, continue.

---

## The Positive Machine, Stated Once

The machine has one mode: positive, forward, additive.

    a(n+1) = a(n) + a(n-1)    [always. phi side.]

The negative sign (phi*psi = -1) never stops the machine.
It triggers the mirror instruction:
    flip x to y, run chi=2 recursions, return to x, continue.

After the mirror crossing, the machine is d=3 steps ahead.
The next crossing is d steps later. And so on.

The machine hits a mirror crossing at every cube step:

    n=3:   phi^3 * psi^3 = -1    [first crossing]
    n=6:   phi^6 * psi^6 = +1    [not a crossing: even power, product = +1]
    n=9:   phi^9 * psi^9 = -1    [second crossing]
    ...

Wait: phi^n * psi^n = (phi*psi)^n = (-1)^n. Crossings at all ODD n. Not just cube steps.

The cube formula L_{n+3} = chi*L_{n+1} + L_n captures the NET effect of going through
d=3 steps, which includes: one even step (no crossing), one odd step (crossing, chi applied),
one even step (return). The chi=2 appears once per d-step cycle.

---

## The Mirror Side

When you flip from x to y (phi to psi):

    the psi machine is identical to the phi machine but contracting.
    the recursion is the same: a(n+1) = a(n) + a(n-1)
    the only difference: psi < 1, so each step shrinks instead of grows.

The mirror side runs the same machine. Same recursion. Same formula.
The flip changes the direction (expanding to contracting), not the rule.

After chi=2 steps on the mirror side, you flip back.
The net result: L_{n+3} = chi*L_{n+1} + L_n. The positive machine, three steps ahead.

---

## The Machine State

Full machine state at step n:

    x = phi^n    [+half, positive, expanding]
    y = psi^n    [-half, mirror, contracting]

Positive step: x_new = x + x_prev.
Mirror step (when sign flips): compute chi=2 steps on y, return to x.

The machine is always on the x side. The y side is the mirror.
The -1 in phi*psi = -1 is the toggle between them.
Every toggle costs chi=2 recursions. One toggle every d=3 steps.

---

## The Termination

The machine runs 291 forward steps. Terminals when:

    x = phi^291    [max: universe boundary]
    y = psi^291    [min: sub-Planck, invisible]

The mirror toggles occurred at every odd n: n = 1, 3, 5, ..., 291.
Total toggles: 146 (the 146 odd numbers from 1 to 291).
Each toggle: chi=2 recursions on y side.
Total mirror recursions: 146 * 2 = 292 = 291 + 1 = phi^291 unit count + the seed.

The machine terminates when the mirror side (y = psi^291) goes sub-Planck.
At that point, the toggle no longer has a real y to flip to.
The -1 instruction becomes vacuous. The machine closes.

---

## Counting the Inverse

Run the machine for n steps. Some steps are positive (phi side). Some are negative
(psi side: reading the square from the other face). Count the negatives.

    total steps:    n
    negative steps: k    [the steps that crossed to the mirror]
    positive steps: n - k

The inverse of the machine at step n:

    inverse(n) = k / n    [negatives divided by total]

This is the other band. The fraction of the machine that lived on the psi side.
You do not need to compute psi^n directly. You count.

Example: 360 steps, 30 of them negative.

    inverse = 30 / 360 = 1/12
    other band = 1/12 of the circle = 30 degrees

The phi band covers 330/360 = 11/12 of the steps.
The psi band covers 30/360 = 1/12 of the steps.
They sum to 1. Always.

---

## The Sum of Negatives IS the Inverse

Sum all negative steps up to n. That sum is the total psi-path length.
That total IS the inverse machine value at step n.

    sum(negative steps, 1..n) = psi-path length at n
                               = inverse of machine at x_n

You do not track irrationals. You count faces. The count is the inverse.

---

## The Large-n Limit

As n grows, the fraction of negative steps:

    k / n -> 1/e    [the derangement limit: !n / n! -> 1/e]

At large n:
    phi band: (1 - 1/e) of steps    [~63.2%]
    psi band: (1/e) of steps        [~36.8%]

The psi band is 1/e of the machine. The phi band is (1-1/e).
Their sum is 1. The circle is covered completely. Nothing remains.

At n=291:
    phi band: (1 - 1/e) * 291 ≈ 184 steps
    psi band: (1/e)     * 291 ≈ 107 steps
    inverse(291) = 107/291 ≈ 1/e    [the other band at universe boundary]

The inverse at the terminus is 1/e. Exact in the limit.

---

## s(n) = p + n

The step function:

    s(n) = p + n    [one step forward: adds p to n]

The inverse step function:

    s^{-1}(n) = p + n    [same formula. same rule.]

The inverse of the machine is the machine. The forward step and the backward step
are the same operation because phi and psi satisfy the same axiom:

    phi^2 = phi + 1    [forward reading]
    psi^2 = psi + 1    [backward reading: same formula, mirror side]

The step law does not distinguish direction. Both roots use x^2 = x+1.
Going forward adds p. Going backward adds p. The formula is the same.

    inverse(s(n)) = s(n) = p + n

The machine is its own inverse. The negative step is not a different rule —
it is the same rule read from the other face. The formula p+n runs in both directions.

The self-inverse property closes the framework:
you can traverse the spiral in either direction using one formula.
Forward: phi path, expanding. Backward: psi path, contracting.
Both: s(n) = p + n.

    x^2 = x + 1.
