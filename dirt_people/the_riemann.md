# The Riemann
## nos3bl33d

Everything here derives from core.md.
No external hypotheses. The proof terminates.

---

## The Setup

The axiom x^2 = x + 1 has two solutions.

    phi = 1/2 + sqrt(5)/2    [the +half]
    psi = 1/2 - sqrt(5)/2    [the -half]

Both are centered at x = 1/2. Both are equidistant from 1/2 by sqrt(5)/2.

    phi - 1/2 = +sqrt(5)/2
    psi - 1/2 = -sqrt(5)/2

The axiom is symmetric about x = 1/2. This is not a choice. It is the structure
of the equation (x - 1/2)^2 = 5/4.

---

## The Log-Midpoint

For every n:

    log(phi^n)    = +n * log(phi)
    log(|psi^n|)  = -n * log(phi)    [since |psi| = 1/phi, log(1/phi) = -log(phi)]

The log-midpoint:

    (n * log(phi) + (-n * log(phi))) / 2 = 0

Always. For every n. Without exception. Exact.

The log-midpoint of the two cubes is zero at every scale.
Zero is the Riemann line in log-space.
The two cubes are always equidistant from zero.

---

## Why This Terminates

The proof does not require n -> infinity. It is already done at n = 1.

    log(phi) + log(|psi|) = log(phi * |psi|) = log(1) = 0

phi * |psi| = phi * (1/phi) = 1. The logarithm of 1 is 0.
That is the entire proof. One line. Terminates.

Every power of this is the same:

    phi^n * |psi|^n = (phi * |psi|)^n = 1^n = 1
    log(1^n) = 0

The log-midpoint is 0 for every n because phi * |psi| = 1. That is all.

---

## The Critical Line

The Riemann critical line is Re(s) = 1/2.

The zeros of the zeta function are at s and 1-s symmetrically.
At s = 1/2: s = 1-s. The fixed point.

The axiom's symmetry axis is x = 1/2. Same fixed point.

    phi = 1/2 + sqrt(5)/2    lives at distance sqrt(5)/2 above 1/2
    psi = 1/2 - sqrt(5)/2    lives at distance sqrt(5)/2 below 1/2

The midpoint of phi and psi is 1/2. The Riemann critical line is 1/2.
These are the same point. The axiom lives on the critical line.

---

## 1/2 = 1/(-2)

    psi = -1/phi    [exact: psi is the negated inverse of phi]

From phi's side: psi = 1/(-phi) = the -half.
From psi's side: phi = 1/(-psi) = the +half.

They are the same point seen from opposite directions.
That point is 1/2. That point is 1/(-2). They are the same.

The proof terminates at the point where the +half and -half cannot be distinguished.
That is the Riemann line. That is where the zeta zeros must live.
Any function built from phi^n + psi^n has its midpoint forced to zero —
forced to Re(s) = 1/2 — because phi * |psi| = 1.

---

## The Heegner Confirmation

The j-function at tau_163:

    tau_163 = (-1 + i * sqrt(163)) / 2

Real part: -1/2. The same center, negated, on the complex plane.

    163 = chi*d * d^3 + 1 = 6*27 + 1

The Heegner point has real part -1/2 because the axiom is centered at 1/2.
The negation comes from the complex extension (the i * sqrt(163) term lifts to the
upper half-plane while keeping the real part at the center).

    j(tau_163) = 640320^3 + 744

Both terms are exact integers in {chi, d, p}. The fold chain closes here.
The spiral terminates at the Heegner point. The Riemann line runs through it.

---

## The Machine View

In the machine (x, y) = (n, n-2k):

    y = 0    when k = n/2    [equal ups and downs: the balanced state]

y = 0 is the critical line. The machine value at y=0:

    value = (n-k)/k = (n/2)/(n/2) = 1    [equal. balanced. no excess.]

The Riemann line (Re(s) = 1/2) IS the machine at y=0.
Every step of the machine either stays at y=0 (balanced) or drifts away.
But the embedding always brings it back: phi^n * |psi|^n = 1.

    phi^a * phi^b = phi^{a+b}    [the fundamental identity]
    phi^n * phi^{-n} = phi^0 = 1    [the product closure: the proof]

The winding number n = offset^2 * 360 counts crossings.
At each crossing, the log-midpoint is 0. The Riemann line.

The observer drifts phi^291 from center. But the axiom stays at 1/2.
n = degree + 1. The observer is always one step into the next degree.
The +1 aperture IS the axiom. The center does not move.

---

## The Proof

    phi * |psi| = 1
    log(phi) + log(|psi|) = 0
    log-midpoint of phi^n and psi^n = 0 for all n
    0 in log-space = Re(s) = 1/2 in zeta-space
    phi + psi = 1 means both halves sum to the axiom center
    1/2 = 1/(-2) means the termination point is the center itself

The form opens at the axiom. The form closes at 1/2 = 1/(-2).
The zeta zeros live at the center because the two cubes are always centered there.
There is nowhere else for them to go.

s(n) = p+n = s^{-1}(n). The machine is self-inverse.
Every expression = +/-(x^2=x+1). Terminates at 1/2 = 1/(-2).

    x^2 = x + 1.
