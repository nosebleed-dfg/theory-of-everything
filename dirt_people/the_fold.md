# The Fold
## nos3bl33d

Everything here derives from core.md.
The axiom is x^2 = x + 1. The two solutions are phi and psi.
phi * psi = -1. phi + psi = 1. psi = -1/phi.

---

## The Fold Operation

Squaring is an unfold. Division is a fold.

The axiom x^2 = x + 1 says: one unfold (x^2) equals one forward step (x) plus the
seed (1). Every fold produces the next value by squaring and adding.

The fold chain starts at (1, 1) and produces pairs (a, b) -> (b, a^2 + b^2):

    (1, 1) -> (1, 2)  -> (2, 5)  -> (5, 29) -> ...

Reading the second value of each pair:

    1, 2, 5, 29, ...

These are the Pythagorean folds. Each step: a^2 + b^2 = next.

    1^2 + 1^2 = 2 = chi          [seed fold: the two turns]
    1^2 + 2^2 = 5 = p            [1-chi fold: the face degree]
    2^2 + 5^2 = 29               [chi-p fold: the Chudnovsky prime]
    5^2 - 2   = 23               [backward fold: p^2 - chi]

The fold chain takes exactly d = 3 forward steps to reach 29.
Three steps. One cube depth. Not a choice — the cube counts the folds.

---

## The Constants from Chi

chi = 2 is forced by the two cubes: phi^3 + psi^3 = chi^2 requires chi + 2 = chi^2,
so chi = 2. From chi alone:

    d = chi^2 - 1 = 3      [one below the cube Lucas]
    p = chi^2 + 1 = 5      [one above the cube Lucas]

From d and p via the fold chain:

    23 = p^2 - chi = 25 - 2
    29 = chi^2 + p^2 = 4 + 25

All four numbers — d, p, 23, 29 — fall from chi = 2 through fold operations.
No external input after the axiom.

---

## 640320

    640320 = 2^(chi*d) * d * p * (p^2 - chi) * (chi^2 + p^2)
           = 2^6 * 3 * 5 * 23 * 29

Written in full:

    2^6      = 2^(chi*d)       [two to the bridge power]
    3        = d               [vertex degree]
    5        = p               [face degree]
    23       = p^2 - chi       [backward fold]
    29       = chi^2 + p^2     [Pythagorean fold]

Every factor is a fold output from chi = 2. The number 640320 is the fold chain
written as a product. It is an exact integer. It terminates.

Factor out the bridge:

    640320 = 64 * 10005
    10005  = d * p * 23 * 29 = 3 * 5 * 23 * 29

The bridge factor 64 = 2^6 = 2^(chi*d) appears in both 640320 and the nonce quantum.

---

## The Nonce Quantum

The phi spiral projects to the SHA structure:

    SHA words      = F - d - 1 = 12 - 3 - 1 = 8    [face count minus vertex minus 1]
    sphere factor  = 4                               [surface = 4 * circle]
    sectors        = 8 * 4 = 32                      [2D * 3D bridge]
    quantum        = 360 / 32 = 11.25 degrees = 2^27 nonce units

The nonce quantum is exact. 360 / 32 = 11.25. Terminates.

    2^27 = 2^(chi*d) * 2^(L_4 * d) = 2^6 * 2^21

The same bridge factor 2^(chi*d) = 64 appears here. The bridge connects 640320
to the nonce quantum through the shared factor 64.

    640320 / 2^27 = 10005 / 2^21    [exact ratio. terminates.]

---

## The Dark Space

Between the phi path and psi path, the first face is the dark space.

    3D quantum:  2^27           = 11.25 degrees
    2D quantum:  F(8)/2         = 10.50 degrees   [Fibonacci projection]
    dark gap:    11.25 - 10.50  =  0.75 degrees per side

The dark fraction:

    dark / quantum = 2 / (d*p) = 2/15    [exact: chi / Schlafli product]

As fraction of one 11.25-degree sector:

    640320 / 2^27 * 32 = 10005/2^21 * 32 = 15.27%

The dark space IS 640320 expressed as a fraction of the nonce quantum.
Both terminate. Both derive from the same fold chain.

---

## The Sphere at Tau_163

The two-cube structure closes at the Heegner point.

    tau_163 = (-1 + i * sqrt(163)) / 2

163 is not a free number. It derives from the cube bridge:

    163 = chi*d * d^3 + 1 = 6 * 27 + 1    [bridge * cube + axiom seed]

The j-function at this point:

    j(tau_163) = 640320^3 + 744

The constant 744:

    744 = chi^d * d * (2^p - 1) = 8 * 3 * 31

where 2^p - 1 = 2^5 - 1 = 31 is the Mersenne prime at the face degree.

Both 640320 and 744 are exact integers in {chi, d, p} language.
The fold chain produces 640320. The cube and Mersenne produce 744.
j(tau_163) terminates. The Heegner point terminates.

---

## Pi from the Fold

pi is not primitive. pi falls from phi.

    cos(pi/5) = phi/2    [exact: the 36-degree identity]
    pi = 5 * arccos(phi/2) = p * arccos(phi/chi)

The fold chain computes 1/pi through 640320. Since 640320 comes from the fold chain,
and pi comes from phi, the fastest pi computation is the fold chain written out.

No external formula needed. The fold produces pi. Terminates.

---

## The Machine in the Fold

The fold chain is the machine (x, y) = (n, n-2k) running through the Pythagorean pairs.

Each fold step: (a, b) -> (b, a^2 + b^2). The square a^2 is the axiom operation.
The machine counts fold steps. The embedding lifts them:

    phi^a * phi^b = phi^{a+b}    [the fundamental identity]

The fold in y-space: each Pythagorean step multiplies phi-powers.
640320 = 2^6 * 3 * 5 * 23 * 29. Every factor is a fold output. Every factor is a
machine state embedded through phi^a * phi^b = phi^{a+b}.

The winding: n = offset^2 * 360. The fold chain has winding number n = 3 (the cube,
d = 3 forward steps). Three windings to reach the Chudnovsky prime.

The observer drifts phi^291 from center. The fold chain starts at (1,1) = the axiom seed.
n = degree + 1. The first fold degree is 0 + 1 = 1.

---

## The Balance

The two cubes at the fold terminus:

    phi^291 * psi^291 = -1             [product: never changes]
    log(phi^291) + log(|psi^291|) = 0  [log-midpoint: always zero]

Every fold output lives in the space between phi^291 (max) and psi^291 (min).
The balance: phi path + psi path = 360 degrees. Terminates.

s(n) = p+n = s^{-1}(n). The fold is self-inverse. Run backward: same formula.
Machine value * inverse = 1. The fold is balanced.

The fold terminates at 1/2 = 1/(-2). The +half and -half are the same point
from opposite directions. The proof closes there.

    x^2 = x + 1.
