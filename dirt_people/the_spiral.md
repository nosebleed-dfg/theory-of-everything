# The Spiral
## nos3bl33d

Everything here derives from core.md.
The path n -> n+1 -> n+2 -> n+3 is a spiral. Not a line.

---

## The Move

At each step n -> n+1 you make two motions simultaneously:

    +1 on the line      [the axiom: x^2 = x+1, the discrete step]
    through the sphere  [the rotation: the n+1 sphere has radius 1 in log-space]

The sphere is exact:

    phi^n * psi^n = (-1)^n    [the n+1 sphere. always unit product. always.]

You cannot take the line step without also going through the sphere.
The sphere is not separate from the move. It IS the move.

The angle through the sphere is not 90 degrees. It is 90 degrees minus an offset.
That offset is the dark space at that step. It accumulates.

---

## The Two Paths

At every step, two paths are available:

    phi path:  forward, expanding, the +half direction
    psi path:  inverse, contracting, the -half direction = 1/(-2) from phi

The phi path IS the n+1 step you take.
The psi path IS every step you did not take.

The faces between the paths at step n = L_n = phi^n + psi^n.
These are exact integers. They terminate.

    L_1 = 1, L_2 = 3, L_3 = 4, L_4 = 7, L_5 = 11, L_6 = 18, L_7 = 29

---

## The Faces Not Taken

At step n, the number of paths not taken = !n (the derangements).

    !1 = 0     (one step, no alternatives)
    !2 = 1     (two steps, one alternative)
    !3 = 2     (three steps, two alternatives)

Sum through the cube (d = 3 steps):

    !1 + !2 + !3 = 0 + 1 + 2 = 3 = d

The alternative faces sum to d. The vertex degree IS the count of all paths
not taken through one cube cycle. Terminates.

As n grows, !n/n! -> 1/e. The phi path takes (1 - 1/e) of the angular space.
The psi path (all not-taken faces) takes 1/e of the angular space.

    phi path:  (1 - 1/e) * 360 = 227.6 degrees
    psi path:  (1/e)     * 360 = 132.4 degrees
    sum:        360 degrees     [always. terminates.]

---

## The Face Sum

Sum all faces along the d-step spiral:

    L_1 + L_2 + L_3 = 1 + 3 + 4 = 8 = 2^d

This is exact. Terminates. It holds because:

    sum L_1..L_n = L_{n+2} - 3    [Lucas identity, exact]
    at n = d = 3: L_5 - 3 = 11 - 3 = 8 = 2^d

Divide by one cube (2^d = 8):

    8 / 8 = 1 -> 1 * 360 = 360 degrees

The spiral always closes. No remainder. Terminates.

---

## Gamma

At each step n, the spiral moves at (90 degrees - offset_n) through the sphere.
The offset at step n is approximately 1/n of a full step.

The accumulated offset across all steps:

    gamma = lim (1 + 1/2 + 1/3 + ... + 1/n - ln(n))

Gamma is the total angular deficit of the two-cube spiral from pure 90-degree rotation.
It is not a free constant imported from outside.
It is the name for the offset that accumulates as the spiral runs.

From d = 3:

    gamma = H_{2^d} - d * ln(2) + d Bernoulli corrections
          = H_8    - 3 * ln(2)  + B_3

Every term controlled by d. Derives from the cube depth. Terminates.

---

## The SHA Structure

The SHA-256 structure falls from the spiral.

The phi spiral has 8 words because:

    SHA words = F - d - 1 = 12 - 3 - 1 = 8

F = 12 faces of the dodecahedron. d = 3 vertex degree. 1 = the identity.
Subtract the cube depth and the seed from the face count: 8 words remain.

The rounds:

    SHA rounds = (SHA words)^2 = 8^2 = 64

The bit counts:

    nonce bits = 2^p = 32
    hash bits  = 2^(p+d) = 2^8 = 256

All exact. All terminate. All fall from d and p.

---

## The Rotations

Every SHA rotation constant is a linear combination of d and p.

    Sigma0: 2, 13, 22    [2 = p-d, 13 = d+2p, 22 = 9d-p]
    Sigma1: 6, 11, 25    [6 = 2d,  11 = 2d+p, 25 = 5p]
    sigma0: 7, 18        [7 = 4d-p, 18 = 6d]
    sigma1: 17, 19       [17 = 4d+p, 19 = 8d-p]

Grand total of all rotations:

    37 + 42 + 25 + 36 = 140 = V * L_4 = 20 * 7

Every rotation is in {d, p} language. Every sum terminates.

---

## The n+3 Step

From the two-cube identity phi^3 = chi*phi + 1:

    phi^(n+3) = chi * phi^(n+1) + phi^n
    L_{n+3}   = chi * L_{n+1}  + L_n    [exact for all n]

The n+3 rule applies chi = 2 at every third step. The cube propagates by multiplying
the previous step by 2 and adding the one before.

Key chain:

    L_1 = 1  ->  L_4 = 7  ->  L_7 = 29

Each jump = +3 = one cube depth. Three jumps reach the Chudnovsky prime.
The spiral covers the Chudnovsky constant in exactly d cube-depth steps. Terminates.

---

## The Termination

The spiral terminates at phi^291. Beyond phi^291 is outside the universe.
The minimum step to reach the maximum = 291.

    291 = p * chi * L_7 + 1 = 5 * 2 * 29 + 1 = 290 + 1

Five Chudnovsky turns (p * chi * L_7 = 290 steps) plus the axiom seed (+1).

The spiral terminates where the +half and -half identify:

    1/2 = 1/(-2)

That is the center of the axiom. That is where the spiral lands.
The form closes. Gamma is the accumulated offset to that point.

    x^2 = x + 1.
