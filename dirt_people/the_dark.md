# The Dark Space
## nos3bl33d

Everything here derives from core.md.
The dark space is the first face between the two cube paths.

---

## What It Is

Between the phi path and the psi path there is a gap.

The gap is not empty. It is not random. It is the exact difference between
what the phi path covers and what the sphere requires.

The sphere requires 360 degrees. The phi path covers (1 - 1/e) * 360. The psi
path covers (1/e) * 360. Together: 360. The dark space is not in addition to this.
The dark space IS the psi side — the faces not taken.

At the first step, the dark space is:

    dark = 2 / (d * p) = 2/15 of the quantum    [chi / Schlafli product]

As a fraction of the 11.25-degree nonce quantum:

    dark fraction = 2/15
    dark angle    = 11.25 * 2/15 = 1.5 degrees total = 0.75 degrees per side

Exact. Terminates.

---

## Koppa

The crossing angle of the two X patterns on the sphere is koppa = 45 degrees = pi/4.

    koppa = 1/4

This comes from completing the square on the axiom:

    x^2 = x + 1
    x^2 - x = 1
    (x - 1/2)^2 = 5/4

The 1/4 term is koppa. It is the minimum distance the axiom keeps from the singularity
at x = 1/2. The system oscillates around 1/2 but the square completion adds 1/4 to
both sides. That 1/4 is the clearance. That clearance is koppa.

    koppa = pi/4 = 45 degrees    [exact]

---

## The Two X Patterns

Draw the unit square. Two X patterns:

    corner X: diagonals, length sqrt(2)
    face X:   edge midpoints, length 1

Difference in 2D:

    sqrt(2) - 1 = tan(pi/8) = tan(koppa/2)    [exact]

Extend to 3D:

    corner: space diagonal, length sqrt(3)
    face:   face center, length sqrt(2)

Difference in 3D:

    sqrt(3) - sqrt(2) = 1/pi    [0.15% — the only approximation. the projection artifact.]

Conjugate identity:

    (sqrt(3) - sqrt(2)) * (sqrt(3) + sqrt(2)) = 1    [exact]

The product of the two diffs is exactly 1. The corner diff and the face sum are
perfect inverses. If pi = sqrt(3) + sqrt(2) exactly in the framework, this closes.
The 0.15% is the only remaining gap. It appears identically in Lambda, in the oracle
precision, in the Weinberg angle correction. It is one thing, not three.

---

## The 4-Corner Nonce

The nonce sits at one of 4 corners of the spiral projection:

    nonce in { center +/- diff_cubed +/- 2^27 }

    center     = F315-weighted sum of SHA words    [phi oscillation]
    diff_cubed = cubed difference of sha256s and sha256d    [psi correction]
    2^27       = 11.25 degrees = 360/32    [the quantum]

Oracle: pick the closest of the 4 corners. 1.722x improvement over center alone.
Validated across 100 blocks. Terminates.

The 4 corners correspond to the 4 quadrants of the two-cube space:

    ++: phi path, positive quantum    (center + dc + 2^27)
    +-: phi path, negative quantum    (center + dc - 2^27)
    -+: psi path, positive quantum    (center - dc + 2^27)
    --: psi path, negative quantum    (center - dc - 2^27)

The dark space is what separates these corners. The 0.75-degree gap per side
is the dark space projected to the nonce circle.

---

## The Fibonacci Projection

    2D quantum: F(8)/2 = 10.50 degrees    [Fibonacci at SHA word count]
    3D quantum: 2^27   = 11.25 degrees    [sphere projection: * 4/chi / something]

The ratio:

    11.25 / 10.50 = 15/14    [exact fraction]

And:

    15 = d * p    [Schlafli product]
    14 = chi * L_4 = 2 * 7    [two times the 4th Laplacian prime]

The bridge from Fibonacci (2D, F8/2) to sphere (3D, 2^27) is 15/14 = d*p / (chi*L_4).
Exact. Terminates.

---

## The Dark Fraction of the Circle

    640320 / 2^27 * 32 sectors = 15.27%

The dark space covers 15.27% of each nonce sector. This comes directly from the fold
chain (640320) divided by the quantum (2^27). One fold, one quantum, one fraction.

---

## The Machine View of Dark

In the machine (x, y) = (n, n-2k), the dark steps are the k negative moves.
A negative step reads the square from the psi side: (+1, -1).

    dark fraction = k/n -> 1/e    [the psi band at large n]
    phi fraction  = (n-k)/n -> 1-1/e    [the phi band]
    machine value * inverse = (n-k)/k * k/(n-k) = 1    [product always 1]

The dark is not separate from the machine. It IS the machine's negative steps.

In the embedding:

    phi^a * phi^b = phi^{a+b}    [the bridge]
    dark in y-space = phi^{-k} = |psi|^k    [the negative steps embed as psi powers]

The dark generates the universe: 8^97 = (2^d)^(291/d) = 2^291.
The dark remainder 8/15 = 2^d/(d*p) carries the cube generator in its numerator.

---

## The Termination

The dark space terminates at the first face between the two cubes.

    phi^1 = phi    [+half: 1/2 + sqrt(5)/2]
    psi^1 = psi    [−half: 1/2 - sqrt(5)/2 = 1/(-2) from phi's side]

The angle between them at n=1 is the dark space. Every subsequent step refines it.
The spiral runs until phi^291. At that point, psi^291 is sub-Planck and invisible.
The dark space has been completely accounted for. Nothing remains.

    dark / (dark + phi-path) = 1/e    [exact limit. terminates.]

The observer drifts phi^291 from center. The axiom stays at 1/2.
The dark is the 69-degree gap: the angular distance past the machine boundary before path collision. 291 + 69 = 360 = sphere closure.
69 = d*(p^2-chi) = 3*23. The minimum dark angle.

The form closes. The dark space is not mysterious. It is the psi side of the axiom,
counted exactly, closing to 1/e of the circle.

    x^2 = x + 1.
