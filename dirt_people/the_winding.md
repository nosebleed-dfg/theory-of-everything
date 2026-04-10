# The Winding
## nos3bl33d

Everything here derives from core.md.
n is the winding number. The number of iterations to return to a.

---

## What n Is

n is not a step count. n is a winding number.

n = how many times the path must circle before the point returns to a.
n = how many times the path crosses itself before closure.
n = how many iterations make one full 360 degrees.

    n iterations * (360/n degrees) = 360    always.

Every form is closed. Every closed form has a winding number n.
The winding number is what you count. Not time. Not distance. Crossings.

---

## The Slice

The inverse of n is the slice:

    slice = 1/n    [one section of the circle]
    slice angle = 360/n degrees

The circle divides into n equal slices. Each slice is the angular distance
between one crossing and the next. The inverse gives you the slice width.

n = 3:   slice = 120 degrees    [the cube: minimum closed form]
n = 5:   slice = 72 degrees     [the pentagon]
n = 8:   slice = 45 degrees     [koppa: 8 SHA words * 45 = 360]
n = 32:  slice = 11.25 degrees  [the nonce quantum: 360/32 = 2^27 units]
n = 291: slice = 360/291 degrees [the spiral: closes after 291 winds]

---

## The Measurable Space

The space between a-to-a and b-to-b is cubed.

Between one visit to a and the next visit to a: one full winding.
That space is 3-dimensional: the path passes through (x, y, z) during one wind.
The 3D form gives it volume. The cube makes it measurable.

    one winding = one face of the cubic space
    n windings  = n faces
    faces counted = n

This is the same face count from the spiral: the faces between the phi path and psi
path are the Lucas numbers L_n. Each winding passes through one face.
The face count IS the winding number.

---

## A to A

The path starts at a. It spirals outward (phi) or inward (psi).
It does NOT return to a after one revolution — because phi is irrational,
the spiral does not close after exactly 360 degrees.

The path must wind n times before the crossing density allows return to a.
That n is the winding number for that orbit.

At n=3 (the cube): the minimum winding that allows closure.
    phi^3 = chi*phi + 1: three winds, one cube depth, the first return.

At n=291: the full spiral closes. phi^291 is the universe boundary.
After 291 winds, the path has covered the full space between phi^291 and psi^291.
The point returns — not to a exactly, but to a modulo the Planck scale.
psi^291 is sub-Planck. The return is exact at observable scales.

---

## The Inverse Gives the Slices

You have n crossings. The inverse 1/n gives you the slice between each crossing.

    n crossings -> 1/n slice width
    sum of all slices = n * (1/n) = 1 = 360 degrees    [always closes]

The slices are the measurable faces of the cubed space between a-to-a.
Count the slices: you have counted the faces.
Count the faces: you have the winding number.
Invert the winding number: you have the slice.

One operation. Forward (n) and inverse (1/n) are the same thing read two ways.
The winding number and the slice are the two faces of the same measurement.

---

## The Embedding Bridge

The winding number n lives in x-space. The embedding lifts it to phi-space:

    phi^a * phi^b = phi^{a+b}    [THE fundamental identity]

n crossings in x-space = phi^n in y-space. The winding number IS the exponent.
The slice 1/n in x-space = phi^{-n} in y-space.

    n crossings -> phi^n in the embedding
    1/n slice -> phi^{-n} = |psi|^n in the embedding
    product: phi^n * phi^{-n} = phi^0 = 1    [always closes]

The winding count and its inverse are the same pair viewed through the embedding.

---

## Framework Winding Numbers

    n=1:    trivial. one circle. no structure.
    n=2:    chi. the minimum non-trivial. the mirror crossing count.
    n=3:    d. the cube. the first closed 3D form. minimum winding with volume.
    n=5:    p. the pentagon. 5 faces of the dodecahedron.
    n=7:    L_4. the fourth Lucas prime. first Mersenne-prime-indexed winding.
    n=8:    2^d. eight SHA words. eight koppa slices = 360.
    n=20:   V. twenty dodecahedron vertices. twenty slices of 18 degrees.
    n=32:   8*4. the nonce sectors. the 2D*3D bridge quantum.
    n=291:  the universe. 291 winds. phi^291 = the unit.

Each is a closed form. Each has a measurable cubed space between a-to-a.
Each has a slice: 360/n degrees. Each has a face count: n.

---

## n as Offset Square Times 360

    n = offset^2 * 360

The winding number is the square of the offset times the full circle (360 degrees).

offset = the (x, y) displacement from the axiom center.
offset^2 = the area of the offset. the square. the axiom operation.
n = that area times 360. the winding count in degrees.

From the axiom: x^2 = x + 1. The square of the offset IS the axiom.

    offset = x:     offset^2 = x^2 = x + 1 = phi + 1 = phi^2
    n = phi^2 * 360 = (phi+1) * 360    [for unit phi offset]

The winding number n is the axiom (x^2) scaled to the circle (360).
The axiom squares the offset. The circle converts it to a count.

The face at offset x:

    face = phi^291 / n = phi^291 / (x^2 * 360)

Minimum face (maximum offset, x = phi^291):

    face = phi^291 / (phi^582 * 360) = 1 / (phi^291 * 360) -> 0    [sub-Planck. terminates.]

---

## n as (x, y) Offset

n is not a scalar. n is an offset: a little x and a little y.

    n = (x_offset, y_offset)

Each combination of x and y offset gives a different number of path crossings
before the point returns to a. Different offsets → different winding numbers.

The crossing count to return to a depends on the ratio x_offset / y_offset:
- Rational ratio p/q: returns after a finite number of crossings.
- Irrational ratio (phi): never exactly returns. Approximates return at ^291.

---

## The Equilateral Path: 45 Degrees

The maximum efficiency path is at 45 degrees = koppa.

    45 degrees: x_offset = y_offset    [equal x and y movement. equilateral.]

This is the fastest distance per move. Equal x and y gives the maximum diagonal
coverage. Any other angle wastes motion on one axis.

    angles < 45:   more x than y. slower coverage.
    angles > 45:   more y than x. slower coverage.
    angle = 45:    equal. maximum. the equilateral case. koppa.

The 45-degree path is the path that covers phi^291 worth of space in the
minimum number of crossings. It is the most efficient encoding of the spiral.

---

## The Minimum Face

At 45 degrees (the equilateral, maximum-efficiency path):

    face = 2^291 / 360^2

2^291: the universe boundary expressed in binary steps.
360^2: the area of the full angular space (circle squared — 2D).

The face at the equilateral path is the minimum measurable unit of the 2D
angular space at universe scale. It is the finest grain: the most path crossings
possible per unit of angular area, achieved at exactly 45 degrees.

    2^291 / 360^2 = 2^291 / (2^6 * 3^4 * 5^2)
                  = 2^285 / (3^4 * 5^2)
                  = 2^285 / 2025

This is the crossing density at maximum efficiency. The number of distinct
path crossings per unit of angular area when the path travels at koppa = 45 degrees.

---

## Why 45 Degrees Is Maximum

The path goes away from the axiom at an angle. The axiom is at x = 1/2 (the center).
The path that departs at 45 degrees is equidistant from both axes. This is koppa:

    koppa = pi/4 = 45 degrees    [from completing the square: the 1/4 term]

At 45 degrees: the path has maximum information per step. Each crossing encodes
one full unit of x AND one full unit of y simultaneously. No redundancy.

At any other angle: the path has less than maximum information per step.
It is biased toward one axis. The crossing density is lower.

The equilateral (45-degree) path is the path of maximum crossing density.
Maximum crossing density = minimum face = 2^291 / 360^2.

---

## 69 Degrees: The Minimum Dark Space

The circle splits at the universe boundary:

    291 degrees:  the phi side    [p * chi * L_7 + 1 = 5*2*29 + 1]
    69 degrees:   the dark side   [360 - 291 = 69]
    sum:          360 degrees     [always]

69 is not arbitrary:

    69 = d * (p^2 - chi) = 3 * 23    [vertex degree * backward fold]

23 = p^2 - chi is the backward fold — the minimum correction in the fold chain.
Scaled by d = 3 (the cube): 3 * 23 = 69 degrees.

The dark side is d times the backward fold projected to the circle.
It cannot be smaller. 23 is the minimum backward fold. d is the minimum cube depth.
69 is the minimum dark angle. The smallest dark point of view possible.

Inverse of 69 in the circle:

    360 / 69 = 360 / (3*23) = 120/23    [the slice count at 69-degree dark resolution]

The dark space at minimum = 69 degrees = 1 backward fold, cubically scaled.
Its inverse in the full circle: 120/23 slices of 69 degrees each cannot tile 360
evenly (120/23 is not an integer). The dark space is incommensurable with the circle.

That incommensurability IS the dark. It is what cannot be counted by the phi ruler.
The 69-degree remainder is the gap the phi path cannot fill.

---

## 8/15 Makes 2^291

The dark fractional remainder:

    8/15 = 2^d / (d*p)

    numerator:   8 = 2^d    [the cube face count. the SHA word count.]
    denominator: 15 = d*p   [the Schlafli product.]

The universe boundary: 291 = d * 97. Exactly 97 cube cycles.

    291 / d = 97    [exact integer. 97 cube depths to reach the universe.]

The numerator 8 = 2^d, raised to 97 cube cycles:

    8^97 = (2^d)^(291/d) = 2^291    [exact.]

The 8 in 8/15 IS the generator of 2^291.
The dark fractional remainder carries 2^d in its numerator.
That generator, cycled through 291/d = 97 cube depths, produces the universe scale.

The dark is not separate from the universe. The dark remainder (8/15) encodes the
generator (2^d = 8) that builds the universe (2^291) through repeated cubing.

    8/15 -> 8^97 = 2^291    [the dark generates the universe. one mechanism.]

    x^2 = x + 1.
