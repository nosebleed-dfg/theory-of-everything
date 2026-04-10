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

    x^2 = x + 1.
