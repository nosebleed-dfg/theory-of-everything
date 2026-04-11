# Pi from the Pythagorean Theorem
## nos3bl33d | (DFG) DeadFoxGroup | x^2 = x + 1 | A Wake In Outerspace

Everything here derives from core.md.
Pi is not assumed. Pi is not an input. Pi falls out of a right triangle.

---

## The Axiom

    a^2 + b^2 = c^2

That's it. A right triangle. Two legs. One hypotenuse. The squares add.

This is not a formula you learn. This is the rule the universe runs on.
Everything below is one move applied over and over. Same move each time.

---

## The Unit Circle Is the Axiom

Set c = 1. The hypotenuse is 1. Then:

    a^2 + b^2 = 1

Every point (a, b) on this curve satisfies the axiom. The curve IS a circle.
The circle is not something separate from the axiom. The circle IS the axiom
with c = 1.

Now name the legs:

    cos(theta) = a    [the horizontal leg]
    sin(theta) = b    [the vertical leg]

These are just names for the two sides of the triangle. The angle theta
is how far around the circle you've gone. That's it.

So:

    cos^2(theta) + sin^2(theta) = 1

This is not a "trig identity." This IS the Pythagorean theorem.
Same equation. Different names for the parts.

---

## What Is Pi?

Pi is the distance from (1, 0) to (-1, 0) going the long way around
the top of the circle. The half-circumference.

We haven't assumed this distance exists, or what number it is.
We're going to COMPUTE it. From the axiom.

Take two points on the circle, very close together. The tiny straight
line between them has length:

    ds^2 = dx^2 + dy^2

That's the axiom again. Pythagorean theorem on a tiny triangle.
The tiny triangle has legs dx and dy and hypotenuse ds.

Now, these points are on the circle x^2 + y^2 = 1. Differentiate:

    2x dx + 2y dy = 0
    dy = -(x/y) dx

Plug into ds^2:

    ds^2 = dx^2 + (x/y)^2 dx^2
         = dx^2 (1 + x^2/y^2)
         = dx^2 (y^2 + x^2) / y^2
         = dx^2 / y^2          [because x^2 + y^2 = 1. the axiom.]

So:

    ds = dx / sqrt(1 - x^2)

Add up all the tiny pieces from x = -1 to x = 1 (the top half):

    pi = integral from -1 to 1 of dx / sqrt(1 - x^2)

That IS pi. Defined. Computed. From the axiom. No assumptions.

The integrand 1/sqrt(1 - x^2) is literally 1/(the other leg).
You're adding up "one over the vertical leg" as the horizontal leg
sweeps across. That's the circumference. That's pi.

---

## The Tangent: Another Pythagorean Ratio

In the right triangle with legs 1 (adjacent) and u (opposite):

    hypotenuse = sqrt(1 + u^2)    [the axiom]
    tan(theta) = u/1 = u          [ratio of the legs]

The angle whose tangent is u: theta = arctan(u).

And here's the key step. If you differentiate arctan:

    d(theta) = du / (1 + u^2)

Why? Because sec^2(theta) = 1 + tan^2(theta), and THAT comes from
dividing cos^2 + sin^2 = 1 by cos^2. The axiom again.

So:

    arctan(t) = integral from 0 to t of du/(1 + u^2)

Every ingredient: Pythagorean.

---

## Pi = 4 * arctan(1)

Take the isosceles right triangle: legs 1 and 1.

    1^2 + 1^2 = 2        [the axiom]
    hypotenuse = sqrt(2)
    tan(theta) = 1/1 = 1
    theta = arctan(1)

This angle is pi/4. Why? Because four copies of a right angle
fill a full turn (this is what "right angle" means geometrically --
the angle that tiles four times around a point). A full turn = 2*pi.
So right angle = pi/2. And the isosceles right triangle splits
that right angle in half: theta = pi/4.

Therefore:

    pi = 4 * arctan(1) = 4 * integral from 0 to 1 of du/(1 + u^2)

And if you expand 1/(1 + u^2) as a geometric series and integrate
term by term:

    pi/4 = 1 - 1/3 + 1/5 - 1/7 + 1/9 - ...

The odd numbers. The reciprocals of the odd numbers, alternating sign.
That series IS pi/4. Derived. Not assumed.

---

## Pi = arctan(1) + arctan(2) + arctan(3)

Now the framework enters.

The axiom x^2 = x + 1 forces chi = 2, d = 3, p = 5.
(See core.md. The derivation is: phi^3 + psi^3 = chi^2, giving chi = 2.
Then d = chi^2 - 1 = 3, p = chi^2 + 1 = 5.)

Watch what happens when you add arctan(1) + arctan(2):

    tan(arctan(1) + arctan(2))
    = (1 + 2) / (1 - 1*2)
    = 3 / (-1)
    = -3

The tangent addition formula is Pythagorean -- it's just combining
two right triangles and reading the new ratio. And the result is -3.

Since arctan(1) + arctan(2) is between pi/2 and pi (both angles
are positive, their sum passes the quarter-turn), and the tangent is -3:

    arctan(1) + arctan(2) = pi - arctan(3)

Rearrange:

    arctan(1) + arctan(2) + arctan(3) = pi

Three right triangles. Legs (1,1), (1,2), (1,3). Their base angles
sum to exactly half a turn.

And the arguments are (1, chi, d) = (1, 2, 3). The unit, the Euler
characteristic, the vertex degree. All forced by the axiom.

---

## The Golden Ratio

Now the deepest path. Apply the axiom to legs 1 and 2:

    1^2 + 2^2 = 5
    diagonal = sqrt(5)

From sqrt(5):

    phi = (1 + sqrt(5)) / 2 = 1.618...
    psi = (1 - sqrt(5)) / 2 = -0.618...

These satisfy:

    phi^2 = phi + 1    [the golden axiom]

That's right: x^2 = x + 1 falls directly out of the right triangle
(1, 2, sqrt(5)). The axiom is the FIRST nontrivial output of the
Pythagorean theorem on integer sides.

---

## cos(pi/5) = phi/2

This is the bridge. The thing that connects the golden ratio to pi.

Start with: 5 * theta = pi, so cos(5*theta) = cos(pi) = -1.

Expand cos(5*theta) using the angle-addition formula repeatedly
(which is Pythagorean -- combining triangles):

    cos(5*theta) = 16c^5 - 20c^3 + 5c

where c = cos(theta). Setting this to -1:

    16c^5 - 20c^3 + 5c + 1 = 0

Factor out (c + 1) [because c = -1 is a solution for theta = pi]:

    (c + 1)(16c^4 - 16c^3 - 4c^2 + 4c + 1) = 0

The quartic is a perfect square:

    (4c^2 - 2c - 1)^2 = 0

Solve 4c^2 - 2c - 1 = 0:

    c = (1 +/- sqrt(5)) / 4 = phi/2 or psi/2

Since theta = pi/5 is in the first quadrant, cos(theta) > 0.
Only phi/2 is positive. So:

    cos(pi/5) = phi/2

And here's the kill shot. Verify by plugging c = phi/2 back into
4c^2 - 2c - 1:

    4*(phi/2)^2 - 2*(phi/2) - 1
    = phi^2 - phi - 1
    = (phi + 1) - phi - 1    [used phi^2 = phi + 1]
    = 0

The proof REDUCES TO THE AXIOM. The entire thing collapses to
phi^2 = phi + 1. That's all it is.

Therefore:

    pi/5 = arccos(phi/2)
    pi = 5 * arccos(phi/2) = p * arccos(phi/chi)

Pi equals p times the arccosine of phi over chi. Every piece forced
by the axiom. Exact. Not approximate. Not a series. Exact.

---

## The Dodecahedron Dihedral

The dodecahedron {5, 3} is the polyhedron forced by the axiom: p = 5
faces per polygon, d = 3 faces meeting at each vertex, chi = 2 as the
Euler characteristic.

Its dihedral angle (the angle between two adjacent faces) is:

    delta = pi - arctan(2) = pi - arctan(chi)

So: delta + arctan(chi) = pi. The dihedral angle and arctan(2) are
supplementary. They add to exactly one half-turn.

Verify: cos(delta) = cos(pi - arctan(2)) = -cos(arctan(2)).
From the (1, 2, sqrt(5)) triangle: cos(arctan(2)) = 1/sqrt(5).
So cos(delta) = -1/sqrt(5) = -1/sqrt(p).

The dihedral angle of the dodecahedron is determined by 1/sqrt(p),
and p comes from the axiom. The dodecahedron's face-meeting angle
is pi minus the arctan of the Euler characteristic.

---

## The 336-Degree System

360 is Babylonian. It's a good number for dividing easily.
But it's not the natural angular quantum of the axiom.

    360 = chi^3 * d^2 * p = 8 * 9 * 5

vs.

    336 = chi^4 * d * L_4 = 16 * 3 * 7

336 is the number of oriented edges of the Klein quartic -- the
genus-3 surface with maximal symmetry. Genus 3 because d = 3.
L_4 = 7 because phi^4 + psi^4 = 7.

The Klein quartic hits the Hurwitz bound:
    max automorphisms for genus g = 84(g - 1)
    at g = d = 3: 84 * 2 = 168 = |PSL(2,7)|
    full group (with orientation-reversing): 336

The half-turn in 336 degrees is 168 = |PSL(2,7)|. That's the
orientation-preserving automorphism group of the maximally symmetric
surface at the axiom's dimension.

The difference:

    360 - 336 = 24 = chi^3 * d = |S_4|

24 is the order of the cube's rotation group. The gap between
the Babylonian circle and the axiom's circle is exactly the cube.

The dodecahedral dihedral angle in 336 degrees:

    delta_336 = 168 - (336/(2*pi)) * arctan(2)

The leading term is 168, the PSL(2,7) order. The dihedral angle
starts at the half-turn and subtracts the arctan(chi) correction.

The byte conversion:

    256/336 = chi^4 / (d * L_4) = 16/21

This converts from geometric angle to computational address.
Exact. Rational. Every factor from the axiom.

    336 * 1024/21 = 2^14

336 grid steps times the byte-to-face ratio gives exactly 2^14.
The byte wrap. The grid period. Terminates.

---

## The Full Picture

Four paths from a^2 + b^2 = c^2 to pi:

    1. CIRCUMFERENCE: unit circle -> arclength integral -> pi

    2. ISOSCELES: (1,1,sqrt(2)) -> arctan(1) -> pi = 4*arctan(1)

    3. GOLDEN: (1,2,sqrt(5)) -> phi -> cos(pi/5) = phi/2 -> pi = 5*arccos(phi/2)

    4. DECOMPOSITION: arctan(1) + arctan(chi) + arctan(d) = pi

All four give the same number. All four use nothing but the axiom.

Pi is not a free constant. It is not a mysterious transcendental
that appears from nowhere. It is the angular content of the geometry
that a^2 + b^2 = c^2 forces into existence.

The Pythagorean theorem gives you the circle.
The circle gives you the circumference.
The circumference gives you pi.
One move. Applied repeatedly. Same move every time.

    a^2 + b^2 = c^2
        -> sqrt(5)
            -> phi, chi=2, d=3, p=5
                -> pi = p * arccos(phi/chi)

---

## The Chain

    step    what happens                              what you get
    ------- ----------------------------------------  ----------------
    0       a^2 + b^2 = c^2                          the axiom
    1       set c = 1                                 unit circle
    2       arclength of the circle                   pi (defined)
    3       (1,1) triangle                            arctan(1) = pi/4
    4       pi = 4 * arctan(1)                        pi (first computation)
    5       (1,2) triangle                            sqrt(5)
    6       phi = (1+sqrt(5))/2                       golden ratio
    7       phi^2 = phi + 1                           golden axiom
    8       cos(pi/5) = phi/2                         the bridge
    9       pi = 5*arccos(phi/2)                      pi (second computation)
    10      arctan(1)+arctan(2)+arctan(3)=pi          pi (third computation)
    11      dodecahedral dihedral = pi - arctan(2)    dihedral from axiom
    12      336 = chi^4 * d * L_4                     natural angular quantum
    ------- ----------------------------------------  ----------------

Every row is one application of a^2 + b^2 = c^2 or its direct consequences.
No external constants. No assumptions. No circular definitions.

Pi is output. The axiom is input. Same as everything else in the framework.
