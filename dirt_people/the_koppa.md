# The Koppa
## nos3bl33d

Everything here derives from core.md.
The crossing angle is koppa. Pi and koppa are the same mechanism.

---

## The Clearance

The axiom:

    x^2 = x + 1
    x^2 - x = 1
    (x - 1/2)^2 = 5/4

Completing the square adds 1/4 to both sides. That 1/4 is koppa.

    koppa = 1/4

The system oscillates around x = 1/2 but never reaches it. The 1/4 is the minimum
clearance. Koppa prevents actual division by zero. The fixed point x = 1/2 is the
place where x = 1/x — the inversion singularity. The axiom keeps koppa distance from it.

As an angle:

    koppa = pi/4 = 45 degrees    [exact]

The oscillation sweeps in 45-degree steps. 8 steps = 360 degrees. 8 is the SHA word count.

---

## The Two X Patterns

Draw the unit square. Two X patterns:

    corner X:  diagonals from corner to corner,      length sqrt(2)
    face X:    lines from edge midpoint to midpoint,  length 1

The difference (exact):

    diff_2D = sqrt(2) - 1 = tan(pi/8) = tan(koppa/2)

This is exact. The 2D corner offset IS the half-angle tangent of koppa. Terminates.

---

## The 3D Extension

Extend to a cube:

    corner: space diagonal,   length sqrt(3)    [corner-to-corner]
    face:   face center path, length sqrt(2)    [face-to-face]

The conjugate identity:

    (sqrt(3) - sqrt(2)) * (sqrt(3) + sqrt(2)) = 3 - 2 = 1    [exact]

The corner diff and the face sum are perfect inverses. Their product is 1. Terminates.

    corner diff * face sum = 1/pi * pi = 1

The 0.15% gap between sqrt(3)-sqrt(2) and 1/pi is the projection artifact. It appears
identically in Lambda, in the oracle correction, in the Weinberg offset. It is one thing,
not three. The framework carries it as a single error from the 2D→3D projection.

---

## The Three Directions

One oscillation. Splits by inversion:

    x oscillates forward       [the phi path: +half]
    y oscillates backward      [the psi path: -half = 1/(-2) from phi's side]
    z does not oscillate       [z is the corner: z = x^3 + y^3]

x and y are perfect inverses: same amplitude, opposite direction. They cancel at the
crossing. The crossing point is x = 1/2 = the singularity. Koppa prevents hitting it.
What remains after x and y cancel is z.

The z value (exact):

    z = phi^3 + psi^3 = chi^2 = 4

Always. For every step of the spiral. z is forced by chi=2. Terminates.
You cannot compute z from x alone. You cannot compute z from y alone. z requires both
cubes together. This is why chi is not free — z must close, and only chi=2 closes it.

---

## The Circle

8 SHA words × koppa (45 degrees) = 360 degrees.

The entire circle is 8 koppa steps. The SHA word count is NOT arbitrary — it is the
number of 45-degree sectors in one full rotation.

The 3D bridge: sphere surface = 4 * circle.

    sectors = 8 * 4 = 32
    quantum  = 360 / 32 = 11.25 degrees = 2^27 nonce units    [exact]

---

## The Nonce Position

The nonce sits at one of 4 corners determined by the two oscillations:

    nonce in { center +/- dc +/- 2^27 }

    center  = F315-weighted sum of sha256s words    [phi oscillation]
    dc      = cubed difference (sha256s - sha256d)  [psi correction]
    2^27    = the quantum                           [the 3D sector width]

The 4 corners are:

    (center + dc + 2^27):  corner X, phi side, positive quantum
    (center + dc - 2^27):  corner X, phi side, negative quantum
    (center - dc + 2^27):  face X,  psi side,  positive quantum
    (center - dc - 2^27):  face X,  psi side,  negative quantum

Oracle: pick the closest corner. 1.722x improvement. Terminates.

---

## The Catalan Structure

The sha256s and sha256d oscillations are exact inverses running in opposite Fibonacci
directions. They cancel by:

    F(n+k) * F(n-k) = F(n)^2 - (-1)^(n-k) * F(k)^2

Left side: the product of the two oscillations.
Right side: a perfect square minus the Fibonacci remainder.

The F(k)^2 term is what survives. That is z. That is the corner offset.
The Catalan identity makes the z-extraction exact. It terminates.

---

## The Machine at Koppa

In the machine (x, y) = (n, n-2k), koppa is the angle of maximum efficiency.

    45 degrees: x_offset = y_offset    [equal movement. equilateral.]
    machine value = (n-k)/k = 1       [at 45 degrees: equal ups and downs]

The 45-degree path IS the machine at balance: k = n/2, y = 0, value = 1.
That is the critical line. That is koppa. That is 1/2.

The embedding at koppa:

    phi^a * phi^b = phi^{a+b}    [the fundamental identity]
    phi^{n/2} * phi^{-n/2} = 1    [balanced: forward and inverse cancel]

The winding number n at koppa = 8 (SHA words). 8 koppa steps = 360.
n = degree + 1. The observer steps through 8 sectors of 45 degrees each.

---

## The Termination

Pi is not the crossing angle in radians. Pi IS the reciprocal of the crossing gap.

    (sqrt(3) - sqrt(2)) * (sqrt(3) + sqrt(2)) = 1
    diff * sum = 1
    diff = 1/sum

The sum is pi (0.15% off). The diff is 1/pi (0.15% off). The product is exactly 1.
Pi is the mechanism that keeps the framework balanced. Koppa is pi/4. Both close.

The clearance at the singularity:

    koppa = 1/4
    koppa = pi/4

Both are the same 1/4 seen in the algebraic and geometric frames. The axiom and the
circle agree. The form closes. Terminates.

    x^2 = x + 1.
