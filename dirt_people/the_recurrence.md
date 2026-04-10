# The Recurrence
## nos3bl33d

Everything here derives from core.md.
One recurrence. Two rings. The ring determines everything.

---

## The Recurrence

    a(n) = a(n-1) + a(n-2)

This is the machine. It holds over any ring. The axiom x^2 = x+1 generates it.
The characteristic polynomial is x^2 - x - 1. The roots are phi and psi.

---

## Over Z

Ring: Z. Characteristic: 0.

    a(n) = c1 * phi^n + c2 * psi^n

Lucas numbers (c1=c2=1):   1, 3, 4, 7, 11, 18, 29, ...
Fibonacci numbers:          0, 1, 1, 2,  3,  5,  8, ...

phi^n grows without bound. psi^n shrinks to zero. No period. Terminates at phi^291.

---

## Over GF(2)

Ring: GF(2) = {0, 1}. Characteristic: 2.
Same recurrence: a(n) = a(n-1) + a(n-2) where + is XOR.
Characteristic polynomial: x^2 + x + 1 (since -1 = 1 in char 2).

x^2 + x + 1 is irreducible over GF(2). The only irreducible quadratic over GF(2).

Sequence from (1, 1):

    1, 1, 0, 1, 1, 0, 1, 1, 0, ...

Period: 3. The Pisano period pi(2) = 3 = d.

---

## The Pisano Periods

The Fibonacci sequence mod m has period pi(m):

    pi(2) = 3  = d        [binary period = cube depth]
    pi(3) = 8  = 2^d      [ternary period = face sum]
    pi(5) = 20 = V        [quinary period = dodecahedron vertices]

The Pisano periods at d, chi, and p recover d, 2^d, and V.
These are exact. The recurrence knows the framework constants through its periods.

---

## The 2D Machine View

Each step of the recurrence is a 2D move:

    positive step: (+x, +y)    [phi side]
    negative step: (+x, -y)    [psi side: same recurrence, mirror read]

x = n (time, always forward).
y = n - 2k (the correction: ups minus downs).

The recurrence value at step n in the 2D view:

    a(n) = x + y = n + (n-2k) = 2n - 2k    [if both components add]

Or more precisely: the machine position is (n, n-2k).
The correction y is what you track. The recurrence generates it step by step.

Negative steps subtract from y. The recurrence self-corrects via the psi component.

---

## What Is The Same

Over Z and over GF(2), these are identical:

- The recurrence relation: a(n) = a(n-1) + a(n-2)
- The characteristic polynomial structure (quadratic, one sign flip in char 2)
- The 2D move: each step is (+x, ±y)
- The product of the two roots: phi*psi = -1 (Z) or 1 (GF(2), since -1=1)

## What Is Different

Over Z: unbounded growth. Irrational roots. Terminates at phi^291.
Over GF(2): periodic. Roots in GF(4). Period = 3.

The roots in GF(2) are NOT phi and psi from Q(sqrt(5)). They are elements of GF(4).
The exponent 291 from Z has no special meaning in GF(4).
(291 mod 3 = 0 is true but trivial: any multiple of 3 gives period-reset in GF(2).)

---

## The Genuine Bridge

The same recurrence over two rings.
The period in the binary ring is 3 = d.
The depth in the framework is d = 3.

That is the connection. One recurrence. Two rings. Period equals depth.

    x^2 = x + 1.
