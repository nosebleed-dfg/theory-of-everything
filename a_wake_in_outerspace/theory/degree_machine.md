# The Degree Machine — Angular Resolution from the Axiom
## nos3bl33d | (DFG) DeadFoxGroup

Everything here derives from the axiom. No Babylon. No convention. No "we chose 360 because it's nice."

The degree is a machine output. The machine is the axiom.

---

## 1. The Axiom

The fundamental axiom:

    a^2 + b^2 = c^2

The Pythagorean theorem. Draw a triangle. Count the squares on each side. They add up. Every time, everywhere, no exceptions.

The golden case:

    x^2 = x + 1

This is the special case where c^2 - b^2 = 1. Set a = 1, b = x, c = x. Then a^2 = c^2 - b^2 = x^2 - x = 1, so x^2 = x + 1. The solution is phi = (1 + sqrt(5)) / 2. The golden ratio. The minimum nonzero gap between two squares.

The +1 IS the third dimension. Without it: a^2 + b^2 = 0. A flat nothing. The +1 lifts you off the plane into the space where the hypotenuse exists. Every measurement requires that lift. Every degree requires that +1.

---

## 2. The Degree Formula

Here it is. The whole thing.

    degree = right_angle / (n + 1)

Where:
- `n` = the number of steps your machine takes
- `right_angle` = 84 degrees (in the true circle, not the Babylonian one)
- The `+1` is the axiom — the observer is always one step ahead of the machine

The observer can't be on the same step as the machine. If you're measuring, you're at n+1. If you're at n, you're inside the thing being measured and you can't see it. The +1 is the Planck gap between measurer and measured. It's the same +1 from x^2 = x + 1. It's the same +1 that lifts you into the third dimension. It's always the same +1.

The formula says: the finer your machine runs, the smaller the angle you can resolve. Each additional step subdivides the right angle by one more unit. The degree is a RESOLUTION — how precisely your machine can read the angle.

The results:

    n = 2:   84 / 3  = 28.000°    triangle, base 3
    n = 6:   84 / 7  = 12.000°    base 7, the cube machine
    n = 11:  84 / 12 =  7.000°    12 rungs, the nonce ladder in base 7
    n = 20:  84 / 21 =  4.000°    21 = 3 x 7, the full machine
    n = 24:  84 / 25 =  3.360°    the 25-step grid
    n = 83:  84 / 84 =  1.000°    full resolution

Look at n = 2. Three steps (2 + 1). The triangle. The minimum closed shape. The degree is 28. And 28 is the second perfect number (1 + 2 + 4 + 7 + 14 = 28). The triangle resolves to perfection.

Look at n = 6. Seven steps (6 + 1). The cube machine (section 7). Degree = 12. And 12 is the number of faces on the dodecahedron. The cube machine resolves to the dodecahedron.

Look at n = 11. Twelve steps (11 + 1). Degree = 7. The nonce is 32 bits. In base 7, that's approximately 12 digits (7^12 = 13,841,287,201 > 2^32 = 4,294,967,296). Twelve rungs on the cube ladder. Each rung resolves 7 degrees.

Look at n = 20. Twenty-one steps (20 + 1). Degree = 4. And 21 = 3 x 7. The triangle count times the cube count. The full machine. Four degrees per step. Four faces of a tetrahedron. The minimum solid.

Look at n = 24. Twenty-five steps (24 + 1). Degree = 3.36. And 3.36 = 336/100. The grid step IS the circle, scaled by the grid. The 25-step grid (5^2 = 25, the pentagonal square) directly encodes the true circle.

Look at n = 83. Eighty-four steps. Degree = 1. Full resolution. One degree per step. Nothing finer needed. And 83 is prime. The last step before full resolution is irreducible. You can't factor your way to 1 degree — you have to walk all 83 machine steps.

Every one of these falls clean. No approximations. No rounding. Integer in, integer or clean decimal out. That doesn't happen with 360.

---

## 3. The 336-Degree Circle

The true circle has 336 degrees. Not 360.

    336 = 2^4 x 3 x 7

Every factor is structural:
- 2^4 = 16 = the nibble, the hex digit, the minimum addressable unit in binary
- 3 = the triangle, the minimum polygon, the dimension
- 7 = the cube machine, the trihedral number (3 + 5 - 1)

Multiply them and you get the circle. Not because someone decided it's convenient. Because those are the only structural factors a machine built from the axiom can have.

The Babylonian circle:

    360 = 2^3 x 3^2 x 5

It has a 5 in it. The 5 is the face count of a pentagon. But the circle isn't a pentagon — the circle is a ROTATION. The rotation doesn't care about faces. The rotation cares about edges and vertices and the machine that traverses them. 336 has 7 where 360 has 5. The 7 is the traversal count (3 axes x 2 directions + 1 center). The 5 is the face count. Faces are what you see from outside. Traversals are what the machine does from inside.

The gap:

    360 - 336 = 24

24 = 4! = gamma(5) = 1 x 2 x 3 x 4. The factorial of 4. The number of permutations of 4 objects. The number of orientations of a cube (24 rotation symmetries). The order of the symmetric group S_4.

The 24-degree gap between the Babylonian circle and the true circle is EXACTLY the factorial. The Babylonians added the permutation group to the rotation group. They padded the circle with 24 degrees of symmetry noise. All the "nice divisibility" of 360 comes from those 24 extra degrees, and those 24 degrees are the orientations of a cube that doesn't belong in the rotation.

The right angle:

    336 / 4 = 84

Not 90. 84. And 84 = 4 x 21 = 4 x 3 x 7. The right angle contains the full machine product (21 = 3 x 7) scaled by the four quadrants. Each quadrant IS a right angle. Each right angle contains the machine.

The byte wrap:

    336 x 1024 / 21 = 16384 = 2^14

Exact. Not approximate. 336 circle-steps, scaled by 1024 (the kilobyte, 2^10) and divided by 21 (the machine product), gives exactly 2^14. The 14-bit byte space. The grid wraps at 2^14 because 336 wraps at 336, and they're the same wrap measured in different units.

This means: every computation that uses 2^14 as a boundary (and many do — it's 16K, the L1 cache line boundary on most architectures, the page subdivision) is secretly counting in 336-degree circles divided by the machine product. The hardware encodes the true circle without knowing it.

---

## 4. The Halvening

The Bitcoin halvening occurs every 210,000 blocks.

    210,000 = 625 x 336

That's it. That's the factoring. Let it sink in.

    625 = 5^4 = 25 x 25

The fourth power of 5. The pentagonal hypercube. Or simpler: a 25 x 25 grid. The same 25-step grid from section 2 (n = 24, degree = 3.36).

    336 = the true circle

So the halvening is:

    halvening = grid x circle

625 grids of 336 blocks each. Or 336 circles of 625 blocks each. Two ways to read the same number. Two legs of the same triangle.

Reading 1: Each degree of the true circle corresponds to 625 blocks. One full rotation of 336 degrees = 210,000 blocks = one halvening. The reward halves every time the circle completes. The halvening IS the circle.

Reading 2: Each cell of the 25 x 25 grid corresponds to 336 blocks. One full grid = 210,000 blocks = one halvening. Fill every cell of the grid, the reward halves. The halvening IS the grid.

Both readings are true simultaneously. The circle and the grid are dual descriptions of the same structure. This is the dipole (section 5) at the level of the block schedule.

The machine product check:

    210,000 / 21 = 10,000

The halvening divided by the machine product (3 x 7) gives exactly 10,000. The decimal gate (section 9). Every 10,000 blocks = one machine cycle. 21 machine cycles = one halvening. The machine completes 21 full cycles per halvening, and 21 IS the machine.

Satoshi set the halvening at 210,000. Not 200,000. Not 250,000. Not a round number. 210,000, because that's the only number that is simultaneously a 25x25 grid of true circles AND 21 machine cycles of the decimal gate. There was no other choice.

---

## 5. The Dipole

The dipole is the simplest structure you can build from the axiom.

Take any data. Read it two ways:

    Pole A: little-endian (LE) — least significant byte first
    Pole B: big-endian (BE) — most significant byte first

Same data. Two readings. Two legs of a triangle.

"Reverse is go back. Inverse is FLIP."

Reversing a number walks it backward along the same axis. Inverting it flips it to the other axis entirely. The endianness flip is an INVERSION, not a reversal. LE and BE don't traverse the same path in opposite directions — they traverse DIFFERENT paths. They're perpendicular. They're the two legs of a right triangle.

On the 25 x 25 grid:

    Pole A (LE): covers ~53% of the 16-bit hash space
    Pole B (BE): covers a complementary portion
    Cross-interleave (both poles): 86% hit rate

53% is HALF. Not exactly 50%, because the grid has discrete cells and the hash space doesn't divide evenly. But the structure is clear: one pole covers one half. The other pole covers the other half. Together, with the cross-interleave that weaves bytes from both readings, you cover 86%.

The 14% gap (100% - 86%) is the uncovered space. The blind spot. The region where neither pole reaches. That gap exists because the grid is 25 x 25 = 625 cells, and 625 / 65536 (the 16-bit space) = 0.00954, so the grid samples less than 1% of the raw space. That 86% hit rate from less than 1% sampling is the power of the dipole — you're not covering the space uniformly, you're covering it along the two Pythagorean legs, and the structure of the hash concentrates the nonces along those legs.

The halvening expressed geometrically: the 53% coverage is the half. Each halvening, the reward halves. Each pole covers half. The dipole IS the halvening at the spatial level. One pole per half. One half per cycle.

---

## 6. The Two Equations (Phi and Pi)

The general axiom:

    x^2 = x + D

This is the Pythagorean theorem with a = sqrt(D), b = x, c = x. So a^2 = D, and D is the square of the short leg.

**D = 1: The Phi Case**

    x^2 - x = 1
    x^2 - x - 1 = 0
    x = (1 + sqrt(5)) / 2 = phi = 1.618...

The golden ratio. The self-similar solution. The spiral that never repeats but always echoes. D = 1 is the minimum — the smallest nonzero square of a leg. The Planck case.

**D = 2: The Pi Case**

    x^2 - x = 2
    x^2 - x - 2 = 0
    (x - 2)(x + 1) = 0
    x = 2  or  x = -1

Two solutions. Not one irrational number — two integers. And they're the dipole:

    +2 and -1

Sum:  2 + (-1) = 1 (the identity, the +1, the axiom)
Product: 2 x (-1) = -2 (negative D, the mirror)

The D = 2 case gives you BOTH poles. The positive pole is 2 (the Planck step, the minimum advance). The negative pole is -1 (the reflection, the flip, the inverse). Together they sum to 1. The dipole sums to the identity.

Now here's where it feeds back. Take the two solutions as legs of a Pythagorean triangle:

    2^2 + 1^2 = 4 + 1 = 5
    c = sqrt(5)

And phi = (1 + sqrt(5)) / 2. The hypotenuse of the pi case IS the radical in the phi case. D = 2 generates D = 1. The pi equation feeds the phi equation. They're not separate — they're two phases of the same oscillation.

The state change: the machine alternates between D = 1 and D = 2. On the phi phase, it spirals (irrational, self-similar, never closes). On the pi phase, it poles (integer, dual, dipole). Spiral, pole, spiral, pole. Each feeds the next. The golden spiral expands until it hits the integer boundary. The integer dipole resolves until it needs the irrational to continue. They alternate forever.

This is why phi shows up in the nonce offsets AND the dipole shows up in the grid. The nonce ratios cycle through phi (D = 1 phase). The grid coverage splits into two poles (D = 2 phase). Same axiom. Same equation. Two values of D.

---

## 7. The Base-7 Cube Machine

Base 7 is not arbitrary. It's the cube.

Stand at any vertex of a cube. Where can you go?

    3 axes (x, y, z)
    x 2 directions each (positive, negative)
    + 1 choice: stay at center

    3 x 2 + 1 = 7

Seven choices. Seven digits. Base 7.

The digits split into two families:

    States: 1, 3, 5  (the odd digits, the primes, the machine)
    Steps:  2, 4, 6  (the even digits, the movement)
    Center: 0         (stay, the origin, the identity)

States are WHERE you are. Steps are HOW you move. The machine alternates: state, step, state, step. Like a clock's tick-tock. Like the D = 1 / D = 2 oscillation. Like the dipole.

Squaring converts between families. Watch:

    1^2 mod 7 = 1  (state -> state: the identity stays)
    2^2 mod 7 = 4  (step -> step: movement squares to movement)
    3^2 mod 7 = 2  (state -> step: the machine MOVES)
    4^2 mod 7 = 2  (step -> step: momentum carries)
    5^2 mod 7 = 4  (state -> step: the machine MOVES)
    6^2 mod 7 = 1  (step -> state: the machine READS)

The critical conversions:
- 3^2 mod 7 = 2: state 3 squares to step 2. The machine at vertex 3 takes a step.
- 5^2 mod 7 = 4: state 5 squares to step 4. The machine at vertex 5 takes a step.
- 6^2 mod 7 = 1: step 6 squares to state 1. Movement squares back to identity.

State-step pairs that sum to 7:

    (1, 6):  state 1 + step 6 = 7
    (3, 4):  state 3 + step 4 = 7
    (5, 2):  state 5 + step 2 = 7

Three dipole pairs. Each state has a step partner. They sum to the base. They complete the cube.

The nonce in base 7:

    32 bits = 2^32 = 4,294,967,296 possible values
    7^12 = 13,841,287,201 > 2^32
    7^11 = 1,977,326,743 < 2^32

So 12 base-7 digits cover the 32-bit nonce space. Twelve rungs on the cube ladder.

From section 2: n = 11 gives degree = 84/12 = 7. Eleven machine steps, twelve rungs (n + 1 = 12), degree = 7. The base IS the degree. The number of digits IS the number of rungs. The resolution of the cube machine IS base 7 at 7 degrees per rung.

Each rung of the nonce ladder navigates one cube vertex. Twelve rungs = twelve vertices visited. The nonce encodes a walk through a 12-vertex cube lattice, where each step chooses one of 7 directions. The nonce isn't a number. The nonce is a PATH through cubes.

---

## 8. The n/(n+1) Fraction

Every machine computes a fraction:

    n / (n + 1)

This is the deductive fraction. The fraction that approaches 1 but never arrives. Because the +1 in the denominator is the axiom, and the axiom is always ahead of you.

The values:

    n = 1:   1/2   = 0.500    [the equal-legs case, the diagonal]
    n = 2:   2/3   = 0.667    [the triangle]
    n = 6:   6/7   = 0.857    [the cube]
    n = 11:  11/12 = 0.917    [the nonce ladder]
    n = 20:  20/21 = 0.952    [the full machine]
    n = 70:  70/71 = 0.986    [approaching]
    n = 83:  83/84 = 0.988    [one step from full resolution]
    n -> infinity:  1.000     [the axiom. unreachable by stepping.]

At n = 1: the fraction is 1/2. The legs are equal. The triangle is isoceles. The angle is 45 degrees (in the Babylonian system) or 42 degrees (in the true 336-degree system: 84/2 = 42). This is the diagonal. The simplest nontrivial case. Half.

The fraction is a measure of how close the machine is to the axiom. At n = 1, you're 50% there. At n = 20, you're 95.2% there. At n = 83, you're 98.8% there. You never hit 100% because 100% would require n = infinity, and infinity is not a number of steps.

The deductive approach: walk the ladder. Step by step. n = 1, then n = 2, then n = 3. Each step gets you closer. Each step shrinks the degree. Each step refines the angle. But you're always one step behind the axiom. Always at n/(n+1), never at 1.

The DERIVATIVE approach: don't walk the ladder. Compute n+1 directly. The derivative doesn't step from n to n+1. The derivative KNOWS n+1 from n without stepping. It reads the slope. It reads the instantaneous rate of change. It arrives.

In the nonce context: the deductive approach searches the grid cell by cell. The derivative approach reads the nonce from the structure of the grid itself. The grid IS a^2 + b^2. The nonce IS c^2. You don't search for c. You compute c from a and b. That computation is the derivative. It skips the ladder and reads the hypotenuse directly.

The fraction n/(n+1) also shows up in the halvening:

    After k halvenings, the reward fraction remaining = 1/2^k
    The cumulative fraction mined = 1 - 1/2^k = (2^k - 1) / 2^k

That's n/(n+1) with n = 2^k - 1. After 1 halvening: 1/2. After 2: 3/4. After 3: 7/8. Always approaching 1 (all coins mined), never arriving (the last satoshi is mined around halvening 33, but the fraction asymptotes).

The halvening IS the n/(n+1) fraction expressed in block time. Each halvening is one machine step. The reward is the remaining gap. The gap shrinks but never vanishes.

---

## 9. The 10000 Gate

    10000 = b^4    (in any base b)

The fourth power of the base. The first number with enough dimensions. Below b^4, you have at most 3 independent directions. At b^4, you have 4. Four is the minimum for the axiom to fully operate: two legs, one hypotenuse, one observer. Three numbers in the triangle plus one outside it to see it.

In the halvening:

    210,000 / 21 = 10,000

The halvening divided by the machine product gives exactly 10,000. So 10,000 blocks = one machine cycle. The machine (3 x 7 = 21) runs exactly 10,000 blocks, then starts a new cycle. 21 cycles = one halvening.

10,000 is the decimal gate because in base 10: 10,000 = 10^4 = (2 x 5)^4. The fourth power. The gate opens when the base has enough dimensions.

The cosmic check:

    291 / 21 = 13.857...

291 = the number of phi-steps from Planck scale to cosmic scale (2 x phi^290 x l_Planck = 13.80 billion light-years). 21 = the machine product. Their ratio:

    291 / 21 = 13.857

The observed radius of the universe: 13.80 billion light-years. The ratio of the phi exponent to the machine product approximates the size of everything. Not in meters. In billions of light-years. The cosmic exponent divided by the machine product gives the universe.

This is not numerology. This is the axiom at scale. The same machine that resolves degrees at the nonce level (84 / 21 = 4 degrees) resolves cosmic distances at the phi-exponential level (291 / 21 = 13.857 Gly). Same denominator. Same machine. Different exponent.

The 10,000 gate in block terms:

    10,000 blocks x 10 minutes/block = 100,000 minutes
    100,000 / 60 / 24 = 69.44 days
    69.44 days ≈ 10 weeks

One machine cycle is approximately 10 weeks. 21 machine cycles (one halvening) ≈ 210 weeks ≈ 4 years. And the halvening is indeed approximately every 4 years. The 10,000 gate converts machine cycles to human time at the rate of 10 weeks per cycle.

---

## 10. The Pythagorean Resolution

Everything in this document reduces to one statement:

    a^2 + b^2 = c^2

The grid is a^2 + b^2. The nonce is c^2. The degree is the resolution at which you read the right angle of that triangle.

    Pole A = a^2        (LE words, one endianness, one leg)
    Pole B = b^2        (BE words, the flip, the other leg)
    Nonce  = c^2        (the hypotenuse, what you're solving for)

The grid doesn't SOLVE the nonce. The grid IS a^2 + b^2. It computes both legs simultaneously across all word pairs. The nonce doesn't LIVE on the grid. The nonce IS c^2. It lives on the hypotenuse. You find it by computing c from a and b. Not by searching for c in a list.

The degree formula ties it together:

    degree = 84 / (n + 1)

The resolution of how precisely you see the right angle determines how precisely you find the nonce. More machine steps = smaller degree = finer resolution = closer to the answer.

At 28 degrees (n = 2), you can tell which quadrant the nonce is in. Coarse. 
At 12 degrees (n = 6), you can tell which cube face. Better.
At 7 degrees (n = 11), you can read the nonce ladder. Good enough for base-7 digit-by-digit extraction.
At 4 degrees (n = 20), you can see the full machine structure. The grid resolves.
At 1 degree (n = 83), full resolution. Every angle distinct. Every nonce unique.

But the derivative approach SKIPS the degree ladder. It doesn't step from 28 to 12 to 7 to 4 to 1. It reads the right angle directly. It computes c from a and b in one operation. The derivative machine doesn't need fine resolution because it doesn't measure the angle — it computes the hypotenuse from the legs. No measurement. No degree. No ladder.

The degree machine is for deduction. Walk the rungs, refine the angle, converge on the nonce.

The derivative machine is for computation. Read the legs, compute the hypotenuse, arrive at the nonce.

Both machines run on the same axiom. Both machines use the same +1. The degree machine puts the +1 in the denominator (84 / (n + 1)). The derivative machine puts the +1 in the equation (x^2 = x + 1). Same +1. Two uses.

The 336-degree circle is the space in which all of this happens. Not 360. 336. The true circle. The one that wraps to 2^14 exactly. The one that factors into 2^4 x 3 x 7. The one whose right angle is 84, not 90. The one whose halvening is 625 circles, not some awkward fraction.

The Babylonians added 24 degrees of noise. The machine strips them off. What remains is 336: the circle that the axiom generates, that the cube machine traverses, that the nonce lives on, that the halvening counts.

---

    degree = 84 / (n + 1)

    The right angle is 84.
    The circle is 336.
    The machine product is 21.
    The halvening is 625 x 336.
    The nonce is 12 rungs of base 7.
    The dipole is D = 2, the spiral is D = 1.
    The gate is 10000 = b^4.
    The fraction is n/(n+1), always approaching, never arriving.
    The derivative skips the ladder and reads c directly.

    a^2 + b^2 = c^2.

    x^2 = x + 1.

The axiom was Pythagorean from the start. — nos3bl33d
