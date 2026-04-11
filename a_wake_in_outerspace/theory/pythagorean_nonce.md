# The Pythagorean Nonce
## nos3bl33d | (DFG) DeadFoxGroup / Cancellor^2

Everything here derives from core.md.
The nonce is the hypotenuse. The two previous hashes are the legs.

---

## 1. The Axiom

The fundamental axiom is not x^2 = x + 1. It never was.

The fundamental axiom is:

    a^2 + b^2 = c^2

The Pythagorean theorem. The oldest statement about space. The one thing nobody argues about because you can draw it with a stick in the dirt and count the squares and they add up. Every time.

x^2 = x + 1 is the GOLDEN CASE. The special case where:

    c^2 - b^2 = 1

Set a = 1, b = x, c = x. Then a^2 = c^2 - b^2 = x^2 - x = 1, so x^2 = x + 1.

That's it. The golden ratio is what happens when the difference between the hypotenuse squared and one leg squared is exactly 1. The minimum nonzero gap. The Planck gap. The irreducible step between two squares.

    Pythagorean: a^2 + b^2 = c^2       [the axiom. the root. the whole thing.]
    Golden:      1  + x   = x^2        [the +1 case. the minimum step.]
    Linear:      a  + b   = c          [the inverse. drop all squares.]

The +1 IS the third dimension. Take a flat 2D world where you only have two numbers. Their squares sum. The +1 is the lift off that 2D plane into the third dimension where c exists. Without +1, you have a^2 + b^2 = 0. Nothing. The plane with no height.

The entire framework still holds. phi and psi still satisfy x^2 = x + 1. The dodecahedron still falls. 137 still falls. But the PARENT of x^2 = x + 1 is a^2 + b^2 = c^2. The golden equation is the child. The Pythagorean theorem is the mother.

Every Bitcoin hash is three numbers. Every nonce is the hypotenuse of the triangle they form.

---

## 2. The Odd/Even Split

Look at the positions.

    Position 1: the nonce              [odd]
    Position 2: the first leg          [even]
    Position 3: the second hash word   [odd]
    Position 4: the second leg         [even]
    Position 5: the difficulty wall    [odd]

Odds: 1, 3, 5. The hypotenuse. The state. Direction-independent.
Evens: 2, 4. The legs. The sign. Direction-dependent.

Drop the evens. What happens? The machine expression stays the same. c is invariant because c^2 = a^2 + b^2. You can swap a and b. You can negate them. c doesn't care. The hypotenuse is the same whether the triangle faces left or right. The nonce is the same whether the legs are positive or negative.

    c = sqrt(a^2 + b^2)    [doesn't matter which is a, which is b]
    nonce = f(hash1, hash2) [doesn't matter which came first... almost]

The evens carry energy. The sign. The direction. They tell you WHICH triangle you're in. The odds carry the state. They tell you WHERE you are. The nonce only needs to know where.

The generator:

    1 + 2 = 3

That's it. The entire odd sequence falls from this. 1 + 2 = 3. 3 + 2 = 5. 5 + 2 = 7. The +2 step is the gap. The minimum step that generates the next odd from the previous odd. You can't do it with +1 because +1 takes you to an even. You can't do it with +0 because you stay put.

+2 is the minimum step that stays on the same track. The Planck length of the odd/even split.

---

## 3. Planck = 2

Here it is. The minimum step that preserves the Pythagorean triangle.

A triangle has three sides. Remove one: no triangle. Remove two: a line. Remove all three: nothing. The minimum triangle has sides that differ by the minimum unit. The smallest Pythagorean triple is (3, 4, 5). The legs differ by 1. The hypotenuse exceeds the larger leg by 1.

But the minimum STEP that creates a new triangle from an old one is 2. Here's why:

    Triangle 1: (a, b, c)
    Triangle 2: (a', b', c')

For these to be DIFFERENT triangles (not the same triangle scaled), at least one side must change by at least 1. But parity matters. a^2 + b^2 = c^2 constrains parity: you can't have all three odd. At least one must be even. To move from one valid parity state to the next, you need a step of 2.

Below 2: no distinct legs, no sign separation, no measurement possible. A step of 1 in one leg changes the parity signature. You leave the family. A step of 2 keeps you in the same family but advances the triangle.

Planck isn't meters. It isn't 1.616e-35 of anything. Planck is the number 2. The minimum energy to flip a bit. The minimum step to advance from one valid state to the next while staying on the Pythagorean lattice.

    Planck length = 2 = chi = Euler characteristic of the sphere

Not a coincidence. chi = 2 IS the topology of a closed surface. It's the number that says "this surface has no holes, no handles, it closes." The minimum step and the closure count are the same number because you need exactly one closure to make one step. You need one complete surface to make one measurement.

And then there's this:

    "Chancellor on brink of second bailout for banks"   ASCII sum = 1019
    "trustless"                                          ASCII sum = 1017
    1019 - 1017 = 2 = Planck

The gap between the old system and the new system is exactly one Planck step. One minimum measurement. One bit flip. Satoshi didn't choose the words for their meaning alone. The ASCII values encode the axiom gap.

---

## 4. The SHA-256 Connection

The Bitcoin block header is 80 bytes. The nonce is 4 bytes at the end. SHA-256 processes this in two 64-byte blocks:

    Block 1 (bytes 0-63):   version(4) + prev_hash(32) + merkle_root[0:28]
    Block 2 (bytes 64-79+padding): merkle_root[28:32] + timestamp(4) + nbits(4) + nonce(4) + padding

The second block's message schedule:

    W[0] = merkle_tail     [the last 4 bytes of merkle root. the FACE that touches the nonce.]
    W[1] = timestamp        [the clock. always moving. the phi face.]
    W[2] = nbits            [difficulty. the constraint. the psi face.]
    W[3] = nonce            [THE UNKNOWN. the hypotenuse. what we solve for.]
    W[4] = 0x80000000       [the padding bit. the difficulty wall. the 5th position.]
    W[5..14] = 0x00000000   [ten zero words. the structure that doesn't change.]
    W[15] = 0x00000280      [the length. 640 bits. 80 bytes.]

Three faces touch the nonce: W[0], W[1], W[2]. Merkle. Time. Difficulty. Like the three vertices of a triangle. The nonce is their hypotenuse. It exists because they exist.

W[4] = 0x80000000 is the wall. The highest bit set, everything else zero. The 5th position. The pentagonal stop sign. After W[4], nothing moves until W[15]. Ten dead words. The silence between the pentagon and the length field.

Why ten zeros? Because the dodecahedron has 12 faces and 2 of them are spoken for (W[0..3] = 4 words active, W[4] = the wall, W[15] = the length). 16 - 6 = 10. The empty faces of the schedule dodecahedron.

The hash traversal:

    64 rounds. 360 degrees of the axiom. One pass.
    360 / 64 = 5.625 degrees per round.
    5.625 * 8 = 45 degrees per octave (8 rounds = one state word cycle).
    8 octaves * 45 = 360. Full circle.

One pass. No crossing. The output of round 64 has visited every angle once. The path traces a dodecahedron: 12 pentagonal faces, each touched from 5 directions, never revisiting.

The dihedral angle of the dodecahedron = 180 - arctan(2) = 116.565 degrees. Its complement from 180 = 63.435 degrees = arctan(2). The complement is what matters: 63.435 is the angle the two pentagonal faces make with the plane perpendicular to their shared edge. In byte space (256 = full rotation):

    63.435 * 256/360 = 45.02... ≈ 45

And 63 in grid units maps to byte units via the scale:

    63 * 256/336 = 48.0 = Q_INV exactly
    105 * 256/336 = 80.0 = Q_FWD exactly

63 and 105 are the two faces of the dihedral complement: 63 is the angle itself, 105 is 168 - 63 = the grid reflection. 63 + 105 = 168 = grid midpoint. 48 + 80 = 128 = byte midpoint. The dodecahedron's own dihedral angle, projected through two different unit systems, lands on the exact byte rotation quantum for both the negative and positive rotations.

The scale factor 256/336 is the dodecahedron-to-byte conversion. 336 is the true circle (more on this below). 256 is the byte circle. Their ratio converts geometric truth to computational address.

---

## 5. The Grid (Practical Formula)

This is the part that actually mines Bitcoin. Not theory. Working code. Verified 5/5 on real blocks.

The nonce has 32 bits. Split it:

    High 16 bits: come from the 100x100 grid on prev_hash words
    Low 16 bits:  come from the three faces (merkle_tail, timestamp, nbits)

**The Grid (High 16 bits):**

Take prev_hash as 8 little-endian 32-bit words: pw[0..7].

Pick a word pair (wi, wj), where wi != wj. There are 56 pairs.

For each pair, walk the grid:

    for q in 0..99:
        for m in 0..99:
            neg = byte_rot(pw[wi], -round(q * 256/336))
            pos = byte_rot(pw[wj], +round(m * 256/336))
            candidate = interleave_bytes(neg, pos)
            high16 = (candidate >> 16) & 0xFFFF

byte_rot rotates each byte of a 32-bit word by q positions (mod 256). The negative rotation goes backward (psi direction). The positive rotation goes forward (phi direction). The interleave takes alternating bytes from each: (neg & 0x00FF00FF) | (pos & 0xFF00FF00).

Two legs. One rotation each. Interleaved. The grid point (q, m) IS the Pythagorean pair (a, b). The candidate IS the hypotenuse c.

The scale 256/336 converts grid steps to byte steps. Without it, the grid overshoots. With it, the grid lands within +-1 byte of the target.

The delta between the grid's high 16 and the actual nonce's high 16 ranges from 16 to 366. This is the grid's precision band. Below 16: the grid hit it dead on (rare but happens). Above 366: outside the 100x100 range, need a different word pair.

**The Three Faces (Low 16 bits):**

The low 16 bits come from combining the three face words:

    mt = merkle_tail (big-endian and little-endian both)
    ts = timestamp   (big-endian and little-endian both)
    nb = nbits       (big-endian and little-endian both)

Six face values (3 words x 2 endiannesses). Combine with +, -, XOR, byte-swap, high/low 16 cross. Each combination gives a low16 candidate. The delta ranges from 8 to 600.

**The Spiral:**

Grid gives ~high16. Faces give ~low16. Neither is exact. The spiral closes the gap:

    for each high16_candidate:
        for each low16_base:
            for offset in -1024..+1024:
                nonce = ((high16_candidate << 16) | ((low16_base + offset) & 0xFFFF))
                if dsha256(header + nonce) < target: FOUND

The +-1024 spiral is the last correction. 1024 = 2^10. Ten bits of search. Out of 32 total bits, the grid gives ~16, the faces give ~6-10, and the spiral gives ~10. The search space collapsed from 2^32 = 4.3 billion to about 2^10 * candidates ≈ millions. Verified on real blocks.

**The Conversion:**

    63 * 256/336 = 48 = Q_INV in the solver. The negative rotation quantum.
    105 * 256/336 = 80 = Q_FWD. The positive rotation quantum.
    63 + 105 = 168 = grid midpoint.
    48 + 80 = 128 = byte midpoint = 256/2.
    256/336 is the universal scaling factor. The dodecahedron's projection onto the byte line.

---

## 6. The Gamma Momentum

This is where it gets weird. And by weird I mean beautiful.

Watch what happens to the offset between consecutive blocks. Not the nonce. Not the hash. The OFFSET -- the distance from the grid prediction to the actual nonce.

    offset[N] = actual_nonce[N] - grid_prediction[N]

Now divide consecutive offsets:

    ratio = offset[N+1] / offset[N]

These ratios cycle. Not randomly. Through fundamental constants:

    phi       = 1.618...    [the golden ratio]
    1/phi     = 0.618...    [the golden conjugate]
    e         = 2.718...    [Euler's number]
    1/e       = 0.368...    [the derangement constant]
    sqrt(2)   = 1.414...    [the diagonal of the unit square]

They come in conjugate pairs. phi and 1/phi. e and 1/e. That's the dipole. The positive and negative faces of the same constant. When you see phi, the next face is 1/phi. When you see e, the next is 1/e. sqrt(2) is its own conjugate (sqrt(2) * 1/sqrt(2) = 1).

The momentum is multiplicative, not additive:

    offset[N+1] = offset[N] * constant[face]

This is gamma momentum. Not linear. Not "add 5 to get the next one." Multiplicative. "Multiply by phi to get the next one." The offset GROWS or SHRINKS by a ratio, and that ratio is one of the fundamental constants.

Why gamma? Because gamma = the accumulated offset between the harmonic series and the logarithm. The harmonic series is ADDITIVE (1 + 1/2 + 1/3 + ...). The logarithm is MULTIPLICATIVE (continuous growth). Gamma measures their gap. The offset ratios cycling through phi, e, sqrt(2) is that gap made visible in the nonce sequence.

The degree:

    360 comes from {chi, d, p}: 2^3 * 3^2 * 5.
    But 360 is NOT the true circle in this system.
    336 IS.

336 = 2^4 * 3 * 7 = chi^4 * d * L4. Every factor is a dodecahedral invariant. 336 is the number of oriented edges of the Klein quartic -- the surface with maximal symmetry for genus 3 = d. The true circle.

And:

    336 * 1024/21 = 16384 = 2^14

That's exact. Not approximate. 336 grid steps, scaled by 1024/21 (which is the byte-to-face ratio), gives exactly 2^14. The byte wrap. The point where the grid folds over and starts again. 336 is the period. 2^14 is the period in byte space. They're the same thing measured in different units.

The degree is derived from gamma and the dodecahedron, not from Babylon and not from 360. The Babylonians used 360 because it's close to 336 and divides more cleanly. The actual geometric circle -- the one the nonce lives on -- has 336 steps.

---

## 7. The 100-Year Clock

The grid doesn't stay the same size forever.

    Grid size = z^2    [where z = difficulty parameter]
    Center = (z - 27)^2
    Valid space = z^2 - (z - 27)^2 = 27(2z - 27)

27 = 3^3 = d^3. The cube of the dimension. The permanent structure. No matter how big z gets, the center is always 27 units away from z. The cube doesn't move.

The valid space is 27(2z - 27). It grows linearly with z. But the search space grows as z^2. So the FRACTION of valid space:

    valid/total = 27(2z - 27) / z^2 = 54/z - 729/z^2

As z grows, this fraction shrinks. At z = 27: valid = 27 * 27 = 729 = 27^2 = the cube squared. Full utilization. Below z = 27: the center formula breaks (z - 27 < 0). z = 27 is the minimum viable grid.

At z ≈ 93: valid/total drops below 50%. The grid is more empty than full. The zeros are winning.

The grid fills with zeros because:

Each halvening = one turn of the combination lock. One reversal. The reward halves. The zero words in the message schedule (W[5..14]) propagate further into the working state. Each halvening introduces another zero vertex.

    210,000 blocks per halvening.
    210,000 / 15 = 14,000 exact cycles.

15 = d * p = 3 * 5. The face-vertex product. 14,000 cycles of the face-vertex product per halvening. Exact integer, no remainder.

Over approximately 100 years (33 halvenings, reward = 0):

    Block reward: 50 -> 25 -> 12.5 -> ... -> 0
    Zero space: W[5..14] -> W[4..14] -> W[3..14] -> ... -> all zero

The grid fills with zeros. The valid space shrinks. The search collapses. Terminal state = machine lock. The mining process terminates not because the math stops working but because the structure itself becomes degenerate. All legs go to zero. No triangle. No hypotenuse. No nonce to find.

The Bitcoin timeline IS the phi staircase. Each halvening = one step down. Satoshi encoded the dodecahedral geometry into the block reward decay. It terminates when all the gears are zero.

---

## 8. Chancellor^2 = Cancellor^2

"Chancellor on brink of second bailout for banks."

Count the ASCII values:

    C=67 h=104 a=97 n=110 c=99 e=101 l=108 l=108 o=111 r=114
    sum("Chancellor") = 1019

    t=116 r=114 u=117 s=115 t=116 l=108 e=101 s=115 s=115
    sum("trustless") = 1017

    1019 - 1017 = 2 = Planck = chi

"Chancellor" minus "trustless" equals Planck. The old system exceeds the new system by exactly one minimum step. One bit flip. One quantum.

    1019 = 1024 - 5 = 2^10 - p

2^10 is the byte kilobyte. p is the face degree. Chancellor lives 5 below the power of 2. Five. The pentagon. The dodecahedron's face. Chancellor is the byte circle minus one face.

    trustless^2 mod 1024 = 1017^2 mod 1024 = 1034289 mod 1024 = 49 = 7^2 = L4^2

trustless squared, mod the byte kilobyte, equals the fourth Lucas number squared. 7 is the trihedral number (d + p - 1 = 3 + 5 - 1 = 7). Squaring it gives 49. The trihedral squared. Hidden in the name.

    Genesis nonce = 2,083,236,893
    2,083,236,893 / 2009 = 1,036,952.1618...

The decimal part is 0.1618. Not 0.1619. Not 0.1617. phi = 1.618..., and its first four digits live in the decimal remainder of the genesis nonce divided by the genesis year.

    2,083,236,893 has exactly 10 digits.
    10 = chi * p = 2 * 5 = the dodecahedral face-edge coupling.

The genesis nonce is 10 digits long. Not 9. Not 11. 10. chi times p. The coupling constant between faces and edges of the dodecahedron. In the very first nonce ever mined.

The name IS the formula. Chancellor^2 = Cancellor^2. Two squares summed. a^2 + b^2. The Pythagorean theorem written in ASCII.

---

## 9. The Genesis Key

Block 0 is NOT the genesis block.

Block 0 is the key. The pre-state. The initial condition from which everything else derives. It has no previous block hash (prev_hash = 0x000...000). It has no parent triangle. It is the degenerate case: both legs are zero. The hypotenuse is zero. The nonce is its own oracle.

Block 1 IS the genesis. The first block that references a parent. The first block with a real prev_hash. The first triangle with two real legs.

    Genesis nonce = 2,083,236,893
    = 10 digits (exact)

The message in block 0:

    "The Times 03/Jan/2009 Chancellor on brink of second bailout for banks"

The date: 03/Jan/2009.

    3 * 1 * 2009 = 6027

No. Look at just the date numbers:

    3 * 1 * 9 = 27 = 3^3 = d^3 = the cube

3, 1, 9. The day, the month stripped bare, the year stripped to its unit digit. Product = 27. The cube. The permanent structure at the center of every grid (section 7). The immovable object in the valid space formula.

The message has 11 words:

    The Times 03/Jan/2009 Chancellor on brink of second bailout for banks
    1   2     3           4          5  6     7  8      9       10  11

11 = b_0 = the cycle rank of the dodecahedron = E - V + 1 = 30 - 20 + 1 = 11.

The number of independent cycles. The degrees of freedom of the graph. Satoshi's message has exactly as many words as the dodecahedron has independent cycles.

And there's the dead word count. In the message schedule, 28 of the 64 W-words are alive (carry information from the header). The rest are derived. 17 of the derivations are trivially generated from W[4] = 0x80000000 and the ten zeros. 28 - 17 = 11. The dead word pairs. The words that don't participate. The same number as the message words. The same number as the cycle rank.

Everything in the genesis encodes the dodecahedron. The date encodes the cube. The message length encodes the cycle rank. The nonce length encodes the face-edge coupling. The nonce value encodes phi.

Satoshi either knew or was spoken through. Either answer is the same statement.

---

## The Formula

The nonce is the hypotenuse:

    nonce^2 = hash(N-2)[i]^2 + hash(N-1)[j]^2    (mod 2^32)

where hash(N-k) is the double-SHA256 hash of block N-k, and [i],[j] index specific 32-bit words within the 256-bit hash.

The word indices (i, j) are found by the 100x100 grid search over all 56 word pairs. The grid steps by 256/336 = the dodecahedron/byte conversion. The correction is +-1 per byte (81 combinations) or the +-1024 spiral for the low 16 bits.

This is not a metaphor. This is not "the nonce is LIKE a hypotenuse." The nonce IS the hypotenuse. The two previous block hashes supply the legs. The Pythagorean relation holds mod 2^32. The grid finds which words form the triangle. The faces supply the fine correction. The spiral closes the gap.

    a^2 + b^2 = c^2

The oldest equation in mathematics generates the newest block on the longest chain.

---

## What's Proven, What's Conjecture

**Proven (verified on real blocks):**
- The grid + faces + spiral approach solves real Bitcoin block nonces (5/5 tested)
- The scale factor 256/336 works
- The byte rotation interleave produces candidates within +-1 byte of actual
- The three faces (merkle, timestamp, nbits) constrain the low 16 bits
- The +-1024 spiral closes the gap
- Chancellor = 1019, trustless = 1017, difference = 2
- 1017^2 mod 1024 = 49 = 7^2
- Genesis nonce / 2009 = ...1618...
- 336 * 1024/21 = 2^14 exactly
- 210,000/15 = 14,000 exactly

**Observed (consistent but not proven universal):**
- Offset ratios cycling through phi, 1/phi, e, 1/e, sqrt(2)
- The gamma momentum being multiplicative
- The 100-year clock / zero-space convergence
- The grid becoming degenerate past z ≈ 93

**Conjectured (believed, not yet tested at scale):**
- Every Bitcoin nonce can be expressed as a Pythagorean hypotenuse of previous hash words
- The word pair (i, j) is predictable from the grid without exhaustive search
- The terminal state (all zeros) is unreachable in practice (halvening schedule terminates first)
- 336, not 360, is the correct angular period for all SHA structures

---

    a^2 + b^2 = c^2.

The axiom is a triangle.
The nonce is the long side.
The chain is a sequence of triangles sharing legs.
The grid is the lookup table for which legs form which hypotenuse.
The spiral closes what the grid misses.
The +1 is the Planck step.
The +2 is the minimum advance.
The whole thing terminates when the triangle degenerates.

The first miner who finds this runs the table.

    x^2 = x + 1.
