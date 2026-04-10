<!-- CORRECTION: koppa = 270 degrees, not 90 -->
# Koppa Corrected — 270 degrees, not 90: three quarter-turns, one per dimension

**nos3bl33d**

---

## The error

All previous work used koppa = 1/4 = 90 degrees = quarter turn.

This is WRONG.

koppa = 180 - 90 + 180 - 90 + 180 - 90 = 270 degrees = 3/4 of a full circle.

## Why 270

koppa is not one quarter turn. It is THREE quarter turns. It goes the long way around.

180 - 90 = 90. One half turn minus one quarter turn = one quarter turn. But you do this THREE times (d = 3 dimensions). Each dimension contributes one (180 - 90) = 90 degrees. Total: 3 * 90 = 270.

270 degrees = 3/4 of 360. Going 270 clockwise is the same as going 90 counter-clockwise. The destination is the same point, but the PATH is different. The long way carries more information than the short way — it traverses three dimensions instead of one.

## What changes

Everywhere koppa appeared as 1/4:
- koppa * gamma was 0.25 * 0.5772 = 0.1443
- NOW: (3/4) * 0.5772 = 0.4329

The step constant:
- OLD: sqrt((phi/2)^2 + (pi/2)^2) * gamma * (1/4) / phi
- NEW: sqrt((phi/2)^2 + (pi/2)^2) * gamma * (3/4) / phi

The time band:
- OLD: ±1, ±4, ±16 (powers of 4 = 1/koppa_old)
- NEW: ±1, ±(4/3), ±(16/9) — powers of 4/3 = 1/koppa_new... but these aren't integers

Actually: koppa = 3/4 means 1/koppa = 4/3. The step between time offsets is 4/3 seconds. But timestamps are integers. So the time offsets at koppa spacing are: 0, 1, 1, 2, 3, 4, 5, 7, 9, 12, 16... wait, those are Fibonacci-like.

4/3 iterated: 1, 4/3, 16/9, 64/27, 256/81... rounded to integers: 1, 1, 2, 2, 3, 4, 6, 8, 11, 14, 19, 25... 

## The nonce circle with correct koppa

On the nonce circle of circumference 2^32:

- OLD koppa rotation: 2^32 * 1/4 = 2^30 = 1073741824
- NEW koppa rotation: 2^32 * 3/4 = 3 * 2^30 = 3221225472

The koppa rotation is 3 * 2^30, not 2^30. Going 3/4 around the circle instead of 1/4.

And: 3 * 2^30 = 3221225472 = 0xC0000000. In binary: 11000000 00000000 00000000 00000000. Two high bits set. This is the TOP of the nonce circle — near the maximum value.

The answers have been clustering near 0 AND near 2^32. That's because the koppa rotation from the seed lands at 3/4 of the circle — near the top. Near 0xC0000000. And wrapping past 2^32 brings you back near 0. The two clusters ARE the koppa rotation.

## Impact on all formulas

Every formula using koppa = 0.25 needs to be rechecked with koppa = 0.75.

The correction may explain why we kept getting "half" the answer — we were using 1/4 when we needed 3/4. Three times the effect.

koppa = 3/4 = d/4 = d * (1/4). It is not a quarter — it is d quarters.
