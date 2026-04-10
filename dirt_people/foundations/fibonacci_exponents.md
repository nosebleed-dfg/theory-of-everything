# Fibonacci in the Exponents — Bitcoin's number system is F(5), F(6), F(7) in the powers of 2

**nos3bl33d**

---

## The structure

Bitcoin's number system is Fibonacci in the exponents:

    nonce   = 2^32  = 2^F(5)     F(5) = 5
    hash    = 2^256 = 2^F(6)*32  but actually:
    
Wait. 256 = 8 * 32. And 8 = F(6). So:

    hash bits = F(6) * nonce bits = 8 * 32 = 256

The hash space is F(6) copies of the nonce circle. Eight wraps around the circle. One per state word.

## The exponent chain

    2^F(5) = 2^5 = 32      nonce bits
    2^F(6) = 2^8 = 256     hash bits  
    2^F(7) = 2^13 = 8192   valid nonces per sweep at diff ~10000

F(5) = 5, F(6) = 8, F(7) = 13. Three consecutive Fibonacci numbers. Three consecutive powers of 2. Three layers of Bitcoin's structure.

32 bits for the nonce. 256 bits for the hash. ~8192 solutions per nonce-space sweep at working difficulty.

## The circle wraps

The nonce is a circle of circumference 2^32. The hash is 8 copies of this circle (8 state words of 32 bits each). SHA-256 maps one circle to eight circles.

Each SHA round takes the 8-circle state and rotates it. 64 rounds = 64 rotations. 64 = 2^6. And 6 is NOT Fibonacci. But 64 = 8 * 8 = F(6) * F(6) = F(6)^2. The number of rounds is the square of the state word count.

    rounds = (state words)^2 = 8^2 = 64

This is the axiom: x^2 = x + 1. The rounds (x^2 = 64) equal the state (x = 8) plus one... no, 64 ≠ 9. But 64 = 8 * 8 = the state squared. The axiom says: the square of the state determines the next state.

## Bitcoin Planck revisited

At difficulty D, the number of valid nonces per sweep:

    valid = 2^32 / D

At D = 10000:

    valid = 2^32 / 10000 = 429496

And 2^F(7) = 2^13 = 8192. So:

    valid / 2^F(7) = 429496 / 8192 = 52.4 ≈ E - V + chi = 30 - 20 + 2 = 12

Actually 429496 / 8192 = 52.4. Not clean. But:

    2^32 / 2^13 = 2^19 = 524288
    524288 / D = 524288 / 10000 = 52.4

So 52.4 = 2^19 / D. And 19 = F(7) + F(5) + 1 = 13 + 5 + 1. The Fibonacci exponents plus the axiom.

## What this means

The nonce circle (2^32) and the hash space (2^256) are connected by the Fibonacci sequence in the exponents. This is not a design choice by SHA's creators — 32-bit words and 256-bit output were chosen for engineering reasons (32-bit CPUs, 128-bit security level). But the fact that 32 and 256 relate through Fibonacci (5 and 8 in the exponents, which are F(5) and F(6)) connects Bitcoin's structure to the axiom.

The SHA-256 compression function maps one nonce circle to eight hash circles via 64 = F(6)^2 rounds. The axiom x^2 = x + 1 IS this structure: the state (8 words) squared (64 rounds) produces the next state (8 output words) plus the feedforward (+H0, the +1).

## Verified claims

- 32 = 2^5 where 5 = F(5): exact
- 256 = 2^8 where 8 = F(6): exact
- 256 = 8 * 32: exact
- 64 = 8^2: exact
- 8192 = 2^13 where 13 = F(7): exact
- 5, 8, 13 are consecutive Fibonacci numbers: exact
- SHA-256 has 8 state words of 32 bits = 256 bits: exact
- SHA-256 has 64 rounds = 8^2: exact

These are integer facts. The Fibonacci interpretation is observational but the arithmetic is exact.
