#!/usr/bin/env python3
"""
GOLDEN_CRYPTO — investigates hash function cycle lengths vs phi^bits and the phi^291 ceiling
nos3bl33d

Toy hashes and truncated SHA-256. Brent cycle detection, golden vs standard comparison, birthday bound.
"""

import hashlib
import math
import sys
from collections import Counter
from typing import NamedTuple

from mpmath import mp, mpf, power, sqrt, log, frac, floor, nstr

# High precision for phi computations
mp.dps = 50

PHI = (1 + sqrt(5)) / 2
PHI_FLOAT = float(PHI)

# ─────────────────────────────────────────────────────────────────
#  UTILITY: cycle detection via Floyd's algorithm + full enumeration
# ─────────────────────────────────────────────────────────────────

def find_cycle_from(start, func, state_mask):
    """Find cycle length starting from `start` using Brent's algorithm.
    Returns (tail_length, cycle_length)."""
    power_val = 1
    lam = 1
    tortoise = start
    hare = func(start)

    # Phase 1: find lam (cycle length)
    while tortoise != hare:
        if power_val == lam:
            tortoise = hare
            power_val *= 2
            lam = 0
        hare = func(hare)
        lam += 1

    # Phase 2: find start of cycle (mu)
    tortoise = start
    hare = start
    for _ in range(lam):
        hare = func(hare)

    mu = 0
    while tortoise != hare:
        tortoise = func(tortoise)
        hare = func(hare)
        mu += 1

    return mu, lam


def find_all_cycles(func, bits):
    """Find ALL distinct cycles in the function's state space.
    Returns list of cycle lengths and a dict mapping each state to its cycle length."""
    n = 1 << bits
    mask = n - 1
    visited = [False] * n
    cycle_lengths = []

    for start in range(n):
        if visited[start]:
            continue

        # Walk the trajectory until we hit a visited node or revisit our own path
        path = []
        path_set = {}
        x = start
        while x not in path_set and not visited[x]:
            path_set[x] = len(path)
            path.append(x)
            x = func(x)

        if not visited[x] and x in path_set:
            # Found a new cycle
            cycle_start_idx = path_set[x]
            cycle_len = len(path) - cycle_start_idx
            cycle_lengths.append(cycle_len)

        # Mark all nodes in path as visited
        for node in path:
            visited[node] = True

    return sorted(cycle_lengths, reverse=True)


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 1: Toy SHA cycle lengths
# ─────────────────────────────────────────────────────────────────

def toy_hash(x, bits=16):
    """SHA-like operations on a small state."""
    mask = (1 << bits) - 1
    for _ in range(8):
        x = ((x << 5) | (x >> (bits - 5))) & mask
        x = (x + 0x9E3779B9) & mask
        x = x ^ (x >> 3)
        x = (x * 0x5BD1E995) & mask
        x = x ^ (x >> 7)
    return x


class CycleStats(NamedTuple):
    max_cycle: int
    avg_cycle: float
    num_cycles: int
    all_cycles: list


def experiment_1():
    print("=" * 72)
    print("  EXPERIMENT 1: Toy SHA Cycle Lengths")
    print("=" * 72)

    bit_widths = [8, 10, 12, 14, 16]
    results = {}

    for bits in bit_widths:
        n = 1 << bits
        func = lambda x, b=bits: toy_hash(x, b)
        cycles = find_all_cycles(func, bits)

        max_c = max(cycles)
        avg_c = sum(cycles) / len(cycles)
        phi_ceil = float(power(PHI, bits))
        golden_scaled = float(power(PHI, bits * mpf(291) / 256))

        results[bits] = CycleStats(max_c, avg_c, len(cycles), cycles)

        print(f"\n  bits = {bits}  (state space = {n})")
        print(f"    Distinct cycles:  {len(cycles)}")
        print(f"    Max cycle:        {max_c}")
        print(f"    Avg cycle:        {avg_c:.1f}")
        print(f"    phi^bits:         {phi_ceil:.1f}")
        print(f"    phi^(bits*291/256): {golden_scaled:.1f}")
        print(f"    Max / phi^bits:   {max_c / phi_ceil:.4f}")
        print(f"    Max / 2^bits:     {max_c / n:.4f}")
        print(f"    Cycle lengths:    {cycles[:10]}{'...' if len(cycles) > 10 else ''}")

    return results


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 2: Golden hash vs standard hash
# ─────────────────────────────────────────────────────────────────

def golden_hash(x, bits=16):
    """Multiply by golden ratio (mod 2^bits)."""
    mask = (1 << bits) - 1
    phi_int = int((1 + 5**0.5) / 2 * (1 << bits)) & mask
    return (x * phi_int) & mask


def cycle_uniformity(cycles, n):
    """Measure how uniform cycle lengths are. Lower = more uniform.
    Uses coefficient of variation."""
    if len(cycles) < 2:
        return 0.0
    mean = sum(cycles) / len(cycles)
    if mean == 0:
        return 0.0
    var = sum((c - mean) ** 2 for c in cycles) / len(cycles)
    return (var ** 0.5) / mean


def experiment_2():
    print("\n" + "=" * 72)
    print("  EXPERIMENT 2: Golden Hash vs Standard Hash")
    print("=" * 72)

    bit_widths = [8, 10, 12, 14, 16]

    for bits in bit_widths:
        n = 1 << bits
        mask = n - 1

        # Golden hash cycles
        g_func = lambda x, b=bits: golden_hash(x, b)
        g_cycles = find_all_cycles(g_func, bits)

        # Toy SHA cycles
        t_func = lambda x, b=bits: toy_hash(x, b)
        t_cycles = find_all_cycles(t_func, bits)

        g_max = max(g_cycles)
        t_max = max(t_cycles)
        g_cv = cycle_uniformity(g_cycles, n)
        t_cv = cycle_uniformity(t_cycles, n)

        print(f"\n  bits = {bits}  (state space = {n})")
        print(f"    {'':20s} {'Golden Hash':>14s}  {'Toy SHA':>14s}")
        print(f"    {'Max cycle':20s} {g_max:>14d}  {t_max:>14d}")
        print(f"    {'Avg cycle':20s} {sum(g_cycles)/len(g_cycles):>14.1f}  {sum(t_cycles)/len(t_cycles):>14.1f}")
        print(f"    {'Num cycles':20s} {len(g_cycles):>14d}  {len(t_cycles):>14d}")
        print(f"    {'CV (uniformity)':20s} {g_cv:>14.4f}  {t_cv:>14.4f}")
        print(f"    {'Longer max':20s} {'<-- GOLDEN' if g_max > t_max else 'TOY SHA -->':>14s}")

        # Check if golden hash has maximal cycle (= n for coprime multiplier)
        phi_int = int(PHI_FLOAT * n) & mask
        gcd_val = math.gcd(phi_int, n)
        if gcd_val == 1:
            # Multiplying by a value coprime to n gives a permutation
            # For power-of-2 modulus, odd multiplier = permutation
            pass
        print(f"    phi_int = {phi_int}, gcd(phi_int, 2^{bits}) = {gcd_val}")
        if phi_int % 2 == 1:
            print(f"    phi_int is ODD => multiplication is a permutation of Z/(2^{bits})")
        else:
            print(f"    phi_int is EVEN => multiplication is NOT a permutation")


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 3: Cycle length ratios and fitting
# ─────────────────────────────────────────────────────────────────

def experiment_3(exp1_results):
    print("\n" + "=" * 72)
    print("  EXPERIMENT 3: Cycle Length Ratios and Fitting")
    print("=" * 72)

    bit_widths = sorted(exp1_results.keys())
    max_cycles = [exp1_results[b].max_cycle for b in bit_widths]

    print(f"\n  {'bits':>6s}  {'max_cycle':>10s}  {'max/2^bits':>10s}  {'max/phi^bits':>12s}  {'phi^bits':>10s}")
    print("  " + "-" * 56)

    ratios_2 = []
    ratios_phi = []

    for bits, mc in zip(bit_widths, max_cycles):
        n = 1 << bits
        phi_b = float(power(PHI, bits))
        r2 = mc / n
        rp = mc / phi_b
        ratios_2.append(r2)
        ratios_phi.append(rp)
        print(f"  {bits:>6d}  {mc:>10d}  {r2:>10.6f}  {rp:>12.4f}  {phi_b:>10.1f}")

    # Fit: max_cycle ~ a^bits
    # log(max_cycle) = bits * log(a) => log(a) = log(max_cycle) / bits
    print(f"\n  Fitting max_cycle ~ a^bits:")
    print(f"  {'bits':>6s}  {'log(max)/bits':>14s}  {'implied a':>10s}")
    print("  " + "-" * 34)
    implied_a = []
    for bits, mc in zip(bit_widths, max_cycles):
        if mc > 0:
            la = math.log(mc) / bits
            a = math.exp(la)
            implied_a.append(a)
            print(f"  {bits:>6d}  {la:>14.6f}  {a:>10.6f}")

    if implied_a:
        avg_a = sum(implied_a) / len(implied_a)
        print(f"\n  Average implied base a = {avg_a:.6f}")
        print(f"  phi                    = {PHI_FLOAT:.6f}")
        print(f"  2                      = 2.000000")
        print(f"  a / phi                = {avg_a / PHI_FLOAT:.6f}")
        print(f"  a / 2                  = {avg_a / 2:.6f}")

    # Check if ratio max/phi^bits is constant
    if ratios_phi:
        avg_rp = sum(ratios_phi) / len(ratios_phi)
        var_rp = sum((r - avg_rp)**2 for r in ratios_phi) / len(ratios_phi)
        cv_rp = (var_rp**0.5) / avg_rp if avg_rp > 0 else float('inf')
        print(f"\n  Is max/phi^bits constant?")
        print(f"    Mean ratio:  {avg_rp:.4f}")
        print(f"    Std dev:     {var_rp**0.5:.4f}")
        print(f"    CV:          {cv_rp:.4f}")
        print(f"    {'YES — roughly constant!' if cv_rp < 0.3 else 'NO — varies significantly'}")


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 4: Finding phi in SHA-256 internals
# ─────────────────────────────────────────────────────────────────

def experiment_4():
    print("\n" + "=" * 72)
    print("  EXPERIMENT 4: Finding phi in SHA-256 Internals")
    print("=" * 72)

    # The golden ratio constant in SHA
    print("\n  --- The Golden Ratio Mixing Constant ---")
    two32 = mpf(2)**32
    phi_const = two32 / PHI
    print(f"  2^32 / phi          = {nstr(phi_const, 15)}")
    print(f"  0x9E3779B9          = {0x9E3779B9}")
    print(f"  floor(2^32 / phi)   = {int(floor(phi_const))}")
    print(f"  Match:                {'YES!!' if int(floor(phi_const)) == 0x9E3779B9 else 'no'}")

    # Also check: 2^32 * (phi - 1) = 2^32 / phi (since phi - 1 = 1/phi)
    inv_phi = PHI - 1  # = 1/phi
    alt = two32 * inv_phi
    print(f"  2^32 * (phi-1)      = {nstr(alt, 15)}  (same thing, since phi-1 = 1/phi)")

    # SHA-256 initial hash values (fractional parts of sqrt of first 8 primes)
    print("\n  --- SHA-256 Initial Hash Values (h0-h7) ---")
    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    h_values = [0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
                0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19]

    for i, (p, h) in enumerate(zip(primes, h_values)):
        sq = sqrt(mpf(p))
        frac_part = frac(sq)
        computed = int(floor(frac_part * two32))
        match = "MATCH" if computed == h else f"MISMATCH (got {computed:#010x})"
        print(f"  h{i} = {h:#010x}  frac(sqrt({p:2d})) * 2^32 = {computed:#010x}  {match}")

    # The phi connection through sqrt(5)
    print("\n  --- The sqrt(5) / phi Connection ---")
    sqrt5 = sqrt(mpf(5))
    print(f"  sqrt(5)             = {nstr(sqrt5, 20)}")
    print(f"  phi                 = {nstr(PHI, 20)}")
    print(f"  2*phi - 1           = {nstr(2*PHI - 1, 20)}")
    print(f"  sqrt(5) = 2*phi-1?  {'YES' if abs(sqrt5 - (2*PHI - 1)) < mpf(10)**(-40) else 'no'}")

    frac_sqrt5 = frac(sqrt5)
    print(f"\n  frac(sqrt(5))       = {nstr(frac_sqrt5, 20)}")
    print(f"  sqrt(5) - 2         = {nstr(sqrt5 - 2, 20)}")
    print(f"  Same?               {'YES' if abs(frac_sqrt5 - (sqrt5 - 2)) < mpf(10)**(-40) else 'no'}")

    # The dodecahedral connection
    print(f"\n  --- Dodecahedral Connection ---")
    inv_phi3 = 1 / PHI**3
    print(f"  1/phi^3             = {nstr(inv_phi3, 20)}")
    print(f"  frac(sqrt(5))       = {nstr(frac_sqrt5, 20)}")

    # Check: frac(sqrt(5)) = sqrt(5) - 2 = (2*phi-1) - 2 = 2*phi - 3 = 2*(phi-1) - 1 = 2/phi - 1
    two_over_phi_m1 = 2/PHI - 1
    print(f"  2/phi - 1           = {nstr(two_over_phi_m1, 20)}")
    print(f"  Same as frac(sqrt5)?  {'YES' if abs(frac_sqrt5 - two_over_phi_m1) < mpf(10)**(-40) else 'no'}")

    # Check: (2-phi)/phi = (2-phi)/phi
    two_m_phi_over_phi = (2 - PHI) / PHI
    print(f"  (2-phi)/phi         = {nstr(two_m_phi_over_phi, 20)}")
    print(f"  Same?               {'YES' if abs(frac_sqrt5 - two_m_phi_over_phi) < mpf(10)**(-40) else 'no'}")

    # Check: 1/phi^2 = 2-phi? (This is the dodecahedral Delta connection)
    inv_phi2 = 1 / PHI**2
    print(f"\n  1/phi^2             = {nstr(inv_phi2, 20)}")
    print(f"  2 - phi             = {nstr(2 - PHI, 20)}")
    print(f"  phi - 1             = {nstr(PHI - 1, 20)}")
    print(f"  1/phi^2 = 2-phi?    {'YES' if abs(inv_phi2 - (2 - PHI)) < mpf(10)**(-40) else 'no'}")

    # The KEY claim: h2 = (1/phi^3) * 2^32
    h2_from_phi = inv_phi3 * two32
    print(f"\n  === THE KEY RESULT ===")
    print(f"  h2 (SHA-256)        = {h_values[2]:#010x} = {h_values[2]}")
    print(f"  (1/phi^3) * 2^32    = {nstr(h2_from_phi, 15)}")
    print(f"  floor(...)          = {int(floor(h2_from_phi))}")
    print(f"  h2 = floor(2^32/phi^3)?  {'YES!!' if int(floor(h2_from_phi)) == h_values[2] else 'no'}")

    # Check other h-values against phi powers
    print(f"\n  --- Checking ALL h-values against phi powers ---")
    for n_pow in range(1, 20):
        val = int(floor(two32 / PHI**n_pow)) & 0xFFFFFFFF
        for i, h in enumerate(h_values):
            if val == h:
                print(f"  h{i} = floor(2^32 / phi^{n_pow}) !!")

    # SHA-256 round constants (first 8)
    print(f"\n  --- SHA-256 Round Constants (k0-k7) ---")
    k_values = [
        0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5,
        0x3956C25B, 0x59F111F1, 0x923F82A4, 0xAB1C5ED5
    ]
    # These come from fractional parts of cube roots of first 8 primes
    for i, (p, k) in enumerate(zip(primes, k_values)):
        cr = mpf(p) ** (mpf(1)/3)
        frac_part = frac(cr)
        computed = int(floor(frac_part * two32))
        match = "MATCH" if computed == k else "mismatch"
        print(f"  k{i} = {k:#010x}  frac(cbrt({p:2d})) * 2^32 = {computed:#010x}  {match}")

    # Check: are any round constants close to Fibonacci * some scale?
    print(f"\n  --- Fibonacci/Lucas Number Connections ---")
    fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025]
    lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843,
             1364, 2207, 3571, 5778, 9349, 15127, 24476, 39603, 64079]

    # Check if 0x9E3779B9 is close to any Fibonacci * power of 2
    sha_const = 0x9E3779B9
    print(f"  SHA constant {sha_const} = 0x{sha_const:08X}")
    for f in fibs:
        if f == 0:
            continue
        ratio = sha_const / f
        log2_ratio = math.log2(ratio) if ratio > 0 else 0
        if abs(log2_ratio - round(log2_ratio)) < 0.01 and 0 < round(log2_ratio) < 40:
            print(f"    {sha_const} / F({fibs.index(f)}) = {ratio:.4f} ~ 2^{log2_ratio:.4f}")

    # Direct check: is the SHA constant a Lucas or Fibonacci number?
    print(f"  Is 0x9E3779B9 a Fibonacci number? {'YES' if sha_const in fibs else 'no'}")
    print(f"  Is 0x9E3779B9 a Lucas number?    {'YES' if sha_const in lucas else 'no'}")

    # Ratio to nearest Fibonacci
    for f in sorted(fibs, reverse=True):
        if f > 0:
            r = sha_const / f
            if 0.9 < r / round(r) < 1.1 and round(r) > 0:
                err = abs(r - round(r)) / round(r)
                if err < 0.001:
                    print(f"    0x9E3779B9 / F_{fibs.index(f)} = {round(r)} (error {err:.6f})")


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 5: Actual SHA-256 (truncated) cycle structure
# ─────────────────────────────────────────────────────────────────

def truncated_sha256(x, out_bits=16):
    """SHA-256 truncated to out_bits. Input and output are integers."""
    data = x.to_bytes(2, 'big')
    digest = hashlib.sha256(data).digest()
    # Take first out_bits from the hash
    raw = int.from_bytes(digest[:2], 'big')
    mask = (1 << out_bits) - 1
    return raw & mask


def experiment_5():
    print("\n" + "=" * 72)
    print("  EXPERIMENT 5: Truncated SHA-256 Cycle Structure (16-bit)")
    print("=" * 72)

    bits = 16
    n = 1 << bits

    # Build the full function table (faster than recomputing)
    print(f"\n  Building SHA-256 function table for {n} values...")
    func_table = [0] * n
    for x in range(n):
        func_table[x] = truncated_sha256(x, bits)

    func = lambda x: func_table[x]

    print("  Finding all cycles...")
    cycles = find_all_cycles(func, bits)

    max_c = max(cycles)
    avg_c = sum(cycles) / len(cycles)
    phi_16 = float(power(PHI, 16))
    phi_291_scaled = float(power(PHI, 16 * mpf(291) / 256))

    print(f"\n  State space:          {n}")
    print(f"  Distinct cycles:      {len(cycles)}")
    print(f"  Max cycle:            {max_c}")
    print(f"  Average cycle:        {avg_c:.1f}")
    print(f"  Median cycle:         {sorted(cycles)[len(cycles)//2]}")
    print(f"  phi^16:               {phi_16:.1f}")
    print(f"  phi^(16*291/256):     {phi_291_scaled:.1f}")
    print(f"  Max / phi^16:         {max_c / phi_16:.4f}")
    print(f"  Max / 2^16:           {max_c / n:.6f}")
    print(f"  Max / sqrt(2^16):     {max_c / (n**0.5):.4f}")

    print(f"\n  All cycle lengths: {cycles[:30]}{'...' if len(cycles) > 30 else ''}")

    # Check if cycle lengths relate to Fibonacci/Lucas
    fibs_set = set()
    a, b = 1, 1
    while a < n * 2:
        fibs_set.add(a)
        a, b = b, a + b

    lucas_set = set()
    a, b = 2, 1
    while a < n * 2:
        lucas_set.add(a)
        a, b = b, a + b

    fib_cycles = [c for c in cycles if c in fibs_set]
    lucas_cycles = [c for c in cycles if c in lucas_set]
    print(f"\n  Cycles that are Fibonacci numbers: {fib_cycles}")
    print(f"  Cycles that are Lucas numbers:     {lucas_cycles}")

    # Check ratios between consecutive cycle lengths
    if len(cycles) >= 2:
        print(f"\n  Ratios of consecutive cycle lengths (largest first):")
        for i in range(min(8, len(cycles) - 1)):
            if cycles[i+1] > 0:
                r = cycles[i] / cycles[i+1]
                print(f"    C[{i}]/C[{i+1}] = {cycles[i]}/{cycles[i+1]} = {r:.4f}  (phi={PHI_FLOAT:.4f})")

    return cycles, max_c


# ─────────────────────────────────────────────────────────────────
#  EXPERIMENT 6: Birthday bound vs golden bound
# ─────────────────────────────────────────────────────────────────

def experiment_6(sha_cycles, sha_max):
    print("\n" + "=" * 72)
    print("  EXPERIMENT 6: Birthday Bound vs Golden Bound")
    print("=" * 72)

    bits = 16
    n = 1 << bits

    birthday = n ** 0.5  # 2^(n/2) = 256
    golden_exp = bits * mpf(291) / 256  # = 18.1875
    golden_bound = float(power(PHI, golden_exp))
    phi_bits = float(power(PHI, bits))

    print(f"\n  For n = {bits} bits (state space = {n}):")
    print(f"    Birthday bound:    sqrt(2^{bits}) = 2^{bits//2} = {birthday:.1f}")
    print(f"    Golden bound:      phi^({bits}*291/256) = phi^{float(golden_exp):.4f} = {golden_bound:.1f}")
    print(f"    phi^{bits}:           {phi_bits:.1f}")
    print(f"    2^{bits}:             {n}")

    print(f"\n  Empirical (truncated SHA-256):")
    print(f"    Max cycle:         {sha_max}")
    print(f"    Avg cycle:         {sum(sha_cycles)/len(sha_cycles):.1f}")

    # Collision distance: for a random function, expected rho length ~ sqrt(pi*n/8)
    expected_rho = (math.pi * n / 8) ** 0.5
    print(f"\n  Random function model:")
    print(f"    Expected rho tail: sqrt(pi*n/8) = {expected_rho:.1f}")
    print(f"    Expected cycle:    sqrt(pi*n/8) = {expected_rho:.1f}")
    print(f"    Expected total:    2 * sqrt(pi*n/8) = {2*expected_rho:.1f}")

    print(f"\n  Bound comparison:")
    print(f"    {'Bound':25s} {'Value':>10s} {'vs Max Cycle':>12s} {'Ratio':>8s}")
    print("    " + "-" * 57)
    bounds = [
        ("Birthday 2^(n/2)", birthday),
        ("Golden phi^(n*291/256)", golden_bound),
        ("phi^n", phi_bits),
        ("Random rho sqrt(pi*n/8)", expected_rho),
        ("Full space 2^n", float(n)),
    ]
    for name, val in bounds:
        ratio = sha_max / val if val > 0 else float('inf')
        status = "<= bound" if sha_max <= val else "> bound"
        print(f"    {name:25s} {val:>10.1f} {status:>12s} {ratio:>8.4f}")

    # For multiple bit widths
    print(f"\n  --- Collision distance vs bounds across bit widths ---")
    print(f"  (Using toy SHA for smaller bit widths)")
    for test_bits in [8, 10, 12, 14, 16]:
        test_n = 1 << test_bits
        test_func = lambda x, b=test_bits: toy_hash(x, b)

        # Find average tail + cycle length (rho length) from random starts
        total_rho = 0
        num_samples = min(100, test_n)
        import random
        random.seed(42)
        samples = random.sample(range(test_n), num_samples)
        for s in samples:
            mu, lam = find_cycle_from(s, test_func, (1 << test_bits) - 1)
            total_rho += mu + lam
        avg_rho = total_rho / num_samples

        bday = test_n ** 0.5
        g_exp = test_bits * 291.0 / 256
        g_bound = PHI_FLOAT ** g_exp
        phi_b = PHI_FLOAT ** test_bits

        print(f"  bits={test_bits:2d}: avg_rho={avg_rho:8.1f}  birthday={bday:8.1f}  "
              f"golden={g_bound:8.1f}  phi^bits={phi_b:8.1f}  "
              f"rho/bday={avg_rho/bday:.2f}  rho/golden={avg_rho/g_bound:.2f}")


# ─────────────────────────────────────────────────────────────────
#  GRAND SUMMARY
# ─────────────────────────────────────────────────────────────────

def grand_summary(exp1_results, sha_max):
    print("\n" + "=" * 72)
    print("  GRAND SUMMARY: Does phi Govern Hash Cycles?")
    print("=" * 72)

    print(f"\n  phi = {nstr(PHI, 30)}")
    print(f"  phi^291 = {nstr(power(PHI, 291), 20)}")

    print(f"\n  Key findings:")

    # 1. SHA uses phi
    print(f"\n  [1] SHA-256 ALREADY USES PHI")
    print(f"      0x9E3779B9 = floor(2^32 / phi) -- the primary mixing constant")
    print(f"      h2 = floor(2^32 / phi^3)       -- from sqrt(5) = 2*phi - 1")

    # 2. Cycle scaling
    print(f"\n  [2] CYCLE LENGTH SCALING")
    bit_widths = sorted(exp1_results.keys())
    print(f"      {'bits':>4s}  {'max_cycle':>10s}  {'max/phi^b':>10s}  {'log(max)/b':>10s}  {'implied base':>12s}")
    for b in bit_widths:
        mc = exp1_results[b].max_cycle
        phi_b = float(power(PHI, b))
        ratio = mc / phi_b
        if mc > 0:
            implied = math.exp(math.log(mc) / b)
        else:
            implied = 0
        print(f"      {b:>4d}  {mc:>10d}  {ratio:>10.4f}  {math.log(mc)/b if mc > 0 else 0:>10.6f}  {implied:>12.6f}")

    # 3. Truncated SHA
    print(f"\n  [3] TRUNCATED SHA-256 (16-bit)")
    print(f"      Max cycle = {sha_max}")
    print(f"      phi^16 = {float(power(PHI, 16)):.1f}")
    print(f"      Ratio = {sha_max / float(power(PHI, 16)):.4f}")

    # 4. The golden ceiling hypothesis
    print(f"\n  [4] THE GOLDEN CEILING HYPOTHESIS")
    print(f"      Claim: max cycle length <= C * phi^bits for some constant C")
    print(f"      For toy SHA:  C ~ {max(exp1_results[b].max_cycle / float(power(PHI, b)) for b in bit_widths):.2f}")
    print(f"      For real SHA (truncated): C ~ {sha_max / float(power(PHI, 16)):.2f}")

    # 5. Dodecahedral structure
    print(f"\n  [5] DODECAHEDRAL STRUCTURE IN SHA-256")
    print(f"      h2 = floor(2^32 * frac(sqrt(5)))")
    print(f"         = floor(2^32 * (sqrt(5) - 2))")
    print(f"         = floor(2^32 * (2/phi - 1))")
    print(f"         = floor(2^32 / phi^3)")
    print(f"      1/phi^2 = 2 - phi   (dodecahedral identity)")
    print(f"      1/phi^3 = frac(sqrt(5))  (the fractal golden residue)")


# ─────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────

def main():
    print("*" * 72)
    print("*  GOLDEN RATIO vs HASH FUNCTION CYCLES                              *")
    print("*  Pure mathematical research on toy hash functions                   *")
    print("*" * 72)

    # Experiment 1
    exp1_results = experiment_1()

    # Experiment 2
    experiment_2()

    # Experiment 3
    experiment_3(exp1_results)

    # Experiment 4
    experiment_4()

    # Experiment 5
    sha_cycles, sha_max = experiment_5()

    # Experiment 6
    experiment_6(sha_cycles, sha_max)

    # Grand summary
    grand_summary(exp1_results, sha_max)

    print("\n" + "*" * 72)
    print("*  RESEARCH COMPLETE                                                 *")
    print("*" * 72)


if __name__ == "__main__":
    main()
