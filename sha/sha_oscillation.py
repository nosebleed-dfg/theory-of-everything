"""
SHA_OSCILLATION — tests Ch/Maj amplification/dampening oscillation of golden signal through SHA-256 rounds
nos3bl33d

Seven tests: Ch amplifier, Maj dampener, oscillation, phase detection, carry growth, resonance search, MI.
"""

import hashlib
import struct
import math
import numpy as np
from typing import List, Tuple
import sys

# ─── Constants ───────────────────────────────────────────────────────────────

PHI = (1 + math.sqrt(5)) / 2          # golden ratio ~1.618
GOLDEN_ANGLE = 137.507764              # degrees
MOD32 = 2**32
MASK32 = MOD32 - 1
N_TRIALS = 10000
N_ROUNDS_OSCILLATION = 20
N_SHA_ROUNDS = 64

# SHA-256 initial hash values
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# SHA-256 round constants
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

np.random.seed(42)


# ─── Primitives ──────────────────────────────────────────────────────────────

def ch(e: int, f: int, g: int) -> int:
    """SHA-256 Ch (choice/select): Ch(e,f,g) = (e AND f) XOR (NOT e AND g)"""
    return ((e & f) ^ (~e & g)) & MASK32


def maj(a: int, b: int, c: int) -> int:
    """SHA-256 Maj (majority): Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)"""
    return ((a & b) ^ (a & c) ^ (b & c)) & MASK32


def rotr(x: int, n: int) -> int:
    """32-bit right rotate."""
    return ((x >> n) | (x << (32 - n))) & MASK32


def sigma0(x: int) -> int:
    """SHA-256 Sigma0: ROTR(2) XOR ROTR(13) XOR ROTR(22)"""
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)


def sigma1(x: int) -> int:
    """SHA-256 Sigma1: ROTR(6) XOR ROTR(11) XOR ROTR(25)"""
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)


def lsigma0(x: int) -> int:
    """SHA-256 lowercase sigma0 for message schedule."""
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)


def lsigma1(x: int) -> int:
    """SHA-256 lowercase sigma1 for message schedule."""
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)


def add32(*args: int) -> int:
    """Modular addition mod 2^32."""
    return sum(args) & MASK32


# ─── Golden Generators ──────────────────────────────────────────────────────

def fibonacci_mod32(n: int) -> List[int]:
    """Generate n Fibonacci numbers mod 2^32."""
    fibs = [1, 1]
    for i in range(2, n):
        fibs.append((fibs[-1] + fibs[-2]) & MASK32)
    return fibs[:n]


def golden_sequence(n: int) -> List[int]:
    """Generate n values from golden angle spacing: int(i * phi * 2^32) mod 2^32."""
    return [int(i * PHI * MOD32) & MASK32 for i in range(n)]


def golden_template(n: int) -> np.ndarray:
    """Binary golden template: bit patterns from golden angle spacing."""
    seq = golden_sequence(n)
    # Extract bit 15 (middle bit) as a simple binary template
    return np.array([(x >> 15) & 1 for x in seq], dtype=np.float64)


def random_u32(n: int) -> List[int]:
    """Generate n random 32-bit unsigned ints."""
    return [int(x) for x in np.random.randint(0, MOD32, size=n, dtype=np.uint64)]


# ─── Correlation Metrics ─────────────────────────────────────────────────────

def bit_correlation(values: List[int], template: np.ndarray) -> float:
    """
    Measure correlation between a sequence of 32-bit values and a golden template.
    Extracts bit 15 from each value and computes Pearson correlation with template.
    """
    n = min(len(values), len(template))
    bits = np.array([(v >> 15) & 1 for v in values[:n]], dtype=np.float64)
    t = template[:n]

    # Handle constant arrays
    if np.std(bits) < 1e-12 or np.std(t) < 1e-12:
        return 0.0

    return float(np.corrcoef(bits, t)[0, 1])


def multi_bit_correlation(values: List[int], template_seq: List[int]) -> float:
    """
    Correlation across all 32 bit positions, averaged.
    More robust than single-bit measurement.
    """
    n = min(len(values), len(template_seq))
    total_corr = 0.0
    valid_bits = 0

    for bit in range(32):
        val_bits = np.array([(v >> bit) & 1 for v in values[:n]], dtype=np.float64)
        tpl_bits = np.array([(t >> bit) & 1 for t in template_seq[:n]], dtype=np.float64)

        if np.std(val_bits) < 1e-12 or np.std(tpl_bits) < 1e-12:
            continue

        corr = np.corrcoef(val_bits, tpl_bits)[0, 1]
        if not np.isnan(corr):
            total_corr += abs(corr)
            valid_bits += 1

    return total_corr / valid_bits if valid_bits > 0 else 0.0


def autocorrelation_at_offset(values: List[int], offset: int) -> float:
    """Autocorrelation of the sequence at a given offset, using all 32 bits."""
    n = len(values) - offset
    if n <= 0:
        return 0.0

    # Flatten to bit sequences across all 32 bits
    x = np.array(values[:n], dtype=np.float64)
    y = np.array(values[offset:offset + n], dtype=np.float64)

    if np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return 0.0

    return float(np.corrcoef(x, y)[0, 1])


def popcount_correlation(values: List[int], template_seq: List[int]) -> float:
    """Correlation using popcount (number of set bits) as the signal."""
    n = min(len(values), len(template_seq))
    v_pop = np.array([bin(v).count('1') for v in values[:n]], dtype=np.float64)
    t_pop = np.array([bin(t).count('1') for t in template_seq[:n]], dtype=np.float64)

    if np.std(v_pop) < 1e-12 or np.std(t_pop) < 1e-12:
        return 0.0

    return float(np.corrcoef(v_pop, t_pop)[0, 1])


# ─── SHA-256 Internal Round (partial, for phase detection) ───────────────────

def sha256_message_schedule(block: bytes) -> List[int]:
    """Expand a 512-bit (64-byte) block into 64 message schedule words."""
    assert len(block) == 64
    W = list(struct.unpack('>16I', block))
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    return W


def sha256_round_states(block: bytes) -> List[Tuple[int, ...]]:
    """
    Run SHA-256 compression on one block, returning the 8-word state after EACH round.
    Returns list of 64 states, each a tuple of (a, b, c, d, e, f, g, h).
    """
    W = sha256_message_schedule(block)

    a, b, c, d, e, f, g, h = H0

    states = []
    for i in range(64):
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))

        h = g
        g = f
        f = e
        e = add32(d, T1)
        d = c
        c = b
        b = a
        a = add32(T1, T2)

        states.append((a, b, c, d, e, f, g, h))

    return states


# ─── Carry Chain Analysis ────────────────────────────────────────────────────

def measure_carry_chain(a: int, b: int) -> int:
    """
    Measure the longest carry chain when adding a + b in 32-bit arithmetic.
    A carry chain is a consecutive sequence of bit positions where a carry propagates.
    """
    carry = 0
    max_chain = 0
    current_chain = 0

    for bit in range(32):
        bit_a = (a >> bit) & 1
        bit_b = (b >> bit) & 1

        total = bit_a + bit_b + carry
        new_carry = total >> 1

        if new_carry:
            current_chain += 1
            max_chain = max(max_chain, current_chain)
        else:
            current_chain = 0

        carry = new_carry

    return max_chain


# ─── Mutual Information ──────────────────────────────────────────────────────

def mutual_information_binary(x: np.ndarray, y: np.ndarray) -> float:
    """
    Compute mutual information between two binary (0/1) sequences.
    I(X;Y) = sum p(x,y) * log2(p(x,y) / (p(x)*p(y)))
    """
    n = len(x)
    if n == 0:
        return 0.0

    # Joint distribution
    p00 = np.sum((x == 0) & (y == 0)) / n
    p01 = np.sum((x == 0) & (y == 1)) / n
    p10 = np.sum((x == 1) & (y == 0)) / n
    p11 = np.sum((x == 1) & (y == 1)) / n

    # Marginals
    px0 = p00 + p01
    px1 = p10 + p11
    py0 = p00 + p10
    py1 = p01 + p11

    mi = 0.0
    for pxy, px, py in [(p00, px0, py0), (p01, px0, py1),
                         (p10, px1, py0), (p11, px1, py1)]:
        if pxy > 0 and px > 0 and py > 0:
            mi += pxy * math.log2(pxy / (px * py))

    return mi


# ─── TEST 1: Ch as Amplifier ─────────────────────────────────────────────────

def test1_ch_amplifier():
    print("=" * 70)
    print("TEST 1: Ch (Choice) as AMPLIFIER")
    print("=" * 70)

    fibs = fibonacci_mod32(N_TRIALS)
    golden = golden_sequence(N_TRIALS)
    rand_f = random_u32(N_TRIALS)
    rand_g = random_u32(N_TRIALS)
    rand_e = random_u32(N_TRIALS)

    # Ch with golden e
    ch_golden = [ch(fibs[i], rand_f[i], rand_g[i]) for i in range(N_TRIALS)]

    # Ch with random e
    ch_random = [ch(rand_e[i], rand_f[i], rand_g[i]) for i in range(N_TRIALS)]

    # Measure correlations against golden template
    corr_golden_multi = multi_bit_correlation(ch_golden, fibs)
    corr_random_multi = multi_bit_correlation(ch_random, rand_e)

    corr_golden_pop = popcount_correlation(ch_golden, fibs)
    corr_random_pop = popcount_correlation(ch_random, rand_e)

    # Also measure how much golden structure survives through Ch
    corr_input_golden = multi_bit_correlation(fibs, golden)
    corr_ch_output_golden = multi_bit_correlation(ch_golden, golden)
    corr_ch_random_golden = multi_bit_correlation(ch_random, golden)

    print(f"\n  Multi-bit correlation (Ch_golden vs golden_e):  {corr_golden_multi:.6f}")
    print(f"  Multi-bit correlation (Ch_random vs random_e):  {corr_random_multi:.6f}")
    print(f"  Popcount correlation (Ch_golden vs golden_e):   {corr_golden_pop:.6f}")
    print(f"  Popcount correlation (Ch_random vs random_e):   {corr_random_pop:.6f}")
    print()
    print(f"  Golden structure in input (fibs):               {corr_input_golden:.6f}")
    print(f"  Golden structure in Ch(golden_e, rand, rand):   {corr_ch_output_golden:.6f}")
    print(f"  Golden structure in Ch(random_e, rand, rand):   {corr_ch_random_golden:.6f}")

    if corr_random_multi > 1e-10:
        amp_factor = corr_golden_multi / corr_random_multi
        print(f"\n  >> AMPLIFICATION FACTOR (multi-bit): {amp_factor:.4f}")
    else:
        print(f"\n  >> Random baseline too low for ratio; golden abs: {corr_golden_multi:.6f}")

    if corr_ch_random_golden > 1e-10:
        struct_ratio = corr_ch_output_golden / corr_ch_random_golden
        print(f"  >> STRUCTURAL RATIO (golden survival): {struct_ratio:.4f}")
        if struct_ratio > 1:
            print("  >> VERDICT: Ch AMPLIFIES golden structure!")
        elif struct_ratio < 1:
            print("  >> VERDICT: Ch DAMPENS golden structure.")
        else:
            print("  >> VERDICT: Ch is NEUTRAL to golden structure.")
    else:
        print(f"  >> Random golden baseline too low; golden abs: {corr_ch_output_golden:.6f}")

    print()
    return corr_golden_multi, corr_random_multi


# ─── TEST 2: Maj as Dampener ─────────────────────────────────────────────────

def test2_maj_dampener():
    print("=" * 70)
    print("TEST 2: Maj (Majority) as DAMPENER")
    print("=" * 70)

    fibs = fibonacci_mod32(N_TRIALS)
    golden = golden_sequence(N_TRIALS)
    rand_b = random_u32(N_TRIALS)
    rand_c = random_u32(N_TRIALS)
    rand_a = random_u32(N_TRIALS)

    # Maj with golden a
    maj_golden = [maj(fibs[i], rand_b[i], rand_c[i]) for i in range(N_TRIALS)]

    # Maj with random a
    maj_random = [maj(rand_a[i], rand_b[i], rand_c[i]) for i in range(N_TRIALS)]

    # Correlations
    corr_golden_multi = multi_bit_correlation(maj_golden, fibs)
    corr_random_multi = multi_bit_correlation(maj_random, rand_a)

    corr_golden_pop = popcount_correlation(maj_golden, fibs)
    corr_random_pop = popcount_correlation(maj_random, rand_a)

    # Golden structure survival
    corr_input_golden = multi_bit_correlation(fibs, golden)
    corr_maj_output_golden = multi_bit_correlation(maj_golden, golden)
    corr_maj_random_golden = multi_bit_correlation(maj_random, golden)

    print(f"\n  Multi-bit correlation (Maj_golden vs golden_a): {corr_golden_multi:.6f}")
    print(f"  Multi-bit correlation (Maj_random vs random_a): {corr_random_multi:.6f}")
    print(f"  Popcount correlation (Maj_golden vs golden_a):  {corr_golden_pop:.6f}")
    print(f"  Popcount correlation (Maj_random vs random_a):  {corr_random_pop:.6f}")
    print()
    print(f"  Golden structure in input (fibs):               {corr_input_golden:.6f}")
    print(f"  Golden structure in Maj(golden_a, rand, rand):  {corr_maj_output_golden:.6f}")
    print(f"  Golden structure in Maj(random_a, rand, rand):  {corr_maj_random_golden:.6f}")

    if corr_random_multi > 1e-10:
        damp_factor = corr_golden_multi / corr_random_multi
        print(f"\n  >> DAMPING FACTOR (multi-bit): {damp_factor:.4f}")
    else:
        print(f"\n  >> Random baseline too low for ratio; golden abs: {corr_golden_multi:.6f}")

    if corr_maj_random_golden > 1e-10:
        struct_ratio = corr_maj_output_golden / corr_maj_random_golden
        print(f"  >> STRUCTURAL RATIO (golden survival): {struct_ratio:.4f}")
        if struct_ratio > 1:
            print("  >> VERDICT: Maj AMPLIFIES golden structure (unexpected!)")
        elif struct_ratio < 1:
            print("  >> VERDICT: Maj DAMPENS golden structure (as predicted).")
        else:
            print("  >> VERDICT: Maj is NEUTRAL to golden structure.")
    else:
        print(f"  >> Random golden baseline too low; golden abs: {corr_maj_output_golden:.6f}")

    # Direct comparison: how much signal survives Ch vs Maj
    print(f"\n  [Comparison] Ch golden survival: see Test 1")
    print(f"  [Comparison] Maj golden survival: {corr_maj_output_golden:.6f}")
    print()
    return corr_golden_multi, corr_random_multi


# ─── TEST 3: The Oscillation ─────────────────────────────────────────────────

def test3_oscillation():
    print("=" * 70)
    print("TEST 3: Ch/Maj OSCILLATION over 20 rounds")
    print("=" * 70)

    n = 2000  # Smaller for multi-round
    fibs = fibonacci_mod32(n)
    golden = golden_sequence(n)

    # Start with golden signal
    signal = list(fibs)
    measurements = []

    for rnd in range(N_ROUNDS_OSCILLATION):
        rand1 = random_u32(n)
        rand2 = random_u32(n)

        if rnd % 2 == 0:
            # Ch round (expected: amplify)
            signal = [ch(signal[i], rand1[i], rand2[i]) for i in range(n)]
            op = "Ch"
        else:
            # Maj round (expected: dampen)
            signal = [maj(signal[i], rand1[i], rand2[i]) for i in range(n)]
            op = "Maj"

        # Measure golden correlation after this round
        corr = multi_bit_correlation(signal, golden)
        pop_corr = popcount_correlation(signal, golden)
        measurements.append((rnd + 1, op, corr, pop_corr))

        print(f"  Round {rnd+1:2d} ({op:3s}): multi-bit={corr:.6f}  popcount={pop_corr:+.6f}")

    # Analyze oscillation pattern
    ch_corrs = [m[2] for m in measurements if m[1] == "Ch"]
    maj_corrs = [m[2] for m in measurements if m[1] == "Maj"]

    avg_ch = np.mean(ch_corrs)
    avg_maj = np.mean(maj_corrs)

    print(f"\n  Average correlation after Ch rounds:  {avg_ch:.6f}")
    print(f"  Average correlation after Maj rounds: {avg_maj:.6f}")

    if avg_ch > avg_maj:
        print("  >> PATTERN: Ch rounds show HIGHER correlation (amplification)")
    elif avg_ch < avg_maj:
        print("  >> PATTERN: Maj rounds show HIGHER correlation (unexpected)")
    else:
        print("  >> PATTERN: No difference between Ch and Maj rounds")

    # Check for oscillation: do correlations alternate?
    corr_seq = [m[2] for m in measurements]
    diffs = [corr_seq[i+1] - corr_seq[i] for i in range(len(corr_seq)-1)]
    sign_changes = sum(1 for i in range(len(diffs)-1) if diffs[i] * diffs[i+1] < 0)
    max_possible = len(diffs) - 1

    print(f"  Sign changes in diff sequence: {sign_changes}/{max_possible}")
    if sign_changes > max_possible * 0.6:
        print("  >> OSCILLATION DETECTED! Signal alternates direction.")
    elif sign_changes > max_possible * 0.3:
        print("  >> WEAK OSCILLATION. Some alternation present.")
    else:
        print("  >> NO OSCILLATION. Monotonic decay or flat.")

    # Decay analysis
    if len(corr_seq) >= 4:
        first_half = np.mean(corr_seq[:len(corr_seq)//2])
        second_half = np.mean(corr_seq[len(corr_seq)//2:])
        print(f"  First-half avg: {first_half:.6f}, Second-half avg: {second_half:.6f}")
        if first_half > second_half * 1.1:
            print("  >> DECAYING envelope detected.")
        else:
            print("  >> No clear decay in envelope.")

    print()
    return measurements


# ─── TEST 4: Phase Detection across 64 SHA-256 rounds ────────────────────────

def test4_phase_detection():
    print("=" * 70)
    print("TEST 4: Phase Detection — signal at EVERY SHA-256 round (1-64)")
    print("=" * 70)

    n_inputs = 500
    fibs = fibonacci_mod32(16)  # 16 words = 512 bits = one block

    # Build a golden-structured block
    golden_block = b''
    for f in fibs:
        golden_block += struct.pack('>I', f & MASK32)
    assert len(golden_block) == 64

    # Build a random block for comparison
    rand_block = bytes(np.random.randint(0, 256, size=64, dtype=np.uint8))

    # Get round states for golden input
    golden_states = sha256_round_states(golden_block)
    random_states = sha256_round_states(rand_block)

    # For each round: measure how structured the state is
    # Use popcount deviation from expected 16 (for 32-bit words) as a proxy
    golden_template_vals = golden_sequence(8)  # 8-word template

    print(f"\n  {'Round':>5s}  {'Golden_e':>10s}  {'Random_e':>10s}  {'Diff':>10s}  {'Signal':>10s}")
    print("  " + "-" * 50)

    round_signals = []
    peak_rounds = []

    for rnd in range(64):
        gs = golden_states[rnd]
        rs = random_states[rnd]

        # Measure multi-bit structure in the 'a' and 'e' registers
        # (a gets Maj output, e gets Ch output)
        golden_a, golden_e = gs[0], gs[4]
        random_a, random_e = rs[0], rs[4]

        # Bit pattern distance from golden template
        golden_dist_e = bin(golden_e ^ golden_template_vals[4 % len(golden_template_vals)]).count('1')
        random_dist_e = bin(random_e ^ golden_template_vals[4 % len(golden_template_vals)]).count('1')

        # Signal = how much closer the golden state is to golden template
        # Negative = golden state is closer, positive = random state is closer
        signal = random_dist_e - golden_dist_e
        round_signals.append(signal)

        if abs(signal) >= 3:
            peak_rounds.append(rnd + 1)

        marker = " **" if abs(signal) >= 3 else ""
        print(f"  {rnd+1:5d}  {golden_dist_e:10d}  {random_dist_e:10d}  {signal:+10d}  {marker}")

    # Also measure using popcount entropy across all 8 state words
    print(f"\n  Full-state analysis (all 8 registers):")
    print(f"  {'Round':>5s}  {'GoldenPC':>10s}  {'RandomPC':>10s}  {'Diff':>10s}")
    print("  " + "-" * 40)

    full_signals = []
    for rnd in range(64):
        gs = golden_states[rnd]
        rs = random_states[rnd]

        golden_pc = sum(bin(w).count('1') for w in gs)
        random_pc = sum(bin(w).count('1') for w in rs)
        diff = golden_pc - random_pc
        full_signals.append(diff)

        marker = " **" if abs(diff) >= 8 else ""
        print(f"  {rnd+1:5d}  {golden_pc:10d}  {random_pc:10d}  {diff:+10d}{marker}")

    # Statistical analysis of round signals
    signals_arr = np.array(round_signals)
    full_arr = np.array(full_signals)

    print(f"\n  Single-register signal stats:")
    print(f"    Mean: {np.mean(signals_arr):+.3f}, Std: {np.std(signals_arr):.3f}")
    print(f"    Max: {np.max(signals_arr):+d} at round {np.argmax(signals_arr)+1}")
    print(f"    Min: {np.min(signals_arr):+d} at round {np.argmin(signals_arr)+1}")
    print(f"    Significant rounds (|signal| >= 3): {peak_rounds}")

    print(f"\n  Full-state signal stats:")
    print(f"    Mean: {np.mean(full_arr):+.3f}, Std: {np.std(full_arr):.3f}")
    print(f"    Max: {np.max(full_arr):+d} at round {np.argmax(full_arr)+1}")
    print(f"    Min: {np.min(full_arr):+d} at round {np.argmin(full_arr)+1}")

    # Check for periodicity in the signal
    if len(round_signals) > 4:
        fft = np.fft.rfft(signals_arr)
        magnitudes = np.abs(fft)
        # Skip DC component
        if len(magnitudes) > 1:
            peak_freq_idx = np.argmax(magnitudes[1:]) + 1
            peak_period = 64.0 / peak_freq_idx if peak_freq_idx > 0 else float('inf')
            print(f"\n  FFT dominant period: {peak_period:.1f} rounds (freq idx {peak_freq_idx})")
            print(f"  FFT magnitude at dominant: {magnitudes[peak_freq_idx]:.3f}")
            print(f"  FFT DC magnitude: {magnitudes[0]:.3f}")

    print()
    return round_signals, full_signals


# ─── TEST 5: Carry Propagation as Growth ──────────────────────────────────────

def test5_carry_propagation():
    print("=" * 70)
    print("TEST 5: Carry Propagation — Growth Analysis")
    print("=" * 70)

    fibs = fibonacci_mod32(N_TRIALS)
    golden = golden_sequence(N_TRIALS)
    rand_a = random_u32(N_TRIALS)
    rand_b = random_u32(N_TRIALS)

    # Measure carry chains: golden + golden
    carries_gg = [measure_carry_chain(fibs[i], golden[i]) for i in range(N_TRIALS)]

    # Measure carry chains: golden + random
    carries_gr = [measure_carry_chain(fibs[i], rand_b[i]) for i in range(N_TRIALS)]

    # Measure carry chains: random + random
    carries_rr = [measure_carry_chain(rand_a[i], rand_b[i]) for i in range(N_TRIALS)]

    avg_gg = np.mean(carries_gg)
    avg_gr = np.mean(carries_gr)
    avg_rr = np.mean(carries_rr)

    std_gg = np.std(carries_gg)
    std_gr = np.std(carries_gr)
    std_rr = np.std(carries_rr)

    max_gg = np.max(carries_gg)
    max_gr = np.max(carries_gr)
    max_rr = np.max(carries_rr)

    print(f"\n  Carry chain lengths:")
    print(f"    Golden + Golden: mean={avg_gg:.3f}, std={std_gg:.3f}, max={max_gg}")
    print(f"    Golden + Random: mean={avg_gr:.3f}, std={std_gr:.3f}, max={max_gr}")
    print(f"    Random + Random: mean={avg_rr:.3f}, std={std_rr:.3f}, max={max_rr}")

    print(f"\n  Golden vs Random growth difference: {avg_gg - avg_rr:+.4f}")
    if avg_gg > avg_rr:
        print("  >> Golden inputs create LONGER carry chains (more growth)")
    elif avg_gg < avg_rr:
        print("  >> Golden inputs create SHORTER carry chains (less growth)")
    else:
        print("  >> No difference in carry chain lengths")

    # Distribution of carry chain lengths
    print(f"\n  Carry chain distribution (golden+golden):")
    hist_gg, _ = np.histogram(carries_gg, bins=range(0, 35))
    for i, count in enumerate(hist_gg):
        if count > 0:
            bar = '#' * (count * 50 // N_TRIALS)
            print(f"    len={i:2d}: {count:5d} {bar}")

    print(f"\n  Carry chain distribution (random+random):")
    hist_rr, _ = np.histogram(carries_rr, bins=range(0, 35))
    for i, count in enumerate(hist_rr):
        if count > 0:
            bar = '#' * (count * 50 // N_TRIALS)
            print(f"    len={i:2d}: {count:5d} {bar}")

    # Statistical significance
    from scipy import stats as sp_stats
    t_stat, p_val = sp_stats.ttest_ind(carries_gg, carries_rr)
    print(f"\n  t-test (golden vs random): t={t_stat:.4f}, p={p_val:.6f}")
    if p_val < 0.05:
        print("  >> STATISTICALLY SIGNIFICANT difference in carry propagation!")
    else:
        print("  >> No statistically significant difference.")

    print()
    return avg_gg, avg_gr, avg_rr


# ─── TEST 6: Resonance Search ────────────────────────────────────────────────

def test6_resonance_search():
    print("=" * 70)
    print("TEST 6: Resonance Search — correlation vs input angle (0-360)")
    print("=" * 70)

    n_per_angle = 200
    angles = range(0, 360)
    correlations = []

    print(f"\n  Scanning {len(list(angles))} angles, {n_per_angle} samples each...")

    for theta in angles:
        # Generate inputs at this angle
        inputs = [int(i * MOD32 * theta / 360) & MASK32 for i in range(n_per_angle)]

        # Hash each input through full SHA-256
        hashes = []
        for val in inputs:
            # Pack as 32-bit big-endian, pad to 64 bytes
            msg = struct.pack('>I', val)
            h = hashlib.sha256(msg).digest()
            # Extract first 32-bit word from hash
            hw = struct.unpack('>I', h[:4])[0]
            hashes.append(hw)

        # Measure autocorrelation at offset 1
        ac = autocorrelation_at_offset(hashes, 1)
        correlations.append((theta, ac))

    # Find peaks
    corr_vals = np.array([c[1] for c in correlations])
    mean_corr = np.mean(corr_vals)
    std_corr = np.std(corr_vals)

    # Significant peaks: > 2 sigma above mean
    threshold = mean_corr + 2 * std_corr
    peaks = [(theta, c) for theta, c in correlations if c > threshold]
    troughs = [(theta, c) for theta, c in correlations if c < mean_corr - 2 * std_corr]

    print(f"\n  Mean autocorrelation: {mean_corr:.6f}")
    print(f"  Std:                  {std_corr:.6f}")
    print(f"  Threshold (2-sigma):  {threshold:.6f}")

    # Check golden angle specifically
    golden_idx = int(round(GOLDEN_ANGLE))
    golden_corr = correlations[golden_idx][1]
    golden_z = (golden_corr - mean_corr) / std_corr if std_corr > 0 else 0

    print(f"\n  Golden angle ({GOLDEN_ANGLE:.1f} deg) autocorrelation: {golden_corr:.6f}")
    print(f"  Golden angle z-score: {golden_z:+.3f}")

    if abs(golden_z) > 2:
        print("  >> SIGNIFICANT! Golden angle shows anomalous correlation.")
    elif abs(golden_z) > 1:
        print("  >> MARGINAL. Golden angle slightly anomalous.")
    else:
        print("  >> Not significant at golden angle.")

    # Top 10 peaks
    sorted_corrs = sorted(correlations, key=lambda x: x[1], reverse=True)
    print(f"\n  Top 10 peaks:")
    for theta, c in sorted_corrs[:10]:
        z = (c - mean_corr) / std_corr if std_corr > 0 else 0
        golden_mark = " <-- GOLDEN!" if abs(theta - GOLDEN_ANGLE) < 1.5 else ""
        print(f"    {theta:3d} deg: {c:+.6f} (z={z:+.3f}){golden_mark}")

    # Bottom 10 troughs
    print(f"\n  Bottom 10 troughs:")
    for theta, c in sorted_corrs[-10:]:
        z = (c - mean_corr) / std_corr if std_corr > 0 else 0
        print(f"    {theta:3d} deg: {c:+.6f} (z={z:+.3f})")

    # Check if ANY angle shows significant resonance
    max_z = max(abs((c - mean_corr) / std_corr) for _, c in correlations) if std_corr > 0 else 0
    print(f"\n  Maximum |z-score| across all angles: {max_z:.3f}")
    if max_z > 3:
        print("  >> STRONG RESONANCE detected at some angle!")
    elif max_z > 2:
        print("  >> MODERATE RESONANCE detected.")
    else:
        print("  >> No significant resonance at any angle (total damping).")

    # ASCII plot of correlation vs angle
    print(f"\n  Autocorrelation vs Input Angle (ASCII plot):")
    print(f"  {'angle':>5s} | correlation")
    print(f"  " + "-" * 60)

    # Bin to every 10 degrees for readable plot
    for theta_start in range(0, 360, 10):
        chunk = corr_vals[theta_start:theta_start+10]
        avg_c = np.mean(chunk)
        # Normalize to bar width
        bar_pos = int((avg_c - mean_corr) / (3 * std_corr + 1e-12) * 25) if std_corr > 0 else 0
        bar_pos = max(-25, min(25, bar_pos))

        line = [' '] * 51
        line[25] = '|'  # center
        if bar_pos > 0:
            for j in range(26, 26 + bar_pos):
                if j < 51:
                    line[j] = '#'
        elif bar_pos < 0:
            for j in range(25 + bar_pos, 25):
                if 0 <= j:
                    line[j] = '#'

        angle_mark = " *" if abs(theta_start - GOLDEN_ANGLE) < 10 else ""
        print(f"  {theta_start:3d}-{theta_start+9:3d} {''.join(line)}{angle_mark}")

    print()
    return correlations


# ─── TEST 7: Bitwise Channel Analysis ────────────────────────────────────────

def test7_bitwise_channels():
    print("=" * 70)
    print("TEST 7: Bitwise Channel Analysis — Mutual Information per bit")
    print("=" * 70)

    n_samples = 2000
    fibs = fibonacci_mod32(n_samples)
    rand_inputs = random_u32(n_samples)

    # Hash golden inputs
    golden_hashes = []
    for val in fibs:
        msg = struct.pack('>I', val)
        h = hashlib.sha256(msg).digest()
        golden_hashes.append(h)

    # Hash random inputs
    random_hashes = []
    for val in rand_inputs:
        msg = struct.pack('>I', val)
        h = hashlib.sha256(msg).digest()
        random_hashes.append(h)

    # Build input bit arrays (use bit 15 of input as the reference signal)
    golden_input_bits = np.array([(f >> 15) & 1 for f in fibs], dtype=np.int32)
    random_input_bits = np.array([(r >> 15) & 1 for r in rand_inputs], dtype=np.int32)

    # Compute mutual information for each of the 256 output bits
    golden_mi = []
    random_mi = []

    for bit_pos in range(256):
        byte_idx = bit_pos // 8
        bit_idx = 7 - (bit_pos % 8)  # MSB first

        golden_output_bits = np.array([(h[byte_idx] >> bit_idx) & 1 for h in golden_hashes], dtype=np.int32)
        random_output_bits = np.array([(h[byte_idx] >> bit_idx) & 1 for h in random_hashes], dtype=np.int32)

        mi_g = mutual_information_binary(golden_input_bits, golden_output_bits)
        mi_r = mutual_information_binary(random_input_bits, random_output_bits)

        golden_mi.append(mi_g)
        random_mi.append(mi_r)

    golden_mi = np.array(golden_mi)
    random_mi = np.array(random_mi)

    # Statistics
    print(f"\n  Golden input MI stats (across 256 output bits):")
    print(f"    Mean: {np.mean(golden_mi):.8f}")
    print(f"    Max:  {np.max(golden_mi):.8f} at bit {np.argmax(golden_mi)}")
    print(f"    Std:  {np.std(golden_mi):.8f}")

    print(f"\n  Random input MI stats (across 256 output bits):")
    print(f"    Mean: {np.mean(random_mi):.8f}")
    print(f"    Max:  {np.max(random_mi):.8f} at bit {np.argmax(random_mi)}")
    print(f"    Std:  {np.std(random_mi):.8f}")

    # Difference: which channels leak MORE golden info?
    mi_diff = golden_mi - random_mi

    print(f"\n  MI difference (golden - random):")
    print(f"    Mean: {np.mean(mi_diff):+.8f}")
    print(f"    Max:  {np.max(mi_diff):+.8f} at bit {np.argmax(mi_diff)}")
    print(f"    Min:  {np.min(mi_diff):+.8f} at bit {np.argmin(mi_diff)}")

    # Top leaky channels
    sorted_idx = np.argsort(mi_diff)[::-1]
    print(f"\n  Top 20 'leaky' channels (golden MI >> random MI):")
    print(f"  {'Bit':>5s}  {'Byte':>5s}  {'Golden_MI':>12s}  {'Random_MI':>12s}  {'Diff':>12s}")
    print("  " + "-" * 55)
    for rank, idx in enumerate(sorted_idx[:20]):
        byte_num = idx // 8
        bit_num = idx % 8
        print(f"  {idx:5d}  {byte_num:3d}.{bit_num:1d}  {golden_mi[idx]:.8f}  {random_mi[idx]:.8f}  {mi_diff[idx]:+.8f}")

    # Are there statistically significant channels?
    mi_diff_mean = np.mean(mi_diff)
    mi_diff_std = np.std(mi_diff)
    if mi_diff_std > 0:
        z_scores = (mi_diff - mi_diff_mean) / mi_diff_std
        sig_channels = np.sum(np.abs(z_scores) > 2)
        print(f"\n  Channels with |z| > 2: {sig_channels}/256")
        print(f"  Expected by chance:    ~10 (5%)")

        if sig_channels > 20:
            print("  >> SIGNIFICANT excess of leaky channels!")
        elif sig_channels > 12:
            print("  >> MARGINAL excess of leaky channels.")
        else:
            print("  >> No significant excess (consistent with chance).")

    # Multi-bit MI: use ALL 32 input bits instead of just bit 15
    print(f"\n  Extended analysis: MI using all 32 input bits as reference...")
    golden_mi_full = np.zeros(256)
    random_mi_full = np.zeros(256)

    for input_bit in range(32):
        g_in = np.array([(f >> input_bit) & 1 for f in fibs], dtype=np.int32)
        r_in = np.array([(r >> input_bit) & 1 for r in rand_inputs], dtype=np.int32)

        for out_bit in range(256):
            byte_idx = out_bit // 8
            bit_idx = 7 - (out_bit % 8)

            g_out = np.array([(h[byte_idx] >> bit_idx) & 1 for h in golden_hashes], dtype=np.int32)
            r_out = np.array([(h[byte_idx] >> bit_idx) & 1 for h in random_hashes], dtype=np.int32)

            golden_mi_full[out_bit] += mutual_information_binary(g_in, g_out)
            random_mi_full[out_bit] += mutual_information_binary(r_in, r_out)

    # Average over input bits
    golden_mi_full /= 32
    random_mi_full /= 32
    mi_diff_full = golden_mi_full - random_mi_full

    print(f"  Full MI (avg over 32 input bits):")
    print(f"    Golden mean: {np.mean(golden_mi_full):.8f}")
    print(f"    Random mean: {np.mean(random_mi_full):.8f}")
    print(f"    Diff mean:   {np.mean(mi_diff_full):+.8f}")
    print(f"    Diff max:    {np.max(mi_diff_full):+.8f} at bit {np.argmax(mi_diff_full)}")

    # Spatial pattern: group by hash word (32-bit blocks)
    print(f"\n  MI by hash word (spatial pattern):")
    for word in range(8):
        word_mi_g = np.mean(golden_mi_full[word*32:(word+1)*32])
        word_mi_r = np.mean(random_mi_full[word*32:(word+1)*32])
        diff = word_mi_g - word_mi_r
        bar = '#' * int(abs(diff) * 100000) if abs(diff) > 0 else ''
        sign = '+' if diff > 0 else '-' if diff < 0 else ' '
        print(f"    Word {word} (H{word}): golden={word_mi_g:.8f} random={word_mi_r:.8f} diff={sign}{abs(diff):.8f} {bar}")

    print()
    return golden_mi, random_mi


# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "#" * 70)
    print("#  SHA-256 OSCILLATION PATTERN ANALYSIS")
    print("#  Hypothesis: Ch amplifies, Maj dampens golden signals")
    print("#" * 70 + "\n")

    # Run all tests
    results = {}

    results['test1'] = test1_ch_amplifier()
    results['test2'] = test2_maj_dampener()
    results['test3'] = test3_oscillation()
    results['test4'] = test4_phase_detection()
    results['test5'] = test5_carry_propagation()
    results['test6'] = test6_resonance_search()
    results['test7'] = test7_bitwise_channels()

    # ─── SYNTHESIS ────────────────────────────────────────────────────────

    print("=" * 70)
    print("SYNTHESIS: Does SHA-256 oscillate between golden amplification and damping?")
    print("=" * 70)

    t1_golden, t1_random = results['test1']
    t2_golden, t2_random = results['test2']

    print(f"\n  1. Ch amplification factor:  {t1_golden:.6f} (vs random {t1_random:.6f})")
    print(f"  2. Maj damping factor:       {t2_golden:.6f} (vs random {t2_random:.6f})")

    if t1_golden > t1_random and t2_golden < t2_random:
        print("  >> CONFIRMED: Ch amplifies, Maj dampens (as hypothesized)")
    elif t1_golden > t1_random:
        print("  >> PARTIAL: Ch amplifies but Maj does not clearly dampen")
    elif t2_golden < t2_random:
        print("  >> PARTIAL: Maj dampens but Ch does not clearly amplify")
    else:
        print("  >> Neither clear amplification nor damping detected")

    # Oscillation verdict
    measurements = results['test3']
    corr_seq = [m[2] for m in measurements]
    diffs = [corr_seq[i+1] - corr_seq[i] for i in range(len(corr_seq)-1)]
    sign_changes = sum(1 for i in range(len(diffs)-1) if diffs[i] * diffs[i+1] < 0)

    print(f"\n  3. Oscillation: {sign_changes} sign changes in 18 possible")

    # Phase verdict
    round_signals, _ = results['test4']
    max_signal = max(round_signals)
    max_round = round_signals.index(max_signal) + 1
    min_signal = min(round_signals)
    min_round = round_signals.index(min_signal) + 1

    print(f"\n  4. Phase: peak signal at round {max_round} ({max_signal:+d})")
    print(f"            trough signal at round {min_round} ({min_signal:+d})")

    # Carry verdict
    avg_gg, avg_gr, avg_rr = results['test5']
    print(f"\n  5. Carry chains: golden={avg_gg:.3f}, random={avg_rr:.3f}")

    print("\n" + "#" * 70)
    print("#  ANALYSIS COMPLETE")
    print("#" * 70 + "\n")


if __name__ == "__main__":
    main()
