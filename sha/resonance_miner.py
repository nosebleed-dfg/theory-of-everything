"""
RESONANCE_MINER — tests whether valid Bitcoin nonces come from Ch/Maj/Sigma constructive interference
nos3bl33d

Metrics: resonance score, phase alignment, golden phase distance, FFT golden frequency dominance.
"""

import struct
import hashlib
import socket
import json
import time
import math
import sys
import os
from typing import List, Tuple, Optional, Dict

# ─── Constants ──────────────────────────────────────────────────────────────

PHI = (1 + 5**0.5) / 2
PI = math.pi
GAMMA = 0.5772156649015329
MASK32 = 0xFFFFFFFF
MOD32 = 2**32
D = 3   # vertex degree (dimension)
P = 5   # face degree (pentagon)
CHI = 2 # Euler characteristic

WALLET = "bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv"
POOL_HOST = "solo.ckpool.org"
POOL_PORT = 3333

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

# ─── SHA-256 Primitives ─────────────────────────────────────────────────────

def rotr(x: int, n: int) -> int:
    return ((x >> n) | (x << (32 - n))) & MASK32

def ch(e: int, f: int, g: int) -> int:
    """Ch(e,f,g) = (e AND f) XOR (NOT e AND g)"""
    return ((e & f) ^ (~e & g)) & MASK32

def maj(a: int, b: int, c: int) -> int:
    """Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)"""
    return ((a & b) ^ (a & c) ^ (b & c)) & MASK32

def sigma0(x: int) -> int:
    """Big Sigma0: ROTR(2) XOR ROTR(13) XOR ROTR(22)"""
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sigma1(x: int) -> int:
    """Big Sigma1: ROTR(6) XOR ROTR(11) XOR ROTR(25)"""
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def lsigma0(x: int) -> int:
    """Small sigma0 for message schedule"""
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def lsigma1(x: int) -> int:
    """Small sigma1 for message schedule"""
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def add32(*args: int) -> int:
    return sum(args) & MASK32

def sha256d(d: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(d).digest()).digest()

def popcount32(x: int) -> int:
    """Count set bits in a 32-bit integer."""
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return ((x * 0x01010101) >> 24) & 0xFF

def leading_zeros(hash_bytes: bytes) -> int:
    """Count leading zero bits in hash (big-endian)."""
    h_int = int.from_bytes(hash_bytes, 'big')
    if h_int == 0:
        return 256
    return 256 - h_int.bit_length()

# ─── SHA-256 Compression with Band Tracking ─────────────────────────────────

def message_schedule(block: bytes) -> List[int]:
    """Expand 64-byte block into 64 message schedule words."""
    assert len(block) == 64, f"Block must be 64 bytes, got {len(block)}"
    W = list(struct.unpack('>16I', block))
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    return W

def sha256_compress_with_bands(block: bytes) -> Tuple[bytes, List[int], List[int], List[int]]:
    """
    Run SHA-256 compression on one 64-byte block.
    Returns (final_hash_of_block, ch_values, maj_values, sigma_values) for all 64 rounds.

    ch_values[i]    = Ch(e,f,g) at round i
    maj_values[i]   = Maj(a,b,c) at round i
    sigma_values[i] = Sigma0(a) + Sigma1(e) at round i (mod 2^32)
    """
    W = message_schedule(block)

    a, b, c, d, e, f, g, h = H0

    ch_vals = []
    maj_vals = []
    sig_vals = []

    for i in range(64):
        # Compute the three bands
        ch_i = ch(e, f, g)
        maj_i = maj(a, b, c)
        sig_i = add32(sigma0(a), sigma1(e))

        ch_vals.append(ch_i)
        maj_vals.append(maj_i)
        sig_vals.append(sig_i)

        # Standard SHA-256 round
        T1 = add32(h, sigma1(e), ch_i, K[i], W[i])
        T2 = add32(sigma0(a), maj_i)

        h = g
        g = f
        f = e
        e = add32(d, T1)
        d = c
        c = b
        b = a
        a = add32(T1, T2)

    # Feedforward: add initial hash values
    final = [
        add32(a, H0[0]), add32(b, H0[1]), add32(c, H0[2]), add32(d, H0[3]),
        add32(e, H0[4]), add32(f, H0[5]), add32(g, H0[6]), add32(h, H0[7]),
    ]

    result = struct.pack('>8I', *final)
    return result, ch_vals, maj_vals, sig_vals


def double_sha256_with_bands(header_80: bytes) -> Tuple[bytes, List[int], List[int], List[int]]:
    """
    Full double-SHA256 on an 80-byte Bitcoin header.

    Bitcoin header = 80 bytes. SHA-256 processes in 64-byte blocks.
    First hash: block1 (bytes 0-63) + block2 (bytes 64-79, padded to 64).
    Second hash: 32-byte first hash, padded to 64.

    Returns bands from the SECOND block of the FIRST hash (where the nonce lives)
    because that's the block directly influenced by the nonce.
    """
    # Pad the 80-byte header to get two 64-byte blocks for the first SHA-256
    # SHA-256 padding: append 0x80, then zeros, then 64-bit big-endian bit length
    msg = header_80 + b'\x80' + b'\x00' * (55 - (80 % 64))
    # Actually, let me do this properly.
    # 80 bytes = 640 bits. We need to pad to a multiple of 512 bits (64 bytes).
    # After 80 bytes, append 0x80 (1 byte), then zeros until 8 bytes before end,
    # then 8-byte big-endian length (640 = 0x280).
    # 80 + 1 = 81. Next multiple of 64 = 128. Zeros needed = 128 - 81 - 8 = 39.
    padded = header_80 + b'\x80' + b'\x00' * 39 + struct.pack('>Q', 640)
    assert len(padded) == 128, f"Padded length should be 128, got {len(padded)}"

    block1 = padded[:64]
    block2 = padded[64:128]

    # First SHA-256: process block1 with H0
    W1 = message_schedule(block1)
    a, b, c, d, e, f, g, h = H0
    for i in range(64):
        ch_i = ch(e, f, g)
        maj_i = maj(a, b, c)
        T1 = add32(h, sigma1(e), ch_i, K[i], W1[i])
        T2 = add32(sigma0(a), maj_i)
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)

    mid_state = [
        add32(a, H0[0]), add32(b, H0[1]), add32(c, H0[2]), add32(d, H0[3]),
        add32(e, H0[4]), add32(f, H0[5]), add32(g, H0[6]), add32(h, H0[7]),
    ]

    # Second block (contains the nonce at bytes 64-67 of the original, plus padding)
    # This is WHERE THE NONCE LIVES — the bands here are what we care about
    W2 = message_schedule(block2)
    a, b, c, d, e, f, g, h = mid_state

    ch_vals = []
    maj_vals = []
    sig_vals = []

    for i in range(64):
        ch_i = ch(e, f, g)
        maj_i = maj(a, b, c)
        sig_i = add32(sigma0(a), sigma1(e))

        ch_vals.append(ch_i)
        maj_vals.append(maj_i)
        sig_vals.append(sig_i)

        T1 = add32(h, sigma1(e), ch_i, K[i], W2[i])
        T2 = add32(sigma0(a), maj_i)
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)

    first_hash_state = [
        add32(a, mid_state[0]), add32(b, mid_state[1]),
        add32(c, mid_state[2]), add32(d, mid_state[3]),
        add32(e, mid_state[4]), add32(f, mid_state[5]),
        add32(g, mid_state[6]), add32(h, mid_state[7]),
    ]
    first_hash = struct.pack('>8I', *first_hash_state)

    # Second SHA-256: hash the 32-byte first hash
    # Pad: 32 bytes + 0x80 + 23 zeros + 8-byte length (256 bits = 0x100)
    second_block = first_hash + b'\x80' + b'\x00' * 23 + struct.pack('>Q', 256)
    assert len(second_block) == 64

    final_hash_bytes = sha256_single_block(second_block)

    return final_hash_bytes, ch_vals, maj_vals, sig_vals


def sha256_single_block(block: bytes) -> bytes:
    """SHA-256 compression of a single 64-byte block starting from H0."""
    W = message_schedule(block)
    a, b, c, d, e, f, g, h = H0
    for i in range(64):
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)
    final = [
        add32(a, H0[0]), add32(b, H0[1]), add32(c, H0[2]), add32(d, H0[3]),
        add32(e, H0[4]), add32(f, H0[5]), add32(g, H0[6]), add32(h, H0[7]),
    ]
    return struct.pack('>8I', *final)


# ─── Resonance Metrics ──────────────────────────────────────────────────────

def bitwise_agreement(a: int, b: int) -> int:
    """Count bits where a and b AGREE (both 0 or both 1). = 32 - popcount(a XOR b)."""
    return 32 - popcount32(a ^ b)

def metric_resonance_score(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    RESONANCE SCORE: How aligned are the three bands across all 64 rounds?

    For each round, count bits where all three bands agree (all 0 or all 1).
    Three-way agreement = (Ch XNOR Maj) AND (Maj XNOR Sigma) per bit.
    Sum across all 64 rounds. Normalize to [0, 1].

    High resonance = bands are in phase = constructive interference.
    """
    total_agreement = 0
    for i in range(64):
        c, m, s = ch_vals[i], maj_vals[i], sig_vals[i]
        # Three-way agreement: bits where all three are the same
        # = bits where (c XOR m) is 0 AND (m XOR s) is 0
        # = NOT(c XOR m) AND NOT(m XOR s)
        # = XNOR(c,m) AND XNOR(m,s)
        cm_agree = ~(c ^ m) & MASK32   # XNOR
        ms_agree = ~(m ^ s) & MASK32   # XNOR
        all_agree = cm_agree & ms_agree
        total_agreement += popcount32(all_agree)

    # Max possible: 64 rounds * 32 bits = 2048
    return total_agreement / 2048.0

def metric_pairwise_agreement(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    PAIRWISE AGREEMENT: Sum of bitwise agreements for all three pairs.
    Softer than three-way: counts even when only two out of three agree.
    """
    total = 0
    for i in range(64):
        c, m, s = ch_vals[i], maj_vals[i], sig_vals[i]
        total += bitwise_agreement(c, m)
        total += bitwise_agreement(m, s)
        total += bitwise_agreement(s, c)
    # Max: 64 * 3 * 32 = 6144
    return total / 6144.0

def metric_phase_alignment(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    PHASE ALIGNMENT: Treat each band value as a phase angle on the nonce circle.
    Compute cos(angle between Ch and Maj) + cos(angle between Maj and Sigma)
    + cos(angle between Sigma and Ch), averaged across rounds.

    Phase angle = value * 2*pi / 2^32 (wraps the 32-bit value onto [0, 2pi))
    cos(angle_between) = cos(phase_a - phase_b)
    """
    total = 0.0
    two_pi_over_mod = 2.0 * PI / MOD32

    for i in range(64):
        c_phase = ch_vals[i] * two_pi_over_mod
        m_phase = maj_vals[i] * two_pi_over_mod
        s_phase = sig_vals[i] * two_pi_over_mod

        total += math.cos(c_phase - m_phase)
        total += math.cos(m_phase - s_phase)
        total += math.cos(s_phase - c_phase)

    # Range: [-192, 192] (64 rounds * 3 pairs * [-1, 1])
    # Normalize to [-1, 1]
    return total / 192.0

def metric_golden_phase(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    GOLDEN PHASE: How close is (Ch + Maj + Sigma) to phi * 2^32 on the circle?
    Average across rounds. Small distance = high golden alignment.
    """
    phi_on_circle = int(PHI * MOD32) & MASK32  # phi * 2^32 mod 2^32

    total_closeness = 0.0
    for i in range(64):
        combined = add32(ch_vals[i], maj_vals[i], sig_vals[i])
        # Circular distance: min of forward and backward distance
        diff = abs(combined - phi_on_circle)
        circ_dist = min(diff, MOD32 - diff)
        # Normalize distance to [0, 1] where 0 = at phi, 1 = at antipode
        closeness = 1.0 - (circ_dist / (MOD32 / 2))
        total_closeness += closeness

    return total_closeness / 64.0

def metric_golden_product(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    GOLDEN PRODUCT: Check if the RATIO Ch/Maj is close to phi for each round.
    Uses the circular distance between Ch and phi*Maj (mod 2^32).
    """
    total_closeness = 0.0
    count = 0

    for i in range(64):
        m = maj_vals[i]
        if m == 0:
            continue
        # phi * Maj mod 2^32
        phi_maj = int(PHI * m) & MASK32
        diff = abs(ch_vals[i] - phi_maj)
        circ_dist = min(diff, MOD32 - diff)
        closeness = 1.0 - (circ_dist / (MOD32 / 2))
        total_closeness += closeness
        count += 1

    return total_closeness / count if count > 0 else 0.0

def metric_cancellation(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    CANCELLATION (the ANTI-metric): Ch XOR Maj XOR Sigma.
    High cancellation was shown to correlate with LOW zeros.
    This is the OPPOSITE of resonance. Included as control/comparison.
    """
    total_set_bits = 0
    for i in range(64):
        xor_all = ch_vals[i] ^ maj_vals[i] ^ sig_vals[i]
        total_set_bits += popcount32(xor_all)
    # Normalize: max = 64 * 32 = 2048
    return total_set_bits / 2048.0

def metric_frequency_golden(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    FREQUENCY METRIC: FFT of the 64-round Ch/Maj/Sigma sequences.
    Check if golden frequency (1/phi of Nyquist) dominates.

    The golden frequency index in a 64-point FFT = 64/phi ≈ 39.5 → indices 39 and 40.
    Also check phi^-1 * 64 ≈ 39.5 and phi^-2 * 64 ≈ 24.4.
    """
    import numpy as np

    ch_arr = np.array(ch_vals, dtype=np.float64)
    maj_arr = np.array(maj_vals, dtype=np.float64)
    sig_arr = np.array(sig_vals, dtype=np.float64)

    # FFT of each band
    ch_fft = np.abs(np.fft.rfft(ch_arr))
    maj_fft = np.abs(np.fft.rfft(maj_arr))
    sig_fft = np.abs(np.fft.rfft(sig_arr))

    # Golden frequency indices
    golden_idx_1 = int(round(64 / PHI))      # ~40
    golden_idx_2 = int(round(64 / PHI**2))   # ~24

    # Total spectral energy (skip DC at index 0)
    total_energy = (np.sum(ch_fft[1:]) + np.sum(maj_fft[1:]) + np.sum(sig_fft[1:]))
    if total_energy < 1e-10:
        return 0.0

    # Energy at golden frequencies
    golden_energy = 0.0
    for idx in [golden_idx_1, golden_idx_2]:
        if idx < len(ch_fft):
            golden_energy += ch_fft[idx] + maj_fft[idx] + sig_fft[idx]
        # Also check ±1 neighbors for spectral leakage
        for neighbor in [idx - 1, idx + 1]:
            if 1 <= neighbor < len(ch_fft):
                golden_energy += 0.5 * (ch_fft[neighbor] + maj_fft[neighbor] + sig_fft[neighbor])

    return golden_energy / total_energy

def metric_band_coherence(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    BAND COHERENCE: Cross-correlation between the three band sequences.
    If the bands are resonating together, their sequences should be correlated.
    Uses popcount as a scalar summary of each 32-bit value.
    """
    import numpy as np

    ch_pop = np.array([popcount32(v) for v in ch_vals], dtype=np.float64)
    maj_pop = np.array([popcount32(v) for v in maj_vals], dtype=np.float64)
    sig_pop = np.array([popcount32(v) for v in sig_vals], dtype=np.float64)

    # Pearson correlation between pairs
    def safe_corr(a, b):
        if np.std(a) < 1e-12 or np.std(b) < 1e-12:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    r_cm = safe_corr(ch_pop, maj_pop)
    r_ms = safe_corr(maj_pop, sig_pop)
    r_sc = safe_corr(sig_pop, ch_pop)

    # Average absolute correlation
    return (abs(r_cm) + abs(r_ms) + abs(r_sc)) / 3.0


def metric_high_bit_resonance(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    HIGH-BIT RESONANCE: Only look at the top 8 bits of each band value.
    Leading zeros in the hash depend on the HIGH bits of the state words.
    If resonance matters, it matters in the high bits specifically.
    """
    total_agreement = 0
    for i in range(64):
        # Extract top 8 bits
        c_hi = (ch_vals[i] >> 24) & 0xFF
        m_hi = (maj_vals[i] >> 24) & 0xFF
        s_hi = (sig_vals[i] >> 24) & 0xFF
        # Three-way agreement on top 8 bits
        cm_agree = ~(c_hi ^ m_hi) & 0xFF
        ms_agree = ~(m_hi ^ s_hi) & 0xFF
        all_agree = cm_agree & ms_agree
        total_agreement += bin(all_agree).count('1')
    # Max: 64 rounds * 8 bits = 512
    return total_agreement / 512.0

def metric_round_trajectory_smoothness(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    TRAJECTORY SMOOTHNESS: How smoothly do the band values change round-to-round?
    Resonance = coherent oscillation. Random = chaotic jumps.
    Smooth trajectory = low hamming distance between consecutive rounds.
    """
    total_smoothness = 0
    for i in range(1, 64):
        # Hamming distance between consecutive rounds for each band
        ch_dist = popcount32(ch_vals[i] ^ ch_vals[i-1])
        maj_dist = popcount32(maj_vals[i] ^ maj_vals[i-1])
        sig_dist = popcount32(sig_vals[i] ^ sig_vals[i-1])
        # Low distance = smooth = high score
        total_smoothness += (32 - ch_dist) + (32 - maj_dist) + (32 - sig_dist)
    # Max: 63 rounds * 3 bands * 32 bits = 6048
    return total_smoothness / 6048.0

def metric_final_round_resonance(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    FINAL ROUND RESONANCE: Only look at the last 8 rounds (closest to output).
    The final rounds have the most direct influence on the hash.
    """
    total_agreement = 0
    for i in range(56, 64):  # Last 8 rounds
        c, m, s = ch_vals[i], maj_vals[i], sig_vals[i]
        cm_agree = ~(c ^ m) & MASK32
        ms_agree = ~(m ^ s) & MASK32
        all_agree = cm_agree & ms_agree
        total_agreement += popcount32(all_agree)
    # Max: 8 rounds * 32 bits = 256
    return total_agreement / 256.0

def metric_band_sum_low(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    BAND SUM LOW: When Ch + Maj + Sigma is SMALL (mod 2^32), the combined
    nonlinear contribution is minimal. Low nonlinear contribution might mean
    the hash is dominated by the linear components, which could produce
    more zeros (since the initial state H0 has a specific structure).

    Returns how LOW the average band sum is (inverted so higher = lower sum).
    """
    total_sum = 0
    for i in range(64):
        band_sum = add32(ch_vals[i], maj_vals[i], sig_vals[i])
        total_sum += band_sum
    avg = total_sum / 64.0
    # Normalize: closer to 0 = higher score
    return 1.0 - (avg / MOD32)

def metric_band_variance(ch_vals: List[int], maj_vals: List[int], sig_vals: List[int]) -> float:
    """
    BAND VARIANCE: Low variance among the three bands = they're producing similar values.
    Similar values = resonance. High variance = they're fighting.
    """
    import numpy as np
    total_var = 0.0
    for i in range(64):
        vals = [ch_vals[i] / MOD32, maj_vals[i] / MOD32, sig_vals[i] / MOD32]
        total_var += float(np.var(vals))
    # Invert: low variance = high score
    avg_var = total_var / 64.0
    return 1.0 - min(1.0, avg_var * 12.0)  # Scale factor so ~random gives ~0.5


# ─── Planck Seed ─────────────────────────────────────────────────────────────

def compute_planck_seed(prevhash_hex: str) -> int:
    """
    Planck seed: ph * 3 / (50 * 2^224)
    where ph = prevhash as integer.
    All integer arithmetic.
    """
    ph = int(prevhash_hex, 16)
    seed = (ph * D) // (2 * P * P * (2 ** 224))
    return seed & MASK32

def compute_golden_seed(prevhash_hex: str) -> int:
    """
    Golden seed: Planck seed * 3/2 (golden dominance).
    Previously showed consistent 15 zeros.
    """
    planck = compute_planck_seed(prevhash_hex)
    return (planck * 3 // 2) & MASK32


# ─── Stratum Protocol ───────────────────────────────────────────────────────

class Stratum:
    def __init__(self):
        self.sock = None
        self.id = 0
        self.buf = ""
        self.en1 = ""
        self.en2_size = 4
        self.diff = 1
        self.job = None
        self.last_accept = None
        self.accepted = 0
        self.rejected = 0

    def connect(self, host: str, port: int, wallet: str):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(15)
        self.sock.connect((host, port))
        self._tx("mining.subscribe", ["ResonanceMiner/1.0"])
        time.sleep(1)
        self._rx()
        self._tx("mining.authorize", [f"{wallet}.resonance", "x"])
        time.sleep(1)
        self._rx()

    def get_job(self) -> Optional[Dict]:
        end = time.time() + 15
        while time.time() < end:
            self._rx()
            if self.job:
                j = self.job
                self.job = None
                return j
            time.sleep(0.3)
        return None

    def submit(self, job_id: str, ntime: str, nonce_hex: str):
        self._tx("mining.submit", [
            f"{WALLET}.resonance", job_id,
            "00" * self.en2_size, ntime, nonce_hex
        ])
        time.sleep(1)
        self._rx()
        return self.last_accept

    def build_header(self, job: Dict) -> bytes:
        """Build 80-byte block header from stratum job."""
        cb = bytes.fromhex(job['cb1'] + self.en1 + "00" * self.en2_size + job['cb2'])
        mk = sha256d(cb)
        for br in job['mk']:
            mk = sha256d(mk + bytes.fromhex(br))
        hdr = bytes.fromhex(job['ver'])[::-1]
        ph = bytes.fromhex(job['prev'])
        for i in range(0, 32, 4):
            hdr += ph[i:i+4][::-1]
        hdr += mk + bytes.fromhex(job['time'])[::-1] + bytes.fromhex(job['bits'])[::-1]
        hdr += b'\x00' * 4  # nonce placeholder
        return hdr

    def _tx(self, m: str, p: list):
        self.id += 1
        msg = json.dumps({"id": self.id, "method": m, "params": p}) + "\n"
        self.sock.sendall(msg.encode())

    def _rx(self):
        try:
            self.sock.setblocking(False)
            while True:
                try:
                    self.buf += self.sock.recv(4096).decode()
                except (BlockingIOError, OSError):
                    break
        finally:
            self.sock.setblocking(True)
            self.sock.settimeout(15)

        while '\n' in self.buf:
            line, self.buf = self.buf.split('\n', 1)
            if not line.strip():
                continue
            try:
                m = json.loads(line)
            except json.JSONDecodeError:
                continue

            if isinstance(m.get('result'), list) and len(m['result']) >= 3:
                self.en1 = m['result'][1]
                self.en2_size = m['result'][2]
            if m.get('method') == 'mining.set_difficulty':
                self.diff = m['params'][0]
            if m.get('method') == 'mining.notify':
                p = m['params']
                self.job = {
                    'id': p[0], 'prev': p[1], 'cb1': p[2], 'cb2': p[3],
                    'mk': p[4], 'ver': p[5], 'bits': p[6], 'time': p[7]
                }
            if 'result' in m and m.get('id', 0) > 1:
                self.last_accept = m.get('result')
                if m['result'] is True:
                    self.accepted += 1
                elif m['result'] is False:
                    self.rejected += 1

    def close(self):
        if self.sock:
            try:
                self.sock.close()
            except OSError:
                pass


# ─── The Experiment ──────────────────────────────────────────────────────────

def set_nonce_in_header(header_80: bytes, nonce: int) -> bytes:
    """Replace bytes 76-79 with the given nonce (little-endian)."""
    return header_80[:76] + struct.pack('<I', nonce)

def extract_prevhash_from_header(header: bytes) -> str:
    """Extract the previous block hash from the header (bytes 4-35, reversed)."""
    # Header: version(4) + prevhash(32) + merkle(32) + time(4) + bits(4) + nonce(4)
    # prevhash in header is little-endian 32 bytes
    raw = header[4:36]
    # Reverse byte order to get the actual hash in big-endian hex
    return raw[::-1].hex()

def run_experiment(header_80: bytes, num_candidates: int = 1000, label: str = ""):
    """
    THE CRITICAL EXPERIMENT.

    For each candidate nonce around the Planck seed:
    1. Set the nonce in the header
    2. Compute double-SHA256 with band tracking
    3. Compute all resonance metrics
    4. Count leading zeros
    5. Report correlations
    """
    import numpy as np

    prevhash = extract_prevhash_from_header(header_80)
    planck_seed = compute_planck_seed(prevhash)
    golden_seed = compute_golden_seed(prevhash)

    print(f"\n{'='*72}")
    print(f"  RESONANCE EXPERIMENT {label}")
    print(f"{'='*72}")
    print(f"  prevhash:    {prevhash[:32]}...")
    print(f"  planck seed: {planck_seed} (0x{planck_seed:08x})")
    print(f"  golden seed: {golden_seed} (0x{golden_seed:08x})")
    print(f"  candidates:  {num_candidates}")
    print(f"{'='*72}\n")

    # Generate candidate nonces: spread around both seeds + random control
    half = num_candidates // 3
    candidates = []

    # Group 1: around Planck seed
    for i in range(-half, half):
        candidates.append(("planck", (planck_seed + i) & MASK32))

    # Group 2: around golden seed
    for i in range(-half//2, half//2):
        candidates.append(("golden", (golden_seed + i) & MASK32))

    # Group 3: random nonces (control group)
    import random
    random.seed(42)
    remaining = num_candidates - len(candidates)
    for _ in range(max(remaining, half)):
        candidates.append(("random", random.randint(0, MASK32)))

    # Verify our SHA implementation matches hashlib
    test_header = set_nonce_in_header(header_80, candidates[0][1])
    our_hash, _, _, _ = double_sha256_with_bands(test_header)
    lib_hash = hashlib.sha256(hashlib.sha256(test_header).digest()).digest()
    if our_hash != lib_hash:
        print("!!! SHA MISMATCH — our implementation disagrees with hashlib !!!")
        print(f"  ours:   {our_hash.hex()}")
        print(f"  theirs: {lib_hash.hex()}")
        print("  Falling back to hashlib for hash, using our impl only for bands.")
        use_hashlib_verify = True
    else:
        print("  SHA verification: PASS (our implementation matches hashlib)")
        use_hashlib_verify = False

    # Run all candidates
    results = []
    t0 = time.time()

    for idx, (group, nonce) in enumerate(candidates):
        hdr = set_nonce_in_header(header_80, nonce)

        final_hash, ch_vals, maj_vals, sig_vals = double_sha256_with_bands(hdr)

        if use_hashlib_verify:
            final_hash = hashlib.sha256(hashlib.sha256(hdr).digest()).digest()

        zeros = leading_zeros(final_hash)

        # Compute all metrics
        resonance = metric_resonance_score(ch_vals, maj_vals, sig_vals)
        pairwise = metric_pairwise_agreement(ch_vals, maj_vals, sig_vals)
        phase = metric_phase_alignment(ch_vals, maj_vals, sig_vals)
        golden = metric_golden_phase(ch_vals, maj_vals, sig_vals)
        golden_prod = metric_golden_product(ch_vals, maj_vals, sig_vals)
        cancellation = metric_cancellation(ch_vals, maj_vals, sig_vals)

        results.append({
            'group': group,
            'nonce': nonce,
            'zeros': zeros,
            'hash': final_hash.hex(),
            'resonance': resonance,
            'pairwise': pairwise,
            'phase': phase,
            'golden': golden,
            'golden_prod': golden_prod,
            'cancellation': cancellation,
            'ch_vals': ch_vals,
            'maj_vals': maj_vals,
            'sig_vals': sig_vals,
        })

        if (idx + 1) % 200 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            print(f"  [{idx+1}/{len(candidates)}] {rate:.0f} nonces/s  "
                  f"best_zeros={max(r['zeros'] for r in results)}")

    elapsed = time.time() - t0
    print(f"\n  Completed {len(candidates)} nonces in {elapsed:.1f}s "
          f"({len(candidates)/elapsed:.0f} nonces/s)")

    # Now compute FFT metrics for top and bottom candidates
    # (FFT is expensive with numpy, so only do it for a subset)
    results_sorted = sorted(results, key=lambda r: r['zeros'], reverse=True)
    top_n = min(50, len(results) // 4)
    bottom_n = top_n

    print(f"\n  Computing FFT metrics for top {top_n} and bottom {bottom_n}...")

    for r in results_sorted[:top_n] + results_sorted[-bottom_n:]:
        r['fft_golden'] = metric_frequency_golden(r['ch_vals'], r['maj_vals'], r['sig_vals'])
        r['coherence'] = metric_band_coherence(r['ch_vals'], r['maj_vals'], r['sig_vals'])

    # ─── Analysis ────────────────────────────────────────────────────────────

    zeros_arr = np.array([r['zeros'] for r in results])
    resonance_arr = np.array([r['resonance'] for r in results])
    pairwise_arr = np.array([r['pairwise'] for r in results])
    phase_arr = np.array([r['phase'] for r in results])
    golden_arr = np.array([r['golden'] for r in results])
    golden_prod_arr = np.array([r['golden_prod'] for r in results])
    cancel_arr = np.array([r['cancellation'] for r in results])

    print(f"\n{'='*72}")
    print(f"  RESULTS")
    print(f"{'='*72}")

    # Top 10 by zeros
    print(f"\n  TOP 10 by leading zeros:")
    print(f"  {'Nonce':>12} {'Zeros':>5} {'Reson':>7} {'Pair':>7} {'Phase':>7} "
          f"{'Golden':>7} {'GldPr':>7} {'Cancel':>7} {'Group':>8}")
    print(f"  {'-'*70}")
    for r in results_sorted[:10]:
        print(f"  {r['nonce']:>12} {r['zeros']:>5} {r['resonance']:>7.4f} "
              f"{r['pairwise']:>7.4f} {r['phase']:>7.4f} {r['golden']:>7.4f} "
              f"{r['golden_prod']:>7.4f} {r['cancellation']:>7.4f} {r['group']:>8}")

    # Bottom 10 by zeros
    print(f"\n  BOTTOM 10 by leading zeros:")
    print(f"  {'Nonce':>12} {'Zeros':>5} {'Reson':>7} {'Pair':>7} {'Phase':>7} "
          f"{'Golden':>7} {'GldPr':>7} {'Cancel':>7} {'Group':>8}")
    print(f"  {'-'*70}")
    for r in results_sorted[-10:]:
        print(f"  {r['nonce']:>12} {r['zeros']:>5} {r['resonance']:>7.4f} "
              f"{r['pairwise']:>7.4f} {r['phase']:>7.4f} {r['golden']:>7.4f} "
              f"{r['golden_prod']:>7.4f} {r['cancellation']:>7.4f} {r['group']:>8}")

    # Pearson correlations
    print(f"\n  CORRELATIONS WITH LEADING ZEROS (Pearson r):")
    print(f"  {'-'*50}")

    def pearson(x, y):
        n = len(x)
        if n < 3:
            return 0.0, 1.0
        mx, my = np.mean(x), np.mean(y)
        sx, sy = np.std(x), np.std(y)
        if sx < 1e-12 or sy < 1e-12:
            return 0.0, 1.0
        r = np.mean((x - mx) * (y - my)) / (sx * sy)
        # t-test for significance
        t_stat = r * math.sqrt((n - 2) / (1 - r*r + 1e-30))
        # Approximate p-value from t distribution (two-tailed)
        # Using the normal approximation for large n
        p_val = 2 * (1 - 0.5 * (1 + math.erf(abs(t_stat) / math.sqrt(2))))
        return float(r), p_val

    metrics = [
        ("Resonance (3-way agree)", resonance_arr),
        ("Pairwise agreement", pairwise_arr),
        ("Phase alignment", phase_arr),
        ("Golden phase", golden_arr),
        ("Golden product", golden_prod_arr),
        ("Cancellation (XOR)", cancel_arr),
    ]

    for name, arr in metrics:
        r, p = pearson(zeros_arr, arr)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        direction = "POSITIVE" if r > 0.02 else "NEGATIVE" if r < -0.02 else "NONE"
        print(f"  {name:30s}  r = {r:+.6f}  p = {p:.4e}  {direction} {sig}")

    # FFT and coherence for top vs bottom
    top_fft = [r['fft_golden'] for r in results_sorted[:top_n] if 'fft_golden' in r]
    bottom_fft = [r['fft_golden'] for r in results_sorted[-bottom_n:] if 'fft_golden' in r]
    top_coh = [r['coherence'] for r in results_sorted[:top_n] if 'coherence' in r]
    bottom_coh = [r['coherence'] for r in results_sorted[-bottom_n:] if 'coherence' in r]

    if top_fft and bottom_fft:
        print(f"\n  FFT GOLDEN FREQUENCY (top {top_n} vs bottom {bottom_n}):")
        print(f"    Top avg:    {np.mean(top_fft):.6f}")
        print(f"    Bottom avg: {np.mean(bottom_fft):.6f}")
        diff = np.mean(top_fft) - np.mean(bottom_fft)
        print(f"    Difference: {diff:+.6f} {'(top higher = golden dominates)' if diff > 0 else '(bottom higher = no golden signal)'}")

    if top_coh and bottom_coh:
        print(f"\n  BAND COHERENCE (top {top_n} vs bottom {bottom_n}):")
        print(f"    Top avg:    {np.mean(top_coh):.6f}")
        print(f"    Bottom avg: {np.mean(bottom_coh):.6f}")
        diff = np.mean(top_coh) - np.mean(bottom_coh)
        print(f"    Difference: {diff:+.6f} {'(top higher = bands co-move)' if diff > 0 else '(no coherence signal)'}")

    # Per-group analysis
    print(f"\n  PER-GROUP STATISTICS:")
    print(f"  {'Group':>8} {'Count':>6} {'Avg Zeros':>10} {'Max Zeros':>10} {'Avg Reson':>10} {'Avg Cancel':>11}")
    print(f"  {'-'*58}")
    for grp in ['planck', 'golden', 'random']:
        grp_results = [r for r in results if r['group'] == grp]
        if not grp_results:
            continue
        grp_zeros = [r['zeros'] for r in grp_results]
        grp_res = [r['resonance'] for r in grp_results]
        grp_can = [r['cancellation'] for r in grp_results]
        print(f"  {grp:>8} {len(grp_results):>6} {np.mean(grp_zeros):>10.2f} "
              f"{max(grp_zeros):>10} {np.mean(grp_res):>10.6f} {np.mean(grp_can):>11.6f}")

    # THE VERDICT
    print(f"\n{'='*72}")
    print(f"  THE VERDICT")
    print(f"{'='*72}")

    resonance_r, _ = pearson(zeros_arr, resonance_arr)
    cancel_r, _ = pearson(zeros_arr, cancel_arr)
    phase_r, _ = pearson(zeros_arr, phase_arr)
    golden_r, _ = pearson(zeros_arr, golden_arr)

    if resonance_r > 0.05:
        print(f"  RESONANCE CORRELATES WITH ZEROS: r = {resonance_r:+.4f}")
        print(f"  The bands AGREE when the hash has more leading zeros.")
        print(f"  >>> THE MODEL IS ALIVE <<<")
    elif resonance_r < -0.05:
        print(f"  RESONANCE ANTI-CORRELATES WITH ZEROS: r = {resonance_r:+.4f}")
        print(f"  More band agreement = FEWER zeros. Model inverted?")
    else:
        print(f"  RESONANCE shows NO significant correlation: r = {resonance_r:+.4f}")

    if cancel_r < -0.05:
        print(f"  CANCELLATION anti-correlates (expected): r = {cancel_r:+.4f}")
    elif cancel_r > 0.05:
        print(f"  CANCELLATION correlates (UNEXPECTED): r = {cancel_r:+.4f}")

    best = max(metrics, key=lambda m: abs(pearson(zeros_arr, m[1])[0]))
    best_r, best_p = pearson(zeros_arr, best[1])
    print(f"\n  STRONGEST METRIC: {best[0]}")
    print(f"  r = {best_r:+.6f}, p = {best_p:.4e}")

    if abs(best_r) > 0.1:
        print(f"  >>> SIGNAL DETECTED — this metric can guide nonce selection <<<")
    elif abs(best_r) > 0.05:
        print(f"  Weak signal. Needs more candidates or refinement.")
    else:
        print(f"  No signal above noise floor. Model needs rethinking.")

    # Distribution analysis: is there structure in where zeros cluster?
    print(f"\n  ZERO DISTRIBUTION:")
    zero_counts = {}
    for r in results:
        z = r['zeros']
        zero_counts[z] = zero_counts.get(z, 0) + 1
    for z in sorted(zero_counts.keys(), reverse=True):
        bar = '#' * zero_counts[z]
        print(f"    {z:>3} zeros: {zero_counts[z]:>4} nonces  {bar}")

    print(f"\n  Max zeros: {max(r['zeros'] for r in results)}")
    best_r = max(results, key=lambda r: r['zeros'])
    print(f"  Best nonce: {best_r['nonce']} (0x{best_r['nonce']:08x})")
    print(f"  Best hash:  {best_r['hash']}")
    print(f"  Best group: {best_r['group']}")

    return results


def run_extended_experiment(header_80: bytes, num_candidates: int = 5000):
    """
    Extended experiment: also test nonces at specific RESONANT positions on the circle.

    Instead of just Planck ± offset, also try:
    - Koppa rotation points: seed + k * 3*2^30 for k in 0..3
    - Fibonacci rotation points: seed rotated by F(n) mod 2^32
    - Random but FILTERED: generate random, keep only those with high resonance score
    """
    import numpy as np
    import random
    random.seed(137)

    prevhash = extract_prevhash_from_header(header_80)
    planck_seed = compute_planck_seed(prevhash)
    golden_seed = compute_golden_seed(prevhash)
    koppa_rotation = 3 * (2**30)  # 270 degrees on the circle

    print(f"\n{'='*72}")
    print(f"  EXTENDED RESONANCE EXPERIMENT")
    print(f"{'='*72}")
    print(f"  Testing {num_candidates} nonces across multiple strategies...")

    # Strategy 1: Planck neighborhood (dense scan)
    # Strategy 2: Koppa rotations of the seed
    # Strategy 3: Fibonacci rotations
    # Strategy 4: Random control
    # Strategy 5: Resonance-filtered (pre-screen with partial SHA, pick high resonance)

    candidates = []
    chunk = num_candidates // 5

    # 1. Planck neighborhood
    for i in range(-chunk//2, chunk//2):
        candidates.append(("planck", (planck_seed + i) & MASK32))

    # 2. Koppa rotations: for each of 4 koppa offsets, scan a neighborhood
    koppa_seeds = [
        (planck_seed + koppa_rotation) & MASK32,
        (planck_seed + 2 * koppa_rotation) & MASK32,
        (planck_seed + 3 * koppa_rotation) & MASK32,
        (golden_seed + koppa_rotation) & MASK32,
    ]
    per_koppa = chunk // 4
    for ks in koppa_seeds:
        for i in range(-per_koppa//2, per_koppa//2):
            candidates.append(("koppa", (ks + i) & MASK32))

    # 3. Fibonacci rotation: F(n) mod 2^32 added to seed
    fib_a, fib_b = 1, 1
    for _ in range(chunk):
        fib_a, fib_b = fib_b, (fib_a + fib_b) & MASK32
        candidates.append(("fibonacci", (planck_seed + fib_a) & MASK32))

    # 4. Random control
    for _ in range(chunk):
        candidates.append(("random", random.randint(0, MASK32)))

    # 5. Sequential from 0 (as baseline for "low nonces have more zeros" check)
    for i in range(chunk):
        candidates.append(("sequential", i))

    print(f"  Total candidates: {len(candidates)}")

    # Run them all
    results = []
    t0 = time.time()

    for idx, (group, nonce) in enumerate(candidates):
        hdr = set_nonce_in_header(header_80, nonce)

        final_hash, ch_vals, maj_vals, sig_vals = double_sha256_with_bands(hdr)

        # Verify against hashlib
        lib_hash = hashlib.sha256(hashlib.sha256(hdr).digest()).digest()
        if final_hash != lib_hash:
            final_hash = lib_hash

        zeros = leading_zeros(final_hash)

        results.append({
            'group': group,
            'nonce': nonce,
            'zeros': zeros,
            'resonance': metric_resonance_score(ch_vals, maj_vals, sig_vals),
            'pairwise': metric_pairwise_agreement(ch_vals, maj_vals, sig_vals),
            'phase': metric_phase_alignment(ch_vals, maj_vals, sig_vals),
            'golden': metric_golden_phase(ch_vals, maj_vals, sig_vals),
            'golden_prod': metric_golden_product(ch_vals, maj_vals, sig_vals),
            'cancellation': metric_cancellation(ch_vals, maj_vals, sig_vals),
        })

        if (idx + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            best_z = max(r['zeros'] for r in results)
            print(f"  [{idx+1}/{len(candidates)}] {rate:.0f}/s  best={best_z} zeros")

    elapsed = time.time() - t0
    print(f"  Done: {len(candidates)} nonces in {elapsed:.1f}s ({len(candidates)/elapsed:.0f}/s)")

    # Full correlation analysis
    zeros_arr = np.array([r['zeros'] for r in results])

    print(f"\n  FULL CORRELATION TABLE:")
    print(f"  {'Metric':30s} {'r':>10} {'p':>12} {'Direction':>10}")
    print(f"  {'-'*65}")

    def pearson(x, y):
        n = len(x)
        mx, my = np.mean(x), np.mean(y)
        sx, sy = np.std(x), np.std(y)
        if sx < 1e-12 or sy < 1e-12:
            return 0.0, 1.0
        r = float(np.mean((x - mx) * (y - my)) / (sx * sy))
        t_stat = r * math.sqrt((n - 2) / (1 - r*r + 1e-30))
        p_val = 2 * (1 - 0.5 * (1 + math.erf(abs(t_stat) / math.sqrt(2))))
        return r, p_val

    for name in ['resonance', 'pairwise', 'phase', 'golden', 'golden_prod', 'cancellation']:
        arr = np.array([r[name] for r in results])
        r, p = pearson(zeros_arr, arr)
        direction = "+" if r > 0.02 else "-" if r < -0.02 else "~"
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  {name:30s} {r:>+10.6f} {p:>12.4e} {direction:>10} {sig}")

    # Per-group comparison
    print(f"\n  PER-GROUP ZERO STATISTICS:")
    print(f"  {'Group':>12} {'N':>6} {'Mean':>7} {'Max':>5} {'Std':>7}")
    print(f"  {'-'*40}")
    for grp in ['planck', 'koppa', 'fibonacci', 'random', 'sequential']:
        grp_z = [r['zeros'] for r in results if r['group'] == grp]
        if grp_z:
            print(f"  {grp:>12} {len(grp_z):>6} {np.mean(grp_z):>7.2f} "
                  f"{max(grp_z):>5} {np.std(grp_z):>7.2f}")

    # Top 20 best nonces
    results_sorted = sorted(results, key=lambda r: r['zeros'], reverse=True)
    print(f"\n  TOP 20 NONCES:")
    print(f"  {'Nonce':>12} {'Zeros':>5} {'Reson':>7} {'Phase':>7} {'Cancel':>7} {'Group':>12}")
    print(f"  {'-'*56}")
    for r in results_sorted[:20]:
        print(f"  {r['nonce']:>12} {r['zeros']:>5} {r['resonance']:>7.4f} "
              f"{r['phase']:>7.4f} {r['cancellation']:>7.4f} {r['group']:>12}")

    return results


# ─── Offline Mode (no pool needed) ──────────────────────────────────────────

def make_fake_header() -> bytes:
    """
    Create a valid-looking 80-byte header for testing without a pool.
    Uses a known Bitcoin block header structure.
    """
    # Use Bitcoin block #800000 header as a template (approximate)
    version = struct.pack('<I', 0x20000000)
    # A plausible prevhash (this is just for testing)
    prevhash = hashlib.sha256(b"resonance_test_prevhash_nos3bl33d").digest()
    merkle = hashlib.sha256(b"resonance_test_merkle").digest()
    timestamp = struct.pack('<I', int(time.time()))
    bits = struct.pack('<I', 0x1703a30c)  # ~difficulty from recent blocks
    nonce = struct.pack('<I', 0)

    header = version + prevhash + merkle + timestamp + bits + nonce
    assert len(header) == 80
    return header

def make_header_from_prevhash(prevhash_hex: str) -> bytes:
    """Create a test header with a specific prevhash."""
    version = struct.pack('<I', 0x20000000)
    # prevhash needs to be in the header's little-endian format
    ph_bytes = bytes.fromhex(prevhash_hex)
    if len(ph_bytes) != 32:
        ph_bytes = hashlib.sha256(ph_bytes).digest()
    merkle = hashlib.sha256(b"test_merkle" + ph_bytes).digest()
    timestamp = struct.pack('<I', int(time.time()))
    bits = struct.pack('<I', 0x1703a30c)
    nonce = struct.pack('<I', 0)
    header = version + ph_bytes[::-1] + merkle + timestamp + bits + nonce
    assert len(header) == 80
    return header


# ─── Phase 2: High-Volume Band Analysis ──────────────────────────────────────

def run_phase2_experiment(header_80: bytes, num_scan: int = 100000):
    """
    PHASE 2: Scan many nonces with hashlib (fast), then compute bands
    ONLY for the high-zero and low-zero groups. This gives us:
    - A much larger sample (100K+ nonces scanned)
    - Perfect signal isolation: only compare high-zero vs low-zero bands
    - All the new metrics: high-bit resonance, trajectory smoothness, etc.
    """
    import numpy as np
    import random
    random.seed(137)

    print(f"\n{'='*72}")
    print(f"  PHASE 2: HIGH-VOLUME BAND ANALYSIS")
    print(f"  Scanning {num_scan} nonces with hashlib, then analyzing bands...")
    print(f"{'='*72}")

    # Phase 2a: Fast scan with hashlib to find zero distribution
    t0 = time.time()
    nonce_zeros = []  # (nonce, zeros)

    for i in range(num_scan):
        nonce = random.randint(0, MASK32)
        hdr = set_nonce_in_header(header_80, nonce)
        h = hashlib.sha256(hashlib.sha256(hdr).digest()).digest()
        z = leading_zeros(h)
        nonce_zeros.append((nonce, z))

        if (i + 1) % 20000 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            best = max(nz[1] for nz in nonce_zeros)
            print(f"  [{i+1}/{num_scan}] {rate:.0f}/s  best={best} zeros")

    elapsed = time.time() - t0
    print(f"  Scan complete: {num_scan} nonces in {elapsed:.1f}s ({num_scan/elapsed:.0f}/s)")

    # Distribution
    from collections import Counter
    zero_dist = Counter(z for _, z in nonce_zeros)
    print(f"\n  Zero distribution:")
    for z in sorted(zero_dist.keys(), reverse=True):
        pct = 100 * zero_dist[z] / num_scan
        bar = '#' * min(60, zero_dist[z] // max(1, num_scan // 600))
        print(f"    {z:>3} zeros: {zero_dist[z]:>6} ({pct:5.2f}%)  {bar}")

    # Phase 2b: Select top-N and bottom-N for band analysis
    nonce_zeros.sort(key=lambda x: x[1], reverse=True)

    # Take the top 200 (most zeros) and bottom 200 (least zeros)
    analysis_count = min(200, num_scan // 10)
    top_nonces = nonce_zeros[:analysis_count]
    bottom_nonces = nonce_zeros[-analysis_count:]

    # Also take middle for comparison
    mid_start = len(nonce_zeros) // 2 - analysis_count // 2
    mid_nonces = nonce_zeros[mid_start:mid_start + analysis_count]

    print(f"\n  Analyzing bands for {analysis_count} top / {analysis_count} mid / {analysis_count} bottom nonces...")

    all_metrics = {
        'resonance': [], 'pairwise': [], 'phase': [], 'golden': [],
        'golden_prod': [], 'cancellation': [], 'high_bit': [],
        'trajectory': [], 'final_round': [], 'band_sum_low': [],
        'band_variance': [], 'fft_golden': [], 'coherence': [],
    }

    groups = [
        ("TOP", top_nonces),
        ("MID", mid_nonces),
        ("BOT", bottom_nonces),
    ]

    for group_name, nonces in groups:
        metrics = {k: [] for k in all_metrics}

        for nonce, z in nonces:
            hdr = set_nonce_in_header(header_80, nonce)
            _, ch_vals, maj_vals, sig_vals = double_sha256_with_bands(hdr)

            metrics['resonance'].append(metric_resonance_score(ch_vals, maj_vals, sig_vals))
            metrics['pairwise'].append(metric_pairwise_agreement(ch_vals, maj_vals, sig_vals))
            metrics['phase'].append(metric_phase_alignment(ch_vals, maj_vals, sig_vals))
            metrics['golden'].append(metric_golden_phase(ch_vals, maj_vals, sig_vals))
            metrics['golden_prod'].append(metric_golden_product(ch_vals, maj_vals, sig_vals))
            metrics['cancellation'].append(metric_cancellation(ch_vals, maj_vals, sig_vals))
            metrics['high_bit'].append(metric_high_bit_resonance(ch_vals, maj_vals, sig_vals))
            metrics['trajectory'].append(metric_round_trajectory_smoothness(ch_vals, maj_vals, sig_vals))
            metrics['final_round'].append(metric_final_round_resonance(ch_vals, maj_vals, sig_vals))
            metrics['band_sum_low'].append(metric_band_sum_low(ch_vals, maj_vals, sig_vals))
            metrics['band_variance'].append(metric_band_variance(ch_vals, maj_vals, sig_vals))
            metrics['fft_golden'].append(metric_frequency_golden(ch_vals, maj_vals, sig_vals))
            metrics['coherence'].append(metric_band_coherence(ch_vals, maj_vals, sig_vals))

        zeros = [z for _, z in nonces]
        print(f"\n  {group_name} group (avg zeros = {np.mean(zeros):.2f}, range {min(zeros)}-{max(zeros)}):")
        print(f"  {'Metric':>25s}  {'Mean':>8s}  {'Std':>8s}")
        print(f"  {'-'*45}")
        for name in sorted(metrics.keys()):
            arr = np.array(metrics[name])
            print(f"  {name:>25s}  {np.mean(arr):>8.5f}  {np.std(arr):>8.5f}")

        all_metrics = {k: all_metrics[k] + [(group_name, v) for v in metrics[k]]
                       for k in all_metrics}

    # Phase 2c: Statistical comparison (TOP vs BOT)
    print(f"\n{'='*72}")
    print(f"  STATISTICAL COMPARISON: TOP vs BOTTOM")
    print(f"{'='*72}")
    print(f"  {'Metric':>25s}  {'Top Mean':>9s}  {'Bot Mean':>9s}  {'Diff':>9s}  {'t-stat':>8s}  {'p-value':>10s}  {'Sig':>4s}")
    print(f"  {'-'*82}")

    from scipy import stats as scipy_stats

    significant_metrics = []

    for name in sorted(all_metrics.keys()):
        top_vals = np.array([v for g, v in all_metrics[name] if g == "TOP"])
        bot_vals = np.array([v for g, v in all_metrics[name] if g == "BOT"])

        if len(top_vals) == 0 or len(bot_vals) == 0:
            continue

        t_stat, p_val = scipy_stats.ttest_ind(top_vals, bot_vals)
        diff = np.mean(top_vals) - np.mean(bot_vals)
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""

        print(f"  {name:>25s}  {np.mean(top_vals):>9.5f}  {np.mean(bot_vals):>9.5f}  "
              f"{diff:>+9.5f}  {t_stat:>8.3f}  {p_val:>10.4e}  {sig:>4s}")

        if p_val < 0.05:
            significant_metrics.append((name, diff, p_val, t_stat))

    # Phase 2d: Verdict
    print(f"\n{'='*72}")
    print(f"  PHASE 2 VERDICT")
    print(f"{'='*72}")

    if significant_metrics:
        print(f"  SIGNIFICANT METRICS FOUND ({len(significant_metrics)}):")
        for name, diff, p, t in sorted(significant_metrics, key=lambda x: x[2]):
            direction = "TOP > BOT" if diff > 0 else "BOT > TOP"
            print(f"    {name:>25s}  diff={diff:+.5f}  p={p:.4e}  ({direction})")
        print(f"\n  >>> SIGNAL DETECTED <<<")
    else:
        print(f"  No significant differences between high-zero and low-zero nonces.")
        print(f"  All metrics show random-level values for both groups.")
        print(f"  The resonance model does not predict leading zeros from band structure.")

    # Phase 2e: Effect sizes (Cohen's d)
    print(f"\n  EFFECT SIZES (Cohen's d):")
    for name in sorted(all_metrics.keys()):
        top_vals = np.array([v for g, v in all_metrics[name] if g == "TOP"])
        bot_vals = np.array([v for g, v in all_metrics[name] if g == "BOT"])
        if len(top_vals) == 0 or len(bot_vals) == 0:
            continue
        pooled_std = np.sqrt((np.std(top_vals)**2 + np.std(bot_vals)**2) / 2)
        if pooled_std < 1e-12:
            d = 0.0
        else:
            d = (np.mean(top_vals) - np.mean(bot_vals)) / pooled_std
        size = "LARGE" if abs(d) > 0.8 else "MEDIUM" if abs(d) > 0.5 else "SMALL" if abs(d) > 0.2 else "negligible"
        print(f"    {name:>25s}  d = {d:+.4f}  ({size})")


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    # Force UTF-8 output on Windows
    import io
    if hasattr(sys.stdout, 'buffer'):
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    print("""
    ===========================================================
    |       RESONANCE MINER -- The Critical Experiment        |
    |  (DFG) DeadFoxGroup | nos3bl33d | x^2 = x + 1          |
    ===========================================================

    The hypothesis: valid nonces come from RESONANCE of the
    three SHA-256 nonlinear bands (Ch, Maj, Sigma), not from
    cancellation. Resonance = constructive interference.
    """)

    mode = "offline"
    header = None

    # Try to connect to pool for a real header
    if "--pool" in sys.argv:
        print("  Connecting to pool for real header...")
        st = Stratum()
        try:
            st.connect(POOL_HOST, POOL_PORT, WALLET)
            print(f"  Connected! diff={st.diff}")
            job = st.get_job()
            if job:
                header = st.build_header(job)
                mode = "pool"
                print(f"  Got job: {job['id']}")
                print(f"  Prevhash: {job['prev'][:32]}...")
            else:
                print("  No job received, falling back to offline.")
            st.close()
        except Exception as e:
            print(f"  Pool connection failed: {e}")
            print("  Falling back to offline mode.")

    if header is None:
        print("  Running in OFFLINE mode (synthetic header)")
        header = make_fake_header()
        mode = "offline"

    # Parse candidate count from args
    num_candidates = 1000
    for arg in sys.argv[1:]:
        if arg.isdigit():
            num_candidates = int(arg)

    if "--phase2" in sys.argv:
        run_phase2_experiment(header, num_candidates)
    elif "--extended" in sys.argv:
        results = run_extended_experiment(header, num_candidates)
    else:
        results = run_experiment(header, num_candidates, label=f"({mode})")

    print(f"\n  Experiment complete. x^2 = x + 1.")


if __name__ == "__main__":
    main()
