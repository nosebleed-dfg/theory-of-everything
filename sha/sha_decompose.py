"""
SHA_DECOMPOSE — decomposes SHA-256 into Ch/Maj/carry components; measures golden signal decay per round
nos3bl33d

Configs A-H toggle each nonlinear op independently.
Metrics: per-bit bias, inter-output correlation, spectral energy, Hamming distance, Cohen's d.
"""

import struct
import os
import sys
import time
import numpy as np
from scipy import stats

# ============================================================
# SHA-256 Constants
# ============================================================

SHA256_H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

SHA256_K = [
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

PHI = (1 + 5**0.5) / 2

# ============================================================
# Bit-level operations
# ============================================================

MASK32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return (x >> n) & MASK32

# ============================================================
# Message schedule
# ============================================================

def message_schedule(block_bytes):
    assert len(block_bytes) == 64
    W = list(struct.unpack('>16I', block_bytes))
    for i in range(16, 64):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ shr(W[i-15], 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ shr(W[i-2], 10)
        W.append((W[i-16] + s0 + W[i-7] + s1) & MASK32)
    return W

# ============================================================
# Modular SHA round
# ============================================================

def sha_round(state, k_i, w_i, use_ch=True, use_maj=True, use_carries=True):
    a, b, c, d, e, f, g, h = state

    S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

    if use_ch:
        ch = ((e & f) ^ ((~e & MASK32) & g)) & MASK32
    else:
        ch = (e ^ f ^ g) & MASK32

    S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

    if use_maj:
        maj = ((a & b) ^ (a & c) ^ (b & c)) & MASK32
    else:
        maj = (a ^ b ^ c) & MASK32

    if use_carries:
        T1 = (h + S1 + ch + k_i + w_i) & MASK32
        T2 = (S0 + maj) & MASK32
        new_e = (d + T1) & MASK32
        new_a = (T1 + T2) & MASK32
    else:
        T1 = (h ^ S1 ^ ch ^ k_i ^ w_i) & MASK32
        T2 = (S0 ^ maj) & MASK32
        new_e = (d ^ T1) & MASK32
        new_a = (T1 ^ T2) & MASK32

    return (new_a, a, b, c, new_e, e, f, g)


def sha_compress(block_bytes, num_rounds, use_ch=True, use_maj=True, use_carries=True):
    W = message_schedule(block_bytes)
    state = tuple(SHA256_H0)
    for i in range(min(num_rounds, 64)):
        state = sha_round(state, SHA256_K[i], W[i],
                          use_ch=use_ch, use_maj=use_maj, use_carries=use_carries)
    return state


def state_to_words(state):
    return np.array(state, dtype=np.uint64)


def state_to_bits(state):
    bits = []
    for word in state:
        for bit_pos in range(31, -1, -1):
            bits.append((word >> bit_pos) & 1)
    return np.array(bits, dtype=np.float64)


def state_to_bitvec(state):
    bits = []
    for word in state:
        for bit_pos in range(31, -1, -1):
            bits.append((word >> bit_pos) & 1)
    return np.array(bits, dtype=np.uint8)


# ============================================================
# Input generation
# ============================================================

def make_golden_inputs(n=5000):
    """Generate n golden-ratio-spaced 32-bit integers, packed as 64-byte blocks."""
    blocks = []
    for i in range(n):
        val = int(i * PHI * (2**32)) % (2**32)
        block = bytearray(64)
        block[0:4] = struct.pack('>I', val)
        block[4] = 0x80
        block[62] = 0x00
        block[63] = 0x20
        blocks.append(bytes(block))
    return blocks


def make_random_inputs(n=5000, seed=42):
    """Generate n random 64-byte blocks with proper SHA padding."""
    rng = np.random.RandomState(seed)
    blocks = []
    for i in range(n):
        val = int(rng.randint(0, 2**31 - 1)) * 2 + int(rng.randint(0, 2))
        block = bytearray(64)
        block[0:4] = struct.pack('>I', val & MASK32)
        block[4] = 0x80
        block[62] = 0x00
        block[63] = 0x20
        blocks.append(bytes(block))
    return blocks


# ============================================================
# Metric computation
# ============================================================

def compute_all_states(blocks, num_rounds, use_ch, use_maj, use_carries):
    """Run SHA on all blocks, return array of states and bit arrays."""
    word_array = np.zeros((len(blocks), 8), dtype=np.uint64)
    bit_matrix = np.zeros((len(blocks), 256), dtype=np.float64)

    for idx, block in enumerate(blocks):
        state = sha_compress(block, num_rounds, use_ch, use_maj, use_carries)
        word_array[idx] = state_to_words(state)
        bit_matrix[idx] = state_to_bits(state)

    return word_array, bit_matrix


def measure_signal(golden_words, golden_bits, random_words, random_bits):
    """
    Compute multiple signal strength metrics comparing golden vs random outputs.

    Returns dict of metric_name -> signal_strength (Cohen's d or equivalent).
    """
    n_golden = golden_words.shape[0]
    n_random = random_words.shape[0]

    metrics = {}

    # === METRIC 1: Bit bias per position ===
    # For each of 256 bit positions, compute the mean across all samples.
    # For truly random outputs, each bit should be ~0.5.
    # Golden inputs might bias specific bits.
    golden_bit_means = np.mean(golden_bits, axis=0)  # shape (256,)
    random_bit_means = np.mean(random_bits, axis=0)
    # Signal = max |difference| across all bit positions, normalized
    # Use chi-squared-like: sum of squared deviations
    golden_bias = np.sum((golden_bit_means - 0.5)**2)
    random_bias = np.sum((random_bit_means - 0.5)**2)
    # Bootstrap std of random bias
    n_boot = 100
    rng = np.random.RandomState(123)
    random_biases = []
    for _ in range(n_boot):
        idx = rng.choice(n_random, n_random, replace=True)
        rb = np.sum((np.mean(random_bits[idx], axis=0) - 0.5)**2)
        random_biases.append(rb)
    std_rb = np.std(random_biases)
    if std_rb > 1e-12:
        metrics['bit_bias'] = abs(golden_bias - random_bias) / std_rb
    else:
        metrics['bit_bias'] = 0.0

    # === METRIC 2: Hamming distance between consecutive outputs ===
    golden_hamming = np.sum(np.abs(golden_bits[1:] - golden_bits[:-1]), axis=1)
    random_hamming = np.sum(np.abs(random_bits[1:] - random_bits[:-1]), axis=1)
    std_rh = np.std(random_hamming)
    if std_rh > 1e-12:
        metrics['hamming_consecutive'] = abs(np.mean(golden_hamming) - np.mean(random_hamming)) / std_rh
    else:
        metrics['hamming_consecutive'] = 0.0

    # === METRIC 3: Word-level spectral analysis ===
    # Take the first word (a) of each output as a sequence of uint32 values.
    # Compute DFT, look for peaks at golden-ratio-related frequencies.
    golden_a = golden_words[:, 0].astype(np.float64)
    random_a = random_words[:, 0].astype(np.float64)

    # Normalize to [0, 1]
    golden_norm = golden_a / (2**32)
    random_norm = random_a / (2**32)

    # DFT
    golden_fft = np.abs(np.fft.rfft(golden_norm - np.mean(golden_norm)))
    random_fft = np.abs(np.fft.rfft(random_norm - np.mean(random_norm)))

    # Look at the energy ratio: total golden spectral energy vs random
    golden_energy = np.sum(golden_fft**2)
    random_energy = np.sum(random_fft**2)

    # More specifically: look at the peak-to-mean ratio
    if len(golden_fft) > 1:
        golden_peak_ratio = np.max(golden_fft[1:]) / (np.mean(golden_fft[1:]) + 1e-30)
        random_peak_ratio = np.max(random_fft[1:]) / (np.mean(random_fft[1:]) + 1e-30)

        # Bootstrap std of random peak ratio
        random_prs = []
        for _ in range(n_boot):
            idx = rng.choice(n_random, n_random, replace=True)
            rn = random_norm[idx]
            rf = np.abs(np.fft.rfft(rn - np.mean(rn)))
            if len(rf) > 1:
                random_prs.append(np.max(rf[1:]) / (np.mean(rf[1:]) + 1e-30))
        std_pr = np.std(random_prs) if random_prs else 1e-12
        if std_pr > 1e-12:
            metrics['spectral_peak'] = abs(golden_peak_ratio - random_peak_ratio) / std_pr
        else:
            metrics['spectral_peak'] = 0.0
    else:
        metrics['spectral_peak'] = 0.0

    # === METRIC 4: Autocorrelation of the OUTPUT SEQUENCE (not per-output bits) ===
    # Treat the sequence of first-word outputs as a time series.
    # Compute lag-1 autocorrelation of this sequence.
    def seq_autocorr(x, lag=1):
        if len(x) < lag + 1:
            return 0.0
        xm = x - np.mean(x)
        c0 = np.sum(xm**2)
        if c0 < 1e-30:
            return 0.0
        return np.sum(xm[:-lag] * xm[lag:]) / c0

    golden_ac1 = seq_autocorr(golden_norm, 1)
    random_ac1 = seq_autocorr(random_norm, 1)

    # Bootstrap for std
    random_acs = []
    for _ in range(n_boot):
        idx = rng.choice(n_random, n_random, replace=True)
        random_acs.append(seq_autocorr(random_norm[idx], 1))
    std_ac = np.std(random_acs)
    if std_ac > 1e-12:
        metrics['seq_autocorr'] = abs(golden_ac1 - random_ac1) / std_ac
    else:
        metrics['seq_autocorr'] = 0.0

    # === METRIC 5: Consecutive difference distribution ===
    # For golden inputs, consecutive inputs differ by ~phi*2^32.
    # After N rounds of SHA, do consecutive OUTPUTS differ by a detectable amount?
    golden_diffs = np.diff(golden_a)  # raw uint32 differences
    random_diffs = np.diff(random_a)
    # Look at the variance: golden diffs might be more uniform (less variance)
    golden_diff_var = np.var(golden_diffs)
    random_diff_var = np.var(random_diffs)

    random_dvars = []
    for _ in range(n_boot):
        idx = rng.choice(n_random, n_random, replace=True)
        random_dvars.append(np.var(np.diff(random_a[idx])))
    std_dv = np.std(random_dvars)
    if std_dv > 1e-12:
        metrics['diff_variance'] = abs(golden_diff_var - random_diff_var) / std_dv
    else:
        metrics['diff_variance'] = 0.0

    # === METRIC 6: KS test on output word values ===
    # Are the distributions of output words different?
    ks_stat, ks_p = stats.ks_2samp(golden_a, random_a)
    metrics['ks_pvalue'] = ks_p
    metrics['ks_stat'] = ks_stat
    # Convert KS to signal: -log10(p) as a measure, capped
    if ks_p > 0:
        metrics['ks_signal'] = min(-np.log10(ks_p), 50)
    else:
        metrics['ks_signal'] = 50

    # === METRIC 7: Bit-level mutual information proxy ===
    # For each bit position, compute the XOR between consecutive outputs.
    # In golden case, the XOR pattern might be biased.
    golden_xor = np.abs(golden_bits[1:] - golden_bits[:-1])  # 0 or 1
    random_xor = np.abs(random_bits[1:] - random_bits[:-1])
    # Mean XOR per bit position
    golden_xor_mean = np.mean(golden_xor, axis=0)
    random_xor_mean = np.mean(random_xor, axis=0)
    # For random, each XOR bit should be ~0.5
    golden_xor_bias = np.sum((golden_xor_mean - 0.5)**2)
    random_xor_bias = np.sum((random_xor_mean - 0.5)**2)

    random_xbs = []
    for _ in range(n_boot):
        idx = rng.choice(n_random - 1, n_random - 1, replace=True)
        rxm = np.mean(random_xor[idx], axis=0)
        random_xbs.append(np.sum((rxm - 0.5)**2))
    std_xb = np.std(random_xbs)
    if std_xb > 1e-12:
        metrics['xor_bias'] = abs(golden_xor_bias - random_xor_bias) / std_xb
    else:
        metrics['xor_bias'] = 0.0

    # === COMPOSITE SIGNAL ===
    # Average of the main Cohen's d metrics (excluding ks which is on different scale)
    cohens = [metrics['bit_bias'], metrics['hamming_consecutive'],
              metrics['spectral_peak'], metrics['seq_autocorr'],
              metrics['diff_variance'], metrics['xor_bias']]
    metrics['composite'] = np.mean(cohens)
    metrics['max_signal'] = np.max(cohens)

    return metrics


# ============================================================
# Configurations
# ============================================================

CONFIGS = {
    'A': {'use_ch': False, 'use_maj': False, 'use_carries': False, 'label': 'skeleton only'},
    'B': {'use_ch': True,  'use_maj': False, 'use_carries': False, 'label': 'skeleton + Ch'},
    'C': {'use_ch': False, 'use_maj': True,  'use_carries': False, 'label': 'skeleton + Maj'},
    'D': {'use_ch': False, 'use_maj': False, 'use_carries': True,  'label': 'skeleton + carries'},
    'E': {'use_ch': True,  'use_maj': True,  'use_carries': False, 'label': 'Ch + Maj, no carries'},
    'F': {'use_ch': True,  'use_maj': False, 'use_carries': True,  'label': 'Ch + carries, no Maj'},
    'G': {'use_ch': False, 'use_maj': True,  'use_carries': True,  'label': 'Maj + carries, no Ch'},
    'H': {'use_ch': True,  'use_maj': True,  'use_carries': True,  'label': 'full SHA-256'},
}

ROUND_COUNTS = [1, 2, 4, 8, 16, 32, 64]

# ============================================================
# Forward-Reverse Residual Analysis
# ============================================================

def build_linear_round_matrix():
    """
    Build the 256x256 GF(2) matrix representing one round of the fully-linearized
    SHA-256 compression (use_ch=False, use_maj=False, use_carries=False).
    """
    M = np.zeros((256, 256), dtype=np.uint8)

    def idx(word, bit):
        return word * 32 + bit

    # Word indices: a=0, b=1, c=2, d=3, e=4, f=5, g=6, h=7

    for b in range(32):
        # T1[b] sources: h[b], S1(e)[b], ch_lin(e,f,g)[b]
        # S1(e)[b] = e[(b+6)%32] ^ e[(b+11)%32] ^ e[(b+25)%32]
        # ch_lin[b] = e[b] ^ f[b] ^ g[b]
        t1_sources = []
        t1_sources.append(idx(7, b))  # h
        t1_sources.append(idx(4, (b+6)%32))   # S1
        t1_sources.append(idx(4, (b+11)%32))  # S1
        t1_sources.append(idx(4, (b+25)%32))  # S1
        t1_sources.append(idx(4, b))   # ch_lin: e
        t1_sources.append(idx(5, b))   # ch_lin: f
        t1_sources.append(idx(6, b))   # ch_lin: g

        t1_counts = {}
        for s in t1_sources:
            t1_counts[s] = t1_counts.get(s, 0) + 1
        t1_active = {s for s, c in t1_counts.items() if c % 2 == 1}

        # T2[b] sources: S0(a)[b], maj_lin(a,b,c)[b]
        # S0(a)[b] = a[(b+2)%32] ^ a[(b+13)%32] ^ a[(b+22)%32]
        # maj_lin[b] = a[b] ^ b[b] ^ c[b]
        t2_sources = []
        t2_sources.append(idx(0, (b+2)%32))   # S0
        t2_sources.append(idx(0, (b+13)%32))  # S0
        t2_sources.append(idx(0, (b+22)%32))  # S0
        t2_sources.append(idx(0, b))   # maj_lin: a
        t2_sources.append(idx(1, b))   # maj_lin: b
        t2_sources.append(idx(2, b))   # maj_lin: c

        t2_counts = {}
        for s in t2_sources:
            t2_counts[s] = t2_counts.get(s, 0) + 1
        t2_active = {s for s, c in t2_counts.items() if c % 2 == 1}

        # new_a = T1 ^ T2 (word 0)
        new_a_counts = {}
        for s in t1_active:
            new_a_counts[s] = new_a_counts.get(s, 0) + 1
        for s in t2_active:
            new_a_counts[s] = new_a_counts.get(s, 0) + 1
        for s, c in new_a_counts.items():
            if c % 2 == 1:
                M[idx(0, b), s] = 1

        # new_e = d ^ T1 (word 4)
        new_e_counts = {}
        new_e_counts[idx(3, b)] = 1
        for s in t1_active:
            new_e_counts[s] = new_e_counts.get(s, 0) + 1
        for s, c in new_e_counts.items():
            if c % 2 == 1:
                M[idx(4, b), s] = 1

    # Shift registers
    for b in range(32):
        M[idx(1, b), idx(0, b)] = 1  # new_b = old_a
        M[idx(2, b), idx(1, b)] = 1  # new_c = old_b
        M[idx(3, b), idx(2, b)] = 1  # new_d = old_c
        M[idx(5, b), idx(4, b)] = 1  # new_f = old_e
        M[idx(6, b), idx(5, b)] = 1  # new_g = old_f
        M[idx(7, b), idx(6, b)] = 1  # new_h = old_g

    return M


def gf2_mat_pow(M, n):
    size = M.shape[0]
    result = np.eye(size, dtype=np.uint8)
    base = M.copy()
    while n > 0:
        if n % 2 == 1:
            result = np.dot(result, base) % 2
        base = np.dot(base, base) % 2
        n //= 2
    return result


def gf2_mat_inv(M):
    n = M.shape[0]
    aug = np.zeros((n, 2*n), dtype=np.uint8)
    aug[:, :n] = M.copy()
    aug[:, n:] = np.eye(n, dtype=np.uint8)
    for col in range(n):
        pivot_row = None
        for row in range(col, n):
            if aug[row, col] == 1:
                pivot_row = row
                break
        if pivot_row is None:
            return None
        if pivot_row != col:
            aug[[col, pivot_row]] = aug[[pivot_row, col]]
        for row in range(n):
            if row != col and aug[row, col] == 1:
                aug[row] = (aug[row] + aug[col]) % 2
    return aug[:, n:]


def forward_reverse_test(golden_blocks):
    print("\n" + "="*70)
    print("FORWARD-REVERSE RESIDUAL ANALYSIS")
    print("="*70)
    print("Forward N rounds (full SHA), invert with linear matrix M^(-N).")
    print("Residual = Hamming(linear_inverse(full_output), linear_inverse(linear_output))")
    print("         = the bits flipped by nonlinearity alone.")
    print()

    M1 = build_linear_round_matrix()
    M1_inv = gf2_mat_inv(M1)

    if M1_inv is None:
        print("Linear round matrix is SINGULAR over GF(2).")
        rank = np.linalg.matrix_rank(M1.astype(np.float64))
        print(f"Rank: {rank}/256 (loses {256-rank} bits per round)")
        print()
        print("Fallback: comparing full SHA output vs linear skeleton output directly.")
        print(f"{'Rounds':>6} | {'Mean Hamming':>12} | {'Std':>8} | {'% of 256':>10}")
        print("-" * 50)

        for nr in [1, 2, 4, 8, 16, 32]:
            diffs = []
            for block in golden_blocks[:1000]:
                sf = sha_compress(block, nr, True, True, True)
                sl = sha_compress(block, nr, False, False, False)
                diff = sum(bin(a ^ b).count('1') for a, b in zip(sf, sl))
                diffs.append(diff)
            diffs = np.array(diffs)
            print(f"{nr:>6} | {np.mean(diffs):>12.2f} | {np.std(diffs):>8.2f} | {100*np.mean(diffs)/256:>9.1f}%")
        print()
        return

    print("Linear round matrix is INVERTIBLE over GF(2).")
    print()

    # Also do: full_SHA vs linear for each individual nonlinear component
    print("--- Full SHA vs Linear Skeleton (direct comparison) ---")
    print(f"{'Rounds':>6} | {'Mean Hamming':>12} | {'Std':>8} | {'% bits':>8} | {'Growth':>8}")
    print("-" * 55)

    prev_mean = None
    for nr in [1, 2, 4, 8, 16, 32]:
        diffs = []
        for block in golden_blocks[:1000]:
            sf = sha_compress(block, nr, True, True, True)
            sl = sha_compress(block, nr, False, False, False)
            diff = sum(bin(a ^ b).count('1') for a, b in zip(sf, sl))
            diffs.append(diff)
        diffs = np.array(diffs)
        mean_d = np.mean(diffs)
        growth = f"{mean_d/prev_mean:.2f}x" if prev_mean and prev_mean > 0 else ""
        prev_mean = mean_d
        print(f"{nr:>6} | {mean_d:>12.2f} | {np.std(diffs):>8.2f} | {100*mean_d/256:>7.1f}% | {growth:>8}")

    print()

    # Per-component contribution to nonlinearity
    print("--- Per-component nonlinear contribution (Hamming from skeleton) ---")
    components = [
        ('B', 'Ch only',     True,  False, False),
        ('C', 'Maj only',    False, True,  False),
        ('D', 'Carries only',False, False, True),
        ('F', 'Ch+Carries',  True,  False, True),
        ('G', 'Maj+Carries', False, True,  True),
        ('H', 'Full SHA',    True,  True,  True),
    ]

    for nr in [1, 2, 4, 8]:
        print(f"\n  Round {nr}:")
        for label_short, label, uc, um, ucar in components:
            diffs = []
            for block in golden_blocks[:1000]:
                sc = sha_compress(block, nr, uc, um, ucar)
                sl = sha_compress(block, nr, False, False, False)
                diff = sum(bin(a ^ b).count('1') for a, b in zip(sc, sl))
                diffs.append(diff)
            print(f"    {label:<15}: mean={np.mean(diffs):>7.2f} bits ({100*np.mean(diffs)/256:>5.1f}%)")

    print()

    # Residual via matrix inverse
    print("--- Residual via M^(-N) inverse ---")
    print(f"{'Rounds':>6} | {'Mean Residual':>14} | {'Std':>8} | {'% bits':>8} | {'Growth':>8}")
    print("-" * 60)

    prev_mean = None
    for nr in [1, 2, 4, 8, 16, 32]:
        M_inv_N = gf2_mat_pow(M1_inv, nr)
        residuals = []
        for block in golden_blocks[:1000]:
            sf = sha_compress(block, nr, True, True, True)
            sl = sha_compress(block, nr, False, False, False)
            bits_full = state_to_bitvec(sf)
            bits_lin = state_to_bitvec(sl)
            rec_full = np.dot(M_inv_N, bits_full) % 2
            rec_lin = np.dot(M_inv_N, bits_lin) % 2
            residual = int(np.sum(rec_full != rec_lin))
            residuals.append(residual)
        residuals = np.array(residuals)
        mean_r = np.mean(residuals)
        growth = f"{mean_r/prev_mean:.2f}x" if prev_mean and prev_mean > 0 else ""
        prev_mean = mean_r
        print(f"{nr:>6} | {mean_r:>14.2f} | {np.std(residuals):>8.2f} | {100*mean_r/256:>7.1f}% | {growth:>8}")

    print()


# ============================================================
# Main experiment
# ============================================================

def run_experiment():
    print("="*70)
    print("SHA-256 DECOMPOSITION: Golden Signal Decay Per Component")
    print("="*70)
    print()

    t_start = time.time()

    N_SAMPLES = 5000
    print(f"Generating {N_SAMPLES} golden + {N_SAMPLES} random inputs...")
    golden_blocks = make_golden_inputs(N_SAMPLES)
    random_blocks = make_random_inputs(N_SAMPLES, seed=42)
    print()

    # Storage: results[config][rounds] = metrics dict
    results = {}

    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        results[cfg_name] = {}
        print(f"Config {cfg_name}: {cfg['label']}")
        sys.stdout.flush()

        for nr in ROUND_COUNTS:
            gw, gb = compute_all_states(golden_blocks, nr,
                                        cfg['use_ch'], cfg['use_maj'], cfg['use_carries'])
            rw, rb = compute_all_states(random_blocks, nr,
                                        cfg['use_ch'], cfg['use_maj'], cfg['use_carries'])

            m = measure_signal(gw, gb, rw, rb)
            results[cfg_name][nr] = m

            print(f"  r={nr:>2}: composite={m['composite']:.3f} max={m['max_signal']:.3f} "
                  f"hamming={m['hamming_consecutive']:.3f} spectral={m['spectral_peak']:.3f} "
                  f"seq_ac={m['seq_autocorr']:.3f} KS_p={m['ks_pvalue']:.2e}")
            sys.stdout.flush()

        print()

    # ============================================================
    # Summary Tables
    # ============================================================

    for metric_name, metric_label in [
        ('composite', 'COMPOSITE SIGNAL (mean of all Cohen\'s d metrics)'),
        ('max_signal', 'MAX SIGNAL (strongest single metric)'),
        ('hamming_consecutive', 'HAMMING DISTANCE (consecutive outputs)'),
        ('spectral_peak', 'SPECTRAL PEAK RATIO'),
        ('seq_autocorr', 'SEQUENCE AUTOCORRELATION'),
        ('bit_bias', 'BIT BIAS'),
        ('xor_bias', 'XOR BIAS (consecutive output bit-flip pattern)'),
        ('diff_variance', 'DIFFERENCE VARIANCE'),
    ]:
        print(f"\n{'='*70}")
        print(f"SIGNAL STRENGTH: {metric_label}")
        print(f"{'='*70}")
        print(f"{'Config':<6} {'Description':<25}", end="")
        for nr in ROUND_COUNTS:
            print(f"{'r='+str(nr):>8}", end="")
        print()
        print("-" * (31 + 8 * len(ROUND_COUNTS)))

        for cfg_name in sorted(CONFIGS.keys()):
            cfg = CONFIGS[cfg_name]
            print(f"{cfg_name:<6} {cfg['label']:<25}", end="")
            for nr in ROUND_COUNTS:
                val = results[cfg_name][nr][metric_name]
                if val >= 10:
                    print(f"{val:>8.1f}", end="")
                elif val >= 1:
                    print(f"{val:>8.2f}", end="")
                else:
                    print(f"{val:>8.3f}", end="")
            print()

    # ============================================================
    # KS p-value table
    # ============================================================

    print(f"\n{'='*70}")
    print("KS TEST p-values (output word distribution, <0.05 = significant)")
    print(f"{'='*70}")
    print(f"{'Config':<6} {'Description':<25}", end="")
    for nr in ROUND_COUNTS:
        print(f"{'r='+str(nr):>10}", end="")
    print()
    print("-" * (31 + 10 * len(ROUND_COUNTS)))

    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        print(f"{cfg_name:<6} {cfg['label']:<25}", end="")
        for nr in ROUND_COUNTS:
            val = results[cfg_name][nr]['ks_pvalue']
            if val < 0.001:
                print(f"{'<0.001':>10}", end="")
            else:
                print(f"{val:>10.4f}", end="")
        print()

    # ============================================================
    # Analysis
    # ============================================================

    print("\n" + "="*70)
    print("ANALYSIS")
    print("="*70)

    def cs(cfg_name, nr):
        return results[cfg_name][nr]['composite']

    # Q1: Steepest decay
    print("\n--- Q1: Which operation DAMPENS the most? (steepest decay) ---")
    print()
    decay_data = []
    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        sig_1 = cs(cfg_name, 1)
        sig_64 = cs(cfg_name, 64)
        # Also measure "area under curve" for smoother comparison
        auc = sum(cs(cfg_name, nr) for nr in ROUND_COUNTS) / len(ROUND_COUNTS)
        decay = sig_1 - sig_64
        rate = decay / sig_1 if sig_1 > 1e-6 else 0
        decay_data.append((cfg_name, sig_1, sig_64, decay, rate, auc))
        print(f"  {cfg_name} ({cfg['label']:<25}): sig_1={sig_1:.3f} sig_64={sig_64:.3f} "
              f"decay={decay:+.3f} rate={rate:.0%} auc={auc:.3f}")

    print()
    print("  Marginal damping vs skeleton (A) at round 64:")
    a_64 = cs('A', 64)
    for cfg_name, op_name in [('B','Ch'), ('C','Maj'), ('D','Carries')]:
        delta = a_64 - cs(cfg_name, 64)
        print(f"    {op_name}: {delta:+.4f} (positive = extra damping)")

    # Q2: Amplification
    print("\n--- Q2: Which operation AMPLIFIES? (signal grows initially) ---")
    print()
    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        signals = [cs(cfg_name, nr) for nr in ROUND_COUNTS]
        found = False
        for i in range(1, len(signals)):
            if signals[i] > signals[0] * 1.1:  # 10% growth threshold
                print(f"  {cfg_name} ({cfg['label']}): AMPLIFIES at round {ROUND_COUNTS[i]} "
                      f"({signals[0]:.3f} -> {signals[i]:.3f}, +{(signals[i]/signals[0]-1)*100:.0f}%)")
                found = True
                break
        if not found:
            peak_idx = np.argmax(signals)
            print(f"  {cfg_name} ({cfg['label']}): peak at round {ROUND_COUNTS[peak_idx]} "
                  f"= {signals[peak_idx]:.3f}")

    # Q3: Noise floor
    print("\n--- Q3: Round where composite signal drops below 0.1 ---")
    print()
    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        noise_round = "never"
        for nr in ROUND_COUNTS:
            if cs(cfg_name, nr) < 0.1:
                noise_round = str(nr)
                break
        print(f"  {cfg_name} ({cfg['label']:<25}): round {noise_round}")

    # Q4: Skeleton persistence
    print("\n--- Q4: Does skeleton (A) maintain signal at all rounds? ---")
    print()
    for nr in ROUND_COUNTS:
        print(f"  r={nr}: composite={cs('A', nr):.4f}", end="")
        if cs('A', nr) > 0.1:
            print("  [DETECTABLE]")
        else:
            print("  [noise]")

    # Q5: Oscillation
    print("\n--- Q5: Does any config show oscillation? ---")
    print()
    for cfg_name in sorted(CONFIGS.keys()):
        cfg = CONFIGS[cfg_name]
        signals = [cs(cfg_name, nr) for nr in ROUND_COUNTS]
        dirs = []
        for i in range(1, len(signals)):
            if signals[i] > signals[i-1] + 0.02:
                dirs.append('+')
            elif signals[i] < signals[i-1] - 0.02:
                dirs.append('-')
            else:
                dirs.append('=')
        pattern = ''.join(dirs)
        # Oscillation = at least one -+ transition (drop then rise)
        osc = '-+' in pattern
        flag = " << OSCILLATES" if osc else ""
        print(f"  {cfg_name} ({cfg['label']:<25}): {pattern}{flag}")
        if osc:
            print(f"      values: {['%.3f'%s for s in signals]}")

    # ============================================================
    # Forward-Reverse
    # ============================================================

    forward_reverse_test(golden_blocks)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed:.1f}s")


if __name__ == '__main__':
    run_experiment()
