"""
GOLDEN_RESONANCE_TEST — tests whether golden structure survives 64 rounds of nonlinear SHA-256
nos3bl33d

char_poly+1 = x^4*(x^2-x-1)*(x^2-x+1). 8 statistical tests, golden vs random inputs.
"""

import hashlib
import struct
import os
import time
import sys
import numpy as np
from scipy import stats
from collections import Counter

# ─── Constants ───────────────────────────────────────────────────────────────

PHI = (1 + 5**0.5) / 2  # 1.6180339887...
NUM_SAMPLES = 10000
HASH_BITS = 256
HASH_BYTES = 32

# Fibonacci offsets to check in autocorrelation
FIB_OFFSETS = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]
# Golden-spaced offsets: round(phi*n) for n=1..20
GOLDEN_OFFSETS = [round(PHI * n) for n in range(1, 21)]

# ─── Utility Functions ──────────────────────────────────────────────────────

def sha256_int(data: bytes) -> int:
    """SHA-256 hash of data, returned as a 256-bit integer."""
    return int.from_bytes(hashlib.sha256(data).digest(), 'big')

def sha256_bytes(data: bytes) -> bytes:
    """SHA-256 hash of data, returned as raw bytes."""
    return hashlib.sha256(data).digest()

def popcount(x: int) -> int:
    """Count the number of set bits in an integer."""
    return bin(x).count('1')

def hamming_distance(a: int, b: int) -> int:
    """Hamming distance between two integers (popcount of XOR)."""
    return popcount(a ^ b)

def int_to_bits(x: int, nbits: int = 256) -> np.ndarray:
    """Convert integer to numpy array of bits (MSB first)."""
    bits = np.zeros(nbits, dtype=np.float64)
    for i in range(nbits):
        bits[nbits - 1 - i] = (x >> i) & 1
    return bits

def rotate_int(x: int, r: int, nbits: int = 256) -> int:
    """Circular right rotation of x by r bits within nbits-width."""
    r = r % nbits
    mask = (1 << nbits) - 1
    return ((x >> r) | (x << (nbits - r))) & mask

# ─── Generate Fibonacci Numbers ─────────────────────────────────────────────

def generate_fibonacci(n: int) -> list:
    """Generate first n Fibonacci numbers: F_1=1, F_2=1, F_3=2, ..."""
    fibs = []
    a, b = 1, 1
    for _ in range(n):
        fibs.append(a)
        a, b = b, a + b
    return fibs

# ─── Generate Input Sets ────────────────────────────────────────────────────

def generate_inputs():
    """Generate three input sets and their SHA-256 outputs."""
    print("Generating input sets...")
    t0 = time.time()

    # SET A: Hash of first 10000 Fibonacci numbers
    # F_n can be enormous, so we encode as variable-length big-endian bytes
    fibs = generate_fibonacci(NUM_SAMPLES)
    set_a = []
    for f in fibs:
        # Convert to bytes: use minimal byte representation
        byte_len = max(1, (f.bit_length() + 7) // 8)
        data = f.to_bytes(byte_len, 'big')
        set_a.append(sha256_int(data))

    print(f"  SET A (Fibonacci): {len(set_a)} hashes [{time.time()-t0:.2f}s]")

    # SET B: Phi-spaced inputs: int(i * phi * 2^32) mod 2^32
    t1 = time.time()
    set_b = []
    for i in range(NUM_SAMPLES):
        val = int(i * PHI * (2**32)) % (2**32)
        data = struct.pack('>I', val)  # 4-byte big-endian uint32
        set_b.append(sha256_int(data))

    print(f"  SET B (phi-spaced): {len(set_b)} hashes [{time.time()-t1:.2f}s]")

    # SET C: Random inputs: hash(random_bytes(32))
    t2 = time.time()
    rng = np.random.default_rng(seed=42)  # reproducible
    set_c = []
    for _ in range(NUM_SAMPLES):
        data = rng.bytes(32)
        set_c.append(sha256_int(data))

    print(f"  SET C (random): {len(set_c)} hashes [{time.time()-t2:.2f}s]")
    print(f"  Total generation time: {time.time()-t0:.2f}s")

    return set_a, set_b, set_c


# ═══════════════════════════════════════════════════════════════════════════
# TEST 1 & 2: Autocorrelation at Golden Offsets
# ═══════════════════════════════════════════════════════════════════════════

def compute_autocorrelation(outputs: list, max_k: int = 500) -> np.ndarray:
    """
    Compute autocorrelation function.
    autocorrelation(k) = mean over i of: popcount(output[i] XOR output[i+k]) / 256
    For random: ~0.5 for all k.
    """
    n = len(outputs)
    ac = np.zeros(max_k)
    for k in range(1, max_k + 1):
        total = 0.0
        count = 0
        for i in range(n - k):
            total += popcount(outputs[i] ^ outputs[i + k]) / HASH_BITS
            count += 1
        ac[k - 1] = total / count if count > 0 else 0.5
    return ac

def test_autocorrelation(set_a, set_b, set_c):
    """TEST 2: Autocorrelation at golden offsets vs non-golden offsets."""
    print("\n" + "="*72)
    print("TEST 2: Autocorrelation at Golden Offsets")
    print("="*72)

    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        print(f"\n  Computing autocorrelation for SET {name}...")
        t0 = time.time()
        ac = compute_autocorrelation(outputs, max_k=500)
        print(f"  Done in {time.time()-t0:.2f}s")

        # Extract autocorrelation at Fibonacci offsets
        fib_ac = [ac[k - 1] for k in FIB_OFFSETS]
        # Extract at golden-spaced offsets
        gold_ac = [ac[k - 1] for k in GOLDEN_OFFSETS if k <= 500]
        # Non-golden offsets: everything that's NOT Fibonacci or golden-spaced
        special = set(FIB_OFFSETS) | set(GOLDEN_OFFSETS)
        non_special_ac = [ac[k - 1] for k in range(1, 501) if k not in special]

        results[name] = {
            'ac': ac,
            'fib_ac': fib_ac,
            'gold_ac': gold_ac,
            'non_special_ac': non_special_ac,
        }

        print(f"  Mean autocorrelation (all k): {np.mean(ac):.6f}")
        print(f"  Mean at Fibonacci offsets:     {np.mean(fib_ac):.6f}")
        print(f"  Mean at golden offsets:        {np.mean(gold_ac):.6f}")
        print(f"  Mean at non-special offsets:   {np.mean(non_special_ac):.6f}")
        print(f"  Std at Fibonacci offsets:      {np.std(fib_ac):.6f}")
        print(f"  Std at non-special offsets:    {np.std(non_special_ac):.6f}")

    # Statistical comparison: are Fibonacci-k autocorrelations systematically
    # different from non-Fibonacci-k?
    print("\n  --- Statistical Comparison ---")
    for name in results:
        fib_vals = results[name]['fib_ac']
        non_vals = results[name]['non_special_ac']

        # Two-sample t-test
        t_stat, p_val = stats.ttest_ind(fib_vals, non_vals, equal_var=False)
        # Mann-Whitney U test (non-parametric)
        u_stat, p_val_mw = stats.mannwhitneyu(fib_vals, non_vals, alternative='two-sided')

        effect_size = (np.mean(fib_vals) - np.mean(non_vals)) / np.std(non_vals) if np.std(non_vals) > 0 else 0

        print(f"\n  SET {name}:")
        print(f"    t-test: t={t_stat:.4f}, p={p_val:.6f}")
        print(f"    Mann-Whitney U: U={u_stat:.1f}, p={p_val_mw:.6f}")
        print(f"    Effect size (Cohen's d): {effect_size:.6f}")
        results[name]['p_ttest'] = p_val
        results[name]['p_mw'] = p_val_mw
        results[name]['effect_size'] = effect_size

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 3: Bit Distribution
# ═══════════════════════════════════════════════════════════════════════════

def test_bit_distribution(set_a, set_b, set_c):
    """TEST 3: Frequency of each bit position being 1."""
    print("\n" + "="*72)
    print("TEST 3: Bit Distribution")
    print("="*72)

    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        # Count frequency of 1 at each bit position
        bit_counts = np.zeros(HASH_BITS)
        for h in outputs:
            for bit in range(HASH_BITS):
                if (h >> bit) & 1:
                    bit_counts[bit] += 1

        frequencies = bit_counts / len(outputs)
        deviations = frequencies - 0.5
        max_dev = np.max(np.abs(deviations))
        mean_dev = np.mean(np.abs(deviations))
        rms_dev = np.sqrt(np.mean(deviations**2))

        # Under null hypothesis, each bit is Bernoulli(0.5, n=10000)
        # std of sample proportion = sqrt(0.25/n) = 0.005
        expected_std = np.sqrt(0.25 / len(outputs))

        # Chi-square test: are the bit frequencies consistent with fair coins?
        chi2_stat = np.sum(deviations**2) / (0.25 / len(outputs))
        # Under null: chi2 with 256 degrees of freedom
        p_val = 1 - stats.chi2.cdf(chi2_stat, df=HASH_BITS)

        results[name] = {
            'max_dev': max_dev,
            'mean_dev': mean_dev,
            'rms_dev': rms_dev,
            'expected_std': expected_std,
            'chi2': chi2_stat,
            'p_val': p_val,
            'frequencies': frequencies,
        }

        print(f"\n  SET {name}:")
        print(f"    Max deviation from 0.5:  {max_dev:.6f}  (expected std: {expected_std:.6f})")
        print(f"    Mean |deviation|:        {mean_dev:.6f}")
        print(f"    RMS deviation:           {rms_dev:.6f}")
        print(f"    Chi-square (256 dof):    {chi2_stat:.2f}")
        print(f"    p-value:                 {p_val:.6f}")

    # Compare max deviations
    print("\n  --- Comparison ---")
    print(f"    Max dev ratio A/C: {results['A (Fibonacci)']['max_dev'] / results['C (random)']['max_dev']:.4f}")
    print(f"    Max dev ratio B/C: {results['B (phi-spaced)']['max_dev'] / results['C (random)']['max_dev']:.4f}")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 4: XOR Differential
# ═══════════════════════════════════════════════════════════════════════════

def test_xor_differential(set_a, set_b, set_c):
    """TEST 4: Hamming weight distribution of consecutive XOR differentials."""
    print("\n" + "="*72)
    print("TEST 4: XOR Differential")
    print("="*72)

    all_hw = {}
    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        # Compute Hamming weights of output[i] XOR output[i+1]
        hw = []
        for i in range(len(outputs) - 1):
            hw.append(popcount(outputs[i] ^ outputs[i + 1]))

        hw = np.array(hw, dtype=np.float64)
        mean_hw = np.mean(hw)
        std_hw = np.std(hw)
        skew_hw = stats.skew(hw)
        kurt_hw = stats.kurtosis(hw)

        # Expected for random: Binomial(256, 0.5) -> mean=128, std=8
        expected_mean = HASH_BITS / 2
        expected_std = np.sqrt(HASH_BITS / 4)

        all_hw[name] = hw

        results[name] = {
            'mean': mean_hw,
            'std': std_hw,
            'skewness': skew_hw,
            'kurtosis': kurt_hw,
            'hw': hw,
        }

        print(f"\n  SET {name}:")
        print(f"    Mean Hamming weight:   {mean_hw:.4f}  (expected: {expected_mean:.1f})")
        print(f"    Std:                   {std_hw:.4f}  (expected: {expected_std:.4f})")
        print(f"    Skewness:              {skew_hw:.6f}  (expected: ~0)")
        print(f"    Kurtosis:              {kurt_hw:.6f}  (expected: ~0)")

        # Check for exact-zero differentials (identical consecutive hashes)
        n_zeros = int(np.sum(hw == 0))
        if n_zeros > 0:
            print(f"    ** WARNING: {n_zeros} zero-weight differential(s) (identical consecutive hashes)")
            print(f"    ** This inflates kurtosis. Likely cause: F_1=F_2=1 (same input)")
            # Report stats with zeros removed
            hw_clean = hw[hw > 0]
            if len(hw_clean) > 0:
                print(f"    ** After removing zeros: skew={stats.skew(hw_clean):.6f}, "
                      f"kurt={stats.kurtosis(hw_clean):.6f}")
                results[name]['skewness_clean'] = stats.skew(hw_clean)
                results[name]['kurtosis_clean'] = stats.kurtosis(hw_clean)
            results[name]['n_zeros'] = n_zeros
        else:
            results[name]['n_zeros'] = 0

    # PROPER comparison: golden sets vs random set using 2-sample KS test
    # This avoids the artifact of comparing discrete Binomial to continuous Normal
    hw_c = all_hw["C (random)"]
    print("\n  --- Golden vs Random (2-sample KS) ---")
    for name in ["A (Fibonacci)", "B (phi-spaced)"]:
        ks_stat, ks_p = stats.ks_2samp(all_hw[name], hw_c)
        # Also compare moments directly
        skew_diff = abs(results[name]['skewness'] - results['C (random)']['skewness'])
        kurt_diff = abs(results[name]['kurtosis'] - results['C (random)']['kurtosis'])

        results[name]['ks_stat'] = ks_stat
        results[name]['ks_p'] = ks_p

        print(f"\n  SET {name} vs C (random):")
        print(f"    2-sample KS: D={ks_stat:.6f}, p={ks_p:.6f}")
        print(f"    |Skewness diff|: {skew_diff:.6f}")
        print(f"    |Kurtosis diff|: {kurt_diff:.6f}")

        if results[name]['kurtosis'] > 1.0:
            print(f"    ** NOTE: Kurtosis {results[name]['kurtosis']:.4f} is extreme (heavy tails)")

    # Give random a self-reference p-value of 1.0
    results['C (random)']['ks_stat'] = 0.0
    results['C (random)']['ks_p'] = 1.0

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 5: Golden Angle Bit Rotation
# ═══════════════════════════════════════════════════════════════════════════

def test_golden_rotation(set_a, set_b, set_c):
    """TEST 5: Rotate output by golden angle (98 bits), XOR with original."""
    print("\n" + "="*72)
    print("TEST 5: Golden Angle Bit Rotation")
    print("="*72)

    golden_rot = round(256 / PHI**2)  # = 98
    print(f"  Golden rotation: {golden_rot} bits (256/phi^2 = {256/PHI**2:.4f})")

    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        matching_bits = []
        for h in outputs:
            rotated = rotate_int(h, golden_rot)
            xored = h ^ rotated
            # Matching bits = 256 - popcount(xor)
            matching = HASH_BITS - popcount(xored)
            matching_bits.append(matching)

        mb = np.array(matching_bits, dtype=np.float64)
        mean_mb = np.mean(mb)
        std_mb = np.std(mb)

        # Expected: 128 matching bits (50%)
        # Test if mean differs from 128
        t_stat, p_val = stats.ttest_1samp(mb, 128.0)

        results[name] = {
            'mean': mean_mb,
            'std': std_mb,
            't_stat': t_stat,
            'p_val': p_val,
        }

        print(f"\n  SET {name}:")
        print(f"    Mean matching bits: {mean_mb:.4f}  (expected: 128.0)")
        print(f"    Std:                {std_mb:.4f}")
        print(f"    Deviation from 128: {mean_mb - 128:.4f}")
        print(f"    t-test: t={t_stat:.4f}, p={p_val:.6f}")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 6: Spectral Analysis (FFT)
# ═══════════════════════════════════════════════════════════════════════════

def test_spectral(set_a, set_b, set_c):
    """TEST 6: FFT of bit vectors, check for peaks at golden frequencies."""
    print("\n" + "="*72)
    print("TEST 6: Spectral Analysis (FFT)")
    print("="*72)

    # Golden frequencies
    k_phi = round(256 / PHI)       # ~158
    k_phi2 = round(256 / PHI**2)   # ~98
    k_phi3 = round(256 / PHI**3)   # ~60
    print(f"  Golden frequency bins: k={k_phi} (1/phi), k={k_phi2} (1/phi^2), k={k_phi3} (1/phi^3)")

    golden_bins = {k_phi, k_phi2, k_phi3}

    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        # Average FFT magnitude over all outputs
        avg_fft = np.zeros(HASH_BITS // 2 + 1)
        for h in outputs:
            bits = int_to_bits(h, HASH_BITS)
            fft_mag = np.abs(np.fft.rfft(bits - 0.5))  # center around 0
            avg_fft += fft_mag

        avg_fft /= len(outputs)

        # Compare magnitudes at golden bins vs non-golden bins
        # Skip DC component (k=0) and Nyquist
        golden_mags = [avg_fft[k] for k in golden_bins if k < len(avg_fft)]
        non_golden_mags = [avg_fft[k] for k in range(1, len(avg_fft) - 1) if k not in golden_bins]

        mean_golden = np.mean(golden_mags)
        mean_non_golden = np.mean(non_golden_mags)
        std_non_golden = np.std(non_golden_mags)

        # Z-score of golden bins relative to non-golden distribution
        z_scores = [(avg_fft[k] - mean_non_golden) / std_non_golden for k in golden_bins if k < len(avg_fft)]

        results[name] = {
            'avg_fft': avg_fft,
            'mean_golden': mean_golden,
            'mean_non_golden': mean_non_golden,
            'z_scores': z_scores,
            'golden_mags': golden_mags,
        }

        print(f"\n  SET {name}:")
        print(f"    Mean FFT mag at golden bins:     {mean_golden:.6f}")
        print(f"    Mean FFT mag at non-golden bins:  {mean_non_golden:.6f}")
        print(f"    Ratio (golden/non-golden):        {mean_golden/mean_non_golden:.6f}")
        for k, z in zip(sorted(golden_bins), z_scores):
            print(f"    k={k}: mag={avg_fft[k]:.6f}, z-score={z:.4f}")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 7: Chi-Square Test
# ═══════════════════════════════════════════════════════════════════════════

def test_chi_square(set_a, set_b, set_c):
    """TEST 7: Chi-square test on consecutive pairs of output bytes as 16-bit values."""
    print("\n" + "="*72)
    print("TEST 7: Chi-Square Test (16-bit pairs)")
    print("="*72)

    results = {}
    for name, outputs in [("A (Fibonacci)", set_a), ("B (phi-spaced)", set_b), ("C (random)", set_c)]:
        # Extract all 16-bit values from each hash output
        # Each 256-bit output has 16 consecutive 16-bit values
        values_16bit = []
        for h in outputs:
            h_bytes = h.to_bytes(HASH_BYTES, 'big')
            for j in range(0, HASH_BYTES, 2):
                val = (h_bytes[j] << 8) | h_bytes[j + 1]
                values_16bit.append(val)

        # With 10000 hashes * 16 pairs = 160000 samples
        # across 65536 bins, most bins will have ~2.4 entries
        # Chi-square may not be reliable with expected < 5
        # Instead: use a coarser binning (256 bins from high byte)
        n_bins = 256
        counts = Counter()
        for v in values_16bit:
            counts[v >> 8] += 1  # Use high byte as bin

        observed = np.array([counts.get(i, 0) for i in range(n_bins)], dtype=np.float64)
        expected = len(values_16bit) / n_bins

        chi2_stat = np.sum((observed - expected)**2 / expected)
        dof = n_bins - 1
        p_val = 1 - stats.chi2.cdf(chi2_stat, df=dof)

        # Also do a full 65536-bin test with exact counts
        full_counts = Counter(values_16bit)
        full_observed = np.array([full_counts.get(i, 0) for i in range(65536)], dtype=np.float64)
        full_expected = len(values_16bit) / 65536
        full_chi2 = np.sum((full_observed - full_expected)**2 / full_expected)
        full_dof = 65535
        full_p = 1 - stats.chi2.cdf(full_chi2, df=full_dof)

        results[name] = {
            'chi2_256': chi2_stat,
            'p_256': p_val,
            'chi2_65536': full_chi2,
            'p_65536': full_p,
            'n_samples': len(values_16bit),
        }

        print(f"\n  SET {name} ({len(values_16bit)} 16-bit samples):")
        print(f"    256-bin chi-square:   {chi2_stat:.2f}  (dof={dof}, p={p_val:.6f})")
        print(f"    65536-bin chi-square: {full_chi2:.2f}  (dof={full_dof}, p={full_p:.6f})")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# TEST 8: Pisano Correlation
# ═══════════════════════════════════════════════════════════════════════════

def test_pisano(set_a):
    """TEST 8: Check for Pisano period correlation in Fibonacci hash outputs."""
    print("\n" + "="*72)
    print("TEST 8: Pisano Correlation")
    print("="*72)

    # Pisano period mod 2^8 = 768
    pisano_8 = 768

    # We only have 10000 samples, so we can check pairs (i, i+768) for i < 10000-768
    n = len(set_a)
    if n < pisano_8 + 100:
        print("  Not enough samples for Pisano test")
        return {'p_val': 1.0}

    # Extract first byte of each hash
    first_bytes = np.array([(h >> 248) & 0xFF for h in set_a], dtype=np.float64)

    # Correlation between first_bytes[i] and first_bytes[i + pisano_8]
    x = first_bytes[:n - pisano_8]
    y = first_bytes[pisano_8:n]

    # Pearson correlation
    corr, p_val = stats.pearsonr(x, y)

    # Also check at non-Pisano offsets for comparison
    control_offsets = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    control_corrs = []
    for offset in control_offsets:
        if n > offset:
            cx = first_bytes[:n - offset]
            cy = first_bytes[offset:n]
            cc, _ = stats.pearsonr(cx, cy)
            control_corrs.append(cc)

    # Also check at Pisano period for 4-bit: 24
    pisano_4 = 24
    x4 = first_bytes[:n - pisano_4]
    y4 = first_bytes[pisano_4:n]
    corr_4, p_4 = stats.pearsonr(x4, y4)

    mean_control = np.mean(control_corrs)
    std_control = np.std(control_corrs)

    results = {
        'pisano_768_corr': corr,
        'pisano_768_p': p_val,
        'pisano_24_corr': corr_4,
        'pisano_24_p': p_4,
        'control_mean': mean_control,
        'control_std': std_control,
        'control_corrs': control_corrs,
        'p_val': p_val,
    }

    print(f"\n  Pisano period pi(2^8) = {pisano_8}")
    print(f"    Correlation at offset {pisano_8}: r={corr:.6f}, p={p_val:.6f}")
    print(f"  Pisano period pi(2^4) = {pisano_4}")
    print(f"    Correlation at offset {pisano_4}: r={corr_4:.6f}, p={p_4:.6f}")
    print(f"  Control offsets mean correlation: {mean_control:.6f} +/- {std_control:.6f}")
    print(f"  Control correlations: {[f'{c:.6f}' for c in control_corrs]}")

    # Is Pisano correlation an outlier?
    if std_control > 0:
        z_score = (corr - mean_control) / std_control
        print(f"  Z-score of Pisano-768 vs controls: {z_score:.4f}")
        results['z_score'] = z_score
    else:
        results['z_score'] = 0.0

    return results


# ═══════════════════════════════════════════════════════════════════════════
# FINAL ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def final_analysis(results_autocorr, results_bits, results_xor,
                   results_rotation, results_spectral, results_chi2,
                   results_pisano):
    """Synthesize all results and determine if golden structure survives."""
    print("\n" + "="*72)
    print("FINAL ANALYSIS: Does Golden Structure Survive SHA-256?")
    print("="*72)

    significant_tests = []
    all_tests = []

    # --- Test 2: Autocorrelation ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        p = results_autocorr[name]['p_ttest']
        ef = results_autocorr[name]['effect_size']
        label = f"T2-Autocorr-{name}"
        sig = p < 0.01
        all_tests.append((label, p, ef, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 3: Bit Distribution ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        p = results_bits[name]['p_val']
        ef = results_bits[name]['max_dev'] / results_bits[name]['expected_std']
        label = f"T3-BitDist-{name}"
        sig = p < 0.01
        all_tests.append((label, p, ef, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 4: XOR Differential ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        p = results_xor[name]['ks_p']
        ef = results_xor[name]['ks_stat']
        label = f"T4-XOR-{name}"
        sig = p < 0.01
        all_tests.append((label, p, ef, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 5: Golden Rotation ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        p = results_rotation[name]['p_val']
        ef = abs(results_rotation[name]['mean'] - 128.0) / results_rotation[name]['std']
        label = f"T5-Rotation-{name}"
        sig = p < 0.01
        all_tests.append((label, p, ef, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 6: Spectral ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        z_scores = results_spectral[name]['z_scores']
        max_z = max(abs(z) for z in z_scores)
        # Convert max z-score to p-value (two-sided, Bonferroni corrected for 3 bins)
        p = 2 * stats.norm.sf(max_z) * 3
        p = min(p, 1.0)
        label = f"T6-Spectral-{name}"
        sig = p < 0.01
        all_tests.append((label, p, max_z, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 7: Chi-Square ---
    for name in ['A (Fibonacci)', 'B (phi-spaced)']:
        p = results_chi2[name]['p_256']
        ef = (results_chi2[name]['chi2_256'] - 255) / np.sqrt(2 * 255)  # normalized
        label = f"T7-ChiSq-{name}"
        sig = p < 0.01
        all_tests.append((label, p, ef, sig))
        if sig:
            significant_tests.append(label)

    # --- Test 8: Pisano ---
    p = results_pisano['pisano_768_p']
    ef = abs(results_pisano['pisano_768_corr'])
    label = "T8-Pisano"
    sig = p < 0.01
    all_tests.append((label, p, ef, sig))
    if sig:
        significant_tests.append(label)

    # Print summary table
    print("\n  {:40s} {:>12s} {:>12s} {:>8s}".format("Test", "p-value", "Effect Size", "Sig?"))
    print("  " + "-"*76)
    for label, p, ef, sig in all_tests:
        marker = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        print(f"  {label:40s} {p:12.6f} {ef:12.6f} {marker:>8s}")

    # Bonferroni correction
    n_tests = len(all_tests)
    bonferroni_threshold = 0.05 / n_tests
    bonferroni_sig = [(label, p) for label, p, ef, sig in all_tests if p < bonferroni_threshold]

    print(f"\n  Total tests: {n_tests}")
    print(f"  Bonferroni threshold (0.05/{n_tests}): {bonferroni_threshold:.6f}")
    print(f"  Tests significant at nominal p<0.01: {len(significant_tests)}")
    print(f"  Tests significant after Bonferroni:  {len(bonferroni_sig)}")

    if bonferroni_sig:
        print(f"\n  Bonferroni-significant tests:")
        for label, p in bonferroni_sig:
            print(f"    {label}: p={p:.8f}")

    # Overall verdict
    print("\n" + "="*72)
    if len(bonferroni_sig) > 0:
        print("VERDICT: GOLDEN STRUCTURE SURVIVES (partially)")
        print(f"  {len(bonferroni_sig)} test(s) significant after multiple-comparison correction.")
        print("  The golden axiom (phi^2=phi+1) in SHA-256's round matrix skeleton")
        print("  leaves a detectable trace through 64 rounds of nonlinear operations.")
    elif len(significant_tests) > 0:
        print("VERDICT: MARGINAL EVIDENCE")
        print(f"  {len(significant_tests)} test(s) nominally significant (p<0.01),")
        print("  but none survive Bonferroni correction.")
        print("  Suggestive but not conclusive.")
    else:
        print("VERDICT: GOLDEN STRUCTURE IS FULLY DAMPED")
        print("  No test shows significance at p<0.01.")
        print("  SHA-256's nonlinear operations (Ch, Maj, carries) fully destroy")
        print("  the golden-ratio structure in the linearized round matrix.")
    print("="*72)


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("="*72)
    print("  GOLDEN RESONANCE TEST")
    print("  Does golden-ratio structure survive SHA-256?")
    print("="*72)
    print(f"  Samples per set: {NUM_SAMPLES}")
    print(f"  phi = {PHI:.15f}")
    print(f"  Golden rotation = {round(256 / PHI**2)} bits")
    print(f"  char_poly + 1 = x^4 * (x^2-x-1) * (x^2-x+1)")
    print(f"  x^2-x-1 = 0  =>  x = phi")
    print()

    t_start = time.time()

    # Generate all inputs
    set_a, set_b, set_c = generate_inputs()

    # TEST 1 is implicit in the generation
    print("\n" + "="*72)
    print("TEST 1: Input Generation Complete")
    print("="*72)
    print(f"  SET A: {len(set_a)} Fibonacci hashes (F_1 to F_{NUM_SAMPLES})")
    print(f"  SET B: {len(set_b)} phi-spaced hashes")
    print(f"  SET C: {len(set_c)} random hashes")
    # Quick sanity: are outputs roughly uniformly distributed?
    for name, outputs in [("A", set_a), ("B", set_b), ("C", set_c)]:
        mean_bits = np.mean([popcount(h) for h in outputs[:1000]])
        print(f"  SET {name} mean popcount (first 1000): {mean_bits:.2f} (expected: 128)")

    # Run all tests
    r_autocorr = test_autocorrelation(set_a, set_b, set_c)
    r_bits = test_bit_distribution(set_a, set_b, set_c)
    r_xor = test_xor_differential(set_a, set_b, set_c)
    r_rotation = test_golden_rotation(set_a, set_b, set_c)
    r_spectral = test_spectral(set_a, set_b, set_c)
    r_chi2 = test_chi_square(set_a, set_b, set_c)
    r_pisano = test_pisano(set_a)

    # Final analysis
    final_analysis(r_autocorr, r_bits, r_xor, r_rotation, r_spectral, r_chi2, r_pisano)

    print(f"\nTotal runtime: {time.time() - t_start:.2f}s")


if __name__ == "__main__":
    main()
