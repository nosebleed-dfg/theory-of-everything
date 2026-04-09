"""
SHA_PREDICTION_TEST — predicts golden leakage per round via gear model, then measures against real SHA-256
nos3bl33d
"""
import struct
import hashlib
import time
import math
import numpy as np

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
MOD32 = 2**32
MASK32 = MOD32 - 1
N_SAMPLES = 3000

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK32
def maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK32
def sigma0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def sigma1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def lsigma0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def lsigma1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def add32(*args): return sum(args) & MASK32

def make_golden_block():
    words = [int(i * PHI * MOD32) & MASK32 for i in range(16)]
    return struct.pack('>16I', *words)

def make_random_block():
    return bytes(np.random.randint(0, 256, 64, dtype=np.uint8))

def sha256_partial(block, num_rounds):
    W = list(struct.unpack('>16I', block))
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    a, b, c, d, e, f, g, h = H0
    for i in range(min(num_rounds, 64)):
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)
    return (a, b, c, d, e, f, g, h)

def count_leaky_channels(num_rounds, n_samples=N_SAMPLES):
    """Count channels where golden MI significantly exceeds random MI."""
    golden_bits = np.zeros((n_samples, 256), dtype=np.int8)
    random_bits = np.zeros((n_samples, 256), dtype=np.int8)

    for s in range(n_samples):
        gs = sha256_partial(make_golden_block(), num_rounds)
        rs = sha256_partial(make_random_block(), num_rounds)
        for bit in range(256):
            word = bit // 32
            pos = bit % 32
            golden_bits[s, bit] = (gs[word] >> pos) & 1
            random_bits[s, bit] = (rs[word] >> pos) & 1

    # For each bit: chi-square test for uniformity
    leaky = 0
    for bit in range(256):
        g_ones = np.sum(golden_bits[:, bit])
        r_ones = np.sum(random_bits[:, bit])
        g_zeros = n_samples - g_ones
        r_zeros = n_samples - r_ones

        # Chi-square: is golden distribution different from random?
        # Expected: n_samples/2 ones and zeros each
        g_chi = (g_ones - n_samples/2)**2 / (n_samples/2) + (g_zeros - n_samples/2)**2 / (n_samples/2)
        r_chi = (r_ones - n_samples/2)**2 / (n_samples/2) + (r_zeros - n_samples/2)**2 / (n_samples/2)

        # A channel "leaks" if golden has significantly more bias than random
        if g_chi > r_chi + 4.0:  # threshold
            leaky += 1

    return leaky

# GEAR MODEL PREDICTIONS
# Band structure: peaks and troughs by round
# Non-trough bands: 6-10(p=5), 25-36(F=12), 45-46(chi=2), 50-59(E/d=10)
# Kill zone: 60-64 (p=5)
# Shear per round: phi - 1/phi = 1
# Leak = shear * non-trough_fraction - kill_zone_damping

print("=" * 70)
print("GEAR MODEL PREDICTIONS vs MEASUREMENT")
print("=" * 70)
print()

test_rounds = [10, 20, 30, 36, 40, 50, 55, 59, 60, 64]

# Predict leakage for each round count
predictions = {}
for nr in test_rounds:
    # Count non-trough rounds up to nr
    non_trough = set()
    if nr >= 6: non_trough.update(range(6, min(nr+1, 11)))    # band 1
    if nr >= 19: non_trough.add(19)
    if nr >= 22: non_trough.add(22)
    if nr >= 25: non_trough.update(range(25, min(nr+1, 37)))   # band 2
    if nr >= 39: non_trough.add(39)
    if nr >= 45: non_trough.update(range(45, min(nr+1, 47)))   # band 3
    if nr >= 50: non_trough.update(range(50, min(nr+1, 60)))   # band 4

    n_peak = len(non_trough)
    n_trough = nr - n_peak

    # Kill zone rounds (60-64)
    kill = max(0, min(nr, 64) - 59) if nr >= 60 else 0

    # Net signal: peak_rounds - phi*kill_rounds (each kill round removes phi units)
    net_signal = n_peak - PHI * kill

    # Predicted leaky channels: proportional to net_signal / nr
    # Scale: at nr=64, we observed 19 leaky channels
    # At nr=64: n_peak=32, kill=5, net=32-8.09=23.91, ratio=23.91/64=0.374
    # Observed 19. So scale = 19/0.374 = 50.8

    if nr > 0:
        ratio = max(0, net_signal) / nr
    else:
        ratio = 0

    # Predict: leaky ~ ratio * 50.8 (calibrated to nr=64)
    predicted = ratio * 50.8
    predictions[nr] = (n_peak, n_trough, kill, net_signal, ratio, predicted)

print(f"{'Rounds':>6s} {'Peaks':>6s} {'Troughs':>7s} {'Kill':>5s} {'NetSig':>7s} {'Ratio':>6s} {'Pred':>6s}")
print("-" * 50)
for nr in test_rounds:
    p = predictions[nr]
    print(f"{nr:6d} {p[0]:6d} {p[1]:7d} {p[2]:5d} {p[3]:7.1f} {p[4]:6.3f} {p[5]:6.1f}")

print()
print("Now measuring (N=3000 per test, this takes a few minutes)...")
print()

# Measure
results = {}
t_start = time.time()
for nr in test_rounds:
    t0 = time.time()
    leaky = count_leaky_channels(nr)
    elapsed = time.time() - t0
    results[nr] = leaky
    pred = predictions[nr][5]
    err = abs(leaky - pred) / max(pred, 1) * 100 if pred > 0 else 0
    print(f"  Rounds={nr:2d}: predicted={pred:5.1f}, measured={leaky:3d}, "
          f"error={err:5.1f}%, time={elapsed:.1f}s")

total_time = time.time() - t_start
print(f"\nTotal measurement time: {total_time:.1f}s")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print(f"{'Rounds':>6s} {'Predicted':>10s} {'Measured':>10s} {'Error':>8s}")
print("-" * 40)
for nr in test_rounds:
    pred = predictions[nr][5]
    meas = results[nr]
    err = abs(meas - pred) / max(pred, 1) * 100
    match = "OK" if err < 50 else "MISS"
    print(f"{nr:6d} {pred:10.1f} {meas:10d} {err:7.1f}% {match}")
