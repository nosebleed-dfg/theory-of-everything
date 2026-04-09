"""
SHA_LEAKAGE_PROPER — proper MI metric for reduced-round SHA-256; input bit 0 vs each of 256 output bits
nos3bl33d
"""
import struct
import math
import numpy as np
import time

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
MOD32 = 2**32
MASK32 = MOD32 - 1
N = 5000

H0_INIT = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
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

def sha_partial(block_words, num_rounds):
    W = list(block_words)
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    a, b, c, d, e, f, g, h = H0_INIT
    for i in range(min(num_rounds, 64)):
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)
    return [a, b, c, d, e, f, g, h]

def mutual_info_binary(x, y):
    n = len(x)
    p00 = np.sum((x==0)&(y==0)) / n
    p01 = np.sum((x==0)&(y==1)) / n
    p10 = np.sum((x==1)&(y==0)) / n
    p11 = np.sum((x==1)&(y==1)) / n
    px0, px1 = p00+p01, p10+p11
    py0, py1 = p00+p10, p01+p11
    mi = 0.0
    for pxy, px, py in [(p00,px0,py0),(p01,px0,py1),(p10,px1,py0),(p11,px1,py1)]:
        if pxy > 0 and px > 0 and py > 0:
            mi += pxy * math.log2(pxy / (px * py))
    return mi

def measure_leakage(num_rounds, n_samples=N):
    """Measure MI between input structure and output bits.
    Golden inputs: vary bit 0 of word 0, keep golden structure elsewhere.
    Random inputs: vary bit 0 of word 0, random elsewhere.
    Leaky channel = golden MI significantly exceeds random MI."""

    golden_input_bits = np.zeros(n_samples, dtype=np.int8)
    random_input_bits = np.zeros(n_samples, dtype=np.int8)
    golden_output = np.zeros((n_samples, 256), dtype=np.int8)
    random_output = np.zeros((n_samples, 256), dtype=np.int8)

    for s in range(n_samples):
        # Golden: golden-angle spaced words, flip bit 0 randomly
        words_g = [int(i * PHI * MOD32) & MASK32 for i in range(16)]
        flip = np.random.randint(2)
        golden_input_bits[s] = flip
        words_g[0] ^= flip  # flip bit 0

        # Random: random words, flip bit 0 randomly
        words_r = [int(np.random.randint(0, 2**31)) | (np.random.randint(0, 2) << 31) for _ in range(16)]
        random_input_bits[s] = flip  # same flip pattern
        words_r[0] ^= flip

        gs = sha_partial(words_g, num_rounds)
        rs = sha_partial(words_r, num_rounds)

        for bit in range(256):
            w, pos = bit // 32, bit % 32
            golden_output[s, bit] = (gs[w] >> pos) & 1
            random_output[s, bit] = (rs[w] >> pos) & 1

    # For each output bit: MI with input bit
    leaky = 0
    golden_mis = []
    random_mis = []

    for bit in range(256):
        g_mi = mutual_info_binary(golden_input_bits, golden_output[:, bit])
        r_mi = mutual_info_binary(random_input_bits, random_output[:, bit])
        golden_mis.append(g_mi)
        random_mis.append(r_mi)

        if g_mi > r_mi + 0.001:  # golden has meaningfully more MI
            leaky += 1

    avg_g = np.mean(golden_mis)
    avg_r = np.mean(random_mis)
    max_g = np.max(golden_mis)
    max_r = np.max(random_mis)

    return leaky, avg_g, avg_r, max_g, max_r

# Test at key round counts
test_rounds = [10, 20, 30, 36, 46, 50, 55, 59, 60, 62, 64]

print("SHA REDUCED-ROUND LEAKAGE (MI metric, N=5000)")
print("=" * 75)
print()
print(f"{'Rnd':>4s} {'Leaky':>6s} {'G_avg_MI':>10s} {'R_avg_MI':>10s} {'G_max':>10s} {'R_max':>10s} {'Time':>6s}")
print("-" * 65)

results = {}
for nr in test_rounds:
    t0 = time.time()
    leaky, g_avg, r_avg, g_max, r_max = measure_leakage(nr)
    elapsed = time.time() - t0
    results[nr] = (leaky, g_avg, r_avg)
    print(f"{nr:4d} {leaky:6d} {g_avg:10.6f} {r_avg:10.6f} {g_max:10.6f} {r_max:10.6f} {elapsed:5.1f}s")

print()
print("ANALYSIS")
print("=" * 75)
print()

# The prediction: leakage peaks before kill zone (round 59)
# and drops after kill zone (rounds 60-64)
for nr in test_rounds:
    leaky = results[nr][0]
    bar = "#" * min(leaky, 80)
    dodec = {10:'', 20:'V', 30:'E', 36:'', 46:'', 50:'', 55:'', 59:'|A5|-1', 60:'|A5|', 62:'', 64:'2^6'}
    label = dodec.get(nr, '')
    print(f"  {nr:2d} rounds: {leaky:3d} leaky  {bar}  {label}")

print()

# Compare golden vs random average MI
print("Golden/Random MI ratio by round:")
for nr in test_rounds:
    _, g, r = results[nr]
    ratio = g/r if r > 0 else 0
    print(f"  {nr:2d}: golden={g:.6f} random={r:.6f} ratio={ratio:.3f}")
