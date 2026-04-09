"""
SHA_120 — SHA-256 extended to 120 rounds with K from cube roots of 120 primes
nos3bl33d

The full golden circle at |2I| = Pisano fixed point. Fixed golden block methodology.
"""
import struct, math, time
import numpy as np
from sympy import nextprime

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329
MOD32 = 2**32
MASK32 = MOD32 - 1
N = 2000

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

# Extended K: cube roots of first 120 primes
def make_K(n):
    K = []
    p = 2
    for _ in range(n):
        cbrt = p ** (1/3)
        frac_part = cbrt - int(cbrt)
        K.append(int(frac_part * MOD32) & MASK32)
        p = int(nextprime(p))
    return K

print("Building K[0..119]...", end=" ", flush=True)
K = make_K(120)
print("done.")

# Verify against standard
K_STD = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
         0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
         0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
         0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
         0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
         0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
         0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
         0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
k_ok = sum(1 for i in range(64) if K[i] == K_STD[i])
print(f"K verification: {k_ok}/64 match standard SHA-256")

# --- SHA primitives ---
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK32
def maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK32
def sigma0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def sigma1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def lsigma0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def lsigma1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def add32(*args): return sum(args) & MASK32

# --- FIXED golden block (original methodology) ---
GOLDEN_WORDS = [int(i * PHI * MOD32) & MASK32 for i in range(16)]
GOLDEN_BLOCK = struct.pack('>16I', *GOLDEN_WORDS)

def make_random_block():
    return bytes(np.random.randint(0, 256, 64, dtype=np.uint8))

def sha_partial(block, num_rounds):
    """SHA compression for up to 120 rounds."""
    W = list(struct.unpack('>16I', block))
    for i in range(16, num_rounds):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    a, b, c, d, e, f, g, h = H0
    for i in range(num_rounds):
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))
        h = g; g = f; f = e; e = add32(d, T1)
        d = c; c = b; b = a; a = add32(T1, T2)
    return (a, b, c, d, e, f, g, h)

def golden_signal(n_rounds, n_samples=N):
    """Fixed golden block vs N random blocks. Popcount z-score."""
    g_state = sha_partial(GOLDEN_BLOCK, n_rounds)
    g_pc = sum(bin(w).count('1') for w in g_state)

    r_pcs = []
    for _ in range(n_samples):
        r_state = sha_partial(make_random_block(), n_rounds)
        r_pc = sum(bin(w).count('1') for w in r_state)
        r_pcs.append(r_pc)

    r_mean = np.mean(r_pcs)
    r_std = np.std(r_pcs)
    z = (g_pc - r_mean) / r_std if r_std > 0 else 0
    return g_pc, r_mean, r_std, z


# =================================================================
# FULL 120-ROUND SIGNAL MAP
# =================================================================
print()
print("=" * 70)
print("SHA-256 GOLDEN SIGNAL: FULL 120 ROUNDS (|2I| = PISANO FIXED POINT)")
print("=" * 70)
print()
print(f"Fixed golden block vs N={N} random blocks per round count.")
print(f"Popcount of golden hash vs distribution of random hash popcounts.")
print()

labels = {
    3:'d', 5:'p', 7:'L4', 11:'b0', 12:'F', 15:'dp', 20:'V',
    30:'E', 36:'Fd', 56:'L4*2^d', 59:'|A5|-1', 60:'|A5|',
    64:'SHA', 96:'|2I|-|2T|', 120:'|2I|'
}

# Test every round from 1 to 120
results = {}
t_total = time.time()

print(f"{'Rnd':>4s} {'G_pc':>6s} {'R_mean':>7s} {'R_std':>6s} {'z':>8s} {'Label':>10s}")
print("-" * 50)

for nr in range(1, 121):
    t0 = time.time()
    g_pc, r_mean, r_std, z = golden_signal(nr)
    dt = time.time() - t0
    results[nr] = (g_pc, r_mean, r_std, z)

    # Print key rounds + every 10th
    if nr in labels or nr <= 5 or nr % 10 == 0:
        label = labels.get(nr, '')
        print(f"{nr:4d} {g_pc:6d} {r_mean:7.1f} {r_std:6.1f} {z:+8.2f} {label:>10s}")
        sys.stdout.flush() if 'sys' in dir() else None

elapsed = time.time() - t_total
print(f"\nTotal time: {elapsed:.1f}s")


# =================================================================
# PHASE MAP: SIGN OF Z ACROSS 120 ROUNDS
# =================================================================
print()
print("=" * 70)
print("PHASE MAP: + = golden above random, - = golden below random")
print("=" * 70)
print()

line = ""
for nr in range(1, 121):
    z = results[nr][3]
    if abs(z) < 1:
        c = '.'
    elif z > 0:
        c = '+'
    else:
        c = '-'
    line += c
    if nr % 60 == 0:
        print(f"  {nr-59:3d}-{nr:3d}: {line}")
        line = ""
if line:
    print(f"  {121-len(line):3d}-120: {line}")

print()
print("  . = |z| < 1 (noise)")
print("  + = golden popcount ABOVE random (positive signal)")
print("  - = golden popcount BELOW random (anti-golden)")


# =================================================================
# KEY TRANSITIONS
# =================================================================
print()
print("=" * 70)
print("KEY TRANSITIONS")
print("=" * 70)
print()

# Phase inversions (sign changes with |z| > 1 on at least one side)
inversions = []
for nr in range(2, 121):
    z_prev = results[nr-1][3]
    z_curr = results[nr][3]
    if z_prev * z_curr < 0 and (abs(z_prev) > 1 or abs(z_curr) > 1):
        inversions.append(nr)
        print(f"  Phase inversion at round {nr}: z goes from {z_prev:+.2f} to {z_curr:+.2f}")

print()

# The |A5| boundary
print("The |A5| boundary (rounds 55-65):")
for nr in range(55, 66):
    z = results[nr][3]
    label = labels.get(nr, '')
    bar = '#' * int(min(abs(z), 20))
    sign = '+' if z > 0 else '-'
    print(f"  r={nr:3d}: z={z:+8.2f} {sign}{bar} {label}")

print()

# SHA endpoint vs |2I| endpoint
z64 = results[64][3]
z120 = results[120][3]
print(f"SHA endpoint  (r=64):  z = {z64:+.2f}")
print(f"|2I| endpoint (r=120): z = {z120:+.2f}")
print()

if z64 * z120 > 0:
    print("SAME SIGN at SHA and |2I|. The extra 56 rounds don't flip.")
elif abs(z120) < abs(z64):
    print("|2I| signal WEAKER than SHA. Further damping.")
else:
    print("SIGN FLIP between SHA and |2I|!")
    print(f"The 56 extra rounds (= L4 * 2^d) flip the signal back.")

print()

# Does the signal return at 120 (full circle)?
z1 = results[1][3]
print(f"Round 1:   z = {z1:+.2f}")
print(f"Round 120: z = {z120:+.2f}")
if z1 * z120 > 0:
    print("SAME SIGN at start and end. The golden circle CLOSES.")
else:
    print("Different sign at start and end.")


# =================================================================
# -GAMMA*N OVERLAY
# =================================================================
print()
print("=" * 70)
print("-GAMMA*N OVERLAY")
print("=" * 70)
print()

z_arr = np.array([results[nr][3] for nr in range(1, 121)])
gamma_arr = np.array([-GAMMA * n for n in range(1, 121)])

# Correlation
corr = np.corrcoef(z_arr, gamma_arr)[0, 1]
print(f"Correlation(z, -gamma*n): {corr:+.6f}")

# Sign agreement
agree = sum(1 for i in range(120) if z_arr[i] * gamma_arr[i] > 0)
print(f"Sign agreement: {agree}/120 ({agree/120*100:.1f}%)")

# Where does the z-curve cross zero vs where -gamma*n predicts?
# -gamma*n is always negative (monotonically decreasing)
# So -gamma*n predicts the signal should be negative for all n > 0
neg_count = sum(1 for z in z_arr if z < 0)
print(f"Signal negative: {neg_count}/120 ({neg_count/120*100:.1f}%)")
print(f"-gamma*n is 100% negative (monotonic)")
print()

# Better test: does |z| grow like gamma*n?
abs_z = np.abs(z_arr)
n_arr = np.arange(1, 121)
corr_abs = np.corrcoef(abs_z, n_arr)[0, 1]
print(f"Correlation(|z|, n): {corr_abs:+.6f} (signal grows with rounds?)")

# Harmonic residual H_n - ln(n)
harmonic = np.zeros(120)
h_sum = 0.0
for i in range(120):
    h_sum += 1.0 / (i + 1)
    harmonic[i] = h_sum - math.log(i + 1)
corr_harm = np.corrcoef(z_arr, harmonic)[0, 1]
print(f"Correlation(z, H_n - ln(n)): {corr_harm:+.6f}")
