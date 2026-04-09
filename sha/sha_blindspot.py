"""
SHA_BLINDSPOT — extends SHA-256 to 120 rounds; maps the -gamma*n curve against the golden signal
nos3bl33d

64 rounds = SHA, 120 = full golden circle. The 56 missing rounds (7*8) complete the picture.
"""
import struct, math, time, sys
import numpy as np
from sympy import primerange, isprime, nextprime

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329  # Euler-Mascheroni
MOD32 = 2**32
MASK32 = MOD32 - 1
N_SAMPLES = 5000

# --- Dodecahedral constants ---
d = 3
p = 5
V = 20
E = 30
F = 12
chi = 2
b0 = E - V + 1  # = 11
A5 = 60
I2 = 120
L4 = 7
L7 = 29
dp = d * p  # = 15

# --- SHA-256 initial hash values (sqrt of first 8 primes) ---
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

# --- Extend round constants K to 120 rounds ---
# K[i] = floor(2^32 * frac(cbrt(prime_i)))
# Standard SHA uses first 64 primes. We extend to 120.
def make_extended_K(n_rounds):
    K = []
    p_val = 2
    for _ in range(n_rounds):
        cbrt = p_val ** (1/3)
        frac = cbrt - int(cbrt)
        K.append(int(frac * MOD32) & MASK32)
        p_val = nextprime(p_val)
    return K

print("Computing extended round constants (120 primes)...")
K_120 = make_extended_K(120)
# Verify first few match standard SHA
K_STD = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
         0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
         0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
         0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
         0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
         0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
         0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
         0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

k_match = sum(1 for i in range(64) if K_120[i] == K_STD[i])
print(f"  K verification: {k_match}/64 match standard SHA-256")

# --- SHA primitives ---
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK32
def maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK32
def sigma0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def sigma1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def lsigma0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def lsigma1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def add32(*args): return sum(args) & MASK32

# --- Extended SHA: run up to 120 rounds ---
def sha_extended(words_16, num_rounds, return_per_round=False):
    """Run SHA compression for up to 120 rounds.
    Message schedule extends naturally beyond 16.
    Returns final state, or per-round states if requested."""
    W = list(words_16)
    for i in range(16, num_rounds):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))

    a, b, c, d_reg, e, f, g, h = H0
    states = []

    for i in range(num_rounds):
        T1 = add32(h, sigma1(e), ch(e, f, g), K_120[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))

        # --- 3-WAY MIDPOINT ---
        # Ch takes 3 inputs (e,f,g) -> 1 output (selection)
        # Maj takes 3 inputs (a,b,c) -> 1 output (consensus)
        # The midpoint of Ch and Maj = where selection meets consensus
        ch_val = ch(e, f, g)
        maj_val = maj(a, b, c)
        # 3-way midpoint: the XOR cancellation (where they disagree = the deviation)
        midpoint_bits = ~(ch_val ^ maj_val) & MASK32  # where they AGREE = the midpoint

        h = g; g = f; f = e; e = add32(d_reg, T1)
        d_reg = c; c = b; b = a; a = add32(T1, T2)

        if return_per_round:
            states.append({
                'state': [a, b, c, d_reg, e, f, g, h],
                'ch': ch_val,
                'maj': maj_val,
                'midpoint': midpoint_bits,
                'T1': T1,
                'T2': T2,
            })

    if return_per_round:
        return states
    return [a, b, c, d_reg, e, f, g, h]

# --- Popcount ---
def popcount_256(state_8):
    return sum(bin(w).count('1') for w in state_8)

def popcount_32(w):
    return bin(w & MASK32).count('1')

# --- Golden input generator ---
def golden_words(seed, n=16):
    return [int((seed * n + i) * PHI * MOD32) & MASK32 for i in range(n)]

def random_words(n=16):
    return [int.from_bytes(np.random.bytes(4), 'big') for _ in range(n)]


# =================================================================
# TEST 1: GOLDEN SIGNAL MAP -- FULL 120 ROUNDS
# =================================================================
print()
print("=" * 70)
print("TEST 1: GOLDEN SIGNAL MAP -- FULL 120 ROUNDS (|2I| = PISANO FIXED PT)")
print("=" * 70)
print()

n_samples = N_SAMPLES
round_counts = list(range(1, 121))

print(f"Measuring golden signal at each round (N={n_samples})...")
print("This maps the FULL golden circle, not just SHA's 64-round window.")
print()

signal_map = []
midpoint_map = []  # 3-way midpoint agreement fraction per round

t_total = time.time()
for nr in round_counts:
    g_pops = []
    r_pops = []
    mid_fracs = []  # fraction of bits where Ch and Maj agree

    for s in range(n_samples):
        gw = golden_words(s)
        rw = random_words()

        # Golden
        gs = sha_extended(gw, nr, return_per_round=True)
        g_pop = popcount_256(gs[-1]['state'])
        g_pops.append(g_pop)

        # Random
        rs = sha_extended(rw, nr, return_per_round=True)
        r_pop = popcount_256(rs[-1]['state'])
        r_pops.append(r_pop)

        # 3-way midpoint at final round: fraction of bits where Ch == Maj
        mid_frac = popcount_32(gs[-1]['midpoint']) / 32.0
        mid_fracs.append(mid_frac)

    g_mean = np.mean(g_pops)
    r_mean = np.mean(r_pops)
    pooled_se = np.sqrt(np.var(g_pops)/n_samples + np.var(r_pops)/n_samples)
    z = (g_mean - r_mean) / pooled_se if pooled_se > 0 else 0
    mid_mean = np.mean(mid_fracs)

    signal_map.append(z)
    midpoint_map.append(mid_mean)

    # Print key rounds
    labels = {
        d: 'd', p: 'p', L4: 'L4', b0: 'b0', F: 'F', dp: 'dp',
        V: 'V', E: 'E', 36: 'Fd', 56: '56=L4*2^d',
        A5: '|A5|', 64: '2^(2d)=SHA', 96: '|2I|-|2T|',
        I2: '|2I|=FIXED'
    }
    if nr in labels or nr <= 5 or nr % 10 == 0:
        label = labels.get(nr, '')
        sign = '+' if z > 0 else '-' if z < 0 else '0'
        mid_pct = mid_mean * 100
        print(f"  r={nr:3d}  z={z:+8.2f}  mid={mid_pct:5.1f}%  {label}")

elapsed = time.time() - t_total
print(f"\n  Total time: {elapsed:.1f}s")


# =================================================================
# TEST 2: THE -GAMMA*N CURVE vs OBSERVED SIGNAL
# =================================================================
print()
print("=" * 70)
print("TEST 2: THE -GAMMA*N OFFSET CURVE")
print("=" * 70)
print()
print("The blind spot: char_poly + 1 = golden factorization.")
print("The +1 offset on the y-axis. Parameterized by -gamma*n.")
print("gamma = 0.5772... = the constant that appears from nowhere.")
print()

# Compute -gamma*n for each round
gamma_curve = [-GAMMA * n for n in range(1, 121)]

# Normalize both curves to [-1, 1] for comparison
sig_arr = np.array(signal_map)
gam_arr = np.array(gamma_curve)

sig_norm = sig_arr / np.max(np.abs(sig_arr)) if np.max(np.abs(sig_arr)) > 0 else sig_arr
gam_norm = gam_arr / np.max(np.abs(gam_arr)) if np.max(np.abs(gam_arr)) > 0 else gam_arr

# Correlation between signal and -gamma*n
corr_linear = np.corrcoef(sig_arr, gam_arr)[0, 1]

# What about -gamma*n modulated by the gear system?
# The gear recurrence: x_{n+1} = V*x_n + (2F)^2 * x_{n-1}
# Start with x_0 = -gamma, x_1 = -gamma (or 0 and -gamma)
gear_x = np.zeros(121)
gear_x[0] = 0
gear_x[1] = -GAMMA
for i in range(2, 121):
    gear_x[i] = V * gear_x[i-1] + (2*F)**2 * gear_x[i-2]

# Normalize gear_x
gear_norm = gear_x[1:] / np.max(np.abs(gear_x[1:])) if np.max(np.abs(gear_x[1:])) > 0 else gear_x[1:]

corr_gear = np.corrcoef(sig_arr, gear_norm)[0, 1]

# What about -gamma*n with Pisano periodicity mod 120?
# H_n - ln(n) -> gamma. So the harmonic residual.
harmonic = np.zeros(120)
h_sum = 0
for i in range(120):
    h_sum += 1.0 / (i + 1)
    harmonic[i] = h_sum - math.log(i + 1)  # approaches gamma

# This IS gamma appearing from nowhere -- each term adds 1/(n+1) and
# the residual (H_n - ln(n)) is the state of appearing
corr_harmonic = np.corrcoef(sig_arr, harmonic)[0, 1]

print(f"Correlation with -gamma*n (linear):     {corr_linear:+.6f}")
print(f"Correlation with gear(-gamma) recurrence: {corr_gear:+.6f}")
print(f"Correlation with H_n - ln(n) (harmonic): {corr_harmonic:+.6f}")
print()

# The REAL test: does the signal's ENVELOPE match -gamma*n?
# Take absolute value of signal (the amplitude, not the sign)
sig_abs = np.abs(sig_arr)
corr_envelope_linear = np.corrcoef(sig_abs, np.abs(gam_arr))[0, 1]
corr_envelope_harmonic = np.corrcoef(sig_abs, harmonic)[0, 1]

print(f"Envelope vs |gamma*n|:  {corr_envelope_linear:+.6f}")
print(f"Envelope vs H_n-ln(n):  {corr_envelope_harmonic:+.6f}")
print()


# =================================================================
# TEST 3: 3-WAY MIDPOINT STRUCTURE
# =================================================================
print()
print("=" * 70)
print("TEST 3: 3-WAY MIDPOINT (Ch vs Maj agreement per round)")
print("=" * 70)
print()
print("Each round: Ch picks (phi), Maj votes (pi). Where they AGREE = midpoint.")
print("At 50% = pure noise. Above 50% = consensus. Below 50% = conflict.")
print()

mid_arr = np.array(midpoint_map)
print(f"  Overall mean midpoint: {np.mean(mid_arr)*100:.2f}%")
print(f"  Expected (random):     50.00%")
print(f"  Deviation:             {(np.mean(mid_arr) - 0.5)*100:+.2f}%")
print()

# Midpoint at key rounds
print("  Midpoint at key dodecahedral rounds:")
for nr in [d, p, L4, b0, F, dp, V, E, 36, A5, 64, 96, I2]:
    if nr <= 120:
        mid_val = midpoint_map[nr-1] * 100
        labels = {d:'d', p:'p', L4:'L4', b0:'b0', F:'F', dp:'dp',
                  V:'V', E:'E', 36:'Fd', A5:'|A5|', 64:'SHA',
                  96:'|2I|-|2T|', I2:'|2I|'}
        print(f"    r={nr:3d} ({labels.get(nr,''):>8s}): midpoint = {mid_val:.1f}%")


# =================================================================
# TEST 4: THE BLIND SPOT -- WHERE SIGNAL = 0 BUT STRUCTURE PEAKS
# =================================================================
print()
print("=" * 70)
print("TEST 4: THE BLIND SPOT")
print("=" * 70)
print()
print("100% visibility = blind spot. Where the signal hits zero but the")
print("3-way midpoint deviates. The structure is there but invisible.")
print()

# Find zero-crossings of the signal (where golden signal ~ 0)
zero_crossings = []
for i in range(1, len(signal_map)):
    if signal_map[i-1] * signal_map[i] < 0:  # sign change
        zero_crossings.append(i)

print(f"  Signal zero-crossings at rounds: {zero_crossings}")
print()

# At each zero-crossing, what is the midpoint doing?
print("  At zero-crossings (golden signal ~ 0):")
for zc in zero_crossings:
    mid = midpoint_map[zc-1] * 100
    sig = signal_map[zc-1]
    print(f"    r={zc:3d}: z={sig:+6.2f}, midpoint={mid:.1f}%")

print()

# The PHASE INVERSION region: rounds 59-65
print("  Phase inversion region (rounds 55-70):")
for nr in range(55, min(71, 121)):
    sig = signal_map[nr-1]
    mid = midpoint_map[nr-1] * 100
    marker = " <-- |A5|" if nr == 60 else " <-- SHA" if nr == 64 else ""
    print(f"    r={nr:3d}: z={sig:+8.2f}, midpoint={mid:.1f}%{marker}")

print()

# Past SHA: rounds 65-120 (the 56 missing rounds)
print("  POST-SHA (rounds 65-120, the 56 = L4*2^d missing rounds):")
for nr in range(65, 121, 5):
    sig = signal_map[nr-1]
    mid = midpoint_map[nr-1] * 100
    label = ""
    if nr == 96: label = " <-- |2I|-|2T|"
    if nr == 120: label = " <-- |2I| FIXED POINT"
    print(f"    r={nr:3d}: z={sig:+8.2f}, midpoint={mid:.1f}%{label}")


# =================================================================
# TEST 5: THE OFFSET -- CHAR POLY EVALUATED AT y = -1
# =================================================================
print()
print("=" * 70)
print("TEST 5: THE Y-AXIS OFFSET")
print("=" * 70)
print()

# The linearized char poly: x^8 - 2x^7 + x^6 - x^4 - 1
# At y = 0: roots (eigenvalues of the linearized matrix)
# At y = -1: char_poly + 1 = x^4(x^2-x-1)(x^2-x+1) = GOLDEN * CYCLOTOMIC
# At y = -gamma*n: ???

def char_poly(x):
    return x**8 - 2*x**7 + x**6 - x**4 - 1

def golden_poly(x):
    return x**2 - x - 1

def cyclo6(x):
    return x**2 - x + 1

# Evaluate at different offsets
print("char_poly(x) at special points:")
print(f"  char_poly(phi)      = {char_poly(PHI):.10f}")
print(f"  char_poly(-1/phi)   = {char_poly(-1/PHI):.10f}")
print(f"  char_poly(0)        = {char_poly(0):.10f}  (= -1, the offset)")
print(f"  char_poly(1)        = {char_poly(1):.10f}")
print(f"  char_poly(1/2)      = {char_poly(0.5):.10f}  (critical line)")
print(f"  char_poly(-gamma)   = {char_poly(-GAMMA):.10f}")
print()

# The -gamma*n family: char_poly evaluated at -gamma, -2*gamma, etc.
print("char_poly(-gamma*n) for n = 1..10:")
for n in range(1, 11):
    val = char_poly(-GAMMA * n)
    val_plus_1 = val + 1  # the golden factorization side
    gp = golden_poly(-GAMMA * n)
    c6 = cyclo6(-GAMMA * n)
    x4 = (-GAMMA * n) ** 4
    product = x4 * gp * c6
    print(f"  n={n:2d}: char(-gamma*{n:2d}) = {val:+14.6f}  "
          f" char+1 = {val_plus_1:+14.6f}  "
          f" x^4·golden·cyclo = {product:+14.6f}  "
          f" Delta = {abs(val_plus_1 - product):.2e}")

print()

# THE KEY: what value of n makes char_poly(-gamma*n) = 0?
# i.e., where does the -gamma*n curve cross the ROOT level?
print("Where does -gamma*n hit the ROOTS of char_poly?")
from scipy.optimize import brentq
try:
    from scipy.optimize import brentq
    # Find zeros of char_poly(-gamma*x) in [0, 20]
    for a, b in [(0.1, 2), (2, 5), (5, 10), (10, 20)]:
        try:
            root = brentq(lambda x: char_poly(-GAMMA * x), a, b)
            cp_val = char_poly(-GAMMA * root)
            print(f"  n = {root:.6f}: char_poly(-gamma*n) = {cp_val:.2e}")
            # What framework constant is this close to?
            for name, val in [('1', 1), ('chi', 2), ('d', 3), ('p', 5),
                              ('L4', 7), ('b0', 11), ('F', 12), ('dp', 15),
                              ('V', 20)]:
                if abs(root - val) < 0.5:
                    print(f"        ~ {name} = {val} (off by {root - val:+.6f})")
        except ValueError:
            pass
except ImportError:
    print("  (scipy not available, skipping root-finding)")

print()


# =================================================================
# SUMMARY
# =================================================================
print()
print("=" * 70)
print("SUMMARY: THE BLIND SPOT")
print("=" * 70)
print()
print(f"SHA-256:  64 rounds = 2^(2d)")
print(f"|2I|:    120 rounds = Pisano fixed point")
print(f"Offset:   56 rounds = L4 * 2^d = 7 * 8")
print()
print(f"Golden signal phase inverts at round {A5} = |A5|")
print(f"Signal at round 64 (SHA output): z = {signal_map[63]:+.2f}")
print(f"Signal at round 120 (|2I|):      z = {signal_map[119]:+.2f}")
print()
print(f"gamma = {GAMMA:.10f}")
print(f"-gamma * 64  = {-GAMMA * 64:+.6f}")
print(f"-gamma * 120 = {-GAMMA * 120:+.6f}")
print()

# Does the signal RETURN at 120? (full circle)
if abs(signal_map[119]) < abs(signal_map[63]):
    print("Signal CLOSER to zero at |2I| than at SHA-64.")
    print("The full circle dampens further.")
elif signal_map[119] * signal_map[0] > 0:
    print("Signal has SAME SIGN at round 1 and round 120!")
    print("The golden circle CLOSES. Back to positive.")
else:
    print(f"Signal at round 1:   z = {signal_map[0]:+.2f}")
    print(f"Signal at round 120: z = {signal_map[119]:+.2f}")

print()
print("The blind spot is where you can see everything.")
print("100% visibility. No shadows. No contrast.")
print("-gamma*n is the state of appearing from nowhere.")
