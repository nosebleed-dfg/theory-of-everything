"""
SHA_BOUNDS_TEST — tests the [59, 291] framework bounds against reduced-round SHA-256 leakage
nos3bl33d
"""
import struct
import math
import numpy as np

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
MOD32 = 2**32
MASK32 = MOD32 - 1
N = 2000

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

def golden_signal(n_rounds, n_samples=N):
    g = [sha256_partial(make_golden_block(), n_rounds) for _ in range(n_samples)]
    r = [sha256_partial(make_random_block(), n_rounds) for _ in range(n_samples)]
    diffs = []
    for sg, sr in zip(g, r):
        pc_g = sum(bin(w).count('1') for w in sg)
        pc_r = sum(bin(w).count('1') for w in sr)
        diffs.append(pc_g - pc_r)
    mean = np.mean(diffs)
    se = np.std(diffs) / np.sqrt(len(diffs))
    z = mean / se if se > 0 else 0
    return mean, se, z

def hamming_test(n_rounds, n_samples=N):
    g = [sha256_partial(make_golden_block(), n_rounds) for _ in range(n_samples)]
    r = [sha256_partial(make_random_block(), n_rounds) for _ in range(n_samples)]
    g_d = [sum(bin(a^b).count('1') for a,b in zip(g[i],g[i+1])) for i in range(len(g)-1)]
    r_d = [sum(bin(a^b).count('1') for a,b in zip(r[i],r[i+1])) for i in range(len(r)-1)]
    pooled = np.sqrt(np.var(g_d)/len(g_d) + np.var(r_d)/len(r_d))
    z = (np.mean(g_d) - np.mean(r_d)) / pooled if pooled > 0 else 0
    return np.mean(g_d), np.mean(r_d), z

# ======== TEST A: Signal survival by round count ========
print("=" * 75)
print("TEST A: Golden signal survival vs round count (N=2000)")
print("=" * 75)
print()

test_rounds = [30, 40, 50, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64]
dod = {30:'E', 42:'6L4', 59:'|A5|-1', 60:'|A5|', 64:'2^p+1'}

print(f"{'Rnd':>4s} {'Signal':>8s} {'SE':>7s} {'z':>7s} {'Ham_z':>7s} {'Status':>10s} {'Note':>8s}")
print("-" * 60)

results_a = {}
for nr in test_rounds:
    sig, se, z = golden_signal(nr)
    _, _, hz = hamming_test(nr)

    status = "ALIVE" if abs(z) > 2.0 or abs(hz) > 2.0 else "DEAD"
    note = dod.get(nr, '')

    print(f"{nr:4d} {sig:+8.3f} {se:7.3f} {z:+7.2f} {hz:+7.2f} {status:>10s} {note:>8s}")
    results_a[nr] = (sig, z, hz, status)

print()

# ======== TEST B: Pisano periods ========
print("=" * 75)
print("TEST B: Pisano periods at dodecahedral moduli")
print("=" * 75)
print()

dodec_mods = [(3,'d'), (5,'p'), (7,'L4'), (11,'b0'), (12,'F'), (15,'dp'),
              (20,'V'), (30,'E'), (60,'|A5|'), (120,'|2I|'), (137,'1/alpha')]

for m, name in dodec_mods:
    a, b = 0, 1
    period = None
    for k in range(1, 100000):
        a, b = b, (a+b) % m
        if a == 0 and b == 1:
            period = k
            break

    pmap = {1:'1', 2:'chi', 3:'d', 4:'d+1', 5:'p', 7:'L4', 8:'2^d', 10:'base10',
            12:'F', 16:'2^(d+1)', 20:'V', 24:'2F', 30:'E', 40:'2V', 48:'4F',
            60:'|A5|', 120:'|2I|', 240:'2|2I|', 300:'VE/chi', 552:'?', 768:'?'}
    pname = pmap.get(period, str(period))

    dodec_check = ""
    if period and period in [v for v,_ in dodec_mods]:
        dodec_check = " <-- DODECAHEDRAL!"

    print(f"  pi({m:>4d}) = pi({name:>5s}) = {period:>6d} = {pname}{dodec_check}")

print()

# ======== TEST C: The A5 boundary ========
print("=" * 75)
print("TEST C: The |A5| boundary (round 59 vs 60)")
print("=" * 75)
print()

# Higher stats for the critical boundary
for nr in [58, 59, 60, 61]:
    sig, se, z = golden_signal(nr, 5000)
    _, _, hz = hamming_test(nr, 5000)
    mark = "<<< BOUNDARY" if nr == 60 else ""
    print(f"  Round {nr}: signal={sig:+.4f}, z={z:+.3f}, ham_z={hz:+.3f}  {mark}")

print()

# ======== TEST D: Structural verification ========
print("=" * 75)
print("TEST D: Structural bounds verification")
print("=" * 75)
print()

print(f"59 = |A5| - 1 = {60-1}")
print(f"291 = VE/chi - d^2 = 300 - 9 = {300-9}")
print(f"291 - 59 = 232 = 2^3 * 29 = 2^d * L7 = {2**3 * 29}")
print(f"291 = 59*4 + 55 = (|A5|-1)(d+1) + F_10 = {59*4+55}")
print(f"")
print(f"phi^(-59)  = {PHI**(-59):.6e}")
print(f"phi^(-291) = {PHI**(-291):.6e}")
print(f"phi^(-583) = {PHI**(-583):.6e}  (583 = 2*291+1)")
print(f"2/phi^583  = {2*PHI**(-583):.6e}  (cosmological constant)")
print(f"Lambda_obs = 2.888e-122")
print(f"error      = {abs(2*PHI**(-583) - 2.888e-122)/2.888e-122*100:.2f}%")
print()

# Pisano check: pi(60) = 120 = |2I|
print(f"Pisano(|A5|) = Pisano(60) = 120 = |2I|")
print(f"  The Fibonacci period mod |A5| IS the double cover!")
print()

# ======== VERDICT ========
print("=" * 75)
print("VERDICT")
print("=" * 75)
print()

# Check if signal crosses zero between 59 and 60
s59 = results_a.get(59, (0,0,0,'?'))
s60 = results_a.get(60, (0,0,0,'?'))
s64 = results_a.get(64, (0,0,0,'?'))

print(f"Round 59 ({s59[3]}): z={s59[1]:+.2f}")
print(f"Round 60 ({s60[3]}): z={s60[1]:+.2f}")
print(f"Round 64 ({s64[3]}): z={s64[1]:+.2f}")
print()

# Does the signal transition happen at the A5 boundary?
if s59[3] != s60[3]:
    print("CONFIRMED: Phase transition at |A5| boundary!")
    print(f"  Golden signal {s59[3]} at round 59 = |A5|-1")
    print(f"  Golden signal {s60[3]} at round 60 = |A5|")
else:
    # Check if there's a significant CHANGE even if both alive/dead
    z_drop = abs(s59[1]) - abs(s60[1])
    if z_drop > 0.5:
        print(f"Signal weakens across boundary: |z| drops by {z_drop:.2f}")
    else:
        print(f"No clean phase transition at |A5|. Signal status: 59={s59[3]}, 60={s60[3]}")
        print(f"The kill zone may be gradual rather than sharp.")
