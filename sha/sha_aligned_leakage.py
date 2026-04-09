"""
SHA_ALIGNED_LEAKAGE — measures aligned mutual information between SHA-256 output and the golden template
nos3bl33d
"""
import struct, math, time
import numpy as np

np.random.seed(42)
PHI = (1 + 5**0.5) / 2
MOD32 = 2**32
MASK32 = MOD32 - 1
N = 10000

H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ch(e,f,g): return ((e&f)^(~e&g))&MASK32
def maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK32
def s0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def s1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ls0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ls1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def a32(*args): return sum(args)&MASK32

def sha_partial(words, nr):
    W = list(words)
    for i in range(16,64):
        W.append(a32(ls1(W[i-2]),W[i-7],ls0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h = H0
    for i in range(min(nr,64)):
        T1=a32(h,s1(e),ch(e,f,g),K[i],W[i])
        T2=a32(s0(a),maj(a,b,c))
        h=g;g=f;f=e;e=a32(d,T1);d=c;c=b;b=a;a=a32(T1,T2)
    return [a,b,c,d,e,f,g,h]

def golden_template_256():
    """The golden template: what a perfectly golden output would look like."""
    t = []
    for bit in range(256):
        # Golden angle spacing at each bit position
        val = int(bit * PHI * 137) % 2  # golden-angle binary pattern
        t.append(val)
    return np.array(t, dtype=np.float64)

def measure_aligned_leakage(nr, n=N):
    """Aligned leakage: correlation between output popcount pattern
    and the golden template, for golden vs random inputs."""

    golden_template = golden_template_256()

    # Method: for each sample, compute output, measure how many bits
    # match the golden template vs how many match random
    golden_scores = []
    random_scores = []

    for s in range(n):
        # Golden input
        gw = [int((s*16+i) * PHI * MOD32) & MASK32 for i in range(16)]
        gs = sha_partial(gw, nr)

        # Random input
        rw = [int.from_bytes(np.random.bytes(4), 'big') for _ in range(16)]
        rs = sha_partial(rw, nr)

        # Extract output bits
        g_bits = np.zeros(256)
        r_bits = np.zeros(256)
        for bit in range(256):
            w, pos = bit//32, bit%32
            g_bits[bit] = (gs[w] >> pos) & 1
            r_bits[bit] = (rs[w] >> pos) & 1

        # Aligned score: correlation with golden template
        g_score = np.corrcoef(g_bits, golden_template)[0,1]
        r_score = np.corrcoef(r_bits, golden_template)[0,1]

        if not np.isnan(g_score): golden_scores.append(abs(g_score))
        if not np.isnan(r_score): random_scores.append(abs(r_score))

    g_mean = np.mean(golden_scores)
    r_mean = np.mean(random_scores)
    g_std = np.std(golden_scores) / np.sqrt(len(golden_scores))
    r_std = np.std(random_scores) / np.sqrt(len(random_scores))

    # Z-score: is golden mean significantly different from random mean?
    pooled_se = np.sqrt(g_std**2 + r_std**2)
    z = (g_mean - r_mean) / pooled_se if pooled_se > 0 else 0

    return g_mean, r_mean, z

# Also measure: popcount bias (simpler, more robust)
def measure_popcount_bias(nr, n=N):
    """Does golden input produce different popcount than random?"""
    g_pcs = []
    r_pcs = []
    for s in range(n):
        gw = [int((s*16+i) * PHI * MOD32) & MASK32 for i in range(16)]
        gs = sha_partial(gw, nr)

        rw = [int.from_bytes(np.random.bytes(4), 'big') for _ in range(16)]
        rs = sha_partial(rw, nr)

        g_pc = sum(bin(w).count('1') for w in gs)
        r_pc = sum(bin(w).count('1') for w in rs)
        g_pcs.append(g_pc)
        r_pcs.append(r_pc)

    g_mean = np.mean(g_pcs)
    r_mean = np.mean(r_pcs)
    pooled_se = np.sqrt(np.var(g_pcs)/n + np.var(r_pcs)/n)
    z = (g_mean - r_mean) / pooled_se if pooled_se > 0 else 0
    return g_mean, r_mean, z

# Run
test_rounds = [10, 20, 30, 40, 50, 55, 59, 60, 62, 64]

print("SHA ALIGNED LEAKAGE TEST (N=10000)")
print("=" * 75)
print()
print("TEST 1: Aligned correlation with golden template")
print(f"{'Rnd':>4s} {'G_corr':>10s} {'R_corr':>10s} {'z-score':>10s} {'Signal':>8s}")
print("-" * 50)

for nr in test_rounds:
    t0 = time.time()
    g, r, z = measure_aligned_leakage(nr)
    dt = time.time() - t0
    sig = "YES" if abs(z) > 2 else "no"
    print(f"{nr:4d} {g:10.6f} {r:10.6f} {z:+10.2f} {sig:>8s}  ({dt:.1f}s)")

print()
print("TEST 2: Popcount bias (golden vs random)")
print(f"{'Rnd':>4s} {'G_pc':>10s} {'R_pc':>10s} {'z-score':>10s} {'Signal':>8s}")
print("-" * 50)

for nr in test_rounds:
    t0 = time.time()
    g, r, z = measure_popcount_bias(nr)
    dt = time.time() - t0
    sig = "YES" if abs(z) > 2 else "no"
    dodec = {30:'E', 59:'|A5|-1', 60:'|A5|', 64:'2^6'}.get(nr, '')
    print(f"{nr:4d} {g:10.3f} {r:10.3f} {z:+10.2f} {sig:>8s}  {dodec} ({dt:.1f}s)")

print()
print("GEAR MODEL PREDICTION:")
print("  Leakage should PEAK at rounds 55-59 (all build bands, no kill)")
print("  Leakage should DROP at rounds 60-64 (kill zone active)")
print("  The SIGN should FLIP at round 60 = |A5| (phase inversion)")
