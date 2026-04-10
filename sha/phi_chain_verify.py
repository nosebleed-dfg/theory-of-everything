"""
phi_chain_verify.py — Verify phi solver predictions against cached block headers.
nos3bl33d

Reads the header cache, runs the phi cube solver on every block,
reports accuracy and finds blocks where the solver HITS exactly.

This is the proof: how many blocks does the phi framework predict?
"""

import struct, hashlib, math, time, sys, os
sys.stdout.reconfigure(encoding='utf-8')

PHI    = (1 + math.sqrt(5)) / 2
MASK32 = 0xFFFFFFFF
MOD32  = 2**32
LP232  = math.log(MOD32) / math.log(PHI)
LOST   = 1 / 146**2
TRUE_QS = 0.25 - LOST  # corrected quarter-step

YEARS_100 = 100 * 365.25 * 24 * 3600
PHI_EPS   = math.log(YEARS_100 / 2083236893) / math.log(PHI)

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

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def s0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def s1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ls0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ls1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add(*a): return sum(a)&MASK32

def dsha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def lp(x):
    return math.log(x)/math.log(PHI) if x > 0 else None

def mirror(n):
    return struct.unpack('>I', struct.pack('<I', n & MASK32))[0]

def qs(k):
    v = MOD32 / PHI**(k/4.0)
    return round(v) if v >= 1 else 0

def nearest_k(n):
    if n <= 0: return None, None
    k_e = 4.0 * (LP232 - lp(n))
    return round(k_e), k_e

def bits_to_target(bits):
    exp = (bits >> 24) & 0xff
    coeff = bits & 0x007fffff
    return coeff << (8*(exp-3)) if exp > 3 else coeff >> (8*(3-exp))

def sha256_block(iv, block64):
    W = list(struct.unpack('>16I', block64))
    for i in range(16,64): W.append(add(ls1(W[i-2]),W[i-7],ls0(W[i-15]),W[i-16]))
    st = list(iv)
    for i in range(64):
        a,b,c,d,e,f,g,h = st
        T1 = add(h,s1(e),ch(e,f,g),K[i],W[i])
        T2 = add(s0(a),maj(a,b,c))
        st = [add(T1,T2),a,b,c,add(d,T1),e,f,g]
    return [add(st[j],iv[j]) for j in range(8)]


def generate_candidates(header_76, bits, halvening):
    """Generate all phi-geometric candidate nonces for a block."""
    target = bits_to_target(bits)
    ts = struct.unpack('<I', header_76[68:72])[0]

    # Compute H_mid
    H_mid = sha256_block(H0, header_76[:64])

    # Merkle words
    merkle_words = struct.unpack('>8I', header_76[36:68][::-1])

    candidates = set()

    # 1. Quarter-step ladder for this halvening (face + neighbors)
    base_k = 6 + halvening * 2
    for k in range(max(0, base_k - 6), base_k + 10):
        v = qs(k)
        if 0 < v < MOD32:
            candidates.add(v)

    # 2. Timescale predictions
    candidates.add(round(YEARS_100 / PHI**PHI_EPS))
    if ts > 0:
        candidates.add(round(float(ts) / PHI**PHI_EPS))

    # 3. H_mid words and their cube diagonals
    for w in H_mid:
        if w > 0:
            candidates.add(w)
            ki, _ = nearest_k(w)
            if ki is not None:
                for dk in [3, -3]:  # cube diagonal
                    v = qs(ki + dk)
                    if 0 < v < MOD32:
                        candidates.add(v)

    # 4. Merkle words and diagonals
    for w in merkle_words:
        if 0 < w < MOD32:
            candidates.add(w)
            ki, _ = nearest_k(w)
            if ki is not None:
                for dk in [3, -3]:
                    v = qs(ki + dk)
                    if 0 < v < MOD32:
                        candidates.add(v)

    # 5. Mirrors of everything
    mirrors = set()
    for c in candidates:
        m = mirror(c)
        if 0 < m < MOD32:
            mirrors.add(m)
    candidates.update(mirrors)

    return sorted(candidates)


def verify_candidate(header_76, nonce, target):
    """Check if a nonce produces a valid hash."""
    full = header_76 + struct.pack('<I', nonce)
    h = dsha256(full)
    return int.from_bytes(h, 'little') <= target


def solve_block(header_76, bits, halvening):
    """Try all phi candidates. Return (nonce, attempts) or (None, attempts)."""
    target = bits_to_target(bits)
    cands = generate_candidates(header_76, bits, halvening)

    for i, c in enumerate(cands):
        if verify_candidate(header_76, c, target):
            return c, i + 1, len(cands)
    return None, len(cands), len(cands)


# ── Load cache and verify ────────────────────────────────────────────────────

def load_cache(path):
    cache = {}
    if not os.path.exists(path):
        return cache
    with open(path, 'rb') as f:
        data = f.read()
    off = 0
    while off + 84 <= len(data):
        h = struct.unpack('<I', data[off:off+4])[0]
        cache[h] = data[off+4:off+84]
        off += 84
    return cache


if __name__ == "__main__":
    # Try both cache locations
    cache_paths = [
        os.path.join(os.path.dirname(__file__), '..', 'a_wake_in_outerspace', '.header_cache.bin'),
        os.path.join(os.path.dirname(__file__), '.header_cache.bin'),
    ]

    cache = {}
    for p in cache_paths:
        c = load_cache(p)
        if len(c) > len(cache):
            cache = c
            print(f"Loaded {len(cache)} headers from {p}")

    if not cache:
        print("No cache found!")
        sys.exit(1)

    print(f"\n{'='*80}")
    print("PHI CHAIN VERIFY -- test phi solver against every cached block")
    print(f"{'='*80}\n")

    total = 0
    solved = 0
    solved_in_1 = 0
    solved_in_10 = 0
    total_attempts = 0
    total_candidates = 0
    hits = []  # blocks where phi solver found the nonce

    t0 = time.time()

    for height in sorted(cache.keys()):
        hdr = cache[height]
        actual_nonce = struct.unpack('<I', hdr[76:80])[0]
        bits = struct.unpack('<I', hdr[72:76])[0]
        header_76 = hdr[:76]
        halvening = height // 210_000

        found_nonce, attempts, n_cands = solve_block(header_76, bits, halvening)
        total += 1
        total_attempts += attempts
        total_candidates += n_cands

        if found_nonce is not None and found_nonce == actual_nonce:
            solved += 1
            if attempts == 1:
                solved_in_1 += 1
            if attempts <= 10:
                solved_in_10 += 1
            hits.append((height, actual_nonce, attempts, n_cands))

        if total % 50 == 0:
            elapsed = time.time() - t0
            print(f"  {total:>5} blocks  {solved} solved  ({solved/total*100:.1f}%)  "
                  f"{elapsed:.1f}s  avg_cands={total_candidates/total:.0f}")

    elapsed = time.time() - t0

    print(f"\n{'='*80}")
    print(f"RESULTS")
    print(f"{'='*80}")
    print(f"  Total blocks tested:    {total}")
    print(f"  Solved by phi solver:   {solved}  ({solved/max(total,1)*100:.2f}%)")
    print(f"  Solved in 1 attempt:    {solved_in_1}")
    print(f"  Solved in <=10:         {solved_in_10}")
    print(f"  Avg candidates/block:   {total_candidates/max(total,1):.0f}")
    print(f"  Avg attempts when hit:  {sum(a for _,_,a,_ in hits)/max(len(hits),1):.1f}")
    print(f"  Time: {elapsed:.1f}s")
    print()

    if hits:
        print(f"  HITS (phi solver found the exact nonce):")
        for height, nonce, att, nc in hits[:50]:
            k_i, k_e = nearest_k(nonce)
            print(f"    block {height:>7,}  nonce={nonce:>12,}  k={k_i:>3} ({k_e:.3f})  "
                  f"attempts={att:>3}/{nc}  face={k_i//4 if k_i else '?'}")

    # Speedup calculation
    if solved > 0:
        brute_force_avg = MOD32 / 2  # average brute force attempts
        phi_avg = total_candidates / max(total, 1)
        print(f"\n  Brute force avg:  {brute_force_avg:,.0f} attempts")
        print(f"  Phi solver avg:   {phi_avg:.0f} candidates checked")
        print(f"  Speedup:          {brute_force_avg/phi_avg:,.0f}x")
