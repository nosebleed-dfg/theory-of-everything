"""
DIGIT MIRROR SOLVER — test nonce recovery via decimal/hex digit mirroring
nos3bl33d

Hypothesis: nonce = T1_round3 - known, where known is computable from block header.
Mirror known's/nonce's digits around center, drop evens, compare residual.

6 approaches tested on Bitcoin blocks 0-49:
  dec_fwd:  mirror(known) decimal, keep odd -> compare to nonce
  dec_rev:  mirror(nonce) decimal, keep odd -> compare to known's odd digits
  hex_fwd:  mirror(known) hex (d->15-d), keep odd -> compare to nonce
  hex_rev:  mirror(nonce) hex, keep odd -> compare to known's odd hex digits
  t1_dec:   mirror(T1) decimal, keep odd -> compare to nonce
  t1_hex:   mirror(T1) hex, keep odd -> compare to nonce
  xor_hex:  mirror(nonce^known) hex, keep odd -> compare to known's odd hex digits
"""

import struct
import time
import math

# ============================================================
# SHA-256 constants and primitives
# ============================================================

MASK32 = 0xFFFFFFFF

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

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def ch(e, f, g):
    return ((e & f) ^ (~e & g)) & MASK32

def maj(a, b, c):
    return ((a & b) ^ (a & c) ^ (b & c)) & MASK32

def sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def lsigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def lsigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def add32(*args):
    return sum(args) & MASK32

# ============================================================
# SHA-256 compression with state tracking
# ============================================================

def sha256_compress(block_bytes, iv):
    """Compress one 64-byte block. Returns (final_hash, W, states).
    states[i] = state BEFORE round i. states[0] = iv copy.
    """
    W = list(struct.unpack('>16I', block_bytes))
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))

    st = list(iv)
    states = [list(st)]
    for i in range(64):
        a, b, c, d, e, f, g, h = st
        T1 = add32(h, sigma1(e), ch(e, f, g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a, b, c))
        st = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
        states.append(list(st))

    final = [add32(st[j], iv[j]) for j in range(8)]
    return final, W, states

# ============================================================
# Digit mirror operations
# ============================================================

def dec_mirror_str(val):
    """Mirror decimal digits: d -> (10-d)%10. Return as string."""
    return ''.join(str((10 - int(c)) % 10) for c in str(val))

def dec_odd_str(s):
    """Keep only odd decimal digits from string."""
    return ''.join(c for c in s if int(c) % 2 == 1)

def hex_mirror_str(val):
    """Mirror hex digits: d -> 15-d. Zero-pad to 8 hex digits. Return as string."""
    s = format(val, '08x')
    return ''.join(format(15 - int(c, 16), 'x') for c in s)

def hex_odd_str(s):
    """Keep only odd hex digits from string."""
    return ''.join(c for c in s if int(c, 16) % 2 == 1)

def str_to_int(s, base=10):
    """Convert digit string to int. Returns 0 for empty string."""
    if not s:
        return 0
    return int(s, base)

# ============================================================
# Header loading and block2 construction
# ============================================================

def load_headers(path, count=50):
    """Load block headers from binary cache. Format: 4B LE height + 80B header."""
    headers = []
    with open(path, 'rb') as f:
        for _ in range(count):
            entry = f.read(84)
            if len(entry) < 84:
                break
            height = struct.unpack_from('<I', entry, 0)[0]
            hdr = entry[4:84]
            headers.append((height, hdr))
    return headers

def build_block2(hdr):
    """64-byte SHA-256 block 2 from 80-byte Bitcoin header."""
    tail = hdr[64:80]
    padding = b'\x80' + b'\x00' * (64 - 16 - 1 - 8) + struct.pack('>Q', 80 * 8)
    return tail + padding

# ============================================================
# Per-block analysis
# ============================================================

def analyze_block(height, hdr):
    """Compute all mirror approaches for one block. Returns result dict."""
    nonce_be = struct.unpack('>I', hdr[76:80])[0]
    nonce_le = struct.unpack('<I', hdr[76:80])[0]

    # SHA-256 block 1 -> H_mid
    h_mid, _, _ = sha256_compress(hdr[0:64], H0)

    # SHA-256 block 2 -> 3 rounds to get state[3]
    block2 = build_block2(hdr)
    _, W, states = sha256_compress(block2, h_mid)

    assert W[3] == nonce_be

    # State entering round 3
    a, b, c, d, e, f, g, h = states[3]
    known = add32(h, sigma1(e), ch(e, f, g), K[3])
    T1 = add32(known, W[3])

    # Precompute hex/dec representations
    n_hex = format(nonce_be, '08x')
    k_hex = format(known, '08x')
    n_dec = str(nonce_be)
    k_dec = str(known)

    # ---- APPROACH 1: dec_fwd — mirror(known) -> odd -> compare to nonce ----
    km_dec = dec_mirror_str(known)
    km_odd = dec_odd_str(km_dec)
    dec_fwd_val = str_to_int(km_odd)
    d_dec_fwd = abs(nonce_be - dec_fwd_val)

    # ---- APPROACH 2: dec_rev — mirror(nonce) -> odd -> compare to known's odd ----
    nm_dec = dec_mirror_str(nonce_be)
    nm_odd = dec_odd_str(nm_dec)
    nm_val = str_to_int(nm_odd)
    k_odd_dec = dec_odd_str(k_dec)
    k_odd_val = str_to_int(k_odd_dec)
    d_dec_rev = abs(nm_val - k_odd_val)

    # ---- APPROACH 3: hex_fwd — mirror(known) hex -> odd -> compare to nonce ----
    km_hex = hex_mirror_str(known)
    km_hodd = hex_odd_str(km_hex)
    hex_fwd_val = str_to_int(km_hodd, 16)
    d_hex_fwd = abs(nonce_be - hex_fwd_val)

    # ---- APPROACH 4: hex_rev — mirror(nonce) hex -> odd -> compare to known's odd hex ----
    nm_hex = hex_mirror_str(nonce_be)
    nm_hodd = hex_odd_str(nm_hex)
    nm_hval = str_to_int(nm_hodd, 16)
    k_hodd = hex_odd_str(k_hex)
    k_hval = str_to_int(k_hodd, 16)
    d_hex_rev = abs(nm_hval - k_hval)

    # ---- APPROACH 5: t1_dec — mirror(T1) decimal -> odd -> compare to nonce ----
    t1m_dec = dec_mirror_str(T1)
    t1m_odd = dec_odd_str(t1m_dec)
    t1_dec_val = str_to_int(t1m_odd)
    d_t1_dec = abs(nonce_be - t1_dec_val)

    # ---- APPROACH 6: t1_hex — mirror(T1) hex -> odd -> compare to nonce ----
    t1m_hex = hex_mirror_str(T1)
    t1m_hodd = hex_odd_str(t1m_hex)
    t1_hex_val = str_to_int(t1m_hodd, 16)
    d_t1_hex = abs(nonce_be - t1_hex_val)

    # ---- APPROACH 7: xor_hex — mirror(nonce^known) hex -> odd -> compare to known odd hex ----
    xor_val = nonce_be ^ known
    xm_hex = hex_mirror_str(xor_val)
    xm_hodd = hex_odd_str(xm_hex)
    xm_hval = str_to_int(xm_hodd, 16)
    d_xor_hex = abs(xm_hval - k_hval)

    deltas = {
        'dec_fwd': d_dec_fwd,
        'dec_rev': d_dec_rev,
        'hex_fwd': d_hex_fwd,
        'hex_rev': d_hex_rev,
        't1_dec':  d_t1_dec,
        't1_hex':  d_t1_hex,
        'xor_hex': d_xor_hex,
    }
    best = min(deltas, key=deltas.get)

    return {
        'height': height,
        'nonce_le': nonce_le,
        'nonce_be': nonce_be,
        'known': known,
        'T1': T1,
        'deltas': deltas,
        'best': best,
        'best_delta': deltas[best],
        # Detail strings for verbose output
        'detail': {
            'dec_fwd': f"mirror({k_dec})={km_dec} -> odd={km_odd} -> {dec_fwd_val}",
            'dec_rev': f"mirror({n_dec})={nm_dec} -> odd={nm_odd}={nm_val}  vs  known_odd={k_odd_dec}={k_odd_val}",
            'hex_fwd': f"mirror({k_hex})={km_hex} -> odd={km_hodd} -> 0x{hex_fwd_val:x}",
            'hex_rev': f"mirror({n_hex})={nm_hex} -> odd={nm_hodd}=0x{nm_hval:x}  vs  known_odd={k_hodd}=0x{k_hval:x}",
            't1_dec': f"mirror({T1})={t1m_dec} -> odd={t1m_odd} -> {t1_dec_val}",
            't1_hex': f"mirror({format(T1,'08x')})={t1m_hex} -> odd={t1m_hodd} -> 0x{t1_hex_val:x}",
            'xor_hex': f"mirror({format(xor_val,'08x')})={xm_hex} -> odd={xm_hodd}=0x{xm_hval:x}  vs  known_odd={k_hodd}=0x{k_hval:x}",
        },
    }

# ============================================================
# Main
# ============================================================

def main():
    t0 = time.perf_counter()

    cache_path = 'C:/Users/funct/Desktop/axiom/a_wake_in_outerspace/.header_cache.bin'
    headers = load_headers(cache_path, 50)
    print(f"Loaded {len(headers)} block headers\n")

    approach_names = ['dec_fwd', 'dec_rev', 'hex_fwd', 'hex_rev', 't1_dec', 't1_hex', 'xor_hex']

    # ============================================================
    # TABLE 1: All deltas, all blocks
    # ============================================================
    print("=" * 140)
    print(f"{'Blk':>4}  {'Nonce(BE)':>10}  {'Known':>10}  "
          + "  ".join(f"{n:>10}" for n in approach_names)
          + f"  {'BEST':>8}")
    print("-" * 140)

    all_results = []
    wins = {n: 0 for n in approach_names}
    exacts = {n: 0 for n in approach_names}

    for height, hdr in headers:
        r = analyze_block(height, hdr)
        all_results.append(r)
        wins[r['best']] += 1
        for n in approach_names:
            if r['deltas'][n] == 0:
                exacts[n] += 1

        row = f"{r['height']:>4}  0x{r['nonce_be']:08x}  0x{r['known']:08x}  "
        row += "  ".join(f"{r['deltas'][n]:>10}" for n in approach_names)
        row += f"  {r['best']:>8}"
        if r['best_delta'] == 0:
            row += "  *** EXACT ***"
        print(row)

    elapsed = time.perf_counter() - t0

    # ============================================================
    # TABLE 2: Summary statistics
    # ============================================================
    print(f"\n{'=' * 80}")
    print(f"SUMMARY -- {len(headers)} blocks in {elapsed:.3f}s")
    print(f"{'=' * 80}\n")

    print(f"{'Approach':>10}  {'Wins':>5}  {'Exact':>5}  {'Avg Delta':>12}  {'Median Delta':>12}  {'Max Delta':>12}  {'Med Bits Saved':>14}")
    print("-" * 85)
    for n in approach_names:
        ds = sorted([r['deltas'][n] for r in all_results])
        avg_d = sum(ds) / len(ds)
        med_d = ds[len(ds) // 2]
        max_d = ds[-1]
        med_bits = 32 - math.log2(2 * med_d + 1) if med_d > 0 else 32
        print(f"{n:>10}  {wins[n]:>5}  {exacts[n]:>5}  {avg_d:>12.0f}  {med_d:>12}  {max_d:>12}  {med_bits:>14.1f}")

    # ============================================================
    # TABLE 3: Search space reduction
    # ============================================================
    print(f"\n{'=' * 80}")
    print(f"SEARCH SPACE REDUCTION")
    print(f"{'=' * 80}\n")
    print(f"How many blocks have delta < threshold:")
    print(f"{'Approach':>10}", end="")
    for bits in [10, 14, 16, 18, 20, 24]:
        print(f"  {'<2^'+str(bits):>6}", end="")
    print()
    print("-" * 55)
    for n in approach_names:
        ds = [r['deltas'][n] for r in all_results]
        print(f"{n:>10}", end="")
        for bits in [10, 14, 16, 18, 20, 24]:
            count = sum(1 for d in ds if d < 2**bits)
            print(f"  {count:>4}/50", end="")
        print()

    # ============================================================
    # TABLE 4: Detailed view of blocks with smallest best-delta
    # ============================================================
    print(f"\n{'=' * 80}")
    print(f"TOP 10 CLOSEST MATCHES (sorted by best delta)")
    print(f"{'=' * 80}")

    top10 = sorted(all_results, key=lambda r: r['best_delta'])[:10]
    for r in top10:
        print(f"\n--- Block {r['height']} ---")
        print(f"  Nonce (LE): 0x{r['nonce_le']:08x} = {r['nonce_le']}")
        print(f"  Nonce (BE): 0x{r['nonce_be']:08x} = {r['nonce_be']}")
        print(f"  Known:      0x{r['known']:08x} = {r['known']}")
        print(f"  T1 round 3: 0x{r['T1']:08x} = {r['T1']}")
        for n in approach_names:
            marker = " <-- BEST" if n == r['best'] else ""
            print(f"  [{n:>7}] delta={r['deltas'][n]:>12}  |  {r['detail'][n]}{marker}")

    # ============================================================
    # TABLE 5: Bit-level view of nonce vs known for hex_rev winners
    # ============================================================
    print(f"\n{'=' * 80}")
    print(f"BIT PATTERN: nonce XOR known for all blocks")
    print(f"{'=' * 80}\n")
    print(f"{'Blk':>4}  {'nonce^known':>10}  {'hamming':>7}  {'popcount(N)':>11}  {'popcount(K)':>11}")
    print("-" * 55)
    for r in all_results:
        xor = r['nonce_be'] ^ r['known']
        hamming = bin(xor).count('1')
        pop_n = bin(r['nonce_be']).count('1')
        pop_k = bin(r['known']).count('1')
        print(f"{r['height']:>4}  0x{xor:08x}  {hamming:>7}  {pop_n:>11}  {pop_k:>11}")


if __name__ == '__main__':
    main()
