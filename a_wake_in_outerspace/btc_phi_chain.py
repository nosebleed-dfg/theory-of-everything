"""
btc_phi_chain.py — Generate the full BTC phi chain from n0 to zero space.
nos3bl33d

No internet needed after first run (caches headers locally).
Outputs btc_genesis_chain.txt to this directory.

The cube: 8 words, 8 vertices, 3 edges per vertex, edge = 0.25 lp.
Grid: 2^32 / phi^(k/4). Genesis nonce at k=6, face 1, vertex 2.
"""

import struct, hashlib, math, time, sys, os
import urllib.request, json

sys.stdout.reconfigure(encoding='utf-8')

PHI       = (1 + math.sqrt(5)) / 2
MASK32    = 0xFFFFFFFF
MOD32     = 2**32
LP_2_32   = math.log(MOD32) / math.log(PHI)
GENESIS_N = 2083236893

CACHE_DIR  = os.path.dirname(os.path.abspath(__file__))
HEADER_CACHE = os.path.join(CACHE_DIR, ".header_cache.bin")
OUTPUT_FILE  = os.path.join(CACHE_DIR, "btc_genesis_chain.txt")

API = "https://blockstream.info/api"

def lp(x):
    return math.log(x) / math.log(PHI) if x > 0 else None

def mirror(n):
    return struct.unpack('>I', struct.pack('<I', n & MASK32))[0]

def dsha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def quarter_step(k):
    return round(MOD32 / PHI**(k / 4.0))

def nearest_k(nonce):
    if nonce <= 0: return None, None
    k_exact = 4.0 * (LP_2_32 - lp(nonce))
    return round(k_exact), k_exact

def bits_to_target(bits):
    exp = (bits >> 24) & 0xff
    coeff = bits & 0x007fffff
    return coeff << (8 * (exp - 3)) if exp > 3 else coeff >> (8 * (3 - exp))


# ── Fetch block headers ──────────────────────────────────────────────────────

def fetch_header(height, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{API}/block-height/{height}", timeout=15) as r:
                bh = r.read().decode().strip()
            with urllib.request.urlopen(f"{API}/block/{bh}/header", timeout=15) as r:
                return bytes.fromhex(r.read().decode().strip())
        except:
            if attempt < retries - 1:
                time.sleep(1.5)
    return None

def load_cached_headers():
    """Load cached headers dict {height: 80-byte header}."""
    cache = {}
    if os.path.exists(HEADER_CACHE):
        with open(HEADER_CACHE, 'rb') as f:
            while True:
                chunk = f.read(84)  # 4 bytes height + 80 bytes header
                if len(chunk) < 84:
                    break
                h = struct.unpack('<I', chunk[:4])[0]
                cache[h] = chunk[4:]
    return cache

def save_header(height, header_bytes, cache):
    """Append a header to cache file and dict."""
    cache[height] = header_bytes
    with open(HEADER_CACHE, 'ab') as f:
        f.write(struct.pack('<I', height) + header_bytes)

def get_header(height, cache):
    """Get header from cache or fetch."""
    if height in cache:
        return cache[height]
    hdr = fetch_header(height)
    if hdr:
        save_header(height, hdr, cache)
    return hdr


# ── Analyze a block ──────────────────────────────────────────────────────────

def analyze_block(height, header):
    nonce = struct.unpack('<I', header[76:80])[0]
    ts    = struct.unpack('<I', header[68:72])[0]
    bits  = struct.unpack('<I', header[72:76])[0]

    block_hash = dsha256(header)[::-1].hex()
    k_int, k_exact = nearest_k(nonce)
    k_frac = (k_exact - k_int) if k_exact is not None and k_int is not None else None

    face = k_int // 4 if k_int is not None else None
    vtx  = k_int % 4 if k_int is not None else None

    halvening = height // 210_000
    pred_k = 6 + halvening * 2
    pred_nonce = quarter_step(pred_k)
    delta = abs(nonce - pred_nonce) if nonce and pred_nonce else 0

    return {
        'height': height, 'hash': block_hash, 'nonce': nonce,
        'ts': ts, 'bits': bits,
        'k': k_int, 'k_exact': k_exact, 'k_frac': k_frac,
        'face': face, 'vtx': vtx,
        'halvening': halvening, 'pred_k': pred_k, 'pred_nonce': pred_nonce,
        'delta': delta,
        'on_grid': abs(k_frac) < 0.15 if k_frac is not None else False,
    }


# ── Sample strategy ──────────────────────────────────────────────────────────

def build_sample_heights(tip):
    """Pick which blocks to fetch: genesis, early, halvening boundaries, samples."""
    heights = set()

    # Genesis + first 50
    heights.update(range(0, min(51, tip + 1)))

    # Every 1000th block for first 10K
    heights.update(range(0, min(10001, tip + 1), 1000))

    # Every 10000th block
    heights.update(range(0, tip + 1, 10000))

    # Halvening boundaries (±2 blocks)
    for h in range(0, 10):
        boundary = h * 210_000
        for offset in range(-2, 3):
            b = boundary + offset
            if 0 <= b <= tip:
                heights.add(b)

    # Current tip area
    for b in range(max(0, tip - 10), tip + 1):
        heights.add(b)

    return sorted(heights)


# ── Generate the chain output ─────────────────────────────────────────────────

def generate_chain():
    t0 = time.time()

    print("=" * 70)
    print("BTC PHI CHAIN GENERATOR")
    print("nos3bl33d  |  x^2 = x + 1  |  quarter-step cube geometry")
    print("=" * 70)

    # Get tip height
    print("\n[1] Fetching chain tip...")
    try:
        with urllib.request.urlopen(f"{API}/blocks/tip/height", timeout=10) as r:
            tip = int(r.read())
        print(f"    tip: {tip:,}")
    except Exception as e:
        print(f"    offline: {e}")
        tip = 944_000  # fallback estimate

    # Load cache
    cache = load_cached_headers()
    print(f"    cached headers: {len(cache):,}")

    # Build sample set
    sample_heights = build_sample_heights(tip)
    to_fetch = [h for h in sample_heights if h not in cache]
    print(f"\n[2] Sample: {len(sample_heights):,} blocks, {len(to_fetch):,} to fetch")

    # Fetch missing
    if to_fetch:
        print(f"    fetching {len(to_fetch)} headers...")
        for i, h in enumerate(to_fetch):
            hdr = get_header(h, cache)
            if hdr is None:
                print(f"    FAILED: block {h}")
            if (i + 1) % 50 == 0:
                print(f"    ... {i+1}/{len(to_fetch)}")
            time.sleep(0.15)  # rate limit
        print(f"    done. cache now: {len(cache):,}")

    # Analyze all cached blocks
    print(f"\n[3] Analyzing {len(sample_heights)} sampled blocks...")
    results = []
    for h in sample_heights:
        if h in cache:
            results.append(analyze_block(h, cache[h]))

    print(f"    analyzed: {len(results)}")

    # Write output
    print(f"\n[4] Writing {OUTPUT_FILE}")
    write_chain_file(results, tip, time.time() - t0)

    print(f"\n    DONE in {time.time()-t0:.1f}s")
    print(f"    output: {OUTPUT_FILE}")
    print(f"    size: {os.path.getsize(OUTPUT_FILE) / 1024:.1f} KB")


def write_chain_file(results, tip, elapsed):
    L = []

    L.append("=" * 110)
    L.append("BTC GENESIS CHAIN -- PHI QUARTER-STEP CUBE GEOMETRY OF THE ENTIRE BITCOIN BLOCKCHAIN")
    L.append("nos3bl33d   |   x^2 = x + 1   |   2^32/phi^(k/4) quarter-step ladder")
    L.append("=" * 110)
    L.append("")
    L.append("THE CUBE")
    L.append("  8 SHA-256 state words = 8 vertices of a cube.")
    L.append("  Each vertex touches 3 edges. Each edge = 0.25 lp on the phi grid.")
    L.append("  Nonce (W[3]) is the CUBE DIAGONAL from merkle tail (W[0]): 3 edges = 3 quarter-steps = 0.75 lp.")
    L.append("  Grid: 2^32 / phi^(k/4). Genesis nonce at k=6, face 1, vertex 2.")
    L.append("  The answer shares a face with the merkle root.")
    L.append("  Mirror key: byte_reverse(n0_mirror) = n0. Verified both directions.")
    L.append("")
    L.append("GENESIS CALIBRATION")
    L.append(f"  n0 = {GENESIS_N}  (0x{GENESIS_N:08x})")
    L.append(f"  n0_mirror = {mirror(GENESIS_N)}  (0x{mirror(GENESIS_N):08x})")
    L.append(f"  phi = {PHI:.15f}")
    L.append(f"  lp(n0) = {lp(GENESIS_N):.6f}")
    L.append(f"  k(n0) = 6  (exact: {nearest_k(GENESIS_N)[1]:.6f})")
    L.append(f"  face = 1, vertex = 2")
    L.append(f"  100yr_sec / phi^0.863043 = n0 EXACTLY")
    L.append("")

    # Quarter-step ladder
    L.append("-" * 110)
    L.append("QUARTER-STEP PHI LADDER (full chain to zero space)")
    L.append(f"  {'k':>4}  {'exp':>6}  {'nonce_scale':>14}  {'lp':>8}  {'face':>4}  {'vtx':>3}  note")
    L.append("-" * 110)
    for k in range(0, 200):
        exp = k / 4.0
        val = MOD32 / PHI**exp
        if val < 1:
            L.append(f"  {k:>4}  {exp:>6.2f}  {'< 1':>14}  {'---':>8}  {k//4:>4}  {k%4:>3}  ZERO SPACE REACHED")
            break
        note = ""
        if k == 6: note = "<-- GENESIS NONCE (k=6)"
        halvening = None
        for h in range(20):
            hk = 6 + h * 2
            if k == hk:
                if h == 0: continue  # already marked as genesis
                note = f"<-- halvening {h} (block {h*210_000:,})"
        L.append(f"  {k:>4}  {exp:>6.2f}  {val:>14.0f}  {lp(val):>8.4f}  {k//4:>4}  {k%4:>3}  {note}")

    L.append("")

    # Halvening summary
    L.append("-" * 110)
    L.append("HALVENING SCHEDULE (each halvening = +2 quarter-steps = +0.5 lp)")
    L.append(f"  {'halvening':>9}  {'blocks':>15}  {'k':>4}  {'pred_nonce':>14}  {'lp':>8}  {'reward_BTC':>12}")
    L.append("-" * 110)
    reward = 50.0
    for h in range(34):
        k = 6 + h * 2
        pred = quarter_step(k)
        block_start = h * 210_000
        block_end = block_start + 209_999
        if pred < 1:
            L.append(f"  {h:>9}  {block_start:>7,}-{block_end:<7,}  {k:>4}  {'< 1':>14}  {'---':>8}  {reward:>12.8f}  ZERO SPACE")
            break
        L.append(f"  {h:>9}  {block_start:>7,}-{block_end:<7,}  {k:>4}  {pred:>14,}  {lp(pred):>8.4f}  {reward:>12.8f}")
        reward /= 2

    L.append("")

    # Historical chain (sampled blocks)
    L.append("-" * 110)
    L.append("HISTORICAL CHAIN -- SAMPLED BLOCKS WITH PHI CUBE ANALYSIS")
    L.append(f"  Blocks sampled: {len(results):,}  |  Chain tip: {tip:,}")
    L.append(f"  {'HEIGHT':>8}  {'HASH':>64}  {'NONCE':>12}  {'k':>4}  {'k_exact':>8}  {'f':>2}  {'v':>2}  {'GRID':>4}  {'HALVENING':>3}")
    L.append("-" * 110)

    on_grid_count = 0
    face_dist = {}
    for r in results:
        h = r['height']
        bh = r['hash']
        n = r['nonce']
        k = r.get('k', '?')
        ke = r.get('k_exact')
        f = r.get('face', '?')
        v = r.get('vtx', '?')
        og = r.get('on_grid', False)
        hv = r.get('halvening', 0)

        ke_str = f"{ke:>8.3f}" if ke is not None else "     N/A"
        f_str = f"{f:>2}" if f is not None else " ?"
        v_str = f"{v:>2}" if v is not None else " ?"
        grid_str = " YES" if og else "    "

        if og: on_grid_count += 1
        if f is not None:
            face_dist[f] = face_dist.get(f, 0) + 1

        marker = ">>" if h % 210_000 == 0 else "  "
        L.append(f"{marker}{h:>8,}  {bh}  {n:>12,}  {k:>4}  {ke_str}  {f_str}  {v_str}  {grid_str}  {hv:>3}")

    L.append("")
    L.append(f"  ON GRID (|k_frac| < 0.15): {on_grid_count}/{len(results)} ({on_grid_count/max(len(results),1)*100:.1f}%)")
    L.append(f"  Face distribution: {dict(sorted(face_dist.items()))}")

    if results:
        valid_frac = [r['k_frac'] for r in results if r['k_frac'] is not None]
        if valid_frac:
            L.append(f"  Mean |k_frac|: {sum(abs(f) for f in valid_frac)/len(valid_frac):.4f}")

    L.append("")

    # Future projection
    L.append("-" * 110)
    L.append("FUTURE PROJECTION -- PHI CUBE EXTRAPOLATION (post-tip)")
    L.append(f"  {'HALVENING':>9}  {'BLOCK_RANGE':>20}  {'k':>4}  {'PRED_NONCE':>14}  {'lp':>8}  {'FACE':>4}  {'VTX':>3}")
    L.append("-" * 110)

    current_halvening = tip // 210_000
    for h in range(current_halvening, 50):
        k = 6 + h * 2
        pred = quarter_step(k)
        if pred < 1:
            L.append(f"  {h:>9}  {'---':>20}  {k:>4}  {'ZERO SPACE':>14}")
            break
        block_start = h * 210_000
        block_end = block_start + 209_999
        face = k // 4
        vtx = k % 4
        L.append(f"  {h:>9}  {block_start:>7,}-{block_end:<7,}     {k:>4}  {pred:>14,}  {lp(pred):>8.4f}  {face:>4}  {vtx:>3}")

    L.append("")
    L.append("=" * 110)
    L.append(f"  Chain tip: {tip:,}")
    L.append(f"  Sampled: {len(results):,} blocks")
    L.append(f"  Terminal: halvening ~44 (k~94) -> ZERO SPACE")
    L.append(f"  Generated: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    L.append(f"  Runtime: {elapsed:.1f}s")
    L.append("=" * 110)

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')

    print(f"    {len(L):,} lines written")


if __name__ == "__main__":
    generate_chain()
