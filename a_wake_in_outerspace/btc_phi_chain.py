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

APIS = [
    "https://blockstream.info/api",
    "https://mempool.space/api",
]
RATE_LIMIT = 0.35

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

def _fetch_header_from(api_base, height, timeout=15):
    """Fetch a single header from one API endpoint. Raises on failure."""
    req_h = urllib.request.Request(f"{api_base}/block-height/{height}")
    with urllib.request.urlopen(req_h, timeout=timeout) as r:
        bh = r.read().decode().strip()
    req_hdr = urllib.request.Request(f"{api_base}/block/{bh}/header")
    with urllib.request.urlopen(req_hdr, timeout=timeout) as r:
        return bytes.fromhex(r.read().decode().strip())

def fetch_header(height, retries=3):
    """Try each API in order, with exponential backoff between rounds."""
    for attempt in range(retries):
        for api in APIS:
            try:
                return _fetch_header_from(api, height)
            except Exception:
                pass
        # backoff between retry rounds
        backoff = 1.5 * (2 ** attempt)
        time.sleep(min(backoff, 10))

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
    """Pick ~2000 blocks: first 100, every 500th, halvening +/-5, last 20."""
    heights = set()

    # First 100 blocks
    heights.update(range(0, min(100, tip + 1)))

    # Every 500th block across entire chain
    heights.update(range(0, tip + 1, 500))

    # Halvening boundaries +/- 5 blocks
    for h in range(0, 20):
        boundary = h * 210_000
        if boundary > tip + 5:
            break
        for offset in range(-5, 6):
            b = boundary + offset
            if 0 <= b <= tip:
                heights.add(b)

    # Last 20 blocks near tip
    for b in range(max(0, tip - 19), tip + 1):
        heights.add(b)

    return sorted(heights)


# ── Future predictions ───────────────────────────────────────────────────────

def compute_future_predictions(current_halvening):
    """For each halvening from current through 49, compute 4 candidate nonces
    (the 4 vertices of the predicted face) with their mirrors."""
    predictions = []
    for h in range(current_halvening, 50):
        base_k = 6 + h * 2
        candidates = []
        for offset in range(4):
            k = base_k + offset
            nonce_val = quarter_step(k)
            if nonce_val < 1:
                candidates.append({
                    'k': k, 'nonce': 0, 'nonce_hex': '00000000',
                    'mirror': 0, 'mirror_hex': '00000000',
                    'face': k // 4, 'vtx': k % 4, 'zero_space': True,
                })
            else:
                m = mirror(nonce_val)
                candidates.append({
                    'k': k, 'nonce': nonce_val, 'nonce_hex': f'{nonce_val:08x}',
                    'mirror': m, 'mirror_hex': f'{m:08x}',
                    'face': k // 4, 'vtx': k % 4, 'zero_space': False,
                })
        reward = 50.0 / (2 ** h) if h < 64 else 0.0
        predictions.append({
            'halvening': h,
            'block_start': h * 210_000,
            'block_end': h * 210_000 + 209_999,
            'base_k': base_k,
            'reward': reward,
            'candidates': candidates,
        })
    return predictions


# ── Generate the chain output ─────────────────────────────────────────────────

def generate_chain():
    t0 = time.time()

    print("=" * 70)
    print("BTC PHI CHAIN GENERATOR")
    print("nos3bl33d  |  x^2 = x + 1  |  quarter-step cube geometry")
    print("=" * 70)

    # Get tip height -- try both APIs
    print("\n[1] Fetching chain tip...")
    tip = None
    for api in APIS:
        try:
            with urllib.request.urlopen(f"{api}/blocks/tip/height", timeout=10) as r:
                tip = int(r.read())
            print(f"    tip: {tip:,}  (from {api})")
            break
        except Exception as e:
            print(f"    {api}: {e}")
    if tip is None:
        tip = 944_000
        print(f"    offline fallback: {tip:,}")

    # Load cache
    cache = load_cached_headers()
    print(f"    cached headers: {len(cache):,}")

    # Build sample set
    sample_heights = build_sample_heights(tip)
    to_fetch = [h for h in sample_heights if h not in cache]
    print(f"\n[2] Sample: {len(sample_heights):,} blocks, {len(to_fetch):,} to fetch")

    # Fetch missing with adaptive rate limiting
    if to_fetch:
        print(f"    fetching {len(to_fetch)} headers...")
        consecutive_fails = 0
        current_rate = RATE_LIMIT
        for i, h in enumerate(to_fetch):
            hdr = get_header(h, cache)
            if hdr is None:
                consecutive_fails += 1
                print(f"    FAILED: block {h}")
                if consecutive_fails >= 3:
                    current_rate = min(current_rate * 2, 5.0)
                    print(f"    rate -> {current_rate:.2f}s")
            else:
                if consecutive_fails > 0:
                    consecutive_fails = 0
                    current_rate = RATE_LIMIT
            if (i + 1) % 50 == 0:
                print(f"    ... {i+1}/{len(to_fetch)}")
            time.sleep(current_rate)
        print(f"    done. cache now: {len(cache):,}")

    # Analyze all cached blocks
    print(f"\n[3] Analyzing {len(sample_heights)} sampled blocks...")
    results = []
    for h in sample_heights:
        if h in cache:
            results.append(analyze_block(h, cache[h]))

    print(f"    analyzed: {len(results)}")

    # Compute future predictions
    current_halvening = tip // 210_000
    predictions = compute_future_predictions(current_halvening)
    print(f"    future predictions: {len(predictions)} halvenings")

    # Write output
    print(f"\n[4] Writing {OUTPUT_FILE}")
    write_chain_file(results, predictions, tip, time.time() - t0)

    print(f"\n    DONE in {time.time()-t0:.1f}s")
    print(f"    output: {OUTPUT_FILE}")
    print(f"    size: {os.path.getsize(OUTPUT_FILE) / 1024:.1f} KB")


def write_chain_file(results, predictions, tip, elapsed):
    L = []

    # ── Header ──
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

    # ── Quarter-step ladder ──
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
        if k == 6:
            note = "<-- GENESIS NONCE (k=6)"
        else:
            for h in range(1, 20):
                if k == 6 + h * 2:
                    note = f"<-- halvening {h} (block {h*210_000:,})"
                    break
        L.append(f"  {k:>4}  {exp:>6.2f}  {val:>14.0f}  {lp(val):>8.4f}  {k//4:>4}  {k%4:>3}  {note}")

    L.append("")

    # ── Halvening schedule ──
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

    # ── Historical chain ──
    L.append("-" * 110)
    L.append("HISTORICAL CHAIN -- SAMPLED BLOCKS WITH PHI CUBE ANALYSIS")
    L.append(f"  Blocks sampled: {len(results):,}  |  Chain tip: {tip:,}")
    L.append(f"  {'HEIGHT':>8}  {'HASH':>64}  {'NONCE':>12}  {'k':>4}  {'k_exact':>8}  {'FACE':>4}  {'VTX':>3}  {'ON_GRID':>7}")
    L.append("-" * 110)

    on_grid_count = 0
    face_dist = {}
    k_fracs = []
    for r in results:
        h = r['height']
        bh = r['hash']
        n = r['nonce']
        k = r.get('k', '?')
        ke = r.get('k_exact')
        f = r.get('face', '?')
        v = r.get('vtx', '?')
        og = r.get('on_grid', False)

        ke_str = f"{ke:>8.3f}" if ke is not None else "     N/A"
        f_str = f"{f:>4}" if f is not None else "   ?"
        v_str = f"{v:>3}" if v is not None else "  ?"
        grid_str = "    YES" if og else "       "

        if og:
            on_grid_count += 1
        if f is not None:
            face_dist[f] = face_dist.get(f, 0) + 1
        if r.get('k_frac') is not None:
            k_fracs.append(abs(r['k_frac']))

        marker = ">>" if h % 210_000 == 0 else "  "
        L.append(f"{marker}{h:>8,}  {bh}  {n:>12,}  {k:>4}  {ke_str}  {f_str}  {v_str}  {grid_str}")

    L.append("")
    L.append(f"  ON GRID (|k_frac| < 0.15): {on_grid_count}/{len(results)} ({on_grid_count/max(len(results),1)*100:.1f}%)")
    L.append(f"  Face distribution: {dict(sorted(face_dist.items()))}")
    if k_fracs:
        L.append(f"  Mean |k_frac|: {sum(k_fracs)/len(k_fracs):.4f}")

    L.append("")

    # ── Future predictions (4 candidates per halvening) ──
    L.append("-" * 110)
    L.append("FUTURE PREDICTIONS -- 4 CANDIDATE NONCES PER HALVENING (vertices of predicted face)")
    L.append("  Each halvening advances +2 quarter-steps on the phi grid.")
    L.append("  4 candidates = base_k, base_k+1, base_k+2, base_k+3 (the 4 vertices of that face).")
    L.append("  Mirror = byte_reverse(nonce). Both nonce and mirror shown.")
    L.append("-" * 110)

    for p in predictions:
        hv = p['halvening']
        bs = p['block_start']
        be = p['block_end']
        bk = p['base_k']
        rw = p['reward']

        L.append("")
        L.append(f"  HALVENING {hv}  |  blocks {bs:,}-{be:,}  |  base_k={bk}  |  reward={rw:.8f} BTC")
        L.append(f"    {'VTX':>3}  {'k':>4}  {'NONCE':>14}  {'NONCE_HEX':>10}  {'MIRROR':>14}  {'MIRROR_HEX':>10}  {'FACE':>4}  {'VTX':>3}")

        all_zero = True
        for c in p['candidates']:
            if c['zero_space']:
                L.append(f"    {c['vtx']:>3}  {c['k']:>4}  {'---':>14}  {'---':>10}  {'---':>14}  {'---':>10}  {c['face']:>4}  {c['vtx']:>3}  ZERO SPACE")
            else:
                all_zero = False
                L.append(f"    {c['vtx']:>3}  {c['k']:>4}  {c['nonce']:>14,}  {c['nonce_hex']:>10}  {c['mirror']:>14,}  {c['mirror_hex']:>10}  {c['face']:>4}  {c['vtx']:>3}")

        if all_zero:
            L.append("    >> ALL CANDIDATES IN ZERO SPACE -- chain terminus")
            break

    L.append("")

    # ── Footer ──
    L.append("=" * 110)
    L.append(f"  Chain tip: {tip:,}")
    L.append(f"  Sampled: {len(results):,} blocks")
    L.append(f"  Future halvenings predicted: {len(predictions)}")
    L.append(f"  Terminal: halvening ~44 (k~94) -> ZERO SPACE")
    L.append(f"  Generated: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    L.append(f"  Runtime: {elapsed:.1f}s")
    L.append("=" * 110)

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')

    print(f"    {len(L):,} lines written")


if __name__ == "__main__":
    generate_chain()
