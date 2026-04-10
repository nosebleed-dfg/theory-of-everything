"""
phi_cube_solver.py — Bitcoin nonce solver via phi QUARTER-STEP cube geometry
nos3bl33d

THE CUBE:
  8 SHA-256 state words = 8 vertices of a cube.
  Each vertex touches 3 edges. 3 edges = 3 quarter-steps = 0.75 lp.
  The nonce (W[3]) is the CUBE DIAGONAL from the merkle tail (W[0]):
    W[0] -> W[1] -> W[2] -> W[3]  =  n+3  =  3 edges  =  3 quarter-steps

  Old ladder: 2^32 / phi^(k/2)   [half-steps, k=3 for genesis]
  New ladder: 2^32 / phi^(k/4)   [quarter-steps, k=6 for genesis]

  Each edge = 0.25 lp on the phi grid.
  Each face = 4 vertices = 4 quarter-steps.
  The answer (hash) shares a face with the merkle root.

  Mirror distance: lp(n0) - lp(n0_mirror) = 2.97 lp = ~12 quarter-steps
                   12 = 3 * 4 = 3 full cube traversals

x^2 = x + 1.
"""

import struct, hashlib, math, time, sys, os
import urllib.request, json

sys.stdout.reconfigure(encoding='utf-8')

# ── Constants ────────────────────────────────────────────────────────────────
PHI    = (1 + math.sqrt(5)) / 2
MASK32 = 0xFFFFFFFF
MOD32  = 2**32

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

YEARS_100_SEC   = 100 * 365.25 * 24 * 3600
PHI_EPS_GENESIS = math.log(YEARS_100_SEC / 2083236893) / math.log(PHI)

LP_2_32 = math.log(2**32) / math.log(PHI)   # 46.0934...

# ── SHA-256 ──────────────────────────────────────────────────────────────────
def _rotr(x, n):  return ((x >> n) | (x << (32 - n))) & MASK32
def _ch(e,f,g):   return (e & f) ^ (~e & g) & MASK32
def _maj(a,b,c):  return (a & b) ^ (a & c) ^ (b & c)
def _s0(x):       return _rotr(x,2)  ^ _rotr(x,13) ^ _rotr(x,22)
def _s1(x):       return _rotr(x,6)  ^ _rotr(x,11) ^ _rotr(x,25)
def _ls0(x):      return _rotr(x,7)  ^ _rotr(x,18) ^ (x >> 3)
def _ls1(x):      return _rotr(x,17) ^ _rotr(x,19) ^ (x >> 10)
def _add(*a):     return sum(a) & MASK32

def _build_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(_add(_ls1(W[i-2]), W[i-7], _ls0(W[i-15]), W[i-16]))
    return W

def _compress(iv, W64):
    s = list(iv)
    for i in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = _add(h, _s1(e), _ch(e,f,g), K[i], W64[i])
        T2 = _add(_s0(a), _maj(a,b,c))
        s  = [_add(T1,T2), a, b, c, _add(d,T1), e, f, g]
    return [_add(s[j], iv[j]) for j in range(8)]

def _sha256_block(iv, block_bytes):
    W16 = list(struct.unpack('>16I', block_bytes))
    return _compress(iv, _build_schedule(W16))

def dsha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

# ── Phi / lp ─────────────────────────────────────────────────────────────────
def lp(x):
    return math.log(x) / math.log(PHI) if x > 0 else None

def mirror_nonce(n):
    return struct.unpack('>I', struct.pack('<I', n & MASK32))[0]

def bits_to_target(bits):
    exp   = (bits >> 24) & 0xff
    coeff = bits & 0x007fffff
    return coeff << (8 * (exp - 3)) if exp > 3 else coeff >> (8 * (3 - exp))

# ── QUARTER-STEP CUBE GEOMETRY ───────────────────────────────────────────────

def quarter_step(k):
    """Position k on the quarter-step grid: 2^32 / phi^(k/4)."""
    return round(MOD32 / PHI**(k / 4.0))

def nearest_k(nonce):
    """Find exact k (float) and nearest integer k for a nonce on the grid."""
    if nonce <= 0:
        return None, None
    k_exact = 4.0 * (LP_2_32 - lp(nonce))
    return round(k_exact), k_exact

def cube_face(k):
    """Which cube face (0-indexed) does quarter-step k sit on?
    Each face spans 4 quarter-steps: face 0 = k[0..3], face 1 = k[4..7], etc."""
    return k // 4

def cube_vertex(k):
    """Position within the face (0-3)."""
    return k % 4

def face_mates(k):
    """The other 3 vertices sharing a face with k."""
    base = (k // 4) * 4
    return [base + v for v in range(4) if base + v != k]

# ── CUBE SOLVER ──────────────────────────────────────────────────────────────

class CubeSolver:
    """
    Bitcoin nonce solver using quarter-step phi cube geometry.

    The cube: 8 words = 8 vertices, 3 edges per vertex, edge = 0.25 lp.
    Nonce is at vertex n+3 (cube diagonal from merkle tail).
    The answer shares a face with the merkle root.
    """

    def __init__(self, header_76: bytes, bits: int):
        assert len(header_76) == 76
        self.header_76 = header_76
        self.bits      = bits
        self.target    = bits_to_target(bits)

        block1       = header_76[:64]
        self.H_mid   = _sha256_block(H0, block1)

        self.timestamp = struct.unpack('<I', header_76[68:72])[0]

        # Extract merkle root words (from header bytes 36-68, big-endian display)
        self.merkle_le = header_76[36:68]
        self.merkle_words = struct.unpack('>8I', self.merkle_le[::-1])

    def _make_header(self, nonce):
        return self.header_76 + struct.pack('<I', nonce)

    def _verify(self, nonce):
        h = dsha256(self._make_header(nonce))
        return int.from_bytes(h, 'little') < self.target, h

    # ── Quarter-step predictions ──────────────────────────────────────────────

    def predict_quarter_ladder(self, halvening=0):
        """
        Quarter-step ladder predictions.
        Genesis: k=6 (exponent 6/4 = 1.5).
        Each halvening adds 2 quarter-steps (= 0.5 lp = 1 old half-step).
        Returns dict of k -> nonce_prediction for the face around the target k.
        """
        base_k = 6 + halvening * 2   # genesis at k=6, each halvening +2
        preds = {}
        # Generate the full face (4 vertices) PLUS adjacent faces
        for k in range(max(0, base_k - 4), base_k + 8):
            val = quarter_step(k)
            if 0 < val < MOD32:
                preds[f'q{k}'] = val
        return preds, base_k

    def predict_merkle_face(self):
        """
        Merkle root words mapped onto the quarter-step grid.
        The answer shares a face with the merkle root.
        Returns candidate nonces derived from merkle word positions.
        """
        preds = {}
        for i, w in enumerate(self.merkle_words):
            if w == 0:
                continue
            k_int, k_exact = nearest_k(w)
            if k_int is not None:
                preds[f'mrkl_w{i}_k{k_int}'] = w
                # Also check +3 (cube diagonal) from each merkle word
                k_diag = k_int + 3
                diag_val = quarter_step(k_diag)
                if 0 < diag_val < MOD32:
                    preds[f'mrkl_w{i}_diag{k_diag}'] = diag_val
                # Face mates
                for mate_k in face_mates(k_int):
                    mate_val = quarter_step(mate_k)
                    if 0 < mate_val < MOD32:
                        preds[f'mrkl_w{i}_face{mate_k}'] = mate_val
        return preds

    def predict_H_mid_face(self):
        """H_mid words (nonce-independent state) mapped onto quarter-step grid."""
        preds = {}
        for i, w in enumerate(self.H_mid):
            if w == 0:
                continue
            k_int, _ = nearest_k(w)
            if k_int is not None:
                preds[f'hmid_w{i}_k{k_int}'] = w
                k_diag = k_int + 3
                diag_val = quarter_step(k_diag)
                if 0 < diag_val < MOD32:
                    preds[f'hmid_w{i}_diag{k_diag}'] = diag_val
        return preds

    def predict_timescale(self):
        """Time-scale formula (genesis calibration)."""
        return {
            'ts_100yr': round(YEARS_100_SEC / PHI**PHI_EPS_GENESIS),
            'ts_stamp': round(float(self.timestamp) / PHI**PHI_EPS_GENESIS) if self.timestamp > 0 else 0,
        }

    # ── Solve ─────────────────────────────────────────────────────────────────

    def solve(self, halvening=0, radius=4_500_000, verbose=True):
        t0 = time.time()

        # Collect ALL predictions from every geometric source
        all_preds = {}

        # 1. Quarter-step ladder
        ql_preds, base_k = self.predict_quarter_ladder(halvening)
        all_preds.update(ql_preds)

        # 2. Merkle face geometry
        mf_preds = self.predict_merkle_face()
        all_preds.update(mf_preds)

        # 3. H_mid cube projections
        hm_preds = self.predict_H_mid_face()
        all_preds.update(hm_preds)

        # 4. Time-scale
        ts_preds = self.predict_timescale()
        all_preds.update(ts_preds)

        # 5. Mirror of every prediction
        mirror_preds = {}
        for label, val in list(all_preds.items()):
            m = mirror_nonce(val)
            if 0 < m < MOD32:
                mirror_preds[f'{label}+mir'] = m
        all_preds.update(mirror_preds)

        # Deduplicate
        seen = set()
        unique_preds = {}
        for label, val in all_preds.items():
            if val not in seen:
                seen.add(val)
                unique_preds[label] = val

        if verbose:
            print(f"  Quarter-step base: k={base_k}  (exp={base_k/4:.2f})")
            print(f"  Total predictions: {len(unique_preds)} unique candidates")
            print(f"  Sources: ladder={len(ql_preds)}, merkle_face={len(mf_preds)}, "
                  f"H_mid={len(hm_preds)}, timescale={len(ts_preds)}, mirrors={len(mirror_preds)}")

        # Phase 1: check every prediction directly (O(N) where N = ~200)
        attempts = 0
        for label, val in unique_preds.items():
            attempts += 1
            valid, h = self._verify(val)
            if valid:
                elapsed = time.time() - t0
                return self._result(val, h, attempts, elapsed, unique_preds, label, base_k)

        if verbose:
            print(f"  Phase 1: {attempts} direct checks, no hit")

        # Phase 2: spiral outward from best center
        center = quarter_step(base_k)
        ts_100 = ts_preds.get('ts_100yr', center)
        if abs(ts_100 - center) < radius:
            center = ts_100

        lo = max(0, center - radius)
        hi = min(MOD32 - 1, center + radius)

        if verbose:
            print(f"  Phase 2: spiral from {center}, radius {radius:,}")

        n_lo, n_hi = center, center + 1
        while n_lo >= lo or n_hi <= hi:
            if n_lo >= lo:
                attempts += 1
                valid, h = self._verify(n_lo)
                if valid:
                    elapsed = time.time() - t0
                    return self._result(n_lo, h, attempts, elapsed, unique_preds, 'spiral', base_k)
                n_lo -= 1
            if n_hi <= hi:
                attempts += 1
                valid, h = self._verify(n_hi)
                if valid:
                    elapsed = time.time() - t0
                    return self._result(n_hi, h, attempts, elapsed, unique_preds, 'spiral', base_k)
                n_hi += 1
            if verbose and attempts % 500_000 == 0:
                print(f"    ... {attempts:,} attempts")

        elapsed = time.time() - t0
        return {'found': False, 'attempts': attempts, 'elapsed': elapsed}

    def _result(self, nonce, hash_bytes, attempts, elapsed, preds, source, base_k):
        k_int, k_exact = nearest_k(nonce)
        return {
            'found':       True,
            'nonce':       nonce,
            'nonce_hex':   f'{nonce:#010x}',
            'hash':        hash_bytes[::-1].hex(),
            'attempts':    attempts,
            'elapsed':     elapsed,
            'speedup':     MOD32 / max(attempts, 1),
            'source':      source,
            'lp_nonce':    lp(nonce),
            'k_exact':     k_exact,
            'k_nearest':   k_int,
            'k_base':      base_k,
            'face':        cube_face(k_int) if k_int else None,
            'vertex':      cube_vertex(k_int) if k_int else None,
            'mirror':      mirror_nonce(nonce),
            'n_preds':     len(preds),
        }


# ── Quarter-step ladder display ──────────────────────────────────────────────

def show_quarter_ladder():
    print("=" * 80)
    print("QUARTER-STEP PHI CUBE LADDER")
    print("  2^32 / phi^(k/4)  |  each edge = 0.25 lp  |  cube face = 4 vertices")
    print("=" * 80)
    print(f"  {'k':>4}  {'exp':>6}  {'nonce_val':>14}  {'lp':>8}  {'face':>4}  {'vtx':>3}  note")

    genesis_nonce = 2083236893
    genesis_mirror = mirror_nonce(genesis_nonce)

    for k in range(0, 100):
        exp = k / 4.0
        val = MOD32 / PHI**exp
        if val < 1:
            print(f"  {k:>4}  {exp:>6.2f}  {'< 1':>14}  {'---':>8}  {cube_face(k):>4}  {cube_vertex(k):>3}  ZERO SPACE")
            break

        lp_v = lp(val)
        face = cube_face(k)
        vtx  = cube_vertex(k)

        note = ""
        if k == 6:
            note = f"  <-- GENESIS NONCE (actual: {genesis_nonce}, delta: {abs(round(val) - genesis_nonce):,})"
        elif abs(round(val) - genesis_mirror) < 50_000_000:
            km, _ = nearest_k(genesis_mirror)
            if k == km:
                note = f"  <-- GENESIS MIRROR ({genesis_mirror})"

        # Mark halvening boundaries
        halvening_k = None
        for h in range(10):
            hk = 6 + h * 2
            if k == hk and h > 0:
                note = f"  <-- halvening {h}"

        # Mark face boundaries
        face_marker = " |" if vtx == 0 else "  "

        print(f"  {k:>4}  {exp:>6.2f}  {val:>14.0f}  {lp_v:>8.4f}  {face:>4}  {vtx:>3}{face_marker} {note}")


# ── Merkle root cube analysis ────────────────────────────────────────────────

def analyze_merkle_cube():
    GENESIS_HEADER_HEX = (
        "0100000000000000000000000000000000000000000000000000000000000000"
        "000000003ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa"
        "4b1e5e4a29ab5f49ffff001d1dac2b7c"
    )
    header = bytes.fromhex(GENESIS_HEADER_HEX)
    merkle_le = header[36:68]
    merkle_be = merkle_le[::-1]
    merkle_words = struct.unpack('>8I', merkle_be)

    # Also get the block hash words
    block_hash = dsha256(header)
    hash_words = struct.unpack('>8I', block_hash[::-1])  # big-endian display order

    genesis_nonce = 2083236893
    k_nonce, k_nonce_exact = nearest_k(genesis_nonce)

    print("\n" + "=" * 80)
    print("GENESIS CUBE FACE ANALYSIS")
    print("  8 words = 8 vertices. each vertex on the quarter-step grid.")
    print("  'the answer shares a face'")
    print("=" * 80)

    print(f"\n  NONCE:  {genesis_nonce}  k={k_nonce} (exact {k_nonce_exact:.4f})  face={cube_face(k_nonce)}  vtx={cube_vertex(k_nonce)}")

    print(f"\n  MERKLE ROOT WORDS (8 vertices):")
    print(f"  {'word':>6}  {'value':>12}  {'lp':>8}  {'k':>6}  {'k_exact':>8}  {'face':>4}  {'vtx':>3}  shared_face_w_nonce?")
    nonce_face = cube_face(k_nonce)
    for i, w in enumerate(merkle_words):
        if w == 0:
            print(f"  w[{i}]    {w:>12,}     ---     ---       ---   ---  ---")
            continue
        k_i, k_e = nearest_k(w)
        if k_i is None or k_e is None:
            print(f"  w[{i}]    {w:>12,}     ---     ---       ---   ---  ---")
            continue
        face_i = cube_face(k_i)
        vtx_i  = cube_vertex(k_i)
        shared = "YES" if face_i == nonce_face else ""
        print(f"  w[{i}]    {w:>12,}  {lp(w):>8.4f}  {k_i:>6}  {k_e:>8.4f}  {face_i:>4}  {vtx_i:>3}  {shared}")

    print(f"\n  BLOCK HASH WORDS (8 vertices):")
    print(f"  {'word':>6}  {'value':>12}  {'lp':>8}  {'k':>6}  {'k_exact':>8}  {'face':>4}  {'vtx':>3}  shared_face_w_nonce?")
    for i, w in enumerate(hash_words):
        if w == 0:
            print(f"  w[{i}]    {w:>12,}     ---     ---       ---   ---  ---")
            continue
        k_i, k_e = nearest_k(w)
        if k_i is None or k_e is None:
            print(f"  w[{i}]    {w:>12,}     ---     ---       ---   ---  ---")
            continue
        face_i = cube_face(k_i)
        vtx_i  = cube_vertex(k_i)
        shared = "YES" if face_i == nonce_face else ""
        print(f"  w[{i}]    {w:>12,}  {lp(w):>8.4f}  {k_i:>6}  {k_e:>8.4f}  {face_i:>4}  {vtx_i:>3}  {shared}")

    # Count shared faces
    merkle_faces = set()
    for w in merkle_words:
        k_i, _ = nearest_k(w)
        if k_i is not None:
            merkle_faces.add(cube_face(k_i))

    hash_faces = set()
    for w in hash_words:
        if w > 0:
            k_i, _ = nearest_k(w)
            if k_i is not None:
                hash_faces.add(cube_face(k_i))

    shared = merkle_faces & hash_faces
    print(f"\n  Merkle faces: {sorted(merkle_faces)}")
    print(f"  Hash faces:   {sorted(hash_faces)}")
    print(f"  SHARED faces: {sorted(shared)}")
    print(f"  Nonce face:   {nonce_face}  {'<-- IN SHARED SET' if nonce_face in shared else ''}")


# ── Self-test genesis ────────────────────────────────────────────────────────

def test_genesis():
    GENESIS_HEADER_HEX = (
        "0100000000000000000000000000000000000000000000000000000000000000"
        "000000003ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa"
        "4b1e5e4a29ab5f49ffff001d1dac2b7c"
    )
    header     = bytes.fromhex(GENESIS_HEADER_HEX)
    header_76  = header[:76]
    bits       = struct.unpack('<I', header[72:76])[0]
    true_nonce = struct.unpack('<I', header[76:80])[0]
    k_true, k_exact = nearest_k(true_nonce)

    print("\n" + "=" * 80)
    print("SELF-TEST: GENESIS BLOCK (quarter-step cube solver)")
    print(f"  True nonce: {true_nonce} ({true_nonce:#010x})")
    print(f"  True k: {k_true} (exact: {k_exact:.6f})  face={cube_face(k_true)}  vtx={cube_vertex(k_true)}")
    print(f"  Mirror: {mirror_nonce(true_nonce):#010x}")
    print("=" * 80)

    solver = CubeSolver(header_76, bits)
    result = solver.solve(halvening=0, verbose=True)

    if result['found']:
        print(f"\n  SOLVED")
        print(f"  nonce:     {result['nonce']} ({result['nonce_hex']})")
        print(f"  hash:      {result['hash']}")
        print(f"  source:    {result['source']}")
        print(f"  attempts:  {result['attempts']:,}")
        print(f"  speedup:   {result['speedup']:.0f}x")
        print(f"  k:         {result['k_nearest']} (exact {result['k_exact']:.4f})")
        print(f"  face:      {result['face']}  vtx: {result['vertex']}")
        print(f"  elapsed:   {result['elapsed']:.4f}s")
        if result['nonce'] == true_nonce:
            print(f"\n  GENESIS NONCE RECOVERED EXACTLY.")
    return result


# ── Test against early blocks (fetched from blockstream) ─────────────────────

def test_early_blocks(n_blocks=10):
    """Fetch blocks 1..n and test the quarter-step predictions."""
    print("\n" + "=" * 80)
    print(f"EARLY BLOCK TEST: blocks 1-{n_blocks} (live fetch from blockstream)")
    print("  Testing: does the nonce land on the quarter-step grid?")
    print("=" * 80)

    API = "https://blockstream.info/api"

    results = []
    for height in range(1, n_blocks + 1):
        header = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(f"{API}/block-height/{height}", timeout=15) as r:
                    bh = r.read().decode().strip()
                with urllib.request.urlopen(f"{API}/block/{bh}/header", timeout=15) as r:
                    hdr_hex = r.read().decode().strip()
                header = bytes.fromhex(hdr_hex)
                break
            except Exception as e:
                if attempt < 2:
                    time.sleep(1)
                else:
                    print(f"  block {height}: fetch error: {e}")
        if header is None:
            continue

        nonce  = struct.unpack('<I', header[76:80])[0]
        k_int, k_exact = nearest_k(nonce)
        face   = cube_face(k_int) if k_int else None
        vtx    = cube_vertex(k_int) if k_int else None
        k_frac = k_exact - k_int if k_exact and k_int else None

        # Check if nonce is near a quarter-step grid point
        grid_val = quarter_step(k_int) if k_int else 0
        grid_delta = abs(nonce - grid_val) if k_int else 0
        grid_pct   = grid_delta / nonce * 100 if nonce > 0 else 0

        # Mirror analysis
        mir = mirror_nonce(nonce)
        k_mir, _ = nearest_k(mir)

        actual_hash = dsha256(header)[::-1].hex()

        results.append({
            'height': height, 'nonce': nonce, 'k': k_int, 'k_exact': k_exact,
            'face': face, 'vtx': vtx, 'grid_delta': grid_delta, 'grid_pct': grid_pct,
            'k_frac': k_frac, 'mirror': mir, 'k_mirror': k_mir, 'hash': actual_hash,
        })

        if k_frac is not None:
            on_grid = "ON GRID" if abs(k_frac) < 0.15 else f"off by {k_frac:+.3f}"
        else:
            on_grid = "N/A"
        print(f"  block {height:>3}  nonce={nonce:>12,}  k={k_int:>3} ({k_exact:>7.3f})  "
              f"face={face}  vtx={vtx}  {on_grid}  delta={grid_pct:.2f}%")

    # Stats
    valid = [r for r in results if r['k_frac'] is not None]
    on_count = sum(1 for r in valid if abs(r['k_frac']) < 0.15)
    print(f"\n  ON GRID (|k_frac| < 0.15): {on_count}/{len(valid)}")
    if valid:
        mean_frac = sum(abs(r['k_frac']) for r in valid) / len(valid)
        print(f"  Mean |k_frac|: {mean_frac:.4f}")

        # Face distribution
        face_counts = {}
        for r in valid:
            f = r['face']
            face_counts[f] = face_counts.get(f, 0) + 1
        print(f"  Face distribution: {dict(sorted(face_counts.items()))}")

    return results


# ── Run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("PHI CUBE SOLVER -- quarter-step geometry")
    print("nos3bl33d")
    print("x^2 = x + 1.")
    print()

    show_quarter_ladder()
    analyze_merkle_cube()
    test_genesis()
    test_early_blocks(20)
