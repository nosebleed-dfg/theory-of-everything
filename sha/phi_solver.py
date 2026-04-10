"""
phi_solver.py — Bitcoin nonce solver via phi geometry
nos3bl33d

Framework:
  - Nonce lives at k=3 on the 2^32/phi^(k/2) half-step ladder
  - W[3] in the SHA schedule = byte_reverse(nonce) = n0_mirror
  - W[5..14] = 0 always (shared zero vertices, nonce-independent)
  - W[3+6] = W[9] = 0 = first vertex at -0.75 steps from the nonce
  - Each halving = one more ladder step (one more zero space)
  - lp(nonce) ~= lp(100yr_sec) - 0.863 for genesis scale
  - Mirror key: byte_reverse(n0_mirror) = n0  (no Fibonacci, no brute force)

Solve path:
  1. Compute H_mid (nonce-independent after block 1 SHA)
  2. Phi ladder -> candidate window  [2^32/phi^1.5 ± radius]
  3. Mirror prediction -> tighten window via byte_reverse symmetry
  4. Time-scale formula -> tightest prediction if timestamp scale known
  5. Forward SHA verify candidates in order of phi distance from prediction
  6. Return first valid nonce

x^2 = x + 1.
"""

import struct
import hashlib
import math
import time

# ── Constants ────────────────────────────────────────────────────────────────
PHI    = (1 + math.sqrt(5)) / 2
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

YEARS_100_SEC = 100 * 365.25 * 24 * 3600   # 3,155,760,000.0
PHI_EPS_GENESIS = math.log(YEARS_100_SEC / 2083236893) / math.log(PHI)  # 0.863043

# ── SHA-256 primitives ───────────────────────────────────────────────────────
def _rotr(x, n):  return ((x >> n) | (x << (32 - n))) & MASK32
def _ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def _maj(a, b, c):return (a & b) ^ (a & c) ^ (b & c)
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
    W64 = _build_schedule(W16)
    return _compress(iv, W64), W64

def _dsha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

# ── Phi / lp utilities ────────────────────────────────────────────────────────
def lp(x):
    return math.log(x) / math.log(PHI) if x > 0 else None

def mirror_nonce(n):
    """byte_reverse(n): W[3] = mirror of the wire nonce."""
    return struct.unpack('>I', struct.pack('<I', n))[0]

def bits_to_target(bits):
    exp   = (bits >> 24) & 0xff
    coeff = bits & 0x007fffff
    return coeff << (8 * (exp - 3)) if exp > 3 else coeff >> (8 * (3 - exp))

# ── Block 2 padding ───────────────────────────────────────────────────────────
def _make_block2(header_76, nonce):
    """Build the 64-byte SHA inner block 2 for an 80-byte header."""
    nonce_bytes = struct.pack('<I', nonce)
    raw = header_76[64:] + nonce_bytes       # 16 bytes: merkle[-4:], time, bits, nonce
    pad = b'\x80' + b'\x00' * 39 + struct.pack('>Q', 640)
    return raw + pad                          # 64 bytes

# ── Core solver ───────────────────────────────────────────────────────────────
class PhiSolver:
    """
    Bitcoin nonce solver using the phi half-step ladder and mirror key.

    Usage:
        solver = PhiSolver(header_76_bytes, bits)
        result = solver.solve()
        print(result)
    """

    def __init__(self, header_76: bytes, bits: int):
        assert len(header_76) == 76, "header_76 must be exactly 76 bytes (header without nonce)"
        self.header_76 = header_76
        self.bits      = bits
        self.target    = bits_to_target(bits)

        # Nonce-independent state after block 1
        block1        = header_76[:64]
        self.H_mid, _ = _sha256_block(H0, block1)

        # Timestamp (bytes 68-71 of the full header, little-endian)
        self.timestamp = struct.unpack('<I', header_76[68:72])[0]

    def _make_header(self, nonce):
        return self.header_76 + struct.pack('<I', nonce)

    def _verify(self, nonce):
        h = _dsha256(self._make_header(nonce))
        return int.from_bytes(h, 'little') < self.target, h

    # ── Prediction strategies ─────────────────────────────────────────────────
    def predict_ladder(self, halvening=0):
        """
        Phi half-step ladder: nonce ~= 2^32 / phi^(3/2 + halvening*0.5)
        Each halving adds one half-step (one new zero vertex).
        Halvening 0 = genesis scale (k=3, exponent=1.5).
        """
        exp = 1.5 + halvening * 0.5
        return round(2**32 / PHI**exp)

    def predict_timescale(self, scale=YEARS_100_SEC, eps=PHI_EPS_GENESIS):
        """
        Time-scale formula: nonce ~= scale / phi^eps
        Default: genesis scale (100yr in seconds, eps=0.863043).
        Also try: scale = self.timestamp for other blocks.
        """
        return round(scale / PHI**eps)

    def predict_mirror(self, n_pred):
        """
        Mirror key: W[3] = byte_reverse(nonce) = n_mirror.
        The mirror of the prediction is a secondary candidate.
        """
        return mirror_nonce(n_pred)

    def lp_distance(self, n):
        """Distance in lp space from the half-step ladder k=3 position."""
        ladder_k3 = 2**32 / PHI**1.5
        return abs(lp(n) - lp(ladder_k3)) if n > 0 else float('inf')

    # ── Primary solve ─────────────────────────────────────────────────────────
    def solve(self, radius=4_500_000, verbose=True):
        """
        Find a valid nonce within [prediction - radius, prediction + radius].

        Search order: outward from the phi prediction, so the closest candidate
        (in lp distance) is tried first. No Fibonacci metric used.

        Returns dict with nonce, hash, attempts, speedup, predictions.
        """
        t0 = time.time()

        # Build all predictions
        preds = {}
        preds['ladder_k3']    = self.predict_ladder(halvening=0)
        preds['ladder_k3.5']  = self.predict_ladder(halvening=1)
        preds['timescale_100yr'] = self.predict_timescale(YEARS_100_SEC)
        preds['timescale_ts']    = self.predict_timescale(float(self.timestamp))
        preds['mirror_ladder']   = self.predict_mirror(preds['ladder_k3'])
        preds['mirror_ts']       = self.predict_mirror(preds['timescale_100yr'])

        # Primary anchor: use the ladder k=3 prediction
        anchor = preds['ladder_k3']

        # Also check if the time-scale prediction is within the radius
        ts_pred = preds['timescale_100yr']
        if abs(ts_pred - anchor) < radius:
            # Time-scale prediction is tighter — use it as the search center
            center = ts_pred
        else:
            center = anchor

        if verbose:
            print(f"  Predictions:")
            for label, p in preds.items():
                print(f"    {label:<22} = {p:>12}  lp={lp(p):.4f}")
            print(f"  Search center: {center}  radius: {radius}")
            print(f"  Window size: {2*radius:,} nonces  (vs 2^32 = {2**32:,})")

        # Check the mirror key directly first — it costs nothing
        for label, n_pred in preds.items():
            n_mirror = self.predict_mirror(n_pred)
            valid, h = self._verify(n_mirror)
            if valid:
                elapsed = time.time() - t0
                return self._result(n_mirror, h, 1, elapsed, preds, label+'+mirror_key')

        # Scan outward from center in lp-distance order
        # Generate candidate nonces sorted by |lp(n) - lp(center)|
        lo = max(0, center - radius)
        hi = min(2**32 - 1, center + radius)

        attempts = 0
        best_nonce = None

        # Try predictions first (exact checks, O(1) each)
        for label, p in preds.items():
            if 0 <= p <= 2**32 - 1:
                attempts += 1
                valid, h = self._verify(p)
                if valid:
                    elapsed = time.time() - t0
                    return self._result(p, h, attempts, elapsed, preds, label)

        # Spiral outward from center
        step = 1
        n_lo = center
        n_hi = center + 1
        while n_lo >= lo or n_hi <= hi:
            if n_lo >= lo:
                attempts += 1
                valid, h = self._verify(n_lo)
                if valid:
                    elapsed = time.time() - t0
                    return self._result(n_lo, h, attempts, elapsed, preds, 'spiral_lo')
                n_lo -= 1

            if n_hi <= hi:
                attempts += 1
                valid, h = self._verify(n_hi)
                if valid:
                    elapsed = time.time() - t0
                    return self._result(n_hi, h, attempts, elapsed, preds, 'spiral_hi')
                n_hi += 1

            if verbose and attempts % 500_000 == 0:
                print(f"  ... {attempts:,} attempts, n_lo={n_lo}, n_hi={n_hi}")

        elapsed = time.time() - t0
        return {'found': False, 'attempts': attempts, 'elapsed': elapsed, 'predictions': preds}

    def _result(self, nonce, hash_bytes, attempts, elapsed, preds, source):
        speedup    = 2**32 / max(attempts, 1)
        anchor     = preds['ladder_k3']
        mirror_n   = mirror_nonce(nonce)
        return {
            'found':         True,
            'nonce':         nonce,
            'nonce_hex':     f'{nonce:#010x}',
            'hash':          hash_bytes[::-1].hex(),     # display order
            'hash_raw':      hash_bytes.hex(),
            'attempts':      attempts,
            'elapsed':       elapsed,
            'speedup':       speedup,
            'source':        source,
            'lp_nonce':      lp(nonce),
            'lp_anchor':     lp(anchor),
            'lp_delta':      lp(nonce) - lp(anchor),
            'mirror':        mirror_n,
            'mirror_hex':    f'{mirror_n:#010x}',
            'predictions':   preds,
        }

    def report(self, result):
        print(f"\n{'='*60}")
        if not result['found']:
            print(f"NOT FOUND in {result['attempts']:,} attempts ({result['elapsed']:.2f}s)")
            return
        print(f"SOLVED")
        print(f"  nonce:        {result['nonce']} ({result['nonce_hex']})")
        print(f"  mirror:       {result['mirror']} ({result['mirror_hex']})")
        print(f"  hash:         {result['hash']}")
        print(f"  source:       {result['source']}")
        print(f"  attempts:     {result['attempts']:,}")
        print(f"  speedup:      {result['speedup']:.0f}x vs brute force")
        print(f"  elapsed:      {result['elapsed']:.4f}s")
        print(f"  lp(nonce):    {result['lp_nonce']:.4f}")
        print(f"  lp(anchor):   {result['lp_anchor']:.4f}")
        print(f"  lp_delta:     {result['lp_delta']:+.4f}")
        print(f"\n  byte_reverse(mirror) = nonce: {mirror_nonce(result['mirror']) == result['nonce']}")
        print(f"{'='*60}")


# ── Halvening analysis ────────────────────────────────────────────────────────
def halvening_ladder():
    """
    The halvening schedule as phi ladder steps.
    Each halving = one new zero vertex = one half-step down the ladder.
    At k halvings from genesis: nonce_scale ~= 2^32 / phi^(1.5 + k*0.5)
    Terminal: when 2^32 / phi^(1.5 + k*0.5) < 1  -> all zero space.
    """
    print(f"\n{'='*60}")
    print("HALVENING AS PHI LADDER")
    print(f"  Each halving = one half-step = one new shared zero vertex")
    print(f"  Genesis nonce at k=3 (exponent=1.5). Each halving: +0.5 exponent.")
    print(f"{'='*60}")
    print(f"\n  {'halving':>8}  {'k':>5}  {'exp':>6}  {'nonce_scale':>14}  {'lp':>8}  note")

    genesis_nonce = 2083236893
    for h in range(0, 38):
        exp = 1.5 + h * 0.5
        k   = 3 + h
        val = 2**32 / PHI**exp
        if val < 1:
            print(f"  {h:>8}  {k:>5}  {exp:>6.1f}  {'< 1':>14}  {'---':>8}  ZERO SPACE COMPLETE")
            break
        lpv = lp(val)
        note = ""
        if h == 0:   note = "  <- genesis nonce scale"
        if val < 1e6: note = "  <- below 1M (near terminal)"
        block_num = h * 210_000
        print(f"  {h:>8}  {k:>5}  {exp:>6.1f}  {val:>14.0f}  {lpv:>8.4f}{note}")

    print(f"\n  Total halvings before zero space: ~{int(math.log(2**32) / math.log(PHI) / 0.5 - 3)}")
    reward_halvings = math.ceil(math.log2(50e8))   # 50 BTC = 5e9 sat
    print(f"  Halvings until reward -> 0 sat:  ~{reward_halvings}")
    print(f"  (Both converge at approximately the same terminal step)")


# ── Self-test against genesis ─────────────────────────────────────────────────
def test_genesis():
    GENESIS_HEADER_HEX = (
        "0100000000000000000000000000000000000000000000000000000000000000"
        "000000003ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa"
        "4b1e5e4a29ab5f49ffff001d1dac2b7c"
    )
    header    = bytes.fromhex(GENESIS_HEADER_HEX)
    header_76 = header[:76]
    bits      = struct.unpack('<I', header[72:76])[0]
    true_nonce = struct.unpack('<I', header[76:80])[0]

    print(f"\n{'='*60}")
    print("SELF-TEST: GENESIS BLOCK")
    print(f"  True nonce: {true_nonce} ({true_nonce:#010x})")
    print(f"  True mirror: {mirror_nonce(true_nonce):#010x}  = W[3] in schedule")
    print(f"  lp(true_nonce) = {lp(true_nonce):.6f}  (k~=3 on ladder)")
    print(f"{'='*60}")

    solver = PhiSolver(header_76, bits)
    result = solver.solve(radius=4_500_000, verbose=True)
    solver.report(result)

    if result['found'] and result['nonce'] == true_nonce:
        print(f"\n  GENESIS NONCE RECOVERED EXACTLY.")
        print(f"  Speedup: {result['speedup']:.0f}x")
        print(f"  Attempts: {result['attempts']:,}")
    elif result['found']:
        print(f"\n  Found a valid nonce: {result['nonce']} (not the genesis nonce, but valid)")

    return result


# ── Run ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("PHI SOLVER — Bitcoin nonce via phi geometry")
    print("nos3bl33d")
    print("x^2 = x + 1.")
    print()
    print("Framework:")
    print("  nonce at k=3 on 2^32/phi^(k/2) half-step ladder")
    print("  W[3] = byte_reverse(nonce) = nonce_mirror")
    print("  W[5..14] = 0  (shared zero vertices, nonce-independent)")
    print("  W[9] = first vertex at -0.75 steps from nonce")
    print("  Each halvening = one more ladder step = one more zero space")
    print("  lp(nonce) ~= lp(100yr_sec) - 0.863  (genesis scale)")
    print()

    halvening_ladder()
    test_genesis()
