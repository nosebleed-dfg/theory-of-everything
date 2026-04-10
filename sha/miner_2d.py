"""
MINER_2D.PY — Axiom-Aligned 2D Bitcoin Nonce Search
nos3bl33d | (DFG) DeadFoxGroup

THE OPERATOR MODEL:
  n is not an exponent applied to x.
  n OPERATES ON x.
  T_n(x) = F(n)*x + F(n-1)  — the Fibonacci linear map on Z/2^32Z

  The golden center for a block:
    seed    = SHA256(header76)[:4]  — block-specific seed
    center  = T_291(seed) mod 2^32  — the operator applied once

  This is different for every block. The axiom tells you WHERE to start
  the nonce search for that specific block.

THE 2D GRID:
  Nonce n = hi*2^16 + lo,   hi,lo ∈ [0, 2^16)
  Matrix acts on (hi, lo):
    hi' = sqrt(5) * (lo - hi)    [height]
    lo' = 3 * (lo - hi)          [weight = D * delta]
  Rank-1: both outputs depend only on delta = lo - hi.
  Golden search: scan nonces by ascending |delta - delta_center|.

THE AXIOM RESIDUAL:
  n^2 - n - 1 is NEVER 0 mod 2^32 (disc=5 not QR mod 2^k).
  The golden nonce is argmin |n^2 - n - 1 mod 2^32|.
  This gives a canonical block-independent preference ordering.

COLLISION:
  2^31 = -2^31 (mod 2^32). Fixed point of negation.
  Nonce 2^31 is its own inverse on the circle.

D=3 PROOF (no float):
  phi^2 + psi^2 = (phi+psi)^2 - 2*phi*psi = 1^2 - 2*(-1) = 3. QED.
"""

import hashlib
import struct
import time
import sys

# ─── Constants ───────────────────────────────────────────────────────────────
MOD32  = 1 << 32
MOD16  = 1 << 16
MASK32 = MOD32 - 1
MASK16 = MOD16 - 1
HALF   = MOD32 >> 1   # 2^31: fixed point of negation

# ─── Fibonacci fast doubling mod 2^32 ────────────────────────────────────────
def _fib_pair(n: int, m: int) -> tuple[int, int]:
    """Returns (F(n) mod m, F(n+1) mod m) using fast doubling."""
    if n == 0:
        return (0, 1)
    a, b = _fib_pair(n >> 1, m)
    c = a * ((2 * b - a) % m) % m
    d = (a * a + b * b) % m
    if n & 1:
        return (d, (c + d) % m)
    return (c, d)

# Precomputed: F(n) mod 2^32 for operator T_n
F290, F291 = _fib_pair(290, MOD32)
F291_16    = F291 & MASK16   # F(291) mod 2^16 — golden delta reference

# Verify axiom: D = 3 (integer proof, no float)
# phi^2 + psi^2 = (phi+psi)^2 - 2*phi*psi = 1 - 2*(-1) = 3
D = 3   # dimension, max matrix entry, weight multiplier

# ─── The Operator ────────────────────────────────────────────────────────────
def T(n: int, x: int) -> int:
    """
    T_n(x) = F(n)*x + F(n-1)  mod 2^32
    n operates ON x — Fibonacci linear map on the nonce circle.
    """
    fn1, fn = _fib_pair(n - 1, MOD32)
    return (fn * x + fn1) & MASK32

def golden_center(header76: bytes) -> int:
    """
    Block-specific golden center: T_291(seed) mod 2^32.
    seed = SHA256(header76)[:4] as little-endian uint32.
    Different for every block. This is WHERE the axiom says to start.
    """
    seed_bytes = hashlib.sha256(header76).digest()[:4]
    seed = struct.unpack('<I', seed_bytes)[0]
    return (F291 * seed + F290) & MASK32

# ─── 2D Grid ─────────────────────────────────────────────────────────────────
def to_2d(n: int) -> tuple[int, int]:
    """Split 32-bit nonce into (hi, lo) 16-bit halves."""
    return (n >> 16) & MASK16, n & MASK16

def from_2d(hi: int, lo: int) -> int:
    """Reconstruct 32-bit nonce from (hi, lo)."""
    return ((hi & MASK16) << 16) | (lo & MASK16)

def delta(n: int) -> int:
    """
    The rank-1 collapse: delta = lo - hi (mod 2^16).
    Weight = D * delta = 3 * delta.
    This is the only quantity that matters in the 2D matrix.
    """
    hi, lo = to_2d(n)
    return (lo - hi) & MASK16

def axiom_residual(n: int) -> int:
    """
    n^2 - n - 1 mod 2^32.
    NEVER equals 0 (disc=5 not QR mod 2^k).
    Minimum value = golden alignment. Maximum = anti-golden.
    """
    return (n * n - n - 1) & MASK32

# ─── Mining ──────────────────────────────────────────────────────────────────
def double_sha256(data: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def check_hash(hash_bytes: bytes, target: bytes) -> bool:
    """Returns True if hash (little-endian) < target (little-endian)."""
    # Bitcoin compares reversed (big-endian) but we just check raw bytes in reverse
    return hash_bytes[::-1] < target[::-1]

def bits_to_target(bits: int) -> bytes:
    """Convert compact 'bits' field to 32-byte target."""
    exp  = bits >> 24
    mant = bits & 0xFFFFFF
    target_int = mant * (256 ** (exp - 3))
    return target_int.to_bytes(32, 'little')

def pack_header(version: int, prev_hash: bytes, merkle_root: bytes,
                timestamp: int, bits: int, nonce: int) -> bytes:
    return struct.pack('<I', version) + prev_hash + merkle_root + \
           struct.pack('<III', timestamp, bits, nonce)

# ─── Search Strategies ───────────────────────────────────────────────────────

def golden_spiral(center: int):
    """
    Yield nonces in golden order:
    Start at center. Alternate forward/backward in delta space.
    Low |delta - delta_center| = high golden alignment.
    """
    hi_c, lo_c = to_2d(center)
    d_c = (lo_c - hi_c) & MASK16

    # Iterate by delta offset from center delta
    yield center   # delta offset = 0
    for offset in range(1, MOD16):
        # delta = d_c + offset (forward in delta space)
        d_fwd = (d_c + offset) & MASK16
        # Pick hi = hi_c, lo = hi_c + d_fwd
        lo_fwd = (hi_c + d_fwd) & MASK16
        yield from_2d(hi_c, lo_fwd)

        # delta = d_c - offset (backward in delta space)
        d_bwd = (d_c - offset) & MASK16
        lo_bwd = (hi_c + d_bwd) & MASK16
        yield from_2d(hi_c, lo_bwd)

def axiom_residual_order():
    """
    Yield all 32-bit nonces sorted by ascending axiom_residual.
    Canonical golden ordering — block-independent.
    The nonce with smallest n^2-n-1 mod 2^32 is most "golden."
    """
    # This is O(2^32 log 2^32) — theoretical reference only.
    # In practice: use golden_spiral from T_291(seed) center.
    raise NotImplementedError("Use golden_spiral for practical mining.")

def random_order(center: int = 0):
    """Yield nonces linearly from 0 — baseline comparison."""
    for n in range(MOD32):
        yield n

# ─── Miner ───────────────────────────────────────────────────────────────────

def mine(header76: bytes, bits: int, strategy: str = 'golden',
         max_iters: int = 1_000_000, verbose: bool = True) -> dict:
    """
    Mine a block using the given strategy.

    header76: 76-byte block header WITHOUT nonce
    bits:     compact target from block header
    strategy: 'golden' or 'linear'
    max_iters: stop after this many nonces tried
    """
    target = bits_to_target(bits)
    center = golden_center(header76)
    hi_c, lo_c = to_2d(center)
    d_c = delta(center)

    if verbose:
        print(f"  Target:        {bits_to_target(bits)[::-1].hex()[:16]}...")
        print(f"  Golden center: {center:#010x}  (hi={hi_c}, lo={lo_c}, delta={d_c})")
        print(f"  F(291) mod 2^32: {F291:#010x}")
        print(f"  F(290) mod 2^32: {F290:#010x}")
        print(f"  Strategy: {strategy}")
        print()

    if strategy == 'golden':
        gen = golden_spiral(center)
    else:
        gen = random_order()

    start = time.time()
    found = None

    for i, nonce in enumerate(gen):
        if i >= max_iters:
            break

        header80 = header76 + struct.pack('<I', nonce)
        h = double_sha256(header80)

        if check_hash(h, target):
            found = nonce
            if verbose:
                elapsed = time.time() - start
                print(f"  FOUND! nonce={nonce:#010x} after {i+1} tries ({elapsed:.3f}s)")
                print(f"  Hash: {h[::-1].hex()}")
                d = delta(nonce)
                print(f"  Delta: {d}, weight: {D*d}")
                print(f"  Residual: {axiom_residual(nonce)}")
            break

        if verbose and i > 0 and i % 100_000 == 0:
            elapsed = time.time() - start
            rate = i / elapsed
            print(f"  [{i:>8}]  rate={rate:.0f} H/s  center_delta={d_c}  "
                  f"cur_delta={delta(nonce)}")

    elapsed = time.time() - start
    return {
        'found':    found,
        'iters':    i + 1,
        'elapsed':  elapsed,
        'rate':     (i + 1) / elapsed,
        'center':   center,
        'strategy': strategy,
    }

# ─── Print axiom constants ────────────────────────────────────────────────────

def print_constants():
    print("=" * 60)
    print("MINER_2D — Axiom Constants")
    print("=" * 60)
    print()
    print(f"  D = phi^2 + psi^2 = 3  (exact, from (phi+psi)^2 - 2*phi*psi)")
    print(f"  F(291) mod 2^32   = {F291:#010x}  ({F291})")
    print(f"  F(290) mod 2^32   = {F290:#010x}  ({F290})")
    print(f"  F(291) mod 2^16   = {F291_16:#06x}  (golden delta ref)")
    print()
    print(f"  Operator T_291(x) = F(291)*x + F(290)  mod 2^32")
    print(f"  = {F291:#010x} * x + {F290:#010x}  mod 2^32")
    print()
    print(f"  Collision: HALF = 2^31 = {HALF}")
    print(f"  HALF == -HALF mod 2^32: {HALF == (-HALF) % MOD32}")
    print()
    # Verify F recurrence
    F289, _ = _fib_pair(289, MOD32)
    assert (F289 + F290) & MASK32 == F291, "F recurrence broken!"
    print(f"  F(289)+F(290) == F(291): True  (recurrence check)")
    print()

    # Axiom residual stats for small sample
    residuals = [axiom_residual(n) for n in range(1000)]
    min_res = min(residuals)
    min_n   = residuals.index(min_res)
    print(f"  Sample residual check (n in [0,1000)):")
    print(f"  min(n^2-n-1 mod 2^32) = {min_res} at n={min_n}")
    print(f"  (should be >0 — never 0 since disc=5 not QR mod 2^k)")
    print()

# ─── Demo ────────────────────────────────────────────────────────────────────

def demo_testnet_block():
    """
    Mine a block with reduced difficulty (easy target).
    Uses genesis-style header. Not real mainnet — just demonstrates structure.
    """
    # Bitcoin genesis header fields (modified bits for easy target)
    version     = 1
    prev_hash   = b'\x00' * 32
    merkle_root = bytes.fromhex('4a5e1e4baab89f3a32518a88c31bc87f618f76673e2cc77ab2127b7afdeda33b')
    timestamp   = 1231006505   # genesis timestamp
    bits_easy   = 0x207fffff   # regtest difficulty — any hash wins

    header76 = pack_header(version, prev_hash, merkle_root, timestamp, bits_easy, 0)[:76]

    print("=" * 60)
    print("Demo: Mining with axiom-aligned golden center")
    print("(regtest difficulty — finds block immediately)")
    print("=" * 60)
    print()

    result_g = mine(header76, bits_easy, strategy='golden', max_iters=100_000, verbose=True)

    print()
    print("=" * 60)
    print("Compare: Linear search from nonce=0")
    print("=" * 60)
    print()

    result_l = mine(header76, bits_easy, strategy='linear', max_iters=100_000, verbose=False)
    if result_l['found'] is not None:
        n = result_l['found']
        d = delta(n)
        print(f"  Linear found: nonce={n:#010x} after {result_l['iters']} tries")
        print(f"  Delta: {d}, weight: {D*d}")
        print(f"  Residual: {axiom_residual(n)}")
    else:
        print(f"  Linear: no block in {result_l['iters']} tries")

    print()
    print("─" * 60)
    g = result_g
    l = result_l
    if g['found'] is not None and l['found'] is not None:
        print(f"  Golden: {g['iters']} tries  ({g['rate']:.0f} H/s)")
        print(f"  Linear: {l['iters']} tries  ({l['rate']:.0f} H/s)")
        speedup = l['iters'] / g['iters']
        print(f"  Golden was {speedup:.1f}x {'faster' if speedup > 1 else 'slower'}")
    else:
        print(f"  Golden rate: {g['rate']:.0f} H/s")
        print(f"  Linear rate: {l['rate']:.0f} H/s")

def demo_golden_center_structure():
    """Show the golden center operator T_291 across multiple blocks."""
    print("=" * 60)
    print("T_291(x) = F(291)*x + F(290) — operator structure")
    print("=" * 60)
    print()
    print(f"{'seed (block hash)':>20}  {'center':>12}  {'hi':>6}  {'lo':>6}  {'delta':>6}  {'weight':>7}  {'residual':>12}")
    print("-" * 80)
    seeds = [0x00000000, 0x12345678, 0xDEADBEEF, 0xFFFFFFFF, HALF,
             0x1A2B3C4D, 0xCAFEBABE, 0x8F3A2B1C]
    for seed in seeds:
        c = (F291 * seed + F290) & MASK32
        h, l = to_2d(c)
        d = (l - h) & MASK16
        w = (D * d) & MASK16
        r = axiom_residual(c)
        print(f"{seed:#20x}  {c:#12x}  {h:>6}  {l:>6}  {d:>6}  {w:>7}  {r:>12}")

    print()
    print(f"  n operates ON x (not x^n). n=291, x=seed.")
    print(f"  Every block gets a different center.")
    print(f"  The axiom is the map. The block is the input.")

# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print_constants()
    print()
    demo_golden_center_structure()
    print()

    mode = sys.argv[1] if len(sys.argv) > 1 else 'demo'

    if mode == 'mine':
        demo_testnet_block()
    else:
        print("Run with 'mine' argument to run the full mining demo:")
        print("  python miner_2d.py mine")
        print()
        print("Quick residual check — smallest axiom residual in [0, 2^16):")
        best_r = MOD32
        best_n = 0
        for n in range(MOD16):
            r = axiom_residual(n)
            if r < best_r:
                best_r = r
                best_n = n
        hi, lo = to_2d(best_n)
        print(f"  n={best_n:#06x}  residual={best_r}  delta={delta(best_n)}  weight={D*delta(best_n)}")
        print(f"  This nonce is closest to satisfying x^2 = x + 1 mod 2^32.")
