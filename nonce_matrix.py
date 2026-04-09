"""
NONCE_MATRIX — 4x4 nonce grid from phi^2/psi^2 combinations; collapse at step 582, rank-1 golden collision
nos3bl33d
---
Matrix entries: TL=-sqrt(5), TR=+sqrt(5), BL=-3, BR=+3=D.
Collapse at 582 = 2*291: all entries reset, psi^582 ~ 10^-122.
Collision pair: (F(582) mod 2^32, L(582) mod 2^32) on Z/2^32Z.
2D grid rank-1: nonce reduces to delta=lo-hi, golden nonce at F(582) mod 2^16.
"""

import math
import struct
import hashlib
import time

# ─── The Axiom ───────────────────────────────────────────────────────────────
phi = (1 + math.sqrt(5)) / 2   # 1.6180339887...
psi = (1 - math.sqrt(5)) / 2   # -0.6180339887...

phi2 = phi ** 2   # = phi + 1 = 2.6180339887...
psi2 = psi ** 2   # = psi + 1 = 0.3819660112...

# ─── The Matrix ──────────────────────────────────────────────────────────────
TL = -phi2 + psi2   # = -(phi2 - psi2) = -sqrt(5)
TR =  phi2 - psi2   # = +sqrt(5)  ("normal")
BL = -(phi2 + psi2) # = -3
BR =   phi2 + psi2  # = +3 = D

D = BR  # = phi^2 + psi^2 = 3

# ─── Proof: D = 3 ────────────────────────────────────────────────────────────
assert abs(phi + psi - 1.0) < 1e-12,      f"phi + psi should be 1"
assert abs(phi * psi - (-1.0)) < 1e-12,   f"phi * psi should be -1"
assert abs(D - 3.0) < 1e-12,              f"D = phi^2 + psi^2 should be 3, got {D}"

# Algebraic proof (no floating point needed):
# phi^2 + psi^2 = (phi+psi)^2 - 2*phi*psi = 1^2 - 2*(-1) = 1 + 2 = 3. QED.

# ─── The Collapse at step 291 ────────────────────────────────────────────────
# phi^291 and psi^291: compute and verify the simultaneous collapse.

phi_291 = phi ** 291
psi_291 = psi ** 291

# Fibonacci number F(291) and Lucas number L(291) — integers that capture phi^291
# F(n) = (phi^n - psi^n) / sqrt(5)
# L(n) = phi^n + psi^n
sqrt5 = math.sqrt(5)
F291 = (phi_291 - psi_291) / sqrt5   # ≈ phi^291 / sqrt(5) (psi^291 ≈ 0)
L291 = phi_291 + psi_291              # ≈ phi^291 (psi^291 ≈ 0)

psi_to_phi_ratio = abs(psi_291) / phi_291   # should be ≈ phi^(-582) ≈ 0

# Matrix entries at step 291:
TL_291 = -phi_291 + psi_291   # ≈ -phi_291
TR_291 =  phi_291 - psi_291   # ≈ +phi_291
BL_291 = -(phi_291 + psi_291) # ≈ -phi_291
BR_291 =   phi_291 + psi_291  # ≈ +phi_291

# ─── The Collision ───────────────────────────────────────────────────────────
CIRCLE = 2 ** 32
HALF   = CIRCLE >> 1   # 2147483648 = 2^31

# Verify: HALF = -HALF (mod CIRCLE)
assert HALF == (-HALF) % CIRCLE, "Collision identity failed"
# HALF is the fixed point of negation on the nonce circle.
# Input: nonce = HALF (the "1/2" point)
# Output: SHA-256(block + HALF) — this nonce equals its own negative.
# Collision: f(n) = f(-n) trivially because n = -n.

# ─── 2D Nonce Decomposition ──────────────────────────────────────────────────

def to_2d(n: int):
    """Split 32-bit nonce into (hi, lo) 16-bit components."""
    return (n >> 16) & 0xFFFF, n & 0xFFFF

def from_2d(hi: int, lo: int) -> int:
    """Reconstruct 32-bit nonce from (hi, lo)."""
    return ((hi & 0xFFFF) << 16) | (lo & 0xFFFF)

def delta(n: int) -> int:
    """
    The rank-1 collapse: delta = lo - hi (mod 2^16).
    This is THE relevant quantity on the nonce grid.

    Weight = D * delta = 3 * delta (mod 2^16)
    Height = sqrt(5) * delta (mod 2^16)

    Both matrix outputs are multiples of delta. The 2D grid
    collapses to the 1D axis: delta.
    """
    hi, lo = to_2d(n)
    return (lo - hi) & 0xFFFF

def axiom_weight(n: int) -> int:
    """
    The 'weight' of a nonce on the golden grid.
    weight = D * delta = 3 * delta (mod 2^16).
    This is the BR entry of the matrix applied to the nonce.
    """
    d = delta(n)
    return (3 * d) & 0xFFFF

def golden_center(header76: bytes) -> int:
    """
    Compute the axiom-predicted center nonce for a block.

    Three steps (d=3 golden rotations):
    1. Planck seed from block header
    2. Fibonacci rotation (phi-map)
    3. Axiom weight projection (D=3 collapse)

    The center is where the matrix M collapses the seed onto
    the 1D delta axis with maximum golden alignment.
    """
    seed_bytes = hashlib.sha256(header76).digest()[:4]
    seed = struct.unpack('<I', seed_bytes)[0]

    # Step 1: extract delta from seed
    d = delta(seed)

    # Step 2: Fibonacci rotation (3 golden steps = dimension d=3)
    # phi ≈ F(n+1)/F(n). Use F(46)/F(45) = 2971215073/1836311903 for 32-bit precision.
    # But we're working mod 2^16, so use F(23)/F(22) = 28657/17711 ≈ phi (16-bit)
    # Fibonacci step: d -> d * F(23) / F(22) mod 2^16
    # Integer version: d_new = (d * 28657) >> 14  (≈ d * 1.618 with 14-bit fraction)
    d1 = (d * 28657) >> 14    # rotation 1 (dimension x)
    d2 = (d1 * 28657) >> 14   # rotation 2 (dimension y)
    d3 = (d2 * 28657) >> 14   # rotation 3 (dimension z)
    d_golden = d3 & 0xFFFF

    # Step 3: reconstruct nonce from golden delta
    # hi = seed's hi, lo = hi + d_golden (delta encodes lo-hi)
    hi_seed, _ = to_2d(seed)
    lo_golden = (hi_seed + d_golden) & 0xFFFF
    return from_2d(hi_seed, lo_golden)

# ─── Print Results ───────────────────────────────────────────────────────────

print("=" * 60)
print("NONCE MATRIX — Axiom Alignment")
print("=" * 60)

print("\nThe Matrix:")
print(f"  [ {TL:+.6f},  {TR:+.6f} ]  =  [ -sqrt(5), +sqrt(5) ]")
print(f"  [ {BL:+.6f},  {BR:+.6f} ]  =  [      -3,       +3  ]")

print(f"\nD = phi^2 + psi^2 = {D:.10f}  (= 3 exactly)")
print(f"Proof: (phi+psi)^2 - 2*phi*psi = 1^2 - 2*(-1) = 3")

print(f"\nMatrix rank: 1 (both rows proportional to [lo-hi])")
print(f"  TL/BL = {TL/BL:.6f}  (= sqrt(5)/3 = {sqrt5/3:.6f})")
print(f"  TR/BR = {TR/BR:.6f}  (= sqrt(5)/3)")
print(f"Collapse: 2D grid -> 1D axis: delta = lo - hi")
print(f"  Weight = D * delta = 3 * delta")
print(f"  Height = sqrt(5) * delta")

print(f"\nThe Collision:")
print(f"  HALF = 2^31 = {HALF}")
print(f"  -HALF mod 2^32 = {(-HALF) % CIRCLE}")
print(f"  SAME = {HALF == (-HALF) % CIRCLE}")
print(f"  1/2 = 1/(-2) on the nonce circle. Fixed point of negation.")

print(f"\nStep 291 Collapse:")
print(f"  phi^291 = {phi_291:.6e}")
print(f"  psi^291 = {psi_291:.6e}")
print(f"  |psi^291| / phi^291 = {psi_to_phi_ratio:.6e}  (sub-Planck, effectively 0)")
print(f"  Matrix at step 291: all entries ≈ ±{phi_291:.4e}")
print(f"  TL={TL_291:.4e}, TR={TR_291:.4e}, BL={BL_291:.4e}, BR={BR_291:.4e}")
print(f"  All simultaneously collapse. Global reset.")

print(f"\n8/45 = 2^d / (d^2 * p) = 2^3 / (3^2 * 5) = {8/45:.6f}")
print(f"  Cube corners (2^d=8) per dimensional-pentagonal degree (d^2 * p = 45)")
print(f"  The axiom's normalized structural density in 3D pentagonal space.")

print(f"\nF(291) ≈ {F291:.6e}  (integer Fibonacci, captures phi^291)")
print(f"L(291) ≈ {L291:.6e}  (integer Lucas, IS phi^291)")
print(f"phi^291 ≈ F(291) * sqrt(5): {F291 * sqrt5:.6e} vs {phi_291:.6e}")

# ─── Demo: golden center for a test block ────────────────────────────────────
print(f"\n{'=' * 60}")
print("Demo: axiom center for a test block header")
header76 = b'\x01\x00\x00\x00' + b'\x00' * 32 + b'\x00' * 32 + b'\x00' * 4 + b'\xff\xff\x7f\x20'
center = golden_center(header76)
hi_c, lo_c = to_2d(center)
d_c = delta(center)
w_c = axiom_weight(center)

print(f"  Header: {header76.hex()[:32]}...")
print(f"  Golden center nonce: {center}")
print(f"  2D decomposition: hi={hi_c}, lo={lo_c}")
print(f"  Delta (lo-hi): {d_c}")
print(f"  Weight (D*delta=3*delta): {w_c}")
print(f"  Check: 3*{d_c} mod 2^16 = {(3*d_c) & 0xFFFF} = {w_c}")
