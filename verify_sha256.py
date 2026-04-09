"""
ADVERSARIAL VERIFICATION: SHA-256 structural claims
1. Characteristic polynomial of linearized round matrix
2. Jacobian determinant at random inputs
3. GF(2) message schedule rank
4. Fibonacci exponents
"""

import numpy as np
from numpy.linalg import det, eigvals
import hashlib
import struct
import os

print("=" * 70)
print("SHA-256 STRUCTURAL VERIFICATION")
print("=" * 70)

# ======================================================================
# CLAIM 1: Characteristic polynomial of linearized round matrix
# ======================================================================

print("\n" + "=" * 70)
print("CLAIM 1: char(M) = x^8 - 2x^7 + x^6 - x^4 - 1")
print("and char(M) + 1 = x^4 * (x^2 - x - 1) * (x^2 - x + 1)")
print("=" * 70)

# The linearized matrix M from the paper
M = np.array([
    [1, 0, 0, 0, 1, 0, 0, 1],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 1],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
], dtype=float)

# Compute characteristic polynomial using numpy
coeffs = np.round(np.poly(M)).astype(int)
print(f"\nCharacteristic polynomial coefficients (highest to lowest degree):")
print(f"  {coeffs}")

# Expected: x^8 - 2x^7 + x^6 + 0x^5 - x^4 + 0x^3 + 0x^2 + 0x - 1
expected_coeffs = [1, -2, 1, 0, -1, 0, 0, 0, -1]
print(f"Expected: {expected_coeffs}")

match = np.array_equal(coeffs, expected_coeffs)
print(f"MATCH: {match}")

# Verify factorization: char(M) + 1 = x^4 * (x^2 - x - 1) * (x^2 - x + 1)
# char(M) + 1 means adding 1 to the constant term
# char(M)(x) + 1 = x^8 - 2x^7 + x^6 - x^4 - 1 + 1 = x^8 - 2x^7 + x^6 - x^4
# = x^4 * (x^4 - 2x^3 + x^2 - 1)

# Check: (x^2 - x - 1)(x^2 - x + 1) = x^4 - x^3 + x^2 - x^3 + x^2 - x - x^2 + x + 1
# = x^4 - 2x^3 + x^2 - 1... wait let me recompute

# (x^2 - x - 1)(x^2 - x + 1) = (x^2 - x)^2 - 1 = x^4 - 2x^3 + x^2 - 1
print("\nVerifying factorization:")
print("  char(M) + 1 = x^8 - 2x^7 + x^6 - x^4")
print("  = x^4 * (x^4 - 2x^3 + x^2 - 1)")
print("  Need: x^4 - 2x^3 + x^2 - 1 = (x^2-x-1)(x^2-x+1)")

# Expand (x^2 - x - 1)(x^2 - x + 1) = ((x^2-x) - 1)((x^2-x) + 1) = (x^2-x)^2 - 1
# = x^4 - 2x^3 + x^2 - 1
from sympy import symbols, expand, factor, Poly
x = symbols('x')
product = (x**2 - x - 1) * (x**2 - x + 1)
expanded = expand(product)
print(f"  (x^2-x-1)(x^2-x+1) = {expanded}")
target = x**4 - 2*x**3 + x**2 - 1
print(f"  x^4 - 2x^3 + x^2 - 1 = {expand(target)}")
print(f"  Match: {expand(expanded - target) == 0}")

# Full check
char_poly = x**8 - 2*x**7 + x**6 - x**4 - 1
char_plus_1 = char_poly + 1
factored_form = x**4 * (x**2 - x - 1) * (x**2 - x + 1)
print(f"\n  char(M) + 1 = {expand(char_plus_1)}")
print(f"  x^4*(x^2-x-1)*(x^2-x+1) = {expand(factored_form)}")
print(f"  FACTORIZATION VERIFIED: {expand(char_plus_1 - factored_form) == 0}")

# Verify the golden ratio is actually a root
phi = (1 + 5**0.5) / 2
psi = (1 - 5**0.5) / 2
print(f"\n  phi^2 - phi - 1 = {phi**2 - phi - 1:.2e} (should be ~0)")
print(f"  psi^2 - psi - 1 = {psi**2 - psi - 1:.2e} (should be ~0)")

# Check eigenvalues of M
eigs = eigvals(M)
print(f"\n  Eigenvalues of M: {np.sort(np.real(eigs))}")
print(f"  Expected roots: 0 (x4), phi={phi:.6f}, psi={psi:.6f}, e^(i*pi/3), e^(-i*pi/3)")

# ======================================================================
# CLAIM 2: Jacobian determinant = -1 every round
# ======================================================================

print("\n" + "=" * 70)
print("CLAIM 2: det(Jacobian) = -1 every round, EVERY input")
print("THIS IS THE ADVERSARIAL TEST")
print("=" * 70)

# SHA-256 constants
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

MOD = 2**32

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def shr(x, n):
    return x >> n

def ch(e, f, g):
    return (e & f) ^ (~e & g) & 0xFFFFFFFF

def maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def sigma0(a):
    return rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

def sigma1(e):
    return rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

def lsigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def lsigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def sha256_round(state, ki, wi):
    """One SHA-256 round. state = [a,b,c,d,e,f,g,h]"""
    a, b, c, d, e, f, g, h = state
    t1 = (h + sigma1(e) + ch(e, f, g) + ki + wi) % MOD
    t2 = (sigma0(a) + maj(a, b, c)) % MOD
    new_a = (t1 + t2) % MOD
    new_e = (d + t1) % MOD
    return [new_a, a, b, c, new_e, e, f, g]

def compute_message_schedule(block_words):
    """Compute full 64-word message schedule from 16 input words"""
    W = list(block_words)
    for i in range(16, 64):
        W.append((lsigma1(W[i-2]) + W[i-7] + lsigma0(W[i-15]) + W[i-16]) % MOD)
    return W

def numerical_jacobian(state, ki, wi, eps=1):
    """Compute 8x8 Jacobian of one SHA-256 round by finite difference"""
    base = sha256_round(state, ki, wi)
    J = np.zeros((8, 8))

    for j in range(8):
        # Perturb state[j] by +eps
        state_plus = list(state)
        state_plus[j] = (state_plus[j] + eps) % MOD

        state_minus = list(state)
        state_minus[j] = (state_minus[j] - eps) % MOD

        out_plus = sha256_round(state_plus, ki, wi)
        out_minus = sha256_round(state_minus, ki, wi)

        for i in range(8):
            # Central difference
            diff = (out_plus[i] - out_minus[i])
            # Handle modular wraparound
            if diff > MOD // 2:
                diff -= MOD
            elif diff < -MOD // 2:
                diff += MOD
            J[i][j] = diff / (2 * eps)

    return J


print("\nTesting Jacobian determinant at RANDOM inputs...")
print("(The claim says det = -1 for ALL inputs at ALL rounds)")

NUM_RANDOM_TESTS = 10
all_det_minus_1 = True

for test_idx in range(NUM_RANDOM_TESTS):
    # Generate random 512-bit message block
    random_block = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]
    W = compute_message_schedule(random_block)

    # Random initial state
    state = [int.from_bytes(os.urandom(4), 'big') for _ in range(8)]

    round_results = []
    for round_num in range(64):
        J = numerical_jacobian(state, K[round_num], W[round_num])
        d = det(J)

        if abs(d - (-1)) > 0.01:
            round_results.append((round_num, d, "FAIL"))
            all_det_minus_1 = False
        else:
            round_results.append((round_num, d, "OK"))

        # Advance state
        state = sha256_round(state, K[round_num], W[round_num])

    failures = [r for r in round_results if r[2] == "FAIL"]
    if failures:
        print(f"  Test {test_idx}: {len(failures)} rounds FAILED det=-1!")
        for rn, dv, _ in failures[:5]:
            print(f"    Round {rn}: det = {dv:.6f}")
    else:
        det_vals = [r[1] for r in round_results]
        print(f"  Test {test_idx}: 64/64 rounds det=-1 (range: [{min(det_vals):.8f}, {max(det_vals):.8f}])")

# Also test with the standard H0 initial values
print("\nTesting with standard SHA-256 initial values H0...")
for test_idx in range(5):
    random_block = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]
    W = compute_message_schedule(random_block)
    state = list(H0)

    round_results = []
    for round_num in range(64):
        J = numerical_jacobian(state, K[round_num], W[round_num])
        d = det(J)

        if abs(d - (-1)) > 0.01:
            round_results.append((round_num, d, "FAIL"))
            all_det_minus_1 = False

        state = sha256_round(state, K[round_num], W[round_num])

    failures_count = sum(1 for r in round_results if r[2] == "FAIL")
    print(f"  H0 test {test_idx}: {64 - failures_count}/64 rounds det=-1")

# ADVERSARIAL: Test with extreme inputs (all zeros, all ones, alternating)
print("\nAdversarial: extreme inputs...")
extreme_states = [
    [0]*8,
    [0xFFFFFFFF]*8,
    [0xAAAAAAAA]*8,
    [0x55555555]*8,
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
]
extreme_messages = [
    [0]*16,
    [0xFFFFFFFF]*16,
    [i for i in range(16)],
]

for si, state_init in enumerate(extreme_states):
    for mi, msg in enumerate(extreme_messages):
        W = compute_message_schedule(msg)
        state = list(state_init)
        fail_count = 0
        fail_rounds = []

        for round_num in range(64):
            J = numerical_jacobian(state, K[round_num], W[round_num])
            d = det(J)

            if abs(d - (-1)) > 0.01:
                fail_count += 1
                fail_rounds.append((round_num, d))
                all_det_minus_1 = False

            state = sha256_round(state, K[round_num], W[round_num])

        if fail_count > 0:
            print(f"  state={si} msg={mi}: {fail_count} FAILURES!")
            for rn, dv in fail_rounds[:3]:
                print(f"    Round {rn}: det = {dv}")
        else:
            print(f"  state={si} msg={mi}: 64/64 OK")

print(f"\nOVERALL det=-1 claim: {'VERIFIED' if all_det_minus_1 else 'BROKEN'}")

# ======================================================================
# CLAIM 3: GF(2) message schedule rank = 512
# ======================================================================

print("\n" + "=" * 70)
print("CLAIM 3: GF(2) message schedule rank = 512")
print("=" * 70)

def gf2_rotr(bit_index, word_index, n, total_bits=32):
    """For a 32-bit word at position word_index, ROTR by n returns the source bit position"""
    return word_index * 32 + (bit_index + n) % 32

def gf2_shr(bit_index, word_index, n, total_bits=32):
    """SHR by n: if bit_index + n >= 32, result is 0 (no source bit)"""
    src = bit_index + n
    if src >= 32:
        return None  # zero
    return word_index * 32 + src

def build_gf2_schedule_matrix():
    """
    Build the 2048x512 GF(2) matrix representing the SHA-256 message schedule.
    Input: 512 bits (16 words x 32 bits)
    Output: 2048 bits (64 words x 32 bits)

    For words 0-15: output = input (identity)
    For words 16-63: W[i] = sigma1(W[i-2]) XOR W[i-7] XOR sigma0(W[i-15]) XOR W[i-16]
    where sigma0(x) = ROTR7(x) XOR ROTR18(x) XOR SHR3(x)
          sigma1(x) = ROTR17(x) XOR ROTR19(x) XOR SHR10(x)
    """
    # Over GF(2), addition = XOR, and all operations are linear.
    # We build a 2048x512 matrix where rows = output bits, cols = input bits (W[0..15])

    # First, we need to track which input bits contribute to each output bit
    # We'll represent each of the 2048 output bits as a 512-dimensional GF(2) vector

    # Initialize: each W[i] (for i=0..63) is a 32-element array of 512-bit vectors
    # W[i][bit] = set of input bit indices that XOR to produce this output bit

    W = []
    for i in range(64):
        word = [set() for _ in range(32)]
        if i < 16:
            # Identity: W[i][bit] depends on input bit i*32+bit
            for b in range(32):
                word[b] = {i * 32 + b}
        W.append(word)

    # Now compute W[16] through W[63]
    for i in range(16, 64):
        for b in range(32):
            sources = set()

            # sigma1(W[i-2]): ROTR17 XOR ROTR19 XOR SHR10
            # ROTR17: bit b of result comes from bit (b+17)%32 of input
            sources ^= W[i-2][(b+17) % 32]
            # ROTR19
            sources ^= W[i-2][(b+19) % 32]
            # SHR10: bit b comes from bit b+10 if b+10 < 32
            if b + 10 < 32:
                sources ^= W[i-2][b+10]

            # W[i-7]
            sources ^= W[i-7][b]

            # sigma0(W[i-15]): ROTR7 XOR ROTR18 XOR SHR3
            sources ^= W[i-15][(b+7) % 32]
            sources ^= W[i-15][(b+18) % 32]
            if b + 3 < 32:
                sources ^= W[i-15][b+3]

            # W[i-16]
            sources ^= W[i-16][b]

            W[i][b] = sources

    # Build the matrix (2048 x 512) over GF(2)
    matrix = np.zeros((2048, 512), dtype=np.uint8)
    for i in range(64):
        for b in range(32):
            row = i * 32 + b
            for src_bit in W[i][b]:
                matrix[row][src_bit] = 1

    return matrix

print("Building GF(2) message schedule matrix (2048 x 512)...")
M_gf2 = build_gf2_schedule_matrix()

# Compute rank over GF(2) using Gaussian elimination
def gf2_rank(matrix):
    """Compute rank of a binary matrix over GF(2)"""
    m = matrix.copy()
    rows, cols = m.shape
    rank = 0
    pivot_col = 0

    for row in range(rows):
        if pivot_col >= cols:
            break

        # Find pivot
        found = False
        for r in range(rank, rows):
            if m[r][pivot_col] == 1:
                # Swap rows
                m[[rank, r]] = m[[r, rank]]
                found = True
                break

        if not found:
            pivot_col += 1
            continue

        # Eliminate
        for r in range(rows):
            if r != rank and m[r][pivot_col] == 1:
                m[r] = m[r] ^ m[rank]

        rank += 1
        pivot_col += 1

    return rank

# Actually, for a 2048x512 matrix we need a smarter approach
# The first 512 rows are the identity (words 0-15 = input),
# so the rank is trivially 512.

print("First 512 rows form identity block (W[0..15] = input)")
identity_check = np.array_equal(M_gf2[:512, :512], np.eye(512, dtype=np.uint8))
print(f"  Identity block verified: {identity_check}")

# But let's verify by computing actual rank of the full matrix
print("Computing full GF(2) rank (this verifies no degeneracy)...")

# Use the transpose for efficiency (512x2048 has fewer rows)
# rank(M) = rank(M^T) and we want rank of the column space of a 2048x512 matrix
# which equals rank of the 512x2048 matrix

# Since we know the first 512 rows are identity, rank >= 512
# And since the matrix has 512 columns, rank <= 512
# Therefore rank = 512 EXACTLY.

print("  Matrix is 2048 x 512")
print("  First 512 rows = I_512 (identity)")
print("  rank >= 512 (from identity block)")
print("  rank <= 512 (column count)")
print("  Therefore rank = 512 EXACTLY")
print("  CLAIM VERIFIED (by construction)")

# Double-check: make sure the identity block is correct
for i in range(16):
    for b in range(32):
        row = i * 32 + b
        expected_col = i * 32 + b
        assert M_gf2[row][expected_col] == 1, f"Missing diagonal at ({row},{expected_col})"
        for c in range(512):
            if c != expected_col:
                assert M_gf2[row][c] == 0, f"Non-zero off-diagonal at ({row},{c})"

print("  Identity block double-checked: CORRECT")

# Also verify that the GF(2) rank is EXACTLY 512 by doing explicit Gaussian elimination
# on a subset to be thorough
print("  Running explicit GF(2) Gaussian elimination on full matrix...")
rank_val = gf2_rank(M_gf2)
print(f"  GF(2) rank = {rank_val}")
print(f"  CLAIM {'VERIFIED' if rank_val == 512 else 'BROKEN'}")

# ======================================================================
# CLAIM 4: Fibonacci in the exponents
# ======================================================================

print("\n" + "=" * 70)
print("CLAIM 4: Fibonacci numbers in SHA-256 parameters")
print("=" * 70)

# Fibonacci sequence
fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
print(f"Fibonacci sequence: {fib}")
print(f"F(5) = {fib[4]} = 5")
print(f"F(6) = {fib[5]} = 8")
print(f"F(7) = {fib[6]} = 13")

print(f"\nSHA-256 word size: 32 bits = 2^5 = 2^F(5)")
print(f"  32 = 2^5: {32 == 2**5}")
print(f"  5 = F(5): {5 == fib[4]}")

print(f"\nSHA-256 hash size: 256 bits = 2^8 = 2^F(6)")
print(f"  256 = 2^8: {256 == 2**8}")
print(f"  8 = F(6): {8 == fib[5]}")

print(f"\n256 = 8 * 32 = F(6) * 2^F(5): {256 == 8 * 32}")

print(f"\n5, 8, 13 are consecutive Fibonacci: F(5)=5, F(6)=8, F(7)=13")
print(f"  F(6) = F(5) + F(4) = 5 + 3 = 8: {fib[5] == fib[4] + fib[3]}")
print(f"  F(7) = F(6) + F(5) = 8 + 5 = 13: {fib[6] == fib[5] + fib[4]}")

print(f"\nSHA-256 rounds: 64 = 8^2 = (state words)^2")
print(f"  64 = 8^2: {64 == 8**2}")

print(f"\nAll Fibonacci-exponent claims: VERIFIED (these are arithmetic identities)")

# ======================================================================
# CLAIM 5: Eigenvalue analysis of linearized matrix
# ======================================================================

print("\n" + "=" * 70)
print("ADDITIONAL: Eigenvalue analysis of linearized M")
print("=" * 70)

eigs = eigvals(M)
print(f"Eigenvalues of M:")
for i, e in enumerate(sorted(eigs, key=lambda x: abs(x))):
    print(f"  lambda_{i} = {e:.6f} (|lambda| = {abs(e):.6f})")

# Check if phi and psi are eigenvalues
phi_val = (1 + 5**0.5) / 2
psi_val = (1 - 5**0.5) / 2
print(f"\nphi = {phi_val:.6f}")
print(f"psi = {psi_val:.6f}")

# Find closest eigenvalues to phi and psi
for target_name, target_val in [("phi", phi_val), ("psi", psi_val)]:
    dists = [(abs(e - target_val), e) for e in eigs]
    closest_dist, closest_eig = min(dists, key=lambda x: x[0])
    print(f"Closest eigenvalue to {target_name}: {closest_eig:.6f} (distance: {closest_dist:.2e})")

# Note: phi and psi are roots of char(M)+1=0, not char(M)=0
# char(M)(phi) = phi^8 - 2*phi^7 + phi^6 - phi^4 - 1
char_at_phi = phi_val**8 - 2*phi_val**7 + phi_val**6 - phi_val**4 - 1
print(f"\nchar(M)(phi) = {char_at_phi:.6f} (should be nonzero; phi is NOT an eigenvalue of M)")
print(f"char(M)(phi) + 1 = {char_at_phi + 1:.2e} (should be ~0; phi IS a root of char(M)+1)")

print("\nNOTE: The golden ratio is a root of char(M)+1, NOT of char(M) itself.")
print("This is correctly stated in the paper: 'The golden polynomial does not divide")
print("char(M) itself: the remainder of char(M) mod (x^2-x-1) is -1.'")

print("\n" + "=" * 70)
print("ALL SHA-256 CLAIMS VERIFICATION COMPLETE")
print("=" * 70)
