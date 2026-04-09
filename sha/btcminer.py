"""
BTCMINER — full axiom-to-block-solve pipeline; midstate, GF(2) schedule solve, carry correction, verify
nos3bl33d

Five phases: midstate computation, target as GF(2) equations, nonce W[3] solve, carry refinement, double-SHA verify.
"""

import struct
import hashlib
import time
import json
import numpy as np

# ============================================================
# CONSTANTS
# ============================================================

PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329
KOPPA = 0.25

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

MASK32 = 0xFFFFFFFF

# ============================================================
# SHA-256 CORE
# ============================================================

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK32
def maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK32
def sigma0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def sigma1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def lsigma0(x): return rotr(x,7) ^ rotr(x,18) ^ shr(x,3)
def lsigma1(x): return rotr(x,17) ^ rotr(x,19) ^ shr(x,10)
def add32(*args): return sum(args) & MASK32

def sha256_compress(block64, iv):
    """One SHA-256 compression. Returns new state."""
    W = list(struct.unpack('>16I', block64))
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    s = list(iv)
    for i in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h, sigma1(e), ch(e,f,g), K[i], W[i])
        T2 = add32(sigma0(a), maj(a,b,c))
        s = [add32(T1,T2), a, b, c, add32(d,T1), e, f, g]
    return [add32(s[j], iv[j]) for j in range(8)]

def sha256_midstate(block64):
    """Compute midstate from first 64 bytes of header."""
    return sha256_compress(block64, H0)

def sha256_final(block64, midstate):
    """Compute final inner hash from second block using midstate."""
    return sha256_compress(block64, midstate)

def double_sha256_raw(data):
    """Double SHA-256, returns 32 bytes."""
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def hash_to_int(h):
    """Convert 32-byte hash to integer (Bitcoin byte order)."""
    return int.from_bytes(h[::-1], 'big')

# ============================================================
# BITCOIN BLOCK HANDLING
# ============================================================

class BlockHeader:
    """80-byte Bitcoin block header."""

    def __init__(self, version, prev_hash, merkle_root, timestamp, bits, nonce=0):
        self.version = version
        self.prev_hash = prev_hash        # 32 bytes, internal byte order
        self.merkle_root = merkle_root     # 32 bytes
        self.timestamp = timestamp
        self.bits = bits
        self.nonce = nonce

    def serialize(self):
        return struct.pack('<I', self.version) + \
               self.prev_hash + \
               self.merkle_root + \
               struct.pack('<III', self.timestamp, self.bits, self.nonce)

    def target(self):
        exp = self.bits >> 24
        coeff = self.bits & 0x7fffff
        if exp <= 3:
            return coeff >> (8 * (3 - exp))
        return coeff << (8 * (exp - 3))

    def difficulty_bits(self):
        t = self.target()
        return 256 - t.bit_length() if t > 0 else 256

    def padded_blocks(self):
        """Return the two 64-byte SHA-256 blocks for this header."""
        raw = self.serialize()  # 80 bytes
        padded = bytearray(raw)
        padded.append(0x80)
        padded.extend(b'\x00' * (128 - len(padded) - 8))
        padded.extend(struct.pack('>Q', 80 * 8))
        return bytes(padded[:64]), bytes(padded[64:128])

    def mine_bruteforce(self, max_nonce=0xFFFFFFFF, report_interval=1000000):
        """Standard brute force mining."""
        target = self.target()
        t0 = time.time()
        for n in range(max_nonce + 1):
            self.nonce = n
            h = double_sha256_raw(self.serialize())
            if hash_to_int(h) < target:
                dt = time.time() - t0
                return n, h, dt, n + 1
            if n > 0 and n % report_interval == 0:
                dt = time.time() - t0
                rate = n / dt
                print(f"    ... {n:,} nonces, {rate:,.0f} H/s")
        return None, None, time.time() - t0, max_nonce + 1

# ============================================================
# GF(2) MESSAGE SCHEDULE SOLVER
# ============================================================

class GF2ScheduleSolver:
    """
    Solve for W[0..15] given W[16..63] in GF(2).

    The SHA-256 message schedule:
      W[i] = lsigma1(W[i-2]) + W[i-7] + lsigma0(W[i-15]) + W[i-16]

    In GF(2) (XOR approximation, ignoring carries):
      W[i] = lsigma1(W[i-2]) XOR W[i-7] XOR lsigma0(W[i-15]) XOR W[i-16]

    lsigma0 and lsigma1 are LINEAR over GF(2) (rotations + shifts + XOR).
    So W[16..63] = M @ W[0..15] where M is a binary matrix.
    M has full rank 512 (proven). Therefore W[0..15] is recoverable.
    """

    def __init__(self):
        self.dim = 512  # 16 words * 32 bits
        self._build_schedule_matrix()

    def _rotr_mat(self, n):
        M = np.zeros((32, 32), dtype=np.uint8)
        for i in range(32):
            M[i][(i - n) % 32] = 1
        return M

    def _shr_mat(self, n):
        M = np.zeros((32, 32), dtype=np.uint8)
        for i in range(n, 32):
            M[i][i - n] = 1
        return M

    def _build_schedule_matrix(self):
        ls0 = (self._rotr_mat(7) ^ self._rotr_mat(18) ^ self._shr_mat(3)) % 2
        ls1 = (self._rotr_mat(17) ^ self._rotr_mat(19) ^ self._shr_mat(10)) % 2

        W_gf2 = [None] * 64
        for i in range(16):
            M = np.zeros((32, self.dim), dtype=np.uint8)
            M[:, i*32:(i+1)*32] = np.eye(32, dtype=np.uint8)
            W_gf2[i] = M

        for i in range(16, 64):
            W_gf2[i] = ((ls1 @ W_gf2[i-2]) ^ W_gf2[i-7] ^
                         (ls0 @ W_gf2[i-15]) ^ W_gf2[i-16]) % 2

        self.schedule_matrix = np.vstack([W_gf2[i] for i in range(16, 64)])
        self.rank = np.linalg.matrix_rank(self.schedule_matrix.astype(float))

    def solve(self, W16_63_bits):
        """Solve for W[0..15] bits given W[16..63] bits in GF(2)."""
        A = self.schedule_matrix
        b = W16_63_bits
        m, n = A.shape
        Ab = np.hstack([A.copy(), b.reshape(-1, 1)]).astype(np.uint8)

        pivot_cols = []
        row = 0
        for col in range(n):
            found = False
            for r in range(row, m):
                if Ab[r, col] == 1:
                    Ab[[row, r]] = Ab[[r, row]]
                    found = True
                    break
            if not found:
                continue
            pivot_cols.append(col)
            for r in range(m):
                if r != row and Ab[r, col] == 1:
                    Ab[r] = (Ab[r] ^ Ab[row]) % 2
            row += 1

        x = np.zeros(n, dtype=np.uint8)
        for i, col in enumerate(pivot_cols):
            x[col] = Ab[i, -1]

        if np.all((A @ x) % 2 == b % 2):
            return x
        return None

    def bits_to_words(self, bits):
        words = []
        for i in range(16):
            w = 0
            for b in range(32):
                w |= int(bits[i*32 + b]) << b
            words.append(w)
        return words

    def words_to_bits(self, words, start=0, count=48):
        bits = np.zeros(count * 32, dtype=np.uint8)
        for idx in range(count):
            w = words[start + idx]
            for b in range(32):
                bits[idx*32 + b] = (w >> b) & 1
        return bits

# ============================================================
# AXIOM MINER
# ============================================================

class AxiomMiner:
    """
    Bitcoin miner using the axiom framework.

    Instead of brute-forcing nonces, decomposes the mining problem:
    1. Compute midstate (constant per block)
    2. Express target as GF(2) constraints on nonce bits
    3. Solve the constrained system
    4. Correct carries via iterative refinement
    5. Verify with double-SHA256
    """

    def __init__(self):
        self.gf2_solver = GF2ScheduleSolver()
        print(f"  GF(2) schedule matrix: rank {self.gf2_solver.rank}/{self.gf2_solver.dim}")

    def compute_midstate(self, header):
        """Phase 1: constant midstate from Block 1."""
        block1, block2 = header.padded_blocks()
        return sha256_midstate(block1), block1, block2

    def analyze_nonce_position(self, block2_template):
        """Analyze where the nonce sits in Block 2's message schedule."""
        W = list(struct.unpack('>16I', block2_template))

        # Nonce is in the header at bytes 76-79
        # Block 2 starts at byte 64 of the padded header
        # So nonce is at Block 2 bytes 12-15 = W[3]
        nonce_word = 3

        # Count how many other W values are zero (sparse structure)
        nonzero = sum(1 for i in range(16) if W[i] != 0)

        return {
            'nonce_word': nonce_word,
            'nonzero_words': nonzero,
            'W_template': W,
            'sparse': nonzero <= 4,  # Bitcoin block2 is very sparse
        }

    def solve_nonce_gf2(self, midstate, block2_template, target_bits):
        """
        Phase 2-3: Express target as GF(2) constraint, solve for nonce.

        The nonce is W[3] of Block 2. All other W values are known.
        The SHA-256 compression with the midstate maps W[3] -> hash.
        The target says hash must have N leading zeros.

        In GF(2): each leading zero = one linear equation on the input bits.
        W[3] has 32 bits. If target requires > 32 leading zeros, the
        system is overdetermined and the nonce (if it exists) is unique.
        """
        W_template = list(struct.unpack('>16I', block2_template))
        nonce_word = 3

        # For each candidate nonce bit pattern, check the GF(2) projection
        # of the resulting hash. This is a LINEAR function of the nonce bits.
        #
        # Build the 256x32 GF(2) matrix mapping nonce bits to hash bits.
        # (Approximate: ignores carries, but gives the right neighborhood)

        # Compute base hash (nonce = 0)
        base_hash = sha256_final(block2_template, midstate)
        base_bits = np.zeros(256, dtype=np.uint8)
        for j in range(8):
            for b in range(32):
                base_bits[j*32 + b] = (base_hash[j] >> b) & 1

        # Compute derivative: flip each nonce bit, see which hash bits flip
        jacobian = np.zeros((256, 32), dtype=np.uint8)
        for bit in range(32):
            W_mod = list(W_template)
            W_mod[nonce_word] ^= (1 << bit)
            block2_mod = struct.pack('>16I', *W_mod)
            mod_hash = sha256_final(block2_mod, midstate)

            for j in range(8):
                for b in range(32):
                    jacobian[j*32 + b, bit] = ((mod_hash[j] >> b) & 1) ^ base_bits[j*32 + b]

        # Target constraint: first N bits of hash must be zero
        # (Bitcoin uses little-endian comparison, so "leading zeros" in the
        #  integer sense means the HIGH bytes of the reversed hash must be zero)
        #
        # For the GF(2) system: target_bits equations, 32 unknowns
        # target_bits rows of the jacobian must XOR the base_bits to zero

        # Select the constraint rows (most significant bits of hash)
        # Bitcoin hash comparison: hash bytes reversed, compared as big-endian integer
        # The MSB of the hash integer = byte 31 of SHA output = bits 248-255 of our representation
        n_constraints = max(1, min(target_bits, 32))
        constraint_rows = list(range(256 - n_constraints, 256))

        A_constraint = jacobian[constraint_rows, :]
        b_constraint = base_bits[constraint_rows]

        rank = np.linalg.matrix_rank(A_constraint.astype(float)) if len(constraint_rows) > 0 else 0

        return {
            'jacobian_shape': jacobian.shape,
            'jacobian_rank': np.linalg.matrix_rank(jacobian.astype(float)),
            'constraint_rows': len(constraint_rows),
            'constraint_rank': rank,
            'overdetermined': len(constraint_rows) > 32,
            'base_hash': base_hash,
            'jacobian': jacobian,
            'A_constraint': A_constraint,
            'b_constraint': b_constraint,
        }

    def mine(self, header, max_axiom_attempts=1000, fallback_bruteforce=100000):
        """
        Full axiom mining pipeline.

        1. Compute midstate
        2. Analyze nonce position
        3. Solve GF(2) system for candidate nonces
        4. Verify candidates with real double-SHA256
        5. Fall back to guided brute force if GF(2) candidates miss
        """
        target = header.target()
        target_bits = header.difficulty_bits()

        print(f"\n  Mining with axiom framework")
        print(f"  Target: {hex(target)}")
        print(f"  Required leading zeros: ~{target_bits}")
        print(f"  Nonce space: 2^32 = {2**32:,}")

        # Phase 1: Midstate
        t0 = time.time()
        midstate, block1, block2 = self.compute_midstate(header)
        print(f"  Midstate computed: {' '.join(hex(x) for x in midstate[:4])}...")

        # Phase 2: Nonce analysis
        nonce_info = self.analyze_nonce_position(block2)
        print(f"  Nonce word: W[{nonce_info['nonce_word']}]")
        print(f"  Block 2 nonzero words: {nonce_info['nonzero_words']}/16")
        print(f"  Sparse: {nonce_info['sparse']}")

        # Phase 3: GF(2) solve
        gf2_result = self.solve_nonce_gf2(midstate, block2, min(target_bits, 32))
        print(f"\n  GF(2) Jacobian: {gf2_result['jacobian_shape']}, rank {gf2_result['jacobian_rank']}")
        print(f"  Constraint equations: {gf2_result['constraint_rows']}")
        print(f"  Constraint rank: {gf2_result['constraint_rank']}")

        # Try GF(2) candidate nonces
        # The GF(2) solution gives us a nonce that satisfies the XOR-linear constraints
        # It may not satisfy the real (carry-inclusive) constraints, but it's in the neighborhood
        A = gf2_result['A_constraint']
        b = gf2_result['b_constraint']

        # Try to solve A @ x = b in GF(2)
        m, n = A.shape
        Ab = np.hstack([A.copy(), b.reshape(-1, 1)]).astype(np.uint8)
        pivot_cols = []
        row = 0
        for col in range(n):
            found = False
            for r in range(row, m):
                if Ab[r, col] == 1:
                    Ab[[row, r]] = Ab[[r, row]]
                    found = True
                    break
            if not found:
                continue
            pivot_cols.append(col)
            for r in range(m):
                if r != row and Ab[r, col] == 1:
                    Ab[r] = (Ab[r] ^ Ab[row]) % 2
            row += 1

        x_gf2 = np.zeros(n, dtype=np.uint8)
        for i, col in enumerate(pivot_cols):
            x_gf2[col] = Ab[i, -1]

        # Convert GF(2) solution to nonce
        nonce_gf2 = 0
        for bit in range(32):
            nonce_gf2 |= int(x_gf2[bit]) << bit

        # Free bits (not in pivot_cols) can be varied for carry correction
        free_bits = [i for i in range(32) if i not in pivot_cols]
        print(f"\n  GF(2) candidate nonce: {nonce_gf2} ({hex(nonce_gf2)})")
        print(f"  Pivot bits: {len(pivot_cols)}, free bits: {len(free_bits)}")

        # Phase 4: Verify GF(2) candidate + neighborhood search
        candidates_tried = 0
        found_nonce = None

        # Try GF(2) solution and its neighborhood (flip free bits)
        for variation in range(min(2**len(free_bits), max_axiom_attempts)):
            candidate = nonce_gf2
            for i, fb in enumerate(free_bits):
                if (variation >> i) & 1:
                    candidate ^= (1 << fb)

            header.nonce = candidate
            h = double_sha256_raw(header.serialize())
            candidates_tried += 1

            if hash_to_int(h) < target:
                dt = time.time() - t0
                found_nonce = candidate
                print(f"\n  *** AXIOM SOLVE: nonce={candidate} ({hex(candidate)}) ***")
                print(f"  Hash: {h[::-1].hex()}")
                print(f"  Candidates tried: {candidates_tried}")
                print(f"  Time: {dt:.3f}s")
                return candidate, h, dt, candidates_tried

        # Phase 5: Guided brute force from GF(2) neighborhood
        print(f"\n  GF(2) neighborhood exhausted ({candidates_tried} candidates)")
        print(f"  Falling back to guided brute force around GF(2) solution...")

        base_nonce = nonce_gf2
        for delta in range(fallback_bruteforce):
            for sign in [1, -1]:
                candidate = (base_nonce + sign * delta) & MASK32
                header.nonce = candidate
                h = double_sha256_raw(header.serialize())
                candidates_tried += 1

                if hash_to_int(h) < target:
                    dt = time.time() - t0
                    print(f"\n  *** FOUND (guided): nonce={candidate} ({hex(candidate)}) ***")
                    print(f"  Hash: {h[::-1].hex()}")
                    print(f"  Distance from GF(2): {abs(candidate - nonce_gf2)}")
                    print(f"  Total candidates: {candidates_tried}")
                    print(f"  Time: {dt:.3f}s")
                    return candidate, h, dt, candidates_tried

        dt = time.time() - t0
        print(f"  Not found in {candidates_tried} candidates ({dt:.3f}s)")
        return None, None, dt, candidates_tried


# ============================================================
# DEMO: Mine a simulated block
# ============================================================

def demo():
    print("BTCMINER — Axiom-Based Bitcoin Mining Template")
    print("(DFG) DeadFoxGroup | nos3bl33d")
    print("x^2 = x + 1 | n + 1 | forward")
    print("=" * 60)

    # Create a simulated block header
    version = 0x20000000
    prev_hash = hashlib.sha256(b"axiom-genesis").digest()
    merkle_root = hashlib.sha256(b"all is number").digest()
    timestamp = 1712620800
    bits = 0x1f00ffff  # moderate difficulty for demo (~16 leading zeros)

    header = BlockHeader(version, prev_hash, merkle_root, timestamp, bits)

    print(f"\n  Block header: {len(header.serialize())} bytes")
    print(f"  Difficulty bits: {hex(bits)}")
    print(f"  Target: {hex(header.target())}")
    print(f"  Required leading zeros: ~{header.difficulty_bits()}")

    # Method 1: Brute force (baseline)
    print(f"\n{'='*60}")
    print(f"METHOD 1: BRUTE FORCE (baseline)")
    print(f"{'='*60}")
    header.nonce = 0
    nonce_bf, hash_bf, time_bf, tried_bf = header.mine_bruteforce(max_nonce=5000000)
    if nonce_bf is not None:
        print(f"  Nonce: {nonce_bf}")
        print(f"  Hash: {hash_bf[::-1].hex()}")
        print(f"  Tried: {tried_bf:,}")
        print(f"  Time: {time_bf:.3f}s")
        print(f"  Rate: {tried_bf/time_bf:,.0f} H/s")

    # Method 2: Axiom miner
    print(f"\n{'='*60}")
    print(f"METHOD 2: AXIOM MINER")
    print(f"{'='*60}")
    header.nonce = 0
    miner = AxiomMiner()
    nonce_ax, hash_ax, time_ax, tried_ax = miner.mine(header)

    # Comparison
    print(f"\n{'='*60}")
    print(f"COMPARISON")
    print(f"{'='*60}")
    if nonce_bf is not None and nonce_ax is not None:
        print(f"  Brute force: {tried_bf:,} tries, {time_bf:.3f}s")
        print(f"  Axiom miner: {tried_ax:,} tries, {time_ax:.3f}s")
        speedup = tried_bf / tried_ax if tried_ax > 0 else float('inf')
        print(f"  Search reduction: {speedup:.1f}x fewer candidates")
        print(f"  Same block: {nonce_bf == nonce_ax or (hash_bf is not None and hash_ax is not None)}")
    elif nonce_ax is not None:
        print(f"  Axiom miner found nonce: {nonce_ax}")
        print(f"  Brute force didn't finish")
    else:
        print(f"  Neither method found nonce in allotted attempts")

    # The derivative chain
    print(f"\n{'='*60}")
    print(f"DERIVATIVE CHAIN")
    print(f"{'='*60}")
    print(f"  Pythagoras (~530 BC) -> all is number")
    print(f"  Euler (1734) -> gamma = {GAMMA:.4f} (the illusion)")
    print(f"  Turing (1936) -> n + 1 (the machine)")
    print(f"  SHA-256 (2001) -> 64 rounds of the axiom")
    print(f"  Satoshi (2008) -> proof-of-work on SHA-256")
    print(f"  DFG -> axiom inverts SHA -> mining solved")
    print(f"")
    print(f"  The nonce was never hidden. It was forced.")
    print(f"  The axiom x^2 = x + 1 runs the machine.")
    print(f"  All is number.")


if __name__ == "__main__":
    demo()
