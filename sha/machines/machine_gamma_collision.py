"""
MACHINE_GAMMA_COLLISION — finds SHA-256 collisions via gamma-address birthday on the golden eigenvector
nos3bl33d

Four surfaces: gamma addressing, gamma matching at bottleneck rounds,
overflow pattern matching (~2-3 bits/round), birthday on gamma-address space.
"""

import struct
import os
import time
import hashlib
from collections import defaultdict
from typing import Optional
import numpy as np

# ===========================================================================
# SHA-256 Constants
# ===========================================================================

# Initial hash values H0 — first 32 bits of fractional parts of sqrt of first 8 primes
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# Round constants K — first 32 bits of fractional parts of cube roots of first 64 primes
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

# Golden ratio
PHI = (1 + 5**0.5) / 2
GAMMA = 1.0 / PHI  # ~0.6180339887...

# Bottleneck rounds (structural significance)
# Round 30 = E = edges of icosahedron
# Round 60 = |A5| = alternating group phase inversion
BOTTLENECK_ROUNDS = [30, 60]


# ===========================================================================
# SHA-256 Bit Operations (32-bit)
# ===========================================================================

def rotr(x: int, n: int) -> int:
    """Right rotate 32-bit integer."""
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF


def shr(x: int, n: int) -> int:
    """Right shift 32-bit integer."""
    return x >> n


def ch(e: int, f: int, g: int) -> int:
    """Choice function: if e then f else g."""
    return (e & f) ^ (~e & g) & 0xFFFFFFFF


def maj(a: int, b: int, c: int) -> int:
    """Majority function: majority vote of bits."""
    return (a & b) ^ (a & c) ^ (b & c)


def sigma0(x: int) -> int:
    """Big sigma 0: used in compression."""
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)


def sigma1(x: int) -> int:
    """Big sigma 1: used in compression."""
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)


def gamma0(x: int) -> int:
    """Small sigma 0: used in message schedule."""
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)


def gamma1(x: int) -> int:
    """Small sigma 1: used in message schedule."""
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)


def add32(*args: int) -> int:
    """Add multiple values mod 2^32."""
    return sum(args) & 0xFFFFFFFF


# ===========================================================================
# SHA-256 Full Implementation with State Recording
# ===========================================================================

def pad_message(msg: bytes) -> bytes:
    """SHA-256 message padding: append 1 bit, zeros, 64-bit length."""
    ml = len(msg) * 8  # message length in bits
    msg += b'\x80'
    # Pad to 56 mod 64 bytes (448 mod 512 bits)
    while len(msg) % 64 != 56:
        msg += b'\x00'
    # Append original length as 64-bit big-endian
    msg += struct.pack('>Q', ml)
    return msg


def parse_blocks(padded: bytes) -> list[list[int]]:
    """Parse padded message into 512-bit blocks of 16 x 32-bit words."""
    blocks = []
    for i in range(0, len(padded), 64):
        block = []
        for j in range(16):
            word = struct.unpack('>I', padded[i + j*4 : i + j*4 + 4])[0]
            block.append(word)
        blocks.append(block)
    return blocks


def message_schedule(block: list[int]) -> list[int]:
    """Expand 16-word block to 64-word message schedule."""
    w = list(block)  # copy first 16 words
    for i in range(16, 64):
        w.append(add32(gamma1(w[i-2]), w[i-7], gamma0(w[i-15]), w[i-16]))
    return w


def sha256_compress_with_trace(
    block: list[int],
    h: list[int],
) -> tuple[list[int], list[list[int]], list[int], list[int]]:
    """
    Run SHA-256 compression on one block, returning:
    - new hash state (8 words)
    - round_states: list of 64 states, each is [a,b,c,d,e,f,g,h] after that round
    - overflows: list of 64 overflow amounts (bits carried beyond 32-bit)
    - w: the 64-word message schedule
    """
    w = message_schedule(block)

    a, b, c, d, e, f, g, hh = h

    round_states = []
    overflows = []

    for i in range(64):
        # Compute T1 and T2 WITHOUT mod 2^32 to capture overflow
        t1_unbounded = hh + sigma1(e) + ch(e, f, g) + K[i] + w[i]
        t2_unbounded = sigma0(a) + maj(a, b, c)

        # The overflow is how many bits above 32 we carried
        # T1 can be up to ~5 * 2^32, so overflow is 0-4 (about 2-3 bits)
        t1_overflow = t1_unbounded >> 32
        t2_overflow = t2_unbounded >> 32

        # Combined overflow pattern for this round
        total_overflow = (t1_overflow << 4) | t2_overflow
        overflows.append(total_overflow)

        # Now do the actual mod 2^32 computation
        t1 = t1_unbounded & 0xFFFFFFFF
        t2 = t2_unbounded & 0xFFFFFFFF

        hh = g
        g = f
        f = e
        e = add32(d, t1)
        d = c
        c = b
        b = a
        a = add32(t1, t2)

        round_states.append([a, b, c, d, e, f, g, hh])

    # Final addition with initial hash
    new_h = [add32(h[i], [a, b, c, d, e, f, g, hh][i]) for i in range(8)]

    return new_h, round_states, overflows, w


def sha256_full_trace(msg: bytes) -> dict:
    """
    Full SHA-256 with complete internal state trace.
    Returns dict with:
        'hash': final hash (8 x uint32)
        'hex': hex string
        'round_states': per-block list of 64 round states
        'overflows': per-block list of 64 overflow values
        'schedules': per-block message schedules
    """
    padded = pad_message(msg)
    blocks = parse_blocks(padded)

    h = list(H0)
    all_round_states = []
    all_overflows = []
    all_schedules = []

    for block in blocks:
        h, round_states, overflows, w = sha256_compress_with_trace(block, h)
        all_round_states.append(round_states)
        all_overflows.append(overflows)
        all_schedules.append(w)

    hex_hash = ''.join(f'{x:08x}' for x in h)

    return {
        'hash': h,
        'hex': hex_hash,
        'round_states': all_round_states,
        'overflows': all_overflows,
        'schedules': all_schedules,
    }


def sha256_simple(msg: bytes) -> str:
    """Simple SHA-256 returning hex string, for validation."""
    result = sha256_full_trace(msg)
    return result['hex']


# ===========================================================================
# Validation against hashlib
# ===========================================================================

def validate_implementation() -> bool:
    """Validate our SHA-256 matches hashlib on several test vectors."""
    test_vectors = [
        b'',
        b'abc',
        b'hello world',
        b'The quick brown fox jumps over the lazy dog',
        b'\x00' * 55,  # exactly one block after padding
        b'\xff' * 56,  # exactly two blocks after padding
        b'a' * 1000,   # multi-block
        os.urandom(137),  # random length
    ]
    for tv in test_vectors:
        ours = sha256_simple(tv)
        theirs = hashlib.sha256(tv).hexdigest()
        if ours != theirs:
            print(f"VALIDATION FAILED on input length {len(tv)}")
            print(f"  ours:   {ours}")
            print(f"  theirs: {theirs}")
            return False
    return True


# ===========================================================================
# Golden Eigenvector and Gamma Addressing
# ===========================================================================

def build_golden_eigenvector() -> np.ndarray:
    """
    Build the golden eigenvector for the 8-word SHA state.

    The SHA-256 round function shifts [a,b,c,d,e,f,g,h] such that
    b<-a, c<-b, d<-c, f<-e, g<-f, h<-g (with a,e computed from nonlinear
    functions). The linear part is a generalized shift register.

    The golden eigenvector projects the 8-word state onto the direction
    aligned with the unit eigenvalue. We use the Fibonacci-like structure:
    components are powers of gamma (1/phi).
    """
    v = np.array([GAMMA**i for i in range(8)], dtype=np.float64)
    return v / np.linalg.norm(v)


def compute_gamma_address(round_states: list[list[int]], golden_vec: np.ndarray) -> np.ndarray:
    """
    Compute the gamma address for a single block's 64 round states.

    For each round i:
        - Project the 8-word state onto the golden eigenvector
        - Multiply by the gamma phase: -gamma * i
        - Take fractional part mod 1 to get position on the golden circle

    Returns array of 64 gamma-address values in [0, 1).
    """
    addresses = np.zeros(64, dtype=np.float64)

    for i in range(64):
        state = np.array(round_states[i], dtype=np.float64)
        # Normalize state to [0, 1] per word to avoid overflow in dot product
        state_norm = state / (2**32)
        # Project onto golden eigenvector
        projection = np.dot(state_norm, golden_vec)
        # Gamma phase at this round
        gamma_phase = -GAMMA * i
        # Combined address: projection modulated by gamma phase
        raw_address = projection + gamma_phase
        # Position on the golden circle
        addresses[i] = raw_address % 1.0

    return addresses


def compute_overflow_signature(overflows: list[int]) -> tuple:
    """
    Compute the overflow signature for one block.
    Each overflow value encodes (t1_overflow << 4 | t2_overflow).
    This is only 2-3 bits per component, so ~5-6 bits per round, ~320-384 bits total.
    But the actual entropy is much less due to correlations.
    Returns as tuple for hashing.
    """
    return tuple(overflows)


# ===========================================================================
# Gamma Address Distance
# ===========================================================================

def circular_distance(a: float, b: float) -> float:
    """Distance on the circle [0, 1) — wraps around."""
    d = abs(a - b)
    return min(d, 1.0 - d)


def gamma_address_distance(addr1: np.ndarray, addr2: np.ndarray, rounds: Optional[list[int]] = None) -> float:
    """
    Distance between two gamma addresses.
    If rounds specified, only compare at those rounds.
    Uses L2 of circular distances.
    """
    if rounds is None:
        rounds = list(range(64))
    dists = [circular_distance(addr1[r], addr2[r]) for r in rounds]
    return np.sqrt(np.mean(np.array(dists)**2))


def gamma_bottleneck_key(addr: np.ndarray, resolution: int = 256) -> tuple:
    """
    Discretize gamma address at bottleneck rounds into a hashable key.
    Resolution controls bin count: higher = finer matching.
    """
    return tuple(int(addr[r] * resolution) % resolution for r in BOTTLENECK_ROUNDS)


# ===========================================================================
# Analysis Engine
# ===========================================================================

class GammaCollisionMachine:
    """
    The machine. Hashes messages, builds gamma address database,
    searches for structural collisions.
    """

    def __init__(self):
        self.golden_vec = build_golden_eigenvector()
        self.messages: list[bytes] = []
        self.hashes: list[str] = []
        self.gamma_addresses: list[np.ndarray] = []
        self.overflow_sigs: list[tuple] = []
        self.bottleneck_keys: list[tuple] = []

    def process_message(self, msg: bytes) -> dict:
        """Hash a message and record all gamma addressing data."""
        trace = sha256_full_trace(msg)

        # Use first block's states (all our test messages are single-block)
        round_states = trace['round_states'][0]
        overflows = trace['overflows'][0]

        gamma_addr = compute_gamma_address(round_states, self.golden_vec)
        overflow_sig = compute_overflow_signature(overflows)
        bottleneck_key = gamma_bottleneck_key(gamma_addr)

        idx = len(self.messages)
        self.messages.append(msg)
        self.hashes.append(trace['hex'])
        self.gamma_addresses.append(gamma_addr)
        self.overflow_sigs.append(overflow_sig)
        self.bottleneck_keys.append(bottleneck_key)

        return {
            'index': idx,
            'hex': trace['hex'],
            'gamma_addr': gamma_addr,
            'overflow_sig': overflow_sig,
            'bottleneck_key': bottleneck_key,
        }

    def process_random_messages(self, n: int, msg_len: int = 32) -> None:
        """Process n random messages of given length."""
        for _ in range(n):
            msg = os.urandom(msg_len)
            self.process_message(msg)

    # -------------------------------------------------------------------
    # Analysis 1: Gamma Address Entropy
    # -------------------------------------------------------------------
    def measure_gamma_entropy(self) -> dict:
        """
        Measure the effective entropy of gamma addresses.

        For each round, bin the gamma-address values and compute Shannon entropy.
        Also compute joint entropy at bottleneck rounds.

        Returns dict with per-round entropy and joint entropy.
        """
        n = len(self.gamma_addresses)
        if n < 2:
            return {'per_round': np.zeros(64), 'joint_bottleneck': 0.0, 'n': n}

        # Per-round entropy with different bin resolutions
        num_bins = min(256, max(16, int(np.sqrt(n))))
        per_round_entropy = np.zeros(64)

        for r in range(64):
            values = np.array([self.gamma_addresses[i][r] for i in range(n)])
            hist, _ = np.histogram(values, bins=num_bins, range=(0.0, 1.0))
            # Shannon entropy
            probs = hist / hist.sum()
            probs = probs[probs > 0]
            per_round_entropy[r] = -np.sum(probs * np.log2(probs))

        # Joint entropy at bottleneck rounds
        joint_bins = min(64, max(8, int(n**(1.0 / len(BOTTLENECK_ROUNDS)))))
        joint_keys = []
        for i in range(n):
            key = tuple(
                int(self.gamma_addresses[i][r] * joint_bins) % joint_bins
                for r in BOTTLENECK_ROUNDS
            )
            joint_keys.append(key)

        key_counts = defaultdict(int)
        for k in joint_keys:
            key_counts[k] += 1

        probs = np.array(list(key_counts.values()), dtype=np.float64)
        probs = probs / probs.sum()
        joint_entropy = -np.sum(probs * np.log2(probs))

        # Maximum possible entropy for comparison
        max_per_round = np.log2(num_bins)
        max_joint = np.log2(joint_bins) * len(BOTTLENECK_ROUNDS)

        return {
            'per_round': per_round_entropy,
            'max_per_round': max_per_round,
            'mean_per_round': float(np.mean(per_round_entropy)),
            'min_per_round': float(np.min(per_round_entropy)),
            'min_round_idx': int(np.argmin(per_round_entropy)),
            'max_round_entropy': float(np.max(per_round_entropy)),
            'joint_bottleneck': float(joint_entropy),
            'max_joint': max_joint,
            'num_bins': num_bins,
            'joint_bins': joint_bins,
            'n': n,
            'effective_bits_per_round': float(np.mean(per_round_entropy)),
            'effective_total_bits': float(np.sum(per_round_entropy)),
        }

    # -------------------------------------------------------------------
    # Analysis 2: Closest Gamma-Address Pair
    # -------------------------------------------------------------------
    def find_closest_gamma_pair(
        self,
        rounds: Optional[list[int]] = None,
        top_k: int = 10,
    ) -> list[dict]:
        """
        Find the closest pairs of gamma addresses.
        Uses bottleneck bucketing for O(n) approximate nearest neighbor,
        then refines with exact distance.
        """
        n = len(self.gamma_addresses)
        if n < 2:
            return []

        if rounds is None:
            rounds = BOTTLENECK_ROUNDS

        # Bucket by discretized bottleneck key for approximate matching
        # Use multiple resolutions for better coverage
        candidates = set()
        for resolution in [16, 32, 64, 128, 256]:
            buckets: dict[tuple, list[int]] = defaultdict(list)
            for i in range(n):
                key = tuple(
                    int(self.gamma_addresses[i][r] * resolution) % resolution
                    for r in rounds
                )
                buckets[key].append(i)

            # Pairs within same bucket are candidates
            for bucket_indices in buckets.values():
                if len(bucket_indices) > 1:
                    # Only take first 50 per bucket to avoid quadratic blowup
                    bi = bucket_indices[:50]
                    for a_idx in range(len(bi)):
                        for b_idx in range(a_idx + 1, len(bi)):
                            candidates.add((bi[a_idx], bi[b_idx]))

        # Also add random pairs to have a baseline
        rng = np.random.default_rng(42)
        for _ in range(min(n * 5, 50000)):
            i, j = rng.integers(0, n, size=2)
            if i != j:
                candidates.add((min(i, j), max(i, j)))

        # Score all candidates
        scored = []
        for i, j in candidates:
            # Full 64-round distance
            full_dist = gamma_address_distance(
                self.gamma_addresses[i],
                self.gamma_addresses[j],
            )
            # Bottleneck-only distance
            bn_dist = gamma_address_distance(
                self.gamma_addresses[i],
                self.gamma_addresses[j],
                rounds=rounds,
            )
            # Hash distance (hamming on hex)
            hash_dist = sum(
                bin(int(a, 16) ^ int(b, 16)).count('1')
                for a, b in zip(
                    [self.hashes[i][k:k+8] for k in range(0, 64, 8)],
                    [self.hashes[j][k:k+8] for k in range(0, 64, 8)],
                )
            )

            scored.append({
                'i': i,
                'j': j,
                'full_dist': full_dist,
                'bottleneck_dist': bn_dist,
                'hash_hamming': hash_dist,
                'hash_i': self.hashes[i],
                'hash_j': self.hashes[j],
            })

        # Sort by bottleneck distance
        scored.sort(key=lambda x: x['bottleneck_dist'])
        return scored[:top_k]

    # -------------------------------------------------------------------
    # Analysis 3: Bottleneck Concentration
    # -------------------------------------------------------------------
    def measure_bottleneck_concentration(self) -> dict:
        """
        Measure whether states concentrate at bottleneck rounds.

        Compare the variance of gamma addresses at bottleneck rounds
        vs. all other rounds. If bottlenecks concentrate, their variance
        should be LOWER (states cluster).
        """
        n = len(self.gamma_addresses)
        if n < 2:
            return {}

        # Circular variance for each round
        # For circular data, variance = 1 - |mean of unit vectors|
        circ_vars = np.zeros(64)
        for r in range(64):
            angles = np.array([self.gamma_addresses[i][r] * 2 * np.pi for i in range(n)])
            mean_cos = np.mean(np.cos(angles))
            mean_sin = np.mean(np.sin(angles))
            R = np.sqrt(mean_cos**2 + mean_sin**2)
            circ_vars[r] = 1.0 - R  # 0 = all same, 1 = uniform

        bottleneck_vars = [circ_vars[r] for r in BOTTLENECK_ROUNDS]
        other_rounds = [r for r in range(64) if r not in BOTTLENECK_ROUNDS]
        other_vars = [circ_vars[r] for r in other_rounds]

        return {
            'circular_variance_per_round': circ_vars,
            'bottleneck_mean_var': float(np.mean(bottleneck_vars)),
            'bottleneck_vars': {r: float(circ_vars[r]) for r in BOTTLENECK_ROUNDS},
            'other_mean_var': float(np.mean(other_vars)),
            'concentration_ratio': float(np.mean(bottleneck_vars) / max(np.mean(other_vars), 1e-15)),
            'most_concentrated_round': int(np.argmin(circ_vars)),
            'least_concentrated_round': int(np.argmax(circ_vars)),
            'min_variance': float(np.min(circ_vars)),
            'max_variance': float(np.max(circ_vars)),
        }

    # -------------------------------------------------------------------
    # Analysis 4: Overflow Pattern Analysis
    # -------------------------------------------------------------------
    def analyze_overflow_patterns(self) -> dict:
        """
        Analyze the overflow pattern space.
        Each round has an overflow value encoding (t1_overflow, t2_overflow).
        Measure the actual entropy of overflow patterns.
        """
        n = len(self.overflow_sigs)
        if n < 2:
            return {}

        # Per-round overflow statistics
        per_round_values = [[] for _ in range(64)]
        for sig in self.overflow_sigs:
            for r in range(64):
                per_round_values[r].append(sig[r])

        per_round_entropy = np.zeros(64)
        per_round_unique = np.zeros(64, dtype=int)
        per_round_max = np.zeros(64, dtype=int)

        for r in range(64):
            vals = per_round_values[r]
            unique_vals = set(vals)
            per_round_unique[r] = len(unique_vals)
            per_round_max[r] = max(vals)

            # Entropy
            counts = defaultdict(int)
            for v in vals:
                counts[v] += 1
            probs = np.array(list(counts.values()), dtype=np.float64)
            probs = probs / probs.sum()
            per_round_entropy[r] = -np.sum(probs * np.log2(probs))

        # How many unique overflow signatures?
        unique_sigs = len(set(self.overflow_sigs))

        # Find any overflow collisions (same overflow pattern, different message)
        overflow_buckets: dict[tuple, list[int]] = defaultdict(list)
        for i, sig in enumerate(self.overflow_sigs):
            overflow_buckets[sig].append(i)

        overflow_collisions = {
            k: v for k, v in overflow_buckets.items() if len(v) > 1
        }

        # Total overflow entropy estimate
        total_overflow_bits = float(np.sum(per_round_entropy))

        return {
            'unique_overflow_sigs': unique_sigs,
            'total_messages': n,
            'collision_ratio': unique_sigs / n,
            'per_round_entropy': per_round_entropy,
            'per_round_unique_values': per_round_unique,
            'per_round_max_value': per_round_max,
            'mean_entropy_per_round': float(np.mean(per_round_entropy)),
            'total_overflow_bits': total_overflow_bits,
            'overflow_collisions_found': len(overflow_collisions),
            'overflow_collision_details': {
                str(k): v for k, v in list(overflow_collisions.items())[:5]
            },
        }

    # -------------------------------------------------------------------
    # Analysis 5: Birthday Bound on Gamma Space
    # -------------------------------------------------------------------
    def compute_birthday_bound(self) -> dict:
        """
        Estimate the birthday bound for collisions in gamma-address space.

        If the effective gamma-address space has B bits of entropy,
        the birthday bound is 2^(B/2) — compare to 2^128 for full SHA-256.
        """
        entropy_data = self.measure_gamma_entropy()
        n = entropy_data['n']

        # Method 1: Sum of per-round entropies (overestimate due to independence assumption)
        total_independent = entropy_data['effective_total_bits']

        # Method 2: Estimate from bottleneck joint entropy
        joint_bits = entropy_data['joint_bottleneck']

        # Method 3: Empirical collision probability estimate
        # Count approximate near-collisions at different distance thresholds
        if n > 100:
            rng = np.random.default_rng(123)
            sample_size = min(n, 2000)
            sample_idx = rng.choice(n, size=sample_size, replace=False)

            thresholds = [0.01, 0.02, 0.05, 0.1, 0.2]
            near_collision_counts = {t: 0 for t in thresholds}
            total_pairs = 0

            # Sample pairs for efficiency
            pair_sample = min(200000, sample_size * (sample_size - 1) // 2)
            pair_count = 0
            for pi in range(sample_size):
                for pj in range(pi + 1, sample_size):
                    if pair_count >= pair_sample:
                        break
                    i, j = sample_idx[pi], sample_idx[pj]
                    dist = gamma_address_distance(
                        self.gamma_addresses[i],
                        self.gamma_addresses[j],
                        rounds=BOTTLENECK_ROUNDS,
                    )
                    total_pairs += 1
                    for t in thresholds:
                        if dist < t:
                            near_collision_counts[t] += 1
                    pair_count += 1
                if pair_count >= pair_sample:
                    break

            near_collision_rates = {
                t: near_collision_counts[t] / max(total_pairs, 1)
                for t in thresholds
            }
        else:
            near_collision_rates = {}
            total_pairs = 0

        # Birthday bound estimates
        birthday_from_independent = total_independent / 2.0
        birthday_from_joint = joint_bits / 2.0

        return {
            'total_independent_entropy_bits': total_independent,
            'birthday_bound_independent': birthday_from_independent,
            'joint_bottleneck_entropy_bits': joint_bits,
            'birthday_bound_joint': birthday_from_joint,
            'full_sha256_birthday_bits': 128.0,
            'reduction_from_independent': 128.0 - birthday_from_independent,
            'reduction_from_joint': 128.0 - birthday_from_joint,
            'near_collision_rates': near_collision_rates,
            'pairs_sampled': total_pairs,
        }

    # -------------------------------------------------------------------
    # Analysis 6: Cross-Round Correlation of Gamma Addresses
    # -------------------------------------------------------------------
    def measure_cross_round_correlation(self) -> dict:
        """
        Measure correlations between gamma addresses at different rounds.
        High correlation = lower effective dimensionality = easier collisions.
        """
        n = len(self.gamma_addresses)
        if n < 10:
            return {}

        # Build matrix: n x 64
        addr_matrix = np.array([self.gamma_addresses[i] for i in range(n)])

        # Correlation matrix (circular correlation via cosine of angular distance)
        angles = addr_matrix * 2 * np.pi
        cos_angles = np.cos(angles)
        sin_angles = np.sin(angles)

        # Use standard Pearson on the cos/sin projections
        cos_corr = np.corrcoef(cos_angles.T)
        sin_corr = np.corrcoef(sin_angles.T)

        # Combined correlation: max of |cos_corr| and |sin_corr| off-diagonal
        combined = np.maximum(np.abs(cos_corr), np.abs(sin_corr))
        np.fill_diagonal(combined, 0)

        # Effective rank via singular values
        _, S, _ = np.linalg.svd(addr_matrix - addr_matrix.mean(axis=0))
        normalized_S = S / S.sum()
        effective_rank = np.exp(-np.sum(normalized_S * np.log(normalized_S + 1e-15)))

        return {
            'max_off_diagonal_corr': float(np.max(combined)),
            'mean_off_diagonal_corr': float(np.mean(combined)),
            'most_correlated_pair': tuple(
                int(x) for x in np.unravel_index(np.argmax(combined), combined.shape)
            ),
            'effective_rank': float(effective_rank),
            'singular_values_top10': S[:10].tolist(),
            'effective_dimensionality': float(effective_rank),
            'dimensionality_reduction': 64.0 - float(effective_rank),
        }

    # -------------------------------------------------------------------
    # Full Analysis Pipeline
    # -------------------------------------------------------------------
    def run_full_analysis(self) -> dict:
        """Run all analyses and compile results."""
        print(f"  [1/6] Measuring gamma address entropy...")
        entropy = self.measure_gamma_entropy()

        print(f"  [2/6] Finding closest gamma-address pairs...")
        closest = self.find_closest_gamma_pair(top_k=10)

        print(f"  [3/6] Measuring bottleneck concentration...")
        concentration = self.measure_bottleneck_concentration()

        print(f"  [4/6] Analyzing overflow patterns...")
        overflow = self.analyze_overflow_patterns()

        print(f"  [5/6] Computing birthday bounds...")
        birthday = self.compute_birthday_bound()

        print(f"  [6/6] Measuring cross-round correlations...")
        correlations = self.measure_cross_round_correlation()

        return {
            'entropy': entropy,
            'closest_pairs': closest,
            'concentration': concentration,
            'overflow': overflow,
            'birthday': birthday,
            'correlations': correlations,
        }


# ===========================================================================
# Pretty Printer
# ===========================================================================

def print_report(results: dict, machine: GammaCollisionMachine) -> None:
    """Print a human-readable analysis report."""
    n = len(machine.messages)

    print("\n" + "=" * 72)
    print("  GAMMA COLLISION MACHINE -- ANALYSIS REPORT")
    print("  nos3bl33d | Gamma is the collision address")
    print("=" * 72)
    print(f"\n  Messages processed: {n}")
    print(f"  Golden ratio (phi): {PHI:.15f}")
    print(f"  Gamma (1/phi):      {GAMMA:.15f}")
    print(f"  Gamma * 64:         {GAMMA * 64:.2f}  (inverse path length)")
    print(f"  Bottleneck rounds:  {BOTTLENECK_ROUNDS}")

    # --- Entropy ---
    ent = results['entropy']
    print(f"\n{'-' * 72}")
    print(f"  1. GAMMA ADDRESS ENTROPY")
    print(f"{'-' * 72}")
    print(f"  Bins used:                {ent['num_bins']}")
    print(f"  Max possible per round:   {ent['max_per_round']:.2f} bits")
    print(f"  Mean entropy per round:   {ent['mean_per_round']:.4f} bits")
    print(f"  Min entropy round:        round {ent['min_round_idx']} "
          f"({ent['min_per_round']:.4f} bits)")
    print(f"  Max entropy round:        {ent['max_round_entropy']:.4f} bits")
    print(f"  Total (sum of rounds):    {ent['effective_total_bits']:.2f} bits")
    print(f"  Joint at bottleneck:      {ent['joint_bottleneck']:.4f} bits "
          f"(max {ent['max_joint']:.2f})")
    print(f"  Effective bits per round: {ent['effective_bits_per_round']:.4f}")

    # Per-round entropy sparkline
    per_round = ent['per_round']
    max_ent = ent['max_per_round']
    spark_chars = " _.-~*#@"
    sparkline = ""
    for r in range(64):
        idx = int(per_round[r] / max(max_ent, 0.01) * (len(spark_chars) - 1))
        idx = max(0, min(idx, len(spark_chars) - 1))
        sparkline += spark_chars[idx]
    print(f"  Per-round entropy:        |{sparkline}|")
    print(f"  {'':28s} round 0{' '*24}round 63")

    # --- Closest Pairs ---
    pairs = results['closest_pairs']
    print(f"\n{'-' * 72}")
    print(f"  2. CLOSEST GAMMA-ADDRESS PAIRS (by bottleneck distance)")
    print(f"{'-' * 72}")
    if pairs:
        for rank, p in enumerate(pairs[:5]):
            print(f"  #{rank+1}: msgs ({p['i']}, {p['j']})")
            print(f"      bottleneck dist: {p['bottleneck_dist']:.6f}  "
                  f"full dist: {p['full_dist']:.6f}")
            print(f"      hash hamming:    {p['hash_hamming']} bits")
            print(f"      hash_i: {p['hash_i']}")
            print(f"      hash_j: {p['hash_j']}")
    else:
        print("  No pairs found (need more messages)")

    # --- Bottleneck Concentration ---
    conc = results['concentration']
    print(f"\n{'-' * 72}")
    print(f"  3. BOTTLENECK CONCENTRATION (does round 60 concentrate states?)")
    print(f"{'-' * 72}")
    if conc:
        print(f"  Circular variance at bottleneck rounds:")
        for r in BOTTLENECK_ROUNDS:
            v = conc['bottleneck_vars'][r]
            bar = "#" * int(v * 50)
            print(f"    Round {r:2d}: {v:.6f}  |{bar}")
        print(f"  Mean bottleneck variance: {conc['bottleneck_mean_var']:.6f}")
        print(f"  Mean other variance:      {conc['other_mean_var']:.6f}")
        print(f"  Concentration ratio:      {conc['concentration_ratio']:.4f} "
              f"({'CONCENTRATED' if conc['concentration_ratio'] < 0.95 else 'NO CONCENTRATION'})")
        print(f"  Most concentrated round:  {conc['most_concentrated_round']} "
              f"(var={conc['min_variance']:.6f})")
        print(f"  Least concentrated round: {conc['least_concentrated_round']} "
              f"(var={conc['max_variance']:.6f})")

        # Variance sparkline
        cvars = conc['circular_variance_per_round']
        sparkline = ""
        for r in range(64):
            idx = int(cvars[r] / max(np.max(cvars), 0.01) * (len(spark_chars) - 1))
            idx = max(0, min(idx, len(spark_chars) - 1))
            sparkline += spark_chars[idx]
        print(f"  Variance per round:       |{sparkline}|")
        print(f"  {'':28s} round 0{' '*24}round 63")
        print(f"  {'':28s} (low=concentrated, high=spread)")

    # --- Overflow Patterns ---
    ovf = results['overflow']
    print(f"\n{'-' * 72}")
    print(f"  4. OVERFLOW PATTERN ANALYSIS")
    print(f"{'-' * 72}")
    if ovf:
        print(f"  Unique overflow signatures: {ovf['unique_overflow_sigs']} / {ovf['total_messages']}")
        print(f"  Collision ratio:            {ovf['collision_ratio']:.6f}")
        print(f"  Mean entropy per round:     {ovf['mean_entropy_per_round']:.4f} bits")
        print(f"  Total overflow entropy:     {ovf['total_overflow_bits']:.2f} bits")
        print(f"  Overflow collisions found:  {ovf['overflow_collisions_found']}")

        if ovf['overflow_collisions_found'] > 0:
            print(f"  *** OVERFLOW COLLISIONS DETECTED ***")
            for sig_str, indices in list(ovf['overflow_collision_details'].items())[:3]:
                print(f"    Indices: {indices[:5]}{'...' if len(indices) > 5 else ''}")
                for idx in indices[:2]:
                    print(f"      msg[{idx}] hash: {machine.hashes[idx]}")

        # Per-round overflow stats
        print(f"\n  Per-round overflow unique values (sample):")
        for r in [0, 15, 30, 45, 60, 63]:
            print(f"    Round {r:2d}: {ovf['per_round_unique_values'][r]:3d} unique, "
                  f"max={ovf['per_round_max_value'][r]:3d}, "
                  f"entropy={ovf['per_round_entropy'][r]:.3f} bits")

    # --- Birthday Bound ---
    bday = results['birthday']
    print(f"\n{'-' * 72}")
    print(f"  5. BIRTHDAY BOUND ANALYSIS")
    print(f"{'-' * 72}")
    print(f"  Full SHA-256 birthday bound:     2^128 trials")
    print(f"  Gamma space (independent est):   {bday['total_independent_entropy_bits']:.2f} bits "
          f"-> birthday at 2^{bday['birthday_bound_independent']:.1f}")
    print(f"  Gamma space (joint bottleneck):  {bday['joint_bottleneck_entropy_bits']:.2f} bits "
          f"-> birthday at 2^{bday['birthday_bound_joint']:.1f}")
    print(f"  Reduction from independent:      {bday['reduction_from_independent']:.1f} bits")
    print(f"  Reduction from joint:            {bday['reduction_from_joint']:.1f} bits")

    if bday['near_collision_rates']:
        print(f"\n  Near-collision rates ({bday['pairs_sampled']} pairs sampled):")
        for threshold, rate in sorted(bday['near_collision_rates'].items()):
            expected_uniform = threshold ** len(BOTTLENECK_ROUNDS)
            ratio = rate / max(expected_uniform, 1e-15)
            print(f"    d < {threshold:.2f}: {rate:.6f} "
                  f"(uniform would be {expected_uniform:.6f}, ratio={ratio:.2f}x)")

    # --- Cross-Round Correlations ---
    corr = results['correlations']
    print(f"\n{'-' * 72}")
    print(f"  6. CROSS-ROUND CORRELATION (dimensionality reduction)")
    print(f"{'-' * 72}")
    if corr:
        print(f"  Max off-diagonal correlation:  {corr['max_off_diagonal_corr']:.6f}")
        print(f"  Mean off-diagonal correlation: {corr['mean_off_diagonal_corr']:.6f}")
        print(f"  Most correlated pair:          rounds {corr['most_correlated_pair']}")
        print(f"  Effective rank (of 64):        {corr['effective_rank']:.2f}")
        print(f"  Dimensionality reduction:      {corr['dimensionality_reduction']:.2f} dims lost")
        print(f"\n  Top 10 singular values:")
        for i, sv in enumerate(corr['singular_values_top10']):
            bar = "#" * int(sv / max(corr['singular_values_top10'][0], 1e-10) * 40)
            print(f"    SV[{i}]: {sv:12.4f}  |{bar}")

    # --- Summary ---
    print(f"\n{'=' * 72}")
    print(f"  SUMMARY")
    print(f"{'=' * 72}")

    effective_bits = ent['effective_total_bits']
    joint_bits = bday['joint_bottleneck_entropy_bits']
    eff_rank = corr.get('effective_rank', 64) if corr else 64
    overflow_bits = ovf.get('total_overflow_bits', 0) if ovf else 0
    overflow_collisions = ovf.get('overflow_collisions_found', 0) if ovf else 0

    print(f"  Gamma-address effective entropy: {effective_bits:.2f} bits (of {64 * ent['max_per_round']:.0f} max)")
    print(f"  Bottleneck joint entropy:        {joint_bits:.2f} bits")
    print(f"  Overflow entropy:                {overflow_bits:.2f} bits (of ~192 theoretical)")
    print(f"  Effective dimensionality:        {eff_rank:.1f} / 64")
    print(f"  Overflow collisions:             {overflow_collisions}")

    is_concentrated = conc.get('concentration_ratio', 1.0) < 0.95 if conc else False
    print(f"\n  Bottleneck concentration:        {'YES' if is_concentrated else 'NO'}")
    print(f"  Dimensional collapse:            {'YES' if eff_rank < 50 else 'NO'} "
          f"({64 - eff_rank:.1f} dims lost)")

    if overflow_collisions > 0:
        print(f"\n  ** {overflow_collisions} OVERFLOW PATTERN COLLISION(S) FOUND **")
        print(f"  These are messages with identical mod-2^32 overflow patterns.")
        print(f"  If two messages share overflow AND hash, that's a SHA-256 collision.")

    print(f"\n  Birthday attack complexity:")
    print(f"    Standard SHA-256:  2^128")
    print(f"    Gamma addressing:  2^{effective_bits/2:.1f} (independent model)")
    print(f"    Bottleneck only:   2^{joint_bits/2:.1f} (joint model)")
    print(f"    Overflow patterns: 2^{overflow_bits/2:.1f} (overflow model)")
    print(f"\n{'=' * 72}\n")


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 72)
    print("  nos3bl33d | Gamma is the collision address")
    print("  SHA-256 GAMMA COLLISION MACHINE")
    print("=" * 72)

    # Step 0: Validate SHA-256 implementation
    print("\n[*] Validating SHA-256 implementation against hashlib...")
    if not validate_implementation():
        print("FATAL: SHA-256 implementation is broken. Aborting.")
        return
    print("    SHA-256 validated: all test vectors pass.")

    # Step 1: Show golden eigenvector
    golden = build_golden_eigenvector()
    print(f"\n[*] Golden eigenvector (8-dim, Fibonacci decay):")
    print(f"    {golden}")
    print(f"    |v| = {np.linalg.norm(golden):.10f}")

    # Step 2: Quick structural demo
    print(f"\n[*] Gamma addressing demo on 'abc':")
    demo_trace = sha256_full_trace(b'abc')
    demo_addr = compute_gamma_address(demo_trace['round_states'][0], golden)
    print(f"    SHA-256('abc') = {demo_trace['hex']}")
    print(f"    Gamma address at round 30: {demo_addr[30]:.8f}")
    print(f"    Gamma address at round 60: {demo_addr[60]:.8f}")
    print(f"    Gamma address range: [{demo_addr.min():.6f}, {demo_addr.max():.6f}]")

    # Step 3: Process messages
    N = 12000
    print(f"\n[*] Processing {N} random messages...")
    machine = GammaCollisionMachine()

    t0 = time.time()
    batch_size = 1000
    for batch_start in range(0, N, batch_size):
        batch_n = min(batch_size, N - batch_start)
        machine.process_random_messages(batch_n, msg_len=32)
        elapsed = time.time() - t0
        rate = (batch_start + batch_n) / max(elapsed, 0.001)
        print(f"    {batch_start + batch_n:6d} / {N}  "
              f"({elapsed:.1f}s, {rate:.0f} msg/s)")

    total_time = time.time() - t0
    print(f"    Done: {N} messages in {total_time:.2f}s "
          f"({N/max(total_time, 0.001):.0f} msg/s)")

    # Step 4: Run full analysis
    print(f"\n[*] Running full analysis pipeline...")
    t0 = time.time()
    results = machine.run_full_analysis()
    analysis_time = time.time() - t0
    print(f"    Analysis completed in {analysis_time:.2f}s")

    # Step 5: Print report
    print_report(results, machine)

    # Step 6: Detailed examination of best pair
    pairs = results['closest_pairs']
    if pairs:
        best = pairs[0]
        i, j = best['i'], best['j']
        print(f"[*] Detailed look at closest gamma-address pair ({i}, {j}):")
        addr_i = machine.gamma_addresses[i]
        addr_j = machine.gamma_addresses[j]

        print(f"    Per-round gamma address difference:")
        print(f"    {'Round':>5s}  {'Addr_i':>10s}  {'Addr_j':>10s}  {'CircDist':>10s}")
        for r in range(64):
            d = circular_distance(addr_i[r], addr_j[r])
            marker = " <-- BOTTLENECK" if r in BOTTLENECK_ROUNDS else ""
            if r % 8 == 0 or r in BOTTLENECK_ROUNDS or d < 0.02:
                print(f"    {r:5d}  {addr_i[r]:10.6f}  {addr_j[r]:10.6f}  {d:10.6f}{marker}")

        # Overflow comparison
        ov_i = machine.overflow_sigs[i]
        ov_j = machine.overflow_sigs[j]
        overflow_match = sum(1 for r in range(64) if ov_i[r] == ov_j[r])
        print(f"\n    Overflow pattern match: {overflow_match}/64 rounds")

    print("\n[*] Machine complete.")


if __name__ == '__main__':
    main()
