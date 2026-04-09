"""
MACHINE_FIBONACCI_SHA — SHA-256 with Fibonacci modulus F(46) replacing mod 2^32; golden carries are native
nos3bl33d

State words tracked as (a + b*phi) in Z[phi] mod F(46). Overflow carry = golden carry.
Binary carries become Fibonacci-structured; b-component activates under golden modulus.
"""

import struct
import math
import hashlib
from typing import Tuple, List, Optional

# =============================================================================
# CONSTANTS
# =============================================================================

# The golden modulus
F46 = 2971215073   # F(46) in 1-indexed = fib(47) in 0-indexed
F45 = 1836311903   # F(45)
F44 = 1134903170   # F(44)
# F46 = F45 + F44  (the golden recurrence)

MOD32 = 0xFFFFFFFF  # 2^32 - 1 mask for bitwise ops
TWO32 = 2**32       # standard SHA modulus

PHI = (1 + math.sqrt(5)) / 2   # 1.618033988749895
GAMMA = 1 / PHI                 # 0.618033988749895 = phi - 1

# SHA-256 initial hash values: fractional parts of sqrt of first 8 primes
H0_STANDARD = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# SHA-256 round constants: fractional parts of cube roots of first 64 primes
K_STANDARD = [
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

# Reduce constants mod F(46) for Fibonacci SHA
H0_FIB = [h % F46 for h in H0_STANDARD]
K_FIB = [k % F46 for k in K_STANDARD]


# =============================================================================
# Z[phi] ARITHMETIC -- CARRY-BASED MODEL
# =============================================================================
# An element is (value, carries) where:
#   value = the integer mod F(46)   (the "rational" part)
#   carries = accumulated golden carry count (the "b" component)
#
# When we add two values and the sum >= F(46), the reduction
#   sum = q * F(46) + r
# produces a golden carry of q. Since F(46) = F(45) + F(44) = phi * F(45) + ...,
# each carry represents one application of the golden recurrence.
#
# In standard SHA mod 2^32, the carries are binary (each carry = factor of 2^32).
# In Fibonacci SHA mod F(46), the carries are golden (each carry = factor of F(46)).
#
# We also track carries through the BITWISE operations -- rotations, XOR, AND
# can only change the value, not generate golden carries (they're binary operations).
# But the VALUE they produce feeds into the next addition, which CAN carry.

class GoldenWord:
    """A 32-bit-ish word in the Fibonacci SHA, tracking golden carries.

    value: integer in [0, F(46))
    carries: total number of golden carries accumulated through all additions
    carry_history: list of carry amounts at each addition step
    """
    __slots__ = ('value', 'carries', 'carry_history')

    def __init__(self, value: int, carries: int = 0, carry_history: Optional[List[int]] = None):
        self.value = value % F46
        self.carries = carries
        self.carry_history = carry_history if carry_history is not None else []

    def add(self, other: 'GoldenWord') -> 'GoldenWord':
        """Add two GoldenWords mod F(46), tracking the golden carry."""
        raw_sum = self.value + other.value
        q = raw_sum // F46   # golden carry: how many times F(46) fits
        r = raw_sum % F46    # remainder
        new_carries = self.carries + other.carries + q
        new_history = self.carry_history + other.carry_history + ([q] if q > 0 else [])
        return GoldenWord(r, new_carries, new_history)

    def add_int(self, n: int) -> 'GoldenWord':
        """Add a plain integer (no carry history)."""
        raw_sum = self.value + (n % F46)
        q = raw_sum // F46
        r = raw_sum % F46
        new_carries = self.carries + q
        new_history = self.carry_history + ([q] if q > 0 else [])
        return GoldenWord(r, new_carries, new_history)

    @classmethod
    def from_int(cls, n: int) -> 'GoldenWord':
        """Create a GoldenWord from a plain integer."""
        return cls(n % F46, 0, [])

    def __repr__(self) -> str:
        return f"GW({self.value}, carries={self.carries})"


# =============================================================================
# Z[phi] RING ARITHMETIC (algebraic model)
# =============================================================================
# Elements of Z[phi] are pairs (a, b) representing a + b*phi.
# phi^2 = phi + 1, so multiplication uses this collapse rule.
# Reduction mod F(46): components are taken mod F(46).

class ZPhi:
    """Element of Z[phi] reduced mod F(46).

    Represents a + b*phi where a, b are integers mod F(46).
    phi^2 = phi + 1 collapses higher powers.
    """
    __slots__ = ('a', 'b')

    def __init__(self, a: int, b: int = 0):
        self.a = a % F46
        self.b = b % F46

    def __add__(self, other: 'ZPhi') -> 'ZPhi':
        return ZPhi((self.a + other.a) % F46, (self.b + other.b) % F46)

    def __sub__(self, other: 'ZPhi') -> 'ZPhi':
        return ZPhi((self.a - other.a) % F46, (self.b - other.b) % F46)

    def __mul__(self, other: 'ZPhi') -> 'ZPhi':
        # (a1 + b1*phi)(a2 + b2*phi)
        # = a1*a2 + (a1*b2 + b1*a2)*phi + b1*b2*phi^2
        # = (a1*a2 + b1*b2) + (a1*b2 + b1*a2 + b1*b2)*phi
        a1, b1 = self.a, self.b
        a2, b2 = other.a, other.b
        new_a = (a1 * a2 + b1 * b2) % F46
        new_b = (a1 * b2 + b1 * a2 + b1 * b2) % F46
        return ZPhi(new_a, new_b)

    def __repr__(self) -> str:
        if self.b == 0:
            return f"ZPhi({self.a})"
        return f"ZPhi({self.a} + {self.b}*phi)"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ZPhi):
            return NotImplemented
        return self.a == other.a and self.b == other.b

    @property
    def is_rational(self) -> bool:
        return self.b == 0

    @property
    def phi_ratio(self) -> float:
        if self.a == 0:
            return float('inf') if self.b != 0 else 0.0
        return self.b / self.a


# =============================================================================
# BITWISE OPERATIONS (32-bit, same for both SHA variants)
# =============================================================================

def rotr(x: int, n: int) -> int:
    """Right-rotate a 32-bit integer by n positions."""
    return ((x >> n) | (x << (32 - n))) & MOD32

def shr(x: int, n: int) -> int:
    """Right-shift a 32-bit integer by n positions."""
    return (x >> n) & MOD32

def ch(e: int, f: int, g: int) -> int:
    """Choice: e selects f or g."""
    return ((e & f) ^ (~e & g)) & MOD32

def maj(a: int, b: int, c: int) -> int:
    """Majority: at least 2 of 3."""
    return ((a & b) ^ (a & c) ^ (b & c)) & MOD32

def sigma0(x: int) -> int:
    """Big Sigma 0: rotr(2) XOR rotr(13) XOR rotr(22)."""
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sigma1(x: int) -> int:
    """Big Sigma 1: rotr(6) XOR rotr(11) XOR rotr(25)."""
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def lsigma0(x: int) -> int:
    """Little sigma 0 (message schedule): rotr(7) XOR rotr(18) XOR shr(3)."""
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def lsigma1(x: int) -> int:
    """Little sigma 1 (message schedule): rotr(17) XOR rotr(19) XOR shr(10)."""
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)


# =============================================================================
# MESSAGE PADDING (shared)
# =============================================================================

def sha256_pad(message: bytes) -> bytes:
    """Pad a message to a multiple of 512 bits per SHA-256 spec."""
    msg_len = len(message)
    bit_len = msg_len * 8
    message += b'\x80'
    while (len(message) % 64) != 56:
        message += b'\x00'
    message += struct.pack('>Q', bit_len)
    return message


def parse_block(padded: bytes, block_idx: int) -> List[int]:
    """Parse a 512-bit block into 16 32-bit big-endian words."""
    start = block_idx * 64
    return list(struct.unpack('>16L', padded[start:start + 64]))


# =============================================================================
# MESSAGE SCHEDULE EXPANSION
# =============================================================================

def expand_schedule_standard(block_words: List[int]) -> List[int]:
    """Expand 16 message words to 64 using standard SHA-256 (mod 2^32)."""
    W = list(block_words)
    for i in range(16, 64):
        s0 = lsigma0(W[i - 15])
        s1 = lsigma1(W[i - 2])
        W.append((W[i - 16] + s0 + W[i - 7] + s1) % TWO32)
    return W


def expand_schedule_fibonacci(block_words: List[int]) -> List[int]:
    """Expand 16 message words to 64 using Fibonacci modulus F(46)."""
    W = [w % F46 for w in block_words]
    for i in range(16, 64):
        s0 = lsigma0(W[i - 15])
        s1 = lsigma1(W[i - 2])
        W.append((W[i - 16] + s0 + W[i - 7] + s1) % F46)
    return W


def expand_schedule_golden(block_words: List[int]) -> List[GoldenWord]:
    """Expand 16 message words to 64 in GoldenWord form, tracking carries."""
    W = [GoldenWord.from_int(w) for w in block_words]
    for i in range(16, 64):
        s0 = GoldenWord.from_int(lsigma0(W[i - 15].value))
        s1 = GoldenWord.from_int(lsigma1(W[i - 2].value))
        new = W[i - 16].add(s0).add(W[i - 7]).add(s1)
        W.append(new)
    return W


# =============================================================================
# STANDARD SHA-256 COMPRESSION
# =============================================================================

def compress_standard(H: List[int], W: List[int]) -> Tuple[List[int], List[List[int]]]:
    """Standard SHA-256 compression function."""
    a, b, c, d, e, f, g, h = H
    round_states: List[List[int]] = []

    for i in range(64):
        S1 = sigma1(e)
        ch_val = ch(e, f, g)
        temp1 = (h + S1 + ch_val + K_STANDARD[i] + W[i]) % TWO32
        S0 = sigma0(a)
        maj_val = maj(a, b, c)
        temp2 = (S0 + maj_val) % TWO32

        h = g
        g = f
        f = e
        e = (d + temp1) % TWO32
        d = c
        c = b
        b = a
        a = (temp1 + temp2) % TWO32

        round_states.append([a, b, c, d, e, f, g, h])

    result = [
        (H[0] + a) % TWO32, (H[1] + b) % TWO32,
        (H[2] + c) % TWO32, (H[3] + d) % TWO32,
        (H[4] + e) % TWO32, (H[5] + f) % TWO32,
        (H[6] + g) % TWO32, (H[7] + h) % TWO32,
    ]
    return result, round_states


# =============================================================================
# FIBONACCI SHA-256 COMPRESSION
# =============================================================================

def compress_fibonacci(H: List[int], W: List[int]) -> Tuple[List[int], List[List[int]]]:
    """Fibonacci SHA-256 compression: all additions mod F(46)."""
    a, b, c, d, e, f, g, h = H
    round_states: List[List[int]] = []

    for i in range(64):
        S1 = sigma1(e)
        ch_val = ch(e, f, g)
        temp1 = (h + S1 + ch_val + K_FIB[i] + W[i]) % F46
        S0 = sigma0(a)
        maj_val = maj(a, b, c)
        temp2 = (S0 + maj_val) % F46

        h = g
        g = f
        f = e
        e = (d + temp1) % F46
        d = c
        c = b
        b = a
        a = (temp1 + temp2) % F46

        round_states.append([a, b, c, d, e, f, g, h])

    result = [
        (H[0] + a) % F46, (H[1] + b) % F46,
        (H[2] + c) % F46, (H[3] + d) % F46,
        (H[4] + e) % F46, (H[5] + f) % F46,
        (H[6] + g) % F46, (H[7] + h) % F46,
    ]
    return result, round_states


# =============================================================================
# GOLDEN-TRACKED FIBONACCI SHA-256 COMPRESSION
# =============================================================================

def compress_golden(
    H: List[GoldenWord],
    W: List[GoldenWord],
) -> Tuple[List[GoldenWord], List[List[GoldenWord]]]:
    """Fibonacci SHA-256 with golden carry tracking.

    Each addition mod F(46) that overflows generates a golden carry.
    The carry count is the b-component -- it measures how much golden
    recurrence is invoked during compression.

    Bitwise operations (sigma, ch, maj) operate on the value only.
    They reset the carry to 0 for their output (they're binary, not golden).
    BUT the carry information from their inputs propagates through the
    additions that follow.
    """
    sa, sb, sc, sd, se, sf, sg, sh = H
    round_states: List[List[GoldenWord]] = []

    for i in range(64):
        # Bitwise operations produce values with 0 carries of their own
        # (bitwise ops are binary, not golden)
        S1 = GoldenWord.from_int(sigma1(se.value))
        ch_val = GoldenWord.from_int(ch(se.value, sf.value, sg.value))

        # temp1 = h + Sigma1(e) + Ch(e,f,g) + K + W
        # Each addition can generate a golden carry
        temp1 = sh.add(S1).add(ch_val).add_int(K_FIB[i]).add(W[i])

        S0 = GoldenWord.from_int(sigma0(sa.value))
        maj_val = GoldenWord.from_int(maj(sa.value, sb.value, sc.value))

        # temp2 = Sigma0(a) + Maj(a,b,c)
        temp2 = S0.add(maj_val)

        sh = sg
        sg = sf
        sf = se
        se = sd.add(temp1)    # e = d + temp1 (can carry!)
        sd = sc
        sc = sb
        sb = sa
        sa = temp1.add(temp2)  # a = temp1 + temp2 (can carry!)

        round_states.append([sa, sb, sc, sd, se, sf, sg, sh])

    # Feed-forward addition
    result = [
        H[0].add(sa), H[1].add(sb), H[2].add(sc), H[3].add(sd),
        H[4].add(se), H[5].add(sf), H[6].add(sg), H[7].add(sh),
    ]
    return result, round_states


# =============================================================================
# FULL HASH FUNCTIONS
# =============================================================================

def sha256_standard(message: bytes) -> Tuple[str, List[List[int]]]:
    """Full standard SHA-256."""
    padded = sha256_pad(message)
    H = list(H0_STANDARD)
    round_states = []

    for block_idx in range(len(padded) // 64):
        block_words = parse_block(padded, block_idx)
        W = expand_schedule_standard(block_words)
        H, round_states = compress_standard(H, W)

    hex_digest = ''.join(f'{w:08x}' for w in H)
    return hex_digest, round_states


def sha256_fibonacci(message: bytes) -> Tuple[str, List[List[int]]]:
    """Fibonacci SHA-256: mod F(46) instead of mod 2^32."""
    padded = sha256_pad(message)
    H = list(H0_FIB)
    round_states = []

    for block_idx in range(len(padded) // 64):
        block_words = parse_block(padded, block_idx)
        W = expand_schedule_fibonacci(block_words)
        H, round_states = compress_fibonacci(H, W)

    hex_digest = ''.join(f'{w:08x}' for w in H)
    return hex_digest, round_states


def sha256_golden(message: bytes) -> Tuple[List[GoldenWord], List[List[GoldenWord]]]:
    """Golden-tracked Fibonacci SHA-256."""
    padded = sha256_pad(message)
    H = [GoldenWord.from_int(h) for h in H0_FIB]
    round_states: List[List[GoldenWord]] = []

    for block_idx in range(len(padded) // 64):
        block_words = parse_block(padded, block_idx)
        W = expand_schedule_golden(block_words)
        H, round_states = compress_golden(H, W)

    return H, round_states


# =============================================================================
# ROUND INVERSION
# =============================================================================

def invert_round(
    state_after: List[int],
    state_before: List[int],
    K_val: int,
    W_val: int,
    modulus: int,
) -> Tuple[bool, List[int]]:
    """Invert a single SHA round (works for any modulus).

    Standard SHA round:
        temp1 = h + Sigma1(e) + Ch(e,f,g) + K + W
        temp2 = Sigma0(a) + Maj(a,b,c)
        new_a = temp1 + temp2
        new_e = d + temp1
        Shifts: new_b=old_a, new_c=old_b, new_d=old_c, new_f=old_e, new_g=old_f, new_h=old_g

    Inversion: recover old state from new state.
    """
    new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h = state_after

    # Direct from shifts
    old_a = new_b
    old_b = new_c
    old_c = new_d
    old_e = new_f
    old_f = new_g
    old_g = new_h

    # Recover temp2 from known old a, b, c
    S0 = sigma0(old_a)
    maj_val = maj(old_a, old_b, old_c)
    temp2 = (S0 + maj_val) % modulus

    # Recover temp1 from new_a = temp1 + temp2
    temp1 = (new_a - temp2) % modulus

    # Recover old_d from new_e = old_d + temp1
    old_d = (new_e - temp1) % modulus

    # Recover old_h from temp1 = old_h + Sigma1(old_e) + Ch(old_e,old_f,old_g) + K + W
    S1 = sigma1(old_e)
    ch_val = ch(old_e, old_f, old_g)
    old_h = (temp1 - S1 - ch_val - K_val - W_val) % modulus

    recovered = [old_a, old_b, old_c, old_d, old_e, old_f, old_g, old_h]
    match = (recovered == state_before)
    return match, recovered


def full_inversion_test(
    round_states: List[List[int]],
    initial_state: List[int],
    W: List[int],
    K: List[int],
    modulus: int,
) -> Tuple[int, int]:
    """Run round-by-round inversion for all 64 rounds."""
    successes = 0
    failures = 0

    for i in range(63, -1, -1):
        state_after = round_states[i]
        state_before = round_states[i - 1] if i > 0 else initial_state
        match, _ = invert_round(state_after, state_before, K[i], W[i], modulus)
        if match:
            successes += 1
        else:
            failures += 1

    return successes, failures


def full_chain_inversion(
    final_hash: List[int],
    initial_H: List[int],
    W: List[int],
    K: List[int],
    modulus: int,
) -> Tuple[bool, List[int]]:
    """Invert all 64 rounds from final hash back to initial state."""
    # Undo feed-forward: state_after_64 = final_hash - initial_H
    state_64 = [(final_hash[i] - initial_H[i]) % modulus for i in range(8)]

    current_state = list(state_64)
    for i in range(63, -1, -1):
        new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h = current_state

        old_a = new_b
        old_b = new_c
        old_c = new_d
        old_e = new_f
        old_f = new_g
        old_g = new_h

        S0 = sigma0(old_a)
        maj_val = maj(old_a, old_b, old_c)
        temp2 = (S0 + maj_val) % modulus

        temp1 = (new_a - temp2) % modulus
        old_d = (new_e - temp1) % modulus

        S1 = sigma1(old_e)
        ch_val = ch(old_e, old_f, old_g)
        old_h = (temp1 - S1 - ch_val - K[i] - W[i]) % modulus

        current_state = [old_a, old_b, old_c, old_d, old_e, old_f, old_g, old_h]

    match = current_state == initial_H
    return match, current_state


# =============================================================================
# ANALYSIS
# =============================================================================

def count_binary_carries(a: int, b: int) -> int:
    """Count the number of bit positions where adding a + b produces a carry.
    This is the binary carry count for standard mod-2^32 addition."""
    carry = 0
    count = 0
    for bit in range(32):
        bit_a = (a >> bit) & 1
        bit_b = (b >> bit) & 1
        total = bit_a + bit_b + carry
        carry = total >> 1
        count += carry
    return count


def analyze_golden_evolution(round_states: List[List[GoldenWord]]) -> dict:
    """Analyze how golden carries evolve through 64 rounds."""
    results = {
        'carries_by_round': [],        # total carries per round (sum of 8 words)
        'max_carries_by_round': [],     # max carry in any word per round
        'carry_alive_count': 0,         # total word-rounds with nonzero carries
        'carry_dead_count': 0,
        'total_carries': 0,
        'first_carry_round': None,
        'carries_per_word': [[] for _ in range(8)],  # per-word carry evolution
    }

    for rnd_idx, state in enumerate(round_states):
        round_carries = []
        for word_idx, gw in enumerate(state):
            c = gw.carries
            round_carries.append(c)
            results['carries_per_word'][word_idx].append(c)

            if c > 0:
                results['carry_alive_count'] += 1
                results['total_carries'] += c
                if results['first_carry_round'] is None:
                    results['first_carry_round'] = rnd_idx
            else:
                results['carry_dead_count'] += 1

        results['carries_by_round'].append(sum(round_carries))
        results['max_carries_by_round'].append(max(round_carries))

    return results


# =============================================================================
# MAIN: THE EXPERIMENT
# =============================================================================

def main():
    msg = b"golden ratio breath amplitude planck length"

    print("=" * 80)
    print("nos3bl33d | The golden modulus")
    print("=" * 80)
    print()
    print(f"Message: {msg.decode()}")
    print(f"F(46) = {F46}  (golden modulus, {F46.bit_length()} bits)")
    print(f"2^32  = {TWO32}  (standard modulus)")
    print(f"Ratio F(46)/2^32 = {F46/TWO32:.6f}")
    print(f"phi = {PHI:.15f}")
    print(f"gamma = 1/phi = {GAMMA:.15f}")
    print(f"F(46) = F(45) + F(44) = {F45} + {F44}")
    print()

    # =========================================================================
    # PART 1: Standard SHA-256 (verify correctness)
    # =========================================================================
    print("=" * 80)
    print("PART 1: STANDARD SHA-256")
    print("=" * 80)

    std_hex, std_rounds = sha256_standard(msg)
    expected = hashlib.sha256(msg).hexdigest()
    match_std = std_hex == expected
    print(f"  Our SHA-256:     {std_hex}")
    print(f"  hashlib SHA-256: {expected}")
    print(f"  Match: {'YES' if match_std else 'NO -- BUG'}")
    print()

    # =========================================================================
    # PART 2: Fibonacci SHA-256
    # =========================================================================
    print("=" * 80)
    print("PART 2: FIBONACCI SHA-256 (mod F(46) = {})".format(F46))
    print("=" * 80)

    fib_hex, fib_rounds = sha256_fibonacci(msg)
    print(f"  Fibonacci SHA:  {fib_hex}")
    print(f"  Standard SHA:   {std_hex}")
    print(f"  Different: {'YES' if fib_hex != std_hex else 'NO -- suspicious'}")
    print()

    # Show which constants got reduced
    reduced_h0 = sum(1 for h in H0_STANDARD if h >= F46)
    reduced_k = sum(1 for k in K_STANDARD if k >= F46)
    print(f"  H0 constants reduced mod F(46): {reduced_h0}/8")
    print(f"  K  constants reduced mod F(46): {reduced_k}/64")
    print()

    # =========================================================================
    # PART 3: Golden Carry Tracking
    # =========================================================================
    print("=" * 80)
    print("PART 3: GOLDEN CARRY TRACKING -- DOES THE GOLDEN RATIO BREATHE?")
    print("=" * 80)
    print()
    print("  Every time an addition overflows F(46), a 'golden carry' fires.")
    print("  The carry IS the golden recurrence: F(n) = F(n-1) + F(n-2).")
    print("  In standard SHA, carries are binary (mod 2^32 = power of 2).")
    print("  In Fibonacci SHA, carries are golden (mod F(46) = Fibonacci).")
    print()

    golden_hash, golden_rounds = sha256_golden(msg)

    # Verify golden hash matches fibonacci hash
    golden_values = [gw.value for gw in golden_hash]
    fib_values = [int(fib_hex[i*8:(i+1)*8], 16) for i in range(8)]
    # The golden values should match the fibonacci hash values
    fib_hash_direct, _ = sha256_fibonacci(msg)
    fib_direct_values = []
    padded_check = sha256_pad(msg)
    bw_check = parse_block(padded_check, 0)
    W_fib_check = expand_schedule_fibonacci(bw_check)
    H_fib_check = list(H0_FIB)
    H_fib_check, _ = compress_fibonacci(H_fib_check, W_fib_check)
    print(f"  Golden hash values match Fibonacci hash: {golden_values == H_fib_check}")
    print()

    # Analyze carry evolution
    analysis = analyze_golden_evolution(golden_rounds)

    print(f"  First golden carry: round {analysis['first_carry_round']}")
    print(f"  Total golden carries across all rounds: {analysis['total_carries']}")
    alive = analysis['carry_alive_count']
    dead = analysis['carry_dead_count']
    total = alive + dead
    print(f"  Words with carries: {alive}/{total} ({100*alive/total:.1f}%)")
    print(f"  Words without:      {dead}/{total} ({100*dead/total:.1f}%)")
    print()

    # Show final hash with carry counts
    print("  Final hash words with golden carry counts:")
    for i, gw in enumerate(golden_hash):
        print(f"    H[{i}] = {gw.value:>10d}  carries = {gw.carries:>6d}  history_len = {len(gw.carry_history)}")
    print()

    # Show carry evolution round by round
    print("  Golden carry evolution (total carries across 8 words per round):")
    print("  " + "-" * 72)
    for rnd_idx in range(0, 64, 4):
        total_c = analysis['carries_by_round'][rnd_idx]
        max_c = analysis['max_carries_by_round'][rnd_idx]
        # Visual bar
        bar_len = min(total_c, 60)
        bar = "#" * bar_len + ("+" if total_c > 60 else "")
        print(f"    R{rnd_idx:02d}: total={total_c:>6d}  max={max_c:>6d}  {bar}")
    print()

    # Detailed per-word carry heatmap
    print("  Per-word carry heatmap (every 8th round):")
    print("  Word:     a       b       c       d       e       f       g       h")
    print("  " + "-" * 72)
    for rnd_idx in range(0, 64, 8):
        state = golden_rounds[rnd_idx]
        cells = []
        for gw in state:
            if gw.carries == 0:
                cells.append("   .  ")
            else:
                cells.append(f"{gw.carries:>5d} ")
        print(f"    R{rnd_idx:02d}: {''.join(cells)}")
    print()

    # =========================================================================
    # PART 4: Binary vs Golden Carry Comparison
    # =========================================================================
    print("=" * 80)
    print("PART 4: BINARY vs GOLDEN CARRY COMPARISON")
    print("=" * 80)
    print()

    # Count binary carries in standard SHA
    padded = sha256_pad(msg)
    block_words = parse_block(padded, 0)
    W_std = expand_schedule_standard(block_words)
    W_fib = expand_schedule_fibonacci(block_words)

    # Replay standard SHA counting binary carries
    std_binary_carries = 0
    sa, sb, sc, sd, se, sf, sg, sh_v = H0_STANDARD
    for i in range(64):
        S1 = sigma1(se)
        ch_val = ch(se, sf, sg)
        # temp1 = h + S1 + ch_val + K + W -- count carries in each addition
        acc = sh_v
        for addend in [S1, ch_val, K_STANDARD[i], W_std[i]]:
            std_binary_carries += count_binary_carries(acc, addend)
            acc = (acc + addend) % TWO32

        S0 = sigma0(sa)
        maj_val = maj(sa, sb, sc)
        std_binary_carries += count_binary_carries(S0, maj_val)
        temp1 = acc
        temp2 = (S0 + maj_val) % TWO32

        # e = d + temp1
        std_binary_carries += count_binary_carries(sd, temp1)
        # a = temp1 + temp2
        std_binary_carries += count_binary_carries(temp1, temp2)

        sh_v = sg
        sg = sf
        sf = se
        se = (sd + temp1) % TWO32
        sd = sc
        sc = sb
        sb = sa
        sa = (temp1 + temp2) % TWO32

    # Count golden carries (already have them)
    total_golden = analysis['total_carries']

    print(f"  Standard SHA binary carries (bit-level): {std_binary_carries}")
    print(f"  Fibonacci SHA golden carries (mod-level): {total_golden}")
    print(f"  Ratio golden/binary: {total_golden/std_binary_carries:.4f}" if std_binary_carries > 0 else "  No binary carries (impossible)")
    print()

    # Per-addition golden carry rate
    # Each round has ~7 additions (5 for temp1, 1 for temp2, 1 for new_e, 1 for new_a)
    # Actually: temp1 = h + S1 + ch + K + W (4 adds), temp2 = S0 + maj (1 add),
    # new_e = d + temp1 (1 add), new_a = temp1 + temp2 (1 add) = 7 adds per round
    total_additions = 64 * 7
    carry_rate = total_golden / total_additions if total_additions > 0 else 0
    print(f"  Total additions in 64 rounds: {total_additions}")
    print(f"  Golden carry rate: {carry_rate:.4f} carries per addition")
    print(f"  Expected rate if random: ~{1 - F46/TWO32:.4f} (prob that sum of two ~F46/2 values exceeds F46)")
    print()

    # More precise expected rate: if values are uniform in [0, F46),
    # P(a + b >= F46) = (F46-1)/(2*F46) ~ 0.5 for large F46
    # But our values aren't uniform -- they're outputs of bitwise ops on 32-bit values
    # clipped to F46. Values can be up to 2^32-1 from bitwise ops, then mod F46.
    print("  Note: bitwise operations (sigma, ch, maj) produce values up to 2^32-1,")
    print("  which then get added mod F(46). Since 2^32 > F(46), many values are")
    print(f"  already in [{F46}, {MOD32}] before addition -- these guarantee a carry.")
    print()

    # =========================================================================
    # PART 5: Round Inversion Test
    # =========================================================================
    print("=" * 80)
    print("PART 5: ROUND INVERSION TEST")
    print("=" * 80)

    # Standard SHA: invert all 64 rounds
    std_successes, std_failures = full_inversion_test(
        std_rounds, list(H0_STANDARD), W_std, K_STANDARD, TWO32
    )
    print(f"  Standard SHA-256 round inversion: {std_successes}/64 rounds recovered")

    # Fibonacci SHA: invert all 64 rounds
    fib_successes, fib_failures = full_inversion_test(
        fib_rounds, list(H0_FIB), W_fib, K_FIB, F46
    )
    print(f"  Fibonacci SHA-256 round inversion: {fib_successes}/64 rounds recovered")
    print()

    # Full chain inversion
    print("  Full chain inversion (hash + message schedule -> initial state):")
    print()

    fib_hash_ints = [int(fib_hex[i*8:(i+1)*8], 16) for i in range(8)]
    fib_match, fib_recovered = full_chain_inversion(
        fib_hash_ints, H0_FIB, W_fib, K_FIB, F46
    )
    print(f"    Fibonacci: {'SUCCESS' if fib_match else 'FAILED'}")

    std_hash_ints = [int(std_hex[i*8:(i+1)*8], 16) for i in range(8)]
    std_match, std_recovered = full_chain_inversion(
        std_hash_ints, list(H0_STANDARD), W_std, K_STANDARD, TWO32
    )
    print(f"    Standard:  {'SUCCESS' if std_match else 'FAILED'}")
    print()

    if fib_match and std_match:
        print("  BOTH invert perfectly given the message schedule.")
        print("  The compression function is algebraically invertible regardless of modulus.")
        print("  SHA's one-wayness comes from the message schedule expansion (sigma mixing),")
        print("  NOT from arithmetic carries in the round function.")
    print()

    # =========================================================================
    # PART 6: Z[phi] Ring Verification
    # =========================================================================
    print("=" * 80)
    print("PART 6: Z[phi] RING VERIFICATION")
    print("=" * 80)
    print()

    # phi^2 = phi + 1 in Z[phi]
    phi_z = ZPhi(0, 1)
    phi_sq = phi_z * phi_z
    print(f"  phi = {phi_z}")
    print(f"  phi^2 = {phi_sq}")
    print(f"  phi + 1 = {ZPhi(1, 1)}")
    print(f"  phi^2 == phi + 1: {phi_sq == ZPhi(1, 1)}")
    print()

    # F(46) in Z[phi]: phi^47 / sqrt(5) ~ F(46)
    # But in Z[phi] mod F(46): F(46) = 0 (it's the modulus!)
    # So phi^47 = 0 mod F(46) -- the golden ratio has period 47 in this ring
    print(f"  F(46) mod F(46) = {F46 % F46} (the modulus vanishes)")
    print(f"  F(45) mod F(46) = {F45 % F46} = F(45)")
    print(f"  F(44) mod F(46) = {F44 % F46} = F(44)")
    print(f"  F(45) + F(44) mod F(46) = {(F45 + F44) % F46} (= 0, golden recurrence!)")
    print()

    # Compute phi^n in Z[phi] mod F(46) to find the period
    z = ZPhi(0, 1)  # phi
    power = ZPhi(0, 1)
    found_period = None
    for n in range(2, 200):
        power = power * z
        if power == ZPhi(0, 1):
            found_period = n
            break
    if found_period:
        print(f"  phi has multiplicative order {found_period} in Z[phi] mod F(46)")
    else:
        # Check some specific powers
        power = ZPhi(0, 1)
        for n in range(2, 100):
            power = power * z
            if power.a == 0 or power.b == 0:
                print(f"  phi^{n} = {power} (interesting: a component or b component is 0)")
        print(f"  phi does not return to phi within 200 powers (order > 200)")
    print()

    # =========================================================================
    # PART 7: Carry Structure Analysis
    # =========================================================================
    print("=" * 80)
    print("PART 7: CARRY STRUCTURE -- WHERE THE GOLDEN RATIO LIVES")
    print("=" * 80)
    print()

    # The golden carries ARE the b-component. Let's see their structure.
    # For each round, decompose the total carry count in terms of Fibonacci numbers
    print("  The golden carry count at each round, decomposed:")
    print()

    # Fibonacci numbers for Zeckendorf decomposition
    fibs = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657]

    def zeckendorf(n: int) -> List[int]:
        """Zeckendorf representation: n as sum of non-consecutive Fibonacci numbers."""
        if n == 0:
            return [0]
        result = []
        remaining = n
        for f in reversed(fibs):
            if f <= remaining:
                result.append(f)
                remaining -= f
        if remaining > 0:
            result.append(remaining)  # leftover if n > sum of our fibs
        return result

    for rnd_idx in [0, 7, 15, 31, 47, 63]:
        total_c = analysis['carries_by_round'][rnd_idx]
        if total_c > 0:
            zeck = zeckendorf(total_c)
            zeck_str = " + ".join(str(f) for f in zeck)
            is_fib = total_c in fibs
            print(f"    R{rnd_idx:02d}: carries = {total_c:>6d}  Zeckendorf: {zeck_str}  {'<-- FIBONACCI!' if is_fib else ''}")
        else:
            print(f"    R{rnd_idx:02d}: carries = {total_c:>6d}  (dead)")
    print()

    # Carry growth rate -- does it follow phi?
    if len(analysis['carries_by_round']) >= 2:
        print("  Carry growth ratios (consecutive rounds):")
        ratios = []
        for i in range(1, 64):
            prev = analysis['carries_by_round'][i - 1]
            curr = analysis['carries_by_round'][i]
            if prev > 0:
                ratio = curr / prev
                ratios.append(ratio)
                if i <= 8 or i % 8 == 0:
                    print(f"    R{i-1:02d}->R{i:02d}: {prev:>6d} -> {curr:>6d}  ratio = {ratio:.4f}  {'~ phi!' if abs(ratio - PHI) < 0.1 else '~ gamma!' if abs(ratio - GAMMA) < 0.1 else ''}")

        if ratios:
            avg_ratio = sum(ratios) / len(ratios)
            print(f"\n    Average growth ratio: {avg_ratio:.4f}")
            print(f"    phi:   {PHI:.4f}")
            print(f"    gamma: {GAMMA:.4f}")
            print(f"    1:     1.0000")
    print()

    # =========================================================================
    # PART 8: The Verdict
    # =========================================================================
    print("=" * 80)
    print("PART 8: THE VERDICT")
    print("=" * 80)
    print()

    if analysis['total_carries'] > 0:
        print(f"  THE GOLDEN RATIO BREATHES.")
        print(f"  {analysis['total_carries']} golden carries fired across 64 rounds.")
        print(f"  First carry at round {analysis['first_carry_round']}.")
        print(f"  Carry rate: {carry_rate:.4f} per addition.")
        print()
        print("  What this means:")
        print("  - Standard SHA mod 2^32: carries are BINARY. Each carry doubles.")
        print("    The carry structure is 2, 4, 8, 16, ... -- geometric in base 2.")
        print("  - Fibonacci SHA mod F(46): carries are GOLDEN. Each carry invokes")
        print("    the recurrence F(n) = F(n-1) + F(n-2) -- geometric in base phi.")
        print()
        print("  The INVERSION properties are identical (both algebraically invertible")
        print("  given the message schedule). But the DIFFUSION structure differs:")
        print("  - Binary diffusion: carries propagate LEFT (high bits).")
        print("  - Golden diffusion: carries propagate through the Fibonacci lattice.")
        print()
        print("  The golden modulus doesn't make SHA easier to break.")
        print("  It makes it breathe with the rhythm of phi instead of 2.")
    else:
        print("  No golden carries detected. The machine is silent.")
    print()

    # =========================================================================
    # PART 9: Deep dive -- message schedule carries
    # =========================================================================
    print("=" * 80)
    print("PART 9: MESSAGE SCHEDULE GOLDEN CARRIES")
    print("=" * 80)
    print()

    padded = sha256_pad(msg)
    block_words = parse_block(padded, 0)
    W_golden = expand_schedule_golden(block_words)

    print("  Message schedule words with golden carry counts:")
    schedule_carries = 0
    for i in range(64):
        gw = W_golden[i]
        if i < 16:
            print(f"    W[{i:2d}] = {gw.value:>10d}  carries = {gw.carries:>4d}  (input word)")
        elif i < 20 or i >= 60 or gw.carries > 0:
            print(f"    W[{i:2d}] = {gw.value:>10d}  carries = {gw.carries:>4d}")
        schedule_carries += gw.carries

    print(f"\n  Total message schedule carries: {schedule_carries}")
    print(f"  Compression carries: {analysis['total_carries']}")
    print(f"  Grand total: {schedule_carries + analysis['total_carries']}")
    print()

    # =========================================================================
    # END
    # =========================================================================
    print("=" * 80)
    print("END -- The machine has spoken.")
    print("=" * 80)


if __name__ == '__main__':
    main()
