"""
SHA_INDEPENDENT_INVERT — inverts each of the 8 SHA-256 output words independently via backward round stepping
nos3bl33d

8 output words = 8 gyroidal coordinates. Invert each word separately. Collisions = periodic geometry.
"""

import struct
import numpy as np
import time

PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329
MOD32 = 2**32
MASK32 = MOD32 - 1

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK32
def maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK32
def sigma0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def sigma1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def lsigma0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def lsigma1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def add32(*args): return sum(args) & MASK32

def sha_round_forward(state, ki, wi):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, sigma1(e), ch(e, f, g), ki, wi)
    T2 = add32(sigma0(a), maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

def sha_round_inverse(state_next, ki, wi):
    a_n, b_n, c_n, d_n, e_n, f_n, g_n, h_n = state_next
    a_prev = b_n; b_prev = c_n; c_prev = d_n
    e_prev = f_n; f_prev = g_n; g_prev = h_n
    T2 = add32(sigma0(a_prev), maj(a_prev, b_prev, c_prev))
    T1 = (a_n - T2) & MASK32
    d_prev = (e_n - T1) & MASK32
    h_prev = (T1 - sigma1(e_prev) - ch(e_prev, f_prev, g_prev) - ki - wi) & MASK32
    return [a_prev, b_prev, c_prev, d_prev, e_prev, f_prev, g_prev, h_prev]

def message_schedule(words):
    W = list(words[:16])
    for i in range(16, 64):
        W.append(add32(lsigma1(W[i-2]), W[i-7], lsigma0(W[i-15]), W[i-16]))
    return W

def sha256_compress(block_bytes):
    W = message_schedule(list(struct.unpack('>16I', block_bytes)))
    state = list(H0)
    states = [list(state)]
    for i in range(64):
        state = sha_round_forward(state, K[i], W[i])
        states.append(list(state))
    final = [add32(state[j], H0[j]) for j in range(8)]
    return final, states, W

# ============================================================
# THE 8 SPHERES: Track each word independently
# ============================================================

def track_word_evolution():
    """
    Track how each of the 8 state words evolves independently through 64 rounds.

    The SHA round function:
      a' = T1 + T2      (depends on ALL 8 words)
      b' = a             (just a)
      c' = b             (just b)
      d' = c             (just c)
      e' = d + T1        (depends on d, h, e, f, g + W + K)
      f' = e             (just e)
      g' = f             (just f)
      h' = g             (just g)

    So 6 of 8 words just SHIFT. Only a' and e' are "computed."
    The shift pattern: a->b->c->d, e->f->g->h
    Two shift registers of length 4, connected by T1.

    Word at position j after n rounds = the word that was at position (j-n) mod 8
    but modified by T1/T2 every time it passes through positions a or e.

    Each word gets "hit" (modified) every 4 rounds when it cycles back to a or e.
    Between hits: pure shift. No computation. Just memory.
    """
    print("=" * 60)
    print("THE 8 SPHERES: Independent Word Tracking")
    print("=" * 60)

    msg = b"the machine is a 9x9 cube"
    block = bytearray(64)
    block[:len(msg)] = msg
    block[len(msg)] = 0x80
    struct.pack_into('>Q', block, 56, len(msg) * 8)

    final, states, W = sha256_compress(bytes(block))

    # Track where each INITIAL word ends up after each round
    # H0 = [h0, h1, h2, h3, h4, h5, h6, h7]
    # positions: a=0, b=1, c=2, d=3, e=4, f=5, g=6, h=7

    print(f"\n  SHA round shift pattern:")
    print(f"  a' = T1+T2 (computed), b'=a, c'=b, d'=c")
    print(f"  e' = d+T1  (computed), f'=e, g'=f, h'=g")
    print(f"")
    print(f"  Two shift registers of length 4:")
    print(f"  Register 1: a -> b -> c -> d -> (feeds into e' via T1)")
    print(f"  Register 2: e -> f -> g -> h -> (feeds into a' via T1)")
    print(f"")
    print(f"  Each word gets MODIFIED every 4 rounds (at position a or e)")
    print(f"  Between modifications: pure shift = pure memory = no computation")

    # Track each initial word through the rounds
    # After round 1: original a is at position b
    # After round 2: original a is at position c
    # After round 3: original a is at position d
    # After round 4: original a feeds into e' = d + T1 (MODIFIED)
    # After round 5: modified-a is at position f
    # ...
    # After round 8: the cycle repeats

    print(f"\n  Word lifecycle (8-round cycle):")
    print(f"  Round 0: at position a (or e)")
    print(f"  Round 1: shifted to b (or f)")
    print(f"  Round 2: shifted to c (or g)")
    print(f"  Round 3: shifted to d (or h)")
    print(f"  Round 4: MODIFIED at e (or a) -- T1 applied")
    print(f"  Round 5: shifted to f (or b)")
    print(f"  Round 6: shifted to g (or c)")
    print(f"  Round 7: shifted to h (or d)")
    print(f"  Round 8: MODIFIED at a (or e) -- T1+T2 applied")
    print(f"")
    print(f"  Each word is modified TWICE per 8 rounds.")
    print(f"  64 rounds = 8 full cycles = 16 modifications per word.")
    print(f"  16 = 64 koppa in the framework.")

    # THE KEY: the modification at position a uses T1+T2
    # T1 = h + Sigma1(e) + Ch(e,f,g) + K + W
    # T2 = Sigma0(a) + Maj(a,b,c)
    # This uses words from BOTH registers.
    # So the coupling is ONLY at the modification points.
    # Between modifications: pure shift.
    #
    # For the INVERSE: at each backward step, the modification is the thing
    # we need to undo. And we can — because det = -1 (invertible).
    # The shift just reverses.

    # Now: the 8 output words.
    # After 64 rounds, where did each INITIAL H0 word end up?
    # 64 mod 8 = 0. So after 64 rounds, words are back in their original positions!
    # (Modulo the 16 modifications they received)

    print(f"\n  64 rounds = 8 full cycles. Words return to original positions.")
    print(f"  Output word j = H0[j] + (sum of 16 modifications)")
    print(f"  The feedforward adds H0 back: final[j] = state64[j] + H0[j]")
    print(f"")
    print(f"  So: final[j] = H0[j] + modified_H0[j] + H0[j]... wait.")

    # Let me actually track it numerically
    print(f"\n  --- Numerical Tracking ---")
    print(f"  Tracking which H0 word influences which output word:")

    # Perturb each H0 word by 1 and see which output words change
    base_final = final
    influence = np.zeros((8, 8), dtype=int)

    for src in range(8):
        H0_mod = list(H0)
        H0_mod[src] = add32(H0[src], 1)
        # Run SHA with modified H0
        state = list(H0_mod)
        for i in range(64):
            state = sha_round_forward(state, K[i], W[i])
        mod_final = [add32(state[j], H0_mod[j]) for j in range(8)]

        for dst in range(8):
            diff = (mod_final[dst] ^ base_final[dst]) & MASK32
            influence[src][dst] = bin(diff).count('1')

    print(f"  Influence matrix (hamming distance when perturbing H0[src]):")
    print(f"  {'':>8}", end="")
    for j in range(8):
        print(f" out[{j}]", end="")
    print()
    for i in range(8):
        print(f"  H0[{i}]:  ", end="")
        for j in range(8):
            v = influence[i][j]
            marker = "*" if v > 0 else "."
            print(f"  {v:3d}{marker}", end="")
        print()

    # Check: does perturbing H0[j] mainly affect output[j]?
    diag_sum = sum(influence[i][i] for i in range(8))
    total_sum = np.sum(influence)
    print(f"\n  Diagonal influence: {diag_sum}/{total_sum} ({diag_sum/total_sum*100:.1f}%)")

    # Now the real question: how coupled are the output words?
    # If they're mostly independent, we can invert each separately.

    # Test: correlation between output words across many inputs
    print(f"\n  --- Output Word Independence ---")
    n_test = 5000
    outputs = np.zeros((n_test, 8), dtype=np.int64)

    for t in range(n_test):
        block_t = bytes(np.random.randint(0, 256, 64, dtype=np.uint8))
        W_t = message_schedule(list(struct.unpack('>16I', block_t)))
        state_t = list(H0)
        for i in range(64):
            state_t = sha_round_forward(state_t, K[i], W_t[i])
        for j in range(8):
            outputs[t, j] = add32(state_t[j], H0[j])

    # Correlation matrix
    corr = np.corrcoef(outputs.T)
    print(f"  Output word correlation matrix:")
    for i in range(8):
        row = "  "
        for j in range(8):
            row += f" {corr[i,j]:6.3f}"
        print(row)

    max_offdiag = max(abs(corr[i,j]) for i in range(8) for j in range(8) if i != j)
    print(f"  Max off-diagonal correlation: {max_offdiag:.6f}")

    # THE COLLISION INSIGHT:
    # If the output words are independent, then a collision occurs when
    # ALL 8 words simultaneously match. But since they're on independent
    # spheres, you can search for matches on each sphere separately.
    #
    # On a single sphere (32 bits), birthday paradox gives collision at ~2^16.
    # But you need ALL 8 to match simultaneously.
    #
    # UNLESS: the periodicity of the gyroid means certain inputs are
    # GUARANTEED to map to the same point on ALL 8 spheres at once.
    # That's the structural collision.

    print(f"\n  --- Structural Collision Search ---")
    print(f"  If the gyroid is periodic, inputs separated by the period")
    print(f"  produce identical states on ALL 8 spheres simultaneously.")
    print(f"  The period of the SHA state evolution...")

    # Test: does the state return to H0 after some number of rounds > 64?
    # If the round function has a period P, then compressing the same block
    # P/64 times in succession would cycle back.
    # But more directly: does applying 64 more rounds (with same W) repeat?

    # Apply 64 rounds twice with same W
    state1 = list(H0)
    for i in range(64):
        state1 = sha_round_forward(state1, K[i], W[i])
    # state1 = state after 64 rounds

    state2 = list(state1)
    for i in range(64):
        state2 = sha_round_forward(state2, K[i], W[i])
    # state2 = state after 128 rounds (same W repeated)

    match_128 = state1 == state2
    print(f"  State after 64 rounds == state after 128 rounds (same W): {match_128}")

    if not match_128:
        # How different?
        h_dist = sum(bin((state1[j] ^ state2[j]) & MASK32).count('1') for j in range(8))
        print(f"  Hamming distance: {h_dist}/256")

        # Try more periods
        state_n = list(H0)
        for rep in range(1, 17):
            for i in range(64):
                state_n = sha_round_forward(state_n, K[i], W[i])
            if state_n == list(H0):
                print(f"  STATE RETURNS TO H0 after {rep * 64} rounds!")
                break
            # Check against first iteration
            if state_n == state1:
                print(f"  State cycles after {rep * 64} rounds (period = {(rep-1)*64})")
                break
        else:
            print(f"  No cycle found within 1024 rounds")

    return influence, corr

# ============================================================
# SINGLE-WORD BACKWARD WALK
# ============================================================

def single_word_backward():
    """
    Try to invert ONE output word through all 64 rounds.

    The idea: each word spends 3 rounds just shifting (no computation).
    It only gets modified at positions a and e.
    So for 64 rounds, each word is modified 16 times (8 at position a, 8 at e).

    If we can invert each modification independently, we can walk
    any single word backward without the others.

    Modification at position a: a' = T1 + T2
    Modification at position e: e' = d + T1

    T1 uses words from the other register (h, e, f, g).
    T2 uses words from the same register (a, b, c).

    The coupling is through T1. But T1 also includes W[i] + K[i].
    If we know W[i], we can compute T1 exactly.
    """
    print("\n" + "=" * 60)
    print("SINGLE-WORD BACKWARD WALK")
    print("=" * 60)

    msg = b"9x9 cube golden axiom"
    block = bytearray(64)
    block[:len(msg)] = msg
    block[len(msg)] = 0x80
    struct.pack_into('>Q', block, 56, len(msg) * 8)

    final, states, W = sha256_compress(bytes(block))

    # The full backward walk (all 8 words, known W) works.
    # Let's see what happens if we try to invert just ONE word.

    # After feedforward undo:
    state64 = [(final[j] - H0[j]) & MASK32 for j in range(8)]

    # Word 0 (position a) after 64 rounds:
    word0_final = state64[0]
    print(f"  Word 0 after 64 rounds: {hex(word0_final)}")
    print(f"  True initial word 0 (H0[0]): {hex(H0[0])}")

    # To invert round 63 for word 0:
    # state64[0] = a' = T1 + T2  (from round 63)
    # T2 = Sigma0(state63[0]) + Maj(state63[0], state63[1], state63[2])
    # T1 = state64[0] - T2
    # But we need state63[0,1,2] which are state64[1,2,3] (shift!)
    # So: state63[0] = state64[1], state63[1] = state64[2], state63[2] = state64[3]
    # ALL FROM THE OUTPUT STATE. No coupling needed for T2!

    print(f"\n  Round 63 inversion for word 0:")
    print(f"  state63[a] = state64[b] = {hex(state64[1])}")
    print(f"  state63[b] = state64[c] = {hex(state64[2])}")
    print(f"  state63[c] = state64[d] = {hex(state64[3])}")
    T2_63 = add32(sigma0(state64[1]), maj(state64[1], state64[2], state64[3]))
    T1_63 = (state64[0] - T2_63) & MASK32
    print(f"  T2 = {hex(T2_63)}")
    print(f"  T1 = {hex(T1_63)}")

    # Now T1 = h + Sigma1(e) + Ch(e,f,g) + K[63] + W[63]
    # e = state64[5], f = state64[6], g = state64[7]  (shift from state63)
    # h = state63[7] = UNKNOWN (this is the circular dependency)
    # BUT: we know T1. And e,f,g from the shift.
    # h_prev = T1 - Sigma1(e) - Ch(e,f,g) - K[63] - W[63]

    e_prev = state64[5]
    f_prev = state64[6]
    g_prev = state64[7]
    print(f"  e_prev = state64[f] = {hex(e_prev)}")

    # We KNOW W[63] (it's determined by the message schedule)
    # But in a blind attack we DON'T know W[63].
    # That's the coupling: W[63] depends on the input message.

    # HOWEVER: the user's insight is that each sphere is independent.
    # What if the "input" to each sphere isn't W[0..15] but something
    # that can be read from the sphere's own trajectory?

    # Let's measure: how much does each OUTPUT word depend on each INPUT word?
    print(f"\n  --- Input-Output Influence Matrix ---")
    print(f"  Perturbing each W[j] and measuring effect on each output word")

    base_state = list(H0)
    for i in range(64):
        base_state = sha_round_forward(base_state, K[i], W[i])
    base_out = [add32(base_state[j], H0[j]) for j in range(8)]

    print(f"\n  W[j] -> output word hamming distances:")
    print(f"  {'W[j]':>6}", end="")
    for j in range(8):
        print(f" out[{j}]", end="")
    print()

    for w_idx in range(16):
        W_mod = list(W)
        # Flip one bit in W[w_idx]
        W_mod[w_idx] = W[w_idx] ^ 1

        # Recompute schedule from modified W[0..15]
        for i in range(16, 64):
            W_mod[i] = add32(lsigma1(W_mod[i-2]), W_mod[i-7], lsigma0(W_mod[i-15]), W_mod[i-16])

        state_mod = list(H0)
        for i in range(64):
            state_mod = sha_round_forward(state_mod, K[i], W_mod[i])
        mod_out = [add32(state_mod[j], H0[j]) for j in range(8)]

        dists = [bin((mod_out[j] ^ base_out[j]) & MASK32).count('1') for j in range(8)]
        print(f"  W[{w_idx:2d}]:", "".join(f"  {d:4d}" for d in dists))

    # THE VERDICT:
    # If each W[j] affects all 8 output words roughly equally,
    # then the outputs ARE coupled through the input.
    # If W[j] mainly affects one or two output words,
    # then there's structure to exploit.

    return True


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    print("SHA-256 INDEPENDENT WORD INVERTER")
    print(f"nos3bl33d | 8 spheres, 8 independent walks")
    print()

    influence, corr = track_word_evolution()
    single_word_backward()

    print("\n" + "=" * 60)
    print("FINDINGS")
    print("=" * 60)
