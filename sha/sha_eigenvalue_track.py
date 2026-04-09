#!/usr/bin/env python3
"""
SHA_EIGENVALUE_TRACK — tracks Jacobian eigenvalue evolution through 64 nonlinear SHA-256 rounds
nos3bl33d

Five tests: Jacobian tracking, chained effective matrix, Lyapunov exponents,
SVD condition growth, golden trace per round.
"""

import numpy as np
from numpy.linalg import eigvals, svd, norm
import sys
import time

# ─── SHA-256 Constants ──────────────────────────────────────────────────────

MOD = 2**32

# Initial hash values (fractional parts of sqrt of first 8 primes)
H_INIT = [
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
]

# Round constants (fractional parts of cube roots of first 64 primes)
K = [
    0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5,
    0x3956C25B, 0x59F111F1, 0x923F82A4, 0xAB1C5ED5,
    0xD807AA98, 0x12835B01, 0x243185BE, 0x550C7DC3,
    0x72BE5D74, 0x80DEB1FE, 0x9BDC06A7, 0xC19BF174,
    0xE49B69C1, 0xEFBE4786, 0x0FC19DC6, 0x240CA1CC,
    0x2DE92C6F, 0x4A7484AA, 0x5CB0A9DC, 0x76F988DA,
    0x983E5152, 0xA831C66D, 0xB00327C8, 0xBF597FC7,
    0xC6E00BF3, 0xD5A79147, 0x06CA6351, 0x14292967,
    0x27B70A85, 0x2E1B2138, 0x4D2C6DFC, 0x53380D13,
    0x650A7354, 0x766A0ABB, 0x81C2C92E, 0x92722C85,
    0xA2BFE8A1, 0xA81A664B, 0xC24B8B70, 0xC76C51A3,
    0xD192E819, 0xD6990624, 0xF40E3585, 0x106AA070,
    0x19A4C116, 0x1E376C08, 0x2748774C, 0x34B0BCB5,
    0x391C0CB3, 0x4ED8AA4A, 0x5B9CCA4F, 0x682E6FF3,
    0x748F82EE, 0x78A5636F, 0x84C87814, 0x8CC70208,
    0x90BEFFFA, 0xA4506CEB, 0xBEF9A3F7, 0xC67178F2,
]

PHI = (1 + 5**0.5) / 2  # 1.6180339887...
LOG_PHI = np.log(PHI)     # 0.48121182...

# ─── SHA-256 Primitives ────────────────────────────────────────────────────

def rotr(x, n):
    """Right rotate 32-bit integer x by n bits."""
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def shr(x, n):
    """Right shift 32-bit integer."""
    return (x >> n) & 0xFFFFFFFF

def add32(*args):
    """Add any number of 32-bit integers mod 2^32."""
    s = 0
    for a in args:
        s = (s + a) & 0xFFFFFFFF
    return s

def ch(e, f, g):
    """SHA-256 Ch function: (e AND f) XOR (NOT e AND g)."""
    return (e & f) ^ ((~e) & g) & 0xFFFFFFFF

def maj(a, b, c):
    """SHA-256 Maj function: (a AND b) XOR (a AND c) XOR (b AND c)."""
    return (a & b) ^ (a & c) ^ (b & c)

def sigma0(a):
    """SHA-256 Sigma0: ROTR(2) XOR ROTR(13) XOR ROTR(22)."""
    return rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

def sigma1(e):
    """SHA-256 Sigma1: ROTR(6) XOR ROTR(11) XOR ROTR(25)."""
    return rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

def small_sigma0(x):
    """SHA-256 message schedule sigma0: ROTR(7) XOR ROTR(18) XOR SHR(3)."""
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def small_sigma1(x):
    """SHA-256 message schedule sigma1: ROTR(17) XOR ROTR(19) XOR SHR(10)."""
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)


# ─── Message Schedule ──────────────────────────────────────────────────────

def compute_message_schedule(message_block):
    """
    Expand 16-word message block into 64-word schedule.
    message_block: list of 16 uint32 words.
    Returns: list of 64 uint32 words.
    """
    W = list(message_block[:16])
    for i in range(16, 64):
        w = add32(small_sigma1(W[i-2]), W[i-7], small_sigma0(W[i-15]), W[i-16])
        W.append(w)
    return W


# ─── Single SHA-256 Round ──────────────────────────────────────────────────

def sha_round(state, k_val, w_val):
    """
    Execute one SHA-256 compression round.
    state: [a, b, c, d, e, f, g, h] as list of 8 uint32.
    k_val: round constant K[i].
    w_val: message schedule word W[i].
    Returns: new state [a', b', c', d', e', f', g', h'].
    """
    a, b, c, d, e, f, g, h = state

    t1 = add32(h, sigma1(e), ch(e, f, g), k_val, w_val)
    t2 = add32(sigma0(a), maj(a, b, c))

    h_new = g
    g_new = f
    f_new = e
    e_new = add32(d, t1)
    d_new = c
    c_new = b
    b_new = a
    a_new = add32(t1, t2)

    return [a_new, b_new, c_new, d_new, e_new, f_new, g_new, h_new]


def run_n_rounds(state, W, n_rounds):
    """
    Run n_rounds of SHA-256 compression from initial state.
    Returns final state after n_rounds.
    """
    s = list(state)
    for i in range(min(n_rounds, 64)):
        s = sha_round(s, K[i], W[i])
    return s


def run_n_rounds_trace(state, W, n_rounds):
    """
    Run n_rounds and return ALL intermediate states.
    Returns list of (n_rounds+1) states: [state_0, state_1, ..., state_n].
    """
    trace = [list(state)]
    s = list(state)
    for i in range(min(n_rounds, 64)):
        s = sha_round(s, K[i], W[i])
        trace.append(list(s))
    return trace


# ─── Jacobian Computation ──────────────────────────────────────────────────

def compute_jacobian(state, W, n_rounds, epsilon=1):
    """
    Compute the 8x8 Jacobian of the n-round SHA transformation
    at the given state, using finite differences with perturbation epsilon.

    J[i][j] = (f_i(x + eps*e_j) - f_i(x)) / eps

    Since we work mod 2^32, the differences are computed as signed
    (interpreting the uint32 difference as a signed 32-bit value).
    """
    base_output = run_n_rounds(state, W, n_rounds)

    J = np.zeros((8, 8), dtype=np.float64)

    for j in range(8):
        # Perturb word j by +epsilon
        perturbed = list(state)
        perturbed[j] = (perturbed[j] + epsilon) & 0xFFFFFFFF
        perturbed_output = run_n_rounds(perturbed, W, n_rounds)

        for i in range(8):
            diff = (perturbed_output[i] - base_output[i]) & 0xFFFFFFFF
            # Interpret as signed 32-bit
            if diff >= 0x80000000:
                diff -= MOD
            J[i][j] = diff / epsilon

    return J


def compute_local_jacobian(state, round_idx, W, epsilon=1):
    """
    Compute the 8x8 Jacobian of a SINGLE round at the given state.
    This is the local linearization at that specific state.
    """
    k_val = K[round_idx]
    w_val = W[round_idx]
    base_output = sha_round(state, k_val, w_val)

    J = np.zeros((8, 8), dtype=np.float64)

    for j in range(8):
        perturbed = list(state)
        perturbed[j] = (perturbed[j] + epsilon) & 0xFFFFFFFF
        perturbed_output = sha_round(perturbed, k_val, w_val)

        for i in range(8):
            diff = (perturbed_output[i] - base_output[i]) & 0xFFFFFFFF
            if diff >= 0x80000000:
                diff -= MOD
            J[i][j] = diff / epsilon

    return J


# ─── Input Generators ──────────────────────────────────────────────────────

def golden_state():
    """State based on golden ratio: int(phi * 2^32) for each word (shifted)."""
    base = int(PHI * (2**32)) & 0xFFFFFFFF  # 0x9E3779B9
    return [(base * (i + 1)) & 0xFFFFFFFF for i in range(8)]

def fibonacci_state(n=8):
    """State from Fibonacci numbers."""
    fibs = [1, 1]
    while len(fibs) < n + 10:
        fibs.append(fibs[-1] + fibs[-2])
    # Use Fibonacci numbers scaled into 32-bit range
    return [(fibs[i + 5] * 0x9E3779B9) & 0xFFFFFFFF for i in range(8)]

def phi_spaced_state():
    """State with phi-spaced values: floor(i * phi * 2^31)."""
    return [int(i * PHI * (2**31)) & 0xFFFFFFFF for i in range(1, 9)]

def random_state(seed=42):
    """Deterministic random state."""
    rng = np.random.RandomState(seed)
    # Use two 16-bit halves to avoid int32 overflow
    hi = rng.randint(0, 0x10000, size=8).astype(np.uint32)
    lo = rng.randint(0, 0x10000, size=8).astype(np.uint32)
    return [int((h << 16) | l) for h, l in zip(hi, lo)]

def zero_message():
    """All-zero 512-bit message (16 zero words)."""
    return compute_message_schedule([0] * 16)

def golden_message():
    """Message with golden-structured words."""
    base = int(PHI * (2**32)) & 0xFFFFFFFF
    msg = [(base * (i + 1)) & 0xFFFFFFFF for i in range(16)]
    return compute_message_schedule(msg)


# ─── Formatting Helpers ────────────────────────────────────────────────────

def banner(title):
    sep = "=" * 78
    print(f"\n{sep}")
    print(f"  {title}")
    print(f"{sep}")

def sub_banner(title):
    print(f"\n--- {title} ---")


# ─── TEST 1: Jacobian Eigenvalue Tracking ──────────────────────────────────

def test1_jacobian_tracking():
    banner("TEST 1: Jacobian Eigenvalue Tracking (N=1..20 rounds)")
    print(f"  Input: golden state [int(phi*2^32) * i for i=1..8]")
    print(f"  Message: all zeros")
    print(f"  phi = {PHI:.10f}")
    print(f"  Linear dominant eigenvalue = 1.647 (from sha_matrix.py)")

    state = golden_state()
    W = zero_message()

    print(f"\n  Initial state: {[f'0x{w:08X}' for w in state]}")

    results = []

    print(f"\n  {'N':>3s}  {'|lam_max|':>12s}  {'lam_max':>20s}  {'|lam_max|/phi':>12s}  {'|lam_max|/1.647':>14s}")
    print(f"  {'---':>3s}  {'--------':>12s}  {'-------':>20s}  {'--------':>12s}  {'---------':>14s}")

    for n in range(1, 21):
        J = compute_jacobian(state, W, n)
        eigs = eigvals(J)

        # Sort by magnitude
        idx = np.argsort(-np.abs(eigs))
        eigs_sorted = eigs[idx]
        lam_max = eigs_sorted[0]
        mag = abs(lam_max)

        ratio_phi = mag / PHI if mag > 0 else 0
        ratio_lin = mag / 1.647 if mag > 0 else 0

        if abs(lam_max.imag) < 1e-6:
            lam_str = f"{lam_max.real:>20.6f}"
        else:
            lam_str = f"{lam_max.real:>10.4f}{lam_max.imag:>+10.4f}i"

        print(f"  {n:3d}  {mag:12.6f}  {lam_str}  {ratio_phi:12.6f}  {ratio_lin:14.6f}")

        results.append({
            'n': n,
            'eigenvalues': eigs_sorted,
            'lam_max': lam_max,
            'mag_max': mag,
        })

    # Summary analysis
    sub_banner("Eigenvalue Trajectory Analysis")
    mags = [r['mag_max'] for r in results]

    # Check if magnitude grows as phi^N
    print(f"\n  Growth rate analysis (log|lam_max| / N):")
    for r in results:
        n = r['n']
        if r['mag_max'] > 0:
            rate = np.log(r['mag_max']) / n
            print(f"    N={n:2d}: log(|lam_max|)/N = {rate:.6f}  "
                  f"(log(phi)={LOG_PHI:.6f}, log(1.647)={np.log(1.647):.6f})")

    # Check eigenvalue spectrum at round 1 vs round 4 vs round 8
    sub_banner("Full Eigenvalue Spectrum at Key Rounds")
    for target_n in [1, 2, 4, 8, 12, 16, 20]:
        r = results[target_n - 1]
        eigs = r['eigenvalues']
        print(f"\n  N={target_n}: eigenvalues (sorted by |.|):")
        for i, e in enumerate(eigs):
            mag = abs(e)
            phi_ratio = mag / PHI if mag > 1e-10 else 0
            if abs(e.imag) < 1e-6:
                print(f"    [{i}] {e.real:>14.6f}  |.| = {mag:12.6f}  |.|/phi = {phi_ratio:.6f}")
            else:
                print(f"    [{i}] {e.real:>9.4f}{e.imag:>+9.4f}i  |.| = {mag:12.6f}  |.|/phi = {phi_ratio:.6f}")

    return results


# ─── TEST 2: Effective Matrix from Chained Local Jacobians ─────────────────

def test2_effective_matrix():
    banner("TEST 2: Effective Matrix -- Lyapunov Spectrum via QR (64 rounds)")
    print(f"  Using QR decomposition at each step to prevent overflow.")
    print(f"  Accumulates log(R_diag) to get Lyapunov exponents.")
    print(f"  Each J_i is the local Jacobian at the actual nonlinear state.")

    state = golden_state()
    W = zero_message()

    # Run full trace to get all intermediate states
    trace = run_n_rounds_trace(state, W, 64)

    results = []

    # QR-based Lyapunov exponent computation:
    # Start with Q = I. At each step:
    #   M = J_local @ Q
    #   Q, R = QR(M)
    #   lyap_sum += log(|R_diag|)
    # Lyapunov exponents = lyap_sum / N

    Q = np.eye(8)
    lyap_sum = np.zeros(8)

    print(f"\n  {'N':>3s}  {'lyap_1':>12s}  {'lyap_2':>12s}  {'lyap_3':>12s}  {'lyap_1/logphi':>13s}  {'lyap_1-lyap_2':>13s}")
    print(f"  {'---':>3s}  {'------':>12s}  {'------':>12s}  {'------':>12s}  {'---------':>13s}  {'---------':>13s}")

    for round_idx in range(64):
        J_local = compute_local_jacobian(trace[round_idx], round_idx, W)
        M = J_local @ Q
        Q, R = np.linalg.qr(M)

        # Accumulate log of absolute diagonal of R
        r_diag = np.abs(np.diag(R))
        # Prevent log(0)
        r_diag = np.maximum(r_diag, 1e-300)
        lyap_sum += np.log(r_diag)

        n = round_idx + 1
        lyap_exponents = lyap_sum / n  # Current Lyapunov exponents
        # Sort descending
        lyap_sorted = np.sort(lyap_exponents)[::-1]

        results.append({
            'n': n,
            'lyapunov_exponents': lyap_sorted.copy(),
            'lyap_sum': lyap_sum.copy(),
        })

        if n <= 20 or n % 4 == 0 or n == 64:
            l1 = lyap_sorted[0]
            l2 = lyap_sorted[1]
            l3 = lyap_sorted[2]
            phi_ratio = l1 / LOG_PHI if abs(LOG_PHI) > 0 else 0
            gap = l1 - l2
            print(f"  {n:3d}  {l1:12.6f}  {l2:12.6f}  {l3:12.6f}  {phi_ratio:13.6f}  {gap:13.6f}")

    # Summary: does the maximal Lyapunov exponent converge?
    sub_banner("Lyapunov Spectrum Convergence (last 10 rounds)")
    print(f"  log(phi)   = {LOG_PHI:.6f}")
    print(f"  log(1.647) = {np.log(1.647):.6f}")
    print(f"\n  {'N':>3s}  {'lyap_1':>12s}  {'l1/log(phi)':>12s}  {'l1/log(1.647)':>14s}  {'Full spectrum':60s}")

    for r in results[-10:]:
        n = r['n']
        l = r['lyapunov_exponents']
        l1 = l[0]
        phi_r = l1 / LOG_PHI
        lin_r = l1 / np.log(1.647)
        spec = '  '.join(f"{x:8.4f}" for x in l)
        print(f"  {n:3d}  {l1:12.6f}  {phi_r:12.6f}  {lin_r:14.6f}  [{spec}]")

    # The full spectrum at round 64
    sub_banner("Full Lyapunov Spectrum at Round 64")
    final = results[-1]['lyapunov_exponents']
    print(f"\n  Lyapunov exponents (descending):")
    for i, le in enumerate(final):
        phi_r = le / LOG_PHI
        is_phi = abs(phi_r - 1.0) < 0.02
        is_neg_phi = abs(phi_r + 1.0) < 0.02
        marker = ""
        if is_phi:
            marker = " <-- MATCHES log(phi)!!"
        elif is_neg_phi:
            marker = " <-- MATCHES -log(phi)!!"
        elif abs(le / np.log(1.647) - 1.0) < 0.02:
            marker = " <-- MATCHES log(1.647)!!"
        print(f"    lambda_{i+1} = {le:12.6f}  (ratio to log(phi): {phi_r:8.4f}){marker}")

    # Entropy rate = sum of positive Lyapunov exponents
    pos_sum = sum(le for le in final if le > 0)
    neg_sum = sum(le for le in final if le < 0)
    print(f"\n  Sum of positive exponents (KS entropy rate): {pos_sum:.6f}")
    print(f"  Sum of negative exponents:                   {neg_sum:.6f}")
    print(f"  Sum of all exponents:                        {sum(final):.6f}")
    print(f"  Positive sum / log(phi): {pos_sum / LOG_PHI:.6f}")

    return results


# ─── TEST 3: Lyapunov Exponents ────────────────────────────────────────────

def compute_lyapunov(state0, state1, W, n_rounds):
    """
    Compute Lyapunov exponent: lambda(N) = (1/N) * log(|delta(N)| / |delta(0)|)
    where delta is the state difference vector.
    """
    # Initial difference
    diff0 = np.array([(state1[i] - state0[i]) & 0xFFFFFFFF for i in range(8)], dtype=np.float64)
    for i in range(8):
        if diff0[i] >= 0x80000000:
            diff0[i] -= MOD
    norm0 = norm(diff0)
    if norm0 == 0:
        return 0.0, []

    s0 = list(state0)
    s1 = list(state1)
    exponents = []

    for r in range(min(n_rounds, 64)):
        s0 = sha_round(s0, K[r], W[r])
        s1 = sha_round(s1, K[r], W[r])

        diff_n = np.array([(s1[i] - s0[i]) & 0xFFFFFFFF for i in range(8)], dtype=np.float64)
        for i in range(8):
            if diff_n[i] >= 0x80000000:
                diff_n[i] -= MOD

        norm_n = norm(diff_n)
        n = r + 1
        if norm_n > 0:
            lam = np.log(norm_n / norm0) / n
        else:
            lam = float('-inf')
        exponents.append(lam)

    return exponents[-1] if exponents else 0.0, exponents


def test3_lyapunov():
    banner("TEST 3: Lyapunov Exponents — Golden vs Random vs Phi-spaced")
    print(f"  log(phi)  = {LOG_PHI:.6f}")
    print(f"  log(1.647) = {np.log(1.647):.6f}")

    W = zero_message()

    input_sets = [
        ("Golden (phi*2^32)", golden_state()),
        ("Fibonacci", fibonacci_state()),
        ("Phi-spaced", phi_spaced_state()),
        ("Random seed=42", random_state(42)),
        ("Random seed=123", random_state(123)),
        ("Random seed=777", random_state(777)),
        ("SHA-256 H_INIT", list(H_INIT)),
    ]

    all_exponents = {}

    for name, state in input_sets:
        sub_banner(f"Input: {name}")
        print(f"  State: {[f'0x{w:08X}' for w in state[:4]]}...")

        # Perturb by +1 in ALL words simultaneously
        perturbed = [(w + 1) & 0xFFFFFFFF for w in state]

        final_lam, exponents = compute_lyapunov(state, perturbed, W, 64)
        all_exponents[name] = exponents

        print(f"\n  {'Round':>5s}  {'lambda(N)':>12s}  {'lambda/log(phi)':>15s}  {'lambda/log(1.647)':>17s}")
        print(f"  {'-----':>5s}  {'--------':>12s}  {'-----------':>15s}  {'-------------':>17s}")

        for r_idx in [0, 1, 2, 3, 4, 7, 11, 15, 23, 31, 47, 63]:
            if r_idx < len(exponents):
                lam = exponents[r_idx]
                if lam > -100:
                    phi_r = lam / LOG_PHI if abs(LOG_PHI) > 0 else 0
                    lin_r = lam / np.log(1.647)
                    print(f"  {r_idx+1:5d}  {lam:12.6f}  {phi_r:15.6f}  {lin_r:17.6f}")
                else:
                    print(f"  {r_idx+1:5d}  {'  -inf':>12s}  {'diverged':>15s}  {'diverged':>17s}")

    # Compare final Lyapunov exponents
    sub_banner("Final Lyapunov Exponents (round 64)")
    print(f"\n  {'Input':30s}  {'lambda':>12s}  {'lambda/log(phi)':>15s}  {'Close to log(phi)?':>18s}")
    print(f"  {'-----':30s}  {'------':>12s}  {'-----------':>15s}  {'----------------':>18s}")

    for name, exponents in all_exponents.items():
        if exponents:
            lam = exponents[-1]
            if lam > -100:
                phi_r = lam / LOG_PHI
                close = abs(phi_r - 1.0) < 0.1
                print(f"  {name:30s}  {lam:12.6f}  {phi_r:15.6f}  {'YES!!' if close else 'no':>18s}")
            else:
                print(f"  {name:30s}  {'  -inf':>12s}  {'  N/A':>15s}  {'collapsed':>18s}")

    # Per-word Lyapunov: perturb ONE word at a time
    sub_banner("Per-Word Lyapunov (golden state, perturb word j only)")
    state = golden_state()
    print(f"\n  {'Word j':>6s}  {'lambda':>12s}  {'lambda/log(phi)':>15s}")
    print(f"  {'------':>6s}  {'------':>12s}  {'-----------':>15s}")

    for j in range(8):
        perturbed = list(state)
        perturbed[j] = (perturbed[j] + 1) & 0xFFFFFFFF
        final_lam, exponents = compute_lyapunov(state, perturbed, W, 64)
        if final_lam > -100:
            print(f"  {j:6d}  {final_lam:12.6f}  {final_lam/LOG_PHI:15.6f}")
        else:
            print(f"  {j:6d}  {'  -inf':>12s}  {'  N/A':>15s}")

    return all_exponents


# ─── TEST 4: SVD Analysis ──────────────────────────────────────────────────

def test4_svd():
    banner("TEST 4: SVD Analysis -- Log-Singular-Values & Condition Number Growth")
    print(f"  Using QR-based accumulation (same as Test 2) to prevent overflow.")
    print(f"  Singular value growth rates = Lyapunov exponents from Test 2.")
    print(f"  Here we focus on: condition number, SV ratios, spectral gaps.")

    state = golden_state()
    W = zero_message()
    trace = run_n_rounds_trace(state, W, 64)

    Q = np.eye(8)
    log_sv_sum = np.zeros(8)  # Accumulated log-singular-values
    results = []

    print(f"\n  {'N':>3s}  {'log_s1/N':>12s}  {'log_s8/N':>12s}  {'log(cond)/N':>12s}  {'gap(1-2)/N':>12s}  {'s1/s2 ratio':>12s}")
    print(f"  {'---':>3s}  {'--------':>12s}  {'--------':>12s}  {'-----------':>12s}  {'----------':>12s}  {'-----------':>12s}")

    for round_idx in range(64):
        J_local = compute_local_jacobian(trace[round_idx], round_idx, W)
        M = J_local @ Q
        Q, R = np.linalg.qr(M)

        # R diagonal gives the local stretching factors (singular values)
        r_diag = np.abs(np.diag(R))
        r_diag = np.maximum(r_diag, 1e-300)
        log_sv_sum += np.log(r_diag)

        n = round_idx + 1
        # Sort descending
        log_svs = np.sort(log_sv_sum / n)[::-1]

        log_s1 = log_svs[0]
        log_s8 = log_svs[-1]
        log_cond = log_s1 - log_s8
        gap_12 = log_svs[0] - log_svs[1] if len(log_svs) > 1 else 0
        sv_ratio_12 = np.exp(gap_12) if gap_12 < 100 else float('inf')

        results.append({
            'n': n,
            'log_svs': log_svs.copy(),
            'log_cond': log_cond,
            'gap_12': gap_12,
        })

        if n <= 20 or n % 4 == 0 or n == 64:
            print(f"  {n:3d}  {log_s1:12.6f}  {log_s8:12.6f}  {log_cond:12.6f}  {gap_12:12.6f}  {sv_ratio_12:12.6f}")

    # Growth rate of condition number
    sub_banner("Condition Number Growth Rate")
    print(f"  log(cond)/N = difference between largest and smallest Lyapunov exponent")
    print(f"  If this converges: the system has a well-defined spectral gap.")
    print(f"\n  Last 10 rounds:")
    for r in results[-10:]:
        n = r['n']
        lc = r['log_cond']
        print(f"    N={n:2d}: log(cond)/N = {lc:.6f}  "
              f"(ratio to log(phi): {lc/LOG_PHI:.4f})")

    # Full log-singular-value spectrum at key rounds
    sub_banner("Full Log-SV Spectrum (per round) at Key Rounds")
    for target_n in [1, 4, 8, 16, 32, 64]:
        r = results[target_n - 1]
        print(f"\n  N={target_n}:")
        for i, ls in enumerate(r['log_svs']):
            phi_r = ls / LOG_PHI
            marker = ""
            if abs(phi_r - 1.0) < 0.02:
                marker = " <-- log(phi)!!"
            elif abs(phi_r + 1.0) < 0.02:
                marker = " <-- -log(phi)!!"
            print(f"    log(sigma_{i+1})/N = {ls:12.6f}  (ratio to log(phi): {phi_r:8.4f}){marker}")

    # Ratios between consecutive log-singular-values at round 64
    sub_banner("Log-SV Ratios at Round 64 (spectral gaps)")
    final = results[-1]['log_svs']
    for i in range(len(final) - 1):
        gap = final[i] - final[i+1]
        gap_phi = gap / LOG_PHI
        close = abs(gap_phi - 1.0) < 0.1
        print(f"    gap_{i+1},{i+2} = {gap:12.6f}  (ratio to log(phi): {gap_phi:8.4f}){'  <-- log(phi)!!' if close else ''}")

    return results


# ─── TEST 5: Golden Trace ──────────────────────────────────────────────────

def test5_golden_trace():
    banner("TEST 5: Golden Trace — trace(J_local) at Each of 64 Rounds")

    state = golden_state()
    W = zero_message()
    trace = run_n_rounds_trace(state, W, 64)

    traces = []
    fib = [1, 1]
    while len(fib) < 20:
        fib.append(fib[-1] + fib[-2])

    print(f"\n  {'Round':>5s}  {'trace(J)':>14s}  {'trace mod 2^32':>14s}  {'trace/8':>12s}  {'Fib?':>6s}")
    print(f"  {'-----':>5s}  {'--------':>14s}  {'-----------':>14s}  {'-----':>12s}  {'----':>6s}")

    for round_idx in range(64):
        J_local = compute_local_jacobian(trace[round_idx], round_idx, W)
        tr = np.trace(J_local)
        tr_mod = int(tr) % MOD if abs(tr) < 1e18 else float('nan')
        tr_over_8 = tr / 8

        # Check if trace relates to Fibonacci
        is_fib = ""
        for f in fib:
            if abs(tr - f) < 0.5:
                is_fib = f"F={f}"
                break

        traces.append(tr)
        print(f"  {round_idx+1:5d}  {tr:14.4f}  {tr_mod:14d}  {tr_over_8:12.4f}  {is_fib:>6s}")

    # Analysis of trace sequence
    sub_banner("Trace Sequence Analysis")

    traces_arr = np.array(traces)
    print(f"\n  Mean trace:   {np.mean(traces_arr):.6f}")
    print(f"  Std trace:    {np.std(traces_arr):.6f}")
    print(f"  Min trace:    {np.min(traces_arr):.6f}")
    print(f"  Max trace:    {np.max(traces_arr):.6f}")
    print(f"  phi * 8 = {PHI * 8:.6f}  (if trace = 8*phi, each eigenvalue ~ phi)")
    print(f"  1.647 * 8 = {1.647 * 8:.6f}")

    # Check if trace is constant, oscillating, or drifting
    diffs = np.diff(traces_arr)
    print(f"\n  Trace differences (consecutive):")
    print(f"    Mean |diff|: {np.mean(np.abs(diffs)):.6f}")
    print(f"    Max |diff|:  {np.max(np.abs(diffs)):.6f}")
    if np.std(traces_arr) < 0.01 * abs(np.mean(traces_arr)):
        print(f"    VERDICT: Trace is nearly CONSTANT at {np.mean(traces_arr):.6f}")
    elif np.std(diffs) < 0.01 * np.mean(np.abs(diffs)):
        print(f"    VERDICT: Trace has REGULAR oscillation")
    else:
        print(f"    VERDICT: Trace is VARYING (not constant, not simple oscillation)")

    # Autocorrelation of trace sequence at Fibonacci offsets
    sub_banner("Trace Autocorrelation at Fibonacci Offsets")
    for offset in [1, 2, 3, 5, 8, 13, 21]:
        if offset < len(traces_arr):
            c = np.corrcoef(traces_arr[:-offset], traces_arr[offset:])[0, 1]
            if not np.isnan(c):
                print(f"    offset {offset:2d}: r = {c:8.5f}")

    # Check: does trace(J_local) at round i relate to K[i]?
    sub_banner("Trace vs Round Constant K[i] Correlation")
    k_arr = np.array(K[:64], dtype=np.float64)
    corr = np.corrcoef(traces_arr, k_arr)[0, 1]
    print(f"    Pearson r(trace, K) = {corr:.6f}")

    return traces


# ─── BONUS: Compare Multiple Input Types for Test 2 ────────────────────────

def bonus_multi_input_effective_matrix():
    banner("BONUS: Lyapunov Spectra -- Multiple Input Types")
    print(f"  Comparing dominant Lyapunov exponent across different inputs.")
    print(f"  If phi is structural: ALL inputs should give the same exponent.")

    W_zero = zero_message()
    W_golden = golden_message()

    configs = [
        ("Golden state + zero msg", golden_state(), W_zero),
        ("Golden state + golden msg", golden_state(), W_golden),
        ("Fibonacci state + zero msg", fibonacci_state(), W_zero),
        ("SHA H_INIT + zero msg", list(H_INIT), W_zero),
        ("Random(42) + zero msg", random_state(42), W_zero),
    ]

    all_lyap1 = []

    for name, state, W in configs:
        sub_banner(f"Config: {name}")

        trace = run_n_rounds_trace(state, W, 64)
        Q = np.eye(8)
        lyap_sum = np.zeros(8)

        for round_idx in range(64):
            J_local = compute_local_jacobian(trace[round_idx], round_idx, W)
            M = J_local @ Q
            Q, R = np.linalg.qr(M)
            r_diag = np.abs(np.diag(R))
            r_diag = np.maximum(r_diag, 1e-300)
            lyap_sum += np.log(r_diag)

        lyap_exponents = np.sort(lyap_sum / 64)[::-1]

        print(f"  Lyapunov spectrum (per round):")
        for i, le in enumerate(lyap_exponents):
            print(f"    lambda_{i+1} = {le:12.6f}  (ratio to log(phi): {le/LOG_PHI:8.4f})")

        l1 = lyap_exponents[0]
        all_lyap1.append((name, l1))
        print(f"\n  Dominant exponent: {l1:.6f}")
        print(f"  l1/log(phi) = {l1/LOG_PHI:.6f}")
        print(f"  l1/log(1.647) = {l1/np.log(1.647):.6f}")

    # Compare all dominant exponents
    sub_banner("Dominant Exponent Comparison Across Inputs")
    print(f"\n  {'Config':40s}  {'lambda_1':>12s}  {'l1/log(phi)':>12s}")
    print(f"  {'------':40s}  {'--------':>12s}  {'-----------':>12s}")
    for name, l1 in all_lyap1:
        print(f"  {name:40s}  {l1:12.6f}  {l1/LOG_PHI:12.6f}")

    lyap_vals = [l for _, l in all_lyap1]
    print(f"\n  Mean lambda_1:   {np.mean(lyap_vals):.6f}")
    print(f"  Std lambda_1:    {np.std(lyap_vals):.6f}")
    print(f"  Spread (max-min): {max(lyap_vals)-min(lyap_vals):.6f}")
    if np.std(lyap_vals) < 0.01 * abs(np.mean(lyap_vals)):
        print(f"  VERDICT: Dominant exponent is INPUT-INDEPENDENT (structural!)")
    else:
        print(f"  VERDICT: Dominant exponent VARIES with input")


# ─── MAIN ──────────────────────────────────────────────────────────────────

def main():
    print("=" * 78)
    print("  SHA-256 EIGENVALUE TRACKING: GOLDEN RATIO TRAJECTORY")
    print("  phi = 1.6180339887... | log(phi) = 0.4812118250...")
    print("  Linear dominant eigenvalue = 1.647 | log(1.647) = 0.4994...")
    print("=" * 78)

    t0 = time.time()

    # Run all 5 tests + bonus
    results1 = test1_jacobian_tracking()
    results2 = test2_effective_matrix()
    results3 = test3_lyapunov()
    results4 = test4_svd()
    results5 = test5_golden_trace()
    bonus_multi_input_effective_matrix()

    elapsed = time.time() - t0

    # ─── GRAND SUMMARY ─────────────────────────────────────────────────────
    banner("GRAND SUMMARY")

    print(f"\n  Total runtime: {elapsed:.1f}s")
    print(f"\n  Reference values:")
    print(f"    phi          = {PHI:.10f}")
    print(f"    log(phi)     = {LOG_PHI:.10f}")
    print(f"    1/phi        = {1/PHI:.10f}")
    print(f"    phi^2        = {PHI**2:.10f}")
    print(f"    log(1.647)   = {np.log(1.647):.10f}")

    # Test 1 summary
    print(f"\n  TEST 1 — Jacobian eigenvalue at N=1:")
    if results1:
        r1 = results1[0]
        print(f"    |lam_max| = {r1['mag_max']:.6f}")
        print(f"    Close to 1.647? {abs(r1['mag_max'] - 1.647) < 0.05}")
        print(f"    Close to phi?   {abs(r1['mag_max'] - PHI) < 0.05}")

    # Test 2 summary
    print(f"\n  TEST 2 -- Lyapunov spectrum (QR) after 64 rounds:")
    if results2:
        r_final = results2[-1]
        l1 = r_final['lyapunov_exponents'][0]
        print(f"    Max Lyapunov exponent = {l1:.6f}")
        print(f"    l1 / log(phi) = {l1/LOG_PHI:.6f}")
        print(f"    Matches log(phi)?   {abs(l1/LOG_PHI - 1.0) < 0.05}")
        print(f"    Matches log(1.647)? {abs(l1/np.log(1.647) - 1.0) < 0.05}")

    # Test 3 summary
    print(f"\n  TEST 3 — Lyapunov exponents:")
    if results3:
        golden_lam = results3.get("Golden (phi*2^32)", [])
        if golden_lam:
            gl = golden_lam[-1]
            print(f"    Golden input lambda = {gl:.6f}")
            print(f"    lambda / log(phi) = {gl/LOG_PHI:.6f}")

    # Test 4 summary
    print(f"\n  TEST 4 -- Condition number growth:")
    if results4:
        r_final = results4[-1]
        print(f"    log(cond)/64 = {r_final['log_cond']:.6f}")
        print(f"    Ratio to log(phi): {r_final['log_cond']/LOG_PHI:.4f}")

    # Test 5 summary
    print(f"\n  TEST 5 — Local Jacobian traces:")
    if results5:
        mean_tr = np.mean(results5)
        print(f"    Mean trace = {mean_tr:.6f}")
        print(f"    Mean trace / 8 = {mean_tr/8:.6f} (avg eigenvalue)")
        print(f"    Compare phi = {PHI:.6f}")

    print(f"\n{'=' * 78}")
    print(f"  END OF ANALYSIS")
    print(f"{'=' * 78}")


if __name__ == "__main__":
    main()
