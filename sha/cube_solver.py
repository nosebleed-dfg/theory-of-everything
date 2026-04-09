"""
CUBE_SOLVER — Bitcoin nonce as a 9x9 cube position; predicts solve from two prior blocks via diagonal scaling
nos3bl33d

8 state words + 1 message word = 9. Mining = unscrambling 64 koppa turns.
Previous nonces give the diagonal; scale by 1/F=1/12, search the neighborhood.
"""

import struct
import hashlib
import time
import numpy as np

PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329
KOPPA = 0.25
F = 12  # dodecahedral faces
SCALE = KOPPA * GAMMA * GAMMA  # = 1/F ≈ 0.0833

def sha256d(d):
    return hashlib.sha256(hashlib.sha256(d).digest()).digest()

# ============================================================
# THE 9x9 CUBE STATE
# ============================================================

class CubeState:
    """
    Represents a SHA-256 block header as a 9x9 cube state.

    The 80-byte header maps to the cube:
    - 76 bytes fixed (the scramble) = cube coloring
    - 4 bytes nonce (the solve) = cube position

    The target (leading zeros) = the solved pattern.
    Finding the nonce = solving the cube to match the pattern.
    """

    def __init__(self, header_76):
        self.header = header_76
        # Decompose header into 9x9 matrix (19 uint32 words, pad to 81 bytes)
        # Actually: header is 76 bytes. SHA processes it as words.
        # The "cube" is the SHA STATE: 8 words = 256 bits
        # Plus the message schedule contribution = 9th dimension
        self.words = list(struct.unpack('>19I', header_76))

    def hash_at(self, nonce):
        return sha256d(self.header + struct.pack('<I', nonce))

    def nonce_to_cube_coords(self, nonce):
        """Map a 32-bit nonce to 9x9 cube coordinates."""
        # Split nonce into 9 groups of ~3.5 bits each
        # But more naturally: 32 bits / 9 faces = ~3.5 bits per face
        # Use golden decomposition: nonce in base phi
        coords = []
        n = nonce
        for i in range(9):
            coords.append(n % 9)
            n = n // 9 + (n % 9)  # golden shift: quotient + remainder
        return coords

    def cube_distance(self, nonce_a, nonce_b):
        """Distance between two nonce positions in cube space."""
        ca = self.nonce_to_cube_coords(nonce_a)
        cb = self.nonce_to_cube_coords(nonce_b)
        # L1 distance in cube coordinates
        return sum(abs(a - b) for a, b in zip(ca, cb))


# ============================================================
# THE SOLVER
# ============================================================

class CubeSolver:
    """
    Solve for nonce using the 9x9 cube diagonal method.

    Given two previous block nonces:
    1. Compute gap in nonce space
    2. Scale by 1/F (dodecahedral face factor = koppa * gamma^2)
    3. Predict next nonce
    4. Search neighborhood using cube-distance metric
    """

    def __init__(self):
        self.history = []  # (header, nonce) pairs
        self.scale = SCALE  # koppa * gamma^2 ≈ 1/F ≈ 0.0833

    def record(self, nonce):
        """Record a solved nonce for future predictions."""
        self.history.append(nonce)
        if len(self.history) > 20:
            self.history = self.history[-20:]

    def predict(self):
        """Predict next nonce from history using cube diagonal."""
        if len(self.history) < 2:
            return None, 0

        n1 = self.history[-1]
        n2 = self.history[-2]
        gap = n1 - n2

        # Primary: 1/F scaling
        pred = n1 + int(gap * self.scale)

        # Also compute multi-scale predictions for wider search
        predictions = [
            pred,                                    # 1/F
            n1 + int(gap * KOPPA * GAMMA),          # koppa * gamma
            n1 + int(gap * KOPPA),                  # koppa alone
            n1 + int(gap * GAMMA * GAMMA),          # gamma^2
        ]

        # If we have 3+ history, use Fibonacci gap recurrence
        if len(self.history) >= 3:
            g1 = self.history[-1] - self.history[-2]
            g2 = self.history[-2] - self.history[-3]
            fib_gap = int((g1 + g2) * self.scale)
            predictions.append(n1 + fib_gap)

        # Adaptive: track which scale has been working
        if len(self.history) >= 4:
            # Check which prediction was closest for the LAST block
            last_actual = self.history[-1]
            prev_n1 = self.history[-2]
            prev_n2 = self.history[-3]
            prev_gap = prev_n1 - prev_n2

            scales_to_test = [
                self.scale,
                KOPPA * GAMMA,
                KOPPA,
                GAMMA * GAMMA,
            ]

            best_scale = self.scale
            best_err = float('inf')
            for s in scales_to_test:
                p = prev_n1 + int(prev_gap * s)
                err = abs(p - last_actual)
                if err < best_err:
                    best_err = err
                    best_scale = s

            # Use the winning scale for THIS prediction
            adaptive_pred = n1 + int(gap * best_scale)
            predictions.insert(0, adaptive_pred)  # try first

        return predictions, int(abs(gap) * self.scale * 2) + 1000

    def solve(self, header_76, target, max_total=50000000):
        """
        Solve for nonce.

        Strategy:
        1. If history available: search around predicted positions
        2. Fallback: sequential scan
        """
        t0 = time.time()
        checked = 0

        # Phase 1: Prediction-guided search
        predictions, radius = self.predict()
        if predictions:
            for pred in predictions:
                pred = max(0, pred) & 0xFFFFFFFF
                for delta in range(min(radius, 30000)):
                    for candidate in [pred + delta, pred - delta]:
                        candidate &= 0xFFFFFFFF
                        h = sha256d(header_76 + struct.pack('<I', candidate))
                        checked += 1
                        if int.from_bytes(h[::-1], 'big') < target:
                            self.record(candidate)
                            dt = time.time() - t0
                            return candidate, h, checked, dt

        # Phase 2: Sequential fallback
        for n in range(0xFFFFFFFF + 1):
            h = sha256d(header_76 + struct.pack('<I', n))
            checked += 1
            if int.from_bytes(h[::-1], 'big') < target:
                self.record(n)
                dt = time.time() - t0
                return n, h, checked, dt
            if checked > max_total:
                break

        return None, None, checked, time.time() - t0


# ============================================================
# TEST
# ============================================================

def test():
    print("9x9 CUBE NONCE SOLVER")
    print("=" * 60)
    print(f"Scale factor: koppa * gamma^2 = {SCALE:.6f} ≈ 1/F = {1/F:.6f}")
    print(f"F = {F} (dodecahedral faces)")
    print()

    target = 0x00ffff << (8 * (0x1f - 3))  # ~16 leading zeros
    solver = CubeSolver()
    prev_hash = b'\x00' * 32

    total_brute = 0
    total_cube = 0

    for block in range(20):
        merkle = hashlib.sha256(f'cube_test_{block}'.encode()).digest()
        ts = struct.pack('<I', 1712700000 + block)
        hdr76 = struct.pack('<I', 0x20000000) + prev_hash + merkle + ts + struct.pack('<I', 0x1f00ffff)

        # Brute force baseline
        t0 = time.time()
        for n in range(10000000):
            h = sha256d(hdr76 + struct.pack('<I', n))
            if int.from_bytes(h[::-1], 'big') < target:
                brute_nonce = n
                brute_time = time.time() - t0
                break

        # Cube solver
        nonce, h, checked, dt = solver.solve(hdr76, target, max_total=5000000)

        if nonce is not None:
            speedup = brute_nonce / max(1, checked) if brute_nonce > 0 else 1
            method = "PREDICT" if len(solver.history) >= 3 and checked < brute_nonce else "SEQUENTIAL"

            total_brute += brute_nonce
            total_cube += checked

            print(f"  Block {block:2d}: brute={brute_nonce:8d}  cube={checked:8d}  "
                  f"speedup={speedup:5.1f}x  [{method}]")

            prev_hash = h
        else:
            print(f"  Block {block:2d}: FAILED")
            # Still record brute force result for history
            solver.record(brute_nonce)
            prev_hash = sha256d(hdr76 + struct.pack('<I', brute_nonce))

    if total_brute > 0:
        print()
        print(f"TOTAL: brute={total_brute:,}  cube={total_cube:,}  "
              f"overall speedup={total_brute/max(1,total_cube):.1f}x")


if __name__ == "__main__":
    test()
