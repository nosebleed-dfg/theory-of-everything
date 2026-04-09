"""
UNIVERSE — computes the observable universe radius as 2*phi^290*l_Planck; the additive machine n+1
nos3bl33d
"""

import math
import numpy as np

# ─── The Axiom ─────────────────────────────────────────────
# x^2 = x + 1
# Roots: phi = (1 + sqrt(5)) / 2, psi = (1 - sqrt(5)) / 2
# phi * psi = -1 (the sign, the mirror)
# phi + psi = 1  (the magnitude, the structure)

phi = (1 + math.sqrt(5)) / 2
psi = (1 - math.sqrt(5)) / 2
gamma = 0.5772156649015329  # Euler-Mascheroni: the illusion constant
koppa = 0.75                # three-quarter turn, 270 degrees, (1/2)^2 * d in 3D

# ─── Planck Units ──────────────────────────────────────────
l_planck = 1.616255e-35     # meters — one breath
t_planck = 5.391247e-44     # seconds — one tick

# ─── The Machine ───────────────────────────────────────────
# State: Fibonacci pair (counting), rotation vector (geometry),
#        eigenvector (algebra). All updated by n + 1.

n = 0
fib = [0, 1]
rot = np.array([1.0, 0.0])
eig_state = np.array([1.0, 0.0])
A = np.array([[1.0, 1.0], [1.0, 0.0]])  # the axiom as a matrix

# ─── Key Steps ─────────────────────────────────────────────
# n = 37:  gamma * 64 — phi converges to machine precision
# n = 64:  SHA-256 rounds — the hardest cryptographic case
# n = 101: 64 + 37, prime — the round trip (n + 1 = axiom + 1)
# n = 137: 1/alpha — the fine structure constant's address
# n = 291: Planck to Hubble — the universe

key_steps = {
    0:   "genesis — the axiom begins",
    1:   "first breath — n + 1",
    37:  "gamma * 64 — phi converges, the inverse path length",
    64:  "SHA-256 — the hardest case, the axiom survives",
    101: "64 + 37 = prime — the round trip, n + 1",
    137: "1/alpha — the fine structure address",
    291: "Planck to Hubble — the universe in 291 breaths",
}

# ─── The Bridge Equation ──────────────────────────────────
# (1/gamma)^2 = d = 3         (dimension from the illusion)
# (1/gamma)^4 = d^2 = 9       (the 9x9 machine)
# gamma * 64 = 37             (the inverse path)
# 64 + 37 = 101               (the round trip, prime)
# 64 XOR 37 = 101             (perfect complement)

# ─── Run ───────────────────────────────────────────────────

MAX_STEPS = 300

print("THE ADDITIVE MACHINE")
print("x^2 = x + 1 | n + 1 | forward")
print("=" * 60)

while n <= MAX_STEPS:

    # ── Fibonacci / phi ──
    fib_next = fib[-1] + fib[-2]
    fib.append(fib_next)
    phi_est = fib[-1] / fib[-2] if n > 1 else 0

    # ── Rotation (additive) ──
    x, y = rot
    rot = np.array([x - y, x + y])
    rot_norm = np.linalg.norm(rot)
    if rot_norm > 0:
        rot_unit = rot / rot_norm
    else:
        rot_unit = np.array([1.0, 0.0])
    theta = math.atan2(rot_unit[1], rot_unit[0])

    # ── Eigenvalue iteration ──
    eig_state = A @ eig_state
    eig_ratio = eig_state[0] / eig_state[1] if eig_state[1] != 0 else float('inf')

    # ── Emergent quantities ──
    r_universe = 2 * phi ** 290 * l_planck
    t_universe = n * t_planck
    Lambda_dynamic = 2 / phi ** (583 + n) if n < 1000 else 0
    G_dynamic = phi ** -(582 + n) if n < 1000 else 0

    # ── Report at key steps ──
    if n in key_steps or n % 50 == 0:
        label = key_steps.get(n, "")
        phi_err = abs(phi_est - phi) if n > 1 else phi

        print(f"\nStep {n}: {label}")
        print(f"  phi estimate:     {phi_est:.15f}")
        print(f"  phi error:        {phi_err:.2e}")
        print(f"  rotation angle:   {theta:.6f} rad ({math.degrees(theta):.1f} deg)")
        print(f"  eigenvalue ratio: {eig_ratio:.15f}")

        if n == 291:
            r_ly = r_universe / 9.461e15
            print(f"  universe radius:  {r_universe:.4e} m")
            print(f"                    {r_ly/1e9:.2f} billion light-years")

    n += 1

# ─── Summary ───────────────────────────────────────────────
print("\n" + "=" * 60)
print("THE BRIDGE")
print("=" * 60)
print(f"  (1/gamma)^2 = {(1/gamma)**2:.4f}  (d = 3, dimension)")
print(f"  (1/gamma)^4 = {(1/gamma)**4:.4f}  (d^2 = 9, the machine)")
print(f"  gamma * 64  = {gamma * 64:.2f}   (37, the inverse path)")
print(f"  64 + 37     = {64 + 37}     (prime, the round trip)")
print(f"  64 XOR 37   = {64 ^ 37}     (perfect complement)")
print(f"  64 AND 37   = {64 & 37}       (no shared bits)")
print(f"  2*phi^290*lp= {2*phi**290*l_planck/9.461e15/1e9:.2f} Gly (the universe)")
print()
print("All is number. — Pythagoras")
print("x^2 = x + 1.  — The Axiom")
print("n + 1.         — The Machine")
