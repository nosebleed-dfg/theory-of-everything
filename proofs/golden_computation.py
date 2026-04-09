"""
GOLDEN_COMPUTATION — Fibonacci matrix M^2=M+I exploration; eigenvalues, exponential squaring, phi embeddings
nos3bl33d
"""

import numpy as np
import sys

np.set_printoptions(linewidth=200, precision=15)

phi = (1 + np.sqrt(5)) / 2
psi = (1 - np.sqrt(5)) / 2  # = -1/phi

print("=" * 80)
print("PART A: Fibonacci Matrix M = [[1,1],[1,0]]")
print("=" * 80)

M = np.array([[1, 1], [1, 0]], dtype=float)
M2 = M @ M
I2 = np.eye(2)
MpI = M + I2

print("M =")
print(M)
print("M^2 =")
print(M2)
print("M + I =")
print(MpI)
print("M^2 - (M+I) =")
print(M2 - MpI)
print("M^2 == M + I:", np.allclose(M2, MpI))
print("det(M) =", np.linalg.det(M))

eigenvalues = np.linalg.eigvals(M)
print("Eigenvalues of M:", eigenvalues)
print("phi =", phi, ", -1/phi =", -1/phi)
print("Match:", np.allclose(sorted(eigenvalues), sorted([phi, psi])))

print()

# Integer matrix operations
def mat_mult_int(A, B):
    return [
        [A[0][0]*B[0][0] + A[0][1]*B[1][0], A[0][0]*B[0][1] + A[0][1]*B[1][1]],
        [A[1][0]*B[0][0] + A[1][1]*B[1][0], A[1][0]*B[0][1] + A[1][1]*B[1][1]]
    ]

def mat_pow_int(M, n):
    if n == 0:
        return [[1,0],[0,1]]
    if n == 1:
        return [row[:] for row in M]
    if n % 2 == 0:
        half = mat_pow_int(M, n // 2)
        return mat_mult_int(half, half)
    else:
        return mat_mult_int(M, mat_pow_int(M, n - 1))

def mat_mod(M, m):
    return [[M[0][0] % m, M[0][1] % m], [M[1][0] % m, M[1][1] % m]]

def mat_mult_mod(A, B, m):
    return mat_mod(mat_mult_int(A, B), m)

def mat_pow_mod(M, n, m):
    if n == 0:
        return [[1,0],[0,1]]
    if n == 1:
        return mat_mod(M, m)
    if n % 2 == 0:
        half = mat_pow_mod(M, n // 2, m)
        return mat_mult_mod(half, half, m)
    else:
        return mat_mult_mod(M, mat_pow_mod(M, n - 1, m), m)

Mi = [[1,1],[1,0]]

print("=" * 80)
print("PART B: M^291 (the ceiling power)")
print("=" * 80)
M291 = mat_pow_int(Mi, 291)
print("M^291 = [[F_292, F_291], [F_291, F_290]]")
print("F_290 =", M291[1][1])
print("F_291 =", M291[1][0])
print("F_292 =", M291[0][0])
print("Digits in F_291:", len(str(M291[1][0])))
print("F_291 mod 137 =", M291[1][0] % 137)
print("F_291 mod 120 =", M291[1][0] % 120)
print("F_291 mod 20 =", M291[1][0] % 20)
print("F_291 mod 12 =", M291[1][0] % 12)
print("F_291 mod 5 =", M291[1][0] % 5)
print("det(M^291) = (-1)^291 =", (-1)**291)

print()

print("=" * 80)
print("PART C: M^583 (the cosmological power = 2*291 + 1)")
print("=" * 80)
M583 = mat_pow_int(Mi, 583)
print("Digits in F_583:", len(str(M583[1][0])))
print("F_583 mod 137 =", M583[1][0] % 137)
print("F_583 mod 120 =", M583[1][0] % 120)
print("F_583 mod 20 =", M583[1][0] % 20)
print("F_583 mod 5 =", M583[1][0] % 5)
print("Note: 583 = 2*291 + 1")

print()

print("=" * 80)
print("PART D: M^291 mod 137")
print("=" * 80)
M291_137 = mat_pow_mod(Mi, 291, 137)
print("M^291 mod 137 =", M291_137)

def pisano_period(m):
    prev, curr = 0, 1
    for i in range(1, 6*m*m):
        prev, curr = curr, (prev + curr) % m
        if prev == 0 and curr == 1:
            return i
    return -1

pi_137 = pisano_period(137)
print("Pisano period pi(137) =", pi_137)
print("291 mod pi(137) =", 291 % pi_137)
print("M^291 mod 137 entries:", M291_137[0][0], M291_137[0][1], M291_137[1][0], M291_137[1][1])
print("Trace mod 137:", (M291_137[0][0] + M291_137[1][1]) % 137)

print()

print("=" * 80)
print("PART E: Dodecahedron adjacency matrix")
print("=" * 80)

dodec_edges = [
    (0,1),(0,4),(0,5),
    (1,2),(1,6),
    (2,3),(2,7),
    (3,4),(3,8),
    (4,9),
    (5,10),(5,14),
    (6,10),(6,11),
    (7,11),(7,12),
    (8,12),(8,13),
    (9,13),(9,14),
    (10,15),
    (11,16),
    (12,17),
    (13,18),
    (14,19),
    (15,16),(15,19),
    (16,17),
    (17,18),
    (18,19)
]

A = np.zeros((20, 20), dtype=float)
for i, j in dodec_edges:
    A[i][j] = 1
    A[j][i] = 1

row_sums = A.sum(axis=1)
print("Row sums (should all be 3):", np.unique(row_sums))
print("Number of edges:", int(A.sum() / 2))

eigvals = np.linalg.eigvalsh(A)
eigvals_sorted = np.sort(eigvals)[::-1]
print("\nEigenvalues of dodecahedron adjacency matrix:")
unique_eigs = []
for e in eigvals_sorted:
    if not any(abs(e - u) < 1e-10 for u in unique_eigs):
        unique_eigs.append(e)
for e in unique_eigs:
    mult = sum(1 for x in eigvals_sorted if abs(x - e) < 1e-10)
    print("  %.8f (multiplicity %d)" % (e, mult))

print("\nsqrt(5) = %.8f" % np.sqrt(5))
print("phi = %.8f" % phi)
print("Eigenvalues are: 3, sqrt(5), 1, -2, -sqrt(5)")

# Check A^k = c*A + d*I + e*J forms
I20 = np.eye(20)
J20 = np.ones((20, 20))

print("\nChecking A^k ~ c*A + d*I (least squares fit):")
for k in range(2, 8):
    Ak = np.linalg.matrix_power(A, k)
    X = np.column_stack([A.flatten(), I20.flatten()])
    y = Ak.flatten()
    coeffs, residual, _, _ = np.linalg.lstsq(X, y, rcond=None)
    c_fit, d_fit = coeffs
    error = np.max(np.abs(Ak - c_fit * A - d_fit * I20))
    print("  A^%d ~ %.4f*A + %.4f*I  (max error: %.4f)" % (k, c_fit, d_fit, error))

print("\nWith J (all-ones) term:")
for k in range(2, 8):
    Ak = np.linalg.matrix_power(A, k)
    X = np.column_stack([A.flatten(), I20.flatten(), J20.flatten()])
    y = Ak.flatten()
    coeffs, residual, _, _ = np.linalg.lstsq(X, y, rcond=None)
    c_fit, d_fit, e_fit = coeffs
    error = np.max(np.abs(Ak - c_fit * A - d_fit * I20 - e_fit * J20))
    print("  A^%d ~ %.4f*A + %.4f*I + %.4f*J  (max error: %.4f)" % (k, c_fit, d_fit, e_fit, error))

# Verify minimal polynomial
P1 = A - 3*I20
P2 = A - 1*I20
P3 = A + 2*I20
P4 = A @ A - 5*I20
result = P1 @ P2 @ P3 @ P4
print("\n|(A-3I)(A-I)(A+2I)(A^2-5I)| max = %.2e" % np.max(np.abs(result)))

# The golden projector
A2m5 = A @ A - 5*I20
eigvals_A2m5 = np.linalg.eigvalsh(A2m5)
print("\nA^2-5I eigenvalues:", sorted(set(np.round(eigvals_A2m5, 6))))
print("This kills the sqrt(5) eigenspaces (golden eigenspaces)")

# Does A satisfy A^2 = A + cI for any c?
print("\nDoes A^2 = A + c*I for any c? Checking per eigenvalue:")
for lam in [3, np.sqrt(5), 1, -2, -np.sqrt(5)]:
    c_needed = lam**2 - lam
    print("  lambda=%.4f: needs c=%.4f" % (lam, c_needed))
print("ANSWER: No single c works.")

# But check: the LAPLACIAN L = 3I - A
L = 3*I20 - A
print("\nLaplacian L = 3I - A")
L_eigs = np.linalg.eigvalsh(L)
L_eigs_sorted = np.sort(L_eigs)
L_unique = []
for e in L_eigs_sorted:
    if not any(abs(e - u) < 1e-10 for u in L_unique):
        L_unique.append(e)
print("Laplacian eigenvalues:", [round(e, 6) for e in L_unique])

# Golden connection
print("\nGolden connection to dodecahedron eigenvalues:")
print("sqrt(5) = phi + 1/phi = %.8f" % (phi + 1/phi))
print("3 - sqrt(5) = 2/phi^2 = %.8f vs %.8f" % (3 - np.sqrt(5), 2/phi**2))
print("3 + sqrt(5) = 2*phi^2 = %.8f vs %.8f" % (3 + np.sqrt(5), 2*phi**2))

print()

print("=" * 80)
print("PART F: Pisano periods")
print("=" * 80)

print("Pisano periods for key moduli:")
for m in [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 40, 60, 120, 137]:
    pp = pisano_period(m)
    divides = "divides 120" if m <= 120 and 120 % m == 0 else ""
    print("  pi(%3d) = %4d  %s" % (m, pp, divides))

pi_120 = pisano_period(120)
print("\nPisano period pi(120) =", pi_120)

print()

print("=" * 80)
print("PART G: M^120 mod 120 and special powers")
print("=" * 80)
M120_120 = mat_pow_mod(Mi, 120, 120)
print("M^120 mod 120 =", M120_120)
print("Is M^120 = I mod 120?", M120_120 == [[1,0],[0,1]])
print("Is M^120 = -I mod 120?", M120_120 == [[119,0],[0,119]])

M_pisano = mat_pow_mod(Mi, pi_120, 120)
print("M^%d mod 120 = %s" % (pi_120, M_pisano))

print("\nSpecial powers:")
for n in [10, 20, 30, 60, 120, 137, 291, 583]:
    Mn = mat_pow_mod(Mi, n, 120)
    print("  M^%3d mod 120 = %s" % (n, Mn))

print()

print("=" * 80)
print("PART H: Spin / Shifted Operator Interpretation")
print("=" * 80)
print("If A^2 = A + I, then:")
print("  (A - 1/2)^2 = A^2 - A + 1/4 = I + 1/4 = (5/4)I")
print("  So B = A - (1/2)I has B^2 = (5/4)I")
print("  B eigenvalues: +/- sqrt(5)/2 = +/- %.10f" % (np.sqrt(5)/2))
print()
print("Compare with quantum spin:")
print("  Spin-1/2: S^2 = (3/4)I, eigenvalues +/- 1/2")
print("  Golden op: B^2 = (5/4)I, eigenvalues +/- sqrt(5)/2")
print("  Spin-s: s(s+1) = j  =>  s = (-1+sqrt(1+4j))/2")
print("  For j=5/4: s = (-1+sqrt(6))/2 = %.6f" % ((-1+np.sqrt(6))/2))
print()
print("Projection decomposition:")
P_plus = (M + (1/phi)*I2) / (phi + 1/phi)
P_minus = (phi*I2 - M) / (phi + 1/phi)
print("  P_+ =", P_plus.tolist())
print("  P_- =", P_minus.tolist())
print("  P_+^2 = P_+:", np.allclose(P_plus @ P_plus, P_plus))
print("  P_-^2 = P_-:", np.allclose(P_minus @ P_minus, P_minus))
print("  P_+*P_- = 0:", np.allclose(P_plus @ P_minus, 0))
print("  P_+ + P_- = I:", np.allclose(P_plus + P_minus, I2))
print("  phi*P_+ + psi*P_- = M:", np.allclose(phi*P_plus + psi*P_minus, M))

print()

print("=" * 80)
print("PART I: Fibonacci Matrix in Modular Group")
print("=" * 80)
S = np.array([[0,-1],[1,0]], dtype=float)
T = np.array([[1,1],[0,1]], dtype=float)
print("PSL(2,Z) generators: S = [[0,-1],[1,0]], T = [[1,1],[0,1]]")
print("S^2 =", (S @ S).tolist(), "= -I")
print("(ST)^3 =", (np.linalg.matrix_power(S @ T, 3)).tolist())
print()
print("M = [[1,1],[1,0]]")
print("M^2 = [[2,1],[1,1]] IS in SL(2,Z)")

print()

print("=" * 80)
print("PART J: 291 and Pisano Connections")
print("=" * 80)
print("291 = VE/chi - d^2 = 300 - 9 = 20*30/2 - 9")
print("291 = 3 * 97")
print("Pisano(3) =", pisano_period(3))
print("Pisano(97) =", pisano_period(97))
print()

# Check: is 291 a Pisano period for any modulus?
print("Is 291 a Pisano period for any modulus up to 500?")
found_291 = False
for m in range(2, 501):
    if pisano_period(m) == 291:
        print("  pi(%d) = 291" % m)
        found_291 = True
if not found_291:
    print("  No -- 291 is not a Pisano period for any m in [2, 500]")

# F_291 mod dodecahedral numbers
print("\nF_291 mod dodecahedral numbers:")
for m in [2, 3, 5, 12, 20, 30, 120, 137]:
    f291_mod = mat_pow_mod(Mi, 291, m)
    pp = pisano_period(m)
    equiv = 291 % pp
    print("  F_291 mod %3d = %3d  (Pisano=%d, 291 mod Pisano = %d)" % (m, f291_mod[1][0], pp, equiv))

print()

print("=" * 80)
print("PART K: The Fusion Matrix (Fibonacci Anyons)")
print("=" * 80)
print("In the Fibonacci anyon category:")
print("  tau x tau = 1 + tau  (fusion rule)")
print("  This IS x^2 = 1 + x at the level of objects")
print()
N_tau = np.array([[0, 1], [1, 1]])
print("The N-matrix (fusion coefficients):")
print("  N_tau =", N_tau.tolist())
print()
print("N_tau^n gives fusion multiplicities for tau^(tensor n):")
for n in range(1, 10):
    Nn = np.linalg.matrix_power(N_tau, n)
    print("  N_tau^%d = %s  (tau^%d = %d*1 + %d*tau)" % (n, Nn[0].tolist(), n, int(Nn[0][0]), int(Nn[0][1])))
print()
N2 = N_tau @ N_tau
print("N_tau^2 =", N2.tolist())
print("N_tau + I =", (N_tau + np.eye(2, dtype=int)).tolist())
print("N_tau^2 == N_tau + I:", np.allclose(N2, N_tau + np.eye(2)))
print()
print("N_tau = [[0,1],[1,1]] and M = [[1,1],[1,0]] are conjugate")
print("Same char poly: x^2 - x - 1 = 0")
print("Same eigenvalues: phi and -1/phi")

print()

print("=" * 80)
print("PART L: All Integer 2x2 Solutions of X^2 = X + I")
print("=" * 80)
print("Trace = 1, det = -1")
count = 0
solutions = []
for a in range(-5, 6):
    d = 1 - a
    for b in range(-5, 6):
        if b == 0:
            continue
        c_num = a*d + 1  # ad - bc = -1 => c = (ad+1)/b
        if c_num % b == 0:
            c = c_num // b
            mat = np.array([[a, b], [c, d]], dtype=float)
            m2 = mat @ mat
            mpi = mat + np.eye(2)
            if np.allclose(m2, mpi):
                count += 1
                solutions.append((a, b, c, d))
                if count <= 15:
                    print("  [[%d,%d],[%d,%d]]" % (a, b, c, d))
print("Found %d solutions in range [-5,5]" % count)

print()

print("=" * 80)
print("PART M: Higher-dimensional solutions of X^2 = X + I")
print("=" * 80)
# Any diagonal matrix with entries in {phi, -1/phi}
# Any matrix conjugate to such a diagonal
# In nxn: spectrum must be subset of {phi, -1/phi}
# So the space of solutions in M_n(R) is parametrized by
# choosing k eigenvalues = phi and (n-k) = -1/phi, plus a basis

# For the dodecahedron: can we find a 20x20 matrix X^2 = X + I
# that commutes with the dodecahedral symmetry group?
# Yes: X = phi*P_golden + (-1/phi)*P_rest
# where P_golden projects onto the sqrt(5) eigenspaces of A

# Let's compute this!
eigvals_full, eigvecs = np.linalg.eigh(A)
print("Dodecahedron eigenspace dimensions:")
for e in unique_eigs:
    mask = np.abs(eigvals_full - e) < 1e-10
    dim = np.sum(mask)
    print("  lambda = %.4f: dim = %d" % (e, dim))

# Project onto sqrt(5) eigenspace
golden_mask = np.abs(np.abs(eigvals_full) - np.sqrt(5)) < 1e-10
P_golden = np.zeros((20, 20))
for i in range(20):
    if golden_mask[i]:
        v = eigvecs[:, i]
        P_golden += np.outer(v, v)

P_rest = I20 - P_golden

# Build X = phi*P_golden + (-1/phi)*P_rest
X_golden = phi * P_golden + psi * P_rest
X2 = X_golden @ X_golden
XpI = X_golden + I20
print("\nGolden operator X = phi*P_golden + psi*P_rest (20x20)")
print("X^2 = X + I:", np.allclose(X2, XpI))
print("Trace(X) = %.4f" % np.trace(X_golden))

golden_dim = np.sum(golden_mask)
rest_dim = 20 - golden_dim
print("Golden eigenspace dim:", golden_dim)
print("Rest dim:", rest_dim)
print("Expected trace: %d*phi + %d*psi = %.4f" % (golden_dim, rest_dim, golden_dim*phi + rest_dim*psi))
print("Actual trace: %.4f" % np.trace(X_golden))

# Does X commute with A?
print("Does X commute with A?", np.allclose(X_golden @ A, A @ X_golden))

print()
print("=" * 80)
print("PART N: The Temperley-Lieb Connection")
print("=" * 80)
print("Temperley-Lieb algebra TL_n(delta):")
print("  Generators: e_1, ..., e_{n-1}")
print("  Relations: e_i^2 = delta * e_i")
print("            e_i e_{i+1} e_i = e_i")
print("            e_i e_j = e_j e_i if |i-j| >= 2")
print()
print("At delta = phi (golden ratio):")
print("  e_i^2 = phi * e_i  (idempotent up to phi)")
print("  This is the FIBONACCI specialization")
print("  The Jones polynomial at q = e^(2*pi*i/5) uses this value")
print()
print("The Jones-Wenzl projector p_n in TL_n(phi):")
print("  Satisfies e_i * p_n = 0 = p_n * e_i for all i")
print("  dim(image of p_n) = Fibonacci number F_{n+1}")
print("  This is WHERE the Fibonacci numbers come from categorically!")
print()
print("Connection to our equation:")
print("  In the Grothendieck ring K_0(TL(phi)):")
print("  [tau] * [tau] = [1] + [tau]")
print("  This IS phi^2 = 1 + phi in the Grothendieck ring")
print("  The functor 'tensor with tau' satisfies F^2 = F + Id")
print("  at the level of the Grothendieck ring (decategorification)")

print()
print("DONE -- all computations complete.")
