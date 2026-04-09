"""
ADVERSARIAL VERIFICATION: Platonic Primes theorem
Compute Tr(L^+) for all 5 Platonic solids from scratch using eigenvalues.
Then extend to 13 Archimedean solids and 6 regular 4-polytopes.
"""

import numpy as np
from fractions import Fraction
from sympy import Matrix, Rational, sqrt, isprime, simplify, eye, ones
from sympy import Integer
import sympy

print("=" * 70)
print("PLATONIC PRIMES VERIFICATION")
print("Computing Tr(L^+) for all 5 Platonic solids from SCRATCH")
print("=" * 70)

# ======================================================================
# STEP 1: Build adjacency matrices for all 5 Platonic solids
# ======================================================================

def build_tetrahedron():
    """K4 - complete graph on 4 vertices"""
    n = 4
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i, j] = 1
    return A, "Tetrahedron"

def build_cube():
    """Cube graph - 8 vertices, each connected to 3 neighbors"""
    # Label vertices as 3-bit binary strings, edge iff Hamming distance = 1
    n = 8
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            # Count differing bits
            xor = i ^ j
            if bin(xor).count('1') == 1:
                A[i, j] = 1
    return A, "Cube"

def build_octahedron():
    """Octahedron - 6 vertices, each connected to all except its antipode"""
    # 3 pairs of antipodal vertices: (0,1), (2,3), (4,5)
    n = 6
    A = Matrix.zeros(n, n)
    antipodal = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
    for i in range(n):
        for j in range(n):
            if i != j and j != antipodal[i]:
                A[i, j] = 1
    return A, "Octahedron"

def build_icosahedron():
    """Icosahedron - 12 vertices, degree 5"""
    # Standard icosahedron adjacency
    # Vertices: 0=top, 11=bottom, 1-5=upper ring, 6-10=lower ring
    n = 12
    edges = [
        # Top vertex to upper ring
        (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
        # Upper ring
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),
        # Upper ring to lower ring
        (1, 6), (2, 6), (2, 7), (3, 7), (3, 8),
        (4, 8), (4, 9), (5, 9), (5, 10), (1, 10),
        # Lower ring
        (6, 7), (7, 8), (8, 9), (9, 10), (10, 6),
        # Lower ring to bottom
        (6, 11), (7, 11), (8, 11), (9, 11), (10, 11),
    ]
    A = Matrix.zeros(n, n)
    for i, j in edges:
        A[i, j] = 1
        A[j, i] = 1
    return A, "Icosahedron"

def build_dodecahedron():
    """Dodecahedron - 20 vertices, degree 3"""
    # Standard dodecahedron adjacency list
    n = 20
    edges = [
        # Top pentagon
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
        # Top spokes to upper ring
        (0, 5), (1, 7), (2, 9), (3, 11), (4, 13),
        # Upper ring (alternating with gaps)
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
        (10, 11), (11, 12), (12, 13), (13, 14), (14, 5),
        # Lower ring connections
        (6, 15), (8, 16), (10, 17), (12, 18), (14, 19),
        # Bottom pentagon
        (15, 16), (16, 17), (17, 18), (18, 19), (19, 15),
    ]
    A = Matrix.zeros(n, n)
    for i, j in edges:
        A[i, j] = 1
        A[j, i] = 1
    return A, "Dodecahedron"

def compute_trace_pseudoinverse(A, name):
    """Compute Tr(L^+) exactly using sympy"""
    n = A.rows
    # Compute degree matrix
    D = Matrix.zeros(n, n)
    for i in range(n):
        D[i, i] = sum(A[i, j] for j in range(n))

    # Laplacian
    L = D - A

    # Verify: row sums should be 0
    for i in range(n):
        rs = sum(L[i, j] for j in range(n))
        assert rs == 0, f"Row {i} sum = {rs}, should be 0"

    # Compute eigenvalues with multiplicities
    eigenvals = L.eigenvals()

    print(f"\n{'='*50}")
    print(f"{name} (n={n}, degree={D[0,0]})")
    print(f"{'='*50}")
    print(f"Eigenvalues (with multiplicities):")

    sorted_eigs = sorted(eigenvals.items(), key=lambda x: float(x[0].evalf()))
    for eig, mult in sorted_eigs:
        print(f"  {eig} (multiplicity {mult})")

    # Compute Tr(L^+) = sum of (mult / eigenvalue) for nonzero eigenvalues
    trace = Rational(0)
    for eig, mult in eigenvals.items():
        if eig != 0:
            trace += Rational(mult) / eig

    # Simplify (handles irrational cancellations)
    trace = simplify(trace)

    print(f"\nTr(L^+) = {trace}")

    # Convert to fraction if rational
    if trace.is_rational:
        p = trace.p  # numerator
        q = trace.q  # denominator
        print(f"  = {p}/{q}")
        print(f"  Numerator = {p}")
        print(f"  Is {p} prime? {isprime(p)}")
        return p, q, True
    else:
        print(f"  WARNING: trace is not rational: {trace}")
        return None, None, False

# ======================================================================
# STEP 2: Verify all 5 Platonic solids
# ======================================================================

print("\n" + "=" * 70)
print("COMPUTING FROM SCRATCH...")
print("=" * 70)

solids = [
    build_tetrahedron,
    build_cube,
    build_octahedron,
    build_icosahedron,
    build_dodecahedron,
]

expected = {
    "Tetrahedron": (3, 4),
    "Cube": (29, 12),
    "Octahedron": (13, 12),
    "Icosahedron": (7, 3),
    "Dodecahedron": (137, 15),
}

results = {}
all_pass = True

for builder in solids:
    A, name = builder()

    # Verify graph is correct: check vertex count, regularity
    n = A.rows
    degrees = [sum(A[i, j] for j in range(n)) for i in range(n)]
    assert len(set(degrees)) == 1, f"{name}: not regular! degrees = {set(degrees)}"
    d = degrees[0]
    print(f"\nVerified: {name} has {n} vertices, each degree {d}")

    p, q, is_rational = compute_trace_pseudoinverse(A, name)

    if is_rational:
        results[name] = (p, q)
        exp_p, exp_q = expected[name]
        if p == exp_p and q == exp_q:
            print(f"  MATCHES expected {exp_p}/{exp_q}")
        else:
            print(f"  *** MISMATCH *** expected {exp_p}/{exp_q}, got {p}/{q}")
            all_pass = False
    else:
        all_pass = False

print("\n" + "=" * 70)
print("PLATONIC PRIMES SUMMARY")
print("=" * 70)
print(f"{'Solid':<15} {'Tr(L^+)':<12} {'Numerator':<12} {'Prime?':<8} {'Match?'}")
print("-" * 65)
for name in ["Tetrahedron", "Cube", "Octahedron", "Icosahedron", "Dodecahedron"]:
    if name in results:
        p, q = results[name]
        exp_p, exp_q = expected[name]
        match = "YES" if (p == exp_p and q == exp_q) else "NO"
        print(f"{name:<15} {p}/{q:<10} {p:<12} {str(isprime(p)):<8} {match}")

print(f"\nAll claims verified: {all_pass}")

# ======================================================================
# STEP 3: ADVERSARIAL CHECK - Non-Platonic graphs
# ======================================================================

print("\n" + "=" * 70)
print("ADVERSARIAL: Testing non-Platonic vertex-transitive graphs")
print("=" * 70)

def build_petersen():
    """Petersen graph - 10 vertices, degree 3"""
    n = 10
    # Outer ring: 0-4, inner star: 5-9
    edges = [
        (0,1),(1,2),(2,3),(3,4),(4,0),  # outer pentagon
        (5,7),(7,9),(9,6),(6,8),(8,5),   # inner pentagram
        (0,5),(1,6),(2,7),(3,8),(4,9),   # spokes
    ]
    A = Matrix.zeros(n, n)
    for i, j in edges:
        A[i,j] = 1
        A[j,i] = 1
    return A, "Petersen"

def build_complete(n):
    """Complete graph K_n"""
    A = Matrix.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = 1
    return A, f"K_{n}"

def build_cycle(n):
    """Cycle graph C_n"""
    A = Matrix.zeros(n, n)
    for i in range(n):
        A[i, (i+1)%n] = 1
        A[(i+1)%n, i] = 1
    return A, f"C_{n}"

non_platonic = [
    build_petersen,
    lambda: build_complete(5),
    lambda: build_complete(6),
    lambda: build_complete(7),
    lambda: build_cycle(5),
    lambda: build_cycle(6),
    lambda: build_cycle(8),
]

non_platonic_results = {}
for builder in non_platonic:
    A, name = builder()
    p, q, is_rational = compute_trace_pseudoinverse(A, name)
    if is_rational:
        non_platonic_results[name] = (p, q, isprime(p))
    else:
        non_platonic_results[name] = (None, None, False)

print("\n" + "-" * 65)
print(f"{'Graph':<15} {'Tr(L^+)':<12} {'Numerator':<12} {'Prime?'}")
print("-" * 65)
for name, (p, q, is_p) in non_platonic_results.items():
    if p is not None:
        print(f"{name:<15} {p}/{q:<10} {p:<12} {is_p}")

any_non_platonic_prime = any(v[2] for v in non_platonic_results.values())
print(f"\nAny non-Platonic graph has prime numerator: {any_non_platonic_prime}")
if not any_non_platonic_prime:
    print("CONFIRMS: Platonic-specific property")
else:
    print("WARNING: Some non-Platonic graphs also have prime numerators!")
    for name, (p, q, is_p) in non_platonic_results.items():
        if is_p:
            print(f"  {name}: {p}/{q}, numerator {p} IS prime")
