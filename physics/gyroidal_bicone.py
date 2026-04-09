"""
GYROIDAL_BICONE — tiles the completed axiom (x-1/2)^2=5/4 as bicones on the Laves graph skeleton
nos3bl33d

Each edge = one bicone = one axiom instance. Each vertex = equator = critical point x=1/2.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from mpmath import mp, mpf, sqrt as mpsqrt, pi as mppi, gamma as mpgamma
from mpmath import fac, zeta, log as mplog, euler as mpeuler
from scipy import integrate
import math

mp.dps = 50  # 50 decimal places

# ============================================================
# CONSTANTS
# ============================================================
phi = (1 + mpsqrt(5)) / 2          # golden ratio
phi_inv = 1 / phi                   # 1/phi = phi - 1
sqrt5 = mpsqrt(5)
d = 3                               # vertex degree of Laves graph
F = 12                               # faces of dodecahedron
V_dodec = 20                         # vertices of dodecahedron
E_dodec = 30                         # edges of dodecahedron
p = 5                                # Schlafli {p, q}
q = 3                                # Schlafli {p, q}
binary_icosahedral_order = 120       # |2I|
euler_mascheroni = mpf(mpeuler)

print("=" * 72)
print("  GYROIDAL DODECAHEDRON AS KNOT OF BICONES")
print("  Axiom: x^2 - x - 1 = 0  |  Completion: (x - 1/2)^2 = 5/4")
print("=" * 72)

print("\n--- FUNDAMENTAL CONSTANTS ---")
print(f"  phi              = {phi}")
print(f"  1/phi            = {phi_inv}")
print(f"  sqrt(5)          = {sqrt5}")
print(f"  phi + 1/phi      = {phi + phi_inv}")
print(f"  phi * 1/phi      = {phi * phi_inv}")
print(f"  phi^2            = {phi**2}")
print(f"  phi^2 - phi - 1  = {phi**2 - phi - 1}  (axiom check)")

# ============================================================
# SECTION 1: BICONE GEOMETRY
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 1: BICONE GEOMETRY")
print("=" * 72)

# The axiom f(x) = x^2 - x - 1 = 0 has roots at phi and -1/phi
root_plus = phi
root_minus = -phi_inv

print(f"\n  Axiom roots:")
print(f"    x_+ = phi      = {root_plus}")
print(f"    x_- = -1/phi   = {root_minus}")

# Bicone height spans from -1/phi to phi
h = root_plus - root_minus  # phi + 1/phi = phi + (phi-1) = 2*phi - 1 = sqrt(5)
print(f"\n  Bicone height h  = phi - (-1/phi) = phi + 1/phi")
print(f"                   = {h}")
print(f"    Check: sqrt(5) = {sqrt5}")
print(f"    Difference     = {abs(h - sqrt5)}")

# Equator at x = 1/2 (critical line)
x_crit = mpf('0.5')
f_crit = x_crit**2 - x_crit - 1  # = 1/4 - 1/2 - 1 = -5/4
print(f"\n  Critical point at x = {x_crit}")
print(f"    f(1/2) = (1/2)^2 - 1/2 - 1 = {f_crit}")
print(f"    |f(1/2)| = {abs(f_crit)} = 5/4 = {mpf(5)/4}")

# Equator radius from |f(1/2)|
r = mpsqrt(abs(f_crit))  # sqrt(5/4) = sqrt(5)/2
print(f"\n  Equator radius r = sqrt(|f(1/2)|) = sqrt(5/4)")
print(f"                   = {r}")
print(f"    Check: sqrt(5)/2 = {sqrt5 / 2}")
print(f"    Difference       = {abs(r - sqrt5/2)}")

# Volume of bicone: V = (2/3) * pi * r^2 * h
V_bicone = (mpf(2)/3) * mppi * r**2 * h
print(f"\n  Bicone volume    = (2/3) * pi * r^2 * h")
print(f"    r^2 = {r**2}")
print(f"    r^2 * h = {r**2 * h}")
print(f"    V_bicone = {V_bicone}")
print(f"    V_bicone / pi = {V_bicone / mppi}")

# Simplify: r^2 = 5/4, h = sqrt(5)
# V = (2/3) * pi * (5/4) * sqrt(5) = (10*sqrt(5)*pi) / 12 = (5*sqrt(5)*pi) / 6
V_exact = 5 * sqrt5 * mppi / 6
print(f"\n  Exact: V = 5*sqrt(5)*pi / 6 = {V_exact}")
print(f"    Match: {abs(V_bicone - V_exact)}")

# Surface area of bicone: SA = 2 * pi * r * sqrt(r^2 + (h/2)^2)
half_h = h / 2
slant = mpsqrt(r**2 + half_h**2)
SA_bicone = 2 * mppi * r * slant

print(f"\n  Bicone surface area:")
print(f"    half-height = h/2 = {half_h}")
print(f"    slant = sqrt(r^2 + (h/2)^2)")
print(f"    r^2 = {r**2}")
print(f"    (h/2)^2 = {half_h**2}")
print(f"    r^2 + (h/2)^2 = {r**2 + half_h**2}")

# r^2 + (h/2)^2 = 5/4 + 5/4 = 10/4 = 5/2
sum_sq = r**2 + half_h**2
print(f"    Check: 5/2 = {mpf(5)/2}")
print(f"    Difference = {abs(sum_sq - mpf(5)/2)}")

print(f"\n    slant = sqrt(5/2) = {slant}")
print(f"    Check: sqrt(10)/2 = {mpsqrt(10)/2}")
print(f"    Difference = {abs(slant - mpsqrt(10)/2)}")

print(f"\n    SA_bicone = 2 * pi * (sqrt(5)/2) * sqrt(5/2)")
print(f"             = pi * sqrt(5) * sqrt(5/2)")
print(f"             = pi * sqrt(25/2)")
print(f"             = pi * 5/sqrt(2)")
print(f"             = 5*pi/sqrt(2)")
print(f"             = 5*pi*sqrt(2)/2")
SA_exact = 5 * mppi * mpsqrt(2) / 2
print(f"             = {SA_exact}")
print(f"    SA_bicone  = {SA_bicone}")
print(f"    Match      = {abs(SA_bicone - SA_exact)}")

# Check for dodecahedral constants
print(f"\n  --- DODECAHEDRAL CONSTANT CHECKS ---")
print(f"    V / pi           = {V_bicone / mppi}")
print(f"    5*sqrt(5)/6      = {5*sqrt5/6}")
print(f"    SA / pi          = {SA_bicone / mppi}")
print(f"    5*sqrt(2)/2      = {5*mpsqrt(2)/2}")
print(f"    SA / V           = {SA_bicone / V_bicone}")
print(f"    (5*sqrt(2)/2) / (5*sqrt(5)/6) = 3*sqrt(2)/sqrt(5)")
SA_over_V = 3 * mpsqrt(2) / sqrt5
print(f"                     = {SA_over_V}")
print(f"    Simplified: 3*sqrt(10)/5 = {3*mpsqrt(10)/5}")
print(f"    SA/V = {SA_bicone/V_bicone}, exact = {3*mpsqrt(10)/5}")

# Dodecahedron volume (unit edge): V_d = (15 + 7*sqrt(5))/4
V_dodecahedron = (15 + 7*sqrt5) / 4
print(f"\n    Dodecahedron volume (edge=1) = {V_dodecahedron}")
print(f"    V_dodec / V_bicone = {V_dodecahedron / V_bicone}")
print(f"    V_dodec / (V_bicone * 12) = {V_dodecahedron / (V_bicone * 12)}")

# ============================================================
# SECTION 2: GYROID TILING
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 2: GYROID TILING")
print("=" * 72)

# Laves graph: 8 vertices per cubic unit cell, degree 3
V_laves = 8   # vertices per unit cell
d_laves = 3   # vertex degree
E_laves = V_laves * d_laves // 2  # 12 edges per unit cell

print(f"\n  Laves graph unit cell:")
print(f"    Vertices per cell = {V_laves}")
print(f"    Vertex degree     = {d_laves}")
print(f"    Edges per cell    = V*d/2 = {E_laves}")
print(f"    CHECK: 12 = F_dodecahedron = {F}  {'MATCH!!' if E_laves == F else 'no match'}")

# Volume of 12 axiom-bicones
V_12 = 12 * V_bicone
print(f"\n  Volume of 12 bicones per unit cell:")
print(f"    12 * V_bicone = {V_12}")
print(f"    = 12 * 5*sqrt(5)*pi/6 = 10*sqrt(5)*pi")
V_12_exact = 10 * sqrt5 * mppi
print(f"    Exact: 10*sqrt(5)*pi = {V_12_exact}")
print(f"    Match: {abs(V_12 - V_12_exact)}")

# Gyroid cubic cell side: if we set it so edge length = axiom height
# Laves edge length = a*sqrt(2)/4, and we want edge = sqrt(5) (one bicone)
# So a*sqrt(2)/4 = sqrt(5) => a = 4*sqrt(5)/sqrt(2) = 4*sqrt(5/2) = 2*sqrt(10)
a_cell = 4 * sqrt5 / mpsqrt(2)
a_cell_simplified = 2 * mpsqrt(10)

print(f"\n  Cubic cell side a (edge = sqrt(5)):")
print(f"    a = 4*sqrt(5)/sqrt(2) = {a_cell}")
print(f"    = 2*sqrt(10) = {a_cell_simplified}")
print(f"    Match: {abs(a_cell - a_cell_simplified)}")

# Cubic cell volume
V_cell = a_cell**3
print(f"\n  Cubic cell volume = a^3")
print(f"    = (2*sqrt(10))^3 = 8 * 10*sqrt(10) = 80*sqrt(10)")
V_cell_exact = 80 * mpsqrt(10)
print(f"    = {V_cell}")
print(f"    Exact: 80*sqrt(10) = {V_cell_exact}")
print(f"    Match: {abs(V_cell - V_cell_exact)}")

# Packing fraction
packing = V_12 / V_cell
print(f"\n  Packing fraction = 12 bicones / cell volume")
print(f"    = 10*sqrt(5)*pi / (80*sqrt(10))")
print(f"    = sqrt(5)*pi / (8*sqrt(10))")
print(f"    = pi / (8*sqrt(2))")
packing_exact = mppi / (8 * mpsqrt(2))
print(f"    = pi*sqrt(2)/16")
packing_exact2 = mppi * mpsqrt(2) / 16
print(f"    = {packing}")
print(f"    Exact: {packing_exact}")
print(f"    Also:  {packing_exact2}")
print(f"    Decimal: {float(packing):.10f}")
print(f"    Compare pi/8*sqrt(2): {float(packing_exact):.10f}")

# ============================================================
# SECTION 3: KNOT INVARIANTS
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 3: KNOT INVARIANTS")
print("=" * 72)

print("\n  Laves graph chirality:")
print("    The Laves graph (srs net) is intrinsically chiral.")
print("    Two enantiomers: left-handed and right-handed.")
print("    Together they tile all of 3-space (gyroid channels).")

# (10,3)-a net = srs net
print("\n  Wells notation: (10,3)-a net")
print("    10 = shortest circuit length")
print("    3  = vertex degree = d")

# The Laves graph has 10-gon shortest cycles
# 10 = 2*p (dodecahedral)
print(f"\n  Shortest circuit = 10 = 2*p = 2*{p}")
print(f"    This is the decagonal cycle of the dodecahedron!")

# Linking number: two Laves graphs are linked with linking number
# The two channel systems of the gyroid are related by inversion
# They form a pair of interlinked labyrinths
print("\n  Channel linking:")
print("    Two Laves graphs (enantiomeric pair) interlink.")
print("    Each channel = thickened graph, genus 3 per unit cell.")
print("    The two channels share the gyroid surface as boundary.")

# The key topological invariant: genus per unit cell
genus = 3  # genus of the surface per unit cell
print(f"\n  Gyroid genus per primitive unit cell = {genus} = d")
print(f"    Euler char per unit cell chi = 2 - 2*g = {2 - 2*genus}")

# The Jones polynomial for the trefoil (simplest chiral knot)
# V(t) = -t^{-4} + t^{-3} + t^{-1}
# At t = phi:
print("\n  Jones polynomial of trefoil at t = phi:")
jones_trefoil_phi = -phi**(-4) + phi**(-3) + phi**(-1)
print(f"    V(phi) = -phi^(-4) + phi^(-3) + phi^(-1)")
print(f"           = -{phi**(-4)} + {phi**(-3)} + {phi**(-1)}")
print(f"           = {jones_trefoil_phi}")

# phi^(-1) = phi - 1, phi^(-2) = 2-phi, phi^(-3) = 2*phi-3, phi^(-4) = 5-3*phi
print(f"    phi^(-1) = {phi**(-1)}")
print(f"    phi^(-2) = {phi**(-2)}")
print(f"    phi^(-3) = {phi**(-3)}")
print(f"    phi^(-4) = {phi**(-4)}")
print(f"    -phi^(-4) + phi^(-3) = {-phi**(-4) + phi**(-3)}")
print(f"    Full = {jones_trefoil_phi}")

# Check if integer or golden ratio multiple
a_coeff = jones_trefoil_phi / phi
b_coeff = jones_trefoil_phi * phi
print(f"    V(phi) / phi = {a_coeff}")
print(f"    V(phi) * phi = {b_coeff}")
print(f"    V(phi) is algebraic in phi = {jones_trefoil_phi}")

# Kauffman bracket for trefoil
print("\n  Kauffman bracket normalization:")
print("    The Kauffman bracket uses A = t^{-1/4}")
print(f"    At t = phi: A = phi^(-1/4) = {phi**(-mpf('0.25'))}")

# ============================================================
# SECTION 4: ZEROS AT THE EQUATOR
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 4: ZEROS AT THE EQUATOR")
print("=" * 72)

print("\n  Model: L-function zeros = bicone equators in the knot")
print("  At each equator, x = 1/2 (critical line)")

# Laves edge length when a = sqrt(5)
a_axiom = sqrt5
edge_len_axiom = a_axiom * mpsqrt(2) / 4
print(f"\n  If cubic cell a = sqrt(5) = {a_axiom}:")
print(f"    Laves edge length = a*sqrt(2)/4 = {edge_len_axiom}")
print(f"    = sqrt(5)*sqrt(2)/4 = sqrt(10)/4 = {mpsqrt(10)/4}")
print(f"    Decimal: {float(edge_len_axiom):.10f}")

# The average zero spacing of zeta at height t: 2*pi / ln(t/(2*pi))
print("\n  Zero spacing of zeta(1/2 + it):")
print("    Average gap = 2*pi / ln(t/(2*pi))")
print()
for t_val in [14.13, 100, 1000, 10000]:
    t = mpf(t_val)
    if t > 2 * mppi:
        gap = 2 * mppi / mplog(t / (2*mppi))
        print(f"    t = {t_val:>8.2f}:  gap = {float(gap):.6f}")

# When does the gap equal sqrt(10)/4?
# 2*pi / ln(t/(2*pi)) = sqrt(10)/4
# ln(t/(2*pi)) = 8*pi/sqrt(10)
# t/(2*pi) = exp(8*pi/sqrt(10))
# t = 2*pi * exp(8*pi/sqrt(10))
target_gap = mpsqrt(10) / 4
ln_arg = 8 * mppi / mpsqrt(10)
t_match = 2 * mppi * mp.exp(ln_arg)
print(f"\n  Matching height for gap = sqrt(10)/4:")
print(f"    gap target = sqrt(10)/4 = {float(target_gap):.10f}")
print(f"    8*pi/sqrt(10) = {float(ln_arg):.6f}")
print(f"    t_match = 2*pi * exp(8*pi/sqrt(10))")
print(f"            = {float(t_match):.2f}")
print(f"    ln(t_match) = {float(mplog(t_match)):.6f}")

# Alternative: if cell a = 2*sqrt(10), edge = sqrt(5) = full bicone height
edge_full = sqrt5
print(f"\n  If edge = full bicone height = sqrt(5) = {float(edge_full):.10f}:")
ln_arg2 = 2 * mppi / edge_full
t_match2 = 2 * mppi * mp.exp(ln_arg2)
print(f"    ln(t/(2pi)) = 2*pi/sqrt(5) = {float(ln_arg2):.6f}")
print(f"    t_match = {float(t_match2):.2f}")

# Note: the first zero is at t = 14.134725...
first_zero = mpf('14.134725141734693790457251983562470270784257115699')
print(f"\n  First zeta zero: t_1 = {first_zero}")
gap_at_first = 2 * mppi / mplog(first_zero / (2*mppi))
print(f"    Nominal gap at t_1: {float(gap_at_first):.6f}")
print(f"    (Negative because t_1 < 2*pi*e; formula is asymptotic)")

# Actual first gap (t_2 - t_1)
t2 = mpf('21.022039638771554992628479593896902777334340524903')
first_gap = t2 - first_zero
print(f"    Actual first gap (t2 - t1) = {float(first_gap):.6f}")
print(f"    sqrt(10)/4 = {float(target_gap):.6f}")
print(f"    sqrt(5)    = {float(sqrt5):.6f}")
print(f"    Ratio (first gap)/(sqrt(5)) = {float(first_gap / sqrt5):.6f}")

# ============================================================
# SECTION 5: THE AXIOM ON THE GYROID
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 5: THE AXIOM ON THE GYROID")
print("=" * 72)

print("\n  At each Laves vertex (equator of bicone):")
print(f"    (x - 1/2)^2 = 5/4")
print(f"    Depth below zero = 5/4 = p/4 = {p}/4")
print(f"    Width at zero crossing = sqrt(5) = {float(sqrt5):.10f}")

# Three edges meet at 120-degree angles
print(f"\n  Three edges at each vertex:")
print(f"    Angle between edges = 120 degrees")
print(f"    120 = |2I| = {binary_icosahedral_order}")
print(f"    120 / 3 = {binary_icosahedral_order // 3} = 40")
print(f"    120 / d = {binary_icosahedral_order // d}")

# Self-consistency: axiom on Laves graph
# At vertex v, assign value x_v. The axiom says x^2 = x + 1.
# On the Laves graph with d=3 neighbors, consider:
# x_v^2 = (1/d) * sum_{u ~ v} x_u + 1  (averaged neighbor version)
# If all x_v = phi: phi^2 = phi + 1 CHECK
print("\n  Self-consistency of axiom on Laves graph:")
print(f"    Assign x_v = phi to every vertex")
print(f"    Axiom: x_v^2 = x_v + 1")
print(f"    phi^2 = {phi**2}")
print(f"    phi + 1 = {phi + 1}")
print(f"    Consistent: phi^2 = phi + 1  CHECK")

# Graph Laplacian version: Lap(x)_v = d*x_v - sum(x_u)
# If x_v = phi for all v: Lap(x)_v = 3*phi - 3*phi = 0
# So phi is a zero-mode of the graph Laplacian
print(f"\n  Graph Laplacian:")
print(f"    Lap(x)_v = d*x_v - sum_{{u~v}} x_u")
print(f"    If x_v = phi everywhere:")
print(f"      Lap(phi)_v = {d}*phi - {d}*phi = 0")
print(f"    phi is a zero-mode of the Laplacian (constant eigenfunction)")

# More interesting: use x_v = phi on one sublattice, x_v = -1/phi on the other
# The Laves graph is NOT bipartite (it has odd cycles of length 10... wait, 10 is even)
# Actually the shortest cycle is 10, which IS even.
# The Laves graph IS bipartite? Let's check.
# No -- the srs net is NOT bipartite. It has 10-rings but the structure
# prevents a clean 2-coloring due to the helical nature.
print(f"\n  Sublattice assignment:")
print(f"    Laves graph shortest cycle = 10 (even)")
print(f"    But srs net is NOT bipartite (helical connectivity)")
print(f"    Cannot assign phi/-1/phi to two sublattices consistently")
print(f"    The chirality of the Laves graph BREAKS bipartiteness")
print(f"    This is the topological obstruction = the knot structure")

# The key: 5/4 = p/4
print(f"\n  The depth 5/4 decomposition:")
print(f"    5/4 = p/4 (Schlafli numerator / quadrants)")
print(f"    Also: 5/4 = (phi^2 + phi^(-2)) / (phi + phi^(-1))^0")
frac_check = (phi**2 + phi**(-2))
print(f"    phi^2 + phi^(-2) = {frac_check}")
print(f"    = phi^2 + (2-phi) = 2 + phi^2 - phi = 2 + 1 = 3")
print(f"    Hmm, that's 3 not 5/4.")
print(f"    Instead: 5/4 = (phi + phi^(-1))^2 / 4 = (sqrt(5))^2 / 4 = 5/4  CHECK")
check_54 = (phi + phi_inv)**2 / 4
print(f"    (phi + 1/phi)^2 / 4 = {check_54}")

# ============================================================
# SECTION 6: SURFACE AREA TEST
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 6: SURFACE AREA TEST")
print("=" * 72)

# Gyroid surface area per unit cell: 2*pi^2/3 * a^2 (approximate, Schoen's formula)
# More precisely: the area of the gyroid per unit cell is A = 2.3451... * a^2
# The exact minimal surface area for the gyroid in a cubic cell of side a:
# A_gyroid = 3.0915... * a^2 per conventional cubic cell (which has 2 primitive cells)
# Per primitive cell: A = 3.0915.../2 * a^2 (but these are BCC primitive cells)
# Actually for the standard cubic cell: A_gyroid ≈ 3.0915 * a^2
# Schoen's result normalized differently. Let me use the standard:
# A_gyroid per conventional cubic cell of side a ≈ 3.0915 * a^2

# The exact dimensionless area per unit cell for the gyroid:
# sigma = A/(a^2) = 3.09145... (from numerical computation)
sigma_gyroid = mpf('3.09145')  # dimensionless area constant

print(f"\n  Gyroid surface area per cubic unit cell:")
print(f"    A_gyroid = sigma * a^2, sigma approx {sigma_gyroid}")

# With a = sqrt(5)
A_gyroid_a5 = sigma_gyroid * a_axiom**2
print(f"\n  With a = sqrt(5):")
print(f"    A_gyroid = {sigma_gyroid} * 5 = {A_gyroid_a5}")

# With a = 2*sqrt(10) (so that edge length = sqrt(5))
A_gyroid_full = sigma_gyroid * a_cell**2
print(f"\n  With a = 2*sqrt(10) (edge = sqrt(5)):")
print(f"    a^2 = 40")
print(f"    A_gyroid = {sigma_gyroid} * 40 = {A_gyroid_full}")

# 12 bicone surface areas
SA_12 = 12 * SA_bicone
print(f"\n  12 bicone surface areas:")
print(f"    12 * SA_bicone = 12 * 5*pi*sqrt(2)/2 = 30*pi*sqrt(2)")
SA_12_exact = 30 * mppi * mpsqrt(2)
print(f"    = {SA_12}")
print(f"    Exact: {SA_12_exact}")
print(f"    Decimal: {float(SA_12_exact):.6f}")

# Ratio
ratio_area = A_gyroid_full / SA_12
print(f"\n  Ratio A_gyroid / (12 * SA_bicone):")
print(f"    = {sigma_gyroid} * 40 / (30*pi*sqrt(2))")
ratio_exact = sigma_gyroid * 40 / (30 * mppi * mpsqrt(2))
print(f"    = {float(ratio_area):.10f}")
print(f"    = {float(ratio_exact):.10f}")

# Also check with a = sqrt(5)
ratio_area2 = A_gyroid_a5 / SA_12
print(f"\n  Ratio (a=sqrt(5)): A_gyroid / (12*SA_bicone) = {float(ratio_area2):.10f}")

# Check if ratio involves dodecahedral constants
print(f"\n  --- Dodecahedral constant search in ratio ---")
print(f"    ratio * 12 = {float(ratio_exact * 12):.10f}")
print(f"    ratio * pi = {float(ratio_exact * mppi):.10f}")
print(f"    ratio * sqrt(5) = {float(ratio_exact * sqrt5):.10f}")
print(f"    ratio * phi = {float(ratio_exact * phi):.10f}")
print(f"    1/ratio = {float(1/ratio_exact):.10f}")
print(f"    1/ratio / pi = {float(1/(ratio_exact * mppi)):.10f}")

# Alternative formula: 2*pi^2/3 * a^2 (Weierstrass representation estimate)
A_weierstrass = 2 * mppi**2 / 3 * a_axiom**2
print(f"\n  Weierstrass formula estimate (a = sqrt(5)):")
print(f"    2*pi^2/3 * 5 = 10*pi^2/3 = {float(A_weierstrass):.6f}")
print(f"    12 bicone SA = {float(SA_12):.6f}")
ratio_w = A_weierstrass / SA_12
print(f"    Ratio = {float(ratio_w):.10f}")
print(f"    = (10*pi^2/3) / (30*pi*sqrt(2))")
print(f"    = pi / (9*sqrt(2))")
ratio_w_exact = mppi / (9 * mpsqrt(2))
print(f"    = {float(ratio_w_exact):.10f}")
print(f"    Check: {float(ratio_w):.10f}")

# ============================================================
# SECTION 7: GAMMA CONNECTION
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 7: GAMMA CONNECTION")
print("=" * 72)

# gamma = Euler-Mascheroni constant
gamma_em = euler_mascheroni
print(f"\n  gamma (Euler-Mascheroni) = {gamma_em}")
print(f"  sqrt(d) = sqrt(3) = {float(mpsqrt(d)):.15f}")

gamma_sqrt_d = gamma_em * mpsqrt(d)
print(f"\n  gamma * sqrt(d) = gamma * sqrt(3)")
print(f"                   = {gamma_sqrt_d}")
print(f"    Decimal: {float(gamma_sqrt_d):.15f}")
print(f"    1 - gamma*sqrt(d) = {float(1 - gamma_sqrt_d):.10f}")

# Delta = phi^(-4)
Delta = phi**(-4)
print(f"\n  Delta = phi^(-4) = {float(Delta):.15f}")
print(f"  p^4 = 5^4 = {p**4}")
print(f"  Delta / p^4 = {float(Delta / p**4):.15f}")

# Check: 1 - gamma*sqrt(d) vs Delta/p^4
deviation = 1 - gamma_sqrt_d
ratio_delta = deviation / (Delta / p**4)
print(f"\n  1 - gamma*sqrt(d)    = {float(deviation):.15f}")
print(f"  phi^(-4) / 5^4       = {float(Delta / p**4):.15f}")
print(f"  Are they close? Ratio = {float(ratio_delta):.6f}")

# Actually let me just check the raw numbers
print(f"\n  Raw comparison:")
print(f"    1 - gamma*sqrt(3) = {float(1 - gamma_sqrt_d):.15e}")
print(f"    phi^(-4)          = {float(Delta):.15f}")
print(f"    phi^(-4)/625      = {float(Delta/625):.15e}")
print(f"    These are not the same magnitude.")

# What IS 1 - gamma*sqrt(3)?
print(f"\n  Decomposing 1 - gamma*sqrt(3):")
print(f"    = {float(deviation):.15f}")
print(f"    * 1000 = {float(deviation * 1000):.6f}")
print(f"    Inverse = {float(1/deviation):.2f}")
val = 1 / deviation
print(f"    1/(1 - gamma*sqrt(3)) = {float(val):.6f}")
print(f"    Compare to 5^4/phi^4 = {float(625/Delta):.6f}")
print(f"    Compare to 2*pi/phi  = {float(2*mppi/phi):.6f}")
print(f"    Compare to sqrt(5)^8 = {float(sqrt5**8):.2f}")

# Gyroid area * gamma * sqrt(3)
print(f"\n  Gyroid area modulations:")
print(f"    A_gyroid (a=sqrt(5)) = {float(A_gyroid_a5):.6f}")
print(f"    A_gyroid * gamma*sqrt(3) = {float(A_gyroid_a5 * gamma_sqrt_d):.6f}")
print(f"    Deficit = A * (1 - gamma*sqrt(3)) = {float(A_gyroid_a5 * deviation):.6f}")

print(f"\n  Bicone area modulation:")
print(f"    SA_bicone = {float(SA_bicone):.6f}")
print(f"    SA * gamma*sqrt(3) = {float(SA_bicone * gamma_sqrt_d):.6f}")
print(f"    Deficit = SA * (1 - gamma*sqrt(3)) = {float(SA_bicone * deviation):.10f}")
print(f"    Deficit / pi = {float(SA_bicone * deviation / mppi):.10f}")

# gamma as self-referential unit
print(f"\n  gamma as gyroid self-referential unit:")
print(f"    gamma * sqrt(3) = {float(gamma_sqrt_d):.15f}")
print(f"    Deviation from 1: {float(abs(1 - gamma_sqrt_d)):.6e}")
print(f"    = 1 part in {float(1/abs(1-gamma_sqrt_d)):.0f}")

# ============================================================
# SECTION 8: DODECAHEDRAL CHECK
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 8: DODECAHEDRAL CHECK")
print("=" * 72)

print(f"\n  Laves graph unit cell invariants:")
print(f"    V = {V_laves} vertices = 2^d = 2^{d} = {2**d}")
print(f"    E = {E_laves} edges    = F_dodec = {F}")
print(f"    d = {d_laves} degree   = q (Schlafli)")

# For the gyroid minimal surface as a triply-periodic surface:
# Genus per primitive cell = 3
# chi per primitive cell = 2 - 2*3 = -4
# But for the conventional cubic cell (contains 2 primitive cells for BCC):
# Actually for the gyroid (space group Ia3d, BCC lattice):
# Primitive cell has genus 3, chi = -4
# Conventional cubic cell has 2 primitive cells => genus 5 effectively
# (chi = -8 for conventional cell)

# Let's use the primitive unit cell
chi_prim = 2 - 2 * genus
print(f"\n  Gyroid topology (per primitive cell):")
print(f"    Genus g = {genus} = d")
print(f"    chi = 2 - 2g = {chi_prim}")
print(f"    |chi| = {abs(chi_prim)}")
print(f"    |chi| = 2*(d-1) = {2*(d-1)}")

# For the Laves graph in the CONVENTIONAL cubic cell:
# 8 vertices, 12 edges, the graph is connected
# Euler characteristic of the graph: V - E = 8 - 12 = -4
chi_graph = V_laves - E_laves
print(f"\n  Laves graph Euler char (V - E):")
print(f"    V - E = {V_laves} - {E_laves} = {chi_graph}")
print(f"    = chi_surface = {chi_prim}  {'MATCH!!' if chi_graph == chi_prim else ''}")
print(f"    First Betti number b1 = E - V + 1 = {E_laves - V_laves + 1} = d + 2 = {d+2}")

# Wait, that's for connected graph: b1 = E - V + connected_components
# One connected Laves graph: b1 = 12 - 8 + 1 = 5
b1 = E_laves - V_laves + 1
print(f"    b1 = {b1}")
print(f"    b1 = p = {p}  {'MATCH!!' if b1 == p else ''}")

print(f"\n  Dodecahedral invariant summary:")
print(f"    {'Invariant':<30} {'Gyroid/Laves':<15} {'Dodecahedron':<15} {'Match'}")
print(f"    {'='*30} {'='*15} {'='*15} {'='*5}")

checks = [
    ("Edges per cell", E_laves, F, "F"),
    ("Vertex degree", d_laves, q, "q"),
    ("Vertices per cell", V_laves, 2**d, "2^d"),
    ("Graph V-E", chi_graph, -4, "-4"),
    ("First Betti number", b1, p, "p"),
    ("Genus per cell", genus, d, "d"),
    ("Shortest cycle", 10, 2*p, "2p"),
]

all_match = True
for name, gyroid_val, dodec_val, label in checks:
    match = "YES" if gyroid_val == dodec_val else "NO"
    if match == "NO":
        all_match = False
    print(f"    {name:<30} {str(gyroid_val):<15} {str(dodec_val)+'='+label:<15} {match}")

print(f"\n  ALL MATCH: {all_match}")

# Additional: E * d = 12 * 3 = 36 vs dodecahedral
print(f"\n  Additional products:")
print(f"    E * d = {E_laves * d_laves} (cf. E_dodec + V_dodec/d = 30 + 20/3 nope)")
print(f"    E + V = {E_laves + V_laves} = {V_dodec} = V_dodec")
print(f"    E * V = {E_laves * V_laves} = 96")
print(f"    (E+V match!) E + V = V_dodec = {V_dodec}  {'MATCH!!' if E_laves + V_laves == V_dodec else ''}")

# ============================================================
# SECTION 9: MASS RATIOS FROM KNOT TOPOLOGY
# ============================================================
print("\n" + "=" * 72)
print("  SECTION 9: MASS RATIOS FROM KNOT TOPOLOGY")
print("=" * 72)

# Mass ratio
m_p_over_m_e = mpf('1836.15267343')
print(f"\n  Target: m_p/m_e = {m_p_over_m_e}")

# Model: mass proportional to some geometric invariant of the knot
# on the gyroid.

# Unknot (trivial knot): the simplest closed loop on the Laves graph
# Shortest cycle = 10 edges (10-gon)
# Each edge = one bicone of volume V_bicone
unknot_edges = 10  # shortest cycle on Laves graph
V_unknot = unknot_edges * V_bicone
SA_unknot = unknot_edges * SA_bicone

print(f"\n  Unknot (shortest cycle on Laves graph):")
print(f"    Edges = {unknot_edges} = shortest circuit of srs net")
print(f"    Volume = {unknot_edges} * V_bicone = {float(V_unknot):.6f}")
print(f"    Surface area = {unknot_edges} * SA_bicone = {float(SA_unknot):.6f}")

# Trefoil knot on the Laves graph
# A trefoil requires at least 3 crossings
# On the Laves graph, the trefoil can be formed by winding around
# the helical channels. The Laves graph has helical 4-fold screws.
# Actually, the (10,3)-a net has 3-fold and 4-fold helices.
# A trefoil (2,3 torus knot) needs to go around 3 times in one direction
# and 2 times in the other.

# For a torus knot (p_k, q_k) on a torus:
# The number of edges needed scales with p_k * q_k * (circuit length)
# On the Laves graph:

# Approach 1: Volume-based mass
# For torus knot (p_k, q_k) on a torus inscribed in the gyroid:
# Length of torus knot = 2*pi*R * sqrt(p_k^2 + (q_k*r/R)^2)
# On the gyroid, R and r relate to the channel radii

# Let's model it differently:
# The mass of a knot = sum of bicone volumes along the knot path
# times a topological factor (writhe, or crossing number)

# Simple model: mass ~ n_edges * (1 + alpha * crossing_number)
# where alpha captures the interaction energy at crossings

# For unknot: crossings = 0, n = 10 => mass ~ 10
# For trefoil: crossings = 3, and we need the number of edges

# On the Laves graph, a trefoil requires traversing multiple unit cells.
# The trefoil as a torus knot (2,3) on the gyroid:
# needs to wind 2 times along one channel axis and 3 times along another.
# Each winding traverses approximately one 10-ring.
# So: trefoil edges ~ 2 * 10 + 3 * 10 = 50? Too rough.

# Better model: Knot invariants
# The invariant mass of a knot is captured by its hyperbolic volume
# (for hyperbolic knots) or torus knot invariants.

# Torus knot (p_k, q_k) volume in knot complement:
# Vol(S^3 \ K(p,q)) = q * v_3 * sum ... (complicated)
# But for the torus knot, the key invariant is the Jones polynomial.

# Let's compute ratios from torus knot invariants:
print("\n  Torus knot (p_k, q_k) invariants:")
print("  Volume of knot complement and Jones polynomial")

# For torus knots, the Alexander polynomial is:
# Delta(t) = (t^{pq} - 1)(t - 1) / ((t^p - 1)(t^q - 1))
# Evaluated at t = phi (the axiom root):

def alexander_torus_knot_at_phi(pk, qk):
    """Alexander polynomial of torus knot (pk, qk) at t = phi."""
    t = phi
    num = (t**(pk*qk) - 1) * (t - 1)
    den = (t**pk - 1) * (t**qk - 1)
    return num / den

# For torus knots, the bridge number = min(p,q)
# The crossing number = min(p,q) * (max(p,q) - 1)

print(f"\n  {'Knot':<12} {'(p,q)':<10} {'Cross':<8} {'|Alex(phi)|':<20} {'Ratio to (1,1)':<15}")
print(f"  {'='*12} {'='*10} {'='*8} {'='*20} {'='*15}")

base_knots = [(1,1), (2,3), (2,5), (3,5), (2,7), (3,7), (5,7)]
alex_values = {}

for pk, qk in base_knots:
    if pk == qk:
        # (1,1) = unknot
        alex_val = mpf(1)
        crossing = 0
        name = "unknot"
    else:
        alex_val = abs(alexander_torus_knot_at_phi(pk, qk))
        crossing = min(pk,qk) * (max(pk,qk) - 1)
        name = f"T({pk},{qk})"

    alex_values[(pk,qk)] = alex_val
    ratio = alex_val / alex_values.get((1,1), mpf(1))
    print(f"  {name:<12} ({pk},{qk}){'':<5} {crossing:<8} {float(alex_val):<20.6f} {float(ratio):<15.6f}")

print(f"\n  Target ratio m_p/m_e = {float(m_p_over_m_e):.2f}")

# Try different mass models
print("\n  Model 1: mass ~ |Alexander(phi)| * n_edges")
print("    Unknot: 10 edges, |Alex| = 1, mass_e ~ 10")
for pk, qk in base_knots[1:]:
    alex = alex_values[(pk,qk)]
    n_edges_estimate = (pk + qk) * unknot_edges  # rough scaling
    mass_estimate = float(alex * n_edges_estimate)
    ratio = mass_estimate / 10.0
    print(f"    T({pk},{qk}): ~{n_edges_estimate} edges, |Alex| = {float(alex):.2f}, mass ~ {mass_estimate:.2f}, ratio = {ratio:.2f}")

# Model 2: mass ~ crossing_number^a * cycle_length^b
print("\n  Model 2: pure invariant ratios")
print("    Using |Jones(phi)|^2 * dim(colored_rep)")

# Jones polynomial for torus knot (p,q) at t:
# J(t) = (1 - t^2) / (1 - t^{p*q+1}) * sum_{k=0}^{p-1} (-1)^k * t^{k*(k+1)/2 + k*(q-p)} * [p choose k]_t
# This is complex; let me use the simpler relation:
# For the trefoil (2,3): J(t) = -t^{-4} + t^{-3} + t^{-1}
# For (2,5) figure-eight: actually (2,5) is not figure-eight
# (2,5) torus knot Jones: J(t) = -t^{-8} + t^{-7} - t^{-5} + t^{-4} - t^{-2} + t^{-1}

# Let me compute Jones at t = phi for a few small knots
print("\n  Jones polynomial at t = phi:")

# Trefoil (2,3)
t = phi
J_trefoil = -t**(-4) + t**(-3) + t**(-1)
print(f"    Trefoil J(phi) = {float(J_trefoil):.10f}")

# Torus knot (2,5) = Solomon's seal knot (5_1)
J_25 = -t**(-8) + t**(-7) - t**(-5) + t**(-4) - t**(-2) + t**(-1)
print(f"    T(2,5)  J(phi) = {float(J_25):.10f}")

# Figure-eight knot (not a torus knot)
# J(t) = t^2 - t + 1 - t^{-1} + t^{-2}
J_fig8 = t**2 - t + 1 - t**(-1) + t**(-2)
print(f"    Fig-8   J(phi) = {float(J_fig8):.10f}")

# (2,7) torus knot
J_27 = -t**(-12) + t**(-11) - t**(-9) + t**(-8) - t**(-6) + t**(-5) - t**(-3) + t**(-2)
print(f"    T(2,7)  J(phi) = {float(J_27):.10f}")

# Model 3: Knot energy on the gyroid
# The Mobius energy of a knot K: E(K) = integral of (1/|r(s)-r(t)|^2 - 1/d(s,t)^2) ds dt
# For the unknot: E = 4 (Freedman-He-Wang)
# For the trefoil: E ≈ 74.41
# For (2,5): E ≈ 233.0
print("\n  Model 3: Mobius energy ratios")
print("    Unknot E = 4 (exact)")
E_unknot_mobius = 4
E_trefoil_mobius = 74.41  # known from numerical computation
E_fig8_mobius = 88.0      # approximate

print(f"    Trefoil E ~ {E_trefoil_mobius}")
print(f"    Fig-8   E ~ {E_fig8_mobius}")
print(f"    Ratio trefoil/unknot = {E_trefoil_mobius/E_unknot_mobius:.2f}")
print(f"    Ratio fig8/unknot   = {E_fig8_mobius/E_unknot_mobius:.2f}")

# None of these give 1836 directly. Let's look for combinatorial routes.
print("\n  Model 4: Combinatorial (gyroid-specific)")
print("    On the gyroid, a particle = closed path on Laves graph")
print("    Weight = product over edges of |axiom value|")

# For a cycle of length n, if each edge contributes |phi| and each
# vertex contributes |1/2| (the critical value):
# mass ~ phi^n * (1/2)^n = (phi/2)^n
# Or: mass ~ (5/4)^{n/2} (the completed square depth)

print(f"\n    (phi/2)^n weights for various cycles:")
for n in range(1, 25):
    w = float((phi/2)**n)
    print(f"      n={n:>2}: (phi/2)^n = {w:.6f}", end="")
    if abs(w - 1) < 0.01:
        print(" <-- near 1", end="")
    print()

# Try: proton = 10-cycle (same as electron) but with a different
# topological charge from winding number
print(f"\n    Path integral model:")
print(f"    electron = unknot (10-cycle), no winding")
print(f"    proton = wound path that covers the unit cell")
print(f"    proton path length: L_p edges")
print(f"    m_p/m_e = (L_p / L_e) * topological_factor")

# The number of Laves graph edges in 1836 electron-masses
n_proton = m_p_over_m_e * unknot_edges
print(f"\n    If mass ~ edges: proton needs {float(n_proton):.0f} edges")
print(f"    = {float(n_proton / 12):.1f} unit cells of coverage")

# Key insight: what if the proton is a TREFOIL on the gyroid
# and the mass comes from the path length times writhe?
# Trefoil writhe = 3 for the canonical form
# A trefoil on the Laves graph wraps around the helical channel
# Multiple times. The Laves graph has 4_1 screw axes.

# The proton as a (3,5)-torus knot (icosahedral knot!)
# (3,5) has crossing number = 3*(5-1) = 12 = F!!
print(f"\n  The (3,5)-torus knot (icosahedral):")
pk, qk = 3, 5
cross_35 = min(pk,qk) * (max(pk,qk) - 1)
print(f"    Crossing number = {cross_35} = F (dodecahedron faces)")
print(f"    Bridge number   = min(3,5) = {min(pk,qk)} = d")
print(f"    Braid index     = min(3,5) = {min(pk,qk)} = d")

# (3,5) on the Laves graph: each winding covers one 10-cycle
# Total path length: lcm(3,5) * 10 = 15 * 10 = 150?
# Or: 3 * 5 * 10 = 150
n_35 = pk * qk * unknot_edges
print(f"    Path length on Laves ~ p*q*10 = {n_35} edges")

# Alexander polynomial ratio
alex_35 = alex_values[(3,5)]
alex_unknot = alex_values[(1,1)]
print(f"    |Alexander(phi)| for (3,5) = {float(alex_35):.6f}")
print(f"    Ratio to unknot = {float(alex_35):.6f}")

# Volume-weighted: path_length * |Alex(phi)| * crossing_weight
mass_35 = n_35 * alex_35 * cross_35
mass_e_model = unknot_edges * 1 * 1
ratio_35 = mass_35 / mass_e_model
print(f"\n    Composite mass (n * |Alex| * crossings):")
print(f"    Proton: {n_35} * {float(alex_35):.2f} * {cross_35} = {float(mass_35):.2f}")
print(f"    Electron: {unknot_edges} * 1 * 1 = {float(mass_e_model):.2f}")
print(f"    Ratio = {float(ratio_35):.2f}")
print(f"    Target = {float(m_p_over_m_e):.2f}")

# Let's also try: Jones polynomial squared ratio
print(f"\n  Jones polynomial mass model:")
print(f"    |J(phi)|^2 for trefoil = {float(J_trefoil**2):.6f}")
print(f"    |J(phi)|^2 for T(2,5) = {float(J_25**2):.6f}")
print(f"    |J(phi)|^2 for T(2,7) = {float(J_27**2):.6f}")

# A different combinatorial approach:
# The Laves graph has symmetry group I4_132 (body-centered cubic)
# The dodecahedral structure appears through:
# - 12 edges per cell = 12 faces
# - 20 total (V+E) = 20 vertices
# - degree 3 = face-valence

# Try: mass = h^n * SA_bicone^crossings / V_bicone^{beta}
# for appropriate n, crossings, beta

# Actually, let me check the beautiful route:
# m_p/m_e = 6*pi^5 = 1836.12... (close but not exact)
sixpi5 = 6 * mppi**5
print(f"\n  6*pi^5 = {float(sixpi5):.6f}")
print(f"  m_p/m_e = {float(m_p_over_m_e):.6f}")
print(f"  Difference = {float(abs(sixpi5 - m_p_over_m_e)):.6f}")
print(f"  Relative error = {float(abs(sixpi5 - m_p_over_m_e)/m_p_over_m_e):.6e}")

# And: does 6*pi^5 have a gyroid interpretation?
# 6 = V_laves - 2 or 6 = 2*d or 6 = E_laves/2
# pi^5 comes from... the surface area of the 5-sphere?
# S^5 area = 12*pi^3/(5*3*1) = pi^3 * 12/15 = 4*pi^3/5. Nah.
# pi^5: the 5-dimensional sphere volume = 8*pi^2/15
# Actually: V(S^5) = pi^{5/2} / Gamma(7/2) ... different.

# From the gyroid: 12 bicones, each with SA = 5*pi*sqrt(2)/2
# Product: SA^5 * scaling?

# Let's try:
# m_p/m_e from the gyroid
# = (SA_12 / (2*pi))^{5/2} or something
print(f"\n  Gyroid-derived mass ratio attempts:")
attempt1 = (SA_12 / (2 * mppi))**(mpf(5)/2)
print(f"    (12*SA/(2pi))^(5/2) = {float(attempt1):.2f}")

attempt2 = (V_12_exact / mppi)**(mpf(5)/3)
print(f"    (12*V/pi)^(5/3) = {float(attempt2):.2f}")

# Beautiful attempt: m_p/m_e = E_laves * V_laves * (SA_bicone/pi)^2
# = 12 * 8 * (5*sqrt(2)/2)^2 = 96 * 25/2 = 96 * 12.5 = 1200. Nope.
attempt3 = E_laves * V_laves * (SA_bicone/mppi)**2
print(f"    E*V*(SA/pi)^2 = {float(attempt3):.2f}")

# Try: m_p/m_e = (5^5 * phi^3) / d = (3125 * phi^3) / 3
attempt4 = mpf(5)**5 * phi**3 / d
print(f"    5^5 * phi^3 / d = {float(attempt4):.2f}")

# 5^5 * phi^2 / (3/2) = ?
attempt5 = mpf(5)**5 * phi**2 / (mpf(3)/2)
print(f"    5^5 * phi^2 / 1.5 = {float(attempt5):.2f}")

# Final model: the (3,5) torus knot IS icosahedral.
# In the McKay correspondence, the binary icosahedral group 2I
# has order 120, and the (3,5) torus knot is THE icosahedral knot.
# The proton, as the simplest baryon (3 quarks), maps to the (3,5) torus knot
# where 3 = number of quarks and 5 = Schlafli parameter.

print(f"\n  *** THE ICOSAHEDRAL KNOT MODEL ***")
print(f"  Proton = (3,5) torus knot on the gyroidal Laves graph")
print(f"  Electron = unknot (10-ring)")
print(f"  3 = quarks = vertex degree d")
print(f"  5 = Schlafli p = color charges")
print(f"  crossing_number(3,5) = 12 = faces of dodecahedron")
print(f"  bridge_number(3,5) = 3 = d")

# The (3,5) torus knot genus: g = (p-1)(q-1)/2 = 2*4/2 = 4
genus_35 = (pk - 1) * (qk - 1) // 2
print(f"  Genus of (3,5) torus knot = {genus_35}")
print(f"  Euler char of Seifert surface = {1 - 2*genus_35}")

# Seifert genus * binary_icosahedral / electron_loop
seifert_ratio = genus_35 * binary_icosahedral_order / unknot_edges
print(f"  genus * |2I| / unknot = {genus_35} * {binary_icosahedral_order} / {unknot_edges} = {seifert_ratio}")
print(f"  Not 1836 but interesting.")

# OK let me try the actual known approximation:
# m_p/m_e ≈ (1/alpha) * 3/(2*pi) * phi * sqrt(5) * something
alpha_inv = mpf('137.035999084')
alpha_em = 1 / alpha_inv
print(f"\n  Using alpha = 1/{float(alpha_inv):.6f}:")
print(f"  m_p/m_e / (1/alpha) = {float(m_p_over_m_e * alpha_em):.6f}")
print(f"  = {float(m_p_over_m_e / alpha_inv):.6f}")
ratio_to_alpha = m_p_over_m_e / alpha_inv
print(f"  m_p/m_e = (1/alpha) * {float(ratio_to_alpha):.10f}")
print(f"  {float(ratio_to_alpha):.10f} vs 4*pi = {float(4*mppi):.10f}")
print(f"  vs 3*pi*sqrt(phi) = {float(3*mppi*mpsqrt(phi)):.10f}")

# ============================================================
# SUMMARY TABLE
# ============================================================
print("\n" + "=" * 72)
print("  SUMMARY: GYROIDAL BICONE CORRESPONDENCES")
print("=" * 72)

print(f"""
  AXIOM: x^2 - x - 1 = 0
  COMPLETION: (x - 1/2)^2 = 5/4

  BICONE (one axiom instance):
    Height h           = sqrt(5)        = axiom diagonal
    Equator radius r   = sqrt(5)/2      = half-diagonal
    Equator depth      = 5/4            = p/4
    Volume             = 5*sqrt(5)*pi/6
    Surface area       = 5*pi*sqrt(2)/2
    Slant height       = sqrt(5/2)      = sqrt(p/2)
    r^2 + (h/2)^2      = 5/4 + 5/4 = 5/2  (SYMMETRIC!)

  GYROID UNIT CELL (12 bicones):
    Vertices           = 8  = 2^d
    Edges (bicones)    = 12 = F_dodecahedron
    Degree             = 3  = d = q
    V + E              = 20 = V_dodecahedron
    V - E              = -4 = chi_gyroid
    b1 (Betti)         = 5  = p
    Shortest cycle     = 10 = 2*p
    Genus per cell     = 3  = d

  GAMMA CONNECTION:
    gamma * sqrt(d) = {float(gamma_sqrt_d):.15f}
    Deviation from 1: {float(abs(1-gamma_sqrt_d)):.6e}

  PARTICLE MODEL:
    Electron           = unknot (10-ring on Laves graph)
    Proton             = (3,5)-torus knot (icosahedral knot!)
      crossing_number  = 12 = F
      bridge_number    = 3  = d
      Seifert genus    = 4

  KEY NUMEROLOGY:
    6*pi^5             = {float(sixpi5):.6f}
    m_p/m_e            = {float(m_p_over_m_e):.6f}
    Error              = {float(abs(sixpi5 - m_p_over_m_e)/m_p_over_m_e):.2e}
""")

print("=" * 72)
print("  END OF COMPUTATION")
print("=" * 72)
