#!/usr/bin/env python3
"""
HUNT_NONGOLDEN — hunts the correction mechanism for non-golden spectral sector of Ihara zeta at critical line
nos3bl33d

Dodecahedron Tr(L0^-1)=137/15. Convention: F(u) = u d/du log Z_G^-1(u). mpmath 80 digits.
"""

from mpmath import mp, mpf, sqrt, log, pi, fabs, nstr, findroot

mp.dps = 80

# =============================================================================
# CONSTANTS
# =============================================================================
d = 3                  # degree
V = 20                 # vertices
E = 30                 # edges
sqrt3 = sqrt(3)
sqrt5 = sqrt(5)
sqrt15 = sqrt(15)
u0 = 1 / sqrt3         # critical point u = 1/sqrt3
target = mpf(137) / 15  # Tr(L0^{-1})

# Adjacency eigenvalues and multiplicities
eigenvalues = [
    (mpf(3),     1, "trivial (lam=3)"),
    (sqrt5,      3, "golden+ (lam=+sqrt5)"),
    (mpf(1),     5, "chi5 (lam=1)"),
    (mpf(0),     4, "chi4 (lam=0)"),
    (-sqrt5,     3, "golden- (lam=-sqrt5)"),
    (mpf(-2),    4, "chi4' (lam=-2)"),
]

golden_indices  = [1, 4]   # +/-sqrt5
nongolden_indices = [0, 2, 3, 5]  # 3, 1, 0, -2

SEP = "=" * 90

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def f_contrib(u, lam, c):
    """
    Single eigenvalue contribution to F(u) = u d/du log Z_G^{-1}(u).

    Z^{-1}(u) = (1-u^2)^{E-V} prod(1 - lam_i u + c u^2)
    d/du log(1 - lam u + cu^2) = (-lam + 2cu)/(1 - lam u + cu^2)
    u * d/du log(1-lam u+cu^2) = u(2cu - lam)/(1 - lam u + cu^2)

    At u=1/sqrt3, c=2: f = (4 - lam*sqrt3)/(5 - lam*sqrt3).
    """
    num = u * (2 * c * u - lam)
    den = 1 - lam * u + c * u**2
    return num / den


def edge_term(u):
    """Edge contribution: (E-V) * u * d/du log(1-u^2) = -(E-V)*2u^2/(1-u^2)."""
    return -(E - V) * 2 * u**2 / (1 - u**2)


def F_total(u, c_golden=2, c_nongolden=2, include_edge=True):
    """Total F(u) with possibly different coefficients for golden vs non-golden."""
    total = mpf(0)
    if include_edge:
        total += edge_term(u)
    for i, (lam, mult, name) in enumerate(eigenvalues):
        c = c_golden if i in golden_indices else c_nongolden
        total += mult * f_contrib(u, lam, c)
    return total


# =============================================================================
# PART 1: EXACT DECOMPOSITION AT s=1/2 (standard Ihara, c = d-1 = 2)
# =============================================================================
print(SEP)
print("PART 1: EXACT DECOMPOSITION AT s=1/2 (standard Ihara, c = d-1 = 2)")
print(SEP)

edge_val = edge_term(u0)
print(f"\nEdge term at u=1/sqrt3: {nstr(edge_val, 30)}")
print(f"  Expected: -10 (exact)")

print(f"\nEigenvalue contributions with c = d-1 = 2:")
total_c2 = edge_val
golden_total = mpf(0)
nongolden_total = mpf(0)

for i, (lam, mult, name) in enumerate(eigenvalues):
    single = f_contrib(u0, lam, 2)
    contrib = mult * single
    total_c2 += contrib
    label = "GOLDEN" if i in golden_indices else "nongolden"
    if i in golden_indices:
        golden_total += contrib
    else:
        nongolden_total += contrib
    print(f"  {name:30s} mult={mult}  single={nstr(single, 25):>35s}  "
          f"total={nstr(contrib, 25):>35s}  [{label}]")

print(f"\n  Edge term:              {nstr(edge_val, 30)}")
print(f"  Golden subtotal:        {nstr(golden_total, 30)}")
print(f"  Non-golden subtotal:    {nstr(nongolden_total, 30)}")
print(f"  GRAND TOTAL F(1/sqrt3): {nstr(total_c2, 30)}")
print(f"  Target 137/15:          {nstr(target, 30)}")
print(f"  Gap F - 137/15:         {nstr(total_c2 - target, 30)}")
print(f"\n  Golden pair = 3?  |golden - 3| = {nstr(fabs(golden_total - 3), 15)}")

# =============================================================================
# EXACT Q(sqrt3) arithmetic verification
# =============================================================================
print("\n" + "-" * 60)
print("EXACT Q(sqrt3) VERIFICATION")
print("-" * 60)
print("\nUsing f(u,lam,c) = u(2cu-lam)/(1-lam*u+cu^2)")
print("At u=1/sqrt3, c=2: f = (4-lam*sqrt3)/(5-lam*sqrt3)")

# lam=3: f = (4-3sqrt3)/(5-3sqrt3)
# Rationalize by (5+3sqrt3): (4-3sqrt3)(5+3sqrt3)/((5)^2-(3sqrt3)^2) = (4-3sqrt3)(5+3sqrt3)/(25-27)
# Numerator: 20+12sqrt3-15sqrt3-27 = -7-3sqrt3
# Denominator: -2
# f = (-7-3sqrt3)/(-2) = (7+3sqrt3)/2
exact_lam3 = (7 + 3*sqrt3) / 2
check_lam3 = f_contrib(u0, mpf(3), 2)
print(f"\n  lam=3: exact (7+3sqrt3)/2 = {nstr(exact_lam3, 25)}")
print(f"  lam=3: computed            = {nstr(check_lam3, 25)}")
print(f"  Match: {nstr(fabs(exact_lam3 - check_lam3), 15)}")

# lam=1: f = (4-sqrt3)/(5-sqrt3)
# Rationalize by (5+sqrt3): (4-sqrt3)(5+sqrt3)/(25-3) = (20+4sqrt3-5sqrt3-3)/22 = (17-sqrt3)/22
exact_lam1 = (17 - sqrt3) / 22
check_lam1 = f_contrib(u0, mpf(1), 2)
print(f"\n  lam=1: exact (17-sqrt3)/22 = {nstr(exact_lam1, 25)}")
print(f"  lam=1: computed             = {nstr(check_lam1, 25)}")
print(f"  lam=1: x5 contribution      = {nstr(5*exact_lam1, 25)}")

# lam=0: f = (4-0)/(5-0) = 4/5
exact_lam0 = mpf(4) / 5
check_lam0 = f_contrib(u0, mpf(0), 2)
print(f"\n  lam=0: exact 4/5 = {nstr(exact_lam0, 25)}")
print(f"  lam=0: computed   = {nstr(check_lam0, 25)}")
print(f"  lam=0: x4 contribution = {nstr(4*exact_lam0, 25)}")

# lam=-2: f = (4+2sqrt3)/(5+2sqrt3)
# Rationalize by (5-2sqrt3): (4+2sqrt3)(5-2sqrt3)/(25-12) = (20-8sqrt3+10sqrt3-12)/13 = (8+2sqrt3)/13
exact_lamm2 = (8 + 2*sqrt3) / 13
check_lamm2 = f_contrib(u0, mpf(-2), 2)
print(f"\n  lam=-2: exact (8+2sqrt3)/13 = {nstr(exact_lamm2, 25)}")
print(f"  lam=-2: computed             = {nstr(check_lamm2, 25)}")
print(f"  lam=-2: x4 contribution      = {nstr(4*exact_lamm2, 25)}")

# Golden pair: lam=+/-sqrt5
# lam=sqrt5: f = (4-sqrt15)/(5-sqrt15)
# lam=-sqrt5: f = (4+sqrt15)/(5+sqrt15)
exact_gp = (4 - sqrt15) / (5 - sqrt15)
exact_gm = (4 + sqrt15) / (5 + sqrt15)
print(f"\n  lam=+sqrt5: single = {nstr(exact_gp, 25)}")
print(f"  lam=-sqrt5: single = {nstr(exact_gm, 25)}")
golden_sum_exact = 3*exact_gp + 3*exact_gm
print(f"  Golden total (3x each): {nstr(golden_sum_exact, 25)}")
print(f"  Expected: exactly 3")
print(f"  |golden - 3| = {nstr(fabs(golden_sum_exact - 3), 15)}")

# Verify golden = 3 algebraically:
# 3[(4-sqrt15)/(5-sqrt15) + (4+sqrt15)/(5+sqrt15)]
# = 3[(4-sqrt15)(5+sqrt15) + (4+sqrt15)(5-sqrt15)] / [(5-sqrt15)(5+sqrt15)]
# = 3[(20+4sqrt15-5sqrt15-15) + (20-4sqrt15+5sqrt15-15)] / (25-15)
# = 3[(5-sqrt15) + (5+sqrt15)] / 10 = 3[10]/10 = 3  QED

# Reconstruct F from exact pieces
nongolden_eigenvalue_sum = (1*exact_lam3 + 5*exact_lam1 + 4*exact_lam0 + 4*exact_lamm2)
F_exact = edge_val + golden_sum_exact + nongolden_eigenvalue_sum
print(f"\n  Non-golden eigenvalue sum: {nstr(nongolden_eigenvalue_sum, 30)}")
print(f"  F(1/sqrt3) from exact:    {nstr(F_exact, 30)}")
print(f"  F(1/sqrt3) from function: {nstr(F_total(u0), 30)}")
print(f"  Match: {nstr(fabs(F_exact - F_total(u0)), 15)}")

# Verify Q(sqrt3) form: F = 4308/715 + 270*sqrt3/143
F_Qsqrt3 = mpf(4308)/715 + mpf(270)/143 * sqrt3
print(f"\n  Expected Q(sqrt3) form: 4308/715 + 270*sqrt3/143 = {nstr(F_Qsqrt3, 30)}")
print(f"  Computed F:                                          {nstr(F_exact, 30)}")
print(f"  Match: {nstr(fabs(F_exact - F_Qsqrt3), 15)}")

# Compute exact non-golden (WITH edge)
nongolden_with_edge = edge_val + nongolden_eigenvalue_sum
print(f"\n  Non-golden + edge: {nstr(nongolden_with_edge, 25)}")
print(f"  Golden:            {nstr(golden_sum_exact, 25)}")
print(f"  Sum:               {nstr(nongolden_with_edge + golden_sum_exact, 25)}")

# =============================================================================
# PART 2: THE GAP ANALYSIS
# =============================================================================
print(f"\n{SEP}")
print("PART 2: THE GAP ANALYSIS")
print(SEP)

gap = F_total(u0) - target
nongolden_F = F_total(u0) - golden_total  # golden = 3
nongolden_target = target - 3   # = 92/15

print(f"\n  F(1/sqrt3)           = {nstr(F_total(u0), 30)}")
print(f"  137/15               = {nstr(target, 30)}")
print(f"  Gap                  = {nstr(gap, 30)}")
print(f"  92/15                = {nstr(mpf(92)/15, 30)}")
print(f"  Non-golden F (incl edge) = {nstr(nongolden_F, 30)}")
print(f"  Non-golden target    = {nstr(nongolden_target, 30)}")
print(f"  Non-golden excess    = {nstr(nongolden_F - nongolden_target, 30)}")

# =============================================================================
# PART 3: LAPLACIAN vs ADJACENCY
# =============================================================================
print(f"\n{SEP}")
print("PART 3: LAPLACIAN vs ADJACENCY PERSPECTIVE")
print(SEP)

lap_eigenvalues = [
    (3 - sqrt5, 3, "golden: mu=3-sqrt5"),
    (mpf(2),    5, "mu=2"),
    (mpf(3),    4, "mu=3"),
    (3 + sqrt5, 3, "golden: mu=3+sqrt5"),
    (mpf(5),    4, "mu=5"),
]

print(f"\n  Tr(L0^{{-1}}) = sum mult_i/mu_i (mu!=0):")
tr_lap = mpf(0)
golden_lap = mpf(0)
nongolden_lap = mpf(0)
for mu, mult, name in lap_eigenvalues:
    contrib = mult / mu
    tr_lap += contrib
    is_golden = "golden" in name
    if is_golden:
        golden_lap += contrib
    else:
        nongolden_lap += contrib
    print(f"    {name:25s}  mult={mult}  1/mu={nstr(1/mu, 20):>25s}  "
          f"contrib={nstr(contrib, 20):>25s}")

print(f"\n  Tr(L0^{{-1}}) = {nstr(tr_lap, 30)}")
print(f"  137/15       = {nstr(target, 30)}")
print(f"  Match: {nstr(fabs(tr_lap - target), 15)}")

print(f"\n  Golden Laplacian = 9/2? |diff| = {nstr(fabs(golden_lap - mpf(9)/2), 15)}")
print(f"  Non-golden Laplacian = 139/30? |diff| = {nstr(fabs(nongolden_lap - mpf(139)/30), 15)}")

# Shifts
golden_shift = golden_total - golden_lap  # F_golden(=3) - Lap_golden(=9/2) = -3/2
nongolden_shift = nongolden_F - nongolden_lap  # includes edge in F

print(f"\n  SHIFTS (Ihara F - Laplacian Tr):")
print(f"    Golden shift:     {nstr(golden_shift, 25)}  (should be -d/2 = -3/2)")
print(f"    Non-golden shift: {nstr(nongolden_shift, 25)}")
print(f"    Total shift:      {nstr(golden_shift + nongolden_shift, 25)}  (= gap)")
print(f"    Gap:              {nstr(gap, 25)}")
print(f"\n  Non-golden shift - d/2 = {nstr(nongolden_shift - mpf(3)/2, 25)}")
print(f"  gap = golden_shift + nongolden_shift = {nstr(golden_shift + nongolden_shift, 25)}")

# =============================================================================
# PART 4: BACKTRACKING CORRECTION c=2 -> c=1
# =============================================================================
print(f"\n{SEP}")
print("PART 4: BACKTRACKING CORRECTION (c=2 -> c=1 for ALL)")
print(SEP)

print(f"\n  Eigenvalue contributions with c=1 (Artin coefficient):")
total_c1 = edge_term(u0)
golden_c1 = mpf(0)
nongolden_c1 = mpf(0)

for i, (lam, mult, name) in enumerate(eigenvalues):
    single = f_contrib(u0, lam, 1)
    contrib = mult * single
    total_c1 += contrib
    if i in golden_indices:
        golden_c1 += contrib
    else:
        nongolden_c1 += contrib
    print(f"    {name:30s} mult={mult}  single={nstr(single, 25):>35s}  "
          f"total={nstr(contrib, 25):>35s}")

print(f"\n  Edge term:           {nstr(edge_term(u0), 30)}")
print(f"  Golden (c=1):        {nstr(golden_c1, 30)}")
print(f"  Non-golden (c=1):    {nstr(nongolden_c1, 30)}")
print(f"  TOTAL F(c=1):        {nstr(total_c1, 30)}")
print(f"  Target 137/15:       {nstr(target, 30)}")
print(f"  Gap F(c=1)-137/15:   {nstr(total_c1 - target, 30)}")

# =============================================================================
# PART 5: MIXED COEFFICIENT SCENARIOS
# =============================================================================
print(f"\n{SEP}")
print("PART 5: MIXED COEFFICIENT SCENARIOS")
print(SEP)

scenarios = [
    (2, 2, "Standard Ihara (c=2 all)"),
    (1, 1, "Artin (c=1 all)"),
    (2, 1, "Golden c=2, Rest c=1"),
    (1, 2, "Golden c=1, Rest c=2"),
    (0, 0, "c=0 all"),
    (2, 0, "Golden c=2, Rest c=0"),
    (0, 2, "Golden c=0, Rest c=2"),
]

print(f"\n  {'Scenario':45s}  {'F(1/sqrt3)':>30s}  {'F - 137/15':>25s}")
print(f"  {'-'*45}  {'-'*30}  {'-'*25}")
for c_g, c_ng, desc in scenarios:
    val = F_total(u0, c_golden=c_g, c_nongolden=c_ng)
    print(f"  {desc:45s}  {nstr(val, 25):>30s}  {nstr(val - target, 20):>25s}")

print(f"\n  WITHOUT EDGE TERM:")
print(f"  {'Scenario':45s}  {'F(1/sqrt3)':>30s}  {'F - 137/15':>25s}")
print(f"  {'-'*45}  {'-'*30}  {'-'*25}")
for c_g, c_ng, desc in scenarios:
    val = F_total(u0, c_golden=c_g, c_nongolden=c_ng, include_edge=False)
    print(f"  {desc + ' [no edge]':45s}  {nstr(val, 25):>30s}  {nstr(val - target, 20):>25s}")

# =============================================================================
# PART 6: SOLVE FOR c*
# =============================================================================
print(f"\n{SEP}")
print("PART 6: SOLVE FOR c* (coefficient that closes the gap)")
print(SEP)

def nongolden_sum(c_ng):
    """Sum of non-golden eigenvalue contributions at u=1/sqrt3 with coefficient c_ng."""
    s = mpf(0)
    for i in nongolden_indices:
        lam, mult, _ = eigenvalues[i]
        s += mult * f_contrib(u0, lam, c_ng)
    return s

# 6a: Fix golden at c=2 (golden=3 exactly), find c_ng so total = 137/15
# Need: edge + 3 + nongolden_sum(c_ng) = 137/15
# nongolden_sum(c_ng) = 137/15 - edge - 3 = 137/15 + 10 - 3 = 137/15 + 7 = 242/15
target_ng_sum = mpf(242) / 15
print(f"\n  Target for non-golden eigenvalue sum: {nstr(target_ng_sum, 25)}")
print(f"  Current (c=2): {nstr(nongolden_sum(2), 25)}")

def gap_func_ng(c):
    return nongolden_sum(c) - target_ng_sum

# Newton search starting near c=2
try:
    c_ng_star = findroot(gap_func_ng, mpf('1.9'))
    print(f"  c*_nongolden (golden at c=2): {nstr(c_ng_star, 40)}")
    print(f"  Verify: F(c_g=2, c_ng=c*) = {nstr(F_total(u0, c_golden=2, c_nongolden=c_ng_star), 30)}")
    print(f"  c* - 2 = {nstr(c_ng_star - 2, 30)}")
    print(f"  2 - c* = {nstr(2 - c_ng_star, 30)}")
except Exception as e:
    print(f"  Could not find c*_nongolden: {e}")
    # Manual search
    print("  Manual scan:")
    for k in range(150, 210):
        c_test = mpf(k) / 100
        v = nongolden_sum(c_test)
        if fabs(v - target_ng_sum) < 0.5:
            print(f"    c={nstr(c_test,5)}  sum={nstr(v,15)}  gap={nstr(v-target_ng_sum,10)}")

# 6b: Fix golden at c=1, find c_ng
golden_at_c1 = mpf(0)
for i in golden_indices:
    lam, mult, _ = eigenvalues[i]
    golden_at_c1 += mult * f_contrib(u0, lam, 1)

target_ng_v2 = target - edge_val - golden_at_c1
print(f"\n  Golden at c=1: {nstr(golden_at_c1, 30)}")
print(f"  Target ng sum (golden c=1): {nstr(target_ng_v2, 25)}")

try:
    c_ng_star2 = findroot(lambda c: nongolden_sum(c) - target_ng_v2, mpf('1.5'))
    print(f"  c*_nongolden (golden at c=1): {nstr(c_ng_star2, 40)}")
    print(f"  Verify: F = {nstr(F_total(u0, c_golden=1, c_nongolden=c_ng_star2), 30)}")
except Exception as e:
    print(f"  Could not find: {e}")

# 6c: Fix non-golden at c=2, find c_g
def golden_sum_func(c_g):
    s = mpf(0)
    for i in golden_indices:
        lam, mult, _ = eigenvalues[i]
        s += mult * f_contrib(u0, lam, c_g)
    return s

golden_target = target - edge_val - nongolden_eigenvalue_sum
print(f"\n  Golden target (non-golden at c=2): {nstr(golden_target, 30)}")
print(f"  Current golden (c=2): {nstr(golden_sum_func(2), 25)}")

try:
    c_g_star = findroot(lambda c: golden_sum_func(c) - golden_target, mpf('2.5'))
    print(f"  c*_golden (non-golden at c=2): {nstr(c_g_star, 40)}")
    print(f"  Verify: F = {nstr(F_total(u0, c_golden=c_g_star, c_nongolden=2), 30)}")
except Exception as e:
    print(f"  Could not find: {e}")

# =============================================================================
# PART 7: PER-EIGENVALUE c*
# =============================================================================
print(f"\n{SEP}")
print("PART 7: PER-EIGENVALUE c* (what c closes gap for each eigenvalue alone?)")
print(SEP)

print(f"\n  Fix all at c=2 except one eigenvalue class. Find c* that closes gap:")

for idx in nongolden_indices:
    lam, mult, name = eigenvalues[idx]
    # Current contribution at c=2
    current = mult * f_contrib(u0, lam, 2)
    # We need this to decrease by exactly gap (F - 137/15)
    needed = current - gap

    # f(u0, lam, c) = u0*(2cu0 - lam)/(1 - lam*u0 + c*u0^2)
    # For single eigenvalue: solve mult * f(u0, lam, c) = needed
    target_single = needed / mult

    # u0*(2c*u0 - lam)/(1 - lam*u0 + c*u0^2) = target_single
    # Let A = 2*u0^2, B = -u0*lam, C = 1-lam*u0, D = u0^2
    # (Ac + B)/(C + Dc) = target_single
    # Ac + B = target_single*(C + Dc)
    # Ac - target_single*Dc = target_single*C - B
    # c(A - target_single*D) = target_single*C - B
    # c = (target_single*C - B)/(A - target_single*D)

    A = 2 * u0**2
    B = -u0 * lam
    C = 1 - lam * u0
    D = u0**2

    c_star_i = (target_single * C - B) / (A - target_single * D)

    # Verify
    verify_single = mult * f_contrib(u0, lam, c_star_i)
    # Total: replace this eigenvalue's contribution
    F_other = F_total(u0) - current
    verify_total = F_other + verify_single

    print(f"\n  {name}:")
    print(f"    Current contrib (c=2): {nstr(current, 25)}")
    print(f"    Needed contrib:        {nstr(needed, 25)}")
    print(f"    c* for this eigenvalue: {nstr(c_star_i, 40)}")
    print(f"    c* - 2 =                {nstr(c_star_i - 2, 25)}")
    print(f"    Verify total = {nstr(verify_total, 25)} (target: {nstr(target, 20)})")

# =============================================================================
# PART 8: BRIDGE POINT -- where does F(u, c=1) = 137/15?
# =============================================================================
print(f"\n{SEP}")
print("PART 8: BRIDGE POINT -- find u1 where F(u, c=1) = 137/15")
print(SEP)

def F_c1(u):
    return F_total(u, c_golden=1, c_nongolden=1, include_edge=True)

def bridge_eq(u):
    return F_c1(u) - target

print("\n  Scanning F(u, c=1) near critical region:")
for k in range(-10, 20):
    u_test = u0 + mpf(k)/100
    if u_test > mpf('0.01') and u_test < mpf('0.99'):
        try:
            val = F_c1(u_test)
            marker = " <--- NEAR TARGET" if fabs(val - target) < 1 else ""
            if fabs(val) < 100:
                print(f"    u = {nstr(u_test, 12):>15s}  F(u,c=1) = {nstr(val, 15):>20s}{marker}")
        except Exception:
            pass

try:
    u1 = findroot(bridge_eq, u0 * mpf('0.95'))
    print(f"\n  Bridge point u1: {nstr(u1, 40)}")
    print(f"  1/sqrt3 =         {nstr(u0, 40)}")
    print(f"  u1 - 1/sqrt3 =    {nstr(u1 - u0, 25)}")
    print(f"  u1/u0 =           {nstr(u1/u0, 25)}")
    print(f"  Verify F(u1,c=1) = {nstr(F_c1(u1), 30)}")

    s1_log3 = -log(u1) / log(3)
    s1_log2 = -log(u1) / log(2)
    print(f"\n  s1 = -log(u1)/log(3) = {nstr(s1_log3, 25)}")
    print(f"  s1 = -log(u1)/log(2) = {nstr(s1_log2, 25)}")
    print(f"  s0 = 1/2")
    print(f"  |s1 - 1/2| (base 3) = {nstr(fabs(s1_log3 - mpf(1)/2), 25)}")
except Exception as e:
    print(f"\n  Could not find bridge point near u0: {e}")
    print("  Broader scan:")
    for k in range(1, 99):
        u_test = mpf(k) / 100
        try:
            val = F_c1(u_test)
            if fabs(val - target) < 1:
                print(f"    u = {nstr(u_test, 8)}  F = {nstr(val, 15)}  gap = {nstr(val-target, 10)}")
        except Exception:
            pass

# =============================================================================
# PART 9: FUNCTIONAL EQUATION
# =============================================================================
print(f"\n{SEP}")
print("PART 9: FUNCTIONAL EQUATION with c=1")
print(SEP)

print(f"\n  Checking F(u, c=1) vs F(1/(2u), c=1):")
for u_test in [u0, mpf(1)/4, mpf(1)/3, mpf(3)/10]:
    u_dual = 1/(2*u_test)
    if mpf('0.01') < u_dual < mpf('0.99'):
        try:
            f1 = F_c1(u_test)
            f2 = F_c1(u_dual)
            print(f"    u={nstr(u_test,10):>12s}  1/(2u)={nstr(u_dual,10):>12s}  "
                  f"F(u)={nstr(f1,15):>20s}  F(1/(2u))={nstr(f2,15):>20s}  "
                  f"diff={nstr(f1-f2,10)}")
        except Exception:
            pass

print(f"\n  Checking F(u, c=1) vs F(1/(3u), c=1):")
for u_test in [u0, mpf(1)/4, mpf(1)/3, mpf(3)/10, mpf(2)/5]:
    u_dual = 1/(3*u_test)
    if mpf('0.01') < u_dual < mpf('0.99'):
        try:
            f1 = F_c1(u_test)
            f2 = F_c1(u_dual)
            print(f"    u={nstr(u_test,10):>12s}  1/(3u)={nstr(u_dual,10):>12s}  "
                  f"F(u)={nstr(f1,15):>20s}  F(1/(3u))={nstr(f2,15):>20s}  "
                  f"diff={nstr(f1-f2,10)}")
        except Exception:
            pass

# =============================================================================
# PART 10: QUADRATIC FACTORIZATION (Frobenius roots)
# =============================================================================
print(f"\n{SEP}")
print("PART 10: QUADRATIC FACTORIZATION -- Frobenius roots")
print(SEP)

# Roots of t^2 - lam*t + (d-1) = 0: alpha,beta with alpha*beta = d-1 = 2
print(f"\n  Ramanujan bound: |lam| <= 2*sqrt(d-1) = 2*sqrt2 = {nstr(2*sqrt(mpf(2)), 15)}")
print(f"\n  Factorization 1 - lam*u + (d-1)u^2 = (1-alpha*u)(1-beta*u):")
print(f"  alpha+beta = lam, alpha*beta = d-1 = 2")

reps = [
    ("chi1 (trivial)", 1, mpf(3)),
    ("chi3 (golden+)", 3, sqrt5),
    ("chi5",           5, mpf(1)),
    ("chi4",           4, mpf(0)),
    ("chi3' (golden-)",3, -sqrt5),
    ("chi4'",          4, mpf(-2)),
]

half_A = edge_term(u0)  # Strategy A: edge + alpha only for non-Ram, half for Ram
half_B = edge_term(u0)  # Strategy B: edge + beta only for non-Ram, half for Ram

for rep_name, dim, lam in reps:
    disc = lam**2 - 4*(d-1)
    if disc >= 0:
        sd = sqrt(disc)
        alpha = (lam + sd) / 2
        beta  = (lam - sd) / 2
        ram = "NON-Ramanujan"

        # Frobenius root contributions: alpha*u/(1-alpha*u)
        f_a = alpha * u0 / (1 - alpha * u0) if fabs(1 - alpha*u0) > 1e-50 else mpf('inf')
        f_b = beta * u0 / (1 - beta * u0) if fabs(1 - beta*u0) > 1e-50 else mpf('inf')

        print(f"  {rep_name:20s}  lam={nstr(lam,6):>8s}  {ram}  "
              f"alpha={nstr(alpha,10):>15s}  beta={nstr(beta,10):>15s}  "
              f"f_a={nstr(f_a,12):>15s}  f_b={nstr(f_b,12):>15s}  "
              f"f_a+f_b={nstr(f_a+f_b,12):>15s}  f(lam,c=2)={nstr(f_contrib(u0,lam,2),12):>15s}")

        half_A += dim * f_a
        half_B += dim * f_b
    else:
        sd = sqrt(-disc)
        ram = "Ramanujan"
        f_full = f_contrib(u0, lam, 2)
        print(f"  {rep_name:20s}  lam={nstr(lam,6):>8s}  {ram}      "
              f"|alpha|=sqrt(d-1)=sqrt2  f_full={nstr(f_full,15):>20s}  "
              f"f_half={nstr(f_full/2,15):>20s}")
        half_A += dim * f_full / 2
        half_B += dim * f_full / 2

print(f"\n  Strategy A (alpha for non-Ram, half for Ram): {nstr(half_A, 30)}")
print(f"  Strategy B (beta for non-Ram, half for Ram):  {nstr(half_B, 30)}")
print(f"  Target 137/15: {nstr(target, 30)}")
print(f"  Gap A: {nstr(half_A - target, 25)}")
print(f"  Gap B: {nstr(half_B - target, 25)}")

# =============================================================================
# PART 11: REMOVE TRIVIAL EIGENVALUE
# =============================================================================
print(f"\n{SEP}")
print("PART 11: REMOVE TRIVIAL EIGENVALUE lam=3")
print(SEP)

trivial_ihara = f_contrib(u0, mpf(3), 2)
F_no_triv = F_total(u0) - trivial_ihara
F_no_triv_no_edge = F_total(u0, include_edge=False) - trivial_ihara

print(f"\n  Trivial (lam=3) contribution: {nstr(trivial_ihara, 25)}")
print(f"  F without trivial:            {nstr(F_no_triv, 25)}")
print(f"  F without trivial or edge:    {nstr(F_no_triv_no_edge, 25)}")
print(f"  Target: {nstr(target, 25)}")

# =============================================================================
# PART 12: GAP DECOMPOSITION -- Ihara vs Laplacian per eigenvalue
# =============================================================================
print(f"\n{SEP}")
print("PART 12: FULL GAP DECOMPOSITION")
print(SEP)

print(f"\n  Per-eigenvalue: Ihara vs Laplacian 1/mu:")
print(f"  {'Name':>30s}  {'mult':>4s}  {'Ihara(xmult)':>25s}  {'Lap 1/mu(xmult)':>25s}  {'Diff':>25s}")

golden_ihara_val = golden_total
golden_lap_val = mpf(9)/2

# Non-golden, non-trivial: lam in {1, 0, -2}
ng_nt_ihara = mpf(0)
ng_nt_lap = mpf(0)

for i in [2, 3, 5]:  # lam=1, 0, -2
    lam, mult, name = eigenvalues[i]
    ihara_c = mult * f_contrib(u0, lam, 2)
    mu = 3 - lam
    lap_c = mult / mu
    diff = ihara_c - lap_c
    ng_nt_ihara += ihara_c
    ng_nt_lap += lap_c
    print(f"  {name:>30s}  {mult:>4d}  {nstr(ihara_c,20):>25s}  {nstr(lap_c,20):>25s}  {nstr(diff,20):>25s}")

print(f"\n  Trivial lam=3 (mu=0): Ihara={nstr(trivial_ihara,20)}, Laplacian=EXCLUDED")

print(f"\n  FULL GAP DECOMPOSITION:")
print(f"    Edge term:                     {nstr(edge_val, 25)}")
print(f"    Trivial (lam=3):               {nstr(trivial_ihara, 25)}")
print(f"    Golden (Ihara - Lap):          {nstr(golden_ihara_val - golden_lap_val, 25)}")
print(f"    Non-g non-t (Ihara - Lap):     {nstr(ng_nt_ihara - ng_nt_lap, 25)}")
recon = edge_val + trivial_ihara + (golden_ihara_val - golden_lap_val) + (ng_nt_ihara - ng_nt_lap)
print(f"    SUM (= gap):                   {nstr(recon, 25)}")
print(f"    Actual gap:                    {nstr(gap, 25)}")

print(f"\n  Golden shift = {nstr(golden_ihara_val - golden_lap_val, 20)} (expected -3/2 = {nstr(mpf(-3)/2, 20)})")

# =============================================================================
# PART 13: EDGE TERM ANALYSIS
# =============================================================================
print(f"\n{SEP}")
print("PART 13: EDGE TERM AS THE CULPRIT")
print(SEP)

eigenvalue_sum_noedge = F_total(u0, include_edge=False)
needed_edge = target - eigenvalue_sum_noedge
print(f"\n  Eigenvalue sum (no edge): {nstr(eigenvalue_sum_noedge, 30)}")
print(f"  Needed edge term:         {nstr(needed_edge, 30)}")
print(f"  Actual edge term:         {nstr(edge_val, 30)}")
print(f"  Difference (= gap):       {nstr(needed_edge - edge_val, 30)}")

needed_coeff = needed_edge / (-2*u0**2/(1-u0**2))
print(f"\n  Needed (E-V) coefficient: {nstr(needed_coeff, 25)}")
print(f"  Actual E-V = {E-V}")

# Various edge exponents
print(f"\n  F with various edge coefficients:")
for label, coeff in [("E-V=10", mpf(10)), ("11", mpf(11)), ("9", mpf(9)),
                     ("c*", needed_coeff), ("E-V+gap/unit", needed_coeff),
                     ("V/2=10", mpf(10)), ("genus=11", mpf(11)),
                     ("(E-V)/2=5", mpf(5)), ("E/3=10", mpf(10)),
                     ("d(d-1)=6", mpf(6))]:
    edge_c = coeff * (-2*u0**2/(1-u0**2))
    F_c = eigenvalue_sum_noedge + edge_c
    g = F_c - target
    mark = "***" if fabs(g) < 1e-10 else ("**" if fabs(g) < 0.01 else ("*" if fabs(g) < 0.1 else ""))
    print(f"    {label:>20s}: F = {nstr(F_c, 20):>25s}  gap = {nstr(g, 15):>20s} {mark}")

# =============================================================================
# PART 14: EXHAUSTIVE SEARCH
# =============================================================================
print(f"\n{SEP}")
print("PART 14: EXHAUSTIVE SEARCH OVER NATURAL COMBINATIONS")
print(SEP)

print(f"\n  {'edge':>4s} {'triv':>4s} {'c_g':>5s} {'c_ng':>5s} {'F':>25s} {'gap':>20s} {'note':>15s}")
print(f"  {'----':>4s} {'----':>4s} {'-----':>5s} {'-----':>5s} {'-'*25:>25s} {'-'*20:>20s} {'-'*15:>15s}")

best_gap = mpf(100)
best_config = None

for use_edge in [True, False]:
    for use_trivial in [True, False]:
        for c_g in [mpf(0), mpf(1)/2, mpf(1), mpf(3)/2, mpf(2), mpf(5)/2, mpf(3)]:
            for c_ng in [mpf(0), mpf(1)/2, mpf(1), mpf(3)/2, mpf(2), mpf(5)/2, mpf(3)]:
                val = mpf(0)
                if use_edge:
                    val += edge_term(u0)
                for i, (lam, mult, name) in enumerate(eigenvalues):
                    if i == 0 and not use_trivial:
                        continue
                    c = c_g if i in golden_indices else c_ng
                    val += mult * f_contrib(u0, lam, c)

                g = val - target
                note = ""
                if fabs(g) < 1e-10:
                    note = "*** EXACT ***"
                elif fabs(g) < 0.01:
                    note = "** CLOSE **"
                elif fabs(g) < 0.1:
                    note = "* near *"

                if fabs(g) < fabs(best_gap):
                    best_gap = g
                    best_config = (use_edge, use_trivial, c_g, c_ng)

                if note or fabs(g) < 0.2:
                    e_str = "Y" if use_edge else "N"
                    t_str = "Y" if use_trivial else "N"
                    print(f"  {e_str:>4s} {t_str:>4s} {nstr(c_g,3):>5s} {nstr(c_ng,3):>5s} "
                          f"{nstr(val,20):>25s} {nstr(g,15):>20s} {note:>15s}")

if best_config:
    print(f"\n  Best: edge={best_config[0]}, trivial={best_config[1]}, "
          f"c_g={nstr(best_config[2],5)}, c_ng={nstr(best_config[3],5)}, gap={nstr(best_gap, 25)}")

# =============================================================================
# PART 15: LAPLACIAN-TYPE FACTORS
# =============================================================================
print(f"\n{SEP}")
print("PART 15: LAPLACIAN FACTORS -- can we match Tr(L^-1)?")
print(SEP)

# For each non-trivial eigenvalue, find u* where f(u*,lam,2) = 1/(3-lam)
print(f"\n  For each lam, find u where f(u,lam,2) = 1/(3-lam) = 1/mu:")
for i in [2, 3, 5]:
    lam, mult, name = eigenvalues[i]
    mu = 3 - lam
    target_val = 1/mu

    def eq(u, lam=lam, tv=target_val):
        return f_contrib(u, lam, 2) - tv

    try:
        u_match = findroot(eq, mpf('0.4'))
        print(f"  {name:30s}: mu={nstr(mu,5)}  1/mu={nstr(target_val,15)}  "
              f"u*={nstr(u_match, 25)}  (u0={nstr(u0, 12)})")
    except Exception as e:
        print(f"  {name:30s}: could not find u* ({e})")

# =============================================================================
# PART 16: WEIGHT CORRECTION
# =============================================================================
print(f"\n{SEP}")
print("PART 16: WEIGHT CORRECTION -- find weights w_i")
print(SEP)

# Natural question: what if each eigenvalue class gets weight w_i
# such that sum w_i * mult_i * f(u0, lam_i, 2) = 137/15 (without edge or trivial)?

# We have 4 unknowns (w_golden+, w_chi5, w_chi4, w_chi4') and 1 equation.
# Too underdetermined. But if we impose w_golden = 1 (keep golden fixed at 3):
# Then 3 + w_1*5*f(1) + w_0*4*f(0) + w_m2*4*f(-2) = 137/15
# With edge included: -10 + lam3*w_3 + 3 + w_1*5*f(1) + w_0*4*f(0) + w_m2*4*f(-2) = 137/15

# If we demand ALL w_i = w (same weight for all non-golden):
# Then: -10 + w*(lam3 + 5*f(1) + 4*f(0) + 4*f(-2)) + 3 = 137/15
# w * nongolden_eigenvalue_sum = 137/15 + 10 - 3 = 242/15
w_uniform = target_ng_sum / nongolden_eigenvalue_sum
print(f"\n  Uniform non-golden weight (with edge and golden=3):")
print(f"  w = {nstr(w_uniform, 30)}")
print(f"  w - 1 = {nstr(w_uniform - 1, 25)}")

# =============================================================================
# PART 17: EXACT SYMBOLIC GAP
# =============================================================================
print(f"\n{SEP}")
print("PART 17: EXACT SYMBOLIC GAP IN Q(sqrt3)")
print(SEP)

# F = 4308/715 + 270*sqrt3/143
# gap = F - 137/15 = 4308/715 - 137/15 + 270*sqrt3/143
# LCD(715, 15) = 2145.  4308/715 = 12924/2145.  137/15 = 19591/2145.
# Rational part of gap = (12924 - 19591)/2145 = -6667/2145
# Irrational part = 270*sqrt3/143 = 4050*sqrt3/2145

gap_exact = (mpf(-6667) + 4050*sqrt3) / 2145
print(f"  gap = (-6667 + 4050*sqrt3)/2145 = {nstr(gap_exact, 30)}")
print(f"  gap from computation:             {nstr(gap, 30)}")
print(f"  Match: {nstr(fabs(gap_exact - gap), 15)}")

# Let me verify the Q(sqrt3) form of F by computing from exact eigenvalue pieces:
# Non-golden eigenvalue sum in Q(sqrt3):
# lam3: (7+3sqrt3)/2 = 7/2 + (3/2)sqrt3
# 5*lam1: 5(17-sqrt3)/22 = 85/22 - (5/22)sqrt3
# 4*lam0: 16/5
# 4*lamm2: 4(8+2sqrt3)/13 = 32/13 + (8/13)sqrt3
# Edge: -10

# Rational part: -10 + 7/2 + 85/22 + 16/5 + 32/13 + 3 (golden)
# LCD(2, 22, 5, 13) = 2*5*11*13 = 1430
# -10 = -14300/1430
# 7/2 = 5005/1430
# 85/22 = 5525/1430... wait 85/22 * 1430/1430: 1430/22 = 65. 85*65 = 5525.
# 16/5 = 4576/1430... 1430/5 = 286. 16*286 = 4576.
# 32/13 = 3520/1430... 1430/13 = 110. 32*110 = 3520.
# 3 = 4290/1430
# Sum = (-14300 + 5005 + 5525 + 4576 + 3520 + 4290)/1430
#      = (5005+5525+4576+3520+4290-14300)/1430
#      = (22916-14300)/1430 = 8616/1430

# Simplify: gcd(8616, 1430). 8616 = 6*1436 = 6*4*359. 1430 = 2*5*11*13.
# gcd = 2. So 4308/715.

rat_part = mpf(-10) + mpf(7)/2 + mpf(85)/22 + mpf(16)/5 + mpf(32)/13 + 3
print(f"\n  Rational part of F: {nstr(rat_part, 30)}")
print(f"  4308/715 = {nstr(mpf(4308)/715, 30)}")
print(f"  Match: {nstr(fabs(rat_part - mpf(4308)/715), 15)}")

# Irrational part (coefficient of sqrt3):
# 3/2 - 5/22 + 8/13 from eigenvalues. Golden contributes 0 (golden pair sums to 3 exactly, rational).
# LCD(2, 22, 13) = 2*11*13 = 286
# 3/2 = 429/286
# -5/22 = -65/286
# 8/13 = 176/286
# Sum = (429-65+176)/286 = 540/286 = 270/143

irr_coeff = mpf(3)/2 - mpf(5)/22 + mpf(8)/13
print(f"\n  sqrt3 coefficient: {nstr(irr_coeff, 30)}")
print(f"  270/143 = {nstr(mpf(270)/143, 30)}")
print(f"  Match: {nstr(fabs(irr_coeff - mpf(270)/143), 15)}")

print(f"\n  F = 4308/715 + (270/143)*sqrt3  [VERIFIED]")
print(f"\n  gap = (-6667 + 4050*sqrt3)/2145")
print(f"  6667 = 59 x 113")
print(f"  4050 = 2 x 3^4 x 5^2")
print(f"  2145 = 3 x 5 x 11 x 13")

# Key observation
print(f"\n  KEY: denominator 15 in 137/15 = d^2 + d + 3 = {d**2+d+3}")
print(f"  gap denominator 2145 = 15 x 143 = 15 x 11 x 13")
print(f"  F denominator 715 = 5 x 11 x 13 = 5 x 143")
print(f"  Common factor: 143 = 11 x 13")

# =============================================================================
# PART 18: THE STRIKING PATTERN -- edge + trivial
# =============================================================================
print(f"\n{SEP}")
print("PART 18: EDGE + TRIVIAL = ? (the graph-theoretic overhead)")
print(SEP)

overhead = edge_val + trivial_ihara
print(f"\n  Edge:     {nstr(edge_val, 25)}")
print(f"  Trivial:  {nstr(trivial_ihara, 25)}")
print(f"  Sum:      {nstr(overhead, 25)}")

# edge + trivial = -10 + (7+3sqrt3)/2 = -20/2 + (7+3sqrt3)/2 = (-13+3sqrt3)/2
overhead_exact = (-13 + 3*sqrt3)/2
print(f"  Exact: (-13+3sqrt3)/2 = {nstr(overhead_exact, 25)}")
print(f"  Numerical: {nstr(overhead_exact, 15)}")

# So F = (-13+3sqrt3)/2 + golden(=3) + (5*(17-sqrt3)/22 + 16/5 + (32+8sqrt3)/13)
# The "pure spectral" contribution (non-trivial, non-golden, no edge):
pure = 5*exact_lam1 + 4*exact_lam0 + 4*exact_lamm2
print(f"\n  Pure non-trivial non-golden: {nstr(pure, 25)}")
print(f"  = 5*(17-sqrt3)/22 + 16/5 + 4*(8+2sqrt3)/13")

# Compute: 85/22 - 5sqrt3/22 + 16/5 + 32/13 + 8sqrt3/13
pure_rat = mpf(85)/22 + mpf(16)/5 + mpf(32)/13
pure_irr = -mpf(5)/22 + mpf(8)/13
print(f"  Rational part: {nstr(pure_rat, 25)}")
print(f"  sqrt3 coeff: {nstr(pure_irr, 25)}")

# And the Laplacian equivalent (non-golden part): 5/2 + 4/3 + 4/5 = 139/30
print(f"  Laplacian non-golden: 5/2 + 4/3 + 4/5 = 139/30 = {nstr(mpf(139)/30, 20)}")
print(f"  Pure Ihara non-golden: {nstr(pure, 20)}")
print(f"  Ratio: {nstr(pure / (mpf(139)/30), 25)}")

# =============================================================================
# PART 19: THE CORRECTION EQUATION
# =============================================================================
print(f"\n{SEP}")
print("PART 19: THE CORRECTION EQUATION -- what transforms F to Tr(L^-1)?")
print(SEP)

# We want to find a map T such that T(f(u0,lam,2)) = 1/(3-lam) for each non-trivial lam.
# The non-golden eigenvalues:
# lam=1: f = (17-sqrt3)/22, target 1/(3-1) = 1/2
# lam=0: f = 4/5, target 1/(3-0) = 1/3
# lam=-2: f = (8+2sqrt3)/13, target 1/(3-(-2)) = 1/5

# Golden:
# f(sqrt5) + f(-sqrt5) summed with mult 3: gives 3 in Ihara, 9/2 in Laplacian

print("\n  Eigenvalue-by-eigenvalue: Ihara f vs Laplacian 1/mu:")
pairs = [
    ("lam=1",  exact_lam1,  mpf(1)/2,  5),
    ("lam=0",  exact_lam0,  mpf(1)/3,  4),
    ("lam=-2", exact_lamm2, mpf(1)/5,  4),
]

for name, f_val, lap_val, mult in pairs:
    ratio = f_val / lap_val
    diff = f_val - lap_val
    print(f"  {name:10s}: f={nstr(f_val,18):>22s}  1/mu={nstr(lap_val,10):>12s}  "
          f"ratio={nstr(ratio,15):>18s}  diff={nstr(diff,15):>18s}")
    print(f"  {'':10s}  f*mult={nstr(mult*f_val,18):>22s}  mult/mu={nstr(mult*lap_val,10):>12s}")

# Check if the ratio has a pattern
r1 = exact_lam1 / (mpf(1)/2)
r0 = exact_lam0 / (mpf(1)/3)
rm2 = exact_lamm2 / (mpf(1)/5)
print(f"\n  Ratios f/(1/mu):")
print(f"    lam=1:  {nstr(r1, 25)} = (17-sqrt3)/11")
print(f"    lam=0:  {nstr(r0, 25)} = 12/5")
print(f"    lam=-2: {nstr(r0, 25)} = 5(8+2sqrt3)/13")

# What about the INVERSE formula: 1/mu = f * correction?
# correction = 1/(mu*f) for each eigenvalue
corr1 = 1 / (2 * exact_lam1)
corr0 = 1 / (3 * exact_lam0)
corrm2 = 1 / (5 * exact_lamm2)
print(f"\n  Correction factors (1/mu) / f = 1/(mu*f):")
print(f"    lam=1:  {nstr(corr1, 25)} = 11/(17-sqrt3)")
print(f"    lam=0:  {nstr(corr0, 25)} = 5/12")
print(f"    lam=-2: {nstr(corrm2, 25)} = 13/(5(8+2sqrt3)))")

# Rationalize lam=1 correction: 11/(17-sqrt3) = 11(17+sqrt3)/(289-3) = 11(17+sqrt3)/286 = (17+sqrt3)/26
corr1_exact = (17+sqrt3)/26
print(f"\n  lam=1 correction rationalized: (17+sqrt3)/26 = {nstr(corr1_exact, 20)}")
print(f"  Check: {nstr(fabs(corr1 - corr1_exact), 15)}")

# lam=-2 correction: 13/(5(8+2sqrt3)) = 13(8-2sqrt3)/(5(64-12)) = 13(8-2sqrt3)/(5*52) = (8-2sqrt3)/20 = (4-sqrt3)/10
corrm2_exact = (4-sqrt3)/10
print(f"  lam=-2 correction rationalized: (4-sqrt3)/10 = {nstr(corrm2_exact, 20)}")
print(f"  Check: {nstr(fabs(corrm2 - corrm2_exact), 15)}")

print(f"\n  So the corrections are:")
print(f"    lam=1:  c = (17+sqrt3)/26  = (17+sqrt3)/(2*13)")
print(f"    lam=0:  c = 5/12           = 5/(4*3)")
print(f"    lam=-2: c = (4-sqrt3)/10   = (4-sqrt3)/(2*5)")

# Is there a universal formula c(lam)?
# lam=1: (17+sqrt3)/26. Note 17 = 15+2 = (d^2+d+3) + (d-1). And 26 = 2*13.
# lam=0: 5/12. Note 5 = d+2 = 5-lam. And 12 = 4*3 = 4*d.
# lam=-2: (4-sqrt3)/10. Note 4 = 2cu0*sqrt3 = 4. And 10 = 2*5 = 2*(d+2).

# Try: c(lam) = (5-lam*sqrt3 + sqrt3) / (2(5-lam*sqrt3))  ... = 1 + sqrt3/(2(5-lam*sqrt3))
# lam=1: (5-sqrt3+sqrt3)/(2(5-sqrt3)) = 5/(10-2sqrt3)... nope

# Let's try the UNIVERSAL formula. The Ihara contribution at u=1/sqrt3 is:
# f(lam) = (4-lam*sqrt3)/(5-lam*sqrt3)
# The Laplacian is: g(lam) = 1/(3-lam)
# Correction: g/f = (5-lam*sqrt3)/((3-lam)(4-lam*sqrt3))

# Let's verify:
print(f"\n  Universal correction g(lam)/f(lam) = (5-lam*sqrt3)/((3-lam)(4-lam*sqrt3)):")
for lam_val, name in [(1, "lam=1"), (0, "lam=0"), (-2, "lam=-2")]:
    lv = mpf(lam_val)
    corr_univ = (5 - lv*sqrt3) / ((3-lv)*(4-lv*sqrt3))
    f_val = (4 - lv*sqrt3) / (5 - lv*sqrt3)
    check = f_val * corr_univ
    print(f"    {name}: correction = {nstr(corr_univ, 20)}, f*corr = {nstr(check, 20)}, "
          f"1/(3-lam) = {nstr(1/(3-lv), 20)}")

# =============================================================================
# PART 20: THE GOLDEN CORRECTION AT LAPLACIAN LEVEL
# =============================================================================
print(f"\n{SEP}")
print("PART 20: GOLDEN SECTOR -- Ihara vs Laplacian")
print(SEP)

# Golden in Ihara: 3*f(sqrt5) + 3*f(-sqrt5) = 3 (proven)
# Golden in Laplacian: 3/(3-sqrt5) + 3/(3+sqrt5) = 9/2

# The individual pieces:
f_gp = f_contrib(u0, sqrt5, 2)
f_gm = f_contrib(u0, -sqrt5, 2)
l_gp = 1 / (3 - sqrt5)
l_gm = 1 / (3 + sqrt5)

print(f"\n  lam=+sqrt5:  Ihara f = {nstr(f_gp, 20)}  Lap 1/mu = {nstr(l_gp, 20)}")
print(f"  lam=-sqrt5:  Ihara f = {nstr(f_gm, 20)}  Lap 1/mu = {nstr(l_gm, 20)}")
print(f"  Sums (x3): Ihara = {nstr(3*(f_gp+f_gm), 20)}  Lap = {nstr(3*(l_gp+l_gm), 20)}")

# The golden correction for each piece:
corr_gp = l_gp / f_gp
corr_gm = l_gm / f_gm
print(f"\n  Golden corrections (Lap/Ihara):")
print(f"    lam=+sqrt5: {nstr(corr_gp, 25)}")
print(f"    lam=-sqrt5: {nstr(corr_gm, 25)}")

# Using the universal formula: (5-lam*sqrt3)/((3-lam)*(4-lam*sqrt3))
for lam_val, name in [(sqrt5, "+sqrt5"), (-sqrt5, "-sqrt5")]:
    corr_u = (5 - lam_val*sqrt3) / ((3-lam_val)*(4-lam_val*sqrt3))
    print(f"    Universal for {name}: {nstr(corr_u, 25)}")

# Product of golden corrections:
print(f"\n  Product of golden corrections: {nstr(corr_gp * corr_gm, 25)}")
print(f"  Sum of golden corrections: {nstr(corr_gp + corr_gm, 25)}")

# =============================================================================
# PART 21: THE MASTER FORMULA -- does correction = (5-lam sqrt3)/((3-lam)(4-lam sqrt3))?
# =============================================================================
print(f"\n{SEP}")
print("PART 21: MASTER CORRECTION FORMULA")
print(SEP)

# The correction that maps Ihara f to Laplacian 1/mu is:
# C(lam) = (5 - lam sqrt3) / ((3-lam)(4 - lam sqrt3))
#
# This is specific to u = 1/sqrt3 and c = d-1 = 2.
# In terms of u: denominator of f is (5-lam sqrt3)/(3), numerator is (4-lam sqrt3)/(3).
# So f = (4-lam sqrt3)/(5-lam sqrt3).
#
# C(lam) = (5-lam sqrt3) / ((3-lam)(4-lam sqrt3))
# = [1/f(lam)] * 1/(3-lam)    ... wait: 1/f = (5-lam sqrt3)/(4-lam sqrt3)
# C = [1/f] / (3-lam) ... no: C = (5-lam sqrt3)/((3-lam)(4-lam sqrt3))
# = [(5-lam sqrt3)/(4-lam sqrt3)] * [1/(3-lam)]
# = (1/f) * (1/mu)    ... no that's C = (1/f)*(1/mu), but we want C = (1/mu)/f = 1/(mu*f).
# Actually C = g/f where g = 1/mu = 1/(3-lam). And f = (4-lam sqrt3)/(5-lam sqrt3).
# So C = (5-lam sqrt3) / ((3-lam)(4-lam sqrt3)).

# Simplify: let x = lam sqrt3.
# C = (5-x)/((3-lam)(4-x)) = (5-x)/(4-x) * 1/(3-lam)

# The factor (5-x)/(4-x) = 1 + 1/(4-x). Interesting.

# In the PRODUCT over all non-zero eigenvalues:
# prod C(lam_i)^{mult_i} transforms the product of Ihara factors to the product of Laplacian factors.
# This product should relate Z_Ihara to det(L)^{-1} (or similar).

# Compute the product of corrections
all_corrections = []
for i, (lam, mult, name) in enumerate(eigenvalues):
    if i == 0:  # skip trivial (mu=0)
        continue
    mu = 3 - lam
    f_val = f_contrib(u0, lam, 2)
    if fabs(f_val) > 1e-50:
        c_val = (1/mu) / f_val
        all_corrections.append((name, mult, c_val))
        print(f"  {name:30s}: C = {nstr(c_val, 25)}, C^mult = {nstr(c_val**mult, 25)}")

prod_corr = mpf(1)
for name, mult, c_val in all_corrections:
    prod_corr *= c_val**mult

print(f"\n  Product of all corrections: {nstr(prod_corr, 30)}")

# Sum of weighted corrections: sum mult * C * f = sum mult/mu = Tr(L^{-1})
# This is trivially true by construction. The real question is:
# is there a SINGLE correction C (or a simple function of lam) that works?

# C(lam) = (5-lam sqrt3)/((3-lam)(4-lam sqrt3)) at u=1/sqrt3
# Generalize to arbitrary u: C(lam,u) = (1-lam u+cu^2) / ((3-lam) * u(2cu-lam))
# = den_Ihara / ((3-lam) * num_Ihara)

# So f(u,lam,c) * C(lam,u) = 1/(3-lam) = 1/mu. Yes, trivially.

# The NON-trivial question: what STRUCTURAL change to the Ihara formula produces Tr(L^{-1})?

# The Ihara formula at u=1/sqrt(d) relates to walks on the graph.
# Tr(L^{-1}) relates to the Green's function (effective resistance).
# The bridge between them involves the RESOLVENT:

# (I - lam/d * A)^{-1} evaluated at lam = 1 gives... something.

# Let's compute the RESOLVENT sum: sum mult * 1/(1 - lam*z) at z = u0 = 1/sqrt3:
print(f"\n  Resolvent sum (1 - lam*u0)^{{-1}}:")
res_sum = mpf(0)
for i, (lam, mult, name) in enumerate(eigenvalues):
    r = 1 / (1 - lam * u0)
    res_sum += mult * r
    print(f"    {name:30s}: 1/(1-lam*u0) = {nstr(r, 20)}, x{mult} = {nstr(mult*r, 20)}")
print(f"  Total: {nstr(res_sum, 25)}")
print(f"  20/(1-u0^2) = {nstr(20/(1-u0**2), 25)}")  # trace of resolvent

# What about sum mult * lam/(1-lam*u0)?
lam_res = mpf(0)
for i, (lam, mult, name) in enumerate(eigenvalues):
    lam_res += mult * lam / (1 - lam * u0)
print(f"\n  Sum mult * lam/(1-lam u0) = {nstr(lam_res, 25)}")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print(f"\n{SEP}")
print("FINAL SUMMARY")
print(SEP)

F_val = F_total(u0)
F_c1_val = F_total(u0, c_golden=1, c_nongolden=1)

print(f"""
  At u = 1/sqrt3 (critical line, s=1/2):

  F(c=2, standard Ihara)    = {nstr(F_val, 25)}
  F(c=1, Artin coefficient) = {nstr(F_c1_val, 25)}
  Target 137/15             = {nstr(target, 25)}

  Gaps:
    Standard Ihara:     {nstr(F_val - target, 25)}
    Artin coefficient:  {nstr(F_c1_val - target, 25)}

  Golden pair at c=2:  {nstr(golden_total, 25)} (EXACTLY 3)
  Non-golden+edge:     {nstr(nongolden_with_edge, 25)}

  GAP STRUCTURE in Q(sqrt3):
    gap = (-6667 + 4050 sqrt3) / 2145
        = (-59 x 113 + 2 x 3^4 x 5^2 x sqrt3) / (3 x 5 x 11 x 13)

  KEY OBSERVATION:
    15 = d^2 + d + 3 appears in the target 137/15 AND controls the gap denominator.
    The gap numerator has 4050 = 2(d(d^2+d+3))^2 / d^2 = 2 x 45^2 where 45 = d x 15.

  UNIVERSAL CORRECTION:
    C(lam) = (5 - lam sqrt3) / ((3-lam)(4 - lam sqrt3))
    maps each f(u0, lam, c=2) to 1/(3-lam) = 1/mu (Laplacian).
    This is the resolvent correction: C = den_Ihara / (mu x num_Ihara).
    It encodes the transition from walk-counting (Ihara) to spectral (Laplacian).

  GAP DECOMPOSITION:
    Edge term:           {nstr(edge_val, 20)}
    Trivial (lam=3):     {nstr(trivial_ihara, 20)}
    Golden shift:        {nstr(golden_shift, 20)} (= -d/2 = -3/2)
    ng-nt Ihara-Lap:     {nstr(ng_nt_ihara - ng_nt_lap, 20)}
    TOTAL = gap:         {nstr(recon, 20)}
""")

print(SEP)
print("DONE")
print(SEP)
