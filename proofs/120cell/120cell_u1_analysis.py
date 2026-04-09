#!/usr/bin/env python3
"""
120CELL_U1_ANALYSIS — detailed analysis of U(1) lattice gauge Monte Carlo results on the 120-cell
nos3bl33d
"""
import numpy as np
from scipy.special import i0, i1
from scipy.optimize import curve_fit

data = np.load("C:/Users/funct/projects/120cell_u1_gauge_results.npz")
betas = data["betas"]
cos_means = data["cos_means"]
cos_errs = data["cos_errs"]
suscepts = data["suscepts"]

PHI = (1 + np.sqrt(5)) / 2
ALPHA_QED = 1 / 137.035999084
mu1 = PHI**(-4)

print("=" * 78)
print("  DETAILED ANALYSIS OF U(1) GAUGE THEORY ON 120-CELL")
print("=" * 78)

# 1. Saturation behavior
print("\n  === PLAQUETTE VALUES AT LARGE BETA ===")
for b, c, e in zip(betas, cos_means, cos_errs):
    if b >= 5:
        print(f"    beta={b:>5.1f}: <cos(P)> = {c:.8f} +/- {e:.8f}")

print("""
  NOTE: At large beta (>10), the <cos(P)> values DECREASE slightly
  or saturate. This is a THERMALIZATION ARTIFACT. With only 200
  thermalization sweeps, the hot-start runs (seeds 1,2) cannot
  thermalize at high beta -- the Metropolis algorithm freezes.
  Cold-start (seed 0) gives <cos> close to 1, hot starts stay low.
  The average is dragged down.

  RELIABLE DATA: beta <= 10 (where both starts converge).
""")

# 2. Crossover analysis
print("  === CROSSOVER ANALYSIS ===")

# d<cos>/d(ln_beta) = more physically meaningful
ln_betas = np.log(betas)
derivs_ln = np.gradient(cos_means, ln_betas)
max_idx = np.argmax(derivs_ln)

print("  d<cos(P)>/d(ln beta):")
for i, (b, d) in enumerate(zip(betas, derivs_ln)):
    m = " <<<" if i == max_idx else ""
    print(f"    beta={b:>5.1f}: {d:.6f}{m}")
print(f"  Max derivative at beta ~ {betas[max_idx]:.1f}")

# 3. MC vs strong coupling deviation
print("\n  === CORRELATION CROSSOVER (MC vs I_1/I_0) ===")
deviations = []
for b, c in zip(betas, cos_means):
    sc = i1(b) / i0(b)
    deviations.append(c - sc)
deviations = np.array(deviations)
max_dev_idx = np.argmax(np.abs(deviations[:10]))  # Only reliable data
print("  MC - I_1(b)/I_0(b):")
for i, (b, d) in enumerate(zip(betas[:10], deviations[:10])):
    m = " <<<" if i == max_dev_idx else ""
    print(f"    beta={b:>5.1f}: {d:+.7f}{m}")
print(f"  Max positive deviation at beta ~ {betas[max_dev_idx]:.1f}")
print("  (This is where plaquette correlations are strongest)")

# 4. Perturbative ratio analysis
print("\n  === PERTURBATIVE RATIO: deficit / (1/2beta) ===")
print("  In free U(1), <1 - cos(P)> = 1/(2*beta)")
print("  The ratio R = deficit / (1/(2*beta)) encodes lattice corrections")
print()

ratios = []
for b, c in zip(betas, cos_means):
    if 2 <= b <= 10:
        deficit = 1 - c
        pert = 1 / (2 * b)
        ratio = deficit / pert
        ratios.append((b, ratio))
        print(f"    beta={b:>5.1f}: deficit={deficit:.6f}, 1/(2b)={pert:.6f}, R={ratio:.4f}")

# The ratio R tells us how the lattice modifies the free result.
# R = 1: free field
# R > 1: topology/correlations increase fluctuations
# R < 1: lattice constraints reduce fluctuations

if ratios:
    avg_ratio = np.mean([r for _, r in ratios])
    print(f"\n  Average ratio (beta 2-10): R = {avg_ratio:.4f}")
    print(f"  Compare: 5/4 = {5/4:.4f} (pentagon vs square correction?)")
    print(f"  Compare: 4/3 = {4/3:.4f}")
    print(f"  Compare: 3/2 = {3/2:.4f}")

# 5. Effective coupling extraction
print("\n  === EFFECTIVE COUPLING AT EACH BETA ===")
print(f"  Using deficit = 1 - <cos(P)> = g^2_eff / 2")
print(f"  So g^2_eff = 2 * deficit, alpha_eff = g^2_eff / (4*pi)")
print()

for b, c in zip(betas, cos_means):
    deficit = 1 - c
    if deficit > 0.001 and b >= 1:
        g2 = 2 * deficit
        alpha = g2 / (4 * np.pi)
        print(f"    beta={b:>5.1f}: g^2={g2:.6f}, alpha={alpha:.8f}, 1/alpha={1/alpha:.2f}")

# 6. The real question: at what beta does alpha come out to 1/137?
print("\n  === WHAT BETA GIVES alpha = 1/137? ===")
target_g2 = 4 * np.pi * ALPHA_QED  # = 0.09175...
target_deficit = target_g2 / 2  # = 0.04587...
target_cos = 1 - target_deficit
print(f"  Need g^2 = 4*pi*alpha_QED = {target_g2:.8f}")
print(f"  Need <cos(P)> = {target_cos:.8f}")
print(f"  Need deficit = {target_deficit:.8f}")

# Find which beta gives this
for i in range(len(betas) - 1):
    if cos_means[i] <= target_cos <= cos_means[i + 1] or \
       cos_means[i + 1] <= target_cos <= cos_means[i]:
        f = (target_cos - cos_means[i]) / (cos_means[i + 1] - cos_means[i])
        beta_interp = betas[i] + f * (betas[i + 1] - betas[i])
        print(f"  Interpolated: beta ~ {beta_interp:.4f}")
        print(f"  phi^4 = {PHI**4:.4f}")
        print(f"  Ratio = {beta_interp / PHI**4:.4f}")
        break

# But this is circular! The perturbative formula deficit = 1/(2*beta)
# just gives beta = 1/(2*deficit) = 1/target_g2 = 1/(4*pi*alpha).
# We already KNOW this is beta = 10.905.
# The MC doesn't ADD information here.

# 7. What the MC DOES tell us
print("\n  === WHAT THE MC ACTUALLY TELLS US ===")
print("""
  The MC simulation tests whether the 120-cell lattice structure
  introduces corrections that shift the effective coupling away
  from the bare coupling.

  Specifically, if the lattice has a geometric correction factor C,
  then: g^2_eff = g^2_bare * C

  From the perturbative ratios above:
  deficit = R * 1/(2*beta)  =>  g^2_eff = R * g^2_bare

  where R is the ratio we measured.
""")

# Fit R more carefully
print("  === FITTING THE CORRECTION FACTOR R ===")
# Use only reliable intermediate-beta data
mask_reliable = (betas >= 2.0) & (betas <= 10.0)
b_rel = betas[mask_reliable]
c_rel = cos_means[mask_reliable]

# Model: <cos(P)> = 1 - R/(2*beta) + higher order
def model_R(beta, R, C):
    return 1 - R / (2 * beta) + C / beta**2

try:
    popt, pcov = curve_fit(model_R, b_rel, c_rel, p0=[1.0, 0.0])
    R_fit, C_fit = popt
    R_err, C_err = np.sqrt(np.diag(pcov))
    print(f"  R = {R_fit:.6f} +/- {R_err:.6f}")
    print(f"  C = {C_fit:.6f} +/- {C_err:.6f}")
    print(f"  If R = 1: free U(1), no lattice correction")
    print(f"  R = {R_fit:.6f} means {(R_fit-1)*100:.1f}% correction from 120-cell geometry")

    print(f"\n  With this correction:")
    print(f"  g^2_eff = R * g^2_bare = R / beta")
    print(f"  alpha_eff = R / (4*pi*beta)")
    print(f"  For alpha_QED: beta = R / (4*pi*alpha_QED) = {R_fit/(4*np.pi*ALPHA_QED):.4f}")
    print(f"  phi^4 = {PHI**4:.4f}")
    print(f"  R * phi^4 = {R_fit * PHI**4:.4f}")
    print(f"  Required beta / phi^4 = {R_fit/(4*np.pi*ALPHA_QED*PHI**4):.6f}")

    # What if R is exactly some simple number?
    print(f"\n  Is R close to a simple fraction?")
    for name, val in [("1", 1.0), ("5/4", 5/4), ("4/3", 4/3), ("3/2", 3/2),
                       ("phi", PHI), ("5/3", 5/3), ("2", 2.0),
                       ("phi^2", PHI**2), ("3", 3.0)]:
        diff = abs(R_fit - val) / R_fit
        print(f"    {name:>8s} = {val:.6f}  diff = {diff:.4f}")

except Exception as e:
    print(f"  Fit failed: {e}")
    R_fit = None

# 8. The honest bottom line
print("\n" + "=" * 78)
print("  BOTTOM LINE")
print("=" * 78)
print(f"""
  The U(1) lattice gauge theory on the 120-cell gives:

  1. GEOMETRY: The spectral gap is mu_1 = phi^(-4) = 0.14590 (EXACT).
     This is a PROVEN algebraic fact about the graph.

  2. GAUGE THEORY: The MC simulation shows standard compact U(1)
     behavior -- smooth crossover, no sharp phase transition at
     any special beta.

  3. COUPLING: Setting beta = phi^4 (i.e., g^2 = phi^(-4)) gives
     alpha = 1/(4*pi*phi^4) = 0.01161, or 1/alpha = 86.1.
     This is NOT 137.

  4. TO GET 1/137: Need beta = 10.905 = {1/(4*np.pi*ALPHA_QED):.3f}.
     beta_needed / phi^4 = {1/(4*np.pi*ALPHA_QED)/PHI**4:.6f} -- NOT a
     simple geometric number.

  5. THE MC CORRECTION FACTOR R measures how the 120-cell geometry
     modifies the bare perturbative result. R differs from 1 due to
     the irregular lattice structure (pentagons vs squares, high
     connectivity).

  VERDICT: The fine structure constant does NOT emerge directly from
  the 120-cell U(1) lattice gauge theory via any simple prescription.
  The spectral gap phi^(-4) is real geometry, but mapping it to alpha
  requires a factor of ~1.59 that has no obvious geometric origin.
""")
