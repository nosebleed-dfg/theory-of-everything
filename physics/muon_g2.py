"""
MUON_G2 — muon g-2 anomaly from bicone framework; Delta_a_mu = alpha^2 * phi^(-4) * (m_mu/m_p)^2 / (4*pi^2)
nos3bl33d

Two cones (QED + hadronic) plus the phi^(-4) growth band middle ring.
"""

from mpmath import mp, mpf, pi, phi, sqrt, power, log, nstr, fabs

mp.dps = 30

# ============================================================
# Physical constants (CODATA 2018 / PDG 2024)
# ============================================================
alpha_inv = mpf("137.035999084")          # inverse fine structure constant
alpha = 1 / alpha_inv                      # alpha = 7.2973525693e-3

m_e   = mpf("0.51099895000")              # electron mass (MeV)
m_mu  = mpf("105.6583755")               # muon mass (MeV)
m_tau = mpf("1776.86")                    # tau mass (MeV)
m_p   = mpf("938.27208816")              # proton mass (MeV)
m_n   = mpf("939.56542052")              # neutron mass (MeV)

# Experimental and SM values (× 10^-11)
a_mu_exp    = mpf("116592061")            # measured: 116592061(41) × 10^-11
a_mu_exp_err = mpf("41")
a_mu_SM     = mpf("116591810")            # SM prediction: 116591810(43) × 10^-11
a_mu_SM_err = mpf("43")

# The anomaly (× 10^-11)
delta_a_mu_exp = a_mu_exp - a_mu_SM       # = 251
delta_a_mu_err = sqrt(a_mu_exp_err**2 + a_mu_SM_err**2)  # combined error

# Golden ratio
PHI = phi  # mpmath's phi = (1+sqrt(5))/2

print("=" * 72)
print("MUON g-2 ANOMALY FROM THE THREE-FACE BICONE FRAMEWORK")
print("=" * 72)

print(f"\n{'PHYSICAL CONSTANTS':=^72}")
print(f"  alpha           = {nstr(alpha, 15)}")
print(f"  alpha^(-1)      = {nstr(alpha_inv, 15)}")
print(f"  phi             = {nstr(PHI, 30)}")
print(f"  m_e             = {nstr(m_e, 12)} MeV")
print(f"  m_mu            = {nstr(m_mu, 12)} MeV")
print(f"  m_tau           = {nstr(m_tau, 8)} MeV")
print(f"  m_p             = {nstr(m_p, 12)} MeV")
print(f"  m_n             = {nstr(m_n, 12)} MeV")

print(f"\n{'EXPERIMENTAL ANOMALY':=^72}")
print(f"  a_mu(exp)       = {a_mu_exp} x 10^-11")
print(f"  a_mu(SM)        = {a_mu_SM} x 10^-11")
print(f"  Delta_a_mu      = {delta_a_mu_exp} +/- {nstr(delta_a_mu_err, 5)} x 10^-11")

# ============================================================
# Framework quantities
# ============================================================
Delta = PHI**(-4)                          # growth band width
mass_ratio = m_mu / m_p                    # muon-to-proton mass ratio

print(f"\n{'FRAMEWORK QUANTITIES':=^72}")
print(f"  phi^(-4)        = {nstr(Delta, 30)}")
print(f"  phi^(-4) exact  = (2 - phi)^2 = {nstr((2 - PHI)**2, 30)}")
print(f"  m_mu/m_p        = {nstr(mass_ratio, 20)}")
print(f"  (m_mu/m_p)^2    = {nstr(mass_ratio**2, 20)}")
print(f"  alpha^2         = {nstr(alpha**2, 20)}")
print(f"  4*pi^2          = {nstr(4*pi**2, 20)}")

# ============================================================
# APPROACH 1: The master formula
# Delta_a_mu = alpha^2 * phi^(-4) * (m_mu/m_p)^2 / (4*pi^2)
# ============================================================
print(f"\n{'APPROACH 1: MASTER FORMULA':=^72}")
print("  Formula: Delta_a_mu = alpha^2 * phi^(-4) * (m_mu/m_p)^2 / (4*pi^2)")
print()

result_1 = alpha**2 * Delta * mass_ratio**2 / (4 * pi**2)

# Convert to units of 10^-11
result_1_units = result_1 * mpf("1e11")

print(f"  Step-by-step:")
print(f"    alpha^2                    = {nstr(alpha**2, 20)}")
print(f"    alpha^2 * phi^(-4)         = {nstr(alpha**2 * Delta, 20)}")
print(f"    ... * (m_mu/m_p)^2         = {nstr(alpha**2 * Delta * mass_ratio**2, 20)}")
print(f"    ... / (4*pi^2)             = {nstr(result_1, 20)}")
print()
print(f"  RESULT: {nstr(result_1, 20)}")
print(f"  RESULT: {nstr(result_1_units, 15)} x 10^-11")
print(f"  TARGET: {delta_a_mu_exp} +/- {nstr(delta_a_mu_err, 5)} x 10^-11")
print()

deviation_1 = result_1_units - delta_a_mu_exp
sigma_1 = deviation_1 / delta_a_mu_err
pct_1 = (deviation_1 / delta_a_mu_exp) * 100

print(f"  Deviation:  {nstr(deviation_1, 10)} x 10^-11")
print(f"  Deviation:  {nstr(sigma_1, 6)} sigma")
print(f"  Deviation:  {nstr(pct_1, 6)}%")
print(f"  WITHIN ERROR BARS: {'YES' if fabs(sigma_1) < 1 else 'NO'} ({nstr(fabs(sigma_1), 4)} sigma)")

# ============================================================
# APPROACH 2: Ring geometry interpretation
# A_eff = Delta * (m_mu/m_p)^2 / 4
# Delta_a_mu = (alpha/pi)^2 * A_eff
# ============================================================
print(f"\n{'APPROACH 2: RING GEOMETRY':=^72}")
print("  Formula: Delta_a_mu = (alpha/pi)^2 * Delta * (m_mu/m_p)^2 / 4")
print()

alpha_over_pi = alpha / pi
A_eff = Delta * mass_ratio**2 / 4
result_2 = alpha_over_pi**2 * A_eff
result_2_units = result_2 * mpf("1e11")

print(f"  alpha/pi        = {nstr(alpha_over_pi, 20)}")
print(f"  (alpha/pi)^2    = {nstr(alpha_over_pi**2, 20)}")
print(f"  A_eff           = {nstr(A_eff, 20)}")
print(f"  RESULT:           {nstr(result_2_units, 15)} x 10^-11")
print(f"  (Same as Approach 1: {nstr(result_1 - result_2, 5)})")

# ============================================================
# APPROACH 3: Using framework proton mass
# m_p/m_e = 6*pi^5 + phi^(-7) + 3*phi^(-21)
# ============================================================
print(f"\n{'APPROACH 3: FRAMEWORK PROTON MASS':=^72}")

mp_me_framework = 6*pi**5 + PHI**(-7) + 3*PHI**(-21)
mp_framework = m_e * mp_me_framework
mass_ratio_fw = m_mu / mp_framework

print(f"  m_p/m_e (framework) = {nstr(mp_me_framework, 20)}")
print(f"  m_p/m_e (measured)  = {nstr(m_p/m_e, 20)}")
print(f"  m_p (framework)     = {nstr(mp_framework, 15)} MeV")
print(f"  m_p (measured)      = {nstr(m_p, 15)} MeV")
fw_mp_err_ppm = fabs(mp_framework - m_p) / m_p * 1e6
print(f"  m_p error           = {nstr(fw_mp_err_ppm, 6)} ppm")

result_3 = alpha**2 * Delta * mass_ratio_fw**2 / (4 * pi**2)
result_3_units = result_3 * mpf("1e11")
deviation_3 = result_3_units - delta_a_mu_exp
sigma_3 = deviation_3 / delta_a_mu_err

print(f"\n  With framework m_p:")
print(f"    m_mu/m_p(fw)    = {nstr(mass_ratio_fw, 20)}")
print(f"    RESULT:           {nstr(result_3_units, 15)} x 10^-11")
print(f"    TARGET:           {delta_a_mu_exp} +/- {nstr(delta_a_mu_err, 5)} x 10^-11")
print(f"    Deviation:        {nstr(deviation_3, 10)} x 10^-11")
print(f"    Deviation:        {nstr(sigma_3, 6)} sigma")

# ============================================================
# APPROACH 4: Corrections to master formula
# ============================================================
print(f"\n{'APPROACH 4: CORRECTIONS':=^72}")

corrections = [
    ("1 + alpha/pi",            1 + alpha/pi),
    ("1 + alpha/(2*pi)",        1 + alpha/(2*pi)),
    ("1 + phi^(-8)",            1 + PHI**(-8)),
    ("1 + phi^(-10)",           1 + PHI**(-10)),
    ("1 + alpha*phi/pi",        1 + alpha*PHI/pi),
    ("1 + 1/(4*pi^2)",          1 + 1/(4*pi**2)),
    ("1 + (m_e/m_mu)^2",        1 + (m_e/m_mu)**2),
    ("1 + alpha^2",             1 + alpha**2),
]

print(f"  Base result: {nstr(result_1_units, 15)} x 10^-11\n")
print(f"  {'Correction':<25} {'Factor':>20} {'Result (x10^-11)':>20} {'Sigma':>10}")
print(f"  {'-'*25} {'-'*20} {'-'*20} {'-'*10}")

for name, factor in corrections:
    r = result_1 * factor * mpf("1e11")
    d = r - delta_a_mu_exp
    s = d / delta_a_mu_err
    mark = " <--" if fabs(s) < fabs(sigma_1) else ""
    print(f"  {name:<25} {nstr(factor, 15):>20} {nstr(r, 12):>20} {nstr(s, 5):>10}{mark}")

# ============================================================
# APPROACH 5: Dimensional analysis survey
# ============================================================
print(f"\n{'APPROACH 5: DIMENSIONAL ANALYSIS SURVEY':=^72}")

target = mpf("2.51e-9")  # 251 x 10^-11

formulas = [
    ("alpha^2 * Delta",
     alpha**2 * Delta),
    ("alpha^3 * Delta",
     alpha**3 * Delta),
    ("alpha^2 * Delta / (4*pi^2)",
     alpha**2 * Delta / (4*pi**2)),
    ("(alpha/pi)^2 * Delta",
     (alpha/pi)**2 * Delta),
    ("alpha^2*Delta*(m_mu/m_tau)^2/(2*pi)",
     alpha**2 * Delta * (m_mu/m_tau)**2 / (2*pi)),
    ("alpha^3/pi * Delta",
     alpha**3 / pi * Delta),
    ("alpha^2*Delta*(m_mu/m_p)^2/(4*pi^2)",
     alpha**2 * Delta * mass_ratio**2 / (4*pi**2)),
    ("alpha^2*Delta*(m_mu/m_p)^2/(2*pi)^2",
     alpha**2 * Delta * mass_ratio**2 / (2*pi)**2),
]

print(f"  Target: {nstr(target, 10)}\n")
print(f"  {'Formula':<45} {'Value':>15} {'Ratio to target':>18}")
print(f"  {'-'*45} {'-'*15} {'-'*18}")

for name, val in formulas:
    ratio = val / target
    mark = " *** MATCH ***" if fabs(ratio - 1) < mpf("0.05") else ""
    print(f"  {name:<45} {nstr(val, 8):>15} {nstr(ratio, 8):>18}{mark}")

# ============================================================
# APPROACH 6: HVP fraction approach
# ============================================================
print(f"\n{'APPROACH 6: HADRONIC VACUUM POLARIZATION':=^72}")

a_HVP = mpf("6845")  # HVP contribution in units of 10^-11

hvp_formulas = [
    ("a_HVP * phi^(-4)",              a_HVP * Delta),
    ("a_HVP * phi^(-4)^2",           a_HVP * Delta**2),
    ("a_HVP * alpha * phi^(-4)/pi",  a_HVP * alpha * Delta / pi),
    ("a_HVP * alpha * phi^(-4)",     a_HVP * alpha * Delta),
    ("a_HVP * phi^(-4) * alpha^2",   a_HVP * Delta * alpha**2),
    ("a_HVP * (m_mu/m_p)^2 * phi^(-4)", a_HVP * mass_ratio**2 * Delta),
]

print(f"  a_mu(HVP) = {a_HVP} x 10^-11\n")
print(f"  {'Formula':<45} {'Value (x10^-11)':>18} {'vs 251':>10}")
print(f"  {'-'*45} {'-'*18} {'-'*10}")

for name, val in hvp_formulas:
    ratio = val / delta_a_mu_exp
    print(f"  {name:<45} {nstr(val, 8):>18} {nstr(ratio, 5):>10}")

# ============================================================
# Muon-electron mass ratio: searching framework expressions
# ============================================================
print(f"\n{'MUON/ELECTRON MASS RATIO SEARCH':=^72}")

m_mu_over_m_e = m_mu / m_e
print(f"  m_mu/m_e (measured)  = {nstr(m_mu_over_m_e, 20)}")
print(f"  Target: 206.7682830...\n")

# Dodecahedral constants
V  = 20   # vertices
E  = 30   # edges
F  = 12   # faces
d  = 3    # dimension
p  = 5    # pentagons
b0 = 11   # Betti number / first zero
chi = 2   # Euler characteristic

attempts = [
    ("V*b0 - F - chi",                        V*b0 - F - chi),
    ("V*b0 - F - 1",                          V*b0 - F - 1),
    ("V*(b0-1) + p",                          V*(b0-1) + p),
    ("207 - sin^2(theta_W)",                  207 - mpf("0.23122")),
    ("V*b0 - F - chi + phi^(-2)",             V*b0 - F - chi + PHI**(-2)),
    ("V*b0 - F - chi + 1 - 1/d",             V*b0 - F - chi + 1 - mpf(1)/d),
    ("V*b0 - F - chi + phi^(-1) - 1/d",      V*b0 - F - chi + PHI**(-1) - mpf(1)/d),
    ("V*b0 - F - phi^(-1) - 1",              V*b0 - F - PHI**(-1) - 1),
    ("V*b0 - F + phi^(-4) - d",              V*b0 - F + PHI**(-4) - d),
    ("d*(E+p*d) + phi^(-1) + phi^(-4)",       d*(E+p*d) + PHI**(-1) + PHI**(-4)),
    ("4*p*b0 - V + d*chi + phi^(-1)",         4*p*b0 - V + d*chi + PHI**(-1)),
    ("6*E + V + d*chi + phi^(-1)",            6*E + V + d*chi + PHI**(-1)),
    ("E*p + V*d + phi^(-1) - 1/d",           E*p + V*d + PHI**(-1) - mpf(1)/d),
]

print(f"  {'Expression':<45} {'Value':>20} {'Error':>12}")
print(f"  {'-'*45} {'-'*20} {'-'*12}")

for name, val in attempts:
    err = val - m_mu_over_m_e
    err_pct = fabs(err / m_mu_over_m_e) * 100
    mark = " <-- CLOSE" if err_pct < mpf("0.1") else ""
    print(f"  {name:<45} {nstr(val, 12):>20} {nstr(err, 8):>12}{mark}")

# 207 - sin^2(theta_W) is intriguing, let's examine
sin2_thetaW = mpf("0.23122")
m_mu_me_approx = 207 - sin2_thetaW
err_sw = fabs(m_mu_me_approx - m_mu_over_m_e) / m_mu_over_m_e * 100
print(f"\n  ** 207 - sin^2(theta_W) = {nstr(m_mu_me_approx, 12)} (error: {nstr(err_sw, 4)}%)")
print(f"     207 = V*b0 - F - 1 = 20*11 - 12 - 1")
print(f"     sin^2(theta_W) = 0.23122")
print(f"     Residual: {nstr(m_mu_me_approx - m_mu_over_m_e, 10)}")

# ============================================================
# Neutron-proton mass difference
# ============================================================
print(f"\n{'NEUTRON-PROTON MASS DIFFERENCE':=^72}")

delta_np_meas = m_n - m_p
print(f"  m_n - m_p (measured) = {nstr(delta_np_meas, 15)} MeV")

np_formulas = [
    ("alpha * m_p * phi^(-4)",             alpha * m_p * Delta),
    ("alpha * m_p / phi",                  alpha * m_p / PHI),
    ("alpha * m_p * d/(2*pi)",             alpha * m_p * d / (2*pi)),
    ("3*alpha*m_p/(2*pi)",                 3*alpha*m_p/(2*pi)),
    ("m_e * chi + alpha*m_p",              m_e * chi + alpha*m_p),
    ("m_e*(d - 1/phi)",                    m_e*(d - 1/PHI)),
    ("(m_e + alpha*m_p)*phi^(-1)",         (m_e + alpha*m_p)*PHI**(-1)),
    ("alpha*m_p + m_e*(phi^(-4)-1)",       alpha*m_p + m_e*(Delta - 1)),
]

print(f"\n  {'Formula':<40} {'Value (MeV)':>15} {'Error (MeV)':>15} {'Error %':>10}")
print(f"  {'-'*40} {'-'*15} {'-'*15} {'-'*10}")

for name, val in np_formulas:
    err = val - delta_np_meas
    err_pct = fabs(err / delta_np_meas) * 100
    mark = " <--" if err_pct < 5 else ""
    print(f"  {name:<40} {nstr(val, 8):>15} {nstr(err, 6):>15} {nstr(err_pct, 5):>10}{mark}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'FINAL SUMMARY':=^72}")
print()
print("  THE MASTER FORMULA:")
print("  Delta_a_mu = alpha^2 * phi^(-4) * (m_mu/m_p)^2 / (4*pi^2)")
print()
print(f"  Framework prediction: {nstr(result_1_units, 15)} x 10^-11")
print(f"  Measured anomaly:     {delta_a_mu_exp} +/- {nstr(delta_a_mu_err, 5)} x 10^-11")
print(f"  Agreement:            {nstr(fabs(pct_1), 4)}% ({nstr(fabs(sigma_1), 4)} sigma)")
print()
print("  INTERPRETATION:")
print("  - alpha^2:       QED vertex squared (two-loop)")
print("  - phi^(-4):      growth band width = spectral gap of bicone equator")
print("  - (m_mu/m_p)^2:  muon-proton mass ratio = hadronic scale coupling")
print("  - 1/(4*pi^2):    ring circumference normalization")
print()
print("  The Standard Model computes the two cones (QED + hadronic).")
print("  The anomaly = the middle ring at sigma = 1/2 that SM misses.")
print("  Ring area = phi^(-4) * (m_mu/m_p)^2 / 4")
print(f"  Effective ring area = {nstr(A_eff, 15)}")
print()

# Pure framework expression check
print("  PURE FRAMEWORK CHECK:")
print(f"    alpha^2 * phi^(-4) * (m_mu/m_p)^2 / (4*pi^2)")
print(f"    = {nstr(result_1, 25)}")
print(f"    = {nstr(result_1_units, 20)} x 10^-11")
print()

# Using framework proton mass
print("  WITH FRAMEWORK PROTON MASS:")
print(f"    m_p/m_e = 6*pi^5 + phi^(-7) + 3*phi^(-21)")
print(f"    Result: {nstr(result_3_units, 15)} x 10^-11")
print(f"    Sigma:  {nstr(fabs(sigma_3), 4)}")
print()

# The key insight
print("  KEY INSIGHT:")
print("  The formula involves exactly the three icosahedral operations:")
print("    ADDITION:        alpha = 1/137 (coupling strength)")
print("    MULTIPLICATION:  (m_mu/m_p)^2 (mass scale ratio)")
print("    GROWTH:          phi^(-4) = Delta (spectral gap)")
print("    ROTATION:        1/(4*pi^2) (circumference normalization)")
print()
print(f"  Prediction accuracy: {nstr(fabs(pct_1), 4)}% of central value")
print(f"  Well within {nstr(delta_a_mu_err, 4)}-unit experimental+theory uncertainty")
print("=" * 72)
