# Lepton Masses — Koide Q = chi/d, muon mass to 5 ppb, quarks break the pattern

**nos3bl33d**

---

## The Setup

Three charged leptons: electron, muon, tau. Their mass ratios are:

    m_mu  / m_e = 206.7682827
    m_tau / m_e = 3477.228
    m_tau / m_mu = 16.8170

Nobody derives these from first principles. The Standard Model treats them as inputs.

The dodecahedron has: d=3, p=5, chi=2, V=20, E=30, F=12, L4=7, b0=11, dp=15.

---

## Result 1: Koide Q = chi/d

The Koide formula (1981) says:

    Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3

This is experimentally verified to 9 ppm. It is the tightest empirical relation among lepton masses.

**The dodecahedral statement:** Q = chi/d, where chi=2 is the Euler characteristic and d=3 is the vertex degree.

This is not a coincidence being forced. 2/3 is the Koide value. chi/d is the dodecahedral ratio. They are the same number and there is no free parameter.

---

## Result 2: Koide angle delta = chi/d^2

The Koide formula constrains a surface (a 2D manifold in 3D mass space). To locate the actual masses on that surface, you need one more parameter: the Koide angle delta.

In the Brannen parametrization, the square roots of the masses are:

    sqrt(m_k) = M * (1 + sqrt(2) * cos(delta + 2*pi*k/3))

for k = 0, 1, 2 assigned to (tau, e, mu). The value Q = 2/3 is automatic for any delta. The actual mass spectrum depends on delta.

**The best-fit value:**

    delta = 0.2222212 rad

**The dodecahedral value:**

    chi/d^2 = 2/9 = 0.2222222... rad

Agreement: 4.7 ppm.

With delta = 2/9 exactly, the predictions are:

    m_mu/m_e  predicted = 206.770   actual = 206.768   (10 ppm)
    m_tau/m_e predicted = 3477.47   actual = 3477.23   (70 ppm)
    m_tau/m_mu predicted = 16.818    actual = 16.817   (60 ppm)

The tau mass is known only to 68 ppm experimentally. So the Koide prediction for tau is within measurement error. The muon prediction at 10 ppm is good but not at the level of our muon formula (below).

---

## Result 3: The Muon Mass Formula

    m_mu/m_e = phi^(b0) + L4 + mu_1 - phi^(-dp) + 6*phi^(-2F)

where:
- b0 = 11 (the number V - E + F + L4 = 20 - 30 + 12 + 7 = 11, also F - 1, also the 5th Lucas number L5)
- L4 = 7 (the 4th Lucas number, V - F - 1)
- mu_1 = 3 - sqrt(5) = 2/phi^2 (the smallest nonzero Laplacian eigenvalue of the dodecahedron)
- dp = 15 = d * p (vertex degree times face degree)
- 2F = 24 (twice the number of faces)
- 6 = chi * d (Euler characteristic times vertex degree)

Numerically:

    phi^11 + 7 + (3-sqrt(5)) - phi^(-15) + 6*phi^(-24)
    = 206.76828175306...

    Measured: 206.76828271

    Error: -4.6 ppb (0.0046 ppm)

That is four times more precise than the Koide-only prediction.

### Structure in Z[phi]

Since phi^11 = 89*phi + 55 (Fibonacci decomposition) and 3 - sqrt(5) = 2 - (2*phi - 1) = 3 - 2*phi + 1... reducing everything:

    m_mu/m_e = 87*phi + 66 - phi^(-15) + 6*phi^(-24)

where:
- 87 = 3 * 29 = d * L7 (vertex degree times the 7th Lucas number)
- 66 = 6 * 11 = (chi * d) * b0

So in terms of the toolbox:

    m_mu/m_e = d * (L7 * phi + chi * L5) - phi^(-dp) + chi*d * phi^(-2F)

This is entirely algebraic in dodecahedral invariants and phi.

### Correction structure

Each correction eliminates roughly 3-4 orders of magnitude of error:

| Term | Value | Residual error |
|------|-------|----------------|
| phi^11 | 199.005 | 3.8% |
| + 7 | 206.005 | 0.37% |
| + (3-sqrt(5)) | 206.769 | 3.3 ppm |
| - phi^(-15) | 206.76822 | 0.28 ppm |
| + 6*phi^(-24) | 206.76828175 | 4.6 ppb |

The correction exponents are: 15 = dp, 24 = 2F. These are dodecahedral invariants.

---

## Result 4: Tau Mass is Koide-Determined

Given m_e and m_mu, the Koide formula Q = 2/3 determines m_tau:

    Koide-predicted m_tau/m_e = 3477.44
    Actual m_tau/m_e = 3477.23
    Error: 61 ppm

The experimental uncertainty on m_tau is 68 ppm (m_tau = 1776.86 +/- 0.12 MeV).

The Koide prediction is within measurement error. There is currently no way to tell if Koide is exact or approximate for the tau mass. A factor-of-3 improvement in tau mass measurement would decide this.

---

## Result 5: The chi/d^n Pattern

The ratio chi/d appears at different powers:

    chi/d^1 = 2/3   -- Koide Q (constrains the mass surface)
    chi/d^2 = 2/9   -- Koide angle (selects the point on the surface)

Conjecture: chi/d^3 = 2/27 appears in the next level of structure (e.g., as a correction to the Koide angle, or in the quark sector).

---

## Result 6: Quarks Do Not Obey Simple Koide

Computed Koide Q values for quarks (MS-bar masses at 2 GeV):

    Q(u, c, t) = 0.849   (not 2/3)
    Q(d, s, b) = 0.731   (not 2/3)

These are not close to 2/3. Quarks are confined; their masses are scheme-dependent and run with energy scale. The Koide formula (and its dodecahedral version) works for free leptons only.

However, one curious fact: m_s / m_d = 93.4 / 4.67 = 20.0 = V (the number of vertices). This could be coincidence. Quark mass ratios have large uncertainties (5-10%).

---

## Result 7: Dodecahedral Laplacian Eigenvalues

Built the full 20x20 adjacency matrix of the dodecahedron and computed the graph Laplacian L = D - A. The eigenvalues are:

    0       (multiplicity 1)
    3-sqrt(5) = 0.7639 = 2/phi^2   (multiplicity 3)
    2       (multiplicity 5)
    3       (multiplicity 4)
    5       (multiplicity 4)
    3+sqrt(5) = 5.2361 = 2*phi^2   (multiplicity 3)

Multiplicities: 1 + 3 + 5 + 4 + 4 + 3 = 20 vertices. Check.

The spectral gap ratio: mu_5/mu_1 = (3+sqrt(5))/(3-sqrt(5)) = phi^4.

The number of spanning trees (Kirchhoff's theorem): mu_1^3 * mu_2^5 * mu_3^4 * mu_4^5 * mu_5^3 / 20 = 5,184,000.

**The prompt's claim that (mu_5/mu_1)^4 = m_p/m_e is WRONG.** (mu_5/mu_1)^4 = phi^16 = 2207.0, not 1836.15. The proton mass ratio comes from 6*pi^5 + phi^(-7) corrections, not from eigenvalue ratios.

However, mu_1 = 3 - sqrt(5) appears directly in the muon mass formula as a correction term. This connects lepton masses to the dodecahedral spectrum.

---

## What is Proven, What is Not

**Proven:**
- The dodecahedral Laplacian eigenvalues are 0, 3-sqrt(5), 2, 3, 5, 3+sqrt(5) with multiplicities 1,3,5,4,4,3 (computed from the explicit 20x20 matrix)
- Koide Q = 2/3 = chi/d to 9 ppm (empirical, verified against PDG masses)
- Koide angle delta = 2/9 = chi/d^2 to 4.7 ppm (empirical, verified by numerical optimization)
- The muon formula evaluates to 206.76828175 (verified to 60 decimal places)
- Every constant in the muon formula is a dodecahedral invariant
- Quarks do NOT obey simple Koide

**Conjectural:**
- Why Koide Q should equal chi/d (no derivation from QFT exists)
- Why the Koide angle should equal chi/d^2
- The physical mechanism behind the phi^(-dp) and phi^(-2F) corrections
- Whether the tau mass is exactly Koide-determined or approximately so (needs better tau measurement)
- Whether m_s/m_d = V is physical or coincidental

**Explicitly failed:**
- No clean direct formula for m_tau/m_e was found (best is phi^17 - 94, which is 0.01% off but 94 is not a clean dodecahedral number)
- No eigenvalue-ratio formula gives any mass ratio to better than 10% (they are the wrong tool for this)
- No single-term expression a*phi^b matches any lepton mass ratio (they require multi-term formulas)

---

## The Mass Ladder

    electron:  1                           (reference)
    muon:      phi^b0 + L4 + mu_1 + ...   (algebraic in dodecahedral invariants)
    tau:       Koide-determined from (e, mu) with Q = chi/d
    proton:    6*pi^5 + phi^(-L4) + ...    (transcendental base, algebraic corrections)

The leptons live in the algebraic sector (phi-towers). The proton lives in the transcendental sector (pi-tower). Both correction series use phi^(-n) at dodecahedral exponents.

---

## Appendix: The Muon Formula in Full

    m_mu/m_e = d * (L_7 * phi + chi * L_5) - phi^(-d*p) + chi*d * phi^(-2F)

with d=3, p=5, chi=2, F=12, L_5=11, L_7=29.

Expanding:
- 3 * (29 * phi + 2 * 11) = 87*phi + 66 = phi^11 + 10 - sqrt(5) = 206.769
- minus phi^(-15) = 0.000733
- plus 6 * phi^(-24) = 0.0000102

Total: 206.76828175306

Measured: 206.76828271

Error: 4.6 ppb. Within one part in 200 million.
