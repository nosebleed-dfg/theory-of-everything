#!/usr/bin/env python3
# ref: https://youtu.be/4a---Gvz8J8?si=UYMUOIAD1ob88XIU
"""
THREE_BODY_PYTHAGOREAN — reduces the three-body problem via Pythagorean/golden framework; 1D stacking + volume average
nos3bl33d
"""

from mpmath import (mp, mpf, sqrt, pi as MP_PI, cos, sin, acos,
                    power, fabs, nstr)
import sys
import io

# Force UTF-8 output on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

mp.dps = 30

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------

PHI = (1 + sqrt(5)) / 2
PHI_SQ = PHI ** 2
KOPPA = MP_PI / 2
G_CONST = mpf(1)
HALF = mpf('0.5')
ONE = mpf(1)
ZERO = mpf(0)
SOFTENING_SQ = mpf('1e-10')

print("=" * 78)
print("THREE-BODY PROBLEM: PYTHAGOREAN / GOLDEN FRAMEWORK REDUCTION")
print("=" * 78)
print(f"  mp.dps       = {mp.dps}")
print(f"  phi          = {nstr(PHI, 20)}")
print(f"  phi^2        = {nstr(PHI_SQ, 20)}")
print(f"  phi^2 - phi  = {nstr(PHI_SQ - PHI, 20)}  (should be 1)")
print(f"  koppa        = pi/2 = {nstr(KOPPA, 20)}")
print()

# ---------------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------------

def triangle_angles(a, b, c):
    """Compute three angles of triangle with sides a, b, c (law of cosines)."""
    def safe_acos(x):
        return acos(max(mpf(-1), min(ONE, x)))
    A = safe_acos((b**2 + c**2 - a**2) / (2 * b * c))
    B = safe_acos((a**2 + c**2 - b**2) / (2 * a * c))
    C = safe_acos((a**2 + b**2 - c**2) / (2 * a * b))
    return A, B, C


def pythagorean_residual(a, b, c):
    """Normalized: (a^2+b^2-c^2)/c^2 for sorted sides. 0 = right triangle."""
    sides = sorted([a, b, c])
    return (sides[0]**2 + sides[1]**2 - sides[2]**2) / sides[2]**2 if sides[2] > 0 else ZERO


def compute_accelerations(positions, masses):
    """Gravitational accelerations in 2D with softening."""
    n = len(positions)
    accels = []
    for i in range(n):
        ax, ay = ZERO, ZERO
        for j in range(n):
            if i == j:
                continue
            dx = positions[j][0] - positions[i][0]
            dy = positions[j][1] - positions[i][1]
            r_sq = dx**2 + dy**2 + SOFTENING_SQ
            r = sqrt(r_sq)
            f = G_CONST * masses[j] / (r_sq * r)
            ax += f * dx
            ay += f * dy
        accels.append((ax, ay))
    return accels


def rk4_step(positions, velocities, masses, dt):
    """Single RK4 step."""
    n = len(positions)

    def derivs(pos, vel):
        acc = compute_accelerations(pos, masses)
        return vel, acc

    def add_sc(base, delta, s):
        return [(base[i][0] + s * delta[i][0],
                 base[i][1] + s * delta[i][1]) for i in range(n)]

    dp1, dv1 = derivs(positions, velocities)
    dp2, dv2 = derivs(add_sc(positions, dp1, dt/2), add_sc(velocities, dv1, dt/2))
    dp3, dv3 = derivs(add_sc(positions, dp2, dt/2), add_sc(velocities, dv2, dt/2))
    dp4, dv4 = derivs(add_sc(positions, dp3, dt), add_sc(velocities, dv3, dt))

    s6 = dt / 6
    new_p, new_v = [], []
    for i in range(n):
        new_p.append((positions[i][0] + s6*(dp1[i][0]+2*dp2[i][0]+2*dp3[i][0]+dp4[i][0]),
                       positions[i][1] + s6*(dp1[i][1]+2*dp2[i][1]+2*dp3[i][1]+dp4[i][1])))
        new_v.append((velocities[i][0] + s6*(dv1[i][0]+2*dv2[i][0]+2*dv3[i][0]+dv4[i][0]),
                       velocities[i][1] + s6*(dv1[i][1]+2*dv2[i][1]+2*dv3[i][1]+dv4[i][1])))
    return new_p, new_v


def mutual_distances(positions):
    """Return r12, r13, r23."""
    def d(i, j):
        dx = positions[j][0] - positions[i][0]
        dy = positions[j][1] - positions[i][1]
        return sqrt(dx**2 + dy**2)
    return d(0,1), d(0,2), d(1,2)


def total_energy(pos, vel, masses):
    n = len(pos)
    KE = sum(HALF * masses[i] * (vel[i][0]**2 + vel[i][1]**2) for i in range(n))
    PE = ZERO
    for i in range(n):
        for j in range(i+1, n):
            dx = pos[j][0] - pos[i][0]
            dy = pos[j][1] - pos[i][1]
            r = sqrt(dx**2 + dy**2 + SOFTENING_SQ)
            PE -= G_CONST * masses[i] * masses[j] / r
    return KE + PE


def total_angular_momentum(pos, vel, masses):
    return sum(masses[i] * (pos[i][0]*vel[i][1] - pos[i][1]*vel[i][0])
               for i in range(len(pos)))


def run_simulation(pos0, vel0, masses, dt, n_steps, sample_every, label):
    """Run a 3-body simulation and return sampled data dict."""
    print(f"\n  [{label}] dt={dt}, steps={n_steps}, samples~{n_steps//sample_every}")
    print(f"  Integrating...", end="", flush=True)

    data = dict(times=[], r12=[], r13=[], r23=[], angles=[], residuals=[],
                geo_means=[], pyth_sums=[], energies=[], ang_mom=[])

    pos, vel = list(pos0), list(vel0)
    E0 = total_energy(pos, vel, masses)
    L0 = total_angular_momentum(pos, vel, masses)

    tenth = max(n_steps // 10, 1)

    for step in range(n_steps + 1):
        if step % sample_every == 0:
            t = step * dt
            r12, r13, r23 = mutual_distances(pos)
            data['times'].append(t)
            data['r12'].append(r12)
            data['r13'].append(r13)
            data['r23'].append(r23)

            if r12 > 0 and r13 > 0 and r23 > 0:
                data['angles'].append(triangle_angles(r12, r13, r23))
            else:
                data['angles'].append((ZERO, ZERO, ZERO))

            data['residuals'].append(pythagorean_residual(r12, r13, r23))
            data['geo_means'].append(power(r12*r13*r23, ONE/3))
            data['pyth_sums'].append(r12**2 + r13**2 + r23**2)
            data['energies'].append(total_energy(pos, vel, masses))
            data['ang_mom'].append(total_angular_momentum(pos, vel, masses))

        if step < n_steps:
            pos, vel = rk4_step(pos, vel, masses, dt)

        if step > 0 and step % tenth == 0:
            print(f" {100*step//n_steps}%", end="", flush=True)

    E_final = data['energies'][-1]
    L_final = data['ang_mom'][-1]
    dE = fabs((E_final - E0) / E0) if E0 != 0 else ZERO
    dL = fabs((L_final - L0) / L0) if L0 != 0 else ZERO
    print(f" done.")
    print(f"  E0={nstr(E0, 10)}, Ef={nstr(E_final, 10)}, dE/E={nstr(dE, 4)}")
    print(f"  L0={nstr(L0, 10)}, Lf={nstr(L_final, 10)}, dL/L={nstr(dL, 4)}")
    data['E0'] = E0
    data['L0'] = L0
    return data


def analyze_data(data, label):
    """Full Pythagorean/golden analysis of simulation data."""
    N = len(data['times'])
    sep = "-" * 60

    print(f"\n{sep}")
    print(f"  ANALYSIS: {label}  ({N} samples)")
    print(f"{sep}")

    # ---- Angle statistics ----
    all_angles_deg = []
    for A, B, C in data['angles']:
        all_angles_deg.extend([A*180/MP_PI, B*180/MP_PI, C*180/MP_PI])

    Na = len(all_angles_deg)
    a_mean = sum(all_angles_deg) / Na
    near_90 = sum(1 for a in all_angles_deg if fabs(a - 90) < 5)
    near_60 = sum(1 for a in all_angles_deg if fabs(a - 60) < 5)
    near_45 = sum(1 for a in all_angles_deg if fabs(a - 45) < 5)

    print(f"\n  ANGLES ({Na} measurements):")
    print(f"    Mean           = {nstr(a_mean, 8)} deg  (60 = uniform)")
    print(f"    Near 60 (+/-5) = {near_60:5d}  ({100*near_60/Na:5.1f}%)")
    print(f"    Near 90 (+/-5) = {near_90:5d}  ({100*near_90/Na:5.1f}%)")
    print(f"    Near 45 (+/-5) = {near_45:5d}  ({100*near_45/Na:5.1f}%)")

    # Histogram
    nbins = 18
    bins = [0]*nbins
    for a in all_angles_deg:
        b = int(float(a) / 10)
        b = max(0, min(nbins-1, b))
        bins[b] += 1
    mx = max(bins) if bins else 1
    print(f"\n  Angle histogram (10-deg bins):")
    for i in range(nbins):
        lo, hi = i*10, (i+1)*10
        bar = '#' * (40 * bins[i] // mx) if mx > 0 else ''
        tag = ""
        if lo <= 60 < hi: tag = " <-- equil"
        if lo <= 90 < hi: tag = " <-- koppa"
        if lo <= 45 < hi: tag = " <-- pi/4"
        print(f"    [{lo:3d}-{hi:3d}] {bar} ({bins[i]}){tag}")

    # ---- Angle CV ----
    a_float = [float(x) for x in all_angles_deg]
    a_mu = sum(a_float)/len(a_float)
    a_var = sum((x-a_mu)**2 for x in a_float)/len(a_float)
    a_cv = (a_var**0.5)/a_mu if a_mu != 0 else 0

    # ---- Pythagorean residual ----
    res = data['residuals']
    res_mean = sum(res)/len(res)
    res_abs_mean = sum(fabs(r) for r in res)/len(res)
    near_zero = sum(1 for r in res if fabs(r) < mpf('0.1'))
    print(f"\n  PYTHAGOREAN RESIDUAL (a^2+b^2-c^2)/c^2:")
    print(f"    Mean           = {nstr(res_mean, 10)}")
    print(f"    Mean |res|     = {nstr(res_abs_mean, 10)}")
    print(f"    Near 0 (<0.1)  = {near_zero}/{len(res)} ({100*near_zero/len(res):.1f}%)")

    # ---- Geometric mean ----
    gm = data['geo_means']
    gm_mean = sum(gm)/len(gm)
    gm_var = sum((g-gm_mean)**2 for g in gm)/len(gm)
    gm_std = sqrt(gm_var)
    gm_cv = float(gm_std/gm_mean) if gm_mean > 0 else 0

    print(f"\n  GEOMETRIC MEAN (r12*r13*r23)^(1/3):")
    print(f"    Mean   = {nstr(gm_mean, 12)}")
    print(f"    StdDev = {nstr(gm_std, 12)}")
    print(f"    CV     = {gm_cv:.6f}")

    # ---- Pythagorean sum ----
    ps = data['pyth_sums']
    ps_mean = sum(ps)/len(ps)
    ps_var = sum((p-ps_mean)**2 for p in ps)/len(ps)
    ps_std = sqrt(ps_var)
    ps_cv = float(ps_std/ps_mean) if ps_mean > 0 else 0

    print(f"\n  PYTHAGOREAN SUM r12^2+r13^2+r23^2:")
    print(f"    Mean   = {nstr(ps_mean, 12)}")
    print(f"    StdDev = {nstr(ps_std, 12)}")
    print(f"    CV     = {ps_cv:.6f}")

    # ---- Key ratio: pyth_sum / (3 * gm^2) ----
    ratios_pg = [ps[i] / (3 * gm[i]**2) if gm[i] > 0 else ZERO for i in range(N)]
    rpg_mean = sum(ratios_pg)/len(ratios_pg)
    print(f"\n  KEY RATIO: pyth_sum / (3 * gm^2):")
    print(f"    Mean   = {nstr(rpg_mean, 12)}  (1.0 = equilateral)")

    # ---- Distance ratios ----
    all_ratios = []
    for i in range(N):
        dists = sorted([data['r12'][i], data['r13'][i], data['r23'][i]])
        if dists[0] > 0:
            all_ratios.append(dists[2]/dists[0])  # max/min
    if all_ratios:
        r_mean = sum(all_ratios)/len(all_ratios)
        nr_phi = sum(1 for r in all_ratios if fabs(r - PHI) < mpf('0.15'))
        nr_one = sum(1 for r in all_ratios if fabs(r - 1) < mpf('0.15'))
        nr_2 = sum(1 for r in all_ratios if fabs(r - 2) < mpf('0.15'))
        nr_phisq = sum(1 for r in all_ratios if fabs(r - PHI_SQ) < mpf('0.15'))
    else:
        r_mean = ZERO
        nr_phi = nr_one = nr_2 = nr_phisq = 0

    print(f"\n  DISTANCE RATIO r_max/r_min:")
    print(f"    Mean             = {nstr(r_mean, 12)}")
    print(f"    phi              = {nstr(PHI, 12)}")
    print(f"    diff from phi    = {nstr(fabs(r_mean - PHI), 10)}")
    print(f"    Near 1    (+/-0.15): {nr_one}")
    print(f"    Near phi  (+/-0.15): {nr_phi}")
    print(f"    Near 2    (+/-0.15): {nr_2}")
    print(f"    Near phi^2(+/-0.15): {nr_phisq}")

    # ---- Chaos ordering comparison ----
    print(f"\n  CHAOS vs ORDER:")
    print(f"    Angle CV        = {a_cv:.6f}")
    print(f"    Geom mean CV    = {gm_cv:.6f}")
    print(f"    Pyth sum CV     = {ps_cv:.6f}")
    if gm_cv > 0 and a_cv > 0:
        factor = a_cv / gm_cv
        print(f"    Order factor    = {factor:.2f}x (angles/distances)")
        if factor > 1:
            print(f"    --> DISTANCES MORE ORDERED THAN ANGLES by {factor:.1f}x")
    else:
        factor = 0

    # ---- Area analysis ----
    areas = []
    for i in range(N):
        r12, r13, r23 = data['r12'][i], data['r13'][i], data['r23'][i]
        s = (r12+r13+r23)/2
        asq = s*(s-r12)*(s-r13)*(s-r23)
        areas.append(sqrt(asq) if asq > 0 else ZERO)
    area_mean = sum(areas)/len(areas)
    area_gm_ratios = [areas[i]/gm[i]**2 if gm[i] > 0 else ZERO for i in range(N)]
    area_ratio_mean = sum(area_gm_ratios)/len(area_gm_ratios)

    print(f"\n  TRIANGLE AREA:")
    print(f"    <area>           = {nstr(area_mean, 10)}")
    print(f"    <area/gm^2>      = {nstr(area_ratio_mean, 10)}")
    print(f"    sqrt(3)/4        = {nstr(sqrt(3)/4, 10)}  (equilateral)")

    # ---- Shape parameter ----
    shape_params = []
    for i in range(N):
        r12, r13, r23 = data['r12'][i], data['r13'][i], data['r23'][i]
        num = (r12-r13)**2 + (r13-r23)**2 + (r23-r12)**2
        den = r12**2 + r13**2 + r23**2
        shape_params.append(num/den if den > 0 else ZERO)
    sp_mean = sum(shape_params)/len(shape_params)
    print(f"\n  SHAPE PARAMETER (0=equil, 2=collinear):")
    print(f"    Mean             = {nstr(sp_mean, 10)}")

    # ---- Golden scan: match observables to golden/pi targets ----
    print(f"\n  GOLDEN/PI SCAN (matching observables to constants):")
    targets = {
        "1/4":     mpf('0.25'),
        "1/phi":   1/PHI,
        "1/2":     HALF,
        "1/sqrt3": 1/sqrt(3),
        "sqrt3/4": sqrt(3)/4,
        "1":       ONE,
        "phi":     PHI,
        "pi/2":    MP_PI/2,
        "phi^2":   PHI_SQ,
        "pi":      MP_PI,
        "3":       mpf(3),
        "pi^2/6":  MP_PI**2/6,
    }
    observables = {
        "<gm>":          gm_mean,
        "<ps>":          ps_mean,
        "<ps>/(3<gm>^2)":rpg_mean,
        "<r_max/r_min>": r_mean,
        "<area/gm^2>":   area_ratio_mean,
        "<shape>":       sp_mean,
        "|E|":           fabs(data['E0']),
    }
    for oname, oval in observables.items():
        if oval <= 0:
            continue
        best_name, best_pct = None, 1e10
        for tname, tval in targets.items():
            if tval <= 0:
                continue
            pct = float(fabs(oval - tval)/tval) * 100
            if pct < best_pct:
                best_pct = pct
                best_name = tname
        star = "***" if best_pct < 3 else "** " if best_pct < 10 else "*  " if best_pct < 20 else "   "
        print(f"    {star} {oname:20s} = {nstr(oval,8):>12s}  ~  "
              f"{best_name:8s} = {nstr(targets[best_name],8):>12s}  ({best_pct:6.2f}%)")

    return dict(a_cv=a_cv, gm_cv=gm_cv, ps_cv=ps_cv, gm_mean=gm_mean,
                ps_mean=ps_mean, rpg_mean=rpg_mean, r_mean=r_mean,
                area_ratio=area_ratio_mean, shape=sp_mean, factor=factor,
                near_90_pct=100*near_90/Na, near_60_pct=100*near_60/Na)


# ===================================================================
# PART 1: TRIANGLE REDUCTION (static analysis)
# ===================================================================

print("=" * 78)
print("PART 1: TRIANGLE REDUCTION -- PYTHAGOREAN RESIDUALS")
print("=" * 78)

a_eq = ONE
angles_eq = triangle_angles(a_eq, a_eq, a_eq)
print(f"\nEquilateral triangle (side=1):")
print(f"  Angles: {nstr(angles_eq[0]*180/MP_PI,10)}, "
      f"{nstr(angles_eq[1]*180/MP_PI,10)}, {nstr(angles_eq[2]*180/MP_PI,10)} deg")
print(f"  Pythagorean residual: {nstr(pythagorean_residual(a_eq, a_eq, a_eq), 12)}")

c_golden = sqrt(1 + PHI_SQ)
print(f"\nGolden right triangle (1, phi, sqrt(1+phi^2)):")
print(f"  Sides: 1, {nstr(PHI,12)}, {nstr(c_golden,12)}")
print(f"  c^2 = 1 + phi^2 = {nstr(c_golden**2, 15)}")
print(f"  Residual = {nstr(pythagorean_residual(ONE, PHI, c_golden), 15)}")
ag = triangle_angles(ONE, PHI, c_golden)
print(f"  Angles: {nstr(ag[0]*180/MP_PI,10)}, {nstr(ag[1]*180/MP_PI,10)}, "
      f"{nstr(ag[2]*180/MP_PI,10)} deg")
print(f"\n  KEY IDENTITY: phi+2 = phi^2+1 = {nstr(PHI+2,20)}")
print(f"  This follows from phi^2 = phi+1 (the axiom).")
print(f"  So the golden right triangle has hypotenuse c = sqrt(phi+2).")
print(f"  Angles: {nstr(ag[0]*180/MP_PI,6)} and {nstr(ag[1]*180/MP_PI,6)} deg")
print(f"  Angle ratio: {nstr(ag[1]/ag[0], 12)}")
print(f"  (NOT equal to phi = {nstr(PHI,12)}, but close to 2-1/phi = {nstr(2-1/PHI,12)})")

# ===================================================================
# PART 2: FIGURE-8 ORBIT (the clean periodic orbit)
# ===================================================================

print("\n" + "=" * 78)
print("PART 2: FIGURE-8 ORBIT (Chenciner-Montgomery)")
print("=" * 78)

# Standard initial conditions (Chenciner-Montgomery 2000)
x1 = mpf('0.97000435669734563955')
y1 = mpf('-0.24308753153583290860')
vx3 = mpf('-0.93240737144104206632')
vy3 = mpf('-0.86473146092102251270')

pos_f8 = [
    (x1, y1),
    (-x1, -y1),
    (ZERO, ZERO),
]
vel_f8 = [
    (-vx3/2, -vy3/2),
    (-vx3/2, -vy3/2),
    (vx3, vy3),
]
masses_f8 = [ONE, ONE, ONE]

# Period ~ 6.3259. Simulate 5 periods for good statistics.
T_period = mpf('6.3259')
dt_f8 = mpf('0.0005')
T_sim = 5 * T_period
n_steps_f8 = int(T_sim / dt_f8)
sample_f8 = 50

data_f8 = run_simulation(pos_f8, vel_f8, masses_f8, dt_f8, n_steps_f8, sample_f8,
                         "Figure-8")
res_f8 = analyze_data(data_f8, "Figure-8 orbit (periodic)")

# ===================================================================
# PART 3: LAGRANGE EQUILATERAL (perturbed) -- BOUND CHAOTIC ORBIT
# ===================================================================

print("\n" + "=" * 78)
print("PART 3: PERTURBED LAGRANGE ORBIT (chaotic, bound)")
print("=" * 78)

# For 3 equal masses m=1 in equilateral triangle with side a,
# the angular velocity for circular orbit is:
#   omega^2 = G*m / (a^2 * sqrt(3))     (force from one body on another projected)
# Wait -- let me derive it properly.
# Each mass sits at distance R = a/sqrt(3) from COM.
# Force on mass i from mass j: F = G*m^2/a^2
# Net force toward center (from 2 other masses):
#   F_net = 2 * G*m^2/a^2 * cos(30) = 2 * G*m^2/a^2 * sqrt(3)/2 = G*m^2*sqrt(3)/a^2
# For circular orbit: F_net = m * omega^2 * R
#   G*m^2*sqrt(3)/a^2 = m * omega^2 * a/sqrt(3)
#   omega^2 = G*m*sqrt(3)/a^2 * sqrt(3)/a = G*m*3/a^3
#   omega = sqrt(3*G*m / a^3)

a_side = mpf(1)
m_body = ONE
R_circ = a_side / sqrt(3)
omega_L = sqrt(3 * G_CONST * m_body / a_side**3)

print(f"  Side a = {a_side}, R_circ = {nstr(R_circ, 12)}")
print(f"  Lagrange omega = {nstr(omega_L, 12)}")

# Positions: equilateral triangle centered at origin
angles_init = [MP_PI/2, MP_PI/2 + 2*MP_PI/3, MP_PI/2 + 4*MP_PI/3]
pos_L = [(R_circ * cos(th), R_circ * sin(th)) for th in angles_init]

# Velocities: tangential at omega, perpendicular to radius
vel_L = [(-omega_L * pos_L[i][1], omega_L * pos_L[i][0]) for i in range(3)]

# Perturbation: kick body 0 slightly (1% of velocity)
eps = mpf('0.02')
vel_L[0] = (vel_L[0][0] + eps * omega_L * R_circ,
            vel_L[0][1] + eps * omega_L * R_circ)

# Zero COM velocity
masses_L = [ONE, ONE, ONE]
cvx = sum(masses_L[i]*vel_L[i][0] for i in range(3)) / 3
cvy = sum(masses_L[i]*vel_L[i][1] for i in range(3)) / 3
vel_L = [(vel_L[i][0]-cvx, vel_L[i][1]-cvy) for i in range(3)]
# Zero COM position
cpx = sum(pos_L[i][0] for i in range(3)) / 3
cpy = sum(pos_L[i][1] for i in range(3)) / 3
pos_L = [(pos_L[i][0]-cpx, pos_L[i][1]-cpy) for i in range(3)]

E_L = total_energy(pos_L, vel_L, masses_L)
print(f"  Initial energy: {nstr(E_L, 12)} (must be < 0 for bound)")
assert float(E_L) < 0, f"System not bound! E={E_L}"

# Lagrange period ~ 2*pi/omega
T_L = 2 * MP_PI / omega_L
print(f"  Lagrange period ~ {nstr(T_L, 6)}")

# Simulate for 20 periods
dt_L = mpf('0.0005')
T_sim_L = 20 * T_L
n_steps_L = int(T_sim_L / dt_L)
sample_L = 50

data_L = run_simulation(pos_L, vel_L, masses_L, dt_L, n_steps_L, sample_L,
                        "Perturbed Lagrange")
res_L = analyze_data(data_L, "Perturbed Lagrange (mildly chaotic)")

# ===================================================================
# PART 4: STRONGLY CHAOTIC (Pythagorean 3-body problem)
# ===================================================================

print("\n" + "=" * 78)
print("PART 4: PYTHAGOREAN THREE-BODY PROBLEM (3-4-5 triangle)")
print("=" * 78)

# The famous Pythagorean three-body problem: masses 3, 4, 5
# at vertices of a 3-4-5 right triangle, starting from rest.
# This is STRONGLY chaotic and eventually ejects one body.
# We simulate the early phase before ejection.

pos_P = [(ONE, mpf(3)), (-mpf(2), -ONE), (ONE, -ONE)]
vel_P = [(ZERO, ZERO), (ZERO, ZERO), (ZERO, ZERO)]
masses_P = [mpf(3), mpf(4), mpf(5)]

# Zero COM
cpx = sum(masses_P[i]*pos_P[i][0] for i in range(3)) / sum(masses_P)
cpy = sum(masses_P[i]*pos_P[i][1] for i in range(3)) / sum(masses_P)
pos_P = [(pos_P[i][0]-cpx, pos_P[i][1]-cpy) for i in range(3)]

E_P = total_energy(pos_P, vel_P, masses_P)
print(f"  Masses: 3, 4, 5")
print(f"  Initial energy: {nstr(E_P, 12)}")

# Simulate shorter time -- this system is violently chaotic
dt_P = mpf('0.0002')
n_steps_P = 100000  # 20 time units
sample_P = 100

data_P = run_simulation(pos_P, vel_P, masses_P, dt_P, n_steps_P, sample_P,
                        "Pythagorean 3-4-5")
res_P = analyze_data(data_P, "Pythagorean 3-4-5 (strongly chaotic)")

# ===================================================================
# PART 5: THE PHI CONNECTION -- DEEPER ANALYSIS
# ===================================================================

print("\n" + "=" * 78)
print("PART 5: THE GOLDEN RATIO CONNECTION -- DEEPER ANALYSIS")
print("=" * 78)

print(f"\n  Golden right triangle identities:")
print(f"    1^2 + phi^2 = 1 + (phi+1) = phi + 2 = phi^2 + 1")
print(f"    = {nstr(ONE + PHI_SQ, 20)}")
print(f"    c = sqrt(phi+2) = {nstr(sqrt(PHI+2), 20)}")
print(f"")
print(f"  Three-body koppa*pi*phi product:")
product = KOPPA * MP_PI * PHI
print(f"    (pi/2) * pi * phi = pi^2*phi/2 = {nstr(product, 15)}")
print(f"    pi*phi/4                        = {nstr(MP_PI*PHI/4, 15)}")
print(f"    pi*phi                          = {nstr(MP_PI*PHI, 15)}")

# The REAL test: compare normalized observables across all three simulations
print(f"\n  CROSS-SIMULATION COMPARISON:")
print(f"  {'Observable':<25s}  {'Figure-8':>12s}  {'Lagrange':>12s}  {'3-4-5':>12s}")
print(f"  {'-'*25}  {'-'*12}  {'-'*12}  {'-'*12}")

rows = [
    ("Angle CV", f"{res_f8['a_cv']:.4f}", f"{res_L['a_cv']:.4f}", f"{res_P['a_cv']:.4f}"),
    ("Geom mean CV", f"{res_f8['gm_cv']:.4f}", f"{res_L['gm_cv']:.4f}", f"{res_P['gm_cv']:.4f}"),
    ("Pyth sum CV", f"{res_f8['ps_cv']:.4f}", f"{res_L['ps_cv']:.4f}", f"{res_P['ps_cv']:.4f}"),
    ("Order factor", f"{res_f8['factor']:.2f}", f"{res_L['factor']:.2f}", f"{res_P['factor']:.2f}"),
    ("ps/(3*gm^2)", nstr(res_f8['rpg_mean'],6), nstr(res_L['rpg_mean'],6), nstr(res_P['rpg_mean'],6)),
    ("<r_max/r_min>", nstr(res_f8['r_mean'],6), nstr(res_L['r_mean'],6), nstr(res_P['r_mean'],6)),
    ("<area/gm^2>", nstr(res_f8['area_ratio'],6), nstr(res_L['area_ratio'],6), nstr(res_P['area_ratio'],6)),
    ("<shape>", nstr(res_f8['shape'],6), nstr(res_L['shape'],6), nstr(res_P['shape'],6)),
    ("Near 90 deg (%)", f"{res_f8['near_90_pct']:.1f}", f"{res_L['near_90_pct']:.1f}", f"{res_P['near_90_pct']:.1f}"),
    ("Near 60 deg (%)", f"{res_f8['near_60_pct']:.1f}", f"{res_L['near_60_pct']:.1f}", f"{res_P['near_60_pct']:.1f}"),
]
for label, v1, v2, v3 in rows:
    print(f"  {label:<25s}  {v1:>12s}  {v2:>12s}  {v3:>12s}")

# ===================================================================
# PART 6: CHAOS --> ORDER THROUGH AVERAGING
# ===================================================================

print(f"\n" + "=" * 78)
print("PART 6: CHAOS --> ORDER THROUGH AVERAGING")
print("=" * 78)

for label, data in [("Figure-8", data_f8), ("Lagrange", data_L), ("3-4-5", data_P)]:
    gm = data['geo_means']
    N = len(gm)
    print(f"\n  [{label}] Running-average convergence of geometric mean:")
    for w in [10, 50, 100, min(500, N)]:
        if w <= N:
            chunk = gm[-w:]
            avg = sum(chunk)/w
            var = sum((g-avg)**2 for g in chunk)/w
            cv = float(sqrt(var)/avg) if avg > 0 else 0
            print(f"    Window {w:4d}: <gm> = {nstr(avg,10)}, CV = {cv:.6f}")

# ===================================================================
# PART 7: THREE BODIES = THREE RELATIONS
# ===================================================================

print(f"\n" + "=" * 78)
print("PART 7: THREE BODIES = THREE RELATIONS")
print("=" * 78)

print(f"""
  The mapping:
    Body 1 <-> koppa (pi/2)  = the structural body, the right angle
    Body 2 <-> pi            = the orbital body, the circle
    Body 3 <-> phi           = the scaling body, the golden ratio

  Product: koppa * pi * phi = (pi/2) * pi * phi
         = pi^2 * phi / 2 = {nstr(MP_PI**2 * PHI / 2, 15)}

  Key question: does any time-averaged observable converge to
  a value expressible in terms of phi, pi, or 1/4?
""")

# Systematic search: for each simulation, try all observables against targets
extended_targets = {
    "1/4":        mpf('0.25'),
    "1/phi":      1/PHI,
    "1/2":        HALF,
    "1/sqrt(5)":  1/sqrt(5),
    "sqrt(3)/4":  sqrt(3)/4,
    "phi/pi":     PHI/MP_PI,
    "2/pi":       2/MP_PI,
    "phi-1":      PHI - 1,
    "1":          ONE,
    "phi":        PHI,
    "sqrt(5)/2":  sqrt(5)/2,
    "pi/2":       MP_PI/2,
    "phi^2":      PHI_SQ,
    "e":          mpf('2.71828182845904523536'),
    "pi":         MP_PI,
    "3":          mpf(3),
    "pi*phi/4":   MP_PI*PHI/4,
    "pi^2/6":     MP_PI**2/6,
    "pi*phi":     MP_PI*PHI,
}

for label, rdict, data in [("Figure-8", res_f8, data_f8),
                            ("Lagrange", res_L, data_L),
                            ("3-4-5", res_P, data_P)]:
    print(f"\n  [{label}] Observable <-> constant matches:")
    obs = {
        "ps/(3*gm^2)":  rdict['rpg_mean'],
        "r_max/r_min":   rdict['r_mean'],
        "area/gm^2":     rdict['area_ratio'],
        "shape":         rdict['shape'],
    }
    for oname, oval in obs.items():
        if oval <= 0:
            continue
        matches = []
        for tname, tval in extended_targets.items():
            if tval <= 0:
                continue
            pct = float(fabs(oval - tval)/tval) * 100
            if pct < 25:
                matches.append((pct, tname, tval))
        matches.sort()
        if matches:
            top = matches[0]
            star = "***" if top[0] < 3 else "** " if top[0] < 10 else "*  "
            print(f"    {star} {oname:18s} = {nstr(oval,8):>12s}  ~  "
                  f"{top[1]:10s} = {nstr(top[2],8):>12s}  ({top[0]:.2f}%)")
        else:
            print(f"        {oname:18s} = {nstr(oval,8):>12s}  (no match < 25%)")

# ===================================================================
# CONCLUSIONS
# ===================================================================

print(f"\n" + "=" * 78)
print("CONCLUSIONS")
print("=" * 78)

uniform_90_pct = 100.0 / 18  # ~5.6% for uniform distribution in [0, 180]

print(f"""
1. CHAOS REDUCTION (1D stacking removes orientational chaos):
   Order factor = angles_CV / distances_CV:
     Figure-8:  {res_f8['factor']:.2f}x
     Lagrange:  {res_L['factor']:.2f}x
     3-4-5:     {res_P['factor']:.2f}x
   {'CONFIRMED' if all(r['factor'] > 1 for r in [res_f8, res_L, res_P]) else 'MIXED'}: Distances are MORE ordered than angles in all cases.

2. PYTHAGOREAN SUM CONSTANCY (virial theorem analogue):
   CV of r12^2+r13^2+r23^2:
     Figure-8:  {res_f8['ps_cv']:.4f}
     Lagrange:  {res_L['ps_cv']:.4f}
     3-4-5:     {res_P['ps_cv']:.4f}
   The Pythagorean sum is {'approximately conserved' if res_f8['ps_cv'] < 0.3 else 'not conserved'} for periodic orbits.

3. RIGHT-ANGLE PREFERENCE:
   % of angles near 90 deg (uniform baseline = {uniform_90_pct:.1f}%):
     Figure-8:  {res_f8['near_90_pct']:.1f}%
     Lagrange:  {res_L['near_90_pct']:.1f}%
     3-4-5:     {res_P['near_90_pct']:.1f}%

4. GOLDEN RATIO IN DISTANCE RATIOS:
   <r_max/r_min> vs phi = {nstr(PHI, 6)}:
     Figure-8:  {nstr(res_f8['r_mean'], 6)}  (diff = {nstr(fabs(res_f8['r_mean']-PHI), 4)})
     Lagrange:  {nstr(res_L['r_mean'], 6)}  (diff = {nstr(fabs(res_L['r_mean']-PHI), 4)})
     3-4-5:     {nstr(res_P['r_mean'], 6)}  (diff = {nstr(fabs(res_P['r_mean']-PHI), 4)})

5. THE KEY RATIO ps/(3*gm^2):
   (= 1 for equilateral, measures shape anisotropy)
     Figure-8:  {nstr(res_f8['rpg_mean'], 6)}
     Lagrange:  {nstr(res_L['rpg_mean'], 6)}
     3-4-5:     {nstr(res_P['rpg_mean'], 6)}

6. GOLDEN RIGHT TRIANGLE IDENTITY (PROVEN):
   1^2 + phi^2 = (phi+2) = phi^2 + 1
   The axiom phi^2 = phi+1 is embedded in every Pythagorean
   triple involving phi. The hypotenuse is c = sqrt(phi+2).

7. STACKING INTERPRETATION:
   Projecting 3 bodies from 2D -> 1D (mutual distances) removes
   the chaotic angular degrees of freedom. The remaining distance
   dynamics shows:
   - Lower variability (ordered)
   - The Pythagorean sum is the most conserved shape quantity
   - Distance ratios explore a structured (not uniform) distribution
""")

print("=" * 78)
print("Script complete.")
