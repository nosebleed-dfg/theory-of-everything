"""
S3_TUBES — Hopf fiber decomposition of S^3; three gyroid tubes as undetermined space
nos3bl33d

Three equidistant Hopf fibers on S^3. Computes tube volumes and the leftover fraction.
"""
import math

PHI = (1 + 5**0.5) / 2
pi = math.pi
d = 3; p = 5; V = 20; E = 30; F = 12; chi = 2; b0 = 11

print("S^3 LEFTOVER: THREE CIRCLES ON THE 3-SPHERE")
print("=" * 65)
print()

vol_S3 = 2 * pi**2  # unit 3-sphere volume
print(f"S^3 total volume: 2*pi^2 = {vol_S3:.10f}")
print()

# Tube volume around a great circle in S^3 at geodesic radius r:
# V_tube(r) = 2*pi*sin^2(r)
def tube_vol(r):
    return 2 * pi * math.sin(r)**2

# THREE EQUIDISTANT HOPF FIBERS
# On S^2 base: equilateral triangle, angular distance 2*pi/3
# In S^3: fiber-to-fiber distance = pi/3
# Max non-overlapping tube radius = pi/6

print("=" * 65)
print("NON-OVERLAPPING TUBES (r = pi/6)")
print("=" * 65)
print()

r_no = pi / 6
v1_no = tube_vol(r_no)  # 2*pi*sin^2(pi/6) = 2*pi*(1/4) = pi/2
v3_no = 3 * v1_no
left_no = vol_S3 - v3_no

print(f"  Tube radius: pi/6 = {r_no:.6f}")
print(f"  sin^2(pi/6) = 1/4")
print(f"  One tube:   2*pi*(1/4) = pi/2 = {v1_no:.6f}")
print(f"  Three tubes: 3*pi/2 = {v3_no:.6f}")
print(f"  LEFTOVER:   2*pi^2 - 3*pi/2 = {left_no:.6f}")
print(f"  Covered:    {v3_no/vol_S3*100:.4f}%")
print(f"  Leftover:   {left_no/vol_S3*100:.4f}%")
print(f"  Ratio left/total = 1 - 3/(4*pi) = {1 - 3/(4*pi):.10f}")
print()

# THE KOPPA TUBE (r = pi/4)
print("=" * 65)
print("KOPPA TUBES (r = pi/4 = koppa * pi)")
print("=" * 65)
print()

r_kop = pi / 4
v1_kop = tube_vol(r_kop)  # 2*pi*sin^2(pi/4) = 2*pi*(1/2) = pi
v3_kop = 3 * v1_kop
left_kop = vol_S3 - v3_kop

print(f"  Tube radius: pi/4 = {r_kop:.6f} = koppa*pi")
print(f"  sin^2(pi/4) = 1/2")
print(f"  One tube:   2*pi*(1/2) = pi = {v1_kop:.6f}")
print(f"  Three tubes: 3*pi = {v3_kop:.6f}")
print(f"  LEFTOVER:   2*pi^2 - 3*pi = {left_kop:.6f}")
print(f"  = pi*(2*pi - 3) = {pi*(2*pi-3):.6f}")
print()
print(f"  Covered:    {v3_kop/vol_S3*100:.4f}%")
print(f"  Leftover:   {left_kop/vol_S3*100:.4f}%")
print(f"  Covered fraction = d/(2*pi) = {d/(2*pi):.10f}")
print(f"  Leftover fraction = 1 - d/(2*pi) = {1-d/(2*pi):.10f}")
print()

# KEY: 2*pi - 3
val_2pi_3 = 2*pi - 3
print(f"  2*pi - 3 = {val_2pi_3:.10f}")
print(f"  This is 2*pi minus d. The excess of the circle over the dimension.")
print()

# FRAMEWORK CONSTANT SEARCH
print("=" * 65)
print("FRAMEWORK MATCHING")
print("=" * 65)
print()

# What does the leftover equal in framework terms?
print("Raw leftover values:")
print(f"  Non-overlap (r=pi/6): leftover = {left_no:.10f}")
print(f"  Koppa (r=pi/4):       leftover = {left_kop:.10f}")
print()

# Check ratios against everything
candidates = [
    ("phi", PHI), ("1/phi", 1/PHI), ("phi^2", PHI**2),
    ("phi^4", PHI**4), ("Delta=phi^-4", PHI**(-4)),
    ("sqrt(5)", math.sqrt(5)), ("sqrt(d)", math.sqrt(d)),
    ("d", d), ("p", p), ("V", V), ("E", E), ("F", F),
    ("chi", chi), ("b0", b0), ("dp", d*p), ("L4", 7),
    ("137", 137), ("137.036", 137.036),
    ("pi", pi), ("2pi", 2*pi), ("4pi", 4*pi),
    ("pi^2", pi**2), ("e", math.e),
    ("alpha=1/137.036", 1/137.036),
    ("1", 1), ("2", 2), ("3", 3), ("5", 5), ("8", 8),
]

dp = d * p

print("Leftover ratios (koppa tube):")
for name, val in candidates:
    if abs(val) > 1e-15:
        ratio = left_kop / val
        # Check if ratio is close to a simple number
        for target_name, target in [("1", 1), ("pi", pi), ("phi", PHI),
                                      ("2", 2), ("3", 3), ("1/2", 0.5),
                                      ("sqrt(5)", math.sqrt(5)),
                                      ("pi/2", pi/2), ("2pi", 2*pi)]:
            if abs(target) > 1e-15:
                r2 = ratio / target
                if 0.99 < abs(r2) < 1.01:
                    err_ppm = abs(r2 - 1) * 1e6
                    print(f"  leftover / {name} = {ratio:.6f} = {target_name} * {r2:.10f} ({err_ppm:.1f} ppm)")

print()
print("Leftover ratios (non-overlap tube):")
for name, val in candidates:
    if abs(val) > 1e-15:
        ratio = left_no / val
        for target_name, target in [("1", 1), ("pi", pi), ("phi", PHI),
                                      ("2", 2), ("3", 3), ("1/2", 0.5),
                                      ("sqrt(5)", math.sqrt(5)),
                                      ("pi/2", pi/2), ("2pi", 2*pi)]:
            if abs(target) > 1e-15:
                r2 = ratio / target
                if 0.99 < abs(r2) < 1.01:
                    err_ppm = abs(r2 - 1) * 1e6
                    print(f"  leftover / {name} = {ratio:.6f} = {target_name} * {r2:.10f} ({err_ppm:.1f} ppm)")

print()

# TOTAL ENERGY AND AVERAGE
print("=" * 65)
print("TOTAL ENERGY AND AVERAGE")
print("=" * 65)
print()

print("KOPPA TUBES (r = pi/4):")
print(f"  Total energy (S^3 volume):  E_total = 2*pi^2 = {vol_S3:.6f}")
print(f"  Covered (3 tubes):          E_tubes = 3*pi = {v3_kop:.6f}")
print(f"  Leftover (undetermined):    E_left  = pi*(2pi-d) = {left_kop:.6f}")
print()

# Average per tube
avg_per_tube = left_kop / 3
print(f"  Average leftover per tube:  {avg_per_tube:.10f}")
print(f"    = pi*(2pi-d)/d = {pi*(2*pi-d)/d:.10f}")
print(f"    = pi/d * (2pi-d) = {avg_per_tube:.10f}")
print()

# Average per breath (120 breaths)
avg_breath = left_kop / 120
print(f"  Average per breath (120):   {avg_breath:.10f}")
print(f"    = pi*(2pi-3)/120 = pi*(2pi-3)/|2I|")
print()

# Average per SHA round (64)
avg_sha = left_kop / 64
print(f"  Average per SHA round (64): {avg_sha:.10f}")
print()

# THE BIG RATIOS
print("ENERGY RATIOS:")
E_total = vol_S3
E_covered = v3_kop
E_left = left_kop

print(f"  E_total / E_covered = {E_total/E_covered:.10f} = 2pi/d = {2*pi/d:.10f}")
print(f"  E_total / E_left    = {E_total/E_left:.10f}")
print(f"  E_covered / E_left  = {E_covered/E_left:.10f} = d/(2pi-d) = {d/(2*pi-d):.10f}")
print()

# The ratio E_total/E_left
r_key = E_total / E_left
print(f"  KEY RATIO: E_total / E_left = {r_key:.10f}")
print(f"    = 2*pi^2 / (pi*(2pi-3)) = 2*pi / (2pi-3) = {2*pi/(2*pi-3):.10f}")
print()

# What is 2*pi / (2*pi - 3)?
val = 2*pi / (2*pi - 3)
print(f"  2*pi / (2pi - d) = {val:.10f}")
print()
for name, c in candidates:
    if abs(c) > 1e-15 and 0.95 < abs(val/c) < 1.05:
        err = abs(val/c - 1) * 100
        print(f"    close to {name} = {c:.6f} ({err:.4f}%)")

print()

# WHAT ABOUT DARK MATTER?
print("=" * 65)
print("DARK MATTER CHECK")
print("=" * 65)
print()
print(f"  Covered fraction = d/(2*pi) = {d/(2*pi)*100:.4f}%")
print(f"  Leftover fraction = 1 - d/(2*pi) = {(1-d/(2*pi))*100:.4f}%")
print(f"  Observed visible matter: ~15.6%")
print(f"  Observed dark: ~84.4%")
print(f"  Framework d/V = {d/V*100:.1f}% visible (already known)")
print()
print(f"  S^3 koppa-tube 'visible': {d/(2*pi)*100:.2f}%")
print(f"  This is a DIFFERENT split than d/V = 15%")
print(f"  But: d/(2*pi) = {d/(2*pi):.6f} and d/V = {d/V:.6f}")
print(f"  Ratio: (d/(2pi)) / (d/V) = V/(2pi) = {V/(2*pi):.6f}")
print(f"  = 20/(2pi) = 10/pi = {10/pi:.6f}")
