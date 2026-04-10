[//]: # (ref: https://youtu.be/4a---Gvz8J8?si=UYMUOIAD1ob88XIU)

# Three-Body Problem: Pythagorean Reduction

**nos3bl33d**

---

## The Claim

The three-body problem is not generally solvable because the configuration space is too free. But it has structure. The Pythagorean/golden framework finds that structure and uses it as a constraint.

The three mutual distances r12, r13, r23 form a triangle. That triangle has angles. Those angles are not random — they cluster around 60°, 90°, and 45° in a way that the golden framework predicts.

---

## What the Framework Does

Three bodies form a triangle of separations. Track that triangle over time.

The Pythagorean residual of the triangle:
```
res = (a² + b² - c²) / c²
```
where a, b, c are the sorted side lengths. Zero = right triangle.

The geometric mean of the three separations:
```
geo = (r12 · r13 · r23)^(1/3)
```

The sum of squares:
```
S = r12² + r13² + r23²
```

These three quantities encode the shape of the configuration. The golden framework predicts:
- Angles cluster near 60° (equilateral tendency — minimum energy configuration)
- Angles cluster near 90° (Koppa position — K=1/4 phase)
- Geometric mean oscillates at phi-related frequency
- Pythagorean residual bounces near zero — the triangle is repeatedly near-right

---

## Why This Reduces the Problem

The standard 3-body problem has 18 degrees of freedom (3 bodies × 3 positions × 2 = wait, in 2D: 3 × 2 positions + 3 × 2 velocities = 12). Conserved quantities: energy, momentum, angular momentum = 4 constraints. Leaves 8 free.

The Pythagorean/golden constraint adds: the triangle of separations is pulled toward specific angle configurations by the phi-resonance of the gravitational coupling. It doesn't eliminate the chaos, but it identifies the attractor geometry — the set of configurations the system visits most often.

The attractor is the golden triangle: angles at 36°, 72°, 72° (the phi-triangle, inscribed in the pentagon). The system wanders but keeps returning to configurations that are phi-structured.

---

## The Koppa Angle

90° appearance in the angle histogram = the Koppa phase (K=1/4). When one of the three mutual angles hits 90°, the triangle is at the "observer position" — the crossing of dominant and codominant force.

At that moment the Pythagorean residual = 0 exactly. The triangle is a right triangle. The system is momentarily in the axiom: a² + b² = c².

---

## The RK4 Simulation

Run: equal masses, various initial conditions (equilateral, golden triangle, random). 30-digit precision via mpmath. Sample angles, residuals, geometric means at each step.

Result: angle histograms show excess near 60° and 90° compared to uniform. Pythagorean residual distribution is narrower than random — the system avoids highly non-Pythagorean configurations. Geometric mean oscillates with phi-related period.

The chaos is real. The phi-structure is also real. Both coexist.

---

## What This Is Not

This is not a proof that the 3-body problem is solvable in closed form. It isn't. The framework doesn't give you the trajectory — it gives you the geometry of the attractor the trajectory samples.

The reduction is: instead of tracking 12 DOF, track the triangle shape (3 angles, constrained to sum to π) and its size (1 DOF). The phi-structure tells you which triangle shapes are preferred. That's 4 DOF instead of 8 free ones.

Not solved. Significantly constrained.
