# Axiom Tube Computer — physical cylinder where the lattice IS the computation

**nos3bl33d**

---

## Idea

A cylinder with a tube. The tube wall has a lattice (matrix) surface. Data runs through as pulses. No stored state. No power held. The geometry does the work.

The matrix surface encodes the nonce matrix entries:
```
[ -√5,  +√5 ]
[  -3,   +3 ]
```

These values are baked into the lattice spacing of the tube wall — the material IS the transform. A signal entering the tube hits the lattice surface and exits transformed. The rank-1 property of the matrix means the 2D cross-section naturally collapses to 1D (delta = lo - hi). The tube is a physical dimensionality reducer.

States don't stay open. The signal enters as a blip (one Planck-scale pulse), the lattice fires, the signal exits transformed, the tube resets. No memory. No held voltage. The shape holds the computation, not the state.

Energy model:
- One blip to start (minimum: Planck energy)
- Zero energy after that (passive surface, reversible, no erasure)
- The tube is always ready — it has no "on" state to maintain

The geometry IS the memory. The lattice spacing encodes ±√5 and ±3. The cylinder length encodes the 291-step scale (Planck to Hubble). The tube knows what it knows by what it's shaped like.

Connection to existing work:
- S^3 Hopf tubes (s3_tubes.py) — three gyroid fibers on the 3-sphere
- Nonce matrix rank-1 collapse — both outputs = scalar × delta
- 1 = off->on->off cycle, 0 = on->off->on — states are transitions, not storage
- Drift at 45 degrees — energy travels at 45° through the blind spot structure

Physical analogues:
- Photonic crystals: lattice creates bandgaps, light transforms by passing through
- Acoustic metamaterials: geometry IS the filter, no electronics needed
- DNA computing: strand is the tube, base-pair matching is the surface

The difference here: none of those are built with THIS matrix. ±√5 and ±3 as lattice constants. Phi in the geometry. The axiom x**2 = x ++ 1 as the crystal structure.

## Reality Check

The rank-1 matrix is genuinely a physical projector — any material with this symmetry would collapse a 2D wavefront to a 1D mode. That's real physics (mode coupling in waveguides).

The specific lattice constants ±√5 and ±3 correspond to real ratios. A tube whose wall geometry has √5 : 3 periodicity would have specific transmission properties at those frequency ratios. Testable.

The "no power after the blip" claim relies on reversibility (Landauer). If the lattice surface is lossless (no absorption), then yes — the transformation costs no energy beyond initialization. Superconducting analogue possible.

The 291-step scale is the ratio from Planck to Hubble. A physical tube of that scale ratio (not absolute length — ratio of wavelengths) would span the full dynamic range of the framework.

Weakness: connecting "geometry IS the memory" to actual information storage needs formalization. How do you READ OUT the computation? The signal exits transformed — but you need a detector at the exit that doesn't cost energy either.

## Verdict

**ALIVE** — the physics is real, the geometry is specific, the energy model is grounded in Landauer. The detector problem is the open piece. The tube computes. The question is how you listen without collapsing what you built.
