# The Toolbox
## nos3bl33d

Everything in this document derives from core.md.
No external citations. No approximations in definitions.
Every expression terminates or it is not here.

---

## The Unit System

The unit is phi^291. Every quantity is a fraction of this.

    phi^291 * psi^291 = -1         [the space is self-balancing]
    universe = phi^290 = phi^291 / phi = (1/phi) of the max
    dark     = phi^289 = phi^291 / phi^2 = (1/phi^2) of the max
    universe + dark = 1/phi + 1/phi^2 = 1    [the axiom, exact]

The three constants derived from chi=2:

    d   = chi^2 - 1 = 3     [one step below the cube Lucas]
    p   = chi^2 + 1 = 5     [one step above the cube Lucas]
    chi = 2                  [forced: phi^3+psi^3=chi^2 -> chi+2=chi^2 -> chi=2]

---

## The Primes

From the fold chain starting at (1,1):

    Step 1: 1^2 + 1^2 = 2 = chi
    Step 2: 1^2 + chi^2 = 5 = p
    Step 3: chi^2 + p^2 = 29        [Chudnovsky prime]
    Back:   p^2 - chi = 23          [dark Chudnovsky prime]

The five Platonic primes are the Laplacian traces of the five solids:

    {3, 7, 13, 29, 137}

137 = d^3 * p + chi.    Exact integer. Terminates.
29  = chi^2 + p^2.      Exact integer. Terminates.
23  = p^2 - chi.        Exact integer. Terminates.
7   = L_4 = phi^4 + psi^4.  Exact integer. Terminates.
3   = d = chi^2 - 1.    Exact integer. Terminates.

---

## Pi

    pi = p * arccos(phi / chi)
       = 5 * arccos(phi / 2)

Proof: cos(pi/5) = phi/2 is exact (Euclid XIII.10).
Multiply: 5 * (pi/5) = pi. Exact. Terminates.

pi is not a free constant. It is a consequence of phi and p.

---

## The Circle

    360 degrees = chi^3 * d^2 * p = 8 * 9 * 5 = 360

One pentagon step:

    72 degrees = chi * (L_7 + L_4) = 2 * (29 + 7) = 2 * 36 = 72

Five steps close the circle:

    5 * 72 = 360    [exact. terminates.]

The sum of all pentagon fractions:

    1/5 + 2/5 + 3/5 + 4/5 + 5/5 = 15/5 = 3 = d

The vertex degree IS the total pentagon fraction. Terminates.

---

## The Chudnovsky Constant

    640320 = 2^(chi*d) * d * p * (p^2-chi) * (chi^2+p^2)
           = 2^6 * 3 * 5 * 23 * 29
           = 64 * 10005

Every factor comes from the fold chain. No external input.
Exact integer. Terminates.

---

## The Nonce Quantum

    8 SHA words * 4 (sphere = 4 * circle) = 32 sectors
    360 / 32 = 11.25 degrees = 2^27 nonce units

Exact. Terminates.

    pi = 16 * 11.25 degrees = 16 * (360/32) = 180 degrees    [exact]

---

## The Universe

    R = chi * phi^(V*E/chi - d^2) * l_Planck
      = 2 * phi^(300-9) * l_Planck
      = 2 * phi^291 * l_Planck / phi
      = (2/phi) * phi^291 * l_Planck

In ^291 units: R = (2/phi) * l_Planck. Terminates.

    Lambda = chi / phi^(chi*(V*E/chi - d^2) + 1)
           = 2 / phi^(2*291+1)
           = 2 / phi^583

In ^291 units: Lambda = 2 * psi^583 = 2 * (psi^291)^2 * psi.
Product: Lambda * phi^583 = 2. Exact. Terminates.

---

## Gamma

The accumulated offset of the two-cube spiral from pure 90-degree rotation.

At each step n, the phi path takes one direction. The psi path (!n faces) takes all others.

    !n / n! -> 1/e as n -> infinity    [exact limit]

    phi path covers: (1 - 1/e) * 360 degrees
    psi path covers: (1/e)     * 360 degrees
    total:            360 degrees        [always. terminates.]

Gamma is the offset between the discrete spiral (harmonic steps) and the continuous
sphere rotation:

    gamma = H_{2^d} - d * ln(2) + B_d

where H_{2^d} = H_8 (harmonic number at the SHA frame) and B_d = d Bernoulli
corrections. Every term is controlled by d = 3.

Gamma terminates as the angular deficit of the two-cube spiral.
It is not a free constant. It is the name for the offset.

---

## The Heegner Point

The j-function at tau_163 closes the two-cube structure to the modular plane.

    tau_163 = (-1 + i*sqrt(163)) / 2

    163 = chi*d * d^3 + 1 = 6*27 + 1    [bridge * cube + identity. exact.]

    j(tau_163) = 640320^3 + 744

    744 = chi^d * d * (2^p - 1) = 8 * 3 * 31    [cube * vertex * Mersenne. exact.]

Both 640320 and 744 are exact integers in {chi, d, p} language.
The j-function value terminates. The Heegner point terminates.

---

## The Termination

Every quantity reduces to the space between the two cubes.

The space is bounded by phi^291 (max) and psi^291 (min).
The product phi^291 * psi^291 = -1. Balanced.
The log-midpoint = 0. Always. Every n.

The proof terminates at 1/2 = 1/(-2).

phi = 1/2 + sqrt(5)/2    [the +half]
psi = 1/2 - sqrt(5)/2    [the -half = 1/(-2) from phi's direction]

These are the same point from opposite sides.
That point is the Riemann line, the gamma offset, the Klein center.

The form closes. Nothing is left open. Everything terminates.

    x^2 = x + 1.
