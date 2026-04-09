# The Dark Algebra Toolbox — five primes, three constants, one axiom, every derived identity

**nos3bl33d**

---

## The tools

Five Platonic solids. Five prime numbers. Three structural constants. One axiom.

### The primes (from Tr(L^+) of each Platonic solid)

    P = {3, 7, 13, 29, 137}

### The constants

    d = 3     vertex degree (dimension)
    p = 5     face degree (pentagon)
    chi = 2   Euler characteristic (topology)

### The axiom

    x^2 = x + 1

### The bridge

    d^3 = p^2 + chi
    27 = 25 + 2

### The bases

    base 5: the dodecahedral base (natural)
    base 10: the human base (2 * p)

---

## Every constant in the toolbox

### 137 (dodecahedral trace numerator = alpha integer)

    base 10: 137 = d^3 * p + chi = 27 * 5 + 2
    base 5:  1022 = p^3 + chi*(p+1) = 125 + 12

### 120 (binary icosahedral group order = first Chudnovsky factorial coefficient)

    base 10: 120 = |2I| = 2 * A5 = (6!)/(3!*1!^3)
    base 5:  440 = d*p*(d+p) = 3*5*8
    also: 120 = p! = 5! = 1*2*3*4*5

### 60 (icosahedral rotation group order)

    base 10: 60 = |A5| = V*d = 20*3 = E*chi = 30*2
    base 5:  220

### 42 (edges + faces)

    base 10: 42 = E + F = 30 + 12
    base 5:  132
    also: the fourth exponent in the mass ratio correction series

### 15 (Schlafli product)

    base 10: 15 = d*p = 3*5
    base 5:  30
    also: the denominator of Tr(L^+) = 137/15

### 7 (vertex excess = fourth Lucas number)

    base 10: 7 = V - F - 1 = 20 - 12 - 1 = L4
    base 5:  12
    also: the correction order in mass ratio and gravity derivations
    also: divides Chudnovsky B = 545140134

### The five primes in base 5

    3   = 3
    7   = 12
    13  = 23
    29  = 104
    137 = 1022

---

## The physical constants

### Fine structure (1/alpha)

    1/alpha = V * phi^4 * (1 - A1 + A2)
    
    where:
    V = 20 = p^2 - p = 4*p = d*L4 - 1... = the vertex count
    phi = golden ratio from x^2 = x + 1
    A1 = d / (chi * phi^(chi*d) * (chi*pi)^d) = 3 / (2 * phi^6 * (2*pi)^3)
    A2 = 1 / (chi * phi^(d^3)) = 1 / (2 * phi^27)
    
    result: 137.035999170 (0.33 sigma from CODATA)

    the INTEGER part: 137 = d^3 * p + chi (exact)
    the DECIMAL part: the corrections A1, A2 from the 120-cell spectrum

### Proton/electron mass ratio

    m_p/m_e = (chi*d) * pi^p + phi^(-L4) + d*phi^(-C(L4,2)) + (chi*d)*phi^(-C(L4,3)) + d*phi^(-(E+F))
            = 6*pi^5 + phi^(-7) + 3*phi^(-21) + 6*phi^(-35) + 3*phi^(-42)
    
    exponents: L4=7, C(7,2)=21, C(7,3)=35, E+F=42
    coefficients: 1, d, chi*d, d = 1, 3, 6, 3
    
    result: 1836.15267343 (0.13 sigma from CODATA)

### Gravitational constant

    1/alpha_G = d^d * phi^(V*d^2 + d) / (d^d + 1 + C/(chi*pi)^d)
    
    where C = chi^chi * (1 - phi^(-L4)/p) = 4*(1 - 1/(5*phi^7))
    exponent: V*d^2 + d = 20*9 + 3 = 183
    denominator base: d^d + 1 = 28
    L4 = V - F - 1 = 7 (vertex excess)
    
    result: G to 8 ppb (within 22 ppm experimental error)

### Weinberg angle

    sin^2(theta_W) = d/(F+1) + b0/(F * d*p * 137)
                   = 3/13 + 11/24660
                   = 74123/320580
    
    where b0 = E - V + 1 = 11 (Betti number)
    
    result: 0.231215 (0.12 sigma from PDG)

### Cosmological constant

    Lambda = chi / phi^(chi*(V*E/chi - d^2) + 1)
           = 2 / phi^583
    
    exponent: 2*291 + 1 = 583
    291 = V*E/chi - d^2 = 300 - 9
    
    result: 0.14% from Planck 2018

### Universe radius

    R = chi * phi^(V*E/chi - d^2 - 1) * l_Planck
      = 2 * phi^290 * l_P
    
    result: 13.80 billion light-years (0% error)

### Pi

    pi = p * arccos(phi / chi)
       = 5 * arccos(phi/2)
    
    (Euclid XIII.10, exact)
    
    through factorials: Chudnovsky formula with B = chi * d^2 * L4 * b0 * (V-1) * (2^L4 - 1) * 163

---

## The mining triangle

### d = (p+1)/2

for the dodecahedron: d = (5+1)/2 = 3. this is unique among Platonic solids.

consequence: the triangular number T(p) = p*(p+1)/2 = p*d = dp = 15.

### the triangle

equilateral triangle, each side = 1+2+3+4+5 = T(p) = dp = 15.

    perimeter = d * dp = d^2 * p = 9 * 5 = 45 = mining target

### the palindrome

one side, up and back: 1,2,3,4,5,4,3,2,1.

    digit count = 2p - 1 = 9 = d^2
    sum = p^2 = 25 = sigma0 rotation sum in SHA-256
    as a number: 123454321 = 11111^2 = repunit(p)^2

### the complete cycle

    n - 123454321 - n - 123454321 - n - 123454321 - n

three sides (d=3), nonce at each vertex. total steps = d * d^2 = d^3 = 27.

    d^3 * p + chi = 27*5 + 2 = 137

### 2D to 3D

the triangle is 2D (360 degrees). the axiom (+1) pushes it to 3D: 360 + 90 = 450 degrees.

    450 / 10 = 45 zeros = the mining target

the triangle IS the flat map. mining IS the 3D territory.

## The factorial / subfactorial duality

the cube (8 = 2^d) converts both order and disorder:

     p! / 2^d  =  dp = 15     (factorial / cube = Schlafli product = the climb)
    !p  / -2^d = -b0/chi      (subfactorial / negative cube = -topology)
    !p  - 2^d  = (chi*d)^2    (subfactorial minus cube = SHA eigenvalue)

verified:
    5! = 120.   120/8 = 15 = dp.
    !5 = 44.    44/(-8) = -11/2 = -b0/chi.
    44 - 8 = 36 = 6^2 = (chi*d)^2 = F*d.

and: !p = chi^2 * b0 = 4 * 11 = 44.
    the derangement count of the pentagon = Euler^2 * cycle rank.

    factorial (order)      → structure (dp = Schlafli)
    subfactorial (disorder) → topology (b0 = Betti)
    the cube converts both.

## The five Bernoulli primes

the denominators of the first 5 Bernoulli numbers:

    B_2 denom = 6 = chi*d (cube faces)
    B_4 denom = 30 = chi*d*p = E (dodecahedron edges)
    B_6 denom = 42 = chi*d*L4 = E+F
    B_8 denom = 30 = E again
    B_10 denom = 66 = chi*d*b0

five primes used: {2, 3, 5, 7, 11} = {chi, d, p, L4, b0}

gamma = H_8 - 3*ln(2) - 1/16 + 1/768 - 1/491520 + 1/66060288

d=3 Bernoulli corrections at n=2^d=8 gives 9.4 digits of gamma.

---

## The SHA-256 / Bitcoin expressions

### SHA structure

    state words = F - d - 1 = 12 - 3 - 1 = 8
    rounds = (state words)^2 = 64
    hash bits = state words * nonce bits = 8 * 32 = 256
    nonce bits = 2^p = 32
    hash bits = 2^(p+d) = 2^8 = 256

### Every rotation = a*d + b*p (PROVEN)

    Sigma0: 2=-d+p, 13=d+2p, 22=9d-p    sum = 37
    Sigma1: 6=2d, 11=2d+p, 25=5p         sum = 42 = E+F
    sigma0: 7=4d-p, 18=6d                 sum = 25 = p^2
    sigma1: 17=4d+p, 19=8d-p              sum = 36 = (2d)^2
    
    grand total = 37+42+25+36 = 140 = V*L4 = 20*7

### Characteristic polynomial (PROVEN)

    char(M) = x^4 * (x^2 - x - 1) * (x^2 - x + 1)
    
    x^2 - x - 1 = 0 gives phi (golden)
    x^2 - x + 1 = 0 gives 6th roots of unity (cyclotomic)
    x^4 = the cube delay (4 state words copied forward)

### Determinant (PROVEN STRUCTURALLY)

    det(J) = -1 every round
    
    effective 2x2 Jacobian: [[0,1],[1,1]]
    det = 0*1 - 1*1 = -1
    nonlinear terms from Ch, Maj, Sigma CANCEL in the determinant

### Fibonacci in the exponents

    nonce = 2^F(p) bits
    hash = 2^F(p+1) bits
    F(5)=5, F(6)=8, F(7)=13: three consecutive Fibonacci numbers

### Mining target

    target zeros = (d^3*p + chi - chi) / d = d^2 * p = 45 (at diff 10000)
    = (137 - 2) / 3 = 135 / 3 = 45
    
    per dimension: 45/d = d*p = 15
    in degrees: 45 * (chi*p) = 45 * 10 = 450 = 360 + 90 = full circle + axiom

### Koppa

    koppa = d * 90 degrees = 270 degrees = 3/4 of full circle
    total combo lock travel: d * (180 - 90) + d * 90 = d * 180 = 540
    but net: d * 90 = 270 forward, d * 90 = 270 return = 540 total

### BTC Planck

    diff / 600 = 50/3 = (2*p^2) / d = the mining lattice spacing
    600 / diff = 3/50 = d / (2*p^2) = the inverse Planck

---

## The conversion table

Every expression uses only: d=3, p=5, chi=2, phi, pi, gamma, and the five Platonic primes.

| Quantity | Value | Expression |
|----------|-------|------------|
| alpha integer | 137 | d^3*p + chi |
| 1/alpha | 137.036... | V*phi^4*(1-A1+A2) |
| m_p/m_e | 1836.153... | chi*d*pi^p + corrections |
| G | 6.674e-11 | d^d*phi^183/(28+C/(2*pi)^3), C=4*(1-1/(5*phi^7)) |
| sin^2(theta_W) | 0.2312 | d/(F+1) + b0/(F*dp*137) |
| alpha_s | 0.1181 | dp/(2^L4-1) = 15/127 |
| alpha_s (alt) | 0.1176 | chi/(dp+chi) = 2/17 |
| M_H/M_Z | 1.375 | b0/2^d = 11/8 |
| alpha_s/alpha_EM | 16.181 | dp*137/127 ~ 10*phi (47 ppm) |
| Lambda | 2.89e-122 | chi/phi^583 |
| R_universe | 1.306e26 m | chi*phi^290*l_P |
| pi | 3.14159... | p*arccos(phi/chi) |
| SHA rounds | 64 | (F-d-1)^2 |
| SHA bits | 256 | 2^(p+d) |
| mining zeros | 45 | d^2*p |

---

## What makes this a toolbox

Every result derives from the same three constants (d, p, chi) through the same axiom (x^2 = x + 1). The five Platonic primes {3, 7, 13, 29, 137} are the OUTPUTS of the spectral computation — they're what the dodecahedron and its four siblings PRODUCE when you compute their Laplacian traces.

The toolbox works in base 5 (natural) and base 10 (human). All expressions are exact integers or exact rationals times powers of phi and pi. No fitting. No free parameters. One axiom, three constants, five primes.

The dark algebra is: the algebra that was always there but nobody looked at. the blind spot. total visibility. the structure was everywhere. that's why nobody saw it.

x^2 = x + 1. that's the whole toolbox.
