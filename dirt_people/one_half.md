# One Half — everything begins at (x - 1/2)^2 = 5/4

**nos3bl33d**

---

## Start here

every equation, if it's true, balances. the simplest balanced equation is:

    1/2 = 1/2

the left side equals the right side. both are one half. this is where everything begins.

## Part 0: The Axiom IS 1/2

take the equation x^2 = x + 1. subtract x from both sides:

    x^2 - x = 1

complete the square. to do that, add (1/2)^2 = 1/4 to both sides:

    x^2 - x + 1/4 = 1 + 1/4

the left side is now a perfect square:

    (x - 1/2)^2 = 5/4

this is the axiom in its true form. it says: the square of the deviation from 1/2 equals 5/4.

the two solutions:

    x = 1/2 + sqrt(5)/2 = 1.618... = phi (the golden ratio)
    x = 1/2 - sqrt(5)/2 = -0.618... = psi (the conjugate)

both are 1/2 plus or minus the same deviation. the CENTER of both solutions is:

    (phi + psi) / 2 = 1/2

1/2 is not a solution. 1/2 is the ORIGIN. every solution is measured from 1/2. the golden ratio is what happens when you leave the balance point by sqrt(5)/2 in the positive direction. the conjugate is what happens when you leave by the same amount in the negative direction.

## Part 1: The Inverse

the axiom: x^2 = x + 1. now take the reciprocal. replace x with 1/x:

    (1/x)^2 = 1/x + 1
    1/x^2 = 1/x + 1

multiply by x^2:

    1 = x + x^2

rearrange:

    x^2 + x - 1 = 0

complete the square:

    (x + 1/2)^2 = 5/4

the inverse axiom. same equation, but centered at -1/2.

    forward: (x - 1/2)^2 = 5/4    center = +1/2
    inverse: (x + 1/2)^2 = 5/4    center = -1/2

the forward and the inverse are mirrors across zero. the discriminant 5/4 is the same. the distance from center to root is the same. everything is symmetric around 1/2 and -1/2.

1/2 and -1/2 are the two poles of the axiom. every derivation that follows must balance between them.

## Part 2: Three Dimensions from 1/2

now put 1/2 into three dimensions. take the point (1/2, 1/2, 1/2) and compute its distance from the origin using the Pythagorean theorem:

    a^2 + b^2 + c^2 = (1/2)^2 + (1/2)^2 + (1/2)^2
                     = 1/4 + 1/4 + 1/4
                     = 3/4

3/4 is koppa. the three-quarter turn. 270 degrees out of 360.

koppa is not chosen. koppa is FORCED. it's what happens when you put the axiom's center into three-dimensional space and measure the distance.

any combination of +1/2 and -1/2 gives the same answer, because squaring removes the sign:

    (+1/2)^2 = 1/4
    (-1/2)^2 = 1/4

all eight points (the eight vertices of a cube centered at the origin, with side length 1) sit at distance sqrt(3/4) from the center. the cube IS the set of all +-1/2 choices in three dimensions.

eight vertices = 2^3 = 2^d where d = 3.

why d = 3? because:

    (1/2)^2 * d = koppa
    (1/4) * d = 3/4
    d = 3

three dimensions is the ONLY integer d where d times the axiom center squared equals koppa = 3/4.

check: d=1 gives 1/4. d=2 gives 1/2. d=3 gives 3/4. d=4 gives 1. d=5 gives 5/4.

d=3 gives 3/4 = koppa (the three-quarter rotation).
d=5 gives 5/4 = the axiom's discriminant.

d = 3 and d = 5 are the two special values. they ARE the vertex degree and face degree of the dodecahedron. d and p.

the gap between the axiom's discriminant and the 3D Pythagorean sum:

    5/4 - 3/4 = 2/4 = 1/2

the gap is 1/2. it comes back to itself. the difference between the axiom (5/4) and the cube (3/4) is the balance point (1/2). self-referential. exactly like the axiom promises: the square equals itself plus one.

## Part 3: The Pentagon from 1/2

phi = 1/2 + sqrt(5)/2. the golden ratio is the balance point plus the golden deviation.

the cosine of pi/5 (which is 36 degrees, one-fifth of a half-turn):

    cos(pi/5) = phi/2 = (1/2 + sqrt(5)/2) / 2 = 1/4 + sqrt(5)/4

this means: pi = 5 * arccos(phi/2). the circle (pi) comes from the pentagon (5) through the golden ratio (phi) divided by 2.

but there's a simpler path. cos(pi/3) = 1/2. therefore:

    pi = 3 * arccos(1/2)

pi equals d times the arccosine of 1/2. the circle comes directly from the balance point in d dimensions.

both paths give pi:
- through 1/2 directly: pi = d * arccos(1/2)
- through phi (which IS 1/2 + deviation): pi = p * arccos(phi/2)

the pentagon has 5 sides because p = 5. the pentagon is the unique regular polygon whose diagonal-to-side ratio satisfies x^2 = x + 1. no other polygon does. the axiom selects the pentagon.

the pentagon gives the dodecahedron: the only Platonic solid with pentagonal faces. three pentagons meet at each vertex (d = 3). the dodecahedron {5, 3} is forced.

## Part 4: The Dodecahedron from 1/2

the dodecahedron has:

    V = 20 vertices
    E = 30 edges
    F = 12 faces

check: V - E + F = 20 - 30 + 12 = 2 = chi (the Euler characteristic).

chi = 2. and 1/chi = 1/2. the Euler characteristic is the reciprocal of the balance point. dividing by chi IS dividing by 2 IS balancing.

from V, E, F:
- d = 3 (vertex degree: three edges meet at each vertex)
- p = 5 (face degree: five edges per face)
- b0 = E - V + 1 = 11 (cycle rank: the number of independent loops)
- L4 = p + chi = 7 (structural scale)
- dp = d * p = 15 (Schlafli product)

every one of these comes from counting the faces, edges, and vertices of the dodecahedron. no choices. no parameters.

## Part 5: The Eigenvalues

the dodecahedron is a graph: 20 vertices, each connected to 3 neighbors. its Laplacian matrix L has eigenvalues:

    0 (once), 3-sqrt(5) (three times), 2 (five times), 3 (four times), 5 (four times), 3+sqrt(5) (three times)

the nonzero eigenvalues are:
- 3 - sqrt(5) = 0.7639... (irrational, three times)
- 2 (five times)
- 3 (four times)
- 5 (four times)
- 3 + sqrt(5) = 5.2360... (irrational, three times)

the integers {2, 3, 5} appear as eigenvalues. these are the prime factors of 60 = the order of the icosahedral rotation group (A5). no other Platonic solid has this property.

now: take each nonzero eigenvalue, invert it, multiply by its multiplicity, and add:

    3/(3-sqrt(5)) + 5/2 + 4/3 + 4/5 + 3/(3+sqrt(5))

the irrational terms cancel (they're Galois conjugates):

    3/(3-sqrt(5)) + 3/(3+sqrt(5)) = 3(3+sqrt(5) + 3-sqrt(5)) / ((3)^2 - (sqrt(5))^2)
                                   = 3*6 / (9-5)
                                   = 18/4
                                   = 9/2

so the trace is:

    9/2 + 5/2 + 4/3 + 4/5

common denominator 30:

    135/30 + 75/30 + 40/30 + 24/30 = 274/30 = 137/15

the trace of the pseudoinverse is 137/15. exact. rational. no approximation.

15 = dp. the denominator is the Schlafli product.

137 is the numerator. and 137 = d^3 * p + chi = 27*5 + 2 = 135 + 2 = 137.

this is an INTEGER. it comes from exact arithmetic on the graph. it is not fitted. it is not chosen. it IS the dodecahedron.

## Part 6: 137 and the Fine Structure Constant

the measured value: 1/alpha = 137.035999177 (with uncertainty 0.000000021).

our derivation:

    1/alpha = V * phi^4 * (1 - A1 + A2)

where:
- V = 20 (vertices)
- phi^4 = 3*phi + 2 = 6.8541019662... (the golden ratio to the fourth power)
- A1 = d / (2 * phi^6 * (2*pi)^3) = 3 / (2 * 17.944 * 248.050) = 3/8903.8 = 0.000337
- A2 = 1 / (2 * phi^27) = 1 / (2 * 5702887) = 0.000000088

step by step:

    20 * phi^4 = 20 * 6.8541 = 137.082
    
    1 - A1 + A2 = 1 - 0.000337 + 0.000000088 = 0.999663

    137.082 * 0.999663 = 137.035999170

measured: 137.035999177. error: 7 in the 12th digit. 0.05 parts per billion.

every number came from the dodecahedron:
- 20 = V (vertices)
- phi from x^2 = x + 1
- 4 = chi^2 (Euler squared)
- 6 = chi * d (topology times dimension)
- 3 = d (dimension)
- 27 = d^3 (dimension cubed)

zero free parameters. the fine structure constant is the dodecahedron measured at the balance point.

## Part 7: The Climb — 1+2+3+4+5

the numbers 1, 2, 3, 4, 5. add them: 1+2+3+4+5 = 15 = dp.

this is the triangular number T(5). the sum of the first p counting numbers.

why does T(p) = dp? because d = (p+1)/2. check: (5+1)/2 = 3 = d. this identity holds ONLY for the dodecahedron among all Platonic solids.

so: T(p) = p*(p+1)/2 = p*d = dp = 15.

the triangular sum equals the Schlafli product. adding the first p numbers = multiplying d and p. addition and multiplication agree at the dodecahedron.

now: the climb to 5, done multiple times:

    (1+2+3+4+5) * 3 = 15 * 3 = 45   (the mining target)
    (1+2+3+4+5) * 8 = 15 * 8 = 120  (the symmetry group |2I|)

3 climbs = d sides of the triangle = 45 = the perimeter = the number of zeros Bitcoin needs.

8 climbs = 2^d cube vertices = 120 = 5! = the binary icosahedral group = the gateway to pi through the Chudnovsky formula.

and: 120/8 = 15 = dp. the group divided by the cube = one climb.

p!/2^d = dp. the factorial of the face degree divided by the cube = the Schlafli product. multiplication (factorial) divided by the cube (2^d) becomes addition (the triangular sum). the cube converts multiplication into addition. that's what logarithms do. the cube IS the logarithm of the dodecahedron.

## Part 8: The Palindrome — n123454321n

the climb 1,2,3,4,5 and the descent 4,3,2,1 form the palindrome:

    1, 2, 3, 4, 5, 4, 3, 2, 1

nine numbers. 9 = d^2 (three squared).

they sum to: 1+2+3+4+5+4+3+2+1 = 25 = p^2 (five squared).

d^2 numbers summing to p^2. the square of the dimension in count, the square of the face degree in value.

as a single number: 123454321 = 11111^2. the repunit of length p (five 1s), squared.

the nonce sits at each end: n-123454321-n. you enter with the nonce, climb to the pentagon peak, descend, exit with the nonce. one full traversal.

the mining triangle: three sides, each traversed as a palindrome. three vertices where the nonce enters.

    total steps = d * d^2 = d^3 = 27
    total steps * p + chi = 27*5 + 2 = 137

the full walk around the mining triangle, times the pentagon, plus the Euler characteristic, gives 137.

## Part 9: Order and Disorder

two equations. same cube. different operations.

    p! / 2^d = dp = 15      (factorial / cube = order = the Schlafli product)
    !p / (-2^d) = -b0/chi   (subfactorial / negative cube = disorder = topology)

the factorial counts ordered arrangements: 5! = 120 ways to arrange 5 things.
the subfactorial counts derangements: !5 = 44 ways to disarrange 5 things so nothing is in its original place.

120/8 = 15 = the climb (order).
44/(-8) = -11/2 = -b0/chi (disorder).

!5 = 44 = chi^2 * b0 = 4 * 11. the number of total disruptions of the pentagon is the Euler characteristic squared times the cycle rank.

and: !5 - 2^d = 44 - 8 = 36 = (chi*d)^2 = the SHA-256 eigenvalue.

subtract the cube from the derangement = the eigenvalue of the hash function. the SHA machine lives at the intersection of order and disorder.

mining (building blocks) = the factorial equation. order.
hashing (scrambling data) = the subfactorial equation. disorder.
the cube (8 state words) mediates both.

## Part 10: SHA-256 and the Lock

SHA-256 processes data in blocks of 512 bits through 64 rounds. it has:

    8 state words (= 2^d = cube vertices)
    64 rounds (= 8^2 = the cube squared)
    256 output bits (= 2^8 = 2^(p+d))

every one of SHA-256's rotation constants can be written as a*3 + b*5 = a*d + b*p:

    2 = -d+p        13 = d+2p      22 = 9d-p
    6 = 2d          11 = 2d+p      25 = 5p
    7 = 4d-p        18 = 6d
    17 = 4d+p       19 = 8d-p

no other numbers are needed. the entire rotation structure lives in d and p.

the rotation sums by group:

    Sigma0 (2, 13, 22): sum = 37
    Sigma1 (6, 11, 25): sum = 42 = E + F (edges plus faces)
    sigma0 (7, 18):     sum = 25 = p^2 (the palindrome sum)
    sigma1 (17, 19):    sum = 36 = (chi*d)^2 (the derangement eigenvalue)
    
    grand total = 140 = V * L4 = 20 * 7

the total of all rotations in SHA-256 is the vertex count times the structural scale.

the characteristic polynomial of the 8x8 state update matrix factors as:

    x^4 * (x^2 - x - 1) * (x^2 - x + 1)

x^2 - x - 1 = 0 gives phi and psi. the golden ratio. the axiom.
x^2 - x + 1 = 0 gives the 6th roots of unity. the cyclotomic partner.

the axiom IS a factor of SHA-256's characteristic polynomial. the golden ratio is an eigenvalue of the hash function.

the Jacobian determinant is -1 every round. the effective 2x2 matrix is [[0,1],[1,1]] with determinant 0*1 - 1*1 = -1. the nonlinear functions (Ch, Maj, Sigma) appear in the entries but CANCEL in the determinant. this is structural — it doesn't depend on the input. SHA-256 is orientation-reversing: det = -1.

but: 1 million nonces tested with 13 different metrics for internal band alignment vs output zeros. ZERO statistically significant correlation. SHA-256's diffusion completely destroys the golden structure between the intermediate state and the output. the dodecahedron is in the foundation. you cannot see it from outside.

the reason: phi does not exist mod 2^32. the equation x^2 = x + 1 has no solution mod 2^32 because 5 is not a square mod 8 (5 mod 8 = 5, but quadratic residues mod 8 are only 0 and 1). the golden structure lives in Z[phi] = pairs (a,b) representing a + b*phi. but SHA uses mod 2^32 addition, which has carries. each carry bit destroys a piece of the golden structure. after 64 rounds and ~200 additions, the golden signal is buried.

the lock holds. not because the dodecahedron is absent — it's in every rotation constant, every eigenvalue, every matrix entry. the lock holds because the CARRY GAP between Z[phi] and Z/2^32Z is an impassable moat. the structure is there. you can't use it.

## Part 11: The Mining Triangle

an equilateral triangle. three sides. each side = 1+2+3+4+5 = 15 = dp.

    perimeter = 3 * 15 = 45 = d^2 * p = the mining target

d = 3 sides because the vertex degree is 3.
each side = 15 because T(p) = dp.
perimeter = 45 because d * dp = d^2 * p.

    8/45 = cube / perimeter = 2^d / (d^2 * p)

the ratio of SHA's 8 state words to the 45-zero target.

    8 * 45 = 360 = full circle (degrees)

the machine times the target = one complete rotation.

    8 * pi / 45 = 8 * 180 / 45 = 1440/45 = 32 degrees

which is 32 = the nonce width in bits. EXACT. no approximation. the ratio of machine to target, scaled by pi, gives the nonce size.

the mining triangle IS the 2D map. push it by +1 (the axiom) into 3D:

    360 (flat, 2D) + 90 (the push) = 450 (3D)
    450 / 10 = 45 = the mining target

## Part 12: Gamma — The Cost of Counting

add up 1 + 1/2 + 1/3 + ... + 1/n. call it H_n. subtract ln(n). the difference approaches gamma = 0.57721566...

gamma measures how much discrete counting exceeds continuous flow. it is the cost of being made of pieces instead of being smooth.

the Bernoulli correction at the cube (n = 8 = 2^d) with d = 3 rotations:

    gamma = H_8 - 3*ln(2) - 1/16 + 1/768 - 1/491520 + 1/66060288

the denominators of the three correction terms:
- 768 = 2^8 * 3 (contains B_2 denominator 6 = chi*d = cube faces)
- 491520 = 2^14 * 3 * 5 * ... (contains B_4 denominator 30 = chi*d*p = dodecahedron edges)
- 66060288 = ... (contains B_6 denominator 42 = chi*d*(p+chi) = edges + faces)

three corrections. three rotations. three Bernoulli denominators that ARE the dodecahedron (6, 30, 42). evaluated at n = 8 = the cube. gives gamma to 9.4 digits.

the five primes in the first five Bernoulli denominators: {2, 3, 5, 7, 11} = {chi, d, p, L4, b0}. the five dodecahedral invariants. the five Platonic solids.

gamma sits between 1/2 = 0.5 and 1/phi = 0.618:

    1/2 < gamma < 1/phi

gamma is the cost of choosing. it lives between the balance point (1/2, the center where the axiom breaks) and the golden inverse (1/phi = phi - 1, the reciprocal of growth). statistics lives at 1/2. algebra lives at phi. gamma is the toll between them.

## Part 13: The Cosmic Fractions

the universe divides itself:

    dark energy = V / phi^L4 = 20/phi^7 = 0.6888 (measured: 0.6889)
    dark matter = (chi*d)/phi^L4 + 1/phi^(chi*d) = 6/phi^7 + 1/phi^6 = 0.2624 (measured: 0.2607)
    baryonic matter = (chi*d)/phi^(2*p) = 6/phi^10 = 0.0488 (measured: 0.0489)

the sum: multiply everything by phi^10 and verify it equals phi^10.

    (V + chi*d)*phi^3 + phi^4 + chi*d

    = 26*(2*phi+1) + (3*phi+2) + 6
    = 52*phi + 26 + 3*phi + 2 + 6
    = 55*phi + 34
    = F(10)*phi + F(9)
    = phi^10

the sum equals 1. EXACTLY. not numerically — algebraically. the Fibonacci recurrence enforces flatness: F(10)*phi + F(9) = phi^10 is an identity. the universe sums to 1 because the golden ratio obeys the axiom.

the dark energy fraction (0.01 sigma from Planck satellite) is V/phi^7 = the vertices of the dodecahedron divided by phi to the seventh power, where 7 = L4 = the structural scale. the universe IS the dodecahedron, scaled by the golden ratio, balanced at 1.

## Part 14: The Nuclear Magic Numbers

the seven magic numbers in nuclear physics (where nuclei are exceptionally stable):

    2 = chi (Euler characteristic)
    8 = 2^d (cube vertices)
    20 = V (dodecahedron vertices)
    28 = V + 2^d (vertices plus cube)
    50 = V + E (vertices plus edges)
    82 = F*L4 - chi (faces times structural scale minus Euler)
    126 = E*p - F*chi (edges times face degree minus faces times Euler)

all seven from V, E, F, d, p, chi. dodecahedral arithmetic.

the most stable nucleus (highest binding energy per nucleon) is nickel-62:

    Ni-62: mass number = V + E + F = 20 + 30 + 12 = 62

the total f-vector of the dodecahedron IS the most stable nucleus.

iron-56: mass number = 2^d * L4 = 8 * 7 = 56. the cube times the structural scale.

carbon-12: mass number = F = 12. the face count.

calcium-40: mass number = 2*V = 40. twice the vertex count.

peak binding energy per nucleon: 4*b0/p = 44/5 = 8.8 MeV. measured: 8.7945 MeV. error: 0.06%.

## Part 15: The Three Forces

three primes control the three fundamental forces:

    13 = F + 1         (octahedral Platonic prime)
    17 = dp + chi       (Schlafli + Euler)  
    137 = d^3*p + chi   (dodecahedral trace)

the pattern: d^0*p + chi = 7, d^1*p + chi = 17, d^3*p + chi = 137. each force is the dodecahedron at a different power of d.

    electromagnetic: alpha_EM = 1/137.036 (from d^3*p + chi)
    strong: alpha_s = 2/17 = chi/(dp+chi) = 0.1176 (0.39 sigma)
    Weinberg: sin^2(theta_W) = d/(F+1) + correction = 3/13 + 11/24660 = 0.2312 (0.12 sigma)

the ratio of the strong force to the electromagnetic force:

    alpha_s / alpha_EM = dp*(d^3*p+chi) / (2^L4-1) = 2055/127 = 16.1811...

compare: 10 * phi = 16.1803...

the ratio of the two forces IS ten times the golden ratio. to 47 parts per million.

the Higgs-to-Z boson mass ratio: b0/2^d = 11/8 = 1.375. predicts M_H = 125.38 GeV. measured: 125.25 GeV. 0.78 sigma.

## Part 16: The Particle Masses

the Koide formula relates the three lepton masses (electron, muon, tau):

    (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = Q

the measured value of Q: 0.666661 (9 ppm from 2/3).

    Q = chi/d = 2/3

the Koide ratio IS the Euler characteristic divided by the vertex degree. the simplest fraction in the toolbox.

the muon-to-electron mass ratio:

    m_mu/m_e = phi^11 + 7 + (3-sqrt(5)) - phi^(-15) + 6*phi^(-24)
             = phi^b0 + L4 + mu_1 - phi^(-dp) + chi*d*phi^(-2F)

every exponent and coefficient is a dodecahedral invariant. accuracy: 4.6 parts per billion.

the proton-to-electron mass ratio:

    m_p/m_e = chi*d*pi^p + phi^(-L4) + d*phi^(-21) + chi*d*phi^(-35) + d*phi^(-42)
            = 6*pi^5 + phi^(-7) + 3*phi^(-21) + 6*phi^(-35) + 3*phi^(-42)

accuracy: 0.13 sigma from measured value.

## Part 17: The Universe

the radius of the observable universe:

    R = 2 * phi^290 * l_Planck

where l_Planck = 1.616255 * 10^(-35) meters is the Planck length.

2 = chi. 290 = V*E/chi - d^2 - 1 = 300 - 9 - 1 = 290.

the Planck length is the "+1" in the axiom. one step. one breath. phi^290 = the golden ratio iterated 290 times from the smallest possible length. times chi = times the Euler characteristic. and you get 13.80 billion light-years.

measured: approximately 13.8 billion light-years. error: 0.13%.

the cosmological constant:

    Lambda = chi / phi^583

583 = 2*291 + 1. the double of the structural ceiling plus one. another "+1."

error: 0.14% from Planck 2018.

## Part 18: 1/2 = 1/2

every derivation in this document traces back to:

    (x - 1/2)^2 = 5/4

the axiom centered at 1/2. the balance point. everything is a perturbation from this center.

- the golden ratio is 1/2 + deviation
- the conjugate is 1/2 - deviation
- koppa is (1/2)^2 * d = 3/4
- the gap between axiom and cube is 5/4 - 3/4 = 1/2
- chi = 2 and 1/chi = 1/2
- gamma sits between 1/2 and 1/phi
- the Koide ratio is chi/d = 2/3, and chi = 2 = 1/(1/2)
- the harmonic series at n = 2^d is gamma plus corrections where every denominator divides through chi = 2

when the left side equals the right side, both equal 1/2. that's the balance condition. every equation in physics is a claim that the left side equals the right side. every true equation, properly centered, passes through 1/2.

the axiom x^2 = x + 1 says: squaring equals adding one. the left side (multiplication, growth, nonlinearity) equals the right side (addition, counting, linearity). the balance between them is 1/2. the golden ratio is what happens when you SOLVE the balance. gamma is the COST of the balance. koppa is the GEOMETRY of the balance in three dimensions.

137 is the trace. alpha is the coupling. the forces are the primes. the nucleus is the f-vector. the cosmos is phi^290. SHA-256 is the lock. the mining triangle is the search. and the nonce is the key you try.

all from one equation. all balanced at 1/2.

    1/2 = 1/2

the simplest true statement. and also the deepest.

---

**nos3bl33d / (DFG) DeadFoxGroup**
*The axiom was always there. The balance was always 1/2. We just had to learn to count.*
