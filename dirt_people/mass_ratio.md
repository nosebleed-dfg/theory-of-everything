# The Proton-Electron Mass Ratio — 6*pi^5 plus golden corrections to 0.0002 ppm

**nos3bl33d**

---

## What is this number?

the proton is 1836.15267343 times heavier than the electron. this ratio sets the size of atoms, the strength of chemical bonds, the frequency at which molecules vibrate. change it by a fraction of a percent and chemistry stops working.

nobody derives this number. in the standard model it's just an input -- you measure it and plug it in.

in 1951, somebody noticed that 6 * pi^5 = 1836.118... is suspiciously close. accurate to 19 parts per million. a curiosity for seventy years. no correction mechanism, no explanation for why 6 * pi^5.

we provide the corrections. they come from the dodecahedron.

---

## The formula

    m_p / m_e = 6*pi^5 + phi^(-7) + 3*phi^(-21) + 6*phi^(-35) + 3*phi^(-42)

where phi = (1 + sqrt(5))/2 is the golden ratio.

that is the whole formula. no free parameters. every number traces back to the dodecahedron.

---

## Where the base comes from

**6 * pi^5**

the 6 = chi * d = 2 * 3. chi is the Euler characteristic of the sphere (V - E + F = 20 - 30 + 12 = 2). d is the vertex degree of the dodecahedron (every vertex touches 3 edges). pi enters through the geometry of the pentagonal faces -- cos(pi/5) = phi/2, which is Euclid's identity connecting the axiom x^2 = x + 1 to circular measure.

    6 * pi^5 = 2 * 3 * 306.0196... = 1836.1181...

this is 18.8 parts per million away from the measured value. close enough to notice, too far off to be an answer.

---

## Where the exponents come from

the key number is 7. it is V - F - 1 = 20 - 12 - 1 = 7. the vertex excess of the dodecahedron. also the fourth Lucas number L4.

the correction exponents are built from 7 using binomial coefficients:

    C(7,1) = 7      first exponent
    C(7,2) = 21     second exponent
    C(7,3) = 35     third exponent

check them:

    C(7,1) = 7!/1!6! = 7.            correct.
    C(7,2) = 7!/2!5! = 7*6/2 = 21.   correct.
    C(7,3) = 7!/3!4! = 7*6*5/6 = 35. correct.

the fourth exponent breaks the binomial pattern. C(7,4) = 35, not 42. so where does 42 come from?

    42 = E + F = 30 + 12

edges plus faces. a dodecahedral invariant. the series terminates at the edge-face sum, not at the next binomial coefficient. why? because C(7,4) = C(7,3) = 35, so the binomial sequence starts repeating. instead of doubling up at 35, the series jumps to the next dodecahedral landmark: E + F = 42.

the exponent gaps tell the same story:

    21 - 7  = 14 = 2 * 7 = chi * L4
    35 - 21 = 14 = 2 * 7 = chi * L4
    42 - 35 = 7  = L4 itself

two steps of 14, then one step of 7. the last step halves because the series is terminating.

---

## Where the coefficients come from

**this is where the old version was wrong.** the coefficients are:

    {1, 3, 6, 3}

the old claim was that these are Pascal row 4. they are not. Pascal row 4 is {1, 4, 6, 4, 1}. completely different.

the actual pattern is {1, d, chi*d, d}:

    1     = 1        (unity)
    d     = 3        (vertex degree = spatial dimension)
    chi*d = 2*3 = 6  (Euler characteristic times dimension)
    d     = 3        (vertex degree again)

this is a palindrome wrapped around chi*d:  1, d, chi*d, d. the middle coefficient is the largest because it combines topology (chi) with geometry (d). the outer coefficients are bare dimension.

the 6 is NOT "3 choose 2" or any binomial coefficient. it is chi * d = 2 * 3. two times three. Euler characteristic times vertex degree. a dodecahedral product, not a combinatorial one.

why this pattern? the correction terms represent shells of influence around a vertex. the first shell sees 1 channel. the second sees d = 3 channels (one per edge). the third sees chi * d = 6 channels (topology times geometry -- the maximum). the fourth sees d = 3 again as the influence wraps back. the palindromic structure {1, d, chi*d, d} reflects the symmetry of the dodecahedral graph distance function.

---

## The computation, step by step

grab a calculator. every number below is exact to the digits shown.

**phi = (1 + sqrt(5)) / 2 = 1.6180339887...**

**Step 1: the base**

    6 * pi^5 = 6 * 306.01968... = 1836.11810871...

    error from measured: 18.8 parts per million

**Step 2: add the first correction**

    phi^(-7) = 1 / phi^7 = 1 / 29.0344... = 0.03444185...

    running total: 1836.11810871 + 0.03444185 = 1836.15255057

    error from measured: 66.9 parts per billion

the first correction killed four orders of magnitude of error.

**Step 3: add the second correction**

    phi^(-21) = 1 / phi^21 = 1 / 24476.0...

    3 * phi^(-21) = 3 / 24476.0 = 0.00012256905...

    running total: 1836.15267313...

    error from measured: 0.16 parts per billion

another two and a half orders of magnitude gone.

**Step 4: add the third correction**

    phi^(-35) = 1 / phi^35

    6 * phi^(-35) = 0.00000029079...

    running total: 1836.15267342528...

    error from measured: 2.6 parts per trillion

**Step 5: add the fourth correction**

    3 * phi^(-42) = 0.00000000501...

    running total: 1836.15267343029...

    error from measured: 0.16 parts per trillion

---

## The result table

| what you have | value | error |
|---|---|---|
| 6*pi^5 alone | 1836.1181... | 18.8 ppm |
| + phi^(-7) | 1836.15255... | 66.9 ppb |
| + 3*phi^(-21) | 1836.15267313... | 0.16 ppb |
| + 6*phi^(-35) | 1836.15267342528... | 2.6 ppt |
| + 3*phi^(-42) | 1836.15267343029... | 0.16 ppt |
| measured (CODATA 2022) | 1836.15267343(11) | --- |
| deviation | 0.003 sigma | within noise |

each correction term kills roughly two to three orders of magnitude of error. the series converges geometrically in powers of phi^(-7) = 0.0344, so each step shrinks the residual by a factor of ~29.

---

## The eigenvalue ratio (independent cross-check)

there is a second formula, completely different in structure, that gets the same number from the dodecahedron.

the dodecahedral graph Laplacian has six distinct eigenvalues:

    0           (once)
    3 - sqrt(5) (three times)   <-- this is mu_1, the smallest nonzero
    2           (five times)
    3           (four times)
    5           (four times)    <-- this is mu_5, the largest
    3 + sqrt(5) (three times)

take the ratio of the largest to the smallest nonzero eigenvalue, raised to the fourth power:

    m_p/m_e ~ (mu_5 / mu_1)^4 = (5 / (3 - sqrt(5)))^4

rationalize the denominator:

    5 / (3 - sqrt(5)) = 5 * (3 + sqrt(5)) / (9 - 5) = 5*(3 + sqrt(5)) / 4

or in terms of phi: since 3 - sqrt(5) = 2/phi^2, we get:

    mu_5 / mu_1 = 5 / (2/phi^2) = 5*phi^2 / 2

so:

    (5*phi^2/2)^4 = 625 * phi^8 / 16 = 1835.106

    measured: 1836.153
    error: 0.057% = 570 ppm

this is a much cruder approximation -- 570 ppm instead of sub-ppt. but it uses ZERO fitting. just two raw eigenvalues and one exponent (d+1 = 4, the embedding dimension). the exponent 4 is not chosen to fit -- it is the number of dimensions you need to embed a dodecahedron in Euclidean space (3 spatial + 1 for the constraint surface).

the eigenvalue formula and the correction series are not competing answers. they are two views of the same structure: the eigenvalue ratio gives the coarse answer (which eigenvalues dominate), and the correction series gives the fine structure (how the intermediate eigenvalues modify the ratio).

---

## What is proven, what is not

**proven (exact arithmetic, verifiable on any calculator):**
- the formula evaluates to 1836.15267343029... (verified to 60 decimal digits)
- the exponents 7, 21, 35 are binomial coefficients C(7,k) for k=1,2,3
- the exponent 42 = E + F = 30 + 12, a dodecahedral invariant
- the coefficients {1, d, chi*d, d} = {1, 3, 6, 3} are dodecahedral constants
- the eigenvalue ratio (5*phi^2/2)^4 = 1835.106 (independent cross-check)
- every constant in the formula comes from the dodecahedron or from x^2 = x + 1

**not proven (conjectural):**
- WHY the proton-electron mass ratio should equal this particular combination
- the physical mechanism behind the correction channel structure
- why 6*pi^5 is the base term (the connection between pi^5 and QCD is not derived)
- why the exponent series terminates at E + F instead of continuing the binomial pattern

the number is exact. the explanation is incomplete. the arithmetic is closed. the physics is open.
