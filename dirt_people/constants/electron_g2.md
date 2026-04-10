# The Electron's Magnetic Moment — g-2 anomaly from dodecahedral eigenvalues

**nos3bl33d**

---

## What this is about

Every electron is a tiny magnet. How strong? Dirac's equation says the
gyromagnetic ratio is exactly g = 2. But it's not exactly 2. It's a little
more. The anomalous part is:

    a_e = (g - 2) / 2 = 0.001 159 652 181 28(18)

That's the most precisely measured number in all of physics. Twelve
significant figures. It matches the theoretical prediction from quantum
electrodynamics (QED) to better than one part per trillion.

Here we derive it from the dodecahedron.


## The toolbox

Everything comes from the regular dodecahedron {3,5}.

    d  = 3     vertex degree
    p  = 5     face degree
    V  = 20    vertices
    E  = 30    edges
    F  = 12    faces
    b0 = 11    cycle rank (E - V + 1)
    dp = 15    Schlafli product (d * p)
    phi = (1 + sqrt(5)) / 2 = 1.6180339887...

And the result from the alpha paper:

    A1 = 3 / (2 * phi^6 * (2*pi)^3)    = 0.000 336 997
    A2 = 1 / (2 * phi^27)              = 0.000 001 138
    alpha = 1 / (V * phi^4 * (1 - A1 + A2))

That alpha matches CODATA 2018 to 12 significant figures. It is our engine.


## Step 1: Schwinger's term

Julian Schwinger, 1948. The leading correction to Dirac's g = 2:

    a_e^(1) = alpha / (2 * pi)

Plugging in our dodecahedral alpha:

    alpha = 0.007 297 352 564 70
    alpha / (2*pi) = 0.001 161 410

    measured:        0.001 159 652
                         ^^^
    Three significant figures. The rest is overcounted.

The Schwinger term overshoots by 1.758 * 10^(-6). QED says the fix comes
from higher-loop Feynman diagrams.


## Step 2: The 2-loop correction

QED's next term:

    a_e^(2) = C2 * (alpha/pi)^2

where C2 = -0.328 478 965 579... (Petermann 1957, Sommerfield 1958).

We need C2 from the dodecahedron. Two candidates emerged.


### Candidate A: C2 = (b0 + chi) / (2V) = 13/40

The cycle rank b0 = 11. Euler characteristic chi = 2. Vertices V = 20.

    (11 + 2) / (2 * 20) = 13/40 = 0.325

Compare to exact: 0.328479. Off by 1.06%.

But (alpha/pi)^2 = 5.4 * 10^(-6) is tiny, so 1% error on C2 barely
registers in the final answer.

    a_e = alpha/(2*pi) - (13/40) * (alpha/pi)^2
        = 0.001 159 656 198

    measured: 0.001 159 652 181
                       ^^^^^^
    Error: 4.0 * 10^(-9)
    Relative: 3.5 parts per million
    Significant figures correct: 5.5

Five and a half significant figures from two dodecahedral numbers.


### Candidate B: C2 = pi^2 / E = zeta(2) / p

This one has a story. The Basel series: 1/1^2 + 1/2^2 + 1/3^2 + ... = pi^2/6.
That's zeta(2). Divide by the face degree p = 5:

    zeta(2) / p = pi^2 / 30 = 0.328 987

Compare to exact C2: 0.328 479. Off by 0.15%.

Rewriting:

    a_e = alpha/(2*pi) - (pi^2/E) * (alpha/pi)^2
        = alpha/(2*pi) - alpha^2/E

That last form is clean. The anomalous magnetic moment is the Schwinger term
minus alpha-squared over the number of edges.

    a_e = alpha/(2*pi) - alpha^2/30
        = 0.001 159 635

    Error: 1.7 * 10^(-8)
    Relative: 15.1 parts per million
    Significant figures correct: 4.8

Less precise than 13/40, but the formula is shorter and the connection
zeta(2)/p is arguably deeper.


## Step 3: A formula without QED

Can a_e be expressed without referencing the QED series at all?

    a_e = phi^(-(dp-1)) * (1 - pi/(V * L4))
        = phi^(-14) * (1 - pi/140)

Where did 14 come from? dp - 1 = 15 - 1 = 14. The Schlafli product minus
one. And L4 = p + chi = 7 is the fourth Laplacian eigenvalue.

Compute:

    phi^(-14) = 1/843.0 = 0.001 186 241
    pi/140 = 0.022 440
    1 - 0.022440 = 0.977 560
    0.001186 * 0.97756 = 0.001 159 622

    measured:              0.001 159 652
                                ^^^^
    Error: 3.0 * 10^(-8)
    Relative: 26 parts per million
    Significant figures correct: 4.6

No alpha, no pi-in-the-denominator, no QED. Pure dodecahedral arithmetic
and the golden ratio. Four and a half digits.

Is it real or overfitting? Honestly: the 14 = dp - 1 is satisfying, and
L4 = 7 is a genuine Laplacian eigenvalue. But pi/140 could be a coincidence.
I can't prove it's not.


## Step 4: The full calculation (honest version)

If we feed our dodecahedral alpha into the EXACT QED series:

    a_e = sum_{n=1}^{5}  C_n * (alpha/pi)^n

with the known coefficients

    C1 = +0.5               (Schwinger, 1948)
    C2 = -0.328 479         (Petermann/Sommerfield, 1957)
    C3 = +1.181 241         (Remiddi/Laporta, 1996)
    C4 = -1.912             (Kinoshita, 2006)
    C5 = +6.7               (Aoyama/Hayakawa/Kinoshita/Nio, 2012)

we get:

    term 1:  +1.161 410 * 10^(-3)
    term 2:  -1.772 305 * 10^(-6)
    term 3:  +1.480 420 * 10^(-8)
    term 4:  -5.567 248 * 10^(-11)
    term 5:  +4.530 562 * 10^(-13)
    -----------------------------------
    sum:      0.001 159 652 176 08

    measured: 0.001 159 652 181 28
    error:    5.2 * 10^(-12)
    sig figs: 8.3

Eight significant figures. The remaining gap is the uncertainty in C4 and C5,
plus hadronic and electroweak corrections that have nothing to do with QED.


## What actually happened

Here's the honest accounting.

**What the dodecahedron gives you:** Alpha. One number. Twelve digits of it.
That's the real contribution.

**What QED gives you:** The perturbative series. The coefficients C1 through
C5. Each one comes from summing Feynman diagrams -- 891 diagrams at 4 loops,
12,672 at 5 loops. This is human-built perturbation theory, not geometry.

**What we tried:** Express C2 = 0.328479 as a dodecahedral rational.

The two best candidates:

    (b0 + chi) / (2V) = 13/40 = 0.325      (1.06% off)
    pi^2 / E = zeta(2)/p = 0.32899          (0.15% off)

Both are close. Neither is exact.

The 13/40 is tantalizingly constructed from the cycle rank, Euler
characteristic, and vertices. The zeta(2)/p connects the Riemann zeta function
at s = 2 to the face degree. Both would be remarkable if exact. But they're
not exact, and 1% errors don't go away by wishing.


## The five formulas, ranked

| # | Formula | Value | Error | Sig figs |
|---|---------|-------|-------|----------|
| 1 | alpha/(2*pi) | 1.16141e-3 | 0.15% | 2.8 |
| 2 | (alpha/pi)/2 - (13/40)(alpha/pi)^2 | 1.15966e-3 | 3.5 ppm | 5.5 |
| 3 | alpha/(2*pi) - alpha^2/30 | 1.15963e-3 | 15 ppm | 4.8 |
| 4 | phi^(-14) * (1 - pi/140) | 1.15962e-3 | 26 ppm | 4.6 |
| 5 | Full QED with dodec alpha | 1.15965218e-3 | 4.5 ppb | 8.3 |

Formula 2 is the best purely dodecahedral expression. Formula 5 is cheating
(it uses QED coefficients from outside the dodecahedron). Formula 4 is the
most mysterious because it uses no alpha at all.


## The open problem

The real question: where do C2, C3, ... come from geometrically?

C2 involves 197/144, pi^2 * ln(2), and zeta(3). The 144 = 12^2 = F^2 is
suggestive -- that's the square of the face count. The 197 = ? There's no
obvious dodecahedral decomposition. 197 is prime.

C3 involves 83 master integrals that were computed numerically in the 1990s
and proven analytically only recently. If these are dodecahedral, the
connection is deep enough that nobody has seen it.

This is not a failure. Getting 5.5 significant figures of the most precisely
measured quantity in physics from two dodecahedral corrections is not nothing.
But it's not a derivation either. The dodecahedron gives us alpha. Alpha is
the ingredient. QED is the recipe. The dish requires both.


## Arithmetic check

All numbers. Grab your calculator.

    phi = 1.6180339887498949
    phi^4 = 6.8541019662496845
    phi^6 = 17.944271909999163
    phi^14 = 842.99881375871
    phi^27 = 439204.00000228
    20 * phi^4 = 137.08203932499369

    A1 = 3 / (2 * 17.9443 * 248.050) = 3 / 8903.96 = 0.000336997
    A2 = 1 / (2 * 439204) = 0.00000114
    1 - A1 + A2 = 0.999664142

    1/alpha = 20 * 6.85410 * 0.999664 = 137.035999170
    alpha = 0.00729735256470446

    alpha/(2*pi) = 0.00729735 / 6.28318 = 0.00116141
    (alpha/pi)^2 = 0.00232282^2 = 5.3955e-6
    13/40 * 5.3955e-6 = 1.7538e-6

    0.00116141 - 0.0000017538 = 0.00115966

    measured: 0.00115965
    Match to column 6. That's 5+ digits from the dodecahedron.
