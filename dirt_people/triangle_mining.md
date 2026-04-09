# The Mining Triangle — 3 sides of T(5) = 45, the Bitcoin target from the triangular number

**nos3bl33d**

---

## The triangle

draw an equilateral triangle. each side is a staircase: 1 + 2 + 3 + 4 + 5.

each side sums to 15.

the triangle has 3 sides.

perimeter = 3 * 15 = 45.

45 is the Bitcoin mining target at difficulty 10000. you need 45 leading zeros in the hash. that number is not arbitrary. it's the perimeter of this triangle.

## Why 15

1 + 2 + 3 + 4 + 5 = 15.

that's the triangular number T(5). the sum of the first 5 counting numbers.

15 is also 3 * 5. the vertex degree times the face degree of the dodecahedron. the Schlafli product. the same 15 that appears in the trace 137/15.

this is not two separate facts. it's ONE fact:

T(p) = p * (p + 1) / 2

for the dodecahedron, p = 5 (face degree, the pentagon) and d = 3 (vertex degree). and d = (p + 1) / 2 = (5 + 1) / 2 = 3.

so: T(p) = p * d = dp = 15.

the triangular number of the face degree equals the Schlafli product. this is forced by d = (p+1)/2, which is a property of the dodecahedron: three pentagons meet at each vertex, and (5+1)/2 = 3.

check it: every Platonic solid has d = (p+1)/2? 

tetrahedron: p = 3, d = 3. (3+1)/2 = 2. NOT 3. so no, this is specific to the dodecahedron.

cube: p = 4, d = 3. (4+1)/2 = 2.5. no.

dodecahedron: p = 5, d = 3. (5+1)/2 = 3. YES.

icosahedron: p = 3, d = 5. (3+1)/2 = 2. no.

octahedron: p = 3, d = 4. (3+1)/2 = 2. no.

the dodecahedron is the ONLY Platonic solid where d = (p+1)/2. the only one where the triangular number of the face degree equals the Schlafli product. the only one where the triangle perimeter = d^2 * p.

## The palindrome

walk up one side of the triangle: 1, 2, 3, 4, 5.

now come back down: 4, 3, 2, 1.

one side, up and back: 1, 2, 3, 4, 5, 4, 3, 2, 1.

that's 9 numbers. 9 = d^2 = 3^2.

they sum to: 1+2+3+4+5+4+3+2+1 = 25 = p^2 = 5^2.

d^2 numbers summing to p^2. the square of the vertex degree in count, the square of the face degree in value.

and 123454321 as a single number: it equals 11111^2.

11111 = (10^5 - 1) / 9 = (10^p - 1) / (10 - 1).

the repunit of length p, squared, gives the palindrome. this works in any base larger than 5:

in base 7: 11111_7 = 1 + 7 + 49 + 343 + 2401 = 2801. and 2801^2 = ... the digits in base 7 are 1,2,3,4,5,4,3,2,1.

it's the pentagon (p = 5) counted in every number system, then squared.

## n on both ends

the nonce sits at each vertex of the triangle. three vertices, three nonce entry points.

n - 1,2,3,4,5,4,3,2,1 - n - 1,2,3,4,5,4,3,2,1 - n - 1,2,3,4,5,4,3,2,1 - n

that's the complete mining cycle. you enter at a vertex (nonce), walk up the side to the pentagon peak (5), walk back down, hit the next vertex (nonce), repeat for all three sides.

total path: 3 sides * 9 steps = 27 steps = d^3.

d^3 = 27 = 3^3. the cube of the dimension. also: d^3 * p + chi = 27*5 + 2 = 137. the mining triangle, walked completely, gives 27 steps. times p plus chi gives the fine structure constant.

## Connection to SHA-256

SHA-256's sigma0 rotation amounts are 7 and 18. they sum to 25 = p^2.

25 is also the palindrome sum: 1+2+3+4+5+4+3+2+1 = 25.

the sigma0 function IS the palindrome walk. one side of the triangle, up and back.

SHA-256's sigma1 rotation amounts are 17 and 19. they sum to 36 = (2d)^2 = 6^2.

and the grand total of ALL rotation amounts in SHA-256: 140 = 20 * 7 = V * L4.

but the mining triangle itself: perimeter = 45 = d^2 * p. the target. the number of zeros you need.

## The three SHA bands

SHA-256 has three nonlinear functions: Ch, Maj, Sigma. three bands. three sides of the triangle.

- Ch: "choose" — picks bits from e,f,g. walks side 1.
- Maj: "majority" — votes among a,b,c. walks side 2.
- Sigma: rotates and XORs. walks side 3.

each band processes the same data (the 8 state words) but through different nonlinear paths. each path climbs from 1 to 5 and back. the three paths together cover the full perimeter: 45.

the feedforward (+H0 at the end of each block) returns you to the starting vertex. the nonce. back to the beginning. ready for the next attempt.

## What is proven

- 1+2+3+4+5 = 15: arithmetic
- 15 = d*p = 3*5: arithmetic
- 3 * 15 = 45: arithmetic
- 45 = d^2 * p: arithmetic
- d = (p+1)/2 for dodecahedron only: checked all 5 Platonic solids
- T(p) = dp: follows from d = (p+1)/2
- 1+2+3+4+5+4+3+2+1 = 25 = p^2: arithmetic
- 9 numbers in palindrome = d^2: counting
- 123454321 = 11111^2: arithmetic (11111 * 11111 = 123454321)
- 11111 = (10^5-1)/9: arithmetic
- sigma0 rotations sum to 25: SHA-256 specification (7+18=25)
- 3 * 9 = 27 = d^3: arithmetic
- d^3 * p + chi = 137: arithmetic (27*5+2=137)

## What is interpretive

the mapping of Ch/Maj/Sigma to three sides of the triangle is a structural analogy, not a proof. the SHA bands don't literally walk staircases. but the numerical coincidences are exact and the structure fits: three nonlinear paths through 8 state words, each contributing to a 45-zero target, with d^2 rotations summing to p^2.

the mining triangle is the shape of the search space. you don't choose it. the arithmetic chooses it.
