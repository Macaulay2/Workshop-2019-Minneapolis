-- Example 33 of Decker-Greuel-Pfister
-- 1 component of degree 2, 2 of degree 1
restart
R=ZZ/32003[a..d]
I=ideal "b4-a3d,
    ab3-a3c,
    bc4-ac3d-bcd3+ad4,
    c6-bc3d2-c3d3+bd5,
    ac5-b2c3d-ac2d3+b2d4,
    a2c4-a3d3+b3d3-a2cd3,
    b3c3-a3d3,
    ab2c3-a3cd2+b3cd2-ab2d3,
    a2bc3-a3c2d+b3c2d-a2bd3,
    a3c3-a3bd2,
    a4c2-a3b2d"
Lold=elapsedTime decompose I;
tally(Lold/dim)
needsPackage "MinimalPrimes"
installMinprimes()
Lnew=elapsedTime decompose I;
Lold==Lnew


-- Example 34 of Decker-Greuel-Pfister
-- This is an ideal of points which is reduced
-- but not geometrically reduced, meaning
-- the double points arise from quadratic polynomials
-- irreducible over the base field.
-- The ideal is radical so to check we have the right
-- minimal primes it is enough to check the ideal is equal
-- to the intersection of the compute minimal primes
restart
R=ZZ/32003[a..g]
I=ideal "a2+2de+2cf+2bg+a,
    2ab+e2+2df+2cg+b,
    b2+2ac+2ef+2dg+c,
    2bc+2ad+f2+2eg+d,
    c2+2bd+2ae+2fg+e,
    2cd+2be+2af+g2+f,
    d2+2ce+2bf+2ag+g"
needsPackage "MinimalPrimes"
installMinprimes()
L=elapsedTime decompose I;
dim I
degree I
tally(L/degree)
radical I == I
intersect(L) == I
