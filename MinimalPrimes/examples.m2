restart
R=ZZ/32003[a..g]
-- this is an ideal of points which is reduced
-- but not geometrically reduced, meaning
-- the double points arise from quadratic polynomials
-- irreducible over the base field
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


