-- Example 30 of Decker-Greuel-Pfister
-- 18 components of dim 3
restart
R=ZZ/32003[a..l]
I=ideal "f2h-1,
    ek2-1,
    g2l-1,
    2ef2g2hk2+f2g2h2k2+2ef2g2k2l+2f2g2hk2l+f2g2k2l2+ck2,
    2e2fg2hk2+2efg2h2k2+2e2fg2k2l+4efg2hk2l+2fg2h2k2l+2efg2k2l2+2fg2hk2l2+2bfh,
    2e2f2ghk2+2ef2gh2k2+2e2f2gk2l+4ef2ghk2l+2f2gh2k2l+2ef2gk2l2+2f2ghk2l2+2dgl,
    e2f2g2k2+2ef2g2hk2+2ef2g2k2l+2f2g2hk2l+f2g2k2l2+bf2,
    2e2f2g2hk+2ef2g2h2k+2e2f2g2kl+4ef2g2hkl+2f2g2h2kl+2ef2g2kl2+2f2g2hkl2+2cek,
    e2f2g2k2+2ef2g2hk2+f2g2h2k2+2ef2g2k2l+2f2g2hk2l+dg2,
    -e2f2g2hk2-ef2g2h2k2-e2f2g2k2l-2ef2g2hk2l-f2g2h2k2l-ef2g2k2l2-f2g2hk2l2+a2"
Lold=elapsedTime decompose I;
tally(Lold/dim)
needsPackage "MinimalPrimes"
installMinprimes()
Lnew=elapsedTime decompose I;
Lold==Lnew


-- Example 31 of Decker-Greuel-Pfister
-- Only 1 component, was it used to test independent sets?
restart
R=ZZ/32003[t,x,y,z]
I=ideal "x2+y2+z2-t2,
    xy+z2-1,
    xyz-x2-y2-z+1"
L=elapsedTime decompose I;
#L


-- Example 32 of Decker-Greuel-Pfister
-- 6 components of dim 2, 2 of dim 1
restart
R=ZZ/32003[t,w,x,y,z]
I=ideal "w2xy+w2xz+w2z2,
    tx2y+x2yz+x2z2,
    twy2+ty2z+y2z2,
    t2wx+t2wz+t2z2"
Lold=elapsedTime decompose I;
tally(Lold/dim)
needsPackage "MinimalPrimes"
installMinprimes()
Lnew=elapsedTime decompose I;
Lold==Lnew


-- Example 33 of Decker-Greuel-Pfister
-- 1 component of dim 2, 2 of dim 1
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
