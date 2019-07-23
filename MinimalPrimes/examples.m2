-- Example 23 of Decker-Greuel-Pfister
-- 2 components of dim 4, 2 of dim 6, 1 of dim 8
restart
R=ZZ/32003[a..d,u..x]
I=ideal "a+b+c+d,
    u+v+w+x,
    3ab+3ac+3bc+3ad+3bd+3cd+2,
    bu+cu+du+av+cv+dv+aw+bw+dw+ax+bx+cx,
    bcu+bdu+cdu+acv+adv+cdv+abw+adw+bdw+abx+acx+bcx,
    abc+abd+acd+bcd,
    bcdu+acdv+abdw+abcx"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))




-- Example 24 of Decker-Greuel-Pfister
-- 2 components of dim 4, 2 of dim 6, 1 of dim 8
restart
R=ZZ/32003[p,q,s..z]
I=ideal "su,vx,qu,xz,
    stx+ux,
    uv3-uvw+ux,
    -pu2v2+pu2w+qtx,
    tx2y-uv2z+uwz"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 25 of Decker-Greuel-Pfister
-- 2 components of dim 1, 1 of dim 3, 1 of dim 5
restart
R=ZZ/32003[a..h]
I=ideal "59ad+59ah+59dh-705d-1199h,
    330acde+330aceh+330cdeh-407acd-1642ade-1410cde-407ach-407cdh-1642aeh-2398ceh-1642deh,
    -483acd-483ach-483cdh+821ad+705cd+821ah+1199ch+821dh,
    13926abcde+13926abceh+13926bcdeh-9404abcd-9239abde-4968acde-13157bcde-9404abch-9404bcdh-9239abeh-4968aceh-13025bceh-9239bdeh-4968cdeh,
    -cde-377cdh-ceh-deh,
    -54acf-54adf+a+d,
    adfg+a+d"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 26 of Decker-Greuel-Pfister
-- 0 dimensional, 16 reduced points, 24 double points
-- default methods take too long
-- we check the answer using the radical (see Example 34)
-- to speed up radical computation we use the Unmixed option
restart
R=ZZ/32003[a..f]
I=ideal "a2+d2+2ce+2bf+a,
    2ab+2de+2cf+b,
    b2+2ac+e2+2df+c,
    2bc+2ad+2ef+d,
    c2+2bd+2ae+f2+e,
    2cd+2be+2af+f"
needsPackage "MinimalPrimes"
installMinprimes()
L=elapsedTime decompose I;
tally(L/degree)
radical(I,Unmixed=>true) == I
intersect(L) == I


-- Example 27 of Decker-Greuel-Pfister
-- 3 components of dim 4
restart
R=ZZ/32003[a..d,t,x..z]
I=ideal "t-b-d,
    x+y+z+t-a-c-d,
    xz+yz+xt+zt-ac-ad-cd,
    xzt-acd"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 28 of Decker-Greuel-Pfister
-- 2 components of dim 7
restart
R=ZZ/32003[a..h,j..l]
I=ideal "a+c+d+e+f+g+h+j-1,
    -c2k-2cdk-d2k-cek-dek-cfk-dfk-cgk-dgk-egk-fgk-chk-dhk-ehk-fhk+c+d,
    -c2l-cdl-cel-cfl-cgl-dgl-egl-fgl+c2+2cd+d2+cg+dg+ch+dh,
    -b+c+e+g+j"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 29 of Decker-Greuel-Pfister
-- One component
restart
R=ZZ/32003[s,t,u,x,y]
I=ideal "s15,t15,u15,u5-s3tx+s2t2x+s2t2y-st3y"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


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
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 31 of Decker-Greuel-Pfister
-- Only 1 component, was it used to test independent sets?
restart
R=ZZ/32003[t,x,y,z]
I=ideal "x2+y2+z2-t2,
    xy+z2-1,
    xyz-x2-y2-z+1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 32 of Decker-Greuel-Pfister
-- 6 components of dim 2, 2 of dim 1
restart
R=ZZ/32003[t,w,x,y,z]
I=ideal "w2xy+w2xz+w2z2,
    tx2y+x2yz+x2z2,
    twy2+ty2z+y2z2,
    t2wx+t2wz+t2z2"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


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
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


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
