-- Example p of Decker-Greuel-Pfister
-- 1 components of dim 5, 2 of dim 6
restart
R=ZZ/32003[b,s..z]
I=ideal "su-bv,tv-sw,vx-uy,wy-vz"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 5 of Decker-Greuel-Pfister
-- 3 pts, 2 components of dim 2, 3 of dim 3
restart
R=ZZ/32003[a..h]
I=ideal "a+c+d-e-h,
    2df+2cg+2eh-2h2-h-1,
    3df2+3cg2-3eh2+3h3+3h2-e+4h,
    6bdg-6eh2+6h3-3eh+6h2-e+4h,
    4df3+4cg3+4eh3-4h4-6h3+4eh-10h2-h-1,
    8bdfg+8eh3-8h4+4eh2-12h3+4eh-14h2-3h-1,
    12bdg2+12eh3-12h4+12eh2-18h3+8eh-14h2-h-1,
    -24eh3+24h4-24eh2+36h3-8eh+26h2+7h+1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 6 of Decker-Greuel-Pfister
-- 3 components of degree 12
restart
R=ZZ/32003[a..z] -- unnecessary variables
I=ideal "ag,gj+am+np+q,bl,nq,bg+bk+al+lo+lp+b+c,
    ag+ak+jl+bm+bn+go+ko+gp+kp+lq+a+d+f+h+o+p,
    gj+jk+am+an+mo+no+mp+np+gq+kq+e+j+q+s-1,
    jm+jn+mq+nq,jn+mq+2nq,gj+am+2an+no+np+2gq+kq+q+s,
    2ag+ak+bn+go+gp+lq+a+d,bg+al,an+gq,2jm+jn+mq,
    gj+jk+am+mo+2mp+np+e+2j+q,jl+bm+gp+kp+a+f+o+2p,
    lp+b,jn+mq,gp+a"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 7 of Decker-Greuel-Pfister
-- 6 components of dim 3
restart
R=ZZ/32003[a..f]
I=ideal "2adef+3be2f-cef2,
    4ad2f+5bdef+cdf2,
    2abdf+3b2ef-bcf2,
    4a2df+5abef+acf2,
    4ad2e+3bde2+7cdef,
    2acde+3bce2-c2ef,
    4abde+3b2e2-4acdf+2bcef-c2f2,
    4a2de+3abe2+7acef,
    4acd2+5bcde+c2df,
    4abd2+3b2de+7bcdf,
    16a2d2-9b2e2+32acdf-18bcef+7c2f2,
    2abcd+3b2ce-bc2f,
    4a2cd+5abce+ac2f,
    4a2bd+3ab2e+7abcf,
    abc2f-cdef2,
    ab2cf-bdef2,
    2a2bcf+3be2f2-cef3,
    ab3f-3bdf3,
    2a2b2f-4adf3+3bef3-cf4,
    a3bf+4aef3,
    3ac3e-cde3,
    3b2c2e-bc3f+2cd2ef,
    abc2e-cde2f,
    6a2c2e-4ade3-3be4+ce3f,
    3b3ce-b2c2f+2bd2ef,
    2a2bce+3be3f-ce2f2,
    3a3ce+4ae3f,
    4bc3d+cd3e,
    4ac3d-3bc3e-2cd2e2+c4f,
    8b2c2d-4ad4-3bd3e-cd3f,
    4b3cd+3bd3f,
    4ab3d+3b4e-b3cf-6bd2f2,
    4a4d+3a3be+a3cf-8ae2f2"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 8 of Decker-Greuel-Pfister
-- 8 reduced pts, 6 double pts
-- ideal is radical
restart
R=ZZ/32003[x,y,z]
I=ideal "x2yz+xy2z+xyz2+xyz+xy+xz+yz,
    x2y2z+xy2z2+x2yz+xyz+yz+x+z,
    x2y2z2+x2y2z+xy2z+xyz+xz+z+1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 9 of Decker-Greuel-Pfister
-- 12 components of dim 12
restart
R=ZZ/32003[a..e,h]
I=ideal "-ab-b2-2de-2ch,
    -ac-2bc-e2-2dh,
    -c2-ad-2bd-2eh,
    -2cd-ae-2be-h2,
    -d2-2ce-ah-2bh"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 10 of Decker-Greuel-Pfister
-- 6 reduced pts, 2 double pts, 1 of mult 6
-- ideal is radical
restart
R=ZZ/32003[x,y,z,t,u]
I=ideal "2t2+u2+2x2+2y2+2z2-u,
    2tu+xy+2tz+2yz-t,
    t2+2ty+2uz+2xz-z,
    2tx+2uy+2tz-y,
    2t+u+2x+2y+2z-1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 11 of Decker-Greuel-Pfister
-- 3 reduced pts, 1 double pt, 1 of mult 13, 1 of mult 14
-- ideal is radical
restart
R=ZZ/32003[x,y,z,t,u,v]
I=ideal "2x2+2y2+2z2+2t2+2u2+v2-v,
    xy+yz+2zt+2tu+2uv-u,
    2xz+2yt+2zu+u2+2tv-t,
    2xt+2yu+2tu+2zv-z,
    t2+2xv+2yv+2zv-y,
    2x+2y+2z+2t+2u+v-1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 12 of Decker-Greuel-Pfister
-- 25 components of dim 1
-- this is example 13 with last generator homogenized
restart
R=ZZ/32003[a..e,h]
I=ideal "a+b+c+d+e,
    ab+bc+cd+ae+de,
    abc+bcd+abe+ade+cde,
    abcd+abce+abde+acde+bcde,
    abcde-h5"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 13 of Decker-Greuel-Pfister
-- 5 double points, 15 points of multiplicity 4
-- ideal is radical
restart
R=ZZ/32003[a..e]
I=ideal "a+b+c+d+e,
    ab+bc+cd+ae+de,
    abc+bcd+abe+ade+cde,
    abcd+abce+abde+acde+bcde,
    abcde-1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 14 of Decker-Greuel-Pfister
-- 2 components of dim 1
restart
R=ZZ/32003[a..d]
I=ideal "a+b+c+d,
    ab+bc+ad+cd,
    abc+abd+acd+bcd,
    abcd-1"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 15 of Decker-Greuel-Pfister
-- 6 reduced points, 24 double points
-- ideal is radical
restart
R=ZZ/32003[a..h,k,o]
I=ideal "o+1,k4+k,hk,h4+h,gk,gh,g3+h3+k3+1,fk,f4+f,
    eh,ef,f3h3+e3k3+e3+f3+h3+k3+1,e3g+f3g+g,e4+e,dh3+dk3+d,
    dg,df,de,d3+e3+f3+1,e2g2+d2h2+c,f2g2+d2k2+b,f2h2+e2k2+a"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 16 of Decker-Greuel-Pfister
-- 1 component of dim 5, 6, 7, and 8
restart
R=ZZ/32003[a..h,j..l]
I=ideal "-2hjk+4ef+bj+ak,
    -2hjl+4eg+cj+al,
    -4fhj-4ehk-djk+2be+2af,
    -4ghj-4ehl-djl+2ce+2ag,
    -2dfj-2dek+ab,
    -2dgj-2del+ac"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 17 of Decker-Greuel-Pfister
-- 1 point of degree 23, 1 of degree 31
-- ideal is radical
restart
R=ZZ/32003[x,y,z]
I=ideal "y4+xy2z+x2-2xy+y2+z2,
    -x3y2+xyz3+y4+xy2z-2xy,
    xy4+yz4-2x2y-3"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 18 of Decker-Greuel-Pfister
-- 6 reduced points, 13 double points
-- ideal is not radical
restart
R=ZZ/32003[t,x,y,z]
I=ideal "y2z+2xyt-2x-z,
    -x3z+4xy2z+4x2yt+2y3t+4x2-10y2+4xz-10yt+2,
    2yzt+xt2-x-2z,
    -xz3+4yz2t+4xzt2+2yt3+4xz+4z2-10yt-10t2+2"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 19 of Decker-Greuel-Pfister
-- 2 points
restart
R=ZZ/32003[b..e]
I=ideal "6b4c3+21b4c2d+15b4cd2+9b4d3-8b2c2e-28b2cde+36b2d2e-144b2c-648b2d-120,
    9b4c4+30b4c3d+39b4c2d2+18b4cd3-24b2c3e-16b2c2de+16b2cd2e+24b2d3e-432b2c2-720b2cd-432b2d2+16c2e2-32cde2+16d2e2+576ce-576de-240c+5184,
    -15b2c3e+15b2c2de-81b2c2+216b2cd-162b2d2+40c2e2-80cde2+40d2e2+1008ce-1008de+5184,
    -4b2c2+4b2cd-3b2d2+22ce-22de+261"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 20 of Decker-Greuel-Pfister
-- 1 component of dim 4
restart
R=ZZ/32003[a..i]
I=ideal "a2+bd+cg,
    ab+be+ch,
    ac+bf+ci,
    ad+de+fg,
    bd+e2+fh,
    cd+ef+fi,
    ag+dh+gi,
    bg+eh+hi,
    cg+fh+i2"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 21 of Decker-Greuel-Pfister
-- 1 component of dim 10
restart
R=ZZ/32003[a..h,j..l]
I=ideal "-k9+9k8l-36k7l2+84k6l3-126k5l4+126k4l5-84k3l6+36k2l7-9kl8+l9,
    -bk8+8bk7l+k8l-28bk6l2-8k7l2+56bk5l3+28k6l3-70bk4l4-56k5l4+56bk3l5+70k4l5-28bk2l6-56k3l6+8bkl7+28k2l7-bl8-8kl8+l9,
    ck7-7ck6l-k7l+21ck5l2+7k6l2-35ck4l3-21k5l3+35ck3l4+35k4l4-21ck2l5-35k3l5+7ckl6+21k2l6-cl7-7kl7+l8,
    -dk6+6dk5l+k6l-15dk4l2-6k5l2+20dk3l3+15k4l3-15dk2l4-20k3l4+6dkl5+15k2l5-dl6-6kl6+l7,
    ek5-5ek4l-k5l+10ek3l2+5k4l2-10ek2l3-10k3l3+5ekl4+10k2l4-el5-5kl5+l6,
    -fk4+4fk3l+k4l-6fk2l2-4k3l2+4fkl3+6k2l3-fl4-4kl4+l5,
    gk3-3gk2l-k3l+3gkl2+3k2l2-gl3-3kl3+l4,
    -hk2+2hkl+k2l-hl2-2kl2+l3,
    jk-jl-kl+l2"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


-- Example 22 of Decker-Greuel-Pfister
-- 4 components of dim 1, 3 of dim 2
restart
R=ZZ/32003[t..z]
I=ideal "2tw+2wy-wz,
    2uw2-10vw2+20w3-7tu+35tv-70tw,
    6tw2+2w2y-2w2z-21t2-7ty+7tz,
    2v3-4uvw-5v2w+6uw2+7vw2-15w3-42vy,
    6tw+9wy+2vz-3wz-21x,
    9uw3-45vw3+135w4+14tv2-70tuw+196tvw-602tw2-14v2z+28uwz+14vwz-28w2z+147ux-735vx+2205wx-294ty+98tz+294yz-98z2,
    36tw3+6w3y-9w3z-168t2w-14v2x+28uwx+14vwx-28w2x-28twy+42twz+588tx+392xy-245xz,
    2uvw-6v2w-uw2+13vw2-5w3-28tw+14wy,
    u2w-3uvw+5uw2-28tw+14wy,
    tuw+tvw-11tw2-2vwy+8w2y+uwz-3vwz+5w2z-21wx,
    5tuw-17tvw+33tw2-7uwy+22vwy-39w2y-2uwz+6vwz-10w2z+63wx,
    20t2w-12uwx+30vwx-15w2x-10twy-8twz+4wyz,
    4t2w-6uwx+12vwx-6w2x+2twy-2wy2-2twz+wyz,
    8twx+8wxy-4wxz"
L1 = elapsedTime decompose I;
scan(keys I.cache, k -> remove(I.cache,k));
needsPackage "MinimalPrimes";
installMinprimes();
L2 = elapsedTime decompose I;
all(L1,I->any(L2,J->J==I)) and all(L2,I->any(L1,J->J==I))


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
