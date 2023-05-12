restart
needsPackage "CellularResolutions"
S = QQ[x,y,z,w];
I = monomialIdeal (y*w, x*y*z, x^2*y, z^4*w);
DeltaT = taylorComplex I

T = chainComplex DeltaT

v1 = newCell({},y*w);
v2 = newCell({},x*y*z);
v3 = newCell({},x^2*y);
v4 = newCell({},z^4*w);
e12 = newCell({v1,v2});
e13 = newCell({v1,v3});
e23 = newCell({v2,v3});
e14 = newCell({v1,v4});
f123 = newCell({e12,e13,e23});
Delta = cellComplex(S, {f123,e14})
C = chainComplex Delta 

prune HH C
prune HH Delta
isMinimal DeltaT 
isMinimal Delta 

Cnonred = chainComplex(Delta, Reduced => false)

H = hullComplex I 
scarf = scarfComplex I
prune HH scarf

H2 = hullComplex monomialIdeal {x^2*z, x*y*z, y^2*z, x^3*y^5, x^4*y^4, x^5*y^3}
isMinimal H2 

scarf2 = scarfComplex monomialIdeal {x*y, y*z, z*w, w*x}
prune HH scarf2 
isMinimal scarf2

S2 = cellComplexSphere(QQ,2)
cells(S2) 
prune HH S2
T3 = cellComplexTorus(QQ,3)
prune HH T3 
Z2 = ZZ/2;
RP3 = cellComplexRPn(Z2,3)
prune HH RP3

facePoset RP

P = convexHull id_(ZZ^4);
Pcellular = cellComplex(S,P)

verts = cells(0,Pcellular);
H = hashTable  {verts#0 => y*w, verts#1 => x*y*z, verts#2 => x^2*y, verts#3 => z^4*w};
relabeledP = relabelCellComplex(Pcellular, H)
cells(0,relabeledP)/(c -> cellLabel c)

R = QQ[a,b];
P1 = convexHull matrix {{5,3},{1,2}};
P2 = convexHull matrix {{3,2},{2,3}};
P3 = convexHull matrix {{2,0},{3,7}};
PC = polyhedralComplex {P1,P2,P3};
v = vertices PC;
H = hashTable for i to (numColumns v - 1) list (v_i => a^(lift(v_i_0,ZZ))*b^(lift(v_i_1,ZZ)))
PCcellular = cellComplex(R,PC, Labels => H)
apply(cells(0,PCcellular),cellLabel)
apply(cells(1,PCcellular),cellLabel)

needsPackage "NormalToricVarieties"
X = toricProjectiveSpace(1) ** toricProjectiveSpace(2);
S = ring X;
B = ideal X;
Sigma = fan X;
P = polytope Sigma;
d = dim P;
F = faces_d P;
m = product(apply(numgens S, i -> S_i));
G = apply(max X, l -> m//product(apply(l,i -> S_i)));
H = hashTable apply(#G, i -> (j := F#i#0#0;((vertices P)_j,G_i)))
C = cellComplex(S,P,Labels => H);
Cres = (chainComplex C)[-1]
assert(betti (res B) == betti Cres);
assert(HH_0 Cres == S^1/B)
assert(HH_1 Cres == 0)
assert(HH_2 Cres == 0)
assert(HH_3 Cres == 0)
for i in {1,2,3} do assert(HH_i Cres == 0)
