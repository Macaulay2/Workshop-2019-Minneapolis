restart
load "nops.m2"

R = QQ[x]
I = ideal(x)
S = R/I
M = matrix{{x,1,0},{x^2,2*x,2},{x^3,3*x^2,6*x}}
kernel M

use R
I = ideal(x^2)
S = R/(radical I)
M = matrix{{x^2,2*x,2},{x^3,3*x^2,6*x},{x^4,4*x^3,12*x^2}}
kernel M


----------------

R = QQ[getSymbol "x",getSymbol "y"]
I = ideal((x^2+y)^2, x^3)
fi = gens I
bx = basis(0,4,R)
flatten (transpose fi * bx)
bd = basis(0,4,R)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))
S = R/radical(I)
M' = sub(M,S)
K = gens kernel M'

R' = makeWA R
dvars = (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R')
bdd = sub(bd, dvars)
bdd * sub(K, R')
first entries oo / (i -> applyNOp(i, first first entries gens I))
oo / (i -> i % radical I)

NoethOps(I)


-- Compute annihilators using Macaulay Matrix.
-- Consider differential monomials up to degree nd
-- and multiply generators with monomials up to degree nx
MacaulayMatrix = (nx, nd, I) -> (
	R := ring I;
	fi := gens I;
	bx := basis(0,nx,R);
	bd := basis(0,nd,R);

	-- macaulay matrix
	M := transpose diff(transpose bd, flatten (transpose fi*bx));
	S := R/radical(I);
	M' := sub(M,S);
	K := gens kernel M';

	-- Return elements in WeylAlgebra for nice formatting
	R' := makeWA R;
	dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
	bdd := sub(bd, dvars);
	use R;
	flatten entries (bdd * sub(K, R'))
)

MacaulayMatrixPD = (var, nx, nd, I) -> (
	R := ring I;
	fi := gens I;
	bx := basis(0,nx,R, Variables => var);
	bd := basis(0,nd,R, Variables => var);

	-- macaulay matrix
	M := transpose diff(transpose bd, flatten (transpose fi*bx));

	S := R/radical(I);
	M' := sub(M,S);
	K := gens kernel M';

	-- Return elements in WeylAlgebra for nice formatting
	R' := makeWA R;
	dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
	bdd := sub(bd, dvars);
	use R;
	flatten entries (bdd * sub(K, R'))
)



-- Tests
-- Tests if output of NoethOps is same as output of MacaulayMatrix
testMM = (nop, mm) -> (
	sub(nop, ring mm) // mm
)





R = QQ[x,y];
I = ideal(x^2-y ,y^2)
benchmark "NoethOps I"
benchmark "MacaulayMatrix(10,5,I)"

R = QQ[x,y];
I = ideal(x^2-y ,y^2);
sort NoethOps(I)
sort MacaulayMatrix(20,5,I)
apply(oo,ooo,testMM)

R = QQ[x,y];
I = ideal(x^4-x*y-y, y^2);
sort NoethOps(I)
sort MacaulayMatrix( 10, 8, I)
apply(oo,ooo,testMM)

R = QQ[x,y];
I = ideal(x^3-y, y^3);
sort NoethOps(I)
sort MacaulayMatrix( 10, 8, I)
apply(oo,ooo,testMM)

R = QQ[x,y,z];
I = ideal(x^2-z, y^2-z, z^2);
sort NoethOps(I)
sort MacaulayMatrix( 5, 5, I)
apply(oo,ooo,testMM)

-- pos dim

R = QQ[x,y,t];
I = ideal(x^2-t*y ,y^2);
sort NoethOps(I)
sort MacaulayMatrixPD({x,y}, 10, 10, I)
apply(oo,ooo,testMM)

R = QQ[x,y,t,s];
I = ideal(x^4-t*x*y-s*y, y^2);
sort NoethOps(I)
sort MacaulayMatrixPD({x,y}, 10, 10,  I)
apply(oo,ooo,testMM)

R = QQ[x,y,z,s,t];
I = ideal(x^2-t*z, y^2-s*z, z^2);
sort NoethOps(I)
sort MacaulayMatrixPD({x,y,z}, 3, 4, I)
apply(oo,ooo,testMM)



-- pos dim, not noether normalized
-- R = QQ[x,y,t];
-- I = ideal(x^2-t*y ,y^2);
-- netList sort NoethOps(I)
-- netList sort MacaulayMatrixPD({x,y}, 10, 10, I)
-- apply(oo,ooo,testMM)


-- zero dim primary, not centered
R = QQ[x,y]
I = ideal((x+1)^3,(y+1)^2)
MacaulayMatrix(10,6,I)


R = QQ[x,t]
I = ideal((x^2 - t)^2)
MacaulayMatrixPD({t}, 10,5, I)





-- pos dim primary, not centered
R = QQ[x,t]
I = ideal((x^2 - t)^2)
MacaulayMatrixPD({x},4,4,I)

-- zero dim, not primary
R = QQ[x,y]
I = intersect( ideal(x,y), ideal(x-1,y-1))
MacaulayMatrix(10,5,I)


R = QQ[x,y]
I = intersect( ideal(x^2,y^2), ideal(x-1,y-1))
MacaulayMatrix(3,3,I)

J1 = ideal(x^2,y^2)
J2 = ideal(x-1, y-1)
MacaulayMatrix(10,5,J1)