restart
load "nops.m2"

----------------

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
	R' := QQ[toList(set gens R - set var)];

	print"macaulay matrix";
	-- macaulay matrix
	M := transpose diff(transpose bd, flatten (transpose fi*bx));

	print "compute radical";
	S := R/radical(I);
	print "sub";
	M' := sub(M,S);
	print "kernel";
	K := gens kernel M';

	-- Return elements in WeylAlgebra for nice formatting
	print "postprocess";
	R' := makeWA R;
	dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
	bdd := sub(bd, dvars);
	use R;
	print "return";
	flatten entries (bdd * sub(K, R'))
)



-- Tests
-- Tests if output of NoethOps is same as output of MacaulayMatrix
testMM = (nop, mm) -> (
	(sub(nop, ring mm) % mm) == 0
)

-- Divide by gcd of coefficients
modConstant = f -> (
	(mon, coe) := coefficients f;
	g := (flatten entries coe) / (i -> sub(i,QQ));
	f // gcd(g)
)


-- Zero dim, centered, primary
R = QQ[x,y];
I = ideal(x^2-y ,y^2);
nops = sort NoethOps(I)
mm = sort MacaulayMatrix(20,5,I)
apply(nops,mm,testMM)

sort MacaulayMatrix(0,0,I)
sort MacaulayMatrix(1,1,I)
sort MacaulayMatrix(2,2,I)
sort MacaulayMatrix(3,3,I)
sort MacaulayMatrix(4,4,I)
sort MacaulayMatrix(5,5,I)

sort MacaulayMatrix(6,10,I)

sort MacaulayMatrix(10,1,I)
sort MacaulayMatrix(10,2,I)
sort MacaulayMatrix(10,3,I)
sort MacaulayMatrix(10,4,I)
sort MacaulayMatrix(10,5,I)



R = QQ[x,y];
I = ideal(x^2-y ,y^2);
nops = sort NoethOps(I)
mm = sort MacaulayMatrix(20,3,I)
apply(nops,mm,testMM)

mm = sort MacaulayMatrix(5,3,I)
mm = sort MacaulayMatrix(2,3,I)
mm = sort MacaulayMatrix(1,3,I)
mm = sort MacaulayMatrix(0,3,I)




R = QQ[x,y,z];
I = ideal(x^2-z, y^2-z, z^2);
nops = sort NoethOps(I)
mm = sort MacaulayMatrix( 10, 4, I)
apply(nops,mm,testMM)
#mm


mm = sort MacaulayMatrix( 10, 4, I)
mm = sort MacaulayMatrix( 5, 4, I)
mm = sort MacaulayMatrix( 2, 4, I)
mm = sort MacaulayMatrix( 1, 4, I)
#oo





-- zero dim, not centered
R = QQ[x,y]
Ic = ideal(x^2-y ,y^2);
I = sub(Ic, {x => x-2, y => y+1})

nops = sort NoethOps(Ic)
mm = sort MacaulayMatrix(10,5,I)
mm / modConstant

mm = sort MacaulayMatrix(10,5,I)
mm = sort MacaulayMatrix(5,5,I)
mm / modConstant

mm = sort MacaulayMatrix(3,3,I)
tex oo
mm = sort MacaulayMatrix(2,3,I)
tex oo
mm = sort MacaulayMatrix(1,3,I)
tex oo
mm = sort MacaulayMatrix(0,3,I)
tex oo




-- pos dim, centered

R = QQ[x,y,t];
I = ideal(x^2-t*y ,y^2);
nop = sort NoethOps(I)
mm = sort MacaulayMatrixPD({x,y}, 10, 10, I)
apply(mm,nop,testMM)

mm = sort MacaulayMatrixPD({x,y}, 10, 3, I)
mm = sort MacaulayMatrixPD({x,y}, 0, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 1, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 2, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 3, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 4, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 5, 5, I)
mm = sort MacaulayMatrixPD({x,y}, 6, 5, I)



R = QQ[x,y,z,s,t];
I = ideal(x^2-t*z, y^2-s*z, z^2);
sort NoethOps(I)
sort MacaulayMatrixPD({x,y,z}, 3, 4, I)
apply(oo,ooo,testMM)



-- pos dim primary, not centered
R = QQ[x,t]
I = ideal(x^2 - t)
mm = MacaulayMatrixPD({x},5,5,I)

I = ideal((x^2-t)^2)
mm = MacaulayMatrixPD({x},5,5,I)

R = QQ[x,y,t];
Ic = ideal(x^2-t*y ,y^2);
I = sub(Ic, {x=>x+t, y => y-t})
mm = MacaulayMatrixPD({x,y}, 5,5, I)
oo / modConstant

MacaulayMatrixPD({x,y}, 5,5, Ic) / modConstant




-- zero dim, not primary
R = QQ[x,y]
I = intersect( ideal(x^2,y), ideal(x-1,y-1))
MacaulayMatrix(10,5,I) / modConstant
MacaulayMatrix(5,5,I) / modConstant
MacaulayMatrix(3,3,I) / modConstant





-- possible numerical approach

-- zero dim:
R = QQ[x,y]
I = intersect( ideal(x^2,y), ideal(x-1,(y-1)^2))
nx = 2;
nd = 2;
bx := basis(0,nx,R)
bd := basis(0,nd,R)
fi = gens I

M = transpose diff(transpose bd, flatten (transpose fi*bx))
-- Evaluate M at (0,0)
M' = sub(M, {x=>0,y=>0})
K = gens kernel M'
bd * K -- corresponds to {1, dx}

-- Evaluate M at (1,1)
M' = sub(M, {x=>1, y=>1})
K = gens kernel M'
bd * K -- corresponds to {1, dy}




-- positive dim
R = QQ[x,y,t];
Ic = ideal(x^2-t*y ,y^2);
I = sub(Ic, {x => x+t, y => y-t})

-- Points on the variety: (-1,1,1), (-2,2,2), (-3,3,3)
nx = 3;
nd = 3;
bx = basis(0,nx,R)
bd = basis(0,nd,R, Variables => {x,y})
fi = gens I
M = transpose diff(transpose bd, flatten (transpose fi*bx))


-- Sub point (-1,1,1)
M' = sub(M, {x => -1, y => 1, t => 1})
K = gens kernel M'
k1 = flatten entries (bd * K) -- corresponds to {1, dy}

-- Sub point (-2,2,2)
M' = sub(M, {x => -2, y => 2, t => 2})
K = gens kernel M'
k2 = flatten entries (bd * K) -- corresponds to {1, dy}

-- Sub point (-3,3,3)
M' = sub(M, {x => -3, y => 3, t => 3})
K = gens kernel M'
k3 = flatten entries (bd * K) -- corresponds to {1, dy}


netList({{"t = 1"} | k1, {"t = 2"} | k2, {"t = 3"} | k3})

-- Interpolate and conclude
-- {1, dx, 1/6 t dx^3 + dx*dy, 1/2 t dx^2 + dy}

NoethOps Ic -- correct answer




-- Sub point (-1,1,1)
M' = sub(M, {x => -1, y => 1})
K = gens kernel M
k1 = flatten entries (bd * K) -- corresponds to {1, dy}

-- Sub point (-2,2,2)
M' = sub(M, {x => -2, y => 2})
K = gens kernel M'
k2 = flatten entries (bd * K) -- corresponds to {1, dy}

-- Sub point (-3,3,3)
M' = sub(M, {x => -3, y => 3})
K = gens kernel M'
k3 = flatten entries (bd * K) -- corresponds to {1, dy}








-- DEBUG


-- Non noether normalized, positive dimension
R = QQ[x,y,z,w]
J = monomialCurveIdeal(R, {1,2,3})
(f,I,t) = noetherNormalization J

var = {x,y}

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)
R' = QQ[toList(set gens R - set var)]

print"macaulay matrix"
-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))

-- pts {x,y,z,w}
p = {0,0,0,1} / (i->i_QQ);
f' = matrixToMap(inverse mapToMatrix(f), R)
q = flatten entries sub(matrix f', matrix {p})
g = map(QQ,R, matrix {q});
K = gens kernel g M;
R' = makeWA R;
dvars = (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
bdd = sub(bd, dvars);
flatten entries (bdd * sub(K, R'))
pointToKernel q
#oo

p = {0,1,0,0};
q = flatten entries sub(matrix f, matrix {p})
pointToKernel q
#oo

p = {1,0,0,0};
q = flatten entries sub(matrix f, matrix {p})
pointToKernel q
#oo

p = {1,1,1,1};
q = flatten entries sub(matrix f, matrix {p})
pointToKernel q
#oo

p = {4, -1, 1/4, -1/16}
q = flatten entries sub(matrix f, matrix {p})
pointToKernel q



pointToKernel = p -> (
	f := map(QQ,R, matrix {p});
	K := gens kernel f M;
	R' := makeWA R;
	dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
	bdd := sub(bd, dvars);
	flatten entries (bdd * sub(K, R'))
)





print "compute radical"
S := R/radical(I)
print "sub"
M' := sub(M,S)





print "kernel"
K := gens kernel M'

-- Return elements in WeylAlgebra for nice formatting
print "postprocess"
R' := makeWA R
dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R')
bdd := sub(bd, dvars)
use R
print "return"
flatten entries (bdd * sub(K, R'))











needsPackage "NoetherNormalization"
-- S(2,1)
R = QQ[x0,x1,x2,y0,y1]
M = matrix{{x0,x1,y0},{x1,x2,y1}}
I = minors(2,M)
(f,J,t) = noetherNormalization I

primaryDecomposition I
MacaulayMatrixPD({x0,x1}, 2,3, J)

-- Points on the variety:
--  (x0, x1, x2, y0, y1)
--  sub(matrix f, matrix{{ 0,  0,  1,  0,  1}})
--  sub(matrix f, matrix{{ 0,  0,  0,  0,  0}})
--  sub(matrix f, matrix{{ 1,  0,  0,  1,  0}})
--  sub(matrix f, matrix{{ 1,  1,  1,  1,  1}})


matrix {{1, 0, 0, 1, 0}}












-- ennen 10

MacaulayMatrixPD({x0,x1,x2},3,4,J) / modConstant


-- X(2,2)
R = QQ[x_0..x_2,y_0..y_2]
I = ideal(det matrix{{x_0,x_1},{x_1,x_2}}, det matrix{{y_0,y_1},{y_1,y_2}}, det matrix{{x_0 + y_0, x_1 + y_1}, {x_1 + y_1, x_2 + y_2}})
(f,J,par) = noetherNormalization I
MacaulayMatrixPD({x_0,x_1,x_2}, 2,3,J)