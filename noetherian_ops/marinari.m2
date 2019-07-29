restart


R = QQ[x,y,z]
I = ideal(x^2-z, y^2-z, z^10)
t_1 = cpuTime();
a = NoethOps0Slow I;
t_2 = cpuTime();
print("Slow: " | toString(t_2-t_1))
b = NoethOps0Fast I;
print("Fast: " | toString(cpuTime() - t_2))



R = QQ[x,y]
I = ideal (x^2-y ,y^2)
t_1 = cpuTime();
a = NoethOps0Slow I;
t_2 = cpuTime();
print("Slow: " | toString(t_2-t_1))
b = NoethOps0Fast I;
print("Fast: " | toString(cpuTime() - t_2))
set a === set b


R = QQ[x,y]
I = ideal (x^4-xy-y, y^2)
t_1 = cpuTime();
a = NoethOps0Slow I;
t_2 = cpuTime();
print("Slow: " | toString(t_2-t_1))
b = NoethOps0Fast I;
print("Fast: " | toString(cpuTime() - t_2))
set a === set b

R = QQ[x,y]
I = ideal (x^3-y, y^3)
t_1 = cpuTime();
a = NoethOps0Slow I;
t_2 = cpuTime();
print("Slow: " | toString(t_2-t_1))
b = NoethOps0Fast I;
print("Fast: " | toString(cpuTime() - t_2))
set a === set b

R = QQ[x,y,z]
I = ideal (x^2-z, y^2-z, z^2)
t_1 = cpuTime();
a = NoethOps0Slow I;
t_2 = cpuTime();
print("Slow: " | toString(t_2-t_1))
b = NoethOps0Fast I;
print("Fast: " | toString(cpuTime() - t_2))
set a === set b




---- Positive dimension -----

restart
load "nops.m2"


R = QQ[x,y];
I = ideal"x2-y ,y2";
benchmark "NoethOps(I)"

R = QQ[x,y];
I = ideal"x4-xy-y, y2";
benchmark "NoethOps(I)"

R = QQ[x,y];
I = ideal"x3-y, y3";
benchmark "NoethOps(I)"

R = QQ[x,y,z];
I = ideal"x2-z, y2-z, z2";
benchmark "NoethOps(I)"

R = QQ[x,y,t];
I = ideal"x2-ty ,y2";
benchmark "NoethOps(I)"

R = QQ[x,y,t,s];
I = ideal"x4-txy-sy, y2";
benchmark "NoethOps(I)"

R = QQ[x,y,z,s,t];
I = ideal"x2-tz, y2-sz, z2+x";
time NoethOps(I)
sanityCheck(oo,I)




R = QQ[x,y,t]
J = ideal"x3,y2,-xt+y3, y6t"

I = J;

nn = NoethOps (I, Verbose=>true)
sanityCheck(oo, I)

(a,b,c) = noetherNormalization I
NoethOps(b, Normalized => true)
sanityCheck(oo, b)
convertNOps(ooo, a)
sanityCheck(oo, I)


-- fails
restart
load "nops.m2"
R = QQ[x_(0,0)..x_(2,2)]
M = genericMatrix(R,3,3)
I = minors(2,M)

time NoethOps(I, Verbose =>true)

-- R = QQ[x,y,z,s,t]
-- J = ideal(x_1*x_2*x_3-x_4,x_3^2*x_4)
-- J = ideal(x^2,(y+1)^2,-x*t+(y+1))
-- S = QQ[x,y,t]
-- J = ideal(x^2 - t*y,y^2)
-- J = ideal (x^4-t*x*y-s*y, y^2)
-- J = ideal(x^2-t*z, y^2-s*z, z^2)
-- S = QQ[w..z]
-- J = monomialCurveIdeal(S,{1,2,3})
-- S = QQ[x,y]
-- I = ideal(x^2-y,y^2)
-- phi = map(S,S,{x+1,y})
-- J = phi(I)
-- S = QQ[x,y,z,w]
-- J = ideal(x^2,z^2,w^2)
-- R = QQ[x,y,z,t,s]
-- J = ideal"x2-tz, y2-sz, z2"

	startTime := cpuTime();

	S := ring J;

	-- Noether normalize, record parameters,
	-- variables & automorphism
	(f, I, par) := noetherNormalization J;
	var := first entries vars ring J - set par;
	mu := degree I;

	-- Define new ring with correct monomial order
	R := QQ[var, par, MonomialOrder => {#var, #par}];
	I = sub(I,R);
	var = var / (i -> sub(i,R));
	par = par / (i -> sub(i,R));


	-- Saturate wrt I (?)
	G := first entries gens gb I;
	tg := G / (i->leadMonomial(i)) / (i -> first entries monomials(i, Variables => par));
	tg = product flatten tg;
	print tg;


	-- Multiply all by t^gamma
	b := first entries basis(0,mu-1,R, Variables => var);
	bt := b;

	-- TODO add a loop here
	bt = bt / (i->i*tg);
	nf := bt / (i->i % I);

	-- check number of coeffs
	beta = set flatten (nf / (i -> first entries monomials(i, Variables => var)));
	if not mu == #beta then error("Need to saturate more... Fix");


	-- kill zeros
	idx := positions(nf, j -> j!=0);
	b = b_idx;
	nf = nf_idx;

	-- record results in a hash table. Keys are standard monomials
	ht := new MutableHashTable from (applyPairs(beta, (i,j) -> i=> 0));

	-- populate the hash table
	for i from 0 to #idx - 1 do (
		print i;
		(mon,coe) = coefficients (nf#i, Variables => var);
		mon = first entries mon;
		coe = flatten entries coe;
		for j from 0 to #mon - 1 do(
			ht#(mon#j) = ht#(mon#j) + (coe#j) * (b#i);
		);
	);
	-- check number of monomials
	res := values ht;

	-- divide by factorials
	res = res / (i -> addFactorials(i, var));

	-- replace x by dx
	W := makeWA (QQ[var,par]);
	WARules := (options(W))#WeylAlgebra / (i -> (i#0)_W => (i#1)_W);
	diffRules := select (WARules, i -> member(sub((i#0),R), var));
	res = res / (i -> sub(i, W)) / (i -> sub(i, diffRules));

	-- Convert to NOps of previous ideal
	fMatrix := mapToMatrix(f);
	nVars := #par + #var;
	zeroMatrix := map(QQ^nVars, nVars, 0);
	f'Matrix := mapToMatrix(inverse f);
	rules := matrix({{f'Matrix, zeroMatrix},{zeroMatrix,transpose fMatrix}}) * (transpose vars W);
	rules = flatten entries rules;
	phi := map(W,W,rules);

	-- return the results
	res / phi

endTime = cpuTime() - startTime



(
NoethOps0Fast = I -> (
	R = ring I;

	error("Broken, don't use");
	-- Sanity check
	if not isPrimary(I) then (
		error "Ideal not primary!";
	);

	if not dim I == 0 then (
		error "Ideal not zero-dimensional!";
	);

	if not minimalPrimes I == {ideal vars R} then (
		error "Ideal does not vanish in origin!";
	);


	mu = degree(I);
	rules = flatten for i from 0 to mu-1 list (
		b = first entries basis(i,R);
		nf = b / (j -> sub(j,R/I));
		-- Remove zeros
		idx = positions(nf, j -> j!=0);
		if #idx == 0 then break;
		b = b_idx;
		nf = nf_idx;
		apply(b,nf, (i,j) -> sub(j,R) => i)
	);

	ht = hashTable(plus, rules);
	values ht
))


(
NoethOps0Slow = I -> (
	R := ring I;

	-- Sanity check
	if not isPrimary(I) then (
		error "Ideal not primary!";
	);

	if not dim I == 0 then (
		error "Ideal not zero-dimensional!";
	);

	if not minimalPrimes I == {ideal vars R} then (
		error "Ideal does not vanish in origin!";
	);

	mu := degree(I);
	b := first entries basis(0,mu-1,R);
	nf := b / (i -> sub(i,R/I));

	-- Remove zeros
	idx := positions(nf, i->i!=0);
	b = b_idx;
	nf = nf_idx;

	beta := set flatten (nf / (i -> first entries monomials(i)));

	ht := new MutableHashTable from (applyPairs(beta, (i,j) -> i=> 0));

	for i from 0 to #idx - 1 do (
		print i;
		(mon,coe) := coefficients (nf#i);
		mon = first entries mon;
		coe = flatten entries coe;
		for j from 0 to #mon - 1 do(
			ht#(mon#j) = ht#(mon#j) + sub(coe#j, R) * (sub(b#i,R));
		);
	);
	values ht
))

(
NoethOps = J -> (
	startTime := cpuTime();

	S := ring J;

	-- Noether normalize, record parameters,
	-- variables & automorphism
	(f, I, par) := noetherNormalization J;
	var := first entries vars ring J - set par:
	mu := degree I;

	-- Define new ring with correct monomial order
	R := QQ[var, par, MonomialOrder => {#var, #par}];
	I := sub(I,R);
	var := var / (i -> sub(i,R));
	par := par / (i -> sub(i,R));


	-- Saturate wrt I (?)
	G := first entries gens gb I;
	tg := G / (i->leadMonomial(i)) / (i -> first entries monomials(i, Variables => par));
	tg = product flatten tg;


	-- Multiply all by t^gamma
	b := first entries basis(0,mu-1,R, Variables => var);
	bt := b;

	-- TODO add a loop here
	bt = bt / (i->i*tg);
	nf := bt / (i->i % I);

	-- check number of coeffs
	beta = set flatten (nf / (i -> first entries monomials(i, Variables => var)));
	if not mu == #beta then error("Need to saturate more... Fix");


	-- kill zeros
	idx := positions(nf, j -> j!=0);
	b = b_idx;
	nf = nf_idx;

	-- record results in a hash table. Keys are standard monomials
	ht := new MutableHashTable from (applyPairs(beta, (i,j) -> i=> 0));

	-- populate the hash table
	for i from 0 to #idx - 1 do (
		print i;
		(mon,coe) = coefficients (nf#i, Variables => var);
		mon = first entries mon;
		coe = flatten entries coe;
		for j from 0 to #mon - 1 do(
			ht#(mon#j) = ht#(mon#j) + (coe#j) * (b#i);
		);
	);
	-- check number of monomials
	res := values ht;

	-- divide by factorials
	res = res / (i -> addFactorials(i, var));

	-- replace x by dx
	W := makeWA (QQ[var,par]);
	WARules := (options(W))#WeylAlgebra / (i -> (i#0)_W => (i#1)_W);
	diffRules := select (WARules, i -> member(sub((i#0),R), var));
	res = res / (i -> sub(i, W)) / (i -> sub(i, diffRules));

	-- Convert to NOps of previous ideal
	fMatrix := mapToMatrix(f);
	nVars := #par + #var;
	zeroMatrix := map(QQ^nVars, nVars, 0);
	f'Matrix := mapToMatrix(inverse f);
	rules := matrix({{f'Matrix, zeroMatrix},{zeroMatrix,transpose fMatrix}}) * (transpose vars W);
	rules = flatten entries rules;
	phi := map(W,W,rules);

	-- return the results
	res / phi
)
)

(

-- f to Matrix
mapToMatrix = f -> (
	S := f.source;
	v := first entries vars S;
	M := first entries f.matrix;
	matrix table(M,v, (i,j) -> coefficient(j,i))
)


matrixToMap = (M, R) -> (
	f = map(R,R,transpose(M*(transpose vars R)))
)

addFactorials = (f, var) -> (
	(coe, mon) := coefficients(f, Variables => var);
	factorials := first entries coe / (i -> first exponents(i)) / (i -> product( i / (j -> j!)));
	res := apply(flatten entries mon, factorials, (i,j) -> i/j);
	sum apply(first entries coe, res, (i,j) -> i*j)
)


-- Given an element N in Weyl algebra and a polynomial
-- f, compute the result of applying N to f.
applyNOp = (N, f) -> (
	WW := ring N;
	diffVars := first entries vars WW;
	n := #diffVars;
	DOs := (options(WW))#WeylAlgebra / (i -> i#1) / (i-> i_WW);
	(mon,coe) = coefficients(N, Variables => DOs);
	moncoe := apply(flatten entries mon, flatten entries coe, (i,j) -> (i,j));
	sum for i in moncoe list (
		exps = (first exponents i#0)_{n//2..(n-1)};
		deriv = product(apply(exps, diffVars_{0..n//2-1}, (i,j) -> j^i));
		diff(sub(deriv, (ring f)), f)*sub((i#1), ring f)
	)
)
)
