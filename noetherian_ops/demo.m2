restart
load "nops.m2"
-- 0 dim

R = QQ[x,y];
I = ideal(x^2,y)
NoethOps(I)

R = QQ[x,y];
I = ideal(x^2-y ,y^2);
NoethOps(I)

R = QQ[x,y];
I = ideal(x^4-x*y-y, y^2);
NoethOps(I)

R = QQ[x,y];
I = ideal(x^3-y, y^3);
NoethOps(I)

R = QQ[x,y,z];
I = ideal(x^2-z, y^2-z, z^2);
NoethOps(I)

-- pos dim

R = QQ[x,y,t];
I = ideal(x^2-t*y ,y^2);
NoethOps(I)

R = QQ[x,y,t,s];
I = ideal(x^4-t*x*y-s*y, y^2);
NoethOps(I)

R = QQ[x,y,z,s,t];
I = ideal(x^2-t*z, y^2-s*z, z^2);
NoethOps(I)


-- Ordering of variables

R = QQ[x,y,z]
I = ideal(x*z-y,y^2,z^2)
NoethOps(I)
sanityCheck(oo,I)

R = QQ[z,y,x]
J = sub(I, R)
NoethOps(J)

-- Centering
R = QQ[x,y,z,s,t];
I = ideal(x^2-t*z, y^2-s*z, z^2+x);
NoethOps(I)
sanityCheck(oo,I)

R = QQ[x,y]
I = ideal(x^2 + t)
radical(I) == I
NoethOps(I)
sanityCheck(oo, I)


R = QQ[x,y]
I = ideal"x-1,y2"
dim I
degree I
isPrimary I
primaryDecomposition I

NoethOps I
NoethOpsExtra(I,4)


NoethOpsExtra = {Normalized => false, Verbose => false} >> opts -> (J,m) -> (
	S := ring J;

	-- Noether normalize, record parameters,
	-- variables & automorphism
	(f,I,par) := if not opts.Normalized then (
		if opts.Verbose then print "Performing Noether normalization";
		noetherNormalization J
	) else (
		if opts.Verbose then print("Assuming in normal position wrt to " | toString((flatten entries vars S)_{-(dim J)..-1}));
		(map(S,S), J, (flatten entries vars S)_{-(dim J)..-1})
	);
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
	if opts.Verbose then print("t^gamma = " | toString(tg));


	-- Multiply all by t^gamma
	b := first entries basis(0,m-1,R, Variables => var);
	bt := b;

	
	-- Normal form of Taylor polynomial
	bt = bt / (i->i*tg);
	nf := bt / (i->i % I);
	
	-- check number of coeffs
	beta := set flatten (nf / (i -> first entries monomials(i, Variables => var)));
	while not mu == #beta do (
		if opts.Verbose then print("Multiply GB by " | toString(tg));
		nf = nf / (i->i*tg) / (i -> i % I);
		beta = set flatten (nf / (i -> first entries monomials(i, Variables => var)));
	);


	-- kill zeros
	idx := positions(nf, j -> j!=0);
	b = b_idx;
	nf = nf_idx;

	-- record results in a hash table. Keys are standard monomials
	ht := new MutableHashTable from (applyPairs(beta, (i,j) -> i=> 0));

	-- populate the hash table
	for i from 0 to #idx - 1 do (
		if opts.Verbose then print("Processing term " | toString(i+1));
		(mon,coe) = coefficients (nf#i, Variables => var);
		mon = first entries mon;
		coe = flatten entries coe;
		for j from 0 to #mon - 1 do(
			ht#(mon#j) = ht#(mon#j) + (coe#j) * (b#i);
		);
	);
	-- check number of monomials
	if opts.Verbose then print"check number of monomials";
	res := values ht;

	-- divide by factorials
	if opts.Verbose then print"divide by factorials";
	res = res / (i -> addFactorials(i, var));

	-- replace x by dx
	if opts.Verbose then print"replace x by dx";
	W := makeWA (QQ[var,par]);
	WARules := (options(W))#WeylAlgebra / (i -> (i#0)_W => (i#1)_W);
	diffRules := select (WARules, i -> member(sub((i#0),R), var));
	res = res / (i -> sub(i, W)) / (i -> sub(i, diffRules));

	-- Convert to NOps of previous ideal
	if opts.Verbose then print"Convert to NOps of previous ideal";
	fMatrix := mapToMatrix(f);
	nVars := #par + #var;
	zeroMatrix := map(QQ^nVars, nVars, 0);
	f'Matrix := mapToMatrix(inverse f);
	rules := matrix({{f'Matrix, zeroMatrix},{zeroMatrix,transpose fMatrix}}) * (transpose vars W);
	rules = flatten entries rules;
	phi := map(W,W,rules);

	-- return the results
	if opts.Verbose then print"return the results";
	use S;
	res / phi
)

