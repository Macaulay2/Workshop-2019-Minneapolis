needsPackage "NoetherNormalization"
needsPackage "Dmodules"


-- Compute Noetherian operators.
-- Caveat: variety of the extended ideal must be the origin!
brokenNoethOps = {Normalized => false, Verbose => false} >> opts -> J -> (
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
	b := first entries basis(0,mu-1,R, Variables => var);
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

-- given a map f, and NOps for f(I), get NOps for I
convertNOps = (nn,f) -> (
	fMatrix := mapToMatrix(f);
	W := ring (nn#0);
	nVars := (#(gens W))//2;
	zeroMatrix := map(QQ^nVars, nVars, 0);
	f'Matrix := mapToMatrix(inverse f);
	rules := matrix({{f'Matrix, zeroMatrix},{zeroMatrix,transpose fMatrix}}) * (transpose vars W);
	rules = flatten entries rules;
	phi := map(W,W,rules);

	-- return the results
	print"return the results";
	nn / phi
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

-- Try to see if gens of I applied with all Noeth Ops
-- vanish on rad(I)
visualCheck = (nops, I) -> (
	foo := table(nops, I_*, (n,i) -> applyNOp(n,i));
	netList applyTable(foo, (i -> i%(radical I)))
)

sanityCheck = (nops, I) -> (
	foo := flatten table(nops, I_*, (n,i) -> applyNOp(n,i));
	all(foo, i -> i%(radical I) == 0)
)


noethOps = method(Options => {DegreeLimit => 5}) 
noethOps (Ideal) := List => opts -> (I) -> (
	R := ring I;
	var := gens R - set support first independentSets I;
	bx := flatten entries basis(0,opts.DegreeLimit,R, Variables => gens R);
	bd := basis(0,opts.DegreeLimit,R, Variables => var);

	elapsedTime M := diff(bd, transpose matrix {flatten (table(bx,I_*,(i,j) -> i*j))});
	M' := M % ideal gens radical(I);
	K := gens trim kernel M';

	-- Return elements in WeylAlgebra for nice formatting
	R' := makeWA (R, SetVariables => false);
	dvars := (options R').WeylAlgebra / (i -> (i#0)_R => (i#1)_R');
	bdd := sub(bd, dvars);
	flatten entries (bdd * sub(K, R'))
)


socleMonomials = method()
socleMonomials(Ideal) := List => (I) -> (
	R := ring I;
	if codim I < dim R then error("Expected Artinian ideal");
	inI := monomialIdeal leadTerm gens gb I;
	irreducibleDecomposition inI / gens / 
						entries / first / 
						product / exponents / first / 
						(i -> (i / (j -> j-1))) / (e -> R_e)
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


TEST ///
R = QQ[x,y,z]
I = ideal(x^2 - y, y^2)
nops = MacaulayMatrixPD({x,y}, 10, 10, I)
assert(sanityCheck(nops, I))
///

TEST ///
R = QQ[x_0..x_3]
S = QQ[s,t]
I0 = ker map(S,R,{s^5,s^3*t^2, s^2*t^3, t^5})
depvars = gens R - set support first independentSets I0
nops = MacaulayMatrixPD(depvars, 10, 10, I0)
assert(sanityCheck(nops, I0))
I1 = ideal(x_0^2, x_1^2, x_2^2)
depvars = gens R - set support first independentSets I1
nops = MacaulayMatrixPD(depvars, 10, 10, I1)
///

end--

load "nops.m2"



--TODO:

-- 1) Figure out bounds for MacaulayMatrix/MacaulayMatrixPD
