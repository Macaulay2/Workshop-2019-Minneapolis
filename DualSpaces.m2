-- -*- coding: utf-8 -*-
newPackage(
     "DualSpaces",
     PackageExports => {"NAGtypes"},
     Version => "0.1", 
     Date => "Aug 1, 2019",
     Authors => {
        {Name => "Robert Krone", 
    	       Email => "krone@math.gatech.edu"},
        {Name => "Justin Chen", 
               Email => "jchen646@gatech.edu"},
        {Name => "Marc Harkonen", 
               Email => "harkonen@gatech.edu"}
     },
     Headline => "numerically compute local dual spaces, Hilbert functions, and Noetherian operators",
     PackageImports => {"Truncations"},
     AuxiliaryFiles => false
)

export {
     "truncatedDual",
     "zeroDimensionalDual",
     "gCorners",
     "socles",
     "localHilbertRegularity",
     "eliminatingDual",
     "innerProduct",
     "reduceSpace",
     "orthogonalInSubspace",
     "DZ",
     "BM",
     "Normalize",
     "ProduceSB",
     "numericalKernel",
     "numericalImage",
     "colReduce",
     "adjointMatrix",

     "noethOps",
     "DependentSet",
     "numNoethOps",
     "noethOpsFromComponents",
     "coordinateChangeOps",
     "sanityCheck",
     "applyNOp"
     }

--TruncDualData private keys
protect \ {
    Igens,
    syl,
    strategy,
    deg,
    dBasis,
    hIgens,
    BMintegrals,
    BMcoefs,
    BMbasis,
    Seeds
    }

-----------------------------------------------------------------------------------------

-- Default tolerance value respectively for exact fields and inexact fields
defaultT := R -> if precision 1_R == infinity then 0 else 1e-6;

shiftGens := (p,Igens) -> (
    R := ring Igens;
    sub(Igens, matrix{gens R + apply(p.Coordinates,c->sub(c,R))})
    )
    

vanishingOps = method(Options => {Tolerance => null})
vanishingOps(PolySpace,Matrix,Boolean,Ideal) := o -> (polys,ops,normalize,P) -> vanishingOps(gens polys,ops,normalize,P,o)
vanishingOps(Matrix,Matrix,Boolean,Ideal) := o -> (polys,ops,normalize,P) -> (
    R := ring polys;
    t := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    M := if normalize then contract(ops,transpose polys) else diff(ops,transpose polys);
    M = sub(M,R/P);
    kern := gens trim kernel M;
    print(M,ops,polys,kern);
    sub(kern,R)
    )

listFactorial = L -> product(L, l->l!)

truncatedDual = method(TypicalValue => DualSpace, Options => {Strategy => BM, Tolerance => null, Normalize=>true})
truncatedDual(Point,Ideal,ZZ) := o -> (p,I,d) -> truncatedDual(p,gens I,d,o)
truncatedDual(Point,Matrix,ZZ) := o -> (p,Igens,d) -> (
    R := ring Igens;
    if d < 0 then return dualSpace(map(R^1,R^0,0),p);
    t := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    R' := diffAlg R;
    P := ideal ((gens R) - p.Coordinates);
    ops := basis(0,d,R);
    dBasis := vanishingOps(idealBasis(Igens,d,false),ops,o.Normalize,P);
    dBasis = sub(ops,vars R')*dBasis;
    dualSpace(dBasis,p)    
    --TDD := initializeDualData(Igens,false,t,Strategy=>o.Strategy);
    --TDD = nextTDD(d,TDD,t);
    --dualSpace(TDD,p)
    )


zeroDimensionalDual = method(TypicalValue => DualSpace, Options => {Strategy => BM, Tolerance => null, Normalize=>true})
zeroDimensionalDual(Point,Ideal) := o -> (p,I) -> zeroDimensionalDual(p,gens I,o)
zeroDimensionalDual(Ideal,Ideal) := o -> (P,I) -> zeroDimensionalDual(P,gens I,o)
zeroDimensionalDual(Point,Matrix) := o -> (p,Igens) -> (
    R := ring Igens;
    P := ideal ((gens R) - p.Coordinates);
    zeroDimensionalDual(P,Igens,o)
    )
zeroDimensionalDual(Ideal,Matrix) := o -> (P,Igens) -> (
    R := ring Igens;
    t := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    --TDD := initializeDualData(Igens,false,t,Strategy=>o.Strategy);
    dBasis := map(R^1,R^0,0);
    d := 0;
    dDim := -1;
    R' := diffAlg R;
    while numcols dBasis != dDim do (
	dDim = numcols dBasis;
	ops := basis(0,d,R);
	dBasis = vanishingOps(idealBasis(Igens,d,false),ops,o.Normalize,P);
	dBasis = sub(ops,vars R')*dBasis;
    	--TDD = nextTDD(d,TDD,t);
	--dBasis = polySpace TDD;
	d = d+1;
	);
    dualSpace(dBasis,origin(R))
    )

--An object that stores the data for an ongoing iterative tuncated dual space computation
TruncDualData = new Type of MutableHashTable
initializeDualData = method(Options => {Strategy => BM})
initializeDualData (Matrix,Boolean,Number) := o -> (Igens,syl,t) -> (
    H := new MutableHashTable;
    R := ring Igens;
    H.Igens = Igens;
    H.syl = syl;
    H.strategy = o.Strategy;
    H.deg = 0;
    h := symbol h;
    S := (coefficientRing R)(monoid[{h}|gens R, MonomialOrder => {Weights => (numgens R+1):1, 1, (options R).MonomialOrder}]); --projectivization of R
    h = S_0;
    H.hIgens = homogenize(sub(Igens,S), h); 
    T := if syl then S else R;
    H.Seeds = dualSpace(matrix{{1_T}},origin(T));
    H.BMmatrix = innerProduct(polySpace if syl then H.hIgens else H.Igens, H.Seeds);
    H.BMintegrals = gens H.Seeds;
    H.BMcoefs = numericalKernel(H.BMmatrix,t);
    H.BMbasis = H.BMcoefs;
    --print(H.BMmatrix,H.BMcoefs);
    H.dBasis = H.BMintegrals * H.BMcoefs;
    if H.syl then H.dBasis = (map(R,S,{1_R} | gens R)) H.dBasis;
    new TruncDualData from H
    )
-- advances the truncated dual computation from whatever was stored in parameters up to degree d
nextTDD = method()
nextTDD (TruncDualData,Number) := (H,t) -> nextTDD(H.deg + 1,H,t)
nextTDD (ZZ,TruncDualData,Number) := (d,H,t) -> (
    R := ring H.Igens;
    S := ring H.hIgens;
    if H.strategy == DZ then (
	Rd := polySpace basis(0,d,R);
	Id := idealBasis(H.Igens,d,H.syl);
	H.dBasis = gens orthogonalInSubspace(Id,Rd,t);
	); 
    if H.strategy == BM then (
	dehomog := map(R,S,{1_R} | gens R);
     	for e from H.deg+1 to d do (
	    (M,E) := BMmatrix H;
	    H.BMmatrix = M; H.BMintegrals = E;
	    H.BMcoefs = numericalKernel(M,t);
	    --print(M,H.BMcoefs);
	    I := basisIndices(last coefficients E*H.BMcoefs, t);
	    H.BMbasis = submatrix(H.BMcoefs, I);
	    H.dBasis = if H.syl then dehomog(E*H.BMcoefs) else H.dBasis | E*H.BMcoefs;
	    if numcols H.BMcoefs == 0 then break;
	    --print (e, numrows M, numcols M, numcols H.dBasis, dim reduceSpace polySpace H.dBasis);
  	    );
    	);
    H.deg = d;
    H
    )
    
polySpace TruncDualData := o-> H -> polySpace(H.dBasis, Reduced=>false)
dualSpace (TruncDualData,Point) := (H,p) -> dualSpace(polySpace H,p)
homogPolySpace = method()
homogPolySpace TruncDualData := H -> polySpace(H.BMintegrals*H.BMcoefs) 


-- this version gets a piece of the eliminating DS from the "usual" truncated DS 
eliminatingDual = method(Options => {Tolerance => null})
eliminatingDual (Point,Ideal,List,ZZ) := o -> (p,I,ind,d) -> eliminatingDual (p,gens I,ind,d,o)
eliminatingDual (Point,Matrix,List,ZZ) := o -> (p,Igens,ind,d) -> (
    R := ring Igens;
    if d < 0 then return dualSpace(map(R^1,R^0,0),p);
    t := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    Igens = shiftGens(p,Igens);
    n := numgens R;
    if not all(ind, i->class i === ZZ) or not all(ind, i -> i>=0 and i<n)
    then error ("expected a list of nonnegative integers in the range [0," | n | "] as 2nd parameter");
    TDD := initializeDualData(Igens,false,t);
    RdBasis := dualSpace(TDD,p);
    dBold := 0;
    while dBold != dim RdBasis do (
	dBold = dim RdBasis;
	TDD = nextTDD(TDD,t);
	RdBasis = truncate(dualSpace(TDD,p),ind,d);
	--print(dualSpace(TDD,p),RdBasis);
	);
    RdBasis
    )

truncate (PolySpace, ZZ) := (L,d) -> (
    R := ring L;
    if not L.Reduced then L = reduceSpace L;
    tGens := select(flatten entries gens L, q -> first degree q <= d);
    if #tGens == 0 then polySpace map(R^1,R^0,0) else polySpace matrix{tGens}
    )
truncate (DualSpace, ZZ) := (L,d) -> dualSpace(truncate(L.Space,d),L.BasePoint)
truncate (PolySpace, List, ZZ) := (L,ind,d) -> (
    R := ring L;
    n := numgens R;
    if not all(ind, i->class i === ZZ) or not all(ind, i -> i>=0 and i<n)
    then error ("expected a list of nonnegative integers in the range [0," | n | "] as 2nd parameter");
    indC := flatten entries submatrix'(matrix{{0..n-1}},ind);
    T := newRing(R,MonomialOrder=>{#ind,#indC},Degrees=>{(#ind:{1,0})|(#indC:{0,1})});
    TtoR := map(R,T, (vars R)_ind | (vars R)_indC);
    varList := new MutableList from n:0;
    scan(#ind, i->(varList#(ind#i) = T_i));
    scan(#indC, i->(varList#(indC#i) = T_(#ind + i)));
    RtoT := map(T,R, matrix{toList varList});
    TL := reduceSpace polySpace RtoT gens L;
    TL = truncate(TL,d);
    polySpace(TtoR gens TL)
    )
truncate (DualSpace, List, ZZ) := (L,ind,d) -> dualSpace(truncate(L.Space,ind,d),L.BasePoint)

gCorners = method(TypicalValue => Matrix, Options => {Strategy => BM, Tolerance => null, ProduceSB => false})
gCorners (Point,Ideal) := o -> (p,I) -> gCorners(p,gens I,o)
gCorners (Point,Matrix) := o -> (p,Igens) -> (
    R := ring Igens;
    t := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    Igens = sub(Igens, matrix{gens R + apply(p.Coordinates,c->sub(c,R))});

    ecart := max apply(flatten entries Igens, g->(gDegree g - lDegree g)); --max ecart of generators
    GCs := {}; -- g-corners (as pairs: monomial, degree in homogenization)
    SBs := {}; -- standard basis elements (if o.ProduceSB)
    finalDegree := max(flatten entries Igens / gDegree);
    d := 0;
    dBasis := dBasisReduced := polySpace map(R^1,R^0,0); -- Sylvester truncated dual space
    TDD := initializeDualData(Igens,true,t,Strategy=>o.Strategy); -- initial parameters for computing truncated duals
    while d <= finalDegree do (
    	TDD = nextTDD(d,TDD,t);
	dBasis = polySpace TDD;
	dBasisReduced = reduceSpace(dBasis,Tolerance=>t);
	-- Find new g-corners based on what monomials are missing from dBasis.
	newGCs := newGCorners(dBasisReduced,GCs,d,ecart);
	GCs = GCs|newGCs;
	-- If o.ProduceSB then compute a standard basis element for each new g-corner.
	if o.ProduceSB and #newGCs > 0 then SBs = SBs|newSBasis(dBasis,newGCs,d,t);
	-- Update stopping degree if there were new g-corners found.
	if #newGCs > 0 then (
	    topLCMdegree := max apply(subsets(#GCs,2),s->(
		    (a,b) := (GCs#(s#0), GCs#(s#1));
		    tdegree := max(a#1 - first degree a#0, b#1 - first degree b#0);
    		    (first degree mylcm {a#0,b#0}) + tdegree
		    ));
	    finalDegree = max(finalDegree,topLCMdegree);
	    );
	<< "-" << "- at degree " << d << ": dim " << dim dBasisReduced << ", new corners " << newGCs/first << endl;
	--print(d, finalDegree, dim dBasisReduced, newGCs/first);
	d = d+1;
	);
    GCs = if o.ProduceSB then SBs else GCs/first;
    sbReduce matrix {GCs}
    )

mylcm = L -> (
    k := #(first exponents L#0);
    exps := apply(k, i->max apply(L, l->(first exponents l)#i));
    (ring L#0)_exps
    )

-- computes s-corners from the g-corners
-- i.e. the maximal monomials not in the ideal generated by the g-corners
socles = method(TypicalValue => Matrix)
socles MonomialIdeal := I -> socles gens I
socles Matrix := gCs -> (
    R := ring gCs;
    n := numgens R;
    G := flatten entries gCs;
    candidates := subsets(G, n) / mylcm;
    S := select(candidates, c -> (
	    c != 0 and
	    all(G, g -> not isDivisible(c,g)) and
	    all(gens R, v -> any(G, g -> isDivisible(v*c,g)))
	    ));
    S = unique S;
    matrix{S}
    )

hilbertFunction DualSpace := L -> (
    if not L.Space.Reduced then L = reduceSpace L;
    tally(flatten entries gens L / first @@ degree)
    )
hilbertFunction(List,DualSpace) := (LL,L) -> (
    h := hilbertFunction L;
    apply(LL, d->(if h#?d then h#d else 0))
    )
hilbertFunction(ZZ,DualSpace) := (d,L) -> first hilbertFunction({d},L)

localHilbertRegularity = method(TypicalValue => ZZ, Options=>{Tolerance => null})
localHilbertRegularity(Point, Ideal) := o -> (p,I) -> localHilbertRegularity(p,gens I,o)
localHilbertRegularity(Point, Matrix) := o -> (p,Igens) -> (
    n := numgens ring Igens;
    gCs := gCorners(p,Igens,o);
    gCLists := apply(flatten entries gCs, l -> (listForm l)#0#0);
    LCMexp := apply(n, i -> max(apply(gCLists, l->l#i)));
    s := max{sum LCMexp - n + 1, 0}; -- an upperbound
    J := ideal gCs;
    HP := hilbertPolynomial(J, Projective=>false);
    while hilbertFunction(s-1,J) == sub(HP,(ring HP)_0 => s-1) do s = s-1;
    s
    )
     

-- Implementation of algorithm from 1996 paper of Bernard Mourrain.
-- M is the main matrix
-- E is the row matrix of all dual basis integrals (including 1)
-- B contains the most recently found generators (as coefficients in terms of E)
BMmatrix = H -> (
    Igens := if H.syl then H.hIgens else H.Igens;
    (M,E,B,Bfull,homogeneous,Seeds) := (H.BMmatrix,H.BMintegrals,H.BMbasis,H.BMcoefs,H.syl,H.Seeds);
    --print(M,E,(numcols B, numrows B),(numcols Bfull, numrows Bfull));
    R := ring Igens;
    n := numgens R;
    m := numcols Igens;
    snew := numcols B;
    offset := if homogeneous then 0 else dim Seeds;
    s := (numcols E - offset)//n; --number of dual space generators
    npairs := subsets(n,2);
    newMEs := apply(snew, i -> (
	    bcol := B_{i};
	    bpoly := (E*bcol)_(0,0);
	    E' := matrix {apply(n, k->(
			subs := matrix{apply(n, l->(if l > k then 0_R else (gens R)#l))};
			(gens R)#k * sub(bpoly,subs)))};
	    M' := innerProduct(polySpace Igens, polySpace E');
	    if not homogeneous then M' = map(R^(s+(numcols Bfull)),R^n,0) || M';
	    for j from 0 to s-1 do (
		w := apply(n,k->(bcol_(offset + j*n + k,0)));
		v := mutableMatrix(R,#npairs,n);
		for i from 0 to #npairs-1 do (
		    v_(i,npairs#i#0) =  w#(npairs#i#1);
		    v_(i,npairs#i#1) = -w#(npairs#i#0);
		    );
		M' = M' || new Matrix from v;
		);
	    (M',E')
	    ));
    
    if not homogeneous then
    M = transpose Bfull || M || map(R^(m + s*(1+#npairs) - numrows M),R^(numcols E),0)
    else M = map(R^(m + s*#npairs),R^0,0);
    (Mnew, Enew) := (matrix{newMEs/first}, matrix{newMEs/last});
    if #newMEs == 0 then (Mnew, Enew) = (map(R^(numrows M),R^0,0),map(R^(numrows E),R^0,0));
    M = M | Mnew;
    E = if homogeneous then Enew else E | Enew;
    (M,E)
    );


newGCorners = method()
newGCorners (PolySpace,List,ZZ,ZZ) := (dBasis,GCs,d,ecart) -> (
    R := ring dBasis;
    mons := flatten entries sort basis(d-ecart,d,R);
    dBasisMons := sort(flatten entries gens dBasis / gLeadMonomial);
    newGCs := {};
    i := 0;
    for m in mons do (
	while i < #dBasisMons and m > dBasisMons#i do i = i+1;
	if i < #dBasisMons and m == dBasisMons#i then (i = i+1; continue);
	if not any(GCs, g->(isDivisible(m,g#0) and gDegree m - gDegree(g#0) <= d - g#1)) then newGCs = append(newGCs, (m,d));
	);
    newGCs
    )

newSBasis = (dBasis,newGCs,d,t) -> (
    R := ring dBasis;
    Rd := sort basis(0,d,R);
    Id := orthogonalInSubspace(dBasis, polySpace Rd, t);
    iBasis := flatten entries gens reduceSpace(Id,Monomials=>Rd);
    new List from apply(newGCs, n->first select(1,iBasis, b->(lLeadMonomial b == n#0)))
    )


-- PolySpace of ideal basis through degree d.
idealBasis = method()
idealBasis(Ideal, ZZ, Boolean) := (I, d, useGDegree) -> idealBasis(gens I, d, useGDegree)
idealBasis(Matrix, ZZ, Boolean) := (igens, d, useGDegree) -> (
     R := ring igens;
     igens = first entries igens;
     genDeg := if useGDegree then gDegree else lDegree;
     p := map(R^1,R^0,0);
     for g in igens do
	 p = p|(matrix{{g}}*basis(0, d - genDeg g, R));
     polySpace p
     )

--lead monomial and lead monomial degree according to ordering associated with
--the ring (local) and reverse ordering (global)
lLeadMonomial = f -> leadMonomial last terms f
gLeadMonomial = f -> leadMonomial first terms f
lDegree = f -> first degree lLeadMonomial f
gDegree = f -> first degree gLeadMonomial f

--remove non-minimal standard basis elements
sbReduce = L -> (
    n:= numcols L;
    goodi := select(n, i->(
     	    all(0..n-1, j->(j == i or not isDivisible(lLeadMonomial L_(0,i), lLeadMonomial L_(0,j))))
	    ));
    L_goodi
    )

contract (PolySpace, PolySpace) := (S, T) -> (
    T = gens T;
    S = gens S;
    cols := for s in flatten entries S list transpose contract(s,T);
    matrix {cols}
    )

diff (PolySpace, PolySpace) := (S, T) -> (
    T = gens T;
    S = gens S;
    cols := for s in flatten entries S list transpose diff(s,T);
    matrix {cols}
    )



-- Matrix of inner products
-- PolySpace generators as rows, DualSpace generators as columns
innerProduct = method()
innerProduct (PolySpace, PolySpace) := (S, T) -> (
    M := last coefficients(gens S | gens T);
    Svec := submatrix(M,0..dim S-1);
    Tvec := submatrix'(M,0..dim S-1);
    (transpose Svec)*Tvec
    )
innerProduct (PolySpace, DualSpace) := (S, L) -> (
    Sshift := polySpace sub(gens S, matrix{(gens ring L) + coordinates L.BasePoint});
    innerProduct(Sshift, L.Space)
    )
innerProduct (RingElement, DualSpace) := (f, L) -> innerProduct(polySpace matrix{{f}}, L)
innerProduct (PolySpace, RingElement) := (S, l) -> innerProduct(S, polySpace matrix{{l}})
innerProduct (RingElement, RingElement) := (f, l) -> (
    M := last coefficients(matrix{{f,l}});
    ((transpose M_{0})*M_{1})_(0,0)
    )

reduceSpace = method(Options => {Monomials => null,Tolerance=>1e-6})
reduceSpace PolySpace := o -> S -> (
    if dim S == 0 then return polySpace(gens S,Reduced=>true);
    (mons,coefs) := coefficients(gens S, Monomials => o.Monomials);
    M := mons*(colReduce(coefs,o.Tolerance));
    polySpace(M,Reduced=>true)
    )
reduceSpace DualSpace := o -> L -> dualSpace(reduceSpace L.Space,L.BasePoint)

orthogonalInSubspace = method()
orthogonalInSubspace (DualSpace, PolySpace, Number) := (D,S,t) -> (
    M := innerProduct(S,D);
    K := numericalKernel(transpose M,t);
    polySpace((gens S)*K, Reduced=>false)
    )
orthogonalInSubspace (PolySpace, PolySpace, Number) := (T,S,t) -> (
    T' := dualSpace(T, origin(ring S));
    orthogonalInSubspace(T',S,t)
    )

isDivisible = (a, b) -> (
     dif := (listForm a)#0#0 - (listForm b)#0#0;
     all(dif, i->(i >= 0))
     )

check PolySpace := S -> (
    assert(numrows gens S == 1);
    M := last coefficients gens S;
    assert(numgens numericalImage M == numcols M);
    )
check DualSpace := D -> check polySpace gens D

---------------------------------------------
-- Numerical Linear Algebra
---------------------------------------------

numericalImage = method(Options => {Tolerance => null})
numericalImage Matrix := o -> M -> (
    R := ultimate(coefficientRing, ring M);
    tol := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    numericalImage(M,tol)
    )
numericalImage (Matrix, Number) := o -> (M, tol) -> (
    R := ultimate(coefficientRing, ring M);
    M = sub(M, R);
    if numcols M == 0 then return M;
    if numrows M == 0 then return map(R^0,R^0,0);
    if precision 1_(ring M) < infinity then (
	(svs, U, Vt) := SVD M;
	cols := positions(svs, sv->(sv > tol));
	submatrix(U,,cols)
	) else (
	gens image M
	)
    )

numericalKernel = method(Options => {Tolerance => null})
numericalKernel Matrix := o -> M -> (
    R := ultimate(coefficientRing, ring M);
    tol := if o.Tolerance === null then defaultT(R) else o.Tolerance;
    numericalKernel(M,tol)
    )
numericalKernel (Matrix, Number) := o -> (M, tol) -> (
    R := ring M;
    M = sub(M, ultimate(coefficientRing, R));
    if numrows M == 0 then return id_(source M);
    if numcols M == 0 then return map(R^0,R^0,0);
    if precision 1_R < infinity then (
	(svs, U, Vt) := SVD M;
	cols := positions(svs, sv->(sv > tol));
	submatrix'(adjointMatrix Vt,,cols)
	) else (
	gens kernel M
	)
    )

-- produces the conjugate transpose
adjointMatrix = method(TypicalValue => Matrix)
adjointMatrix Matrix := M -> (
    M' := mutableMatrix transpose M;
    for i from 0 to (numrows M')-1 do (
	for j from 0 to (numcols M')-1 do M'_(i,j) = conjugate(M'_(i,j));
	);
    matrix M'
    )

--performs Gaussian reduction on M
colReduce = method(TypicalValue => Matrix)
colReduce (Matrix, Number) := (M, tol) -> (
    M = new MutableMatrix from sub(transpose M, ultimate(coefficientRing, ring M));
    (m,n) := (numrows M, numcols M);
    i := 0; --row of pivot
    for j from 0 to n-1 do (
	if i == m then break;
	a := i + maxPosition apply(i..m-1, l->(abs M_(l,j)));
	c := M_(a,j);
	if abs c <= tol then continue;
	rowSwap(M,a,i);
	for l from 0 to n-1 do M_(i,l) = M_(i,l)/c; --rowMult(M,i,1/c); is bugged
	for k from 0 to m-1 do rowAdd(M,k,-M_(k,j),i);
	i = i+1;
	);
    M = (transpose new Matrix from M)_{0..i-1};
    if tol > 0 then clean(tol,M) else M
    )

--a list of column indices for a basis of the column space of M
basisIndices = (M, tol) -> (
    M = new MutableMatrix from sub(M, ultimate(coefficientRing, ring M));
    (m,n) := (numrows M, numcols M);
    i := 0; --row of pivot
    I := new MutableList;
    for j from 0 to n-1 do (
	if i == m then break;
	a := i + maxPosition apply(i..m-1, l->(abs M_(l,j)));
	c := M_(a,j);
	if abs c <= tol then continue;
	I#(#I) = j;
	rowSwap(M,a,i);
	for l from 0 to n-1 do M_(i,l) = M_(i,l)/c; --rowMult(M,i,1/c); is bugged
	for k from 0 to m-1 do rowAdd(M,k,-M_(k,j),i);
	i = i+1;
	);
    new List from I
    )





--
-----  Noetherian operator code
--





-- Compute Noetherian operators.
-- Caveat: variety of the extended ideal must be the origin!
--brokenNoethOps = {Normalized => false, Verbose => false} >> opts -> J -> (
--    S := ring J;--

--    -- Noether normalize, record parameters,
--    -- variables & automorphism
--    (f,I,par) := if not opts.Normalized then (
--        if opts.Verbose then print "Performing Noether normalization";
--        noetherNormalization J
--    ) else (
--        if opts.Verbose then print("Assuming in normal position wrt to " | toString((flatten entries vars S)_{-(dim J)..-1}));
--        (map(S,S), J, (flatten entries vars S)_{-(dim J)..-1})
--    );
--    var := first entries vars ring J - set par;
--    mu := degree I;--

--    -- Define new ring with correct monomial order
--    R := QQ[var, par, MonomialOrder => {#var, #par}];
--    I = sub(I,R);
--    var = var / (i -> sub(i,R));
--    par = par / (i -> sub(i,R));--
--

--    -- Saturate wrt I (?)
--    G := first entries gens gb I;
--    tg := G / (i->leadMonomial(i)) / (i -> first entries monomials(i, Variables => par));
--    tg = product flatten tg;
--    if opts.Verbose then print("t^gamma = " | toString(tg));--
--

--    -- Multiply all by t^gamma
--    b := first entries basis(0,mu-1,R, Variables => var);
--    bt := b;--

--    
--    -- Normal form of Taylor polynomial
--    bt = bt / (i->i*tg);
--    nf := bt / (i->i % I);
--    
--    -- check number of coeffs
--    beta := set flatten (nf / (i -> first entries monomials(i, Variables => var)));
--    while not mu == #beta do (
--        if opts.Verbose then print("Multiply GB by " | toString(tg));
--        nf = nf / (i->i*tg) / (i -> i % I);
--        beta = set flatten (nf / (i -> first entries monomials(i, Variables => var)));
--    );--
--

--    -- kill zeros
--    idx := positions(nf, j -> j!=0);
--    b = b_idx;
--    nf = nf_idx;--

--    -- record results in a hash table. Keys are standard monomials
--    ht := new MutableHashTable from (applyPairs(beta, (i,j) -> i=> 0));--

--    -- populate the hash table
--    for i from 0 to #idx - 1 do (
--        if opts.Verbose then print("Processing term " | toString(i+1));
--        (mon,coe) = coefficients (nf#i, Variables => var);
--        mon = first entries mon;
--        coe = flatten entries coe;
--        for j from 0 to #mon - 1 do(
--            ht#(mon#j) = ht#(mon#j) + (coe#j) * (b#i);
--        );
--    );
--    -- check number of monomials
--    if opts.Verbose then print"check number of monomials";
--    res := values ht;--

--    -- divide by factorials
--    if opts.Verbose then print"divide by factorials";
--    res = res / (i -> addFactorials(i, var));--

--    -- replace x by dx
--    if opts.Verbose then print"replace x by dx";
--    W := makeWA (QQ[var,par]);
--    WARules := (options(W))#WeylAlgebra / (i -> (i#0)_W => (i#1)_W);
--    diffRules := select (WARules, i -> member(sub((i#0),R), var));
--    res = res / (i -> sub(i, W)) / (i -> sub(i, diffRules));--

--    -- Convert to NOps of previous ideal
--    if opts.Verbose then print"Convert to NOps of previous ideal";
--    fMatrix := mapToMatrix(f);
--    nVars := #par + #var;
--    zeroMatrix := map(QQ^nVars, nVars, 0);
--    f'Matrix := mapToMatrix(inverse f);
--    rules := matrix({{f'Matrix, zeroMatrix},{zeroMatrix,transpose fMatrix}}) * (transpose vars W);
--    rules = flatten entries rules;
--    phi := map(W,W,rules);--

--    -- return the results
--    if opts.Verbose then print"return the results";
--    use S;
--    res / phi
--)

diffAlg = method()
diffAlg(Ring) := R -> (
    diffVars := apply(gens R, i -> value("symbol d" | toString(i)) );
    R[diffVars]
    )

-- addFactorials = (f, var) -> (
--     (coe, mon) := coefficients(f, Variables => var);
--     factorials := first entries coe / (i -> first exponents(i)) / (i -> product( i / (j -> j!)));
--     res := apply(flatten entries mon, factorials, (i,j) -> i/j);
--     sum apply(first entries coe, res, (i,j) -> i*j)
-- )


-- Given an element N in Weyl algebra and a polynomial
-- f, compute the result of applying N to f.
applyNOp = (N, f) -> (
    m := map(ring f, ring N, vars ring f);
    (a,b) := coefficients matrix{{N}};
    diff(m a,f)*sub(b,ring f)
    )

-- Try to see if gens of I applied with all Noeth Ops
-- vanish on rad(I)
visualCheck = (nops, I) -> (
    foo := table(nops, I_*, (n,i) -> applyNOp(n,i));
    netList applyTable(foo, (i -> i%(radical I)))
)

sanityCheck = (nops, I) -> (
    foo := flatten table(nops, I_*, (n,i) -> applyNOp(n,i));
    all(foo, i -> sub(i,ring I)%(radical I) == 0)
)


noethOps = method(Options => {DegreeLimit => 5, DependentSet => null}) 
noethOps (Ideal, Ideal) := List => opts -> (I, P) -> (
    R := ring I;
    depVars := if opts.DependentSet === null then gens R - set support first independentSets P
            else opts.DependentSet;
    indVars := gens R - set depVars;
    S := ((coefficientRing R)[indVars])[depVars];
    dS := zeroDimensionalDual(sub(P,S),sub(I,S),Normalize=>false);
    R' := diffAlg R;
    return flatten entries sub(gens dS,R');
--    bx := flatten entries basis(0,opts.DegreeLimit,R, Variables => gens R);
--    bd := basis(0,opts.DegreeLimit,R, Variables => depVars);

--    elapsedTime M := diff(bd, transpose gens idealBasis(I,opts.DegreeLimit,false));
--    elapsedTime M' := sub(M,R/P);
--    elapsedTime K := gens trim kernel M';

    -- Return elements in WeylAlgebra for nice formatting
--    bdd := sub(bd, vars R');
--    flatten entries (bdd * sub(K, R'))
)
noethOps (Ideal) := List => opts -> (I) -> noethOps(I, ideal gens radical I, opts)
noethOps (Ideal, Point) := List => opts -> (I, p) -> (
    P := ideal ((gens ring I) - p.Coordinates);
    noethOps(I,P,opts)
)

approxKer = method(Options => {Tolerance => null})
approxKer(Matrix) := Matrix => o -> A -> (
    tol := if o.Tolerance === null then defaultT(ring A) else o.Tolerance;
    d := numcols A;
    (S,U,Vh) := SVD A;
    n := #select(S, s -> clean(tol, s) == 0);
    K := transpose Vh^{d-n..d-1};
    if K == 0 then K else conjugate K
)


numNoethOps = method(Options => options noethOps ++ options approxKer)
numNoethOps (Ideal, Matrix) := List => opts -> (I, p) -> (
    R := ring I;
    var := gens R - set support first independentSets I;
    bx := flatten entries basis(0,opts.DegreeLimit,R, Variables => gens R);
    bd := basis(0,opts.DegreeLimit,R, Variables => var);

    elapsedTime M := diff(bd, transpose matrix {flatten (table(bx,I_*,(i,j) -> i*j))});
    elapsedTime M' := sub(M,p);
    elapsedTime K := colReduce(numericalKernel M', opts.Tolerance);

    -- Return elements in WeylAlgebra for nice formatting
    R' := diffAlg R;
    bdd := sub(bd, vars R');
    flatten entries (bdd * sub(K, R'))
)

conjugate(Matrix) := Matrix => M -> (
    matrix table(numrows M, numcols M, (i,j) -> conjugate(M_(i,j)))
)

--socleMonomials = method()
--socleMonomials(Ideal) := List => (I) -> (
--    R := ring I;
--    if codim I < dim R then error("Expected Artinian ideal");
--    inI := monomialIdeal leadTerm gens gb I;
--    irreducibleDecomposition inI / gens / 
--                        entries / first / 
--                        product / exponents / first / 
--                        (i -> (i / (j -> j-1))) / (e -> R_e)
--)

coordinateChangeOps = method()
coordinateChangeOps(RingElement, RingMap) := RingElement => (D, f) -> (
    R := f.target;
    WA := ring D;
    A := f.matrix // vars R;
    A' := inverse A;
    (a,b) := coefficients D;
    b = sub(f sub(b,R), WA);
    

    psi := transpose (sub(transpose A',WA) * (transpose vars WA));
    a = (map(WA,WA,psi)) a;
    (a*b)_(0,0)
)


noethOpsFromComponents = method()
noethOpsFromComponents(HashTable) := List => H -> (
    nops := flatten values H;
    R := ring first nops;
    nops = unique (nops / (f -> sub(f, R)));
    Ps := apply(nops, D -> select(keys H, P -> any(H#P, D' -> D == sub(D',ring D))));
    
    mults := Ps / (Lp -> 
        if set Lp === set keys H then 1_R else (
            J := intersect(keys H - set Lp);
            (sub(gens J,R) * random(R^(#J_*), R^1))_(0,0)
        )
    );

    apply(mults, nops, (i,j) -> i*j)
)

testMM = (nop, mm) -> (
    (sub(nop, ring mm) % mm) == 0
)

-- Divide by gcd of coefficients
modConstant = f -> (
    (mon, coe) := coefficients f;
    g := (flatten entries coe) / (i -> sub(i,QQ));
    f // gcd(g)
)

beginDocumentation()

doc ///
     Key
     	  DualSpaces
     Headline
     	  numerically compute local dual space and Hilbert functions
     Description
     	  Text
	       The @EM "DualSpaces"@ package includes algorithms for computing local dual 
	       spaces of polynomial ideals, and related local combinatorial data about its scheme structure.  These 
	       techniques are numerically stable, and can be used with floating point arithmetic over the complex numbers.  
	       They provide a viable alternative in this setting to purely symbolic methods such as standard bases.  
	       In particular, these methods can be used to compute initial ideals, local Hilbert functions and Hilbert regularity.

	       Methods for computing and manipulating local dual spaces:

	       @UL {
		   {TO "truncatedDual"},
		   {TO "zeroDimensionalDual"},
		   {TO "eliminatingDual"},
		   {TO "localHilbertRegularity"},
		   {TO "gCorners"},
		   {TO "socles"},
		   {TO "innerProduct"},
                   {TO "reduceSpace"}
		   }@

	       Auxiliary numerical linear algebra methods:

	       @UL {
		   {TO "numericalKernel"},
		   {TO "numericalImage"},
		   {TO "colReduce"},
		   {TO "adjointMatrix"}
		   }@

	       The algorithm used for computing truncated dual spaces is that of B. Mourrain ("Isolated points, duality and residues." 
	       J. Pure Appl. Algebra, 117/118:469–493, 1997).  To compute the initial ideal and Hilbert regularity of positive dimensional
	       ideals we use the algorithm of R. Krone ("Numerical algorithms for dual bases of positive-dimensional ideals." Journal of
               Algebra and Its Applications, 12(06):1350018, 2013.).  This package depends on the package @TO NAGtypes@.
///


doc ///
     Key
          truncatedDual
	  (truncatedDual,Point,Ideal,ZZ)
	  (truncatedDual,Point,Matrix,ZZ)
     Headline
          truncated dual space of a polynomial ideal
     Usage
          S = truncatedDual(p, I, d)
     Inputs
     	  p:Point
	  I:Ideal
               or a one-row @TO Matrix@ of generators
	  d:ZZ
     Outputs
          S:DualSpace
     Description
          Text
	       Computes a basis for the local dual space of a polynomial ideal localized at point p, truncated at degree d.
	       Elements are expressed as elements of the polynomial ring of the ideal although this is an abuse of notation.
	       They are really elements of the dual ring.
	  Example
	       R = CC[x,y];
	       I = ideal{x^2, y*x}
	       truncatedDual(origin(R),I,3)
	  Text
	       The functionals in the dual at a point p are expressed in coordinates centered at p.
	  Example
	       p = point matrix{{0., 1}}
	       truncatedDual(p,I,3)
	  Text
	       Over inexact fields, the computation accounts for the possibility of small numerical error in the point p.
	       The optional argument @TO "Tolerance (DualSpaces)"@ can be set to adjust the tolerance of the numerical computations.
	       Higher degree dual computations generally require higher accuracy in the input and larger tolerance value to complete correctly.
	       
	       In this example, the point q is slightly away from the variety of I, but an appropriate @TT "Tolerance"@ value can overcome the error. 
	  Example
	       q = point matrix{{0. + 1e-10, 1}}
	       tol = 1e-6;
	       S = truncatedDual(q,I,3, Tolerance => tol)
	       (m,c) = coefficients gens S;
	       m*clean(tol, c)
	  Text
	       See also @TO zeroDimensionalDual@.
///

TEST ///
R = CC[x,y]
I1 = ideal{x^2,x*y}
D1 = truncatedDual(origin R, I1, 4)
assert(hilbertFunction({0,1,2,3,4}, D1) == {1,2,1,1,1})
///

doc ///
     Key
          zeroDimensionalDual
	  (zeroDimensionalDual,Point,Ideal)
	  (zeroDimensionalDual,Point,Matrix)
     Headline
          dual space of a zero-dimensional polynomial ideal
     Usage
          S = zeroDimensionalDual(p, I)
     Inputs
     	  p:Point
	  I:Ideal
              or a one-row @TO Matrix@ of generators
     Outputs
          S:DualSpace
     Description
          Text
	       Computes a reduced basis of the dual space of a zero-dimensional ideal.  It does not check if the ideal is
	       zero-dimensional and if not then termination will fail.
	       Elements are expressed as elements of the polynomial ring of the ideal although this is an abuse of notation.
	       They are really elements of the dual ring.
	  Example
	       R = QQ[a,b];
	       I = ideal{a^3,b^3}
	       D = zeroDimensionalDual(origin(R), I)
	       dim D
	  Text
	       The dimension of the dual space at p is the multiplicity of the solution at p.
	  Example
	       S = CC[x,y];
	       J = ideal{(y-2)^2,y-x^2}
	       p = point matrix{{1.4142136,2}};
	       D = zeroDimensionalDual(p, J)
	       dim D
	  Text
	       See also @TO truncatedDual@.
     Caveat
	  The computation will not terminate if I is not locally zero-dimensional at the chosen point.  This is not checked.
///

TEST ///
R = CC[x,y]
I2 = ideal{x^2,y^2}
D2 = zeroDimensionalDual(origin R, I2)
assert(hilbertFunction({0,1,2,3,4}, D2) == {1,2,1,0,0})
D2' = zeroDimensionalDual(point matrix{{1,1}}, I2)
assert(dim D2' == 0)
///

doc ///
     Key 
          gCorners
	  (gCorners,Point,Ideal)
	  (gCorners,Point,Matrix)
	  [gCorners,ProduceSB]
	  ProduceSB
     Headline
          generators of the initial ideal of a polynomial ideal
     Usage
          G = gCorners(p, I)
     Inputs
     	  p:Point
	  I:Ideal
              or a one-row @TO Matrix@ of generators
     Outputs
          G:Matrix
	       generators of the initial ideal in a one-row matrix
     Description
          Text
	       Computes the generators of the initial ideal of an ideal, with respect to a local order.  The ring of the
	       ideal should be given a (global) monomial order and the local order will be taken to be the reverse order.
	       The point p is moved to the origin, so the monomial generators represent terms of the Taylor expansion at p.
	  Example
	       R = CC[x,y];
	       I = ideal{x^2-y^2}
	       p = point matrix{{1,1}};
	       gCorners(p, I)
	  Text
	       If the optional argument @TT "ProduceSB"@ is set to true, the output is instead a matrix of elements of the ideal
	       with the p translated to the origin such that the lead terms generate the inital ideal, i.e. a standard basis.
	       Note that the coordinates of the standard basis elements are translated to be centered at the point p.
	  Example
	       S = gCorners(p, I, ProduceSB=>true)
	  Example
	       R = CC[x,y,z];
	       J = ideal{z*(x*y-4), x-y}
	       q = point matrix{{1.4142136, 1.4142136, 0}};
	       gCorners(q, J, Tolerance=>1e-5)
	       gCorners(q, J, ProduceSB=>true)
///

TEST ///
R = CC[x,y]
M = matrix {{x^2-x*y^2,x^3}}
--M = matrix {{x*y, y^2}}
p = point matrix{{0_CC,0_CC}}
q = point matrix{{0_CC,1_CC}}
assert(numcols gCorners(p,M) == 2)
assert(numcols gCorners(q,M) == 1)
LDZ = reduceSpace truncatedDual(p,M,5,Strategy=>DZ)
LBM = reduceSpace truncatedDual(p,M,5,Strategy=>BM)
assert(dim LDZ == dim LBM)
///

doc ///
     Key
          socles
          (socles,MonomialIdeal)
	  (socles,Matrix)
     Headline
          socle corners of a monomial ideal
     Usage
          S = socles I
	  S = socles G
     Inputs
          I:MonomialIdeal
	  G:Matrix
	       one row of g-corners
     Outputs
          S:Matrix
	       socle corners in a one-row matrix
     Description
          Text
	       Computes the maximal monomials which are not in the monomial ideal, i.e. the "socle-corners".
	  Example
	       R = CC[x,y,z];
	       I = monomialIdeal{x^2,y^2,z^2}
	       socles I
	       socles I^2
	       G = vars R
	       socles G
///

TEST ///
R = CC[x,y]
G1 = matrix{{x^2,x*y^2,y^4}}
assert(socles G1 == matrix {{x*y, y^3}})
G2 = matrix{{x*y^2,x^2*y^2,y^4}}
assert(socles G2 == matrix {{y^3}})
G3 = matrix{{x*y^2,x^2*y^2}}
assert(socles G3 == matrix {{}})
///

doc ///
     Key
          localHilbertRegularity
	  (localHilbertRegularity,Point,Matrix)
	  (localHilbertRegularity,Point,Ideal)
	  
     Headline
          regularity of the local Hilbert function of a polynomial ideal
     Usage
          d = localHilbertRegularity(p,I)
     Inputs
          p:Point
	  I:Ideal
              or a one-row @TO Matrix@ of generators
     Outputs
          d:ZZ
     Description
          Text
	       The gCorners of the ideal are computed in order to find the Hilbert polynomial, which is
	       compared to the Hilbert function to find the degree of regularity, which is the degree at
	       which the two become equal.
	  Example
	       R = CC[x,y];
	       I = ideal{x^2,x*y}
	       d = localHilbertRegularity(origin R, I)
	       D = truncatedDual(origin R, I, 3)
	       L = hilbertFunction({0,1,2,3}, D)
	  Text
	       See also @TO gCorners@.
///


doc ///
     Key
          eliminatingDual
	  (eliminatingDual,Point,Matrix,List,ZZ)
	  (eliminatingDual,Point,Ideal,List,ZZ)
     Headline
          eliminating dual space of a polynomial ideal
     Usage
          S = eliminatingDual(p, I, v, d)
     Inputs
     	  p:Point
	  I:Ideal
              or a one-row @TO Matrix@ of generators
	  v:List
	       a list of the integers designating which variables to bound
	  d:ZZ
	       the degree bound for the designated variables
     Outputs
          S:DualSpace
     Description
          Text
	       Given a list of variable indices, compute the a basis for all dual elements
	       orthogonal to I which have total degree in the variables on the list bounded by d.
	  Example
	       R = CC[x,y];
	       I = ideal{x^2-y^3}
	       --bound the x degree to 2
	       eliminatingDual(origin R, I, {0}, 2)
	  Text
	       This function generalizes @TO truncatedDual@ in that if v includes all the variables
	       in the ring, then its behavior is the same.
	  Example
	       eliminatingDual(origin R, I, {0,1}, 2)
	  Text
	       See also @TO truncatedDual@.
     Caveat
	  The space of dual elements satisying the conditions is not in general of finite dimension.
	  If the dimension is infinite, this function will not terminate.  This is not checked.  To ensure
	  termination, the local dimension of I at p should not exceed the length of v, and certain genericity
	  constraints on the coordinates must be met.
///

doc ///
     Key
          innerProduct
	  (innerProduct,PolySpace,DualSpace)
	  (innerProduct,PolySpace,PolySpace)
	  (innerProduct,RingElement,DualSpace)
	  (innerProduct,PolySpace,RingElement)
	  (innerProduct,RingElement,RingElement)
     Headline
          Applies dual space functionals to polynomials
     Usage
          M = innerProduct(S, D)
     Inputs
	  S:PolySpace
	  D:DualSpace
     Outputs
          M:Matrix
	       containing the values of the generators of D applied to the generators of S
     Description
          Text
	       The dual space represents functionals from the polynomial ring to the base field.
	       Given a polySpace S with n generators f_1,...,f_n and a dualSpace D with m generators
	       p_1,...,p_m, innerProduct returns a nxm matrix M over the base field whose entries are p_j(f_i).
	  Text
	       A dual functional is applied to a polynomial by taking the standard inner product of their coefficient
	       vectors.  In other words, the functional represented by the monomial a acts on monomials in the
	       polynomial ring as a(a) = 1 and a(b) = 0 for all other monomials b.
	  Example
	       R = CC[x,y];
	       S = polySpace matrix{{x+y,2*x+y^2}};
	       D = dualSpace(matrix{{1,x,y}}, origin R);
	       M = innerProduct(S, D)
	  Text
	       @TT "innerProduct"@ can also be called with one or both inputs @ofClass RingElement@.  If both arguments
	       are single elements, the output is also a ring element rather than a matrix.
	  Example
	       innerProduct(S, 1+x)
	       innerProduct(x, D)
	       innerProduct(x, 1+x)
///

TEST ///
R = CC[x,y]
S = polySpace(basis(3,R))
P = innerProduct(S,S)
assert(all((0,0)..(3,3), i->(P_i == if i#0 == i#1 then 1 else 0)))
///

doc ///
     Key
          orthogonalInSubspace
	  (orthogonalInSubspace,DualSpace,PolySpace,Number)
	  (orthogonalInSubspace,PolySpace,PolySpace,Number)
     Headline
          Orthogonal of a space
     Usage
          S = orthogonalInSubspace(D, T, tol)
     Inputs
	  D:DualSpace
	       or @ofClass PolySpace@ a space of which to find the orthogonal
	  T:PolySpace
	       ambient space
	  tol:Number
	       a positive number, the numerical tolerance
     Outputs
          S:PolySpace
     Description
          Text
	       Computes the subspace of polynomial space T which is orthogonal to the dual space (or polynomial space) D.
	  Example
	       R = CC[x,y];
	       T = polySpace matrix{{1,x,y}};
	       D = dualSpace(matrix{{x-y}}, origin R);
	       S = orthogonalInSubspace(D, T, 1e-6)
///

doc ///
     Key
          reduceSpace
	  (reduceSpace,DualSpace)
	  (reduceSpace,PolySpace)
	  [reduceSpace,Monomials]
     Headline
          reduce the generators of a space
     Usage
          S = reduceSpace T
     Inputs
     	  T:DualSpace
	       or @ofClass PolySpace@
     Outputs
          S:DualSpace
	       or @ofClass PolySpace@
     Description
          Text
	       Reduces the generators of a DualSpace or PolySpace so that the new generators are linearly independent, and each has
	       a distinct lead monomial.  This is achieved by Gaussian reduction.
	  Example
	       R = CC[x,y];
	       T = polySpace matrix{{x,y,x-y+1e-10}}
	       S = reduceSpace T
	       S = reduceSpace(T, Tolerance=>1e-12)
///

doc ///
     Key
          (truncate,DualSpace,ZZ)
	  (truncate,PolySpace,ZZ)
     Headline
          truncate a polynomial space or dual space
     Usage
          S = truncate(T, d)
     Inputs
     	  T:DualSpace
	       or @ofClass PolySpace@
	  d:ZZ
	       the degree bound for the designated variables
     Outputs
          S:DualSpace
     Description
          Text
	       Truncates a dual space or polynomial space T so that the total degree of all terms does not exceed d.
	  Example
	       R = CC[x,y];
	       I = ideal {x-y};
	       D = truncatedDual(origin R, I, 5)
	       truncate(D, 3)
///

doc ///
     Key
	  (truncate,DualSpace,List,ZZ)
	  (truncate,PolySpace,List,ZZ)
     Headline
          truncate a polynomial space or dual space
     Usage
          S = truncate(T, v, d)
     Inputs
     	  T:DualSpace
	       or @ofClass PolySpace@
	  v:List
	       a list of the integers designating which variables to bound
	  d:ZZ
	       the degree bound for the designated variables
     Outputs
          S:DualSpace
     Description
          Text
	       Truncates a dual space or polynomial space T, so that the total degree of the specified variables is bounded by d.
	  Example
	       R = CC[x,y];
	       I = ideal {x,y};
	       D = zeroDimensionalDual(origin R, I^3)
	       --truncate the x degree to 1
	       truncate(D, {0}, 1)
///

doc ///
     Key
	  (hilbertFunction,ZZ,DualSpace)
	  (hilbertFunction,DualSpace)
	  (hilbertFunction,List,DualSpace)
     Usage
          k = hilbertFunction(d,S)
	  T = hilbertFunction(S)
	  K = hilbertFunction(D,S)
     Inputs
     	  d:
		an integer, or a @TO List@ of arguments
	  S:DualSpace
		a dual space
     Outputs
          T:Tally
	      of dimensions for each degree
	  k:ZZ
	      dimensions
	  K:List
	      of dimensions
     Description
          Text
	      counts the dimension of the dual space in each degree.  A single degree can be
	      specified, or a list of degrees.  If no degree is specified, then a Tally
	      is returned pairing each degree with its dimension.
	  Example
	       R = CC[x,y];
	       I = (ideal {x,y})^4;
	       D = zeroDimensionalDual(origin R, I)
	       hilbertFunction D
///

doc ///
     Key
          "Tolerance (DualSpaces)"
	  [localHilbertRegularity,Tolerance]
	  [eliminatingDual,Tolerance]
	  [gCorners,Tolerance]
	  [zeroDimensionalDual,Tolerance]
	  [truncatedDual,Tolerance]
	  [reduceSpace,Tolerance]
     Headline
          optional argument for numerical tolernace
     Description
          Text
	       Many of the numerical operations require a tolerance value, below which numbers are considered to be numerically zero.
	       {\tt Tolerance} should be given a non-negative real number.
	       Over exact fields the default value is zero, while for inexact fields the default is 1e-6.
	       
	       
	       
	       See also @TO Tolerance@.
///

doc ///
     Key
          "Strategy (DualSpaces)"
	  [gCorners,Strategy]
	  [truncatedDual,Strategy]
	  [zeroDimensionalDual,Strategy]
     Headline
          optional argument for dual space algorithm
     Description
          Text
	       The two available algorithm choices are

	       @UL {
		   {TOH "BM"},
		   {TOH "DZ"}
		   }@

	       Both should produce roughly 
	       the same output, but there may be differences in performance.  Generally @TT "BM"@ should be more efficient and is the default, 
	       but @TT "DZ"@ is left as an option to the user.

///

doc ///
     Key
	  DZ
     Headline
          Macaulay matrix algorithm for dual space computation
     Description
          Text
	       This algorithm is a strategy for computing the truncated dual space of an ideal I at degree d.
	       A matrix is formed with a column for each monomial in the ring of I of degree at most d, and a row for each
	       monomial multiple of a generator of I which has any terms of degree d or less, storing its coefficients.
	       Any vector in the kernel of this matrix is the coeffeicient vector of an element of the dual space.
	  
	       See: B.H. Dayton and Z. Zeng. Computing the multiplicity structure in solving polynomial systems. In M. Kauers,
	       editor, @EM "Proceedings of the 2005 International Symposium on Symbolic and Algebraic Computation"@, pages 116-123. ACM, 2005.
///

doc ///
     Key
	  BM
     Headline
          Mourrain algorithm for dual space computation
     Description
          Text
	       The Mourrain algorithm is a strategy for computing the truncated dual space of an ideal I at degree d.
	  
	       See: B. Mourrain. Isolated points, duality and residues. @EM "J. Pure Appl. Algebra"@, 117/118:469-493, 1997. 
	       Algorithms for algebra (Eindhoven, 1996).
///

doc ///
     Key
          numericalImage
	  (numericalImage,Matrix,Number)
	  (numericalImage,Matrix)
	  [numericalImage,Tolerance]
     Headline
          Image of a matrix
     Usage
          V = numericalImage(M, tol)
     Inputs
	  M:Matrix
	  tol:Number
	       a positive number, the numerical tolerance
     Outputs
          V:Matrix
     Description
          Text
	       Computes the image of a matrix M numerically using singular value decomposition.
	  Example
	       M = matrix {{1., 0, 1}, {0, 1, 1}, {1, 0, 1}}
	       numericalImage(M, 0.01)
	  Text
	       Singular values less than the tolerance are treated as zero.
	  Example
	       M = matrix {{0.999, 2}, {1, 2}}
	       numericalImage(M, 0.01)
///

TEST ///
M = matrix {{0.999, 2}, {1, 2}}
Mimage = numericalImage(M, 0.01)
assert(numcols Mimage == 1)
///

doc ///
     Key
          numericalKernel
	  (numericalKernel,Matrix,Number)
	  (numericalKernel,Matrix)
	  [numericalKernel,Tolerance]
     Headline
          Kernel of a matrix
     Usage
          V = numericalKernel(M, tol)
     Inputs
	  M:Matrix
	  tol:Number
	       a positive number, the numerical tolerance
     Outputs
          V:Matrix
     Description
          Text
	       Computes the kernel of a matrix M numerically using singular value decomposition.
	  Example
	       M = matrix {{1., 1, 1}}
	       numericalKernel(M, 0.01)
	  Text
	       Singular values less than the tolerance are treated as zero.
	  Example
	       M = matrix {{1., 1}, {1.001, 1}}
	       numericalKernel(M, 0.01)
///

doc ///
     Key
          colReduce
	  (colReduce,Matrix,Number)
     Headline
          Column reduces a matrix
     Usage
          N = colReduce(M, tol)
     Inputs
	  M:Matrix
	  tol:Number
	       a positive value, the numerical tolerance
     Outputs
          N:Matrix
	       in reduced column echelon form
     Description
          Text
	       Performs Gaussian column reduction on a matrix M, retaining only the linearly independent columns.
	  Example
	       M = matrix {{1., 2, 3}, {2, 4, 0}}
	       colReduce(M, 0.01) 
	  Text
	       Entries with absolute value below the tolerance are treated as zero and not used as pivots.
	  Example
	       N = matrix {{0.001, 0, 0}, {1, 1, 3}, {2, 2, 5.999}}
	       colReduce(N, 0.01)
///

TEST ///
N = matrix {{0.001, 0, 0}, {1, 1, 3}, {2, 2, 5.999}}
N = colReduce(N, 0.01)
assert(numcols N == 1)
///

doc ///
     Key
          adjointMatrix
	  (adjointMatrix,Matrix)
     Headline
          Conjugate transpose of a complex matrix
     Usage
          N = adjointMatrix M
     Inputs
	  M:Matrix
     Outputs
          N:Matrix
     Description
          Text
	       Returns the conjugate transpose of a matrix with complex entries.
	  Example
	       M = matrix {{1+ii,2*ii},{0,1}}
	       adjointMatrix M
///




-------------- Noetherian operators tests



TEST ///
R = QQ[x,y,z]
I = ideal(x^2 - y, y^2)
nops = noethOps(I, DegreeLimit => 10)
assert(sanityCheck(nops, I))
///

TEST ///
R = QQ[x_0..x_3]
S = QQ[s,t]
I0 = ker map(S,R,{s^5,s^3*t^2, s^2*t^3, t^5})
nops = noethOps(I0, DegreeLimit => 10)
assert(sanityCheck(nops, I0))
I1 = ideal(x_0^2, x_1^2, x_2^2)
nops = noethOps(I1, DegreeLimit => 10)
assert(sanityCheck(nops,I1))
///


TEST ///
R = QQ[x,y]
I = ideal((x-1)^2,(x-1)*(y+1),(y+1)^3)
J = ideal((x)^2,(x)*(y),(y)^3)
Ps = associatedPrimes I
nopsI = noethOps(I, first Ps)
W = ring nopsI
assert(set flatten entries gens nopsI === set{W_1^2, W_0, W_1, 1_W})
nopsJ = noethOps(J, ideal gens R)
W = ring nopsJ
assert((sort flatten entries gens nopsJ) == sort{W_1^2, W_0, W_1, 1_W})
///


TEST /// -- Linear coordinate change test
R = QQ[x,y]
I = ideal(x^2*(y-x))
f = map(R,R,{2*x+y,x+y})
J = f I
NI = noethOps I
NJ = noethOps J
WI = ring first NI
WJ = ring first NJ
convertedNI = NI / (i-> sub(coordinateChangeOps(i,f), WJ))
assert sanityCheck(convertedNI,J)
assert sanityCheck(NJ,J)
assert sanityCheck(NI,I)
assert(set NI === set{1_WI, WI_0*WI_2 - WI_1*WI_2})
assert(set NJ === set{1_WJ, WJ_0*WJ_2})
assert(set convertedNI === set{1_WJ, WJ_0*WJ_2 - WJ_0*WJ_3})
///

TEST /// -- numNoethOps test
R = CC[x,y,t]
I = ideal(x^2-t*y, y^2)
nv = numericalIrreducibleDecomposition(I, Software => BERTINI)
Wsets = select(flatten values nv, x -> class x === WitnessSet )
pt = first (Wsets / sample / matrix)
l = numNoethOps(I, pt, Tolerance => 1e-2)
assert(toString(l/terms/(m -> m/exponents)) === "{{{{0, 0, 0, 0, 0, 0}}}, {{{0, 0, 0, 1, 0, 0}}}, {{{0, 0, 0, 2, 0, 0}}, {{0, 0, 0, 0, 1, 0}}}, {{{0, 0, 0, 3, 0, 0}}, {{0, 0, 0, 1, 1, 0}}}}")
///

end


restart
debug loadPackage("DualSpaces", Reload => true)
installPackage("DualSpaces", RemakeAllDocumentation => true)
--loadPackage "DualSpaces"


needsPackage "NumericalAlgebraicGeometry"


R = CC[x,y]
M = matrix {{x^2-x*y^2,x^3}}
--M = matrix {{x*y}}
--M = matrix {{x^9 - y}}
p = point matrix{{0_CC,0_CC}}
q = point matrix{{1_CC,0_CC}}
L = reduceSpace truncatedDual(M,p,6,Strategy=>DZ)
L = reduceSpace truncatedDual(M,p,6,Strategy=>BM)
--shiftDual(L,q,8)
G = matrix{{x^2,x*y^2,y^4}}
socles G

hilbertFunction(toList(0..8),L)
dualInfo(M,p,Truncate=>8)
standardBasis(M)
dualHilbert(M,Truncate=>25)
dualBasis(M)
dualInfo(M)
dualInfo(M, Strategy=>DZ)
dualInfo(M,Point=>{0.01,0.01})

-- small example
restart
loadPackage ("DualSpaces", Reload => true)
R = QQ[x,y,z, MonomialOrder => {Weights=>{-1,-1,-1}}, Global => false]
M = matrix {{x-y^2, y-z^2}}
RH = QQ[t,x,y,z,MonomialOrder => {Weights=>{-1,-1,-1,-1}}, Global => false]
MH = homogenize(sub(M,RH),t)	  
DH = flatten entries gens dualBasis(MH,4)
DHhat = sub(matrix{select(DH, Q->first first last listForm Q >= 1)}, {t=>1})
(
     sort flatten entries gens dualBasis(M,3)
     ) == (
     sort unique flatten entries sub(DHhat,R)
     )
