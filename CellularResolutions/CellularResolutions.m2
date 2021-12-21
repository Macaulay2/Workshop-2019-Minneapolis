-- -*- coding: utf-8 -*-
newPackage(
    "CellularResolutions",
    Version => "0.9",
    Date => "July 22, 2019",
    Authors => {
        {Name => "Jay Yang", Email => "yangj306@mcmaster.ca"},
        {Name => "Aleksandra Sobieska", Email => "asobieska@math.wisc.edu"}
        },
    Headline => "A package for cellular resolutions of monomial ideals",
    AuxiliaryFiles => true, -- set to true if package comes with auxiliary files
    PackageExports => {"Polyhedra", "SimplicialComplexes", "Posets"}
    )

export {"CellComplex",
        "Cell",
        "cellComplex",
        "isCycle",
        "isSimplex",
	"boundaryCells",
	"relabelCellComplex",
	"InferLabels",
	"cells",
        "cellLabel",
        "newCell",
        "newSimplexCell",
        "isFree",
        "isMinimal",
	"Reduced",
        "CellDimension",
	"maxCells",
        "cellComplexSphere",
        "cellComplexRPn",
        "cellComplexTorus",
	"taylorComplex",
	"hullComplex"
        }
protect labelRing
protect cellDimension
protect CellDimension

CellComplex = new Type of HashTable
CellComplex.synonym = "cell complex"
--Note, the mutable hash table means that equality works "Correctly"
Cell = new Type of MutableHashTable
Cell.synonym = "cell"

--returns a hashtable of lists of cells indexed by dimension
cellsFromMaxCells := lst -> (
    pendingCells := set lst;
    finishedCells := {};
    while #pendingCells !=0 do (
        c := (elements pendingCells)#0;
        pendingCells = pendingCells + ((set ((boundary c)/first)) - finishedCells) - {c};
        finishedCells = append(finishedCells,c);
        );
    partition(dim,finishedCells)
    )

--Private constructor, creates the cache
mkCellComplex := (labelRingVal, cellsVal, maxCellsVal) -> (
    new CellComplex from {
        symbol labelRing => labelRingVal,
        symbol cells => cellsVal,
        cache => new CacheTable from (
            if maxCellsVal === null
            then {}
            else {symbol maxCells => maxCellsVal})
	}
    )

--Adds a single -1 cell by default
--TODO: create an option to make a void complex
cellComplex = method()
cellComplex(Ring,List) := (R,maxCells) -> (
    mkCellComplex(R, (maxAndAllCells maxCells)#1, (maxAndAllCells maxCells)#0)
    )

cellComplex(SimplicialComplex) := (C) -> (
    S := ring C;
    Cfaces := applyValues(new HashTable from faces C,flatten @@ entries);
    --cells indexes Cells by monomials corresponding to faces of the simplicial complex
    cells := new MutableHashTable from {};
    for i from 0 to dim C do (
        for simplex in Cfaces#i do (
            bd := if i==0
                  then {}
                  else for x in gens S list (if simplex%x==0 then cells#(simplex//x) else continue);
            cells#simplex = newCell bd
            );
        );
    cellComplex(S,values cells)
    )

maxAndAllCells = method()
maxAndAllCells(List) := (lst) -> (
    if #lst == 0 then return (new HashTable,new HashTable);
    bdfn := c -> set boundaryCells c;
    maxcells := set lst;
    bdcells := sum (maxcells/bdfn);
    allcells := maxcells + bdcells;
    maxcells = maxcells - bdcells;
    while #bdcells != 0 do (
	bdcells = sum (bdcells/bdfn);
	allcells = allcells + bdcells;
	maxcells = maxcells - bdcells;
	);
    (partition(dim,toList maxcells), partition(dim,toList allcells))
    )

maxCells = method()
maxCells(CellComplex) := (cacheValue (symbol maxCells)) (cellComplex ->
    (
        lst := flatten values cells cellComplex;
        if #lst == 0 then return new HashTable;
        bdfn := c -> set boundaryCells c;
        maxcells := set lst;
        bdcells := sum (maxcells/bdfn);
        maxcells = maxcells - bdcells;
        while #bdcells != 0 do (
            bdcells = sum (bdcells/bdfn);
	    maxcells = maxcells - bdcells;
	    );
        partition(dim,toList maxcells)
        ))

--Define dimension for cell
dim(Cell) := (cell) -> cell.cellDimension

--Define dimension for cell complex 
dim(CellComplex) := (cellComplex) -> max keys cellComplex.cells

--Define ring for cell complex 
ring(CellComplex) := (cellComplex) -> cellComplex.labelRing


cellLabel = method()
cellLabel(Cell) := (cell) -> cell.label

--Make a cell, internal function
makeCell := (lst, l, d) -> (
    bdim := -1;
    for cell in lst do ( 
	if bdim < 0 
	then bdim = dim cell#0
	else assert(bdim == dim cell#0)
	);
    n := max(bdim + 1,d);
    new Cell from {
	symbol cellDimension => n, 
	symbol boundary => lst, -- could verify that it's a list
	symbol label => l
    	}
    );

chainToVirtualTally := (lst) -> (
    if lst == {}
    then new VirtualTally from {}
    else sum(lst, (cell,deg) -> new VirtualTally from {cell => deg})
    )

boundary(Cell) := (cell) -> cell.boundary
boundaryCells = method()
boundaryCells(Cell) := (cell) -> apply(boundary(cell), c -> first c)
--Boundary function, returns the boundary as a VirtualTally
boundaryTally := (cell) -> chainToVirtualTally cell.boundary


--Check if a chain, represented by a list is a boundary
isCycle = method()
isCycle(List) := (lst) ->
    ((sum(lst,l -> (
                c := l#0;
                deg := l#1;
                if deg>0
                then sum(deg,i -> boundaryTally c)
                else - sum(deg,i -> boundaryTally c)))) ? 0) == symbol ==


--Figure out an orientation automatically
inferOrientation := (lst) -> (
    if #lst == 2 and (dim first lst) == 0 then (
        ret := {(lst#0,1),(lst#1,-1)};
        if not isCycle ret then error "The given list of cells do not form a cycle";
        return ret
        );
    boundaryChain := new VirtualTally from {};
    remainingCells := lst;
    --the "?" based comparison is a work arround for "==" not working correctly for VirtualTally==ZZ
    while (boundaryChain ? 0) != symbol == or #remainingCells!=0 list (
        if (boundaryChain ? 0) == symbol ==
        then (
            if remainingCells =!= lst then error "The orientation on the cycle is non-unique";
            nextCell := last remainingCells;
            remainingCells = drop(remainingCells,-1);
            boundaryChain = boundaryTally nextCell;
            (nextCell,1)
            )
        else (
            c := (keys boundaryChain)#0;
            nextElems := select(remainingCells,c2 -> (boundaryTally c2)#?c);
            if #nextElems==0 then error "The given list of cells do not form a cycle";
            newBoundaryComponent := boundaryTally (nextElems#0);
            remainingCells = delete(nextElems#0,remainingCells);--Inefficient
            --check sign equality
            if (boundaryChain#c)*(newBoundaryComponent#c)<0
            then (
                boundaryChain = boundaryChain + boundaryTally (nextElems#0);
                (nextElems#0,1)
                )
            else (
                boundaryChain = boundaryChain - boundaryTally (nextElems#0);
                (nextElems#0,-1)
                )
            )
        )
    )

--Convert it to a submodule of R^1 if possible
toModule := (R,x) -> (
    if instance(x,Module) then return x;
    if instance(x,Ideal) then return module x;
    if instance(x,RingElement) then return image matrix {{x}};
    if instance(x,Number) then return image matrix {{x_R}};
    error "Expected a Module, Ideal, RingElement, or Number"
    )

inferLabel := boundary -> (
    if boundary == {} then return 1;
    if instance(boundary#0,Sequence) then return boundary/first//inferLabel;
    if all(boundary/cellLabel,b -> instance(b,RingElement) or instance(b,Number))
    then boundary/cellLabel//lcm
    else (
        rings := select(boundary/cellLabel, b -> not (instance(b,RingElement) or instance(b,Number)))/ring;
        nonNumberRings := select(rings,r -> ancestor(Number,r));
        R := if nonNumberRings==={} then rings#0 else nonNumberRings#0;
        boundary/cellLabel/(x -> toModule(R,x))//intersect
        )
    )

--Attach a cell
newCell = method(Options => {CellDimension=>null})
newCell(List,Thing) := opt -> (boundary,label) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return newCell(inferOrientation boundary,label,CellDimension=>opt.CellDimension);
    if not isCycle boundary then error "Expected the boundary to be a cycle";
    cd := if opt.CellDimension=!=null then opt.CellDimension else 0;
    c := makeCell(boundary,label,cd);
    if opt.CellDimension=!=null and dim c > cd then error "Incorrect CellDimesion optional parameter";
    c
    )
newCell(List) := opt -> cells -> newCell(cells,inferLabel cells,CellDimension=>opt.CellDimension);



isSimplexBoundary := (lst) -> (
    if #lst==0 then return true;
    bdim := dim first lst#0;
    all(lst,isSimplex @@ first) and
    all(lst,i -> dim first i == bdim) and
    (#lst == bdim+2) and
    (length lst == length unique (lst/first)) and
    (isCycle lst)
    )

isSimplex = method();
isSimplex(Cell) := cell ->
     isSimplexBoundary boundary cell

newSimplexCell = method();
newSimplexCell(List) := (boundary) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return newSimplexCell inferOrientation boundary;
    if not isSimplexBoundary boundary then error "The given boundary is not a valid boundary for a simplex";
    newCell boundary
    )
newSimplexCell(List,Thing) := (boundary,label) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return newSimplexCell(inferOrientation boundary,label);
    if not isSimplexBoundary boundary then error "The given boundary is not a valid boundary for a simplex";
    newCell(boundary,label)
    )

--Relabel function 
relabelCellComplex = method(Options=>{InferLabels=>true});
relabelCellComplex(CellComplex,HashTable) := o -> (C,T) -> (
    dimC := dim C;
    R := ring C;
    tablecellsbydim := for i to dimC list select(keys T, c -> dim c == i);
    relabeledcells := new MutableHashTable;
    for c in cells(0,C) do relabeledcells#c = (
	if any(tablecellsbydim#0, cell -> cell === c) then newCell({},T#c)
	else newCell({},cellLabel c)
	);
    for i from 1 to dimC do (
	for c in cells(i,C) do (
	    newbd := for b in boundaryCells c list relabeledcells#b;
	    newlabel := if any(tablecellsbydim#i, cell -> cell === c) then T#c 
	    else if not o.InferLabels then cellLabel(c)
	    else inferLabel(newbd);
	    relabeledcells#c = newCell(newbd, newlabel);
	    );
	);
    cellComplex(R, flatten values relabeledcells)
    )

RingMap ** CellComplex := (f,c) -> (
    if source f =!= ring c then error "source ring should match label ring";
    R := source f;
    S := target f;
    allCells := flatten values cells(c);
    -- ht := hashTable apply(allCells, c -> (c,f ** toModule(R,cellLabel c)));
    ht := hashTable apply(allCells, c -> (c,f(cellLabel c)));
    tempCellComplex := relabelCellComplex(c,ht);
    cellComplex(S,flatten values cells tempCellComplex) 
    )

--Get list of cells 
cells = method();
cells(CellComplex) := (cellComplex) -> cellComplex.cells
cells(ZZ,CellComplex) := (r,cellComplex) -> (
    if cellComplex.cells#?r 
    then cellComplex.cells#r
    else {}
    )

-- skeleton = method();
skeleton(ZZ,CellComplex) := (n,cellComplex) -> (
    c := new HashTable from select(pairs cellComplex.cells, (k,v) -> k<=n);
    mkCellComplex(cellComplex.labelRing,c,null)
    )

--take a hash table of RingElements/Matrices, and make a matrix, or 0
sparseBlockMatrix := (ht) -> (
    ks := keys ht;
    if ks === {} then return 0;
    rows := max (ks/first) + 1;
    columns := max (ks/(p->p#1)) + 1;
    maybeHt := p -> (
        if ht#?p then ht#p else 0
        );
    matrix apply(rows,i -> apply(columns, j -> maybeHt(i,j))))

--Create chain complex from cell complex 
boundary(ZZ,CellComplex) := (r,cellComplex) -> (
    R := cellComplex.labelRing;
    t := r-1;
    rCells := cells(r,cellComplex);
    tCells := cells(t,cellComplex);
    --We define these tables in two steps so that the ordering of the modules in domain and codomain
    --is consistent, especially between calls to boundary with different values for "r".
    domainModules := apply(toList rCells, c-> (c,toModule(R,cellLabel c)));
    codomainModules := apply(toList tCells, c -> (c,toModule(R,cellLabel c)));
    domainModulesTable :=
        new HashTable from domainModules;
    codomainModulesTable :=
        new HashTable from codomainModules;
    domain := fold((a,b) -> a ++ b, R^0, apply(domainModules,last));
    codomain := if t==-1 then R^1 else fold((a,b) -> a ++ b, R^0, apply(codomainModules,last));
    tCellsIndexed := new HashTable from toList apply(pairs(tCells),reverse);
    i := 0;
    L := flatten for F in rCells list (
	l := if t==-1
             then (0,i) => inducedMap(codomain,domainModulesTable#F)
             else for p in pairs boundaryTally F list(
                 (cell,deg) := p;
                 if dim cell < dim F - 1 then continue;
                 (tCellsIndexed#cell,i) => deg_R*inducedMap(codomainModulesTable#cell,domainModulesTable#F));
	i = i+1;
	l
	);
    map(codomain,domain,sparseBlockMatrix new HashTable from L)
    );

chainComplex(CellComplex) := {Reduced=>true} >> o -> (cellComplex) -> (
    if not cellComplex.cache.?chainComplex then (
        cellComplex.cache.chainComplex =
            (chainComplex apply(max((dim cellComplex) + 1,1), r -> boundary(r,cellComplex)))[1]
        );
    if not o.Reduced then (
	Ccopy := chainComplex cellComplex.cache.chainComplex;
	Ccopy_(-1) = 0*Ccopy_(-1);
	Ccopy
	)
    else cellComplex.cache.chainComplex
    );

--Get homology directly from cell complex
homology(ZZ,CellComplex) := opts -> (i,cellComplex) -> (
    homology_i chainComplex(cellComplex)
    );

homology(CellComplex) := opts -> (cellComplex) -> (
    homology chainComplex(cellComplex)
    );

--Get cohomology directly from cell complex
cohomology(ZZ,CellComplex) := opts -> (i,cellComplex) -> ( 
    cohomology_i Hom(chainComplex(cellComplex),cellComplex.labelRing^1)
    );

----------
---Here there be polyhedra 
----------

faces(Polyhedron) := opts -> (P) -> Polyhedra$faces P
faces(ZZ,Polyhedron) := opts -> (r,P) -> Polyhedra$faces(r,P)
vertices(Polyhedron) := (P) -> Polyhedra$vertices P

faces(PolyhedralComplex) := opts -> (PC) -> Polyhedra$faces PC 
faces(ZZ,PolyhedralComplex) := opts -> (r,PC) -> Polyhedra$faces(r,PC) 
vertices(PolyhedralComplex) := (PC) -> Polyhedra$vertices PC

cellComplex(Ring,Polyhedron) := (R,P) -> (
    if not isCompact P then error "The given polyhedron is not compact.";
    Pdim := dim P;
    Pfaces := applyPairs(faces P, (i,lst) -> (Pdim-i,apply(lst,first)));
    cells := new MutableHashTable;
    for i from 0 to Pdim do (
        for face in Pfaces#i do (
            bd := if i!=0
                  then for f in Pfaces#(i-1) list (if isSubset(f,face) then cells#f else continue)
                  else {};
            cells#face = newCell bd
            );
        );
    cellComplex(R,flatten values cells)
    );

cellComplex(Ring,PolyhedralComplex) := (R,P) -> (
    Pdim := dim P;
    Pfaces := applyPairs(faces P, (i,lst) -> (Pdim-i-1,apply(lst,first)));
    cells := new MutableHashTable;
    for i from 0 to Pdim do (
        for face in Pfaces#i do (
	    bd := if i!=0
                  then for f in Pfaces#(i-1) list (if isSubset(f,face) then cells#f else continue)
                  else {};
            cells#face = newCell bd
            );
        );
    cellComplex(R,flatten values cells)
   );

-------------
-- Posets
-------------

facePoset(CellComplex) := (cellComplex) -> (
    G := flatten values cells cellComplex;
    contain := (a,b) -> member(a,boundaryCells b) or a === b;-- a contained or equal b
    P := poset(G,contain);
    rel := allRelations P;
    M := transitiveClosure(G,rel);
    poset(G,rel,M)
    )

-------------
-- Minimality
-------------

isFree = method(TypicalValue => Boolean);
--check if all the labels are free modules
isFree(CellComplex) := (cellComplex) -> (
    R := cellComplex.labelRing;
    all(flatten values cells cellComplex,c -> isFreeModule prune toModule(R,cellLabel c))
    )

isCellMinimal := (R,cell) -> (
    label := toModule(R,cellLabel cell);
    all(boundary cell, c -> toModule(R,cellLabel first c) != label)
    )

isMinimal = method(TypicalValue => Boolean)
--Check if a labeled cell complex is minimal, Note: we assume the cell complex is free (see isFree)
isMinimal(CellComplex) := (cellComplex) -> (
    R := cellComplex.labelRing;
    all(flatten values cells cellComplex,c -> isCellMinimal(R,c))
    )

---------
-- Output
---------

net(Cell) := (cell) -> (
    "Cell of dimension " | (dim cell)
    )

net(CellComplex) := (cellComplex) -> (
    d := dim cellComplex;
    nTotalCells := #(flatten values cells cellComplex);
    if nTotalCells == 0
    then "empty CellComplex"
    else "CellComplex of dimension " | d | " with " | nTotalCells | " total cells"
    ); 


------------------------
-- Common cell complexes
------------------------

cellComplexSphere = method();
cellComplexSphere(Ring,ZZ) := (R,n) -> (
    if n<0 then error "cellComplexSphere expects a non-negative integer";
    v := newSimplexCell {};
    if n==0 then (
        w := newSimplexCell {};
        cellComplex(R,{v,w})
        )
    else(
        c := newCell({(v,0)},CellDimension=>n);
        cellComplex (R,{c})
        )
    )

cellComplexRPn = method();
cellComplexRPn(Ring,ZZ) := (R,n) -> (
    if n<0 then error "cellComplexRPn expects a non-negative integer";
    t := newSimplexCell {};
    if n==0 then return cellComplex(R,{t});
    for i from 1 to n do(
        attachingDegree := if even i then 2 else 0;
        t = newCell {(t,attachingDegree)};
        );
    cellComplex(R,{t})
    )

cellComplexTorus = method();
cellComplexTorus(Ring,ZZ) := (R,n) -> (
    if n<0 then error "cellComplexTorus expects a non-negative integer";
    v := newSimplexCell {};
    if n==0 then return cellComplex(R,{v});
    cells := new MutableHashTable;
    for s in subsets(n) do (
        k := #s;
        cells#s = newCell apply(subsets(s,k-1), s' -> (cells#s',0))
        );
    cellComplex(R,{cells#(toList (0..(n-1)))})
    )


----------------------------
-- Specific chain complexes
----------------------------

taylorComplex = method();
taylorComplex(MonomialIdeal) := (I) -> (
    gensI := I_*;
    r := #gensI;
    if r == 0 then error "taylorComplex expects a non-zero monomialIdeal";
    cells := new MutableHashTable;
    for i to r-1 do cells#{i} = newSimplexCell({},gensI#i);
    for k from 2 to r do (
	for s in subsets(r,k) do (
	    bd := for t in subsets(s,k-1) list cells#t;
	    cells#s = newSimplexCell(bd, lcm(gensI_s));
	    );
	);
    cellComplex(ring I, {cells#(splice {0..(r-1)})})
    )

hullComplex = method();
hullComplex(MonomialIdeal) := (I) -> ( 
    gensI := I_*;
    R := ring I;
    n := #(gens R);
    t := (n+1)! + 1;
    expvecs := flatten (gensI/exponents);
    verts := for a in expvecs list for i from 0 to (n-1) list t^(a#i);
    P := convexHull(transpose matrix verts, id_(ZZ^n));
    Pdim := dim P;
    Pfaces := new HashTable from select(pairs faces P, (k,v) -> (k <= Pdim and k > 0)); --weird inequalities bc codim
    Pfaces = applyValues(Pfaces, v -> apply(v,p -> if p#1 == {} then p#0)); -- selecting compact faces
    Pfaces = applyValues(Pfaces, v -> delete( , v));
    Pfaces = applyPairs(Pfaces, (d,lst) -> (Pdim-d,lst)); --flipping from codim to dim
    cells := new MutableHashTable;
    for v in Pfaces#0 do cells#v = newCell({},gensI#(v#0));
    for i from 1 to Pdim-1 do (
        for face in Pfaces#i do (
            bd := for f in Pfaces#(i-1) list (if isSubset(f,face) then cells#f else continue);
            cells#face = newCell bd
            );
        );
    cellComplex(R,flatten values cells)
    )

----------------------------

----------------------------


load "./CellularResolutions/doc.m2"
load "./CellularResolutions/tests.m2"

end

restart
installPackage("CellularResolutions")
loadPackage("CellularResolutions", Reload => true)
check(CellularResolutions)
viewHelp CellularResolutions

