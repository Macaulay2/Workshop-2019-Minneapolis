-- -*- coding: utf-8 -*-
newPackage(
    "CellularResolutions",
    Version => "0.1",
    Date => "July 22, 2019",
    Authors => {
        {Name => "Jay Yang", Email => "jkyang@umn.edu"},
        {Name => "Aleksandra Sobieska", Email => "ola@math.tamu.edu"}
        },
    Headline => "A package for cellular resolutions of monomial ideals",
    AuxiliaryFiles => true, -- set to true if package comes with auxiliary files
    PackageExports => {"Polyhedra", "SimplicialComplexes"}
    )

export {"CellComplex",
        "Cell",
        "cellComplex",
        "makeCell",
        "attach",
        "isCycle",
        "attachSimplex",
        "isSimplex",
	"cells",
        "cellLabel",
        "newCell",
        "newSimplexCell",
        "neg1Cell",
        "isFree",
        "isMinimal"
--        "skeleton"
        }
protect labelRing
protect cellDimension

CellComplex = new Type of HashTable
--Note, the mutable hash table means that equality works "Correctly"
Cell = new Type of MutableHashTable

--Creates a single negative 1 dimensional cell 
neg1Cell = new Cell from {
    symbol cellDimension => -1,
    symbol boundary => {},
    symbol label => 1 
    };

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

--Adds a single -1 cell by default
--TODO: create an option to make a void complex
cellComplex = method()
cellComplex(Ring,List) := (R,maxCells) -> (
    new CellComplex from {
    	symbol labelRing => R,
        symbol cells => cellsFromMaxCells maxCells
	}
    )
cellComplex(Ring,SimplicialComplex) := (R,C) -> (
    S := ring C;
    Cfaces := applyValues(new HashTable from faces C,flatten @@ entries);
    cells := new MutableHashTable from {1_S => neg1Cell};
    for i from 0 to dim C do (
        for simplex in Cfaces#i do (
            bd := for x in gens S list (if simplex%x==0 then cells#(simplex//x) else continue);
            cells#simplex = newCell bd
            );
        );
    cellComplex(R,values cells)
    )

--Define dimension for cell
dim(Cell) := (cell) -> cell.cellDimension

--Define dimension for cell complex 
dim(CellComplex) := (cellComplex) -> max keys cellComplex.cells

--Define ring for cell complex 
ring(CellComplex) := (cellComplex) -> cellComplex.labelRing

cellLabel = method()
cellLabel(Cell) := (cell) -> cell.label

--Make cell 
makeCell := (lst, l) -> (
    bdim := -1;
    for cell in lst do ( 
	if bdim < 0 
	then bdim = dim cell#0
	else assert(bdim == dim cell#0)
	);
    n := bdim + 1;
    bd := if lst == {} then {(neg1Cell,1)} else lst;
    new Cell from {
	symbol cellDimension => n, 
	symbol boundary => bd, -- could verify that it's a list
	symbol label => l
    	}
    );

chainToVirtualTally := (lst) -> (
    if lst == {}
    then new VirtualTally from {}
    else sum(lst, (cell,deg) -> new VirtualTally from {cell => deg})
    )

boundary(Cell) := (cell) -> cell.boundary
--Boundary function, returns the boundary as a VirtualTally
boundaryTally := (cell) -> chainToVirtualTally cell.boundary

--Check if a chain, represented by a list is a boundary
isCycle = method()
isCycle(List) := (lst) ->
    (sum(lst,(c,deg) -> if deg>0
                        then sum(deg,i -> boundaryTally c)
                        else - sum(deg,i -> boundaryTally c)) ? 0) == symbol ==


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
newCell = method()
newCell(List,Thing) := (boundary,label) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return newCell(inferOrientation boundary,label);
    if not isCycle boundary then error "Expected the boundary to be a cycle";
    c := makeCell(boundary,label);
    c
    )
newCell(List) := (cells) ->
    newCell(cells,inferLabel cells)

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



--Get list of cells 
cells = method();
cells(CellComplex) := (cellComplex) -> cellComplex.cells
cells(ZZ,CellComplex) := (r,cellComplex) -> (
    if cellComplex.cells#?r 
    then cellComplex.cells#r
    else {}
    )

--TODO polyhedra also defines skeleton
-- skeleton = method();
skeleton(ZZ,CellComplex) := (n,cellComplex) -> (
    c := new HashTable from select(pairs cellComplex.cells, (k,v) -> k<=n);
    new CellComplex from {
        symbol labelRing => cellComplex.labelRing,
        symbol cells => c
	}
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
    domainModules :=
        new HashTable from apply(toList rCells, c-> (c,toModule(R,cellLabel c)));
    codomainModules :=
        new HashTable from apply(toList tCells,c -> (c,toModule(R,cellLabel c)));
    domain := fold((a,b) -> a ++ b, R^0, values domainModules);
    codomain := fold((a,b) -> a ++ b, R^0, values codomainModules);
    tCellsIndexed := new HashTable from toList apply(pairs(tCells),reverse);
    i := 0;
    L := flatten for F in rCells list (
	l := apply(pairs boundaryTally F,
            (cell,deg) -> (tCellsIndexed#cell,i) => deg_R*inducedMap(codomainModules#cell,domainModules#F));
	i = i+1;
	l
	);
    map(codomain,domain,sparseBlockMatrix new HashTable from L)
    );

chainComplex(CellComplex) := (cellComplex) -> (
    (chainComplex apply((dim cellComplex) + 1, r -> boundary(r,cellComplex)))[1]
    );

--Get homology directly from cell complex
homology(ZZ,CellComplex) := opts -> (i,cellComplex) -> (
    homology_i chainComplex cellComplex
    );

homology(CellComplex) := opts -> (cellComplex) -> (
    homology chainComplex cellComplex
    );

--Get cohomology directly from cell complex
cohomology(ZZ,CellComplex) := opts -> (i,cellComplex) -> ( 
    cohomology_i Hom(chainComplex cellComplex,cellComplex.labelRing^1)
    );

----------
---Here there be polyhedra 
----------

faces(Polyhedron) := opts -> (P) -> Polyhedra$faces P
faces(ZZ,Polyhedron) := opts -> (r,P) ->Polyhedra$faces(r,P)
vertices(Polyhedron) := (P) -> Polyhedra$vertices P

cellComplex(Ring,Polyhedron) := (R,P) -> (
    if not isCompact P  then error "The given polyhedron is not compact.";
    Pdim := dim P;
    Pfaces := applyKeys(faces P, i -> Pdim-i); --flips from codim to dim
    Pfaces = applyValues(Pfaces, lst -> if lst != {} then apply(lst, lst -> first lst) else {});
    cells := new MutableHashTable from {{} => neg1Cell};
    for i from 0 to Pdim do (
        for face in Pfaces#i do (
	    bd := for f in Pfaces#(i-1) list (if isSubset(f,face) then cells#f else continue);
            cells#face = newCell bd
            );
        );
    cellComplex(R,flatten values cells)    
    );

-------------
-- Minimality
-------------

isFree = method(TypicalValue => Boolean);
--check if all the labels are free modules
isFree(CellComplex) := (cellComplex) -> (
    R := cellComplex.labelRing;
    all(cells cellComplex,c -> isFree toModule(R,cellLabel c))
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

net(CellComplex) := (cellComplex) -> (
    d := dim cellComplex;
    nMaxCells := #(cells(d,cellComplex));
    nTotalCells := #(cells cellComplex);
    "CellComplex of dimension " | d | " with " | nMaxCells | " maximal cells and " | nTotalCells | " total cells"
    );

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

