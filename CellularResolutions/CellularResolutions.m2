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
	AuxiliaryFiles => false -- set to true if package comes with auxiliary files
    	)

export {"CellComplex",
        "Cell",
        "cellComplex",
        "makeCell",
        "attach",
        "isCycle",
        "attachSimplex",
        "isSimplex"
        }
protect labelRing
protect cells 
protect boundary 
protect cellDimension
protect label

CellComplex = new Type of HashTable
--Note, the mutable hash table means that equality works "Correctly"
Cell = new Type of MutableHashTable

--TODO: Add an option to allow/disallow a single -1 cell
cellComplex = method()
cellComplex(Ring) := (R) -> (
    new CellComplex from {
    	symbol labelRing => R,
    	symbol cells => new MutableHashTable from {}}
    )

--Define dimension for cell
dim(Cell) := (cell) -> cell.cellDimension

--Define dimension for cell complex 
dim(CellComplex) := (cellComplex) -> max keys cellComplex.cells

--Make cell 
makeCell := (lst, l) -> (
    bdim := -1;
    for cell in lst do ( 
	if bdim < 0 
	then bdim = dim cell#0
	else assert(bdim == dim cell#0)
	);
    n := bdim + 1;
    new Cell from {
	symbol cellDimension => n, --FIGURE IT OUT
	symbol boundary => lst, -- could verify that it's a list
	symbol label => l
    	}
    );

chainToVirtualTally := (lst) -> (
    sum(lst, (cell,deg) -> new VirtualTally from {cell => deg})
    )

boundary := (cell) -> cell.boundary
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
            nextElems := select(lst,c2 -> ((boundaryTally c2)#c)!=0);
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

--Attach a cell
attach = method()
attach(CellComplex,List,Thing) := (baseComplex,boundary,label) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return attach(baseComplex,inferOrientation boundary,label);
    c := makeCell(boundary,label);
    if not isCycle boundary then error "Expected the boundary to be a cycle";
    n := dim c;
    if baseComplex.cells#?n 
    then(
        i := #baseComplex.cells#n;
        baseComplex.cells#n#i=c;
	)
    else (
	baseComplex.cells#n = new MutableList from {c};
	);
    c
    )
attach(CellComplex,List) := (baseComplex,cells) -> attach(baseComplex,cells,1)

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

attachSimplex = method();
attachSimplex(CellComplex,List) := (baseComplex,boundary) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return attachSimplex(baseComplex,inferOrientation boundary);
    if not isSimplexBoundary boundary then error "The given boundary is not a valid boundary for a simplex";
    attach(baseComplex,boundary)
    )

--Get list of cells 
cells := (cellcomplex) -> toList cellcomplex.cells

----------------------------

-- THINGS TO DO:
-- - verify that boundary is a cycle 

----------------------------


beginDocumentation()

doc ///
    Key
        CellularResolutions
    Headline
        A package for cellular resolutions of monomial ideals
    Description
        Text
            Cellular Resolutions
///

doc ///
    Key
        CellComplex
    Headline
        the class of all cell complexes
    Description
        Text
            A cell complex in this context is the combinatorial data of a
            CW-complex, i.e. a collection of cells in various dimensions along
            with their boundary expressed as a sequence of cells along with an
            orientation such that the boundary is a cycle.
///

doc ///
    Key
        cellComplex
        (cellComplex,Ring)
    Headline
        create an empty cell complex
    Usage
        cellComplex R
    Inputs
        R : Ring
            that specifies the base ring for the labels
    Outputs
        : CellComplex
            an empty cell complex
    Description
        Text
            This constructs an empty cell complex to which cells can be attached
        Example
            C = cellComplex QQ;
            R = QQ[x,y];
            C = cellComplex R;
    SeeAlso
        (attach,CellComplex,List)
        (attach,CellComplex,List,Thing)
///

doc ///
    Key
        attach
        (attach,CellComplex,List,Thing)
        (attach,CellComplex,List)
    Headline
        attach a cell to a cell complex
    Usage
        attach(C,boundary,label)
        attach(C,boundary)
    Inputs
        C : CellComplex
            the complex to attach a cell to
        boundary : List
            that gives the boundary of the cell to attach as a list of pairs of cells and attaching degree
        label : Thing
            that gives a label to associate to the cell, default of 1
    Outputs
        : Cell
            that was attached
///

TEST ///
C = cellComplex(QQ);
v1 = attachSimplex(C,{});
v2 = attachSimplex(C,{});
assert(isSimplex v1);
assert(isSimplex v2);
assert(dim C==0);
assert(dim v1==0);
assert(dim v2==0);
l1 = attachSimplex(C,{(v1,1),(v2,-1)});
l2 = attachSimplex(C,{v1,v2});
assert(isSimplex l1);
assert(isSimplex l2);
assert(dim C==1);
assert(dim l1==1);
assert(dim l2==1);
f1 = attach(C,{l1,l2});
assert(dim C==2);
assert(dim f1==2);
///


end

restart
installPackage("CellularResolutions")
loadPackage("CellularResolutions", Reload => true)
check(CellularResolutions)
viewHelp CellularResolutions

C = cellComplex()
a = attach(C,{});
b = attach(C,{});

C = cellComplex(R);
a = attach0Cell(C,"a");
b = attach0Cell(C,"b");

C = cellComplex({"a","b"})


a = make0Cell("a");
b = make0Cell("b");
C = cellComplex({a,b});

l1 = attach(C,{a,b});
l2 = attach(C,{a,b});
attach(C,{l1,l2});


