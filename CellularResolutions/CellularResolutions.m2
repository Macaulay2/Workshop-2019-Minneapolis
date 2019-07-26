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
    AuxiliaryFiles => false, -- set to true if package comes with auxiliary files
    PackageExports => {"SimplicialComplexes"}
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
        "neg1Cell"
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
	symbol cellDimension => n, --FIGURE IT OUT
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

----------------------------

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
        Cell
    Headline
        the class of all cells in cell complexes
    Description
        Text
            This class represents a single cell in a cell complex.
    SeeAlso
        CellComplex
///

doc ///
    Key
        cellComplex
        (cellComplex,Ring,List)
    Headline
        create an cell complex
    Usage
        cellComplex(R,maxCells)
    Inputs
        R : Ring
            that specifies the base ring for interpreting the labels
        naxCells : List
            that specifies the maximal cells for the complex
    Outputs
        : CellComplex
            a cell complex
    Description
        Text
            This function constructs a cell complex from a list of cells
        Example
            c = newCell {}
            C = cellComplex(QQ,{c});
            R = QQ[x,y];
            C = cellComplex(R,{c});
    SeeAlso
        (newCell,List)
        (newCell,List,Thing)
///

doc ///
    Key
        newCell
        (newCell,List,Thing)
        (newCell,List)
    Headline
        creates a new cell
    Usage
        newCell(boundary,label)
        newCell(boundary)
    Inputs
        boundary : List
            that gives the boundary of the new cell either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was created
    Description
        Text
            This function creates a new cell, to be added to a cell complex,
            if given a list of cells without any orientation information, it
            attempts to infer the orientation.
        Text
            If given an empty list for the boundary, the function creates a
            0-cell (a vertex).
        Text
            If not given a label, and the labels on the boundary are monomials
            or monomial ideals, from the same ring, then the label is the lcm
            of the labels of the boundary
        Example
            R = QQ[x,y]
            a = newCell({},x);
            b = newCell({},y);
            c1 = newCell {(a,1),(a,-1)};
            c2 = newCell {a,a};
            c3 = newCell {a,b};
            C = cellComplex(R,{c1,c2,c3});
    Caveat
        This function does not check that there is a valid map from the boundary
        of a n-cell to the given boundary. It only checks that the boundary
        forms a cycle in homology.
    SeeAlso
        cellComplex
        newSimplexCell
///

doc ///
    Key
        newSimplexCell
        (newSimplexCell,List,Thing)
        (newSimplexCell,List)
    Headline
        create a new cell
    Usage
        newSimplexCell(boundary,label)
        newSimplexCell(boundary)
    Inputs
        boundary : List
            that gives the boundary of the new cell either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was created
    Description
        Text
            This function will only create simplices, and it will verify that
            the new cell is a simplex, as such does not have the caveat of
            @TO newCell@. Otherwise it has the same behavior. This is
            particularly useful in constructing \Delta-complexes.
    SeeAlso
        cellComplex
        newCell
///

doc /// 
    Key 
    	cellLabel 
	(cellLabel,Cell)
    Headline 
    	return the label of a cell 
    Usage 
    	cellLabel C
    Inputs
    	C : Cell 
    Outputs
    	: Thing 
	    the label of the cell
    SeeAlso
    	Cell
///

doc ///
    Key
        (dim,CellComplex)
    Headline
        compute the dimension of a cell complex
    Usage
        dim C
    Inputs
        C : CellComplex
            the complex to compute the dimension of
    Outputs
        : ZZ
            the dimension of the complex
    SeeAlso
        (dim,Cell)
///

doc ///
    Key
        (dim,Cell)
    Headline
        compute the dimension of a cell
    Usage
        dim C
    Inputs
        C : Cell
            the cell to compute the dimension of
    Outputs
        : ZZ
            the dimension of the cell
    SeeAlso
        (dim,CellComplex)
///

doc ///
    Key 
    	(cells,CellComplex)
    Headline 
    	return the cells of a cell complex
    Usage
    	cells C
    Inputs 
    	C : CellComplex  
	    the CellComplex whose cells are to be returned 
    Outputs 
    	: MutableHashTable 
	    the cells of C 
    SeeAlso 
    	(cells,ZZ,CellComplex)
    	
///

doc /// 
    Key
    	(cells,ZZ,CellComplex)
    Headline 
    	return the cells of a cell complex 
    Usage 
    	cells(ZZ,C)
    Inputs 
    	r : ZZ
	    the dimension of the cells to be returned 
    	C : CellComplex 
	    the CellComplex whose r-cells are to be returned
    Outputs
    	: MutableList 
	    the r-cells of C
///

doc /// 
    Key 
    	boundary 
	(boundary, ZZ, CellComplex) 
    Headline 
    	compute the boundary map of a cell complex from r-faces to (r-1)-faces 
    Usage 
    	boundary(ZZ,CellComplex) 
    Inputs 
    	r : ZZ
	    the source cell-dimension 
	C : CellComplex 
	    the CellComplex 
    Outputs 
    	: Matrix 
	    the boundary map from r-faces to (r-1)-faces of C
    SeeAlso 
    (chainComplex,CellComplex)
    (boundary,SimplicialComplex)
    (chainComplex,SimplicialComplex) 
///

doc ///
    Key
        (chainComplex,CellComplex)
    Headline
        compute the cellular chain complex for a cell complex
    Usage
        chainComplex C
    Inputs
        C : CellComplex
            the complex to compute the chain complex
    Outputs
        :ChainComplex
            the dimension of the complex
    Description
        Text
            This constructs the cellular chain complex for a cell complex,
            taking into account the labels on the cells
        Example
            R = QQ[x]
            a = newSimplexCell({},x);
            b1 = newCell {a,a};
            b2 = newCell {a,a};
            C = cellComplex(R,{b1,b2});
            chainComplex C
    SeeAlso
        (boundary,ZZ,CellComplex)
	(boundary,SimplicialComplex)
        (chainComplex,SimplicialComplex)
///

TEST ///
assert(dim cellComplex(QQ,{}) === -infinity);
assert(dim cellComplex(QQ,{neg1Cell}) === -1);
///


TEST ///
v1 = newSimplexCell {};
v2 = newSimplexCell {};
assert(isSimplex v1);
assert(isSimplex v2);
assert(dim v1===0);
assert(dim v2===0);
assert(v1 =!= v2);
assert(cellLabel v1 === 1)
l1 = newSimplexCell {(v1,1),(v2,-1)};
l2 = newSimplexCell {v1,v2};
assert(isSimplex l1);
assert(isSimplex l2);
assert(l1=!=l2);
assert(dim l1==1);
assert(dim l2==1);
C = cellComplex(QQ,{l1,l2});
CchainComplex = chainComplex C;
assert(HH_0(CchainComplex)==0);
assert(prune HH_1(CchainComplex)==QQ^1);
assert(HH_2(CchainComplex)==0);
f1 = newCell {l1,l2};
C = cellComplex(QQ,{f1});
assert(dim C==2);
assert(dim f1==2);
assert(cellLabel f1 == 1);
delneg1 = boundary(-1,C);
del0 = boundary(0,C);
del1 = boundary(1,C);
del2 = boundary(2,C);
del100 = boundary(100,C);
assert(delneg1 == map(QQ^0,QQ^1,0));
assert(del0 == map(QQ^1,QQ^2, {{1,1}}));
--This doesn't work due to issues about the ordering of the cells
--assert(del1 == map(QQ^2,QQ^2, {{1,1},{-1,-1}}));
assert(del100 == map(QQ^0,QQ^0,{}));
CchainComplex = chainComplex C;
assert(HH_0(CchainComplex)==0);
assert(HH_1(CchainComplex)==0);
assert(HH_2(CchainComplex)==0);
///

TEST /// 
a = newCell {};
b1 = newCell {(a,1),(a,-1)};
b2 = newCell {(a,1),(a,-1)};
D = cellComplex(QQ[x],{b1,b2});
assert(dim D == 1);
assert(isSimplex a);
assert(not isSimplex b1);
assert(not isSimplex b2);
DchainComplex = chainComplex D;
assert(HH_0(DchainComplex)==0);
R = ring D;
assert(prune HH_1(DchainComplex)==R^2);
assert(HH_2(DchainComplex)==0);
///

--Koszul Complex via Taylor resolutions
TEST ///
R = QQ[x,y,z];
vx = newSimplexCell({},x);
vy = newSimplexCell({},y);
vz = newSimplexCell({},z);
lxy = newSimplexCell({vx,vy});
lyz = newSimplexCell({vy,vz});
lxz = newSimplexCell({vx,vz});
fxyz = newSimplexCell({lxy,lyz,lxz});
assert(cellLabel fxyz === x*y*z);
D = cellComplex(R,{fxyz});
C = (chainComplex D)[-1];
assert(HH_0(C)==cokernel matrix {{x,y,z}});
assert(C.dd^2==0);
///

--Monomial ideal labels
TEST ///
R = QQ[x,y,z];
vx = newSimplexCell({},ideal(x));
vy = newSimplexCell({},ideal(y));
vz = newSimplexCell({},ideal(z));
lxy = newSimplexCell({vx,vy});
lyz = newSimplexCell({vy,vz});
lxz = newSimplexCell({vx,vz});
fxyz = newSimplexCell({lxy,lyz,lxz});
D = cellComplex(R,{fxyz});
C = (chainComplex D)[-1];
assert(HH_0(C)==R^1/module ideal(x,y,z))
assert(HH_1(C)==0)
assert(C.dd^2==0);
///

--Non principal labels
TEST ///
R = QQ[x,y,z];
vx = newSimplexCell({},ideal(x,y));
vy = newSimplexCell({},ideal(y,z));
vz = newSimplexCell({},ideal(x,z));
lxy = newSimplexCell {vx,vy};
lyz = newSimplexCell {vy,vz};
lxz = newSimplexCell {vx,vz};
fxyz = newSimplexCell {lxy,lyz,lxz};
D = cellComplex(R,{fxyz});
C = (chainComplex D)[-1];
assert(C.dd^2==0);
///

TEST ///
R = QQ[w,x,y,z]
C = simplicialComplex monomialIdeal(w*x,w*y);
D = cellComplex(QQ,C);
assert(dim D==2);
assert(#cells(2,D)==1);
assert(#cells(1,D)==4);
assert(#cells(0,D)==4);
assert(#cells(-1,D)==1);
///



end

restart
installPackage("CellularResolutions")
loadPackage("CellularResolutions", Reload => true)
check(CellularResolutions)
viewHelp CellularResolutions

