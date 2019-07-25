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
        "cellLabel"
        }
protect labelRing
protect cellDimension

CellComplex = new Type of HashTable
--Note, the mutable hash table means that equality works "Correctly"
Cell = new Type of MutableHashTable

--Creates a single negative 1 dimensional cell 
neg1Cell := new Cell from {
    symbol cellDimension => -1,
    symbol boundary => {},
    symbol label => 1 
    };

--Adds a single -1 cell by default
--TODO: create an option to make a void complex
cellComplex = method()
cellComplex(Ring) := (R) -> (
    new CellComplex from {
    	symbol labelRing => R,
    	symbol cells => new MutableHashTable from 
	    {-1 => new MutableList from {neg1Cell}}
	}
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
attachSimplex(CellComplex,List,Thing) := (baseComplex,boundary,label) -> (
    if #boundary!=0 and instance(boundary#0,Cell)
    then return attachSimplex(baseComplex,inferOrientation boundary);
    if not isSimplexBoundary boundary then error "The given boundary is not a valid boundary for a simplex";
    attach(baseComplex,boundary,label)
    )



--Get list of cells 
cells = method();
cells(CellComplex) := (cellComplex) -> cellComplex.cells
cells(ZZ,CellComplex) := (r,cellComplex) -> (
    if cellComplex.cells#?r 
    then cellComplex.cells#r
    else {}
    )

--Create chain complex from cell complex 
boundary(ZZ,CellComplex) := (r,cellComplex) -> (
    R := cellComplex.labelRing;
    t := r-1;
    rCells := cells(r,cellComplex);
    tCells := cells(t,cellComplex);
    domain := R^(#rCells);
    codomain := R^(#tCells);
    tCellsIndexed := new HashTable from toList apply(pairs(tCells),reverse);
    i := 0;
    L := flatten for F in rCells list (
	l := apply(pairs boundaryTally F, (cell,deg) -> (tCellsIndexed#cell,i) => deg_R*(cellLabel F/cellLabel cell)_R);
	i = i+1;
	l
	);
    map(codomain,domain,L)
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
            that gives the boundary of the cell to attach either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was attached
    Description
        Text
            This function adds cells to a cell complex, if given a list of cells
            without any orientation information, then it attempts to infer the
            orientation.
        Text
            If given an empty list for the boundary, the function adds a 0-cell
            (a vertex).
        Text
            If not given a label, and the labels on the boundary are monomials
            or monomial ideals, from the same ring, then the label is the lcm
            of the labels of the boundary
        Example
            C = cellComplex(QQ[x,y]);
            a = attach(C,{},x);
            b = attach(C,{},y);
            c1 = attach(C,{(a,1),(a,-1)});
            c2 = attach(C,{a,a});
            c3 = attach(C,{a,b});
    Caveat
        This function does not check that there is a valid map from the boundary
        of a n-cell to the given boundary. It only checks that the boundary
        forms a cycle in homology.
    SeeAlso
        cellComplex
        attachSimplex
///

doc ///
    Key
        attachSimplex
        (attachSimplex,CellComplex,List,Thing)
        (attachSimplex,CellComplex,List)
    Headline
        attach a simplex to a cell complex
    Usage
        attachSimplex(C,boundary,label)
        attachSimplex(C,boundary)
    Inputs
        C : CellComplex
            the complex to attach a cell to
        boundary : List
            that gives the boundary of the cell to attach either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was attached
    Description
        Text
            This funciton will only attach simplicies, and it will verify that
            the attached cell is a simplex, as such does not have the caveat of
            @TO attach@. Otherwise it has the same behavior. This is
            particularly useful in constructing \Delta-complexes.
    SeeAlso
        cellComplex
        attach
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
        :ZZ
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
        :ZZ
            the dimension of the cell
    SeeAlso
        (dim,CellComplex)
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
            C = cellComplex(R);
            a = attachSimplex(C,{},x);
            b1 = attach(C,{a,a});
            b2 = attach(C,{a,a});
            chainComplex C
    SeeAlso
        (boundary,ZZ,CellComplex)
        (chainComplex,SimplicialComplex)
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
assert(cellLabel v1 == 1)
l1 = attachSimplex(C,{(v1,1),(v2,-1)});
l2 = attachSimplex(C,{v1,v2});
assert(isSimplex l1);
assert(isSimplex l2);
assert(dim C==1);
assert(dim l1==1);
assert(dim l2==1);
CchainComplex = chainComplex C;
assert(HH_0(CchainComplex)==0);
assert(prune HH_1(CchainComplex)==QQ^1);
assert(HH_2(CchainComplex)==0);
f1 = attach(C,{l1,l2});
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
assert(del1 == map(QQ^2,QQ^2, {{1,1},{-1,-1}}));
assert(del100 == map(QQ^0,QQ^0,{}));
CchainComplex = chainComplex C;
assert(HH_0(CchainComplex)==0);
assert(HH_1(CchainComplex)==0);
assert(HH_2(CchainComplex)==0);
///

TEST /// 
D = cellComplex(QQ[x]);
a = attach(D,{});
b1 = attach(D,{(a,1),(a,-1)});
b2 = attach(D,{(a,1),(a,-1)});
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


