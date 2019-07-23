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

export {"CellComplex","Cell","cellComplex","makeCell","attach"}
protect labelRing
protect cells 
protect boundary 
protect cellDimension 

CellComplex = new Type of HashTable

Cell = new Type of HashTable


cellComplex = method()
cellComplex(Ring) := (R) -> (
    new CellComplex from {
    	symbol labelRing => R,
    	symbol cells => new MutableHashTable from {}}
    )

--Define dimension for cell
dim(Cell) := (cell) -> cell.cellDimension

--Define dimension for cell complex 
dim(CellComplex) := (cellComplex) -> max keys CellComplex.cells

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

--Boundary function, which returns a hashtable with the cells in the boundary and their degrees, which is probably \pm 1
boundary := (cell) -> HashTable tally cell.boundary

--Check if something (a list) is a boundary
isCycle = method()
isCycle() := (lst) -> 

--Attach a cell
attach = method()
attach(CellComplex,List,Thing) := (baseComplex,boundary,label) -> (
    c := makeCell boundary;
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

--Get list of cells 
cells := (cellcomplex) -> toList cellcomplex.cells

----------------------------

-- THINGS TO DO:
-- - verify that boundary is a cycle 

----------------------------


beginDocumentation()
document { 
	Key => CellularResolutions,
	Headline => "A package for cellular resolutions of monomial ideals",
	EM "Cellular Resolutions", "TODO"
	}
       
end

restart
loadPackage("CellularResolutions", Reload => true)


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


