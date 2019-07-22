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

export {"CellComplex"}

CellComplex = new Type of HashTable

Cell = new Type of HashTable

cellComplex = method()
cellComplex(Ring) := (R) -> (
    labelRing = symbol "labelRing";
    h := new HashTable from {
	labelRing => R,
	cells => new MutableHashTable from {}};
    new CellComplex from h
    )

--Make cell 
makeCell := (lst) -> (
    bdim := -1;
    for cell in lst do ( 
	if bdim < 0 
	then bdim = dim cell
	else assert(bdim == dim cell)
	);
    n := bdim + 1;
    h := new HashTable from {
	cellDimension => n, --FIGURE IT OUT
	boundary => lst -- could verify that it's a list
    	};
    new Cell from h
    );

--Attach a cell
attach = method()
attach(CellComplex,List,RingElement) := (baseComplex,boundary,label) -> (
    bdim := -1;
    for cell in boundary do ( 
	if bdim < 0 
	then bdim = dim cell
	else assert(bdim == dim cell)
	);
    n := bdim + 1;
    c := makeCell boundary;
    if baseComplex.cells#?n 
    then(
        i = #baseComplex.cells#n;
        baseComplex.cells#n#i=c;
	)
    else (
	baseComplex.cell#n = new MutableList from {c};
	);
    c
    )
attach(CellComplex,List) := (baseComplex,cells) -> (
    )
attach(CellComplex,List,Number) := (baseComplex,cells) -> (
    )




beginDocumentation()
document { 
	Key => CellularResolutions,
	Headline => "A package for cellular resolutions of monomial ideals",
	EM "Cellular Resolutions", "TODO"
	}
       
end

restart
loadPackage("CellularResolutions")


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


