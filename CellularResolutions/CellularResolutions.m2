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

CellComplex = new Type of MutableHashTable

Cell = new Type of HashTable

cellComplex = method()
cellComplex(Ring) := (R) -> (
    labelRing = symbol "labelRing";
    h := new HashTable from { labelRing => R
                            , cells => new MutableHashTable from {}};
    new CellComplex from h
    )

--Attach a cell
attach = method()
attach(CellComplex,List,RingElement) := (baseComplex,c) -> (
    if baseComplex.cells#?n then(
        i = #baseComplex.cells;
        baseComplex.cells#n#i=c;
        )
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


