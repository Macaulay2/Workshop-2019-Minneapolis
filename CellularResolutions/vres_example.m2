needsPackage "CellResolutions"
needsPackage "NormalToricVarieties"



mkShortVResMonomial := (tv,I,d) -> (
    R = ring tv;
    F = faces polytope fan tv;
    C = cellComplex(ring tv);
    n := dim tv
    ngens := #rays tv
    coneAsIdeal := c -> (
        cSet = set c;
        ideal for i from 0 to ngens-1 list if cSet#?i then R_i else continue 
        )
    
    irrComponents = applyPairs(orbits tv,
        (cd,cones) -> cd => apply(cones,coneAsIdeal))
    --add 0 cells
    cellsTable = new MutableHashTable
    for c from 0 to #(irrComponents#0)-1 do (
        cellsTable#{c} = attach(C,{},intersect((irrComponent#0#c)^d,I));
        );
    for i from 1 to n do (
        for face in F#(n-i) do (
            for s  bv subsets(face,(#face))
            attach(C,{...},);
            )
        )
    
    )

restart
needsPackage "CellularResolutions"
needsPackage "NormalToricVarieties"
tv = toricProjectiveSpace(2)**toricProjectiveSpace(1);
R = ring tv

v013 = newCell({},R_2*R_4)
v014 = newCell({},R_2*R_3)
v023 = newCell({},R_1*R_4)
v024 = newCell({},R_1*R_3)
v123 = newCell({},R_0*R_4)
v124 = newCell({},R_0*R_3)
l01 = newCell({v013,v014});
l02 = newCell({v023,v024});
l03 = newCell({v013,v023});
l04 = newCell({v014,v024});
l12 = newCell({v123,v124});
l13 = newCell({v013,v123});
l14 = newCell({v014,v124});
l23 = newCell({v023,v123});
l24 = newCell({v024,v124});
f0=newCell({l01,l02,l03,l04});
f1=newCell({l01,l12,l13,l14});
f2=newCell({l02,l12,l23,l24});
f3=newCell({l03,l13,l23});
f4=newCell({l04,l14,l24});
p = newCell({f0,f1,f2,f3,f4});
C = cellComplex(R,{p});
CchainComplex = chainComplex C;
prune (CchainComplex[-1])
prune HH prune (CchainComplex[-1])
