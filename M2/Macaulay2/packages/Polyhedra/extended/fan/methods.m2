-- PURPOSE : Computing the stellar subdivision
--   INPUT : '(F,r)', where 'F' is a Fan and 'r' is a ray
--  OUTPUT : A fan, which is the stellar subdivision
stellarSubdivision = method()
stellarSubdivision (Fan,Matrix) := Fan => (F,r) -> (
   -- Checking for input errors
   if numColumns r != 1 or numRows r != ambDim F then error("The ray must be given by a one column matrix in the ambient dimension of the fan");
   divider := (C,r) -> (
      if dim C != 1 then (
         raysC := rays C;
         linC := linealitySpace C;
         flatten apply(faces(1,C), f -> (
            conef := coneFromVData(raysC_f, linC); 
            if not contains(conef,r) then coneFromVData {conef,r} else divider(conef,r)) 
         )
      )
      else {C}
   );
   L := flatten apply(values getProperty(F, honestMaxObjects), C -> if contains(C,r) then divider(C,r) else {C});
   fan L
)



-- PURPOSE : Computes the coarsest common refinement of a given set of rays
--   INPUT : 'M'  a Matrix
--  OUTPUT : 'F'  a Fan, the coarsest common refinement of the rays in 'M'
ccRefinement = method(TypicalValue => Fan)
ccRefinement Matrix := M -> (
     -- Checking for input errors
     M = chkZZQQ(M,"rays");
     -- Extracting data
     n := numRows M;
     m := numColumns M;
     -- Generating all cones generated by 'n' rays in 'M'
     nCones := apply(subsets(m,n), e -> coneFromVData M_e);
     -- Selecting those cones that are 'n' dimensional and do not contain any 
     -- of the others
     nConesfd := select(nCones, C -> dim C == n);
     nConesfd = inclMinCones nConesfd;
     refCones := {};
     while nConesfd != {} do (
	  newCones := {};
	  -- scan through the 'n' dimensional cones and check for each of the cones generated by
	  -- 'n' rays if their intersection is 'n' dimensional and if the first one is not contained 
	  -- in the latter. If true, then their intersection will be saved in the list 'newCones'.
	  -- If false for every cone generated by 'n' rays, then the 'n' dimensional cone will be 
	  -- appended to the list 'refCones'
	  refCones = refCones | (flatten apply(nConesfd, C1 -> (
			 toBeAdded := flatten apply(nCones, C2 -> (
				   C := intersection(C2,C1);
				   if dim C == n and (not contains(C2,C1)) then C
				   else {}));
			 if toBeAdded == {} then C1
			 else (
			      newCones = newCones | toBeAdded;
			      {}))));
	  -- now, the new intersections will be the 'n' dimensional cones and the same procedure 
	  -- starts over again if this list is not empty
	  nConesfd = unique newCones);
     -- Compute the fan generated by the 'refCones'
     fan refCones);



-- PURPOSE : Computes the minimal non-faces of a fan
--   INPUT : 'Phi' a Fan
--  OUTPUT : 'S' a list of minimal non-faces of the fan 'Phi'
minimalNonFaces = method();
minimalNonFaces Fan := Phi -> (
   F := faces Phi;
   S := set subsets toList(0..numColumns rays Phi -1); -- all possible subsets of rays
      for k in keys F do ( 
         S = S - set F#k; -- find the non-faces
         );
    -- for each subset, determine whether it contains another subset 
    S = select(toList S, 
         s -> position(toList S, t -> isSubset(t, s) and t != s) === null
         );
   S    
   )



-- PURPOSE : Computes the Stanley–Reisner ring of a fan
--   INPUT : 'Phi' a Fan
--  OUTPUT : 'R/I' a QuotientRing
stanleyReisnerRing = method();
stanleyReisnerRing Fan := Phi -> (
   S := minimalNonFaces Phi;
   x := getSymbol "x";
   R := QQ[x_0..x_(numColumns rays Phi -1)];
   I := ideal for s in S list product apply(s, i -> R_i);
   R / I
   )
