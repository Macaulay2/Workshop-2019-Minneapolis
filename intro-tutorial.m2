-- emacs (is great!) 
-- C = Ctrl, M = Alt (usually) 
---------------------------
-- C-x C-f : open/create a file in a buffer
-- F12     : start M2
-- C-x o   : jump between the buffers (using a mouse is OK too!)
-- F11     : execute the current line (or the highlighted block) 
-- C-x C-s : save the file (in the current buffer)
-- M-/     : autocomplete (not as frequently used as it should!)

-- RINGS
  -- fields
  KK = ZZ/911 -- pick your favorite prime  
  KK = QQ
  -- polynomial rings
  R1 = KK[x,y,z]
  R2 = KK[a,b][x,y,z]
  R3 = KK[a,b] ** KK[x,y,z] -- tensor product
  describe R1
  describe R2
  describe R3
  x -- is "global"
 
  -- more rings 
  R4 = KK[x,y,z]/(y-x^2,z-x^3)
  y*x
  R5 = frac R4 -- this one is a field
  f = z-y*x 
  use R1 -- but I wanted it in R1!!!
  f = z-y*x
  help Ring
  viewHelp Ring
  
-- IDEALS
  -- generators
  R = R1
  gens R
  I = ideal {y-x^2,z-x^3}
  gens I
  I_*
  I_0
  I_1
  -- Groebner bases
  gb I  
  gens gb I
  groebnerBasis I
  R = ZZ/32003 [x,y,z] 
  J = ideal {random(40,R), random(41,R), random(42,R)}  
  time G1 = gens gb J; -- ";" supresses output 
  time G2 = groebnerBasis(J, Strategy=>"F4");
  G1 == G2
  G1 === G2
  
-- POWERED by GB
  -- numerical invariants: dim, degree, hilbertSeries/Polynomial, betti, regularity
  dim I
  codim I
  degree I
  hilbertSeries I 
  betti I
  -- decomposition, etc.
  x = symbol x -- one way to "clear" x
  R = QQ[a,b]**QQ[c]**QQ[x_1..x_4]
  v = matrix {{a},{b}}
  A = genericMatrix(R,x_1,2,2)
  I = ideal(A*v - c*v) -- do you recognize this problem?
  decompose I  -- also check out the package "MinimalPrimes" 
  primaryDecomposition I
  ass I
  -- intersection, quotient, saturation
  I + ideal random({1,0,0},R)
  I : b     
  J = intersect(I,ideal{a,b^3}) 
  primaryDecomposition J
  J : b
  J : b^2
  saturate (J, b)  
  -- resolution
  C = res I
  C#dd  
  -- elimination
  restart -- one way to get a clean sheet
  R = QQ[x,y] 
  monoms = flatten entries basis(4,R)
  S = QQ[apply(monoms,m->z_m)]
  vars S
  RS = QQ[gens R, gens S]
  vars RS
  I = ideal apply(monoms,m->z_m-sub(m,RS)) -- the graph of the Veronese
  J = eliminate({x,y}, I) -- the ideal defining the image
  toString J
  toExternalString J
  tex J
  -- packages
  needsPackage "Dmodules"
  viewHelp Dmodules
  W = QQ[x,y,Dx,Dy, WeylAlgebra => {x=>Dx,y=>Dy}]
  I = ideal (x*Dx+2*y*Dy-3, Dx^2-Dy)
  isHolonomic I
  charIdeal I
  holonomicRank I
  bFunction (I,{1,0})

-- FLOATING POINT...
  -- RR and CC
  1.2
  R = CC[x,y,z]
  M = random(R^1,R^{-3,-3,-3}) -- 3 cubics
  I = ideal M
  gens gb I -- don't do that!!!
  degree I -- 27... or not?
  -- numerical linear algebra
  A = random(CC^2,CC^2)  
  B = random(CC^2,CC^1)
  solve(A,B) -- solve Ax = B for x
  eigenvalues A 
  eigenvectors A
  LUdecomposition A
  SVD A  

-- NONLINEAR numerical algebra (numerical AG)
  -- homotopy continuation
  needsPackage "NumericalAlgebraicGeometry"
  F1 = {x+y+z-1, x^2+y^2+z^2-1, x^3+y^3+z^3-ii}  
  (F0,sols0) = totalDegreeStartSystem F1
  track(F0,F1,sols0) -- track a homotopy (1-t)*F0+t*F1 for t in [0,1]
  -- blackbox solvers (of 0-dimensional systems) 
  solveSystem F1 
  solveSystem I_*
  #oo
  needsPackage "MonodromySolver"
  viewHelp MonodromySolver
  R = CC[x,y]
  F = {x^3*y^5 + y^2 + x^2*y, x*y + x^2 - 1}
  # solveSystem F
  sparseMonodromySolve polySystem F
  #oo
  needsPackage "PHCpack"
  mixedVolume F
 