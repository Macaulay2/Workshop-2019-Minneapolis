restart
load "nops.m2"

R = QQ[x,y,z,w]
J = monomialCurveIdeal(R, {1,2,3})
J = J^2
(f,I,t) = noetherNormalization J

nx = 2;
nd = 3;
var = {x,y};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))



f' =inverse mapToMatrix f
-- Points
p = matrix{{0,0,0,1}}; -- vanishes on the original ideal
sub(J,p)

p' = p * transpose f' -- vanishes on the new ideal
sub(I,p')

M' = sub(M, p')
K = gens kernel M'
transpose (bd * K)


--------------------------------------------------------


p = matrix{{1,0,0,0}}
sub(J,p) -- sanity check, should get zero ideal

p' = p * transpose f'
sub(I,p') -- sanity check, should get zero again

M' = sub(M,p')
K = gens kernel M'
transpose (bd * K)


--------------------------------------------------------


p = matrix{{1,1,1,1}}
sub(J,p) -- sanity check, should get zero ideal

p' = p * transpose f'
sub(I,p') -- sanity check, should get zero again

M' = sub(M,p')
K = gens kernel M'
transpose (bd * K)


--------------------------------------------------------


p = matrix{{4, -1, 1/4, -1/16}}
sub(J,p) -- sanity check, should get zero ideal

p' = p * transpose f'
sub(I,p') -- sanity check, should get zero again

M' = sub(M,p')
K = gens kernel M'
transpose (bd * K)


--------------------------------------------------------
--------------------------------------------------------
---------------------Friday demo------------------------
--------------------------------------------------------
--------------------------------------------------------

-- Prototype
R = QQ[x,y,z]
line1 = ideal(x^2 - z*y, y^2)
line2 = ideal(x+y+z,x-y+z)

I = intersect(line1,line2)
pd = primaryDecomposition I

MacaulayMatrixPD({x,y}, 10,5,pd#1)

nx = 5;
nd = 5;
var = {x,y};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))


-- line1 nops
p = matrix{{0,0,1}}
M' = sub (M, p);
K = gens kernel M';
transpose (bd * K)

p = matrix{{0,0,2}}
M' = sub (M, p);
K = gens kernel M';
transpose (bd * K)

p = matrix{{0,0,3}}
M' = sub (M, p);
K = gens kernel M';
transpose (bd * K)

p = matrix{{0,0,4}}
M' = sub (M, p);
K = gens kernel M';
transpose (bd * K)




-- Test with embedded components
-- Fat line + fat point in 2d
R = QQ[x,y]
line = ideal((x^2-y)^2)
MacaulayMatrixPD({x},10,3,line) -- 1, dx

point = ideal((x+1),(y-1))
MacaulayMatrix(10,5,point) -- 1

I = line*point

pd = primaryDecomposition I

MacaulayMatrixPD({x}, 10,5, pd#0) -- 1, dx
MacaulayMatrix(10,5,pd#1) -- dx^2, dx, 1

-- Evaluate
nx = 10;
nd = 10;
var = {x};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))


-- Fat point
p = matrix{{-1,1}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)


-- Points on the line
p = matrix{{0,0}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{-2,4}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{3,9}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

--------------------------------------------------------
--------------------------------------------------------

-- Point embedded in line embedden in surface
R = QQ[x,y,z]
surf = ideal(x^2 + y^2 - 1)
line = surf + ideal(x^2 + z^2 - 1)
point = ideal(x-1,y,z)

I = surf * line * point

pd = primaryDecomposition(I)

-- Cylinder x^2 + y^2 - 1
-- Prime ideal
pd#0
dim pd#0
degree pd#0
MacaulayMatrixPD({x}, 10, 5, pd#0)


-- Line1, ideal((x2 + z2 - 1)^2, y-z)
pd#1
dim pd#1
degree pd#1
MacaulayMatrixPD({x,y}, 10, 5, pd#1) -- dx, 1

-- Line2, ideal((x2 + z2 - 1)^2, y+z)
pd#2
dim pd#2
degree pd#2
MacaulayMatrixPD({x,y}, 10, 5, pd#2) -- dx, 1

-- Point, (1,0,0)
pd#3
dim pd#3
degree pd#3
MacaulayMatrix(10, 5, pd#3) -- 1, dz, dy, dx, dx*dz, dy^2, dx*dy, dx^2, dx*dy



-- Evaluate
nx = 6;
nd = 6;
var = {x,y,z};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))



-- cylinder
p = matrix{{1,0,1}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{3/5,-4/5,-1}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

-- line1
p = matrix{{1,0,0}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{-1,0,0}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,1,0}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{3/5,4/5,4/5}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

-- line2
p = matrix{{3/5,4/5,-4/5}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)


-- point
p = matrix{{1,0,0}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)


--------------------------------------------------------
--------------------------------------------------------


R = QQ[x,y,z]
surface = ideal(x^2, x*y)
line = ideal(x^2 - z*y, y^2)

I = surface*line
pd = primaryDecomposition I

MacaulayMatrixPD({x,y},10,10,pd#1)

nx = 5;
nd = 5;
var = {x,y};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))


-- Surface nops
p = matrix{{0,2,-2}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

-- Line nops
p = matrix{{0,0,1}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,0,2}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,0,3}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,0,4}}
M' = sub (M, p)
K = gens kernel M'
transpose (bd * K)









--------------------------------------------------------
--------------------------------------------------------
---------------------End Friday demo--------------------
--------------------------------------------------------
--------------------------------------------------------






-- Test with embedded components
restart
load "nops.m2"
R = QQ[x,y,z]
I1 = ideal"x3y+y3z+z3x" -- surface
I2 = I1 + ideal"x2-y" -- two lines
J = I1 * I2 -- surface with two embedded lines

(f, I, t) = noetherNormalization(J)

nx = 4;
nd = 5;
var = {x};

fi = gens I
bx = basis(0,nx,R, Variables => var)
bd = basis(0,nd,R, Variables => var)

-- macaulay matrix
M = transpose diff(transpose bd, flatten (transpose fi*bx))



f' =inverse mapToMatrix f
-- Points on hypersurface 
p = matrix{{0,0,0}}; -- vanishes on the original ideal
sub(J,p)
p' = p * transpose f' -- vanishes on the new ideal
sub(I,p')
M' = sub(M, p')
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,0,1}}; -- vanishes on one of the lines
sub(J,p)
p' = p * transpose f' -- vanishes on the new ideal
sub(I,p')
M' = sub(M, p')
K = gens kernel M'
transpose (bd * K)

p = matrix{{0,1,0}}; -- vanishes on the original ideal
sub(J,p)
p' = p * transpose f' -- vanishes on the new ideal
sub(I,p')
M' = sub(M, p')
K = gens kernel M'
transpose (bd * K)

p = matrix{{1,0,0}}; -- vanishes on the original ideal
sub(J,p)
p' = p * transpose f' -- vanishes on the new ideal
sub(I,p')
M' = sub(M, p')
K = gens kernel M'
transpose (bd * K)


