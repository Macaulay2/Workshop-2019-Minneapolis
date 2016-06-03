-- Test cayley_polytope/c2.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1,1},{0,1,0,0,1,0,1,0,1},{0,0,1,0,0,1,1,0,0},{1,1,1,0,0,0,0,0,0},{0,0,0,1,1,1,1,0,0}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1},{0,0,0,0,0,1,1,1,1},{0,0,0,1,1,0,0,0,1},{0,0,1,0,1,0,0,1,0},{0,1,0,1,0,0,1,0,1}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test gc_closure/3-0.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1},{3/2,-1/2,1/2,1/2,1/2,1/2,1/2,1/2},{1/2,1/2,3/2,-1/2,1/2,1/2,1/2,1/2},{1/2,1/2,1/2,1/2,3/2,-1/2,1/2,1/2},{1/2,1/2,1/2,1/2,1/2,1/2,3/2,-1/2}};
desiredHB = matrix {{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3},{-1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3},{1,0,1,1,1,1,1,2,-1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,0,1,1,1,1,1,2,1,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2},{1,1,0,1,1,1,2,1,1,0,1,1,1,2,-1,0,0,0,1,1,1,1,1,2,2,2,3,0,1,1,1,2,1,1,0,1,1,1,2,1,1,1,1,2,2,1,1,2,2,1,1,2,2,0,0,1,1,1,1,2,2,2,2,3,3,0,0,1,1,1,1,2,2,2,2,3,3,1,1,2,2,1,1,2,2,0,0,1,1,1,1,2,2,2,2,3,3,0,0,1,1,1,1,2,2,2,2,3,3,1,1,2,2,1,1,2,2,1,1,2,2},{1,1,1,0,1,2,1,1,1,1,0,1,2,1,1,0,1,2,-1,0,1,2,3,0,1,2,1,1,0,1,2,1,1,1,1,0,1,2,1,1,1,1,2,1,2,1,2,1,2,1,2,1,2,1,2,0,1,2,3,0,1,2,3,1,2,1,2,0,1,2,3,0,1,2,3,1,2,1,2,1,2,1,2,1,2,1,2,0,1,2,3,0,1,2,3,1,2,1,2,0,1,2,3,0,1,2,3,1,2,1,2,1,2,1,2,1,2,1,2,1,2}};
desiredHB = sort promote(desiredHB, QQ);
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test SCHLEGEL_VERTEX_COLORS/1.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1},{-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1},{-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test VISUAL_DUAL_FACE_LATTICE/1.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1},{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1},{-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1},{-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test zonotope/6.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1},{0,0,2,0,2,2,0,0,2,2,2,2,2,0,0,0},{0,4,4,0,0,4,4,0,0,4,0,4,0,4,0,4},{0,0,8,8,0,8,0,8,0,0,8,0,8,8,0,8}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},{0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4},{0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test zonotope/9.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1,1,1,1,1},{0,2,0,2,0,2,0,2,0,2,0,2},{0,0,2,2,4,4,0,0,2,2,4,4},{0,0,2,2,0,0,-2,-2,-4,-4,-2,-2}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},{0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4},{-2,-1,0,-3,-2,-1,0,1,-4,-3,-2,-1,0,1,2,-3,-2,-1,0,1,-2,-1,0,-2,-1,0,-3,-2,-1,0,1,-4,-3,-2,-1,0,1,2,-3,-2,-1,0,1,-2,-1,0,-2,-1,0,-3,-2,-1,0,1,-4,-3,-2,-1,0,1,2,-3,-2,-1,0,1,-2,-1,0}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test random_edge_epl/1.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1},{1,0,0,1,0,1,0,1},{3/4,1,0,1/4,0,1/4,1,3/4},{3/16,1/4,0,1/16,1,15/16,3/4,13/16}};
desiredHB = matrix {{1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,13,13,14,14,15,15,16,16,16,16},{0,0,0,0,1,2,0,0,1,1,2,2,3,3,3,3,0,0,1,1,2,2,3,3,3,4,4,4,4,4,4,4,4,1,1,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,13,13,14,14,15,15,16,16,16,16},{0,0,1,2,1,1,3,3,2,2,2,2,1,1,2,2,4,4,3,3,3,3,3,3,3,1,1,1,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,2,2,3,3,4,4,2,2,3,3,4,4,2,2,3,3,4,4,2,2,3,3,4,4,3,3,4,4,3,3,4,4,3,3,4,4,8,8,3,3,4,4,4,4,4,4,4,4,4,4,12,12},{0,1,1,1,1,1,1,2,1,2,1,2,1,2,1,2,1,3,1,3,1,3,1,2,3,1,2,3,1,3,1,2,3,1,4,1,4,1,4,1,2,3,4,1,4,1,4,1,5,1,5,1,5,1,5,1,6,1,6,1,6,1,7,1,7,1,7,1,8,1,8,1,9,1,9,1,10,1,10,2,9,1,11,1,11,1,12,1,13,1,14,1,15,3,13}};
desiredHB = sort promote(desiredHB, QQ);
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

-- Test schlegel_params/2.poly
-- Checking hilbertBasis
TEST ///
raysC = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
desiredHB = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1},{-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1},{-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1}};
desiredHB = sort desiredHB;
C = posHull(raysC)
computedHB = sort matrix {hilbertBasis C};
assert(desiredHB == computedHB);
///

