-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1},{0,0,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{0},{0},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 2, vertices: 1, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3},{4}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0},{0,-1}};
eqrhsPd = matrix {{-3},{-4}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 11, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1,0,1,0},{0,0,1,1,0,0,1,1,0,0,1},{1,1,1,1,0,0,0,0,0,0,0},{0,0,0,0,1,1,1,1,0,0,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0},{-1,0,0,0},{0,-1,0,0},{1,1,-1,-1},{0,0,0,-1},{0,1,0,0},{1,0,0,0},{0,0,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 5, vertices: 11, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1,0,1,0},{0,0,1,1,0,0,1,1,0,0,1},{1,1,1,1,0,0,0,0,0,0,0},{0,0,0,0,1,1,1,1,0,0,0},{0,0,0,0,0,0,0,0,1,1,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0,0},{-1,0,0,0,0},{0,-1,0,0,0},{1,1,-1,-1,0},{0,0,0,-1,0},{0,1,0,0,0},{1,0,0,0,0},{0,0,1,1,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{0},{1},{1},{1}};
eqlhsPd = matrix {{0,0,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 7, facets: 7
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,1,1,3},{0,0,1,1,0,0,4},{1,1,1,1,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,1,0}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,-3,-2,-2},{0,-1,-4,-4,-4},{1,0,2,2,2},{0,1,3,4,4},{0,0,0,0,-1},{0,0,0,-1,0},{0,0,1,1,1}};
ineqrhsPd = matrix {{-3},{-4},{3},{4},{0},{0},{1}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0},{0,0,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 5, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,1,3},{0,0,1,1,0,4},{1,1,1,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,-3,-2,0},{0,-1,-4,-4,0},{1,0,2,2,0},{0,1,3,4,0},{0,0,0,-1,0},{0,0,1,1,0}};
ineqrhsPd = matrix {{-3},{-4},{3},{4},{0},{1}};
eqlhsPd = matrix {{0,0,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 2, vertices: 1, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0},{0,-1}};
eqrhsPd = matrix {{-1},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 4, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,0,0,0},{1,0,0,1,1,0},{0,1,0,1,0,1},{0,0,1,0,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0,0},{0,0,1,0},{-1,0,0,0},{-1,-1,-1,0},{0,1,0,0},{0,0,-1,0},{1,1,1,0},{1,0,0,0}};
ineqrhsPd = matrix {{0},{1},{0},{-1},{1},{0},{2},{1}};
eqlhsPd = matrix {{-1,-1,-1,-1}};
eqrhsPd = matrix {{-2}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 10, facets: 30
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,-1,0,0,0,1,-1},{0,1,0,0,0,-1,0,0,1,-1},{0,0,1,0,0,0,-1,0,1,-1},{0,0,0,1,0,0,0,-1,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,-1,-1},{1,0,-1,-1},{1,-1,1,-1},{1,-1,0,-1},{1,-1,-1,1},{1,-1,-1,0},{-1,-1,1,1},{-1,-1,1,0},{0,-1,1,-1},{-1,0,1,-1},{-1,0,-1,1},{0,-1,-1,1},{-1,-1,0,1},{0,1,-1,-1},{-1,1,-1,0},{-1,1,-1,1},{-1,1,0,-1},{-1,1,1,-1},{-1,0,1,1},{-1,1,0,1},{-1,1,1,0},{0,1,1,-1},{0,1,-1,1},{0,-1,1,1},{1,-1,0,1},{1,-1,1,0},{1,0,-1,1},{1,1,-1,0},{1,0,1,-1},{1,1,0,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 11, facets: 26
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,0,-1,0,0,0,0,1},{0,1,0,0,0,0,-1,0,0,0,1},{0,0,1,0,0,0,0,-1,0,0,1},{0,0,0,1,0,0,0,0,-1,0,1},{0,0,0,0,1,0,0,0,0,-1,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,-1,-1,-1},{1,-1,1,-1,-1},{1,-1,-1,1,-1},{1,-1,-1,-1,-1},{1,-1,-1,-1,1},{-1,-1,1,1,-1},{-1,-1,1,-1,-1},{-1,-1,1,-1,1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,1},{-1,-1,-1,1,1},{-1,-1,-1,1,-1},{-1,1,-1,-1,1},{-1,1,-1,-1,-1},{-1,1,-1,1,-1},{-1,1,1,-1,-1},{-1,1,1,1,-1},{-1,1,1,-1,1},{-1,1,-1,1,1},{-1,-1,1,1,1},{1,-1,-1,1,1},{1,-1,1,-1,1},{1,-1,1,1,-1},{1,1,-1,-1,1},{1,1,-1,1,-1},{1,1,1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 6, vertices: 13, facets: 102
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,0,0,-1,0,0,0,0,0,1},{0,1,0,0,0,0,0,-1,0,0,0,0,1},{0,0,1,0,0,0,0,0,-1,0,0,0,1},{0,0,0,1,0,0,0,0,0,-1,0,0,1},{0,0,0,0,1,0,0,0,0,0,-1,0,1},{0,0,0,0,0,1,0,0,0,0,0,-1,1}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,-1,-1,-1},{1,1,-1,1,-1,-1},{1,1,-1,-1,1,-1},{1,1,-1,-1,-1,-1},{1,1,-1,-1,-1,1},{1,-1,1,1,-1,-1},{1,-1,1,-1,1,-1},{1,-1,1,-1,-1,-1},{1,-1,1,-1,-1,1},{1,-1,-1,1,1,-1},{1,-1,-1,1,-1,-1},{1,-1,-1,1,-1,1},{1,-1,-1,-1,1,-1},{1,-1,-1,-1,1,1},{1,-1,-1,-1,-1,-1},{1,-1,-1,-1,-1,1},{-1,-1,1,1,1,-1},{-1,-1,1,1,-1,-1},{-1,-1,1,1,-1,1},{-1,-1,1,-1,1,-1},{-1,-1,1,-1,1,1},{-1,-1,1,-1,-1,-1},{-1,-1,1,-1,-1,1},{-1,-1,-1,-1,1,-1},{-1,-1,-1,-1,1,1},{-1,-1,-1,-1,-1,1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,1,-1,1},{-1,-1,-1,1,-1,-1},{-1,-1,-1,1,1,1},{-1,-1,-1,1,1,-1},{-1,1,-1,-1,-1,1},{-1,1,-1,-1,-1,-1},{-1,1,-1,-1,1,1},{-1,1,-1,-1,1,-1},{-1,1,-1,1,-1,1},{-1,1,-1,1,-1,-1},{-1,1,-1,1,1,-1},{-1,1,1,-1,-1,1},{-1,1,1,-1,-1,-1},{-1,1,1,-1,1,-1},{-1,1,1,1,-1,-1},{-1,-1,0,1,1,1},{-1,-1,1,0,1,1},{-1,-1,1,1,0,1},{-1,-1,1,1,1,0},{-1,0,1,1,1,-1},{-1,0,1,1,-1,1},{-1,0,1,-1,1,1},{-1,0,-1,1,1,1},{-1,1,-1,0,1,1},{-1,1,-1,1,0,1},{-1,1,-1,1,1,0},{-1,1,0,-1,1,1},{-1,1,1,-1,0,1},{-1,1,1,-1,1,0},{-1,1,0,1,-1,1},{-1,1,1,0,-1,1},{-1,1,1,1,-1,0},{-1,1,0,1,1,-1},{-1,1,1,0,1,-1},{-1,1,1,1,0,-1},{0,1,1,1,-1,-1},{0,1,1,-1,1,-1},{0,1,1,-1,-1,1},{0,1,-1,1,1,-1},{0,1,-1,1,-1,1},{0,1,-1,-1,1,1},{0,-1,1,1,1,-1},{0,-1,1,1,-1,1},{0,-1,1,-1,1,1},{0,-1,-1,1,1,1},{1,-1,-1,0,1,1},{1,-1,-1,1,0,1},{1,-1,-1,1,1,0},{1,-1,0,-1,1,1},{1,-1,1,-1,0,1},{1,-1,1,-1,1,0},{1,-1,0,1,-1,1},{1,-1,1,0,-1,1},{1,-1,1,1,-1,0},{1,-1,0,1,1,-1},{1,-1,1,0,1,-1},{1,-1,1,1,0,-1},{1,0,-1,-1,1,1},{1,1,-1,-1,0,1},{1,1,-1,-1,1,0},{1,0,-1,1,-1,1},{1,1,-1,0,-1,1},{1,1,-1,1,-1,0},{1,0,-1,1,1,-1},{1,1,-1,0,1,-1},{1,1,-1,1,0,-1},{1,0,1,-1,-1,1},{1,1,0,-1,-1,1},{1,1,1,-1,-1,0},{1,0,1,-1,1,-1},{1,1,0,-1,1,-1},{1,1,1,-1,0,-1},{1,0,1,1,-1,-1},{1,1,0,1,-1,-1},{1,1,1,0,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^6, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,-1,0,0,1,-1},{0,1,0,0,-1,0,1,-1},{0,0,1,0,0,-1,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,1},{1,1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 7, vertices: 19, facets: 11
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1},{0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1},{0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0},{0,1,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0},{1,1,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0,1},{1,1,1,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0},{1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,0,0,1,0},{-1,0,0,0,0,0,0},{-1,-1,0,0,-1,-1,0},{0,-1,0,0,0,0,0},{0,0,0,0,1,0,0},{0,1,0,0,0,0,0},{0,0,0,0,-1,-1,0},{1,0,0,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,-1,0,0,0},{1,1,1,1,1,1,0}};
ineqrhsPd = matrix {{1},{0},{-2},{0},{1},{1},{-1},{1},{0},{0},{3}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 7, vertices: 35, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1},{0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1},{0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0},{1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,0,0,1,0,0},{1,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,1,0,0,0},{1,1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,0,0,1,0},{0,-1,0,0,0,0,0},{-1,0,0,0,0,0,0},{-1,-1,-1,-1,-1,-1,0},{0,0,0,0,1,0,0},{0,0,0,-1,0,0,0},{0,0,1,0,0,0,0},{0,0,0,0,-1,0,0},{0,0,0,0,0,-1,0},{1,1,1,1,1,1,0},{0,1,0,0,0,0,0},{1,0,0,0,0,0,0}};
ineqrhsPd = matrix {{0},{1},{1},{0},{0},{-2},{1},{0},{1},{0},{0},{3},{1},{1}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 7, vertices: 30, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1},{0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0},{0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0},{1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,1,0,1},{1,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0},{1,1,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,0,0}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,0,0,1,0},{0,0,0,0,-1,0,0},{0,0,0,1,0,0,0},{-1,0,0,0,0,0,0},{-1,-1,-1,-1,-1,-1,0},{0,0,1,0,0,0,0},{0,-1,0,0,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,0,-1,0},{0,1,0,0,0,0,0},{0,0,-1,-1,-1,-1,0},{1,0,0,0,0,0,0},{0,0,-1,0,0,0,0},{1,1,1,1,1,1,0},{0,0,0,-1,0,0,0},{1,1,1,1,0,0,0}};
ineqrhsPd = matrix {{1},{0},{1},{0},{-2},{1},{0},{1},{0},{1},{-1},{1},{0},{3},{0},{2}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,1/2},{0,0,11/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{11,1},{0,-1},{-11,1}};
ineqrhsPd = matrix {{11},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3/2,-1/2,1/2,1/2,1/2,1/2,1/2,1/2},{1/2,1/2,3/2,-1/2,1/2,1/2,1/2,1/2},{1/2,1/2,1/2,1/2,3/2,-1/2,1/2,1/2},{1/2,1/2,1/2,1/2,1/2,1/2,3/2,-1/2}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,1,1,-1},{-1,1,1,1},{-1,1,-1,-1},{-1,1,-1,1},{-1,-1,-1,-1},{-1,-1,-1,1},{-1,-1,1,1},{-1,-1,1,-1},{1,-1,-1,1},{1,-1,-1,-1},{1,-1,1,1},{1,-1,1,-1},{1,1,-1,1},{1,1,-1,-1},{1,1,1,1},{1,1,1,-1}};
ineqrhsPd = matrix {{1},{2},{0},{1},{-1},{0},{1},{0},{1},{0},{2},{1},{2},{1},{3},{2}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1,0,1,0,1},{1,1,0,0,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{0,-1,0},{-1,0,0},{1,0,0},{0,1,0},{0,0,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1},{0,0,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,0},{1,0},{0,1}};
ineqrhsPd = matrix {{0},{0},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,1/2},{0,0,5}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-10,1},{10,1}};
ineqrhsPd = matrix {{0},{0},{10}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 24, facets: 32
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,2/3,2/3,1/3,1/2,1/2,2/3,2/3,1/3,1/3,1/3,2/3,1/3,2/3,1/2,1/3,1/3,0,2/3,1/3,2/3,1/2,1/2,1/2},{1/2,1/3,2/3,2/3,0,1/2,1/3,2/3,1/3,2/3,1/3,1/3,1/3,1/3,1/2,2/3,1/3,1/2,2/3,2/3,2/3,1,1/2,1/2},{1/2,2/3,2/3,1/3,1/2,0,1/3,1/3,1/3,1/3,1/3,1/3,2/3,2/3,1/2,2/3,2/3,1/2,1/3,2/3,2/3,1/2,1,1/2},{1/2,2/3,1/3,2/3,1/2,1/2,2/3,1/3,2/3,1/3,1/3,1/3,1/3,1/3,0,1/3,2/3,1/2,2/3,2/3,2/3,1/2,1/2,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,-1,-1},{-1,0,-1,-1},{-1,-1,0,-1},{-1,-1,-1,0},{1,0,-1,-1},{1,-1,0,-1},{1,-1,-1,0},{0,1,-1,-1},{0,-1,1,-1},{0,-1,-1,1},{-1,1,0,-1},{-1,1,-1,0},{-1,0,1,-1},{-1,0,-1,1},{-1,-1,1,0},{-1,-1,0,1},{1,1,0,-1},{1,1,-1,0},{1,0,1,-1},{1,0,-1,1},{1,-1,1,0},{1,-1,0,1},{0,1,1,-1},{0,1,-1,1},{0,-1,1,1},{-1,1,1,0},{-1,1,0,1},{-1,0,1,1},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1}};
ineqrhsPd = matrix {{-1},{-1},{-1},{-1},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{2},{2},{2},{2}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3/2,-1/2,1/2,1/2},{1/2,1/2,3/2,-1/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1},{-1,1},{1,1},{1,-1}};
ineqrhsPd = matrix {{0},{1},{2},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{0},{1},{0},{1},{0},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 2, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,-1},{1/2,-1/2}};
eqrhsPd = matrix {{0},{1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 0, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^0, QQ^1, 0);
raysP = map(QQ^0, QQ^0, 0);
linealityP = map(QQ^0, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^1, QQ^0, 0);
ineqrhsPd = matrix {{1}};
eqlhsPd = map(QQ^0, QQ^0, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 2, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = matrix {{0},{1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1,3/2,5/2,5/2},{0,1/2,5/2,5/2,0,1/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,0},{-4,2},{0,1},{2,1},{1,0}};
ineqrhsPd = matrix {{0},{0},{1},{5/2},{11/2},{5/2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{5/3,-1/3,-1,1,1/3},{1/3,-5/3,1,-1,5/3}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-2},{-4,-1},{-1,2},{1,1},{2,-1}};
ineqrhsPd = matrix {{3},{3},{3},{2},{3}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 2, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0},{0,-1}};
eqrhsPd = matrix {{-1},{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1},{-1},{-1}};
raysP = matrix {{1,0,0},{0,1,0},{0,0,1}};
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,15}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1},{-1}};
ineqrhsPd = matrix {{15},{0}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,0}};
raysP = matrix {{0},{1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1}};
ineqrhsPd = matrix {{0},{1},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: -1, ambientDim: 1, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^1, QQ^0, 0);
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^1, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: -1, ambientDim: 2, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^2, QQ^0, 0);
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^2, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: -1, ambientDim: 1, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^1, QQ^0, 0);
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^1, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 2, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0},{0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = matrix {{1},{1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{1,-1}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: -1, ambientDim: 2, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^2, QQ^0, 0);
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^2, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 8, ambientDim: 9, vertices: 12, facets: 36
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,2,2,3,3,0,0,0,0,0,0},{2,3,1,3,1,2,0,0,0,0,0,0},{3,2,3,1,2,1,0,0,0,0,0,0},{1,0,0,0,0,0,1,0,0,0,0,0},{0,1,0,0,0,0,0,1,0,0,0,0},{0,0,1,0,0,0,0,0,1,0,0,0},{0,0,0,1,0,0,0,0,0,1,0,0},{0,0,0,0,1,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0,0,0,0,1}};
raysP = map(QQ^9, QQ^0, 0);
linealityP = map(QQ^9, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,1,1,1,1,1,0},{-5,1,7,-18,-12,-12,0,0,0},{1,-5,1,0,0,0,0,0,0},{-1,-7,5,0,0,-6,0,0,0},{-1,-1,2,-3,0,-3,0,0,0},{-7,-1,5,-6,0,0,0,0,0},{-5,1,1,0,0,0,0,0,0},{-1,-1,1,0,0,0,0,0,0},{-1,1,-1,0,0,0,0,0,0},{-7,5,-1,0,-6,0,0,0,0},{-2,1,1,-3,-3,0,0,0,0},{-5,7,1,-12,-18,0,-12,0,0},{-1,2,-1,0,-3,0,-3,0,0},{-1,5,-7,0,0,0,-6,0,0},{1,1,-5,0,0,0,0,0,0},{1,-1,-1,0,0,0,0,0,0},{5,-1,-7,6,6,6,6,6,0},{5,-7,-1,0,0,0,0,-6,0},{2,-1,-1,3,3,3,3,0,0},{1,1,-2,3,3,3,0,3,0},{7,1,-5,18,18,18,6,6,0},{1,-2,1,0,0,-3,0,-3,0},{7,-5,1,12,12,0,12,-6,0},{1,-5,7,-12,0,-18,0,-12,0},{-1,-1,5,-12,-6,-12,0,-6,0},{1,-1,1,0,2,-2,2,-2,0},{1,7,-5,12,0,12,-6,12,0},{-1,1,1,-4,-4,-2,-2,0,0},{-1,5,-1,0,-6,6,-6,6,0},{1,1,-1,4,2,4,0,2,0},{5,-1,-1,12,12,6,6,0,0},{0,0,0,0,0,0,0,-1,0},{0,0,0,0,0,0,-1,0,0},{0,0,0,-1,0,0,0,0,0},{0,0,0,0,-1,0,0,0,0},{0,0,0,0,0,-1,0,0,0}};
ineqrhsPd = matrix {{1},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{6},{0},{3},{3},{18},{0},{12},{0},{0},{2},{12},{0},{6},{4},{12},{0},{0},{0},{0},{0}};
eqlhsPd = matrix {{0,0,0,-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 10, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3330883453830711/4503599627370496,853624150927223/2251799813685248,-2356553346050111/18014398509481984,8539723551892617/9007199254740992,-851543628576541/1125899906842624,5316837853882195/9007199254740992,-6538305985494021/9007199254740992,-1749267935005635/4503599627370496,7424579411304029/9007199254740992,-6887392800778413/18014398509481984},{-2298030817950255/18014398509481984,7698990389472477/9007199254740992,5054138723166979/9007199254740992,-4611252879029533/72057594037927936,7431314160210763/36028797018963968,3306083656494553/4503599627370496,-5119504026856145/9007199254740992,3819013085252451/4503599627370496,4143898238679627/9007199254740992,-5765394822723937/9007199254740992},{-5952332624486941/9007199254740992,-1596589986555531/4503599627370496,-7361860108925987/9007199254740992,-1402721951416725/4503599627370496,-1397989901669917/2251799813685248,3023262654016433/9007199254740992,-6977355308041419/18014398509481984,3248144300909891/9007199254740992,-1486048987021867/4503599627370496,-6002567455236229/9007199254740992}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{51894724291563534050210877879505/44579928701577317537651253200841,1,-184914068801775237924440155455677/89159857403154635075302506401682},{1,15701373045711447922845966302445/15404771123625084781539131855594,-26933806422179042513723576438601/15404771123625084781539131855594},{21859749371173593885550183008062/2404747912716946201299872316723,24507876296948960756894875847431/2404747912716946201299872316723,-1},{1,55228066624563312478254429271949/6621599294598900589652414097938,7629217969444610664921388076327/6621599294598900589652414097938},{-1,-120142383507582197243199728437381/8084188145793725279203688369944,115580771581373778595914316892651/4042094072896862639601844184972},{-307888710630698160800930201114164/175792615486297400307821108914777,445367236960222658287647234346636/175792615486297400307821108914777,-1},{-29609553812244423504592728106921/7297571300091599151751392189772,-1,-214868202830604470995803376144337/14595142600183198303502784379544},{-151189885071826254692495399928433/57582004016524483838099502969478,-1,-171320946886386780877737144372595/57582004016524483838099502969478},{-53700167453247199213158679994849/3307137770655639630104724597440,1,98822397682530192523626532018309/18189257738606017965575985285920},{-1,116803526287883759434174588839635/34580975135927942243574654312938,-36007901243263107943040099596727/34580975135927942243574654312938},{1,-13850085013344855522773355887391/6818528060217462587852947510909,-197756724530638510637354434024943/13637056120434925175705895021818},{1036135618271528606489034912963057307981344407552/581537443638615074354745762480279589568458618937,-4327653905728583935255115385385967640039770816512/581537443638615074354745762480279589568458618937,5910231771534213223026164067443076355213980860416/581537443638615074354745762480279589568458618937},{1054880121932942658506222715788325980197860933632/534035911503027033021959037060208471395088331547,-1281627460961420112686943019512565527775396495360/178011970501009011007319679020069490465029443849,762141108438165580312953028897043765493543993344/178011970501009011007319679020069490465029443849},{231975478494261133288862020400920/46126130931291503571272729738481,-169100958277387977102335510026312/15375376977097167857090909912827,-1},{10937629307262539547157488890951/2680744756838196321854157910496,1,9595609196525214890978161704275/9382606648933687126489552686736},{119597613426197722851154758944867/25570123189749367929383738810624,1,-6911323054769947339261622162401/2324556653613578902671248982784}};
ineqrhsPd = matrix {{844822554186341186269726002874502975268217280911/401540300577253773676892664419639444784411574272},{259488750129937520885302663064968598812539952375/138753842984171417584403220243956156918816309248},{271004736937658571397159379076252143550481633717/21660043607264033902846299161076252312355209216},{423448316443658173490576444001806454714909083095/59642064231504695724840812927973958165511274496},{-34630449549038960368831490461075377896727262675/18203973360494801078813096920473943068599386112},{3907886527366795718565487502314225730876766478645/1583399115197147713009168145831525123113336438784},{788920134410405211829833976643041257432385113137/65730678775604303964119646763830352431571533824},{3762001182033049909826590173319852249022590795111/1037305167328264277638954820505228106547542884352},{93325951612417297252353420723381015737315530589/10239641796716496655605735175750105074541527040},{896172313316093959034626764131762762974160675273/311477733472546935902002230399753906396904554496},{81236788406770564894193769464202668062514900579/7676980107802658855434829539982379117986185216},{-1},{1},{800941153401737356903768825435157970937740883067/138489084049471416744915206981758066949183504384},{589174034570303102921008420525643417063558720297/169022015231606766094385546579987883038787764224},{1219874388473548040026806220057750201529214825065/230315194538345863983135500535513255463257899008}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 19, facets: 34
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3/32,7/64,3/4,-3/16,-5/32,-3/4,-1/4,1/64,-7/8,-1,5/8,5/32,-3/4,3/64,-7/32,-3/8,3/8,-1,-1/32},{5/8,5/8,3/8,-1/2,1,1/2,-3/32,-1/4,7/256,-1/4,-3/4,-7/128,5/8,-3/8,-7/128,-7/8,5/8,3/16,1},{-3/4,-3/4,-5/8,-3/4,1/4,1/4,1,1,-1/2,-3/16,3/8,-1,7/256,7/8,1,-3/8,3/4,5/64,5/16}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{344/37,-1,365/37},{1,-29/8,31/8},{1,-604/449,-917/898},{1,-11/7,-15/14},{-1,-239/56,229/112},{0,1,-87/32},{-327/232,1,-1425/464},{-249/88,-1,-363/64},{-1,1330/597,-216/199},{-1,35528/30447,-216/199},{0,0,1},{-1,19/16,2437/1536},{-3/2,57/32,1},{-5/4,1,13/8},{-1013/704,-1,-515/352},{-3447/544,1,-28/17},{-373/128,57/32,1},{-17109/6208,283/194,-1},{-53/16,1,3},{-979/408,-1,28/17},{-73/58,-226/145,1},{-199/144,-1,-487/288},{-1,40/27,2},{-592/547,-868/547,1},{-1,-17/10,29/20},{1,145/24,-2},{0,8/3,-1},{-1,740/417,681/278},{1,6/5,451/80},{32/21,143/42,-1},{16,37/2,1},{1,170/151,-3471/1208},{872/629,-656/629,-1},{115/13,1,29/13}};
ineqrhsPd = matrix {{3037/296},{307/64},{8987/7184},{157/112},{2995/896},{341/128},{20759/7424},{14881/2816},{40361/19104},{1412371/974304},{1},{10693/6144},{145/64},{59/32},{11057/5632},{3479/544},{1703/512},{9161/3104},{239/64},{955/408},{3389/2320},{519/256},{1847/864},{6211/4376},{119/64},{517/96},{29/12},{4289/1668},{1713/320},{341/112},{293/16},{28683/9664},{377/296},{73/13}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1},{0,0,2,0,2,2,0,0,2,2,2,2,2,0,0,0},{0,4,4,0,0,4,4,0,0,4,0,4,0,4,0,4},{0,0,8,8,0,8,0,8,0,0,8,0,8,8,0,8}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{0},{8},{4},{2},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 72, facets: 68
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,4,5,-1,1,4,3,2,5,0,-5,0,-2,1,-3,4,-7,-3,-2,-4,5,-9,-4,2,5,-11,2,-1,-10,-12,2,-7,0,3,-2,-4,0,0,-9,-12,0,-12,-3,-6,2,-1,-7,-6,-7,-9,-4,-4,-1,-11,-9,-12,-9,-11,-3,-5,-8,-6,-7,-10,-5,-3,-8,-4,-7,-3,-6,-9},{-2,0,1,-4,-1,3,4,-3,-1,6,-3,-3,-3,1,-4,2,-3,7,-1,-4,2,-1,7,-2,0,0,5,7,-1,1,4,-1,-2,0,6,-1,0,4,6,4,2,2,-3,-4,4,6,6,7,1,5,6,1,3,3,-2,3,5,1,-2,6,4,-3,5,3,4,4,2,5,3,2,0,-1},{0,-1,0,0,-3,0,0,2,2,1,-1,4,-3,-5,2,-3,1,0,-5,2,3,-2,0,5,5,2,-3,-2,2,-1,4,-4,7,7,3,-6,9,6,0,0,8,2,6,4,-6,-5,-2,2,-6,2,-4,-8,-8,3,5,-3,-3,5,9,5,5,7,-5,-5,7,8,7,-7,-7,10,10,8}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{2,-1,-1},{5/4,-7/4,-1},{4/3,-2,-1},{1,-1,-1},{7/5,-1,-8/5},{7,10,1},{5/3,3,-1},{1,3,-3},{6/5,-9/5,-1},{1,-3/2,-3/2},{1,-3/2,-5/3},{1,4,1},{1,-9/2,3/2},{-1,9,-3},{0,1,0},{-1,6,-3/2},{-1,3/2,-6},{1,11/2,3/2},{-1,-13/4,-3/2},{-1,2,2},{1,-9,3},{-1,-3,3},{-1,9/2,-3/2},{-1,-3/2,-3/2},{-4/3,-1,-2},{-1,-4,-1},{-4/3,2,1},{-5/3,-3,1},{-1,3,-1},{-1,1,1},{-5/4,7/4,1},{-2,3,-1},{-2,1,1},{-1,0,0},{-6,-3,1},{-1,-1,0},{-9/2,-3,1},{-4,-11/2,-1},{-1,-1,-1},{-1,3/2,-1},{-7,-10,-1},{-7/5,1,8/5},{-1,-5/2,-1},{-6/5,9/5,1},{-1,3/2,3/2},{-1,3/2,-3/2},{-1,5,-1},{-1,3/2,5/3},{1,-2,-2},{0,-1,0},{-1,-11/2,-3/2},{1,-6,3/2},{1,-3/2,6},{1,13/4,3/2},{1,3/2,3/2},{1,-5,1},{1,-3,1},{1,5/2,1},{1,-3/2,3/2},{4/3,1,2},{9/2,3,-1},{1,-3/2,1},{2,-3,1},{1,1,1},{1,1,0},{4,11/2,1},{6,3,-1},{1,0,0}};
ineqrhsPd = matrix {{9},{6},{20/3},{5},{42/5},{61},{64/3},{32},{6},{7},{47/6},{25},{39/2},{70},{7},{46},{107/2},{71/2},{65/4},{27},{42},{36},{37},{29/2},{61/3},{18},{24},{26},{27},{17},{22},{36},{28},{12},{68},{11},{103/2},{87/2},{12},{39/2},{78},{122/5},{27/2},{109/5},{43/2},{22},{39},{68/3},{10},{4},{23},{24},{54},{22},{15},{19},{13},{16},{27/2},{18},{27},{10},{15},{10},{7},{34},{33},{5}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2,-1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2},{-1,-1,1,1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1},{-2,-2,-2,-2,2,2,2,2,2,-2,-2,-2,-2,2,2,2},{-4,-4,-4,-4,-4,4,4,4,4,4,4,4,4,-4,-4,-4}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-2,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,1,0,0},{2,0,0,0}};
ineqrhsPd = matrix {{4},{1},{1},{2},{4},{2},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,1,-1,0,0},{1,-1,0,0,1,-1},{0,0,-1,1,-1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0},{-1,-1,0},{-1,0,0},{0,1,0},{1,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = matrix {{-1,-1,-1}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,2,0,2,0,2,0,2,0,2,0,2},{0,0,2,2,4,4,0,0,2,2,4,4},{0,0,2,2,0,0,-2,-2,-4,-4,-2,-2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,-1},{0,-1,-1},{-1,0,0},{0,-1,0},{0,-1,1},{0,1,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{6},{2},{0},{0},{0},{4},{4},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-7/2,7/2,-3/2,3/2},{-5,5,-2,2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-6,4},{14,-10},{6,-4},{-14,10}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,-1,1},{-3/2,3/2,-7/2,7/2},{-2,2,-5,5}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{2,-1,0},{-5,1,0},{-2,1,0},{5,-1,0}};
ineqrhsPd = matrix {{3/2},{3/2},{3/2},{3/2}};
eqlhsPd = matrix {{1/3,4/3,-1}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,2,2,0,1,0},{-1,-1,0,1,1,0},{0,-1,-2,-1,-2,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0},{-1,-1,0},{-1,0,0},{0,1,0},{1,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{0},{0},{1},{2},{2}};
eqlhsPd = matrix {{-1,-1,-1}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0,0},{1,0,0,0,0},{0,0,1,0,2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{2,2,1},{-1,0,0},{0,-1,0}};
ineqrhsPd = matrix {{0},{2},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{1,0,0,0},{0,0,1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,0,0},{0,-1,0},{0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 27, facets: 9
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-4/29,-1/2,0,1/19,1/12,19/30,-7/12,1/2,2/3,1/6,4/29,3/19,2/3,2/29,0,0,1/2,7/12,-1/19,-2/29,-19/30,-1/12,-1/2,-2/3,-1/6,-2/3,-3/19},{2/29,-1/2,0,3/19,1/6,2/3,2/3,1/2,-19/30,-1/12,-2/29,-1/19,7/12,4/29,0,0,-1/2,-2/3,-3/19,-4/29,-2/3,-1/6,1/2,19/30,1/12,-7/12,1/19},{-5/29,2,1,3/19,1/4,19/10,-7/4,2,-19/10,-1/4,-5/29,-3/19,7/4,5/29,0,-1,-2,-7/4,3/19,5/29,19/10,1/4,-2,-19/10,-1/4,7/4,-3/19},{8/29,-1/2,0,5/19,1/4,-1/2,-1/2,-1/2,-1/2,1/4,8/29,5/19,-1/2,8/29,1/3,0,-1/2,-1/2,5/19,8/29,-1/2,1/4,-1/2,-1/2,1/4,-1/2,5/19}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-3,3,1,2},{3,-3,1,2},{2,-1,1,3},{-2,1,1,3},{3,3,-1,2},{-3,-3,-1,2},{-1,-2,-1,3},{1,2,-1,3},{0,0,0,-2}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 4, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0},{0,1},{0,1},{1,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0,0},{-1,0,0,0}};
ineqrhsPd = matrix {{1},{0}};
eqlhsPd = matrix {{-1,-1,0,0},{-1/2,1/2,-1,0},{1/3,-1/3,-1/3,-1}};
eqrhsPd = matrix {{-1},{-1/2},{-2/3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,1,-1,-1,-1,-1},{-1,-1,1,1,1,1,-1,-1},{-1,1,1,-1,1,-1,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,2,2,1,0},{0,0,2,2,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{0},{1},{0},{2},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,2,1},{0,0,1,1,1,2},{0,0,0,3,3,3}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-3,0,1},{0,0,-1},{0,0,1},{0,-3,1},{3/2,3/2,-1}};
ineqrhsPd = matrix {{0},{0},{3},{0},{3/2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,5/2,0},{0,0,5/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{5/2},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 9, vertices: 6, facets: 9
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1,0,0},{0,1,0,0,1,0},{0,0,1,0,0,1},{0,1,1,0,0,0},{1,0,0,0,0,1},{0,0,0,1,1,0},{0,0,0,0,1,1},{0,0,1,1,0,0},{1,1,0,0,0,0}};
raysP = map(QQ^9, QQ^0, 0);
linealityP = map(QQ^9, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1,0,-1,-1,0,0,0,0},{-1,0,0,0,0,0,0,0,0},{0,0,0,0,-1,0,0,0,0},{1,0,0,1,0,0,0,0,0},{1,1,0,0,0,0,0,0,0},{0,0,0,1,1,0,0,0,0},{0,1,0,0,1,0,0,0,0},{0,-1,0,0,0,0,0,0,0},{0,0,0,-1,0,0,0,0,0}};
ineqrhsPd = matrix {{-1},{0},{0},{1},{1},{1},{1},{0},{0}};
eqlhsPd = matrix {{-1,-1,-1,0,0,0,0,0,0},{0,0,0,-1,-1,-1,0,0,0},{-2/3,1/3,1/3,-2/3,1/3,1/3,-1,0,0},{1/7,-4/7,3/7,1/7,-4/7,3/7,-2/7,-1,0},{1/5,1/5,-2/5,1/5,1/5,-2/5,-2/5,-2/5,-1}};
eqrhsPd = matrix {{-1},{-1},{-1/3},{-3/7},{-3/5}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,5/2,0},{0,0,5/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{5/2},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 2, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1},{0},{0}};
raysP = matrix {{1},{0},{0}};
linealityP = matrix {{0,0},{1,0},{0,1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1/2,1/2,-1/2,1/2},{-1/2,-1/2,1/2,1/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-2},{-2,0},{0,2},{2,0}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,0,-1,-1},{0,1,1,-1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-2},{-1,1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 2, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,3},{2,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{2,1}};
ineqrhsPd = matrix {{0},{5}};
eqlhsPd = matrix {{-1,-1}};
eqrhsPd = matrix {{-2}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,0,-1,-1},{0,1,1,-1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-2},{-1,1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,0,-1},{0,1,1,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-2,1},{1,-2},{1,0},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,4}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1},{-1}};
ineqrhsPd = matrix {{4},{-1}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1,1},{-1,0,0,0,0},{0,-1,0,0,0},{0,0,-1,0,0},{0,0,0,-1,0},{0,0,0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1},{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{0},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0,0,0},{0,0,1,-1,0,0,0,0},{0,0,0,0,1,-1,0,0},{0,0,0,0,0,0,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,1,1,-1},{-1,1,1,1},{-1,1,-1,-1},{-1,1,-1,1},{-1,-1,-1,-1},{-1,-1,-1,1},{-1,-1,1,1},{-1,-1,1,-1},{1,-1,-1,1},{1,-1,-1,-1},{1,-1,1,1},{1,-1,1,-1},{1,1,-1,1},{1,1,-1,-1},{1,1,1,1},{1,1,1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,2,1,0},{0,0,1,3,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-2,1},{-1,0},{0,-1},{1,-1},{2,1}};
ineqrhsPd = matrix {{1},{0},{0},{1},{5}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3,5/3,1},{3,1,5/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-9,-4},{3,-2},{-1,4}};
ineqrhsPd = matrix {{-19},{3},{9}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 5, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3,1,5/3},{3,5/2,1},{0,0,0},{2,4,1},{3,5,2}};
raysP = matrix {{0},{0},{1},{0},{0}};
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,0,0,-2},{0,0,-1,0,0},{0,-2,0,0,1},{0,4/15,0,0,1/15}};
ineqrhsPd = matrix {{-3},{0},{0},{1}};
eqlhsPd = matrix {{-27/20,7/5,0,-1,0},{-540/1913,560/1913,0,1513/1913,-1}};
eqrhsPd = matrix {{-37/20},{-2653/1913}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,0,-1},{0,-1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 20
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,910872158600853/1125899906842624,-910872158600853/1125899906842624,1/2,-1/2,1/2,-1/2,910872158600853/1125899906842624,-910872158600853/1125899906842624,0,0},{910872158600853/1125899906842624,910872158600853/1125899906842624,1/2,1/2,0,0,0,0,-1/2,-1/2,-910872158600853/1125899906842624,-910872158600853/1125899906842624},{1/2,-1/2,0,0,910872158600853/1125899906842624,910872158600853/1125899906842624,-910872158600853/1125899906842624,-910872158600853/1125899906842624,0,0,1/2,-1/2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{391725578400125525240351555584/829688089314177501862572327609,-1125899906842624/910872158600853,0},{-391725578400125525240351555584/829688089314177501862572327609,-1125899906842624/910872158600853,0},{-1125899906842624/1473822112022165,-1125899906842624/1473822112022165,-1125899906842624/1473822112022165},{0,-391725578400125525240351555584/829688089314177501862572327609,-1125899906842624/910872158600853},{1125899906842624/1473822112022165,-1125899906842624/1473822112022165,-1125899906842624/1473822112022165},{1125899906842624/910872158600853,0,391725578400125525240351555584/829688089314177501862572327609},{1125899906842624/910872158600853,0,-391725578400125525240351555584/829688089314177501862572327609},{1125899906842624/1473822112022165,1125899906842624/1473822112022165,-1125899906842624/1473822112022165},{1125899906842624/1473822112022165,1125899906842624/1473822112022165,1125899906842624/1473822112022165},{391725578400125525240351555584/829688089314177501862572327609,1125899906842624/910872158600853,0},{0,391725578400125525240351555584/829688089314177501862572327609,1125899906842624/910872158600853},{-391725578400125525240351555584/829688089314177501862572327609,1125899906842624/910872158600853,0},{-1125899906842624/1473822112022165,1125899906842624/1473822112022165,1125899906842624/1473822112022165},{-1125899906842624/1473822112022165,1125899906842624/1473822112022165,-1125899906842624/1473822112022165},{0,391725578400125525240351555584/829688089314177501862572327609,-1125899906842624/910872158600853},{-1125899906842624/910872158600853,0,-391725578400125525240351555584/829688089314177501862572327609},{-1125899906842624/910872158600853,0,391725578400125525240351555584/829688089314177501862572327609},{-1125899906842624/1473822112022165,-1125899906842624/1473822112022165,1125899906842624/1473822112022165},{1125899906842624/1473822112022165,-1125899906842624/1473822112022165,1125899906842624/1473822112022165},{0,-391725578400125525240351555584/829688089314177501862572327609,1125899906842624/910872158600853}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,-5200308914369309/9007199254740992},{-5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,5200308914369309/9007199254740992,-5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,5200308914369309/9007199254740992},{-5200308914369309/9007199254740992,-5200308914369309/9007199254740992,-5200308914369309/9007199254740992,-5200308914369309/9007199254740992,5200308914369309/9007199254740992,5200308914369309/9007199254740992,5200308914369309/9007199254740992,5200308914369309/9007199254740992}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-9007199254740992/5200308914369309},{-9007199254740992/5200308914369309,0,0},{0,-9007199254740992/5200308914369309,0},{0,0,9007199254740992/5200308914369309},{0,9007199254740992/5200308914369309,0},{9007199254740992/5200308914369309,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,-2,2,-2,2,-2,2,-2},{-2,-2,2,2,-2,-2,2,2},{-2,-2,-2,-2,2,2,2,2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{2},{2},{2},{2},{2},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1,0,1,0,1},{3/4,1,0,1/4,0,1/4,1,3/4},{3/16,1/4,0,1/16,1,15/16,3/4,13/16}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{1/4,-1,0},{1/4,1,0},{0,1/4,-1},{0,1/4,1}};
ineqrhsPd = matrix {{0},{1},{0},{1},{0},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0},{0,1,1,-1,0,0},{0,0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,1},{-2,-1,1},{-2,-1,-1},{0,1,-1},{1,1,-1},{1,1,1},{1,-1,1},{1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,0,1,1,0,0,-4},{1,0,1,-2,1,1,-3,1},{0,0,-1,1,1,1,1,-3}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1,-2},{1,0,0},{1,0,-1},{0,1,0},{0,0,1},{-1,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{1,0,0,0},{0,0,1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,0,0},{0,-1,0},{0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 34, facets: 19
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-593934448380038187629003561336244223315615416320/581295835759584558808312033228316509435714702879,-779973149963902010794106307206482951439877406720/1148181310580490847889443245285958332842937431117,-834911336348711971508916753479457864880196943872/1033279699092212080961238485377359834950104619125,-451351638156153469675449537298017167288402182144/661074666395213201954148779083788761647422101087,-11875211446039377815154684977628597209263677374464/28146830865538587708960983240939498061618358464803,-1066404108435779557262026998475926048840047656960/2663035629934843887354436756168725949136781115679,-286667771360306766448084578481195066846882037760/1070022436834922326869526480917927343294110940591,-101857867966402892552331775300708148248320147456/313078221326677069773728886384385920072137617303,-247647508815704899421824604945096989996604194816/4715350869261984712862393959738846535578964133343,-231354719769776038518069933565989679675543126016/448263694500983649316778464851002479369071878087,1173995576268867024101222385236247095587712794624/1063363920455119185533336766028459859555986812289,10919588709302286410864254853031789456928047890432/201086830318867823857028347393140421252269033251,7737382246988530702401719335728566856035340386304/302903591957417203012304489947529736347874762775,-208965238551936641108002108112096334628087922688/969015032110259826352173485095262980905288420413,-859290202182082173684458899037542304891078180864/5007207875503957243795218944333844953080912307477,-10197621225091634689632545083717720162533638144/34796618359701895270882545993358473009022808147,-164860432541172988012500288848677826630810337280/369527502856320347152283004429443120261713839807,-4788459856255095265270714418435572229295055044608/17794019300086801999024109977398818324480077290739,-121768151542727502417327058247524448316084453376/214725024102500251121085398955994353664526185507,-4695507374188243411505217745040521455703542464512/8239742805158896400827962926559115961918687751871,-85837091223178700317543387027914465534606835712/134535282588166625496077801426867389035530698279,-580435442581073966343298317053478125353981444096/601721064087771031599502474452759076861582344741,-4399086184007617290603389456832361923188004749312/5587846899827315228018162025780624481889822155377,-8946398221587728817626161242544128840325987827712/11324914804645713522973305338138304500811523904485,-2380243102799262291488880640371567860996900913152/2500718746193551593444075481086648309079292861499,-2101645623657309528239208136369819998420526432256/2575644869312600520837378718018929479143740845201,-1478893849222614588217813826634995473583454552064/1537865163739176575553792502864197330470575328727,-1650966813993801114766861815173619324152071585792/1477194555638248894623263976557986753178854006673},{-22818873475658756052383747193627849323869896704/581295835759584558808312033228316509435714702879,976998916289555033383810311596041748929302233088/1148181310580490847889443245285958332842937431117,7586383608487890572942933955429067934453511225344/9299517291829908728651146368396238514550941572125,1714126851978983796318065960745923606019303276544/1983223999185639605862446337251366284942266303261,9175789268262439599570487075487982327846387318784/28146830865538587708960983240939498061618358464803,2118822138322137696317075292329089245204883439616/2663035629934843887354436756168725949136781115679,-1035135536219048502230860487599639901808833331200/1070022436834922326869526480917927343294110940591,-2151310265872916459076542960878097013057267957760/2191547549286739488416102204690701440504963321121,-4629211660552619529874925514650714422269913858048/4715350869261984712862393959738846535578964133343,-30783432736003875979750701542967690016264814592/64037670642997664188111209264428925624153125441,-293529744116052320987137638985152910776380424192/1063363920455119185533336766028459859555986812289,1818696859734868465859909813538812606321256300544/201086830318867823857028347393140421252269033251,-8797965056266683176252475579552797508633440026624/2726132327616754827110740409527767627130872864975,950884402381782865481050691257735479637701558272/969015032110259826352173485095262980905288420413,4910313066152228029763916645661983662521781846016/5007207875503957243795218944333844953080912307477,200489647855657583570731698405731449071592275968/521949275395528429063238189900377095135342122205,-264064968691283956425736981033314113190331678720/369527502856320347152283004429443120261713839807,3237718166968952441032952974617210505523697287168/17794019300086801999024109977398818324480077290739,-690405230055092879748857330835757135958038806528/1503075168717501757847597792691960475651683298549,-3750862582827464504403021450241614596972344246272/8239742805158896400827962926559115961918687751871,-103079032796038331763922255286326302856768913408/134535282588166625496077801426867389035530698279,-399722854162767122647551157340747551288997183488/6618931704965481347594527218980349845477405792151,-4191403943173474617574437652497846669758179573760/5587846899827315228018162025780624481889822155377,-8544545583914319080093878921340161755887043608576/11324914804645713522973305338138304500811523904485,-1570702742277753177347654238016256171368196866048/2500718746193551593444075481086648309079292861499,3067136313122161735618404786018771451266004418560/7726934607937801562512136154056788437431222535603,8017251242771651453532827175785619515969990819840/13840786473652589179984132525777775974235177958543,-210009175232972847106184266082692189666278375424/1477194555638248894623263976557986753178854006673},{292015572736985309319079896593271243862274211840/581295835759584558808312033228316509435714702879,406448713537773738248476345065189171331470458880/1148181310580490847889443245285958332842937431117,-2419017046828958004959018541626998047413981151232/9299517291829908728651146368396238514550941572125,605841361880050239986422753184209142956366495744/1983223999185639605862446337251366284942266303261,26376580871057029152502082943918053696927447384064/28146830865538587708960983240939498061618358464803,1829625818811552520540966823773129831638487793664/2663035629934843887354436756168725949136781115679,355357895910219463637798838165319385714427691008/1070022436834922326869526480917927343294110940591,-400465054874798630522403667777820192393459138560/2191547549286739488416102204690701440504963321121,2282559012109495953069810963444182641650457116672/4715350869261984712862393959738846535578964133343,361986275197950306925085973513804365362903384064/448263694500983649316778464851002479369071878087,4294441423989885021499182554914770621993064595456/3190091761365357556600010298085379578667960436867,-122506588092093177848048518135056155612152856576/201086830318867823857028347393140421252269033251,470849944371399163319599630506760597200342876160/109045293104670193084429616381110705085234914599,127886593915226235025779439534246475699425640448/969015032110259826352173485095262980905288420413,-1505216801932142153130468954670392960833913094144/5007207875503957243795218944333844953080912307477,-613844172024470845003066640846993936060801613824/521949275395528429063238189900377095135342122205,-254819809966568204935607306891143779378352095232/369527502856320347152283004429443120261713839807,127281167770811709473005771620679312072431894528/124433701399208405587581188653138589681678862173,1197126379723621998890849251275759964883179798528/1503075168717501757847597792691960475651683298549,6562719208971136758614765103246540865130763124736/8239742805158896400827962926559115961918687751871,42417100986611552576579593815379552102260932608/134535282588166625496077801426867389035530698279,3814345438454256446950114483520119135859390808064/6618931704965481347594527218980349845477405792151,-2796281770744980650596513001246509959046460604416/5587846899827315228018162025780624481889822155377,-1110906015150491879853805652350360710054221447168/2264982960929142704594661067627660900162304780897,-162999463660340098562820645803440974937427279872/2500718746193551593444075481086648309079292861499,-7039471309045169588769443884563022662595116007424/7726934607937801562512136154056788437431222535603,-7720896340977146716264447239622474994632773074944/13840786473652589179984132525777775974235177958543,319771745064582453583006035182977561718765387776/1477194555638248894623263976557986753178854006673}};
raysP = matrix {{1,1,1,1,1,1},{251112891175366054364274580161127/1742729142036330334599966324909420,-34025299559200432416283693893587/344714124034586031771472840661500,-166072625277962181545693392506793/1856806596362728607273125232693112,128406741685128773958560033805701/869050463315085459389094472964712,6841838086192401798777624145777/204686637756790519442275767776536,2131537604361729521355148099871/57850743173374257358441466703168},{-28909925018099237201050958950312/118822441502477068268179522152915,-11929325663257510348818223910588/86178531008646507942868210165375,24333199372073205901955950761380/232100824545341075909140654086639,-7193367473035672524891855261596/325893923743157047270910427361767,-1842628407698070068212035118957/4498607423226165262467599291792,-47544614283236920077152578537201/115701486346748514716882933406336}};
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-4741371086824893/9007199254740992,-7647497920718555/144115188075855872,7643335962122821/9007199254740992},{-1078387572567719/2251799813685248,-1969043944499017/2251799813685248,1398252369236353/18014398509481984},{-4545080756734037/18014398509481984,8712378752251247/9007199254740992,-3907903642428347/144115188075855872},{-3812755930250863/9007199254740992,-6613900702997815/9007199254740992,597505859304803/1125899906842624},{-528415934779725/562949953421312,-4929203444735165/36028797018963968,-5702589786781073/18014398509481984},{-7844636121648267/18014398509481984,-5639193008728887/18014398509481984,-475143625477287/562949953421312},{-8197550089758337/9007199254740992,1703653972479579/4503599627370496,3046350937582833/18014398509481984},{-3615904728084257/72057594037927936,4243949184108969/9007199254740992,7931849752919691/9007199254740992},{-3034309338522589/18014398509481984,-5763735044348673/9007199254740992,6753313147092495/9007199254740992},{-5164178031720289/18014398509481984,-7125750769402285/36028797018963968,8443277472466113/9007199254740992},{-8820475294483955/36028797018963968,7356267017525943/9007199254740992,-2353318133783607/4503599627370496},{-5322465915783559/36028797018963968,-8008288277152163/9007199254740992,-1951052092499811/4503599627370496},{-6143462652402181/9007199254740992,1656098959042419/4503599627370496,5693579717943775/9007199254740992},{-781114321112795/2251799813685248,-6104924401936087/18014398509481984,-7877174595240857/9007199254740992},{-6725670116260403/9007199254740992,-4552485229651157/9007199254740992,7789705549160769/18014398509481984},{-2639153864504731/18014398509481984,8908792931407259/9007199254740992,-4720801535434569/288230376151711744},{-3349280896922509/36028797018963968,-560106062785759/562949953421312,5461605447630899/144115188075855872},{-4940367677971645/36028797018963968,2144252301938057/2251799813685248,4914915834536135/18014398509481984},{0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0,1,0},{0,0,1,0,0,1},{0,0,0,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{0,0,-1},{0,0,1},{0,-1,0},{1,1,0}};
ineqrhsPd = matrix {{0},{0},{1},{0},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0},{0,0,1,-1,0,0},{0,0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,3,0,0,1},{0,0,3,0,1},{0,0,0,3,-5}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-5,0,-1},{-1,0,0},{1,1,1},{5,5,-1},{0,-1,0},{0,-5,-1}};
ineqrhsPd = matrix {{0},{0},{3},{15},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 10, facets: 12
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,1},{-1,-1,-1,-1,1,1,1,1,-1,1},{0,0,0,0,0,0,0,0,2,-2}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,1},{0,-1,0,0},{0,1,0,1},{0,0,-1,0},{0,0,1,1},{-1,0,0,-1},{1,0,0,0},{0,-1,0,-1},{0,1,0,0},{0,0,-1,-1},{0,0,1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,3,0},{0,0,3}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{3},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,3,0,1,1},{0,0,3,1,1},{0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,-1},{-1,0,1},{1,1,1},{1,1,-1},{0,-1,1},{0,-1,-1}};
ineqrhsPd = matrix {{0},{0},{3},{3},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2,-1/2,3/4,2,100},{2/3,2/3,2/3,3,200},{3/4,-3/4,1/2,4,300}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-3/2,0},{28/25,27/25,-28/25},{28/3,-19,28/3},{121/20,-753/50,401/50},{1,3294/2941,-3166/2941},{-299/100,-3/2,299/150}};
ineqrhsPd = matrix {{-1},{1},{-1},{-1},{3100/2941},{-1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0,0,0},{0,0,1,-1,0,0,0,0},{0,0,0,0,1,-1,0,0},{0,0,0,0,0,0,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,1,1,1},{1,-1,1,1},{-1,-1,1,1},{1,1,-1,1},{-1,1,-1,1},{1,-1,-1,1},{-1,-1,-1,1},{1,1,1,-1},{-1,1,1,-1},{1,-1,1,-1},{-1,-1,1,-1},{1,1,-1,-1},{-1,1,-1,-1},{1,-1,-1,-1},{-1,-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0},{0,0,1,-1,0,0},{0,0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 10, facets: 7
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1,1,1,1,-1,-1,-1,-1},{1,0,1,-1,0,1,1,-3,-1,1},{-1,-1,1,1,0,0,1,1,-1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1,-1},{1,0,0},{1,0,-1},{0,1,0},{0,0,-1},{0,0,1},{-1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,1,0,-1,-1},{1,1,0,-1,-1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1},{0,-1},{-1,0},{-1,1},{0,1},{1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0},{0,0,1,-1,0,0},{0,0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,-1,1,-1,1,-1,-1},{1,-1,1,1,-1,-1,1,-1},{1,-1,-1,-1,1,1,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 10, facets: 9
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,1,0,0,1,1,0,1},{1,1,1,0,0,0,0,1,0,1},{0,0,1,0,1,1,0,0,0,0},{0,0,1,0,1,1,0,1,1,1},{0,1,1,1,0,1,0,0,1,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0,0},{-1,0,-1,0,-1},{-1,0,0,-1,0},{1,0,1,0,-1},{0,0,1,-1,0},{1,0,0,0,0},{0,0,0,0,1},{1,-1,0,1,0},{-1,1,0,0,0}};
ineqrhsPd = matrix {{0},{-1},{-1},{1},{0},{1},{1},{1},{0}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,-1,0,0},{0,1,0,0,-1,0},{0,0,1,0,0,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 10, facets: 9
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,1,0,1,1,1,0,1,1},{1,0,0,0,1,1,0,0,1,1},{0,1,0,0,0,0,0,1,0,1},{1,1,0,1,0,1,0,1,0,1},{0,0,0,1,0,1,1,1,1,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0,0},{-1,0,-1,0,-1},{-1,0,0,-1,0},{1,0,1,0,-1},{0,0,1,-1,0},{1,0,0,0,0},{0,0,0,0,1},{1,-1,0,1,0},{-1,1,0,0,0}};
ineqrhsPd = matrix {{0},{-1},{-1},{1},{0},{1},{1},{1},{0}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 24, facets: 24
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,1,0,0,1,1,0,1,0,-1,-1,0,0,1,0,0,-1,0,-1,0,0,-1},{0,0,0,-1,0,0,0,-1,0,0,-1,0,0,-1,0,-1,0,-1,1,1,1,1,1,1},{-1,1,1,0,0,0,0,0,1,-1,1,1,0,0,-1,-1,-1,0,0,0,1,-1,0,0},{-1,1,-1,-1,-1,1,-1,0,0,0,0,0,1,0,0,0,1,1,1,-1,0,0,0,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0},{1,1,0,1},{1,0,0,1},{0,-1,0,0},{1,0,1,1},{-1,-1,-1,0},{-1,0,0,-1},{0,0,1,0},{0,0,0,1},{-1,-1,0,0},{-1,0,0,0},{-2,-1,-1,-1},{-1,0,-1,0},{1,1,1,1},{0,1,0,0},{-1,-1,0,-1},{-1,-1,-1,-1},{-1,0,-1,-1},{0,0,0,-1},{1,1,1,0},{1,0,1,0},{2,1,1,1},{1,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 24, facets: 24
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,0,0,0,0,0,0,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1},{0,-1/2,-1/2,-1/2,-1/2,1/2,1/2,1/2,1/2,-1,0,0,0,0,1,-1/2,-1/2,-1/2,-1/2,1/2,1/2,1/2,1/2,0},{0,-1/2,-1/2,1/2,1/2,-1/2,-1/2,1/2,1/2,0,-1,0,0,1,0,-1/2,-1/2,1/2,1/2,-1/2,-1/2,1/2,1/2,0},{0,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,0,0,-1,1,0,0,-1/2,1/2,-1/2,1/2,-1/2,1/2,-1/2,1/2,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,1,0},{0,0,1,1},{0,1,0,1},{0,1,-1,0},{0,0,-1,1},{0,-1,1,0},{0,-1,0,1},{0,-1,-1,0},{-1,1,0,0},{-1,0,0,1},{-1,0,-1,0},{-1,-1,0,0},{-1,0,0,-1},{-1,0,1,0},{0,-1,0,-1},{0,0,-1,-1},{0,1,0,-1},{0,0,1,-1},{1,0,0,-1},{1,-1,0,0},{1,0,-1,0},{1,0,0,1},{1,0,1,0},{1,1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 24, facets: 24
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,-1,-1},{-1,0,-1,-1,0,1,1,1,-1,0,0,1,0,0,1,0,-1,0,0,1,0,0,0,-1},{-1,-1,0,1,1,-1,0,1,0,0,0,0,1,-1,0,1,0,-1,0,0,1,0,-1,0},{0,-1,-1,0,-1,0,-1,0,0,-1,1,0,0,0,1,1,1,1,1,0,0,-1,0,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{-1/2,-1/2,-1/2,-1/2},{-1/2,-1/2,-1/2,1/2},{-1/2,-1/2,1/2,-1/2},{-1/2,-1/2,1/2,1/2},{-1/2,1/2,-1/2,-1/2},{-1/2,1/2,-1/2,1/2},{-1/2,1/2,1/2,-1/2},{-1/2,1/2,1/2,1/2},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1},{0,0,0,1},{0,0,1,0},{0,1,0,0},{1/2,-1/2,-1/2,-1/2},{1/2,-1/2,-1/2,1/2},{1/2,-1/2,1/2,-1/2},{1/2,-1/2,1/2,1/2},{1/2,1/2,-1/2,-1/2},{1/2,1/2,-1/2,1/2},{1/2,1/2,1/2,-1/2},{1/2,1/2,1/2,1/2},{1,0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 7, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,9007199254740992/7286977268806823,-1,29472066294954575717671549009920/77158871278881291906588116670641,0,0,29472066294954575717671549009920/77158871278881291906588116670641},{6223821613304664/8566355544790271,0,-6223821613304664/8566355544790271,-90705693242519399071460246945792/77158871278881291906588116670641,0,0,90705693242519399071460246945792/77158871278881291906588116670641},{0,0,0,0,1,-1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,1},{-347922205179541/1125899906842624,-8566355544790271/9007199254740992,1},{7286977268806823/9007199254740992,-5294298886396511/9007199254740992,1},{7286977268806823/9007199254740992,5294298886396511/9007199254740992,1},{-347922205179541/1125899906842624,8566355544790271/9007199254740992,1},{-1,0,-1},{-347922205179541/1125899906842624,-8566355544790271/9007199254740992,-1},{7286977268806823/9007199254740992,-5294298886396511/9007199254740992,-1},{7286977268806823/9007199254740992,5294298886396511/9007199254740992,-1},{-347922205179541/1125899906842624,8566355544790271/9007199254740992,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{1}};
raysP = matrix {{-1,0},{0,-1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0},{0,1}};
ineqrhsPd = matrix {{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 4, vertices: 1, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1/6},{-1/6},{-1/6},{-1/6}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^4, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = matrix {{0,-1,-2,-3},{0,-1/14,-8/7,11/14},{0,-10/9,2/9,2/9},{-1,0,0,0}};
eqrhsPd = matrix {{1},{1/14},{1/9},{1/6}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 7, vertices: 35, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0},{1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0},{1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,1},{0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,1,1,0,1,1,0,1,1,1,0,1,1},{0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1},{0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,1,0,1,1,1,0,1,1,1,1}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0,0,0,0,0},{-1,0,0,0,0,0,0},{0,1,0,0,0,0,0},{0,-1,0,0,0,0,0},{0,0,1,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,1,0,0,0},{0,0,0,-1,0,0,0},{0,0,0,0,1,0,0},{0,0,0,0,-1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,-1,0},{0,0,0,0,0,0,1},{0,0,0,0,0,0,-1}};
ineqrhsPd = matrix {{1},{0},{1},{0},{1},{0},{1},{0},{1},{0},{1},{0},{1},{0}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-4}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 6, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-5,-5,0,0,5,6},{0,6,-5,5,-5,0},{-5,0,6,-5,0,5},{6,5,0,0,-5,-5},{0,-5,5,-5,6,0},{5,0,-5,6,0,-5}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-211,21,20,-1176/5,-1,0},{1160/121,105/11,-1,1050/121,1055/121,0},{-1,-1281/5,-22,-21,-232,0},{1171/121,1176/121,1,1281/121,116/11,0},{1,-210,22,21,-1171/5,0},{-1276/5,-21,-20,-231,1,0}};
ineqrhsPd = matrix {{10},{1105/121},{-11},{1226/121},{11},{-10}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 6, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-5,-5,0,0,5,6},{0,6,-5,5,-5,0},{-5,0,6,-5,0,5},{6,5,0,0,-5,-5},{0,-5,5,-5,6,0},{5,0,-5,6,0,-5}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-211,21,20,-1176/5,-1,0},{1160/121,105/11,-1,1050/121,1055/121,0},{-1,-1281/5,-22,-21,-232,0},{1171/121,1176/121,1,1281/121,116/11,0},{1,-210,22,21,-1171/5,0},{-1276/5,-21,-20,-231,1,0}};
ineqrhsPd = matrix {{10},{1105/121},{-11},{1226/121},{11},{-10}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 4, vertices: 6, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,1,1,1},{1,1,1,0,0,0},{1,1,2,1,1,2},{1,2,1,1,2,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{1,-1,0,0},{0,0,-1,0},{-1,1,0,0},{0,0,1,1}};
ineqrhsPd = matrix {{-1},{1},{-1},{1},{3}};
eqlhsPd = matrix {{-1,-1,0,0}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 5, vertices: 10, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,1,1,0,0,0,1,1},{0,0,1,0,1,0,1,1,0,0},{0,1,1,0,0,1,0,0,0,1},{1,1,0,0,0,0,0,1,1,0},{1,0,0,1,0,1,1,0,0,0}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0,0,0},{0,0,1,0,0},{-1,0,0,0,0},{-1,-1,-1,-1,0},{0,0,0,1,0},{0,1,0,0,0},{1,0,0,0,0},{0,0,-1,0,0},{0,0,0,-1,0},{1,1,1,1,0}};
ineqrhsPd = matrix {{0},{1},{0},{-1},{1},{1},{1},{0},{0},{2}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-2}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,-1,0},{1,0,0,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1},{-1,1},{1,1},{1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 5, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1},{0,0,1,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0},{0,-1,0,0,0},{0,1,0,0,0},{1,0,0,0,0}};
ineqrhsPd = matrix {{0},{0},{1},{1}};
eqlhsPd = matrix {{0,0,-1,0,0},{0,0,0,-1,0},{0,0,0,0,-1}};
eqrhsPd = matrix {{0},{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 6, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-5,-5,0,0,5,6},{0,6,-5,5,-5,0},{-5,0,6,-5,0,5},{6,5,0,0,-5,-5},{0,-5,5,-5,6,0},{5,0,-5,6,0,-5}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-211,21,20,-1176/5,-1,0},{1160/121,105/11,-1,1050/121,1055/121,0},{-1,-1281/5,-22,-21,-232,0},{1171/121,1176/121,1,1281/121,116/11,0},{1,-210,22,21,-1171/5,0},{-1276/5,-21,-20,-231,1,0}};
ineqrhsPd = matrix {{10},{1105/121},{-11},{1226/121},{11},{-10}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1},{0,1,0},{1,0,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,0},{-1,0,0},{0,-1,0}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = matrix {{-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 4, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^4, QQ^0, 0);
raysP = map(QQ^4, QQ^0, 0);
linealityP = matrix {{-1},{-1},{1},{1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^4, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = matrix {{0,-1,0,-1},{0,-1/2,-1,1/2},{-1,1/3,-1/3,-1/3},{0,0,0,0}};
eqrhsPd = matrix {{1},{1/2},{2/3},{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,0,0,1},{0,-1,1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1},{-1,1},{1,1},{1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1},{0,1,0},{1,0,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,0},{-1,0,0},{0,-1,0}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = matrix {{-1,-1,-1}};
eqrhsPd = matrix {{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-1,1,0,-1,1,0,-1,-1,1,1,0},{-1,0,0,1,-1,-1,-1,1,0,1,0,1},{-1,-1,-1,-1,0,0,1,0,1,0,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,1},{1,-1,1},{0,-1,0},{-1,1,-1},{-1,0,0},{-1,-1,-1},{-1,-1,1},{-1,1,1},{1,-1,-1},{0,0,-1},{1,1,-1},{1,1,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{2},{1},{2},{1},{2},{2},{2},{2},{1},{2},{2},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 6, vertices: 64, facets: 12
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0,0},{1,0,0,0,0,0},{0,-1,0,0,0,0},{0,1,0,0,0,0},{0,0,-1,0,0,0},{0,0,1,0,0,0},{0,0,0,-1,0,0},{0,0,0,1,0,0},{0,0,0,0,-1,0},{0,0,0,0,1,0},{0,0,0,0,0,-1},{0,0,0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^6, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,1,-1,0},{-1,-1,1,1,2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{1,0},{0,-1},{-1,0},{-1,1}};
ineqrhsPd = matrix {{2},{1},{1},{1},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0},{0,0,1,-1,0,0},{0,0,0,0,1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 4, vertices: 24, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4},{2,2,3,3,4,4,1,1,3,3,4,4,1,1,2,2,4,4,1,1,2,2,3,3},{3,4,2,4,2,3,3,4,1,4,1,3,2,4,1,4,1,2,2,3,1,3,1,2},{4,3,4,2,3,2,4,3,4,1,3,1,4,2,4,1,2,1,3,2,3,1,2,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1,0},{0,-1,-1,0},{-1,-1,-1,0},{-1,0,-1,0},{-1,0,0,0},{-1,-1,0,0},{0,1,1,0},{0,0,1,0},{0,1,0,0},{0,-1,0,0},{1,0,1,0},{1,1,1,0},{1,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{-1},{-3},{-6},{-3},{-1},{-3},{7},{4},{4},{-1},{7},{9},{7},{4}};
eqlhsPd = matrix {{-1,-1,-1,-1}};
eqrhsPd = matrix {{-10}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 48, facets: 26
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,3,-3,3,-3,3,-3,3,-3,3,-3,3,-3,3,-3,3,-3},{2,2,-2,-2,2,2,-2,-2,3,3,-3,-3,3,3,-3,-3,1,1,-1,-1,1,1,-1,-1,3,3,-3,-3,3,3,-3,-3,1,1,-1,-1,1,1,-1,-1,2,2,-2,-2,2,2,-2,-2},{3,3,3,3,-3,-3,-3,-3,2,2,2,2,-2,-2,-2,-2,3,3,3,3,-3,-3,-3,-3,1,1,1,1,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,1,1,1,1,-1,-1,-1,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,-1},{1,0,-1},{1,-1,-1},{0,0,-1},{0,1,-1},{0,-1,-1},{-1,1,-1},{-1,-1,-1},{-1,0,-1},{-1,0,0},{-1,-1,0},{-1,0,1},{-1,1,0},{-1,-1,1},{-1,1,1},{0,-1,1},{0,1,1},{0,-1,0},{0,0,1},{0,1,0},{1,-1,0},{1,-1,1},{1,0,1},{1,1,1},{1,1,0},{1,0,0}};
ineqrhsPd = matrix {{6},{5},{6},{3},{5},{5},{6},{6},{5},{3},{5},{5},{5},{6},{6},{5},{5},{3},{3},{3},{5},{6},{5},{6},{5},{3}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 5, vertices: 120, facets: 30
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5},{2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,1,1,1,1,1,1,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,1,1,1,1,1,1,2,2,2,2,2,2,4,4,4,4,4,4,5,5,5,5,5,5,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,5,5,5,5,5,5,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4},{3,3,4,4,5,5,2,2,4,4,5,5,2,2,3,3,5,5,2,2,3,3,4,4,3,3,4,4,5,5,1,1,4,4,5,5,1,1,3,3,5,5,1,1,3,3,4,4,2,2,4,4,5,5,1,1,4,4,5,5,1,1,2,2,5,5,1,1,2,2,4,4,2,2,3,3,5,5,1,1,3,3,5,5,1,1,2,2,5,5,1,1,2,2,3,3,2,2,3,3,4,4,1,1,3,3,4,4,1,1,2,2,4,4,1,1,2,2,3,3},{4,5,3,5,3,4,4,5,2,5,2,4,3,5,2,5,2,3,3,4,2,4,2,3,4,5,3,5,3,4,4,5,1,5,1,4,3,5,1,5,1,3,3,4,1,4,1,3,4,5,2,5,2,4,4,5,1,5,1,4,2,5,1,5,1,2,2,4,1,4,1,2,3,5,2,5,2,3,3,5,1,5,1,3,2,5,1,5,1,2,2,3,1,3,1,2,3,4,2,4,2,3,3,4,1,4,1,3,2,4,1,4,1,2,2,3,1,3,1,2},{5,4,5,3,4,3,5,4,5,2,4,2,5,3,5,2,3,2,4,3,4,2,3,2,5,4,5,3,4,3,5,4,5,1,4,1,5,3,5,1,3,1,4,3,4,1,3,1,5,4,5,2,4,2,5,4,5,1,4,1,5,2,5,1,2,1,4,2,4,1,2,1,5,3,5,2,3,2,5,3,5,1,3,1,5,2,5,1,2,1,3,2,3,1,2,1,4,3,4,2,3,2,4,3,4,1,3,1,4,2,4,1,2,1,3,2,3,1,2,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1,0},{0,0,-1,-1,0},{0,-1,-1,-1,0},{0,-1,0,-1,0},{-1,-1,-1,-1,0},{-1,0,-1,-1,0},{-1,-1,0,-1,0},{-1,0,0,-1,0},{-1,0,0,0,0},{-1,-1,0,0,0},{-1,0,-1,0,0},{0,1,1,1,0},{-1,-1,-1,0,0},{0,0,1,1,0},{0,1,0,1,0},{0,1,1,0,0},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{0,-1,0,0,0},{0,-1,-1,0,0},{1,0,1,1,0},{1,0,0,1,0},{1,0,1,0,0},{0,0,-1,0,0},{1,1,0,1,0},{1,1,1,1,0},{1,1,1,0,0},{1,1,0,0,0},{1,0,0,0,0}};
ineqrhsPd = matrix {{-1},{-3},{-6},{-3},{-10},{-6},{-6},{-3},{-1},{-3},{-3},{12},{-6},{9},{9},{9},{5},{5},{5},{-1},{-3},{12},{9},{9},{-1},{12},{14},{12},{9},{5}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-15}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^4, QQ^0, 0);
raysP = matrix {{-1,-1,-1,-1},{-1,0,0,-1},{-1,-1,0,0},{0,0,0,0}};
linealityP = matrix {{1},{1},{1},{1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1,0,0},{1,0,-1,0},{0,1,0,-1},{0,0,1,-1}};
ineqrhsPd = matrix {{0},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1,-1},{0,1,-1,-1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{0,1},{-2,-1},{1,1},{1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 2, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2},{1/2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0},{0,-1}};
eqrhsPd = matrix {{-1/2},{-1/2}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 3, vertices: 2, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0},{-1},{-1}};
raysP = matrix {{1},{0},{0}};
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0}};
ineqrhsPd = matrix {{0}};
eqlhsPd = matrix {{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{1},{1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: -1, ambientDim: 4, vertices: 0, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^4, QQ^0, 0);
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^4, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 3, vertices: 2, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0},{1},{1}};
raysP = matrix {{1},{0},{0}};
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0}};
ineqrhsPd = matrix {{0}};
eqlhsPd = matrix {{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{-1},{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 7, vertices: 7, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3,37/6,1,79/48},{0,43/6,2,65/48},{4,0,0,31/24},{0,31/3,0,0},{65/2,0,0,0},{0,0,31/2,0},{1,-52/3,-7,-53/12}};
raysP = matrix {{1,1,0},{1,0,0},{0,1,0},{2,1,0},{3,27/2,1},{0,0,1},{-2,1,0}};
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,-1,0,0,0},{0,0,0,0,-1,0,0},{0,0,0,0,0,-1,0}};
ineqrhsPd = matrix {{0},{0},{0},{0},{0}};
eqlhsPd = matrix {{-1,1,1,0,0,0,0},{-1/3,-5/3,4/3,0,0,0,-1},{15/17,7/17,8/17,-1,0,0,-6/17},{87/26,-1/26,44/13,43/13,-1,1,45/13}};
eqrhsPd = matrix {{1},{10/3},{71/17},{-71/13}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 7, vertices: 4, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{0},{0},{4},{43/2},{0},{3}};
raysP = matrix {{1,1,0},{0,1,0},{1,0,0},{1,2,0},{27/2,3,1},{0,0,1},{1,-2,0}};
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,0,0,-1,0}};
ineqrhsPd = matrix {{0},{0},{0}};
eqlhsPd = matrix {{-1,1,1,0,0,0,0},{-1/3,-5/3,4/3,0,0,0,-1},{15/17,7/17,8/17,-1,0,0,-6/17},{87/26,-1/26,44/13,43/13,-1,1,45/13}};
eqrhsPd = matrix {{-1},{-10/3},{-71/17},{71/13}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 6, vertices: 4, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0},{1},{0},{5},{0},{1}};
raysP = matrix {{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}};
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0,0},{0,0,-1,0,0,0},{0,0,0,0,-1,0}};
ineqrhsPd = matrix {{0},{0},{0}};
eqlhsPd = matrix {{-1,1,0,0,0,0},{0,0,-1,1,0,0},{0,0,0,0,-1,1}};
eqrhsPd = matrix {{1},{5},{1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0},{0,0,1}};
raysP = matrix {{1},{1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{-1,1},{1,-1},{0,-1}};
ineqrhsPd = matrix {{0},{1},{1},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 6, vertices: 4, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{0},{5},{0},{1},{0}};
raysP = matrix {{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}};
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1,0,0,0,0},{0,0,0,-1,0,0},{0,0,0,0,0,-1}};
ineqrhsPd = matrix {{0},{0},{0}};
eqlhsPd = matrix {{-1,1,0,0,0,0},{0,0,-1,1,0,0},{0,0,0,0,-1,1}};
eqrhsPd = matrix {{-1},{-5},{-1}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 4, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{29,7,41/2,40,40,40,40,40},{7,7,21/2,21/2,7,21/2,7,245/24},{22,0,0,39/2,0,0,11,165/8},{2,2,3,3,2,3,2,35/12}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,1,1,10},{1,-3,1,0},{1,0,0,0},{0,0,-1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{20},{30},{40},{0},{-2},{3}};
eqlhsPd = matrix {{0,-1,0,7/2}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,1,-1,-2,-1,1},{0,1,1,0,-1,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1},{-1,-1},{-1,1},{0,-1},{1,-1},{1,1}};
ineqrhsPd = matrix {{1},{2},{2},{1},{2},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-2,-1},{0,0,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,-1},{0,1},{-1,-1}};
ineqrhsPd = matrix {{0},{0},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 6, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,1,-1,-2,-1,1},{0,1,1,0,-1,-1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1},{-1,-1},{-1,1},{0,-1},{1,-1},{1,1}};
ineqrhsPd = matrix {{1},{2},{2},{1},{2},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-1,-2},{0,1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{0,-1},{-1,1}};
ineqrhsPd = matrix {{0},{0},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{0,0,1,0},{0,0,0,2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{0,-1,0},{0,0,-1},{2,2,1}};
ineqrhsPd = matrix {{0},{0},{0},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0,0,0},{0,0,1,-1,0,0,0,0},{0,0,0,0,1,-1,0,0},{0,0,0,0,0,0,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,1,1,1},{1,-1,1,1},{-1,-1,1,1},{1,1,-1,1},{-1,1,-1,1},{1,-1,-1,1},{-1,-1,-1,1},{1,1,1,-1},{-1,1,1,-1},{1,-1,1,-1},{-1,-1,1,-1},{1,1,-1,-1},{-1,1,-1,-1},{1,-1,-1,-1},{-1,-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{5/8,-5/8,0,0,0,0,0,0},{0,0,5/8,-5/8,0,0,0,0},{0,0,0,0,5/8,-5/8,0,0},{0,0,0,0,0,0,5/8,-5/8}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,1,1,1},{1,-1,1,1},{-1,-1,1,1},{1,1,-1,1},{-1,1,-1,1},{1,-1,-1,1},{-1,-1,-1,1},{1,1,1,-1},{-1,1,1,-1},{1,-1,1,-1},{-1,-1,1,-1},{1,1,-1,-1},{-1,1,-1,-1},{1,-1,-1,-1},{-1,-1,-1,-1}};
ineqrhsPd = matrix {{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8},{5/8}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1,-1,0,0},{0,1,0,0,-1,0},{1,0,0,0,0,-1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,1,-1},{-1,1,1},{1,1,1},{-1,-1,1},{1,-1,1},{1,1,-1},{-1,-1,-1},{1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 24, facets: 26
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1/4,-1/4,-1/4,-1/4,1/4,1/4,1/4,1/4,1/2,1/2,1/2,1/2,1/4,1/4,1/4,1/4,-1/4,-1/4,-1/4,-1/4,-1/2,-1/2,-1/2,-1/2},{-1/2,-1/2,-1/4,1/4,-1/2,-1/2,-1/4,1/2,-1/4,1/4,1/4,-1/4,1/2,1/4,1/4,-1/4,1/2,1/2,1/4,-1/4,1/4,1/4,-1/4,-1/4},{1/4,-1/4,1/2,1/2,-1/4,1/4,1/2,-1/4,1/4,1/4,-1/4,-1/4,1/4,1/2,-1/2,-1/2,1/4,-1/4,-1/2,-1/2,-1/4,1/4,1/4,-1/4}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1},{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},{0,0,2},{0,0,-2},{4/3,4/3,0},{4/3,-4/3,0},{4/3,0,4/3},{4/3,0,-4/3},{-4/3,4/3,0},{-4/3,-4/3,0},{-4/3,0,4/3},{-4/3,0,-4/3},{0,4/3,4/3},{0,4/3,-4/3},{0,-4/3,4/3},{0,-4/3,-4/3}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 3, vertices: 1, facets: 0
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2},{0},{0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^0, QQ^3, 0);
ineqrhsPd = map(QQ^0, QQ^1, 0);
eqlhsPd = matrix {{-1,0,0},{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{-1/2},{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-1/2,0,1/2,1/2,1/2,1/2,0,-1/2,-1/2,0,-1/2},{-1/2,0,1/2,-1/2,0,1/2,0,1/2,1/2,0,-1/2,-1/2},{1/2,1/2,1/2,0,-1/2,0,1/2,-1/2,0,-1/2,-1/2,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1,1},{1,-1,-1},{-1,1,-1},{1,1,-1},{1,-1,1},{-1,1,1},{1,1,1},{-1,-1,-1},{0,0,1},{0,1,0},{1,0,0},{-1,0,0},{0,-1,0},{0,0,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1/2},{1/2},{1/2},{1/2},{1/2},{1/2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1},{-1,-1,1,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0},{1,0},{0,-1},{0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 10, facets: 7
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1,1,1/2},{0,0,1,1,0,0,1,1,1/2,1},{0,0,0,0,1,1,1,1/2,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0},{1,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{5/2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,1,1,0,0},{1,0,0,1,0,1},{0,0,1,1,1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0},{1,-1,-1},{-1,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{0},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,-1},{1,0,-1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1},{1,-1},{-1,1},{1,1}};
ineqrhsPd = matrix {{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 3, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,0},{0,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0},{-1,0,0}};
ineqrhsPd = matrix {{1},{0}};
eqlhsPd = matrix {{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 3, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,0},{0,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0},{-1,0,0}};
ineqrhsPd = matrix {{1},{0}};
eqlhsPd = matrix {{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1,0,1,1,2/3,1,1,0,0,1/3},{0,1,1,0,0,1,0,1/3,0,1,2/3,1},{0,0,0,1,1,1,0,0,1/3,2/3,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0},{1,-1,-1},{-1,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{2/3},{5/3}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,0,1,1},{1,1/2,0,0,0,1,1/2,1},{0,0,0,1/2,1,1,1,1/2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0},{1/2,-1,-1},{1/2,1,-1},{1/2,-1,1},{1/2,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{0},{1},{1},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 16, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,1/2,1,1,1,1/2,1,1,1/2,1,1,1,1/2},{0,1,0,1,0,1/4,0,3/4,1,1,0,0,1/4,1,3/4,1},{0,0,1,1,0,0,1/4,0,0,1/4,3/4,1,1,3/4,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0},{1/2,-1,-1},{1/2,1,-1},{1/2,-1,1},{1/2,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{1/4},{5/4},{5/4},{9/4}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 12, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2,0,0,0,0},{1/2,-1/2,0,0,1/2,-1/2,0,0,1/2,1/2,-1/2,-1/2},{0,0,1/2,-1/2,0,0,1/2,-1/2,1/2,-1/2,1/2,-1/2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1},{2,0,0},{-2,0,0},{0,2,0},{0,-2,0},{0,0,2},{0,0,-2}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 3, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1/2,1},{0,0},{0,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0},{-1,0,0}};
ineqrhsPd = matrix {{1},{-1/2}};
eqlhsPd = matrix {{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{1,0,0,0},{0,0,1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,0,0},{0,-1,0},{0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 6, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{3/2,-3/2,0,0,0,0},{0,0,3/2,-3/2,0,0},{0,0,0,0,3/2,-3/2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
ineqrhsPd = matrix {{3/2},{3/2},{3/2},{3/2},{3/2},{3/2},{3/2},{3/2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 24, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1/3,-1/3,0,0,-1/3,0,0,1/3,1/3,2/3,2/3,2/3,2/3,1/3,1/3,0,0,0,0,-1/3,-2/3,-2/3,-2/3,-2/3},{-2/3,0,-2/3,-1/3,2/3,1/3,2/3,-2/3,2/3,0,1/3,0,-1/3,0,0,2/3,1/3,-1/3,-2/3,0,1/3,0,0,-1/3},{0,-2/3,-1/3,-2/3,0,-2/3,-1/3,0,0,1/3,0,-1/3,0,2/3,-2/3,1/3,2/3,2/3,1/3,2/3,0,1/3,-1/3,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,1,-1},{-1,1,1},{-1,-1,1},{-1,-1,-1},{1,-1,1},{1,-1,-1},{1,1,1},{1,1,-1},{0,0,-1},{-1,0,0},{0,-1,0},{0,0,1},{0,1,0},{1,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{2/3},{2/3},{2/3},{2/3},{2/3},{2/3}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 10, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-468/2581,-1716/5077,-156/359,-1/2,-156/359,-1716/5077,-468/2581,0,15600/2029},{156/329,1092/2581,1716/5077,468/2513,0,-468/2513,-1716/5077,-1092/2581,-156/329,0},{173/329,1333/2581,2581/5077,1265/2513,1/2,1265/2513,2581/5077,1333/2581,173/329,1873/2029}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,3,3},{-1,5/3,11/6},{-5/3,1,11/6},{-3,-1,3},{-3,1,3},{-5/3,-1,11/6},{-1,-5/3,11/6},{-1,-3,3},{1,-100,100},{1,100,100}};
ineqrhsPd = matrix {{3},{11/6},{11/6},{3},{3},{11/6},{11/6},{3},{100},{100}};
eqlhsPd = matrix {{-34/3,0,658/3}};
eqrhsPd = matrix {{346/3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{0,0,-1},{0,-1,0},{0,1,0},{1,0,0},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0},{1,0,0,1,0},{0,0,1,1,2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1,-1},{0,-1,0},{-1,0,0},{-1,2,1},{2,-1,1},{1,1,-1}};
ineqrhsPd = matrix {{-1},{0},{0},{2},{2},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{1,0,0,0},{0,0,1,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1},{-1,0,0},{0,-1,0},{0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2},{2}};
raysP = matrix {{1,1},{-1,1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0},{-1,1},{-1,-1}};
ineqrhsPd = matrix {{1},{0},{-4}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,2,2},{0,2,0,2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{0},{0},{2},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,1}};
raysP = matrix {{1},{0}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1},{-1,1},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1},{0,1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 12, facets: 7
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,0,0,0,0,0,1,1,1,1},{0,0,0,0,1,1,1,1,0,0,0,0},{0,0,2,2,0,0,2,2,0,0,2,2},{0,2,0,2,0,2,0,2,0,2,0,2}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-1},{-1,0,0,0},{0,0,-1,0},{0,0,0,1},{0,0,1,0},{0,-1,0,0},{1,1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{2},{2},{0},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,1},{2,2},{2,2}};
raysP = matrix {{1,0,0},{0,0,0},{0,1,1},{0,-1,1}};
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,0,0},{-1,1,0,0},{0,-1,0,0},{0,0,-1,1},{0,0,-1,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{-4}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,1}};
raysP = matrix {{1},{0}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1},{-1,1},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2},{2}};
raysP = matrix {{1,1},{-1,1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0},{-1,1},{-1,-1}};
ineqrhsPd = matrix {{1},{0},{-4}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,2,2},{0,2,0,2}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-1},{-1,0},{0,1},{1,0}};
ineqrhsPd = matrix {{0},{0},{2},{2}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1},{0,1},{2,2},{2,2}};
raysP = matrix {{1,0,0},{0,0,0},{0,1,1},{0,-1,1}};
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,1,0,0},{-1,1,0,0},{0,-1,0,0},{0,0,-1,1},{0,0,-1,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{-4}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1},{0,1,0}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1},{-1}};
ineqrhsPd = matrix {{1},{1}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1},{-1}};
ineqrhsPd = matrix {{1},{1}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 9, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1/2,0,1/2,0,1,1/2,1/2,0},{0,1/2,0,-1/2,0,0,-1/2,1/2,0},{1/2,0,0,0,1/2,0,1/2,1/2,1},{1/2,0,0,0,-1/2,0,1/2,-1/2,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,-1,0,0},{0,0,-1,-1},{1,-1,1,-1},{-1,1,0,0},{0,0,-1,1}};
ineqrhsPd = matrix {{1},{0},{0},{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0},{0,0,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,1,-1,1,-1,1,-1},{2,0,0,2,-2,0,-2,0},{0,-2,0,2,0,2,-2,0}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,0,0},{-1,0,0},{1,1,0},{-1,-1,0},{1,0,1},{-1,0,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0},{0,0,1}};
raysP = map(QQ^2, QQ^0, 0);
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1},{-1,0},{0,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0},{0,0,1,0},{0,0,0,2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0},{0,-1,0},{0,0,-1},{2,2,1}};
ineqrhsPd = matrix {{0},{0},{0},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{1,0,0,0},{0,-1,0,0},{0,1,0,0},{0,0,-1,0},{0,0,1,0},{0,0,0,-1},{0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 3, vertices: 10, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,-468/2581,-1716/5077,-156/359,-1/2,-156/359,-1716/5077,-468/2581,0,15600/2029},{156/329,1092/2581,1716/5077,468/2513,0,-468/2513,-1716/5077,-1092/2581,-156/329,0},{173/329,1333/2581,2581/5077,1265/2513,1/2,1265/2513,2581/5077,1333/2581,173/329,1873/2029}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{2029/32900,1,0},{2029/32900,-1,0},{-278/987,-1,0},{-1787/3290,-1,0},{-3103/1974,-1,0},{-936/329,1,0},{-936/329,-1,0},{-3103/1974,1,0},{-1787/3290,1,0},{-278/987,1,0}};
ineqrhsPd = matrix {{156/329},{156/329},{156/329},{858/1645},{286/329},{468/329},{468/329},{286/329},{858/1645},{156/329}};
eqlhsPd = matrix {{-34/3,0,658/3}};
eqrhsPd = matrix {{346/3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 4, vertices: 6, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,3/2},{3,2},{4,5/2},{5,3}};
raysP = matrix {{1,-1,0,0},{-1,1,0,0},{2,-2,1,0},{-2,2,-1,0}};
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,0,0},{-1,-1,0,0},{7,1,-3,0}};
ineqrhsPd = matrix {{5},{-7/2},{5}};
eqlhsPd = matrix {{7/3,7/3,-1,-1}};
eqrhsPd = matrix {{8/3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1}};
raysP = matrix {{1}};
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0},{-1}};
ineqrhsPd = matrix {{1},{-1}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,3/2},{3,2},{4,5/2},{5,3}};
raysP = matrix {{1,-1,0},{3/2,2,0},{2,9/2,1},{3,-7/2,-1}};
linealityP = matrix {{1},{-1},{2},{-2}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,-1,-1},{30/13,2,-15/13,-1},{0,0,0,0},{-2,-2,1,1},{0,-2,0,1},{-7/3,-7/3,1,1}};
ineqrhsPd = matrix {{-2},{41/26},{1},{-1},{-1},{-8/3}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 2, ambientDim: 2, vertices: 3, facets: 3
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0},{0}};
raysP = matrix {{0,1},{1,1}};
linealityP = map(QQ^2, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0},{-1,0},{1,-1}};
ineqrhsPd = matrix {{1},{0},{0}};
eqlhsPd = map(QQ^0, QQ^2, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 1, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1}};
eqrhsPd = matrix {{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1},{1}};
ineqrhsPd = matrix {{1},{1}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 72, facets: 37
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-3/2,-3/2,-3/2,-3/2,-3/2,-3/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2,1/2,1/2,1/2,1/2,1/2,3/2,3/2,3/2,3/2,3/2,3/2,-3/2,-3/2,-3/2,-3/2,-3/2,-3/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2,1/2,1/2,1/2,1/2,1/2,3/2,3/2,3/2,3/2,3/2,3/2,-3/2,-3/2,-3/2,-3/2,-3/2,-3/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2,1/2,1/2,1/2,1/2,1/2,3/2,3/2,3/2,3/2,3/2,3/2},{-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2,-5/2,-3/2,-1/2,1/2,3/2,5/2},{19/18,11/18,7/18,7/18,11/18,19/18,5/6,7/18,1/6,1/6,7/18,5/6,5/6,7/18,1/6,1/6,7/18,5/6,19/18,11/18,7/18,7/18,11/18,19/18,17/18,1/2,5/18,5/18,1/2,17/18,13/18,5/18,1/18,1/18,5/18,13/18,13/18,5/18,1/18,1/18,5/18,13/18,17/18,1/2,5/18,5/18,1/2,17/18,19/18,11/18,7/18,7/18,11/18,19/18,5/6,7/18,1/6,1/6,7/18,5/6,5/6,7/18,1/6,1/6,7/18,5/6,19/18,11/18,7/18,7/18,11/18,19/18}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,2,2,-9},{1,2,-2,-9},{0,0,-1,0},{2,0,4,-18},{2,0,-4,-18},{1,-2,2,-9},{1,-2,-2,-9},{-1,2,2,-9},{-1,2,-2,-9},{-2,0,4,-18},{-2,0,-4,-18},{-1,-2,2,-9},{-1,-2,-2,-9},{-1,0,0,0},{-1,-2,-4,-9},{-2,-4,0,-18},{-1,-2,4,-9},{-1,0,-4,-9},{-2,0,0,-18},{-1,0,4,-9},{-1,2,-4,-9},{-2,4,0,-18},{-1,2,4,-9},{0,-1,0,0},{1,-2,-4,-9},{2,-4,0,-18},{1,-2,4,-9},{1,0,-4,-9},{2,0,0,-18},{1,0,4,-9},{1,2,-4,-9},{2,4,0,-18},{0,0,1,0},{0,0,0,1},{1,2,4,-9},{0,1,0,0},{1,0,0,0}};
ineqrhsPd = matrix {{3/2},{3/2},{5/2},{1},{1},{3/2},{3/2},{3/2},{3/2},{1},{1},{3/2},{3/2},{1},{9/2},{1},{9/2},{7/2},{-1},{7/2},{9/2},{1},{9/2},{3/2},{9/2},{1},{9/2},{7/2},{-1},{7/2},{9/2},{1},{5/2},{19/18},{9/2},{3/2},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 4, facets: 4
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,0,1},{0,0,2,2},{0,3,3,3}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,1},{-1,0,0},{2,-1,0},{0,3/2,-1}};
ineqrhsPd = matrix {{3},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 6, vertices: 7, facets: 7
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,0,0,-1},{0,1,0,0,0,0,-1},{0,0,1,0,0,0,-1},{0,0,0,1,0,0,-1},{0,0,0,0,1,0,-1},{0,0,0,0,0,1,-1}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1,1,1},{1,1,1,1,1,-6},{1,-6,1,1,1,1},{1,1,-6,1,1,1},{1,1,1,-6,1,1},{1,1,1,1,-6,1},{-6,1,1,1,1,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^6, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,5/8,0,0,0},{0,0,5/8,0,0},{0,0,0,5/8,0},{0,0,0,0,5/8}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{8/5,8/5,8/5,8/5},{-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 0, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = map(QQ^0, QQ^1, 0);
raysP = map(QQ^0, QQ^0, 0);
linealityP = map(QQ^0, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = map(QQ^1, QQ^0, 0);
ineqrhsPd = matrix {{1}};
eqlhsPd = map(QQ^0, QQ^0, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 8, facets: 6
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1},{0,0,1,1,0,0,1,1},{0,0,0,0,1,1,1,1}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{0,-1,0},{-1,0,0},{0,0,1},{1,0,0},{0,1,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 9, facets: 9
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0,1,0,1,1/2},{0,0,1,1,0,0,1,1,1/2},{0,0,0,0,1,1,1,1,3/2}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-1},{0,-1,0},{-1,0,0},{-1,0,1},{0,-1,1},{1,0,0},{0,1,0},{1,0,1},{0,1,1}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1},{1},{1},{2},{2}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 9, ambientDim: 16, vertices: 24, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0},{0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1},{0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0},{0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,1,0,0},{0,1,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0},{1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0},{0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1},{0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0},{0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0},{0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0},{1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0},{0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0},{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
raysP = map(QQ^16, QQ^0, 0);
linealityP = map(QQ^16, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0},{-1,-1,-1,0,-1,-1,-1,0,-1,-1,-1,0,0,0,0,0},{-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0},{0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0},{0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0},{0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0},{0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0},{1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0},{1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0}};
ineqrhsPd = matrix {{0},{-2},{0},{0},{1},{1},{0},{0},{1},{1},{0},{0},{0},{0},{1},{1}};
eqlhsPd = matrix {{-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,-1,-1,-1,-1,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,0},{-3/4,1/4,1/4,1/4,-3/4,1/4,1/4,1/4,-3/4,1/4,1/4,1/4,-1,0,0,0},{1/13,-9/13,4/13,4/13,1/13,-9/13,4/13,4/13,1/13,-9/13,4/13,4/13,-3/13,-1,0,0},{1/10,1/10,-3/5,2/5,1/10,1/10,-3/5,2/5,1/10,1/10,-3/5,2/5,-3/10,-3/10,-1,0},{1/7,1/7,1/7,-3/7,1/7,1/7,1/7,-3/7,1/7,1/7,1/7,-3/7,-3/7,-3/7,-3/7,-1}};
eqrhsPd = matrix {{-1},{-1},{-1},{-1/4},{-4/13},{-2/5},{-4/7}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 9, ambientDim: 16, vertices: 12, facets: 64
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,0,0,0,1,0,1,0,0,0,0,0},{0,0,1,0,0,1,0,0,0,0,1,0},{0,1,0,0,0,0,0,0,1,0,0,1},{0,0,0,1,0,0,0,1,0,1,0,0},{0,1,0,0,0,1,0,1,0,0,0,0},{1,0,0,1,0,0,0,0,0,0,0,1},{0,0,1,0,0,0,1,0,0,1,0,0},{0,0,0,0,1,0,0,0,1,0,1,0},{0,0,1,1,0,0,0,0,1,0,0,0},{0,1,0,0,1,0,0,0,0,1,0,0},{1,0,0,0,0,0,0,1,0,0,1,0},{0,0,0,0,0,1,1,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,1,1,1},{0,0,0,0,0,0,1,1,1,0,0,0},{0,0,0,1,1,1,0,0,0,0,0,0},{1,1,1,0,0,0,0,0,0,0,0,0}};
raysP = map(QQ^16, QQ^0, 0);
linealityP = map(QQ^16, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1,-1,0,-1,-1,-1,0,-1,-1,-1,0,0,0,0,0},{-1,-1,-2,0,-1,-2,-1,0,0,-1,-1,0,0,0,0,0},{-1,-1,0,0,-1,-2,-1,0,-2,-1,-1,0,0,0,0,0},{0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0},{-1,0,-1,0,-2,-1,-1,0,-1,-1,-2,0,0,0,0,0},{0,1,1,0,-1,-1,0,0,-1,0,-1,0,0,0,0,0},{-1,-2,-1,0,0,-1,-1,0,-1,-1,-2,0,0,0,0,0},{0,-1,-1,0,1,-1,0,0,1,0,-1,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0},{-2,-1,-1,0,-1,-1,0,0,-1,-2,-1,0,0,0,0,0},{-1,0,1,0,0,-1,1,0,-1,-1,0,0,0,0,0,0},{-1,-1,0,0,1,0,1,0,0,-1,-1,0,0,0,0,0},{-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{-1,1,0,0,-1,0,-1,0,0,1,-1,0,0,0,0,0},{-1,0,-1,0,0,-1,-1,0,1,1,0,0,0,0,0,0},{-2,-1,-1,0,-1,-1,-2,0,-1,0,-1,0,0,0,0,0},{0,1,1,0,1,1,0,0,1,0,1,0,0,0,0,0},{1,2,1,0,0,1,-1,0,1,1,0,0,0,0,0,0},{1,1,0,0,1,0,-1,0,2,1,1,0,0,0,0,0},{0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0},{-1,1,0,0,1,2,1,0,0,1,1,0,0,0,0,0},{-1,0,1,0,0,1,1,0,1,1,2,0,0,0,0,0},{0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0},{0,1,-1,0,1,1,0,0,1,2,1,0,0,0,0,0},{-1,0,-1,0,0,1,-1,0,-1,1,0,0,0,0,0,0},{0,1,1,0,-1,1,0,0,1,2,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0},{-1,-1,0,0,-1,0,-1,0,0,1,1,0,0,0,0,0},{1,0,1,0,2,1,1,0,1,-1,0,0,0,0,0,0},{1,1,2,0,1,0,1,0,0,-1,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0},{2,1,1,0,1,1,0,0,1,0,-1,0,0,0,0,0},{1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0},{1,-1,0,0,1,0,-1,0,0,-1,-1,0,0,0,0,0},{1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},{1,1,0,0,-1,0,-1,0,0,-1,-1,0,0,0,0,0},{2,1,1,0,1,-1,0,0,1,0,1,0,0,0,0,0},{1,0,1,0,0,-1,-1,0,-1,-1,0,0,0,0,0,0},{1,0,-1,0,0,-1,-1,0,1,-1,0,0,0,0,0,0},{0,-1,-1,0,-1,-1,-2,0,-1,-2,-1,0,0,0,0,0},{0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0},{0,1,1,0,1,1,2,0,-1,0,1,0,0,0,0,0},{-1,0,-1,0,0,1,1,0,-1,-1,0,0,0,0,0,0},{1,1,0,0,1,2,1,0,0,1,-1,0,0,0,0,0},{1,0,-1,0,2,1,1,0,1,1,0,0,0,0,0,0},{0,-1,-1,0,1,1,0,0,-1,0,-1,0,0,0,0,0},{1,2,1,0,0,1,1,0,-1,1,0,0,0,0,0,0},{0,1,-1,0,-1,1,0,0,-1,0,-1,0,0,0,0,0},{1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0},{0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0},{0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0},{-1,-1,-2,0,-1,0,-1,0,-2,-1,-1,0,0,0,0,0},{0,-1,1,0,1,1,2,0,1,0,1,0,0,0,0,0},{0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0},{-1,-1,0,0,-1,0,1,0,0,-1,1,0,0,0,0,0},{1,0,1,0,0,1,1,0,1,1,0,0,0,0,0,0},{1,-1,0,0,1,0,1,0,2,1,1,0,0,0,0,0},{0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{1,1,2,0,-1,0,1,0,0,1,1,0,0,0,0,0},{0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0},{1,0,1,0,0,-1,1,0,1,1,2,0,0,0,0,0},{0,-1,1,0,-1,-1,0,0,-1,0,1,0,0,0,0,0},{0,-1,-1,0,-1,-1,0,0,1,0,1,0,0,0,0,0},{-1,-2,-1,0,-2,-1,-1,0,-1,-1,0,0,0,0,0,0}};
ineqrhsPd = matrix {{-2},{-2},{-2},{0},{-2},{0},{-2},{0},{0},{-2},{0},{0},{0},{0},{0},{-2},{2},{2},{2},{0},{2},{2},{1},{2},{0},{2},{1},{0},{2},{2},{0},{2},{1},{0},{1},{0},{2},{0},{0},{-2},{1},{2},{0},{2},{2},{0},{2},{0},{2},{0},{0},{-2},{2},{1},{0},{2},{2},{0},{2},{0},{2},{0},{0},{-2}};
eqlhsPd = matrix {{-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,-1,-1,-1,-1,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,0},{-3/4,1/4,1/4,1/4,-3/4,1/4,1/4,1/4,-3/4,1/4,1/4,1/4,-1,0,0,0},{1/13,-9/13,4/13,4/13,1/13,-9/13,4/13,4/13,1/13,-9/13,4/13,4/13,-3/13,-1,0,0},{1/10,1/10,-3/5,2/5,1/10,1/10,-3/5,2/5,1/10,1/10,-3/5,2/5,-3/10,-3/10,-1,0},{1/7,1/7,1/7,-3/7,1/7,1/7,1/7,-3/7,1/7,1/7,1/7,-3/7,-3/7,-3/7,-3/7,-1}};
eqrhsPd = matrix {{-1},{-1},{-1},{-1/4},{-4/13},{-2/5},{-4/7}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 4, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{2,4},{4,2},{4,2},{2,4}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,-2,0,0},{0,0,0,-2}};
ineqrhsPd = matrix {{-4},{-4}};
eqlhsPd = matrix {{-1,-1,-1,-1},{1,-1,1,-1},{1,1,-1,-1}};
eqrhsPd = matrix {{-12},{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 8, vertices: 22, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{4,4/3,20/3,4,4,4/3,4/3,4/3,20/3,4/3,4,4,4/3,4,4,20/3,4,20/3,4,4/3,20/3,28/3},{4/3,4,4,4,4/3,20/3,4,20/3,4/3,4,20/3,20/3,28/3,20/3,4,4,4/3,4/3,4/3,4,4,4/3},{4/3,4,4/3,4/3,4,4,20/3,20/3,4/3,4,4/3,4,4,4/3,20/3,4,20/3,4,20/3,28/3,4,4/3},{28/3,20/3,4,20/3,20/3,4,4,4/3,20/3,20/3,4,4/3,4/3,4,4/3,4/3,4,4,4,4/3,4/3,4},{20/3,28/3,4,20/3,20/3,4,4,20/3,4/3,20/3,4,4/3,4,4/3,4/3,4/3,4,4,4/3,4,4,4/3},{4,4/3,4/3,4/3,4,4,20/3,4/3,20/3,4,4/3,4,4/3,4,20/3,4,20/3,4,28/3,20/3,4/3,4},{4,4/3,4,4,4/3,20/3,4,4/3,20/3,4,20/3,20/3,20/3,28/3,4,4,4/3,4/3,4,4/3,4/3,4},{4/3,4,20/3,4,4,4/3,4/3,20/3,4/3,4/3,4,4,4,4/3,4,20/3,4,20/3,4/3,4,28/3,20/3}};
raysP = map(QQ^8, QQ^0, 0);
linealityP = map(QQ^8, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,-2,-2,0,0,0,0},{0,-2,0,-2,0,0,0,0},{0,-2,0,0,0,-2,0,0},{0,0,-2,0,-2,0,-4,-2},{0,-2,0,0,-2,-4,0,-2},{0,-2,-2,-4,0,0,0,-2},{0,0,0,-2,0,-2,-2,-4},{0,0,0,-2,0,0,0,-2},{-2,0,0,0,-4,-2,-2,0},{0,0,0,0,0,-2,0,-2},{-2,0,-4,-2,0,0,-2,0},{0,0,0,0,0,0,-2,-2},{-2,-4,0,-2,0,-2,0,0},{-4,-2,-2,0,-2,0,0,0}};
ineqrhsPd = matrix {{-32/3},{-32/3},{-32/3},{-104/3},{-104/3},{-104/3},{-104/3},{-32/3},{-104/3},{-32/3},{-104/3},{-32/3},{-104/3},{-104/3}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1,-1},{1,-1,1,-1,1,-1,1,-1},{1,1,-1,-1,1,1,-1,-1},{1,1,1,1,-1,-1,-1,-1}};
eqrhsPd = matrix {{-32},{0},{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 8, vertices: 132, facets: 20
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,9553577508788659/4503599627370496,38824863540916338353528385288599/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,9553577508788659/4503599627370496,1592262918131443/2251799813685248,1592262918131443/2251799813685248,44765443775022755566307833752773/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,24483034171371086679852902234519/10141204801825835211973625643008,53166692910461591153103775185303/20282409603651670423947251286016,9553577508788659/4503599627370496,38824863540916338353528385288599/20282409603651670423947251286016,1592262918131443/2251799813685248,9553577508788659/4503599627370496,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,24483034171371085553952995391895/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,24483034171371085553952995391895/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,34483405944844081/18014398509481984,48966068342742172439602104088983/20282409603651670423947251286016,47221509289895627/18014398509481984,48966068342742172439602104088983/20282409603651670423947251286016,34483405944844081/18014398509481984,47221509289895627/18014398509481984,34483405944844081/18014398509481984,34483405944844081/18014398509481984,1592262918131443/1125899906842624,1592262918131443/1125899906842624,1592262918131443/1125899906842624,1592262918131443/1125899906842624,12738103345051545/4503599627370496,53166692910461588901303961500055/20282409603651670423947251286016,53166692910461588901303961500055/20282409603651670423947251286016,24483034171371085553952995391895/10141204801825835211973625643008,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,21745302599792537/9007199254740992},{3730904090310553/18014398509481984,1592262918131443/2251799813685248,3730904090310553/18014398509481984,1592262918131443/2251799813685248,21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,9553577508788659/4503599627370496,38824863540916338353528385288599/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,53166692910461591153103775185303/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,34483405944844083/18014398509481984,48966068342742174691401917774231/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,12738103345051545/4503599627370496,47221509289895627/18014398509481984,9553577508788659/4503599627370496,38824863540916338353528385288599/20282409603651670423947251286016,48966068342742172439602104088983/20282409603651670423947251286016,53166692910461588901303961500055/20282409603651670423947251286016,19107155017577317/9007199254740992,47221509289895627/18014398509481984,21745302599792537/9007199254740992,12738103345051545/9007199254740992,34483405944844081/18014398509481984,34483405944844083/18014398509481984,12738103345051545/9007199254740992,34483405944844081/18014398509481984,19107155017577317/9007199254740992,21745302599792537/9007199254740992,1592262918131443/2251799813685248,34483405944844081/18014398509481984,9553577508788659/4503599627370496,47221509289895627/18014398509481984,19107155017577317/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,34483405944844081/18014398509481984,19107155017577317/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,34483405944844083/18014398509481984,21745302599792537/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,39759701109274521/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,34483405944844083/18014398509481984,6369051672525773/9007199254740992,1592262918131443/2251799813685248,9553577508788659/4503599627370496,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,6369051672525773/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,12738103345051545/9007199254740992,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984},{9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,9553577508788659/4503599627370496,9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,38824863540916340605328198973847/20282409603651670423947251286016,6369051672525773/4503599627370496,24483034171371087805752809077143/10141204801825835211973625643008,34483405944844083/18014398509481984,48966068342742174691401917774231/20282409603651670423947251286016,6369051672525773/4503599627370496,53166692910461591153103775185303/20282409603651670423947251286016,34483405944844083/18014398509481984,48966068342742174691401917774231/20282409603651670423947251286016,6369051672525773/4503599627370496,53166692910461591153103775185303/20282409603651670423947251286016,34483405944844083/18014398509481984,47221509289895627/18014398509481984,34483405944844083/18014398509481984,6369051672525773/4503599627370496,21745302599792537/9007199254740992,47221509289895627/18014398509481984,12738103345051545/4503599627370496,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,38824863540916338353528385288599/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,19107155017577317/9007199254740992,19107155017577317/9007199254740992,19107155017577317/9007199254740992,53166692910461588901303961500055/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,19107155017577317/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,44765443775022757818107647438021/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,12738103345051545/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,12738103345051545/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992},{21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,34483405944844083/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/9007199254740992,34483405944844083/18014398509481984,6369051672525773/9007199254740992,39759701109274521/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,19107155017577317/9007199254740992,34483405944844081/18014398509481984,19107155017577317/9007199254740992,34483405944844081/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,21745302599792537/9007199254740992,9553577508788659/4503599627370496,47221509289895627/18014398509481984,1592262918131443/2251799813685248,9553577508788659/4503599627370496,34483405944844081/18014398509481984,12738103345051545/9007199254740992,19107155017577317/9007199254740992,47221509289895627/18014398509481984,34483405944844083/18014398509481984,34483405944844081/18014398509481984,21745302599792537/9007199254740992,12738103345051545/9007199254740992,19107155017577317/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,48966068342742174691401917774231/20282409603651670423947251286016,34483405944844083/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,9553577508788659/4503599627370496,47221509289895627/18014398509481984,48966068342742172439602104088983/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,12738103345051545/4503599627370496,53166692910461591153103775185303/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,53166692910461588901303961500055/20282409603651670423947251286016,9553577508788659/4503599627370496,24483034171371086679852902234519/10141204801825835211973625643008,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,12738103345051545/9007199254740992,6369051672525773/9007199254740992,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,3730904090310553/18014398509481984},{1592262918131443/2251799813685248,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,24483034171371085553952995391895/20282409603651670423947251286016,9553577508788659/4503599627370496,9553577508788659/4503599627370496,53166692910461591153103775185303/20282409603651670423947251286016,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,44765443775022755566307833752773/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,38824863540916338353528385288599/20282409603651670423947251286016,1592262918131443/2251799813685248,9553577508788659/4503599627370496,9553577508788659/4503599627370496,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,24483034171371085553952995391895/20282409603651670423947251286016,24483034171371085553952995391895/20282409603651670423947251286016,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,24483034171371085553952995391895/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,4200624567719415541598157031017/20282409603651670423947251286016,48966068342742172439602104088983/20282409603651670423947251286016,34483405944844081/18014398509481984,34483405944844081/18014398509481984,34483405944844081/18014398509481984,48966068342742172439602104088983/20282409603651670423947251286016,34483405944844081/18014398509481984,47221509289895627/18014398509481984,47221509289895627/18014398509481984,12738103345051545/4503599627370496,53166692910461588901303961500055/20282409603651670423947251286016,53166692910461588901303961500055/20282409603651670423947251286016,24483034171371085553952995391895/10141204801825835211973625643008,1592262918131443/1125899906842624,1592262918131443/1125899906842624,1592262918131443/1125899906842624,1592262918131443/1125899906842624,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,19107155017577317/9007199254740992,38824863540916336101728571603351/20282409603651670423947251286016,38824863540916336101728571603351/20282409603651670423947251286016,19107155017577317/9007199254740992,19107155017577317/9007199254740992,38824863540916336101728571603351/20282409603651670423947251286016,21745302599792537/9007199254740992},{3730904090310553/18014398509481984,3730904090310553/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,21745302599792537/18014398509481984,9553577508788659/4503599627370496,53166692910461588901303961500055/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,38824863540916338353528385288599/20282409603651670423947251286016,53166692910461591153103775185303/20282409603651670423947251286016,12738103345051545/4503599627370496,47221509289895627/18014398509481984,9553577508788659/4503599627370496,38824863540916338353528385288599/20282409603651670423947251286016,48966068342742172439602104088983/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,34483405944844083/18014398509481984,48966068342742174691401917774231/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,19107155017577317/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,21745302599792537/9007199254740992,34483405944844081/18014398509481984,34483405944844083/18014398509481984,47221509289895627/18014398509481984,19107155017577317/9007199254740992,34483405944844081/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,34483405944844083/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,39759701109274521/18014398509481984,21745302599792537/9007199254740992,6369051672525773/9007199254740992,34483405944844083/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,19107155017577317/9007199254740992,19107155017577317/9007199254740992,34483405944844081/18014398509481984,34483405944844081/18014398509481984,21745302599792537/9007199254740992,47221509289895627/18014398509481984,9553577508788659/4503599627370496,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,9553577508788659/4503599627370496,1592262918131443/2251799813685248,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984},{9553577508788659/4503599627370496,9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,38824863540916340605328198973847/20282409603651670423947251286016,9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,9553577508788659/4503599627370496,38824863540916340605328198973847/20282409603651670423947251286016,24483034171371087805752809077143/10141204801825835211973625643008,6369051672525773/4503599627370496,48966068342742174691401917774231/20282409603651670423947251286016,34483405944844083/18014398509481984,53166692910461591153103775185303/20282409603651670423947251286016,6369051672525773/4503599627370496,48966068342742174691401917774231/20282409603651670423947251286016,34483405944844083/18014398509481984,53166692910461591153103775185303/20282409603651670423947251286016,6369051672525773/4503599627370496,47221509289895627/18014398509481984,34483405944844083/18014398509481984,47221509289895627/18014398509481984,12738103345051545/4503599627370496,21745302599792537/9007199254740992,34483405944844083/18014398509481984,6369051672525773/4503599627370496,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,6369051672525773/9007199254740992,6369051672525773/9007199254740992,24483034171371086679852902234519/10141204801825835211973625643008,38824863540916338353528385288599/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,19107155017577317/9007199254740992,53166692910461588901303961500055/20282409603651670423947251286016,19107155017577317/9007199254740992,19107155017577317/9007199254740992,6369051672525773/9007199254740992,19107155017577317/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,24483034171371086679852902234519/10141204801825835211973625643008,38824863540916338353528385288599/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,44765443775022757818107647438021/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,24483034171371087805752809077143/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,6369051672525773/9007199254740992,6369051672525773/9007199254740992,4200624567719417793397970716265/20282409603651670423947251286016,12738103345051545/9007199254740992},{21745302599792537/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,12738103345051545/9007199254740992,3730904090310553/18014398509481984,1592262918131443/2251799813685248,3730904090310553/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,1592262918131443/2251799813685248,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,1592262918131443/2251799813685248,6369051672525773/9007199254740992,12738103345051545/9007199254740992,1592262918131443/2251799813685248,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,1592262918131443/2251799813685248,12738103345051545/9007199254740992,12738103345051545/9007199254740992,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,1592262918131443/2251799813685248,21745302599792537/9007199254740992,34483405944844081/18014398509481984,1592262918131443/2251799813685248,1592262918131443/2251799813685248,34483405944844081/18014398509481984,9553577508788659/4503599627370496,47221509289895627/18014398509481984,19107155017577317/9007199254740992,19107155017577317/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,21745302599792537/9007199254740992,34483405944844083/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,39759701109274521/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,34483405944844083/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,9553577508788659/4503599627370496,1592262918131443/2251799813685248,19107155017577317/9007199254740992,47221509289895627/18014398509481984,34483405944844081/18014398509481984,12738103345051545/9007199254740992,34483405944844083/18014398509481984,34483405944844081/18014398509481984,12738103345051545/9007199254740992,21745302599792537/9007199254740992,19107155017577317/9007199254740992,53166692910461588901303961500055/20282409603651670423947251286016,48966068342742172439602104088983/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,9553577508788659/4503599627370496,47221509289895627/18014398509481984,12738103345051545/4503599627370496,34483405944844083/18014398509481984,12738103345051545/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,48966068342742174691401917774231/20282409603651670423947251286016,12738103345051545/9007199254740992,38824863540916338353528385288599/20282409603651670423947251286016,53166692910461591153103775185303/20282409603651670423947251286016,38824863540916338353528385288599/20282409603651670423947251286016,9553577508788659/4503599627370496,24483034171371086679852902234519/10141204801825835211973625643008,12738103345051545/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,3730904090310553/18014398509481984,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,12738103345051545/9007199254740992,3730904090310553/18014398509481984,3730904090310553/18014398509481984,6369051672525773/9007199254740992,6369051672525773/9007199254740992,6369051672525773/9007199254740992,3730904090310553/18014398509481984,6369051672525773/9007199254740992,3730904090310553/18014398509481984,12738103345051545/9007199254740992,21745302599792537/18014398509481984,21745302599792537/18014398509481984,12738103345051545/9007199254740992,3730904090310553/18014398509481984}};
raysP = map(QQ^8, QQ^0, 0);
linealityP = map(QQ^8, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0,-9007199254740993/9007199254740992,-2,-21745302599792537/9007199254740992,-2,-9007199254740991/9007199254740992},{0,0,0,0,-2,-15376250927266765/4503599627370496,-21745302599792537/6369051672525772,-2},{0,-1592262918131443/2251799813685248,-1,-1592262918131443/2251799813685248,0,0,0,0},{0,-2,-21745302599792537/6369051672525772,-15376250927266765/4503599627370496,-2,0,0,0},{0,-9007199254740991/9007199254740992,-2,-21745302599792537/9007199254740992,-2,-9007199254740993/9007199254740992,0,0},{-1,0,0,0,-1,-2,-3844062731816691/1592262918131443,-2},{-2,0,0,0,0,-2,-21745302599792537/6369051672525772,-15376250927266765/4503599627370496},{-1,0,0,0,0,0,-6369051672525773/6369051672525772,-6369051672525773/4503599627370496},{0,0,-2,-3844062731816691/1125899906842624,-21745302599792537/6369051672525773,-2,0,0},{0,0,-1,-2,-15376250927266765/6369051672525773,-2,-1,0},{-2,-9007199254740993/9007199254740992,0,0,0,-9007199254740991/9007199254740992,-2,-21745302599792537/9007199254740992},{-21745302599792537/6369051672525773,-2,0,0,0,0,-2,-3844062731816691/1125899906842624},{-1,-6369051672525773/9007199254740992,0,0,0,0,0,-6369051672525773/9007199254740992},{0,0,0,-2,-21745302599792537/6369051672525773,-3844062731816691/1125899906842624,-2,0},{-15376250927266765/6369051672525773,-2,-1,0,0,0,-1,-2},{-21745302599792537/6369051672525773,-3844062731816691/1125899906842624,-2,0,0,0,0,-2},{-1,-6369051672525773/4503599627370496,-6369051672525773/6369051672525772,0,0,0,0,0},{-2,-21745302599792537/9007199254740992,-2,-9007199254740991/9007199254740992,0,0,0,-9007199254740993/9007199254740992},{-2,-15376250927266765/4503599627370496,-21745302599792537/6369051672525772,-2,0,0,0,0},{-1,-2,-3844062731816691/1592262918131443,-2,-1,0,0,0}};
ineqrhsPd = matrix {{-1390735264807440175688666356357489/162259276829213363391578010288128},{-1263141805282976595730656393454296059875634256959/129179714808990451881365335705143949221903204352},{-34624238973196921219667310126139/20282409603651670423947251286016},{-1263141805282976595730656393454296059875634256959/129179714808990451881365335705143949221903204352},{-1390735264807440175688666356357489/162259276829213363391578010288128},{-276802023886324979020895840658481896763905387605/32294928702247612970341333926285987305475801088},{-1263141805282976595730656393454296059875634256959/129179714808990451881365335705143949221903204352},{-311867419475353320469282520724398592111640124479/129179714808990451881365335705143949221903204352},{-1263141805282976718866540485255142372683213511743/129179714808990472163774939356814373169154490368},{-1107208095545300078044330995351112253391586093803/129179714808990472163774939356814373169154490368},{-1390735264807440175688666356357489/162259276829213363391578010288128},{-1263141805282976718866540485255142372683213511743/129179714808990472163774939356814373169154490368},{-138496955892787688609573330815109/81129638414606681695789005144064},{-1263141805282976718866540485255142372683213511743/129179714808990472163774939356814373169154490368},{-1107208095545300078044330995351112253391586093803/129179714808990472163774939356814373169154490368},{-1263141805282976718866540485255142372683213511743/129179714808990472163774939356814373169154490368},{-311867419475353320469282520724398592111640124479/129179714808990451881365335705143949221903204352},{-1390735264807440175688666356357489/162259276829213363391578010288128},{-1263141805282976595730656393454296059875634256959/129179714808990451881365335705143949221903204352},{-276802023886324979020895840658481896763905387605/32294928702247612970341333926285987305475801088}};
eqlhsPd = matrix {{-1,-1,-1,-1,-1,-1,-1,-1},{-1,-6369051672525773/9007199254740992,0,6369051672525773/9007199254740992,1,6369051672525773/9007199254740992,0,-6369051672525773/9007199254740992},{0,-1592262918131443/2251799813685248,-1,-1592262918131443/2251799813685248,0,1592262918131443/2251799813685248,1,1592262918131443/2251799813685248}};
eqrhsPd = matrix {{-38214310035154635/4503599627370496},{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 32, facets: 10
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0},{1,0,0,0,0},{0,-1,0,0,0},{0,1,0,0,0},{0,0,-1,0,0},{0,0,1,0,0},{0,0,0,-1,0},{0,0,0,1,0},{0,0,0,0,-1},{0,0,0,0,1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^5, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 7, ambientDim: 7, vertices: 8, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,1}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1,1,1,1},{-1,0,0,0,0,0,0},{0,-1,0,0,0,0,0},{0,0,-1,0,0,0,0},{0,0,0,-1,0,0,0},{0,0,0,0,-1,0,0},{0,0,0,0,0,-1,0},{0,0,0,0,0,0,-1}};
ineqrhsPd = matrix {{1},{0},{0},{0},{0},{0},{0},{0}};
eqlhsPd = map(QQ^0, QQ^7, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 1, ambientDim: 1, vertices: 2, facets: 2
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,2}};
raysP = map(QQ^1, QQ^0, 0);
linealityP = map(QQ^1, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1},{1}};
ineqrhsPd = matrix {{0},{2}};
eqlhsPd = map(QQ^0, QQ^1, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 5, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,1,0,1,0},{0,0,1,1,0},{0,0,0,0,1},{0,0,0,0,0},{0,0,0,0,0}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0},{0,0,-1,0,0},{0,-1,0,0,0},{0,1,1,0,0},{1,0,1,0,0}};
ineqrhsPd = matrix {{0},{0},{0},{1},{1}};
eqlhsPd = matrix {{0,0,0,-1,0},{0,0,0,0,-1}};
eqrhsPd = matrix {{0},{0}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0,0,0},{0,0,1,-1,0,0,0,0},{0,0,0,0,1,-1,0,0},{0,0,0,0,0,0,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-1,-1,-1},{-1,1,-1,-1},{-1,-1,-1,1},{1,1,1,1},{1,-1,-1,1},{-1,1,1,-1},{1,-1,-1,-1},{-1,1,1,1},{-1,-1,1,1},{-1,1,-1,1},{1,-1,1,-1},{1,-1,1,1},{1,1,-1,1},{1,1,-1,-1},{-1,-1,1,-1},{1,1,1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 8, facets: 16
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1,-1,0,0,0,0,0,0},{0,0,1,-1,0,0,0,0},{0,0,0,0,1,-1,0,0},{0,0,0,0,0,0,1,-1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1,1,1,1},{-1,1,1,1},{1,-1,1,1},{-1,-1,1,1},{1,1,-1,1},{-1,1,-1,1},{1,-1,-1,1},{-1,-1,-1,1},{1,1,1,-1},{-1,1,1,-1},{1,-1,1,-1},{-1,-1,1,-1},{1,1,-1,-1},{-1,1,-1,-1},{1,-1,-1,-1},{-1,-1,-1,-1}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 6, ambientDim: 6, vertices: 64, facets: 12
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^6, QQ^0, 0);
linealityP = map(QQ^6, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0,0,0},{0,0,-1,0,0,0},{0,-1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,0,-1,0},{0,0,0,-1,0,0},{0,0,0,1,0,0},{0,1,0,0,0,0},{0,0,0,0,0,-1},{1,0,0,0,0,0},{0,0,0,0,0,1},{0,0,0,0,1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^6, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 3, ambientDim: 3, vertices: 5, facets: 5
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{55},{1250},{2800}};
raysP = matrix {{-1,0,1,0},{0,-1,0,1},{-80,0,140,5000}};
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{0,0,0},{140,0,-1},{80,0,-1},{140,5000,-1},{80,5000,-1}};
ineqrhsPd = matrix {{1},{4900},{1600},{6254900},{6251600}};
eqlhsPd = map(QQ^0, QQ^3, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 7, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},{-5/2,-1,-1,1/2,-2,-1/2,-1/2,1,-1,1/2,1/2,2,-1/2,1,1,5/2},{-5/4,1/4,1/4,7/4,-5/4,1/4,1/4,7/4,-1/4,5/4,5/4,11/4,-1/4,5/4,5/4,11/4},{-1,-1,1/2,1/2,-1/2,-1/2,1,1,-1/2,-1/2,1,1,0,0,3/2,3/2}};
raysP = map(QQ^7, QQ^0, 0);
linealityP = map(QQ^7, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-4121/1206,64/67,-1/804,913/804,-200/201,-257/201,1},{-3/1732,324/433,-23169/3464,1537/3464,-675/866,677/866,-1},{192/347,-1732/1041,108/347,236/347,-85/347,-171/347,-1},{4121/1206,-64/67,1/804,-913/804,200/201,257/201,-1},{-2739/124,-708/31,-1537/248,19521/248,1475/62,351/62,1},{2739/124,708/31,1537/248,-19521/248,-1475/62,-351/62,-1},{-192/347,1732/1041,-108/347,-236/347,85/347,171/347,1},{3/1732,-324/433,23169/3464,-1537/3464,675/866,-677/866,1}};
ineqrhsPd = matrix {{10645/2412},{25875/3464},{9775/4164},{14065/2412},{25825/248},{23595/248},{14935/4164},{23545/3464}};
eqlhsPd = matrix {{3/4,3/4,1/4,3/4,-1,0,0},{15/44,15/44,-3/22,1/11,6/11,-1,0},{-603/2192,1041/2192,433/2192,31/2192,115/548,43/274,-1}};
eqrhsPd = matrix {{0},{-3/4},{-145/1096}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 12, facets: 14
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{0,0,1,1,1,1,0,0,-1,-1,-1,-1},{-1,1,-1,1,0,0,1,-1,0,1,0,-1},{-1,-1,0,0,1,-1,1,1,1,0,-1,0},{0,0,0,0,0,0,0,0,0,0,0,0}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = matrix {{-1},{-1},{-1},{1}};
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{1/2,1/2,1/2,3/2},{0,0,1,1},{0,1,0,1},{-1/2,1/2,1/2,1/2},{1,0,0,1},{1/2,-1/2,1/2,1/2},{-1,0,0,-1},{-1/2,-1/2,1/2,-1/2},{1/2,1/2,-1/2,1/2},{0,0,-1,-1},{-1/2,-1/2,-1/2,-3/2},{1/2,-1/2,-1/2,-1/2},{0,-1,0,-1},{-1/2,1/2,-1/2,-1/2}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 4, ambientDim: 4, vertices: 16, facets: 8
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1},{-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1}};
raysP = map(QQ^4, QQ^0, 0);
linealityP = map(QQ^4, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,0,0,0},{0,0,-1,0},{0,-1,0,0},{1,0,0,0},{0,0,0,1},{0,0,0,-1},{0,1,0,0},{0,0,1,0}};
ineqrhsPd = matrix {{1},{1},{1},{1},{1},{1},{1},{1}};
eqlhsPd = map(QQ^0, QQ^4, 0);
eqrhsPd = map(QQ^0, QQ^1, 0);
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 0, ambientDim: 3, vertices: 1, facets: 1
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{1},{2},{3}};
raysP = map(QQ^3, QQ^0, 0);
linealityP = map(QQ^3, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqlhsPd = matrix {{-1,-2,-3}};
ineqrhsPd = matrix {{1}};
eqlhsPd = matrix {{-1,0,0},{0,-1,0},{0,0,-1}};
eqrhsPd = matrix {{-1},{-2},{-3}};
Pd = polyhedronFromHData(ineqlhsPd, ineqrhsPd, eqlhsPd, eqrhsPd);
assert(Pd == P)
assert(isEmpty Pd === isEmpty P)
assert(isCompact Pd === isCompact P)
assert(isLatticePolytope Pd === isLatticePolytope P)
assert(isNormal Pd === isNormal P)
assert(numColumns vertices Pd == numColumns vertices P)
assert(numColumns rays Pd == numColumns rays P)
assert(numColumns linealitySpace Pd == numColumns linealitySpace P)
facetsP = facets P;
facetsPd = facets Pd;
assert(numRows (facetsPd#0) == numRows (facetsP#0))
assert(numRows (facetsPd#1) == numRows (facetsP#1))
hyperplanesP = hyperplanes P;
hyperplanesPd = hyperplanes Pd;
assert(numRows (hyperplanesPd#0) == numRows (hyperplanesP#0))
assert(numRows (hyperplanesPd#1) == numRows (hyperplanesP#1))
///

-- Test dim: 5, ambientDim: 5, vertices: 74, facets: 1106
-- Checking representation vs dual representation
TEST ///
verticesP = matrix {{6730902641522717314961998953235/104776506286137515566857289990144,-390280997608381835232697072991/1864397174641101581339368357888,-28661454104902942433939732686875/599213099946195397458372881022976,-15406128800245753276724909349157/30099827352902650724266352312320,-1348940276325668500074957593851/2697085334304636491828689895424,-75141446950690149587759747394937/165616557448512116164343317397504,-4127179246096640907403890265885/41904101308994666459112740487168,2163445037703968508482316572317/8272174879337969639020747030528,-450191691851794455064904918875/9911153521006773458158718287872,2279981434851121779512083749439/17676604320812387688659128680448,-27824252277803128204504832190717/58693645000177121121239226646528,432519296535220584799240474927507/1077260866563456139894848957186048,3765428054130167415452303705125/222125502036585245242532474912768,-10521689461745455177822871802093/85753133669755553209037369638912,-1615106097963036509945373835895/10989631586614384549803648352256,-153859104600774931051039606447/20406035371046271441470712971264,4095357961966490213093023541807/10600428188431546152758393438208,-37556954633350929381541664662061/120997085920833975034182430097408,296286671065291413809687610328321/1053538385342145802784220758671360,-10037572016285418472076457697283/34633413130182945354295273324544,997853003413383634343654437239/28057515132355876009082573815808,1582119868326493885513507218639/8407088162979660786601711108096,-81569316783071114122199918965583/121058810132428827263998000889856,1120789335983075619818326397155/4134944818141064063042084405248,7914014844818656873551204740075/22266590606652141677541010178048,25236149475136229335454327319069/79666744713521896378976557858816,34436819126181817978168993557249/82379817115275568750484537737216,1075237021510564600709996702909/17773044687234827937178204504064,-23175127494363607504616048277949/210988502108587902922675355385856,-262546879473029603707919740754503/1622551957453061927871757233946624,-19280515772323983437962092252821/35728275856999840664389622431744,-218131885003462022022715316659/469355616747114736836851793920,-32956926220711635100663405748243/449227965645187879958530102394880,-995867185257037157829295479027/13278017017308155570887822671872,-34890800420277220429041632861991/74325346191872223061463382622208,4881578075620478161802803758307/43604575931687752643748567187456,-20580810229325881370831876823683/230674452417758719947424114147328,-18559593442819017269644306113629/108009397080680772287604176453632,88458460007120463175689674713969/220978085925343540969004637093888,-37800542459003165953439084708129/195948843722016143991203635920896,758858976223275895833211232655/2363840771607767955197225074688,66614867809015005947115586555/1881417544855809094232601264128,9691421936173397280519865426369/48606240382070450486122702176256,-234129344117733669655521557354443/488479183229913383369921964015616,572772271277166369007899321289/4213292614681816007099019165696,-17725253852189936798796711715609/73625179114097004834733902790656,-2467413306791997788510344280955/8988565927180644761690706018304,-85229655145351274060232228811983/280192885126197765903285898706944,-628404180535349928731333609701/3709936664406618205468571467776,-642280552501046512786290404317/6459758645495573539114734059520,-327315657767793243998783766913/3134175500792059494205132636160,3323722793616802310017784382157/8990617413053070246706828804096,1407830498052662806322035048261/2881077377630623370955825610752,35711980245977968636464561209/10093496863672162882410132799488,-1803712978445146309490912595779/48179111502188368715192056414208,-4272000456244429145286834506875/6149444231630612191866909622272,-1899152396436366206979095160857/14913336393539997612264159969280,-79885633539449837090851846233/323468803175407180235067621376,-167630647380197786990945515144791/657217633371141849077965278674944,-749163246137013580167389662235/1658390490567448669885812441088,-5658485928621560047787578053135/20374235221048913858078802182144,21899128469814271135939613863441/45568617888980620247557090574336,8995083771771060735599094535787/24917620553257582720908074680320,-4659018797769907867302859829157/45788183250861295204600509366272,6182761421693475070225869587551/60832267043382222623616446496768,-95181007750031517664768387896297/141956724170309861478622472175616,-2532845396902464879676112717701/9542449889035840280145808064512,-3617307262047392943548475040639/70268006390052649658190697332736,-74121539530598398873606680775621/192839276433825988720292456824832,-49879226967660407713626026523023/584973119586817781864237061111808,-49522201761274099182019521066591/78952935354379874649010377588736,-14421542476081339814983748532481/53054540356761081813775547367424,667312808131112830515958873725/1695014187846971091419750465536,-551749439197412304463264245357805/2671042971747034365646309863456768},{7965064963571047990661974269879/13097063285767189445857161248768,2229088198582871639444364821867/29830354794257625301429893726208,-1151744064633648628088591931673/2340676171664825771321769066496,2145096958654462396472974670829/6019965470580530144853270462464,-10525723290363199551782593673797/86306730697748367738518076653568,-31881212010323691645997977395683/82808278724256058082171658698752,20882301837660656381405878009049/83808202617989332918225480974336,-3247489166253543566073686098855/6204131159503477229265560272896,-48189951258451377090286032216733/2537255301377734005288631881695232,193422173033721230934552429571343/2262605353063985624148368471097344,1664452692166777604995330538197/3668352812511070070077451665408,-349178708316418480086063324531239/1077260866563456139894848957186048,-4137667214039165408457690077963/13882843877286577827658279682048,36628341422630448695173972988305/171506267339511106418074739277824,-104270031989255405946445474888353/703336421543320611187433494544384,-4287408942114899427190338576911/81624141484185085765882851885056,1044531660056480550698210323335/10158743680580231729726793711616,122870349026601083185477506377683/241994171841667950068364860194816,-3064238193395914605517733535229/13169229816776822534802759483392,7121989971065341008516212447887/25975059847637209015721454993408,1607187210814660852458271891653/28057515132355876009082573815808,5974129575653224600489116374775/33628352651918643146406844432384,-27576546548265815473564379749289/242117620264857654527996001779712,10895075389827495465482092520049/33079558545128512504336675241984,29490743568262358155021760113843/66799771819956425032623030534144,22882169707515382444261727501521/79666744713521896378976557858816,-1027137756604630684155139079641/20594954278818892187621134434304,-33787938972684823589097355832555/284368714995757246994851272065024,709700154888250295073392607591/52747125527146975730668838846464,-228753373493680613220504777671/1584523395950255788937262923776,-11609384448259671752805283511119/44660344821249800830487028039680,4951449166715010318244791285653/9011627841544602947267554443264,-5749845985147907118554383737101/14038373926412121248704065699840,-607520463481970650216388406289/3319504254327038892721955667968,1583274583421674849338627977345/5308953299419444504390241615872,-31391449607091359924774996275781/54505719914609690804685708984320,-11555679008077558368969397222637/57668613104439679986856028536832,-2608701412841068154625078007003/6750587317542548267975261028352,1481586222944295658850823100921/55244521481335885242251159273472,-3007045330251727214901715433149/48987210930504035997800908980224,-108240639898642410653474213424931/378214523457242872831556011950080,2075036137626593180757888594613/7525670179423236376930405056512,-5061705390633347433157125601825/48606240382070450486122702176256,-4515290081362007265787106032939/15264974475934793230310061375488,248186171944610231946489491591/1296397727594404925261236666368,382747169834834998861276615505/18406294778524251208683475697664,-479593785294655215245696025383417/2301072877358245058992820740685824,5135289724553075263003159276727/70048221281549441475821474676736,-13279245757276656800060650422041/35244398311862872951951428943872,30325921334999442226877079934151/103356138327929176625835744952320,-107839038677592866986291244713069/401174464101383615258256977428480,8878203957537922206357026747511/35962469652212280986827315216384,674059909643007544673679137861/5762154755261246741911651221504,6669803676586111336491149098487/40373987454688651529640531197952,-13773573518448295337824773687363/192716446008753474860768225656832,-231463636747097253361586597425/12298888463261224383733819244544,52000938440254358608783646009/5592501147577499104599059988480,-23562677087222900245816100187885/66877175056515434513600230719488,-21660451933343641272278053192615/82152204171392731134745659834368,94148550224331818453108899944373/424547965585266859490767984918528,79665672352385144871400198763971/325987763536782621729260834914304,-4046169211321087859962148351051/22784308944490310123778545287168,-4718857893936659579984888877001/24917620553257582720908074680320,-40123600378327465056037029412677/366305466006890361636804074930176,3966006532539263535186966481861/15208066760845555655904111624192,-61744958048349796754460953041/8872295260644366342413904510976,1462184541536367936666149952673/4771224944517920140072904032256,-517523027417861842895135851573499/2248576204481684789062102314647552,29230891278611954387812151423/48209819108456497180073114206208,-5879415085666183789715223121129/9140204993544027841628704079872,-820259240608481307307693487755/19738233838594968662252594397184,3928915770255929627217575343551/53054540356761081813775547367424,-191861865212527600876268409985853/1735694528355298397613824476708864,534768793649473304007534010009/869480134032237749233824825344},{-3342652359699167139022539790471/13097063285767189445857161248768,-1534665979565667194785895009701/3728794349282203162678736715776,8333638905581478487649145193153/37450818746637212341148305063936,227099350507606258171259080649/1881239209556415670266647019520,-2549189918856456235279214270633/21576682674437091934629519163392,-3391055148694122261161421520801/165616557448512116164343317397504,4126594342041522852890080789749/13968033769664888819704246829056,-25553532659770934830940661851423/99266098552055635668248964366336,5611796467550617950615919653593/13214871361342364610878291050496,9117296580989591680997009490247/17676604320812387688659128680448,-1436234466833322011675678050943/7336705625022140140154903330816,6702309251531395125465372072179/16832201040054002185857014956032,-1191375422275600760421122131231/6941421938643288913829139841024,-3714838494133717175410872105647/10719141708719444151129671204864,3601528570793678479734188322871/10989631586614384549803648352256,4124308810249996578507623462813/5101508842761567860367678242816,-8276048895725174667756567107057/15238115520870347594590190567424,7835751045453304933163940987397/15124635740104246879272803762176,8779360110100023196809004789615/26338459633553645069605518966784,-271312425918379852558854083675/811720620238662781741295468544,-654489203874942421078944987163/7014378783088969002270643453952,173707664538688098341627228479/4203544081489830393300855554048,40011078672628179149395582213523/121058810132428827263998000889856,6476870628947590114352801910269/16539779272564256252168337620992,16348225413394055150040310271585/66799771819956425032623030534144,1412423539468898929430913210631/19916686178380474094744139464704,-15560074540843728627741726644121/41189908557637784375242268868608,11413196293167785450985319588221/35546089374469655874356409008128,-311185536419451160233993844703/1198798307435158539333382701056,-183154853190352416769378505685/3169046791900511577874525847552,9527964545234817714684521427921/89320689642499601660974056079360,1579767994401995619705563030509/120155037887261372630234059243520,-58983201425072442429570396645023/112306991411296969989632525598720,7546869778526812900451591245475/53112068069232622283551290687488,-10854321791226946298691499424689/74325346191872223061463382622208,-627983529694430921828919399057/6813214989326211350585713623040,-3601608239468940029746861235859/7208576638054959998357003567104,53825733865033525595133408767513/108009397080680772287604176453632,-363638222449790908339937182019/1649090193472712992306004754432,-2772896998307541263493490388279/12246802732626008999450227245056,78704501027233941700112017267317/236384077160776795519722507468800,-4441638733411078661282566836787/7525670179423236376930405056512,-1753983517695115260870924084057/4418749125642768226011154743296,-1983015490753760259194503607515/30529948951869586460620122750976,29795778558961735941298899667/2808861743121210671399346110464,-821730522190183876889418301579/1150393423657765700542717231104,43622254775655200789644696661961/71908527417445158093525648146432,-2229190949124792244296451655047/35024110640774720737910737338368,75585774919127201135954373455/1958022128436826275108412719104,-99741781664673517904615660281/1291951729099114707822946811904,6483423107591433836871080431083/12536702003168237976820530544640,-108826859649811573355633862181/2247654353263267561676707201024,-3708115436847593006902357490845/23048619021044986967646604886016,3870532237408438547643378079637/20186993727344325764820265598976,23458458620045304127588238297139/48179111502188368715192056414208,-88963458112695015782697589019/878492033090087455980987088896,-3871980035763002169953081979227/11185002295154998209198119976960,-1112963538455023576784955779837/66877175056515434513600230719488,12936122378001387600986836735115/20538051042848182783686414958592,-13388686370365297202034685008483/106136991396316714872691996229632,687461546063083194339580772613/1273389701315557116129925136384,-25805768785268132047843931364717/182274471555922480990228362297344,3035961123884614377973000809937/33223494071010110294544099573760,16833355465486492538105681170713/183152733003445180818402037465088,-1964970727894378894803910954259/5069355586948518551968037208064,-711541546661849495453097131897/8872295260644366342413904510976,2753469344034912002988429924169/14313674833553760420218712096768,163965830150859070435754838739267/281072025560210598632762789330944,-851477847073082923708235893731/24104909554228248590036557103104,8631081558121231969359032165547/73121639948352222733029632638976,12855502019727594208123394296261/78952935354379874649010377588736,46112584465630507570974453164301/212218161427044327255102189469696,12901935126394094592232124130935/27120227005551537462716007448576,-3082081120485313123433067687901/10433761608386852990805897904128},{-12308668566877594927456986443287/52388253143068757783428644995072,12092521140967957947599219546835/29830354794257625301429893726208,-14083213712596794813343921883225/37450818746637212341148305063936,92003761259524493894298939875499/481597237646442411588261636997120,7679610560756775980194537637027/43153365348874183869259038326784,340123136112772266489599412873/20702069681064014520542914674688,-579058088656906051659917785793/2619006331812166653694546280448,1252597485851790063258510952783/3102065579751738614632780136448,1341426339432260279014904680505/19822307042013546916317436575744,1659762896389895632197456817775/70706417283249550754636514721792,4478680683340412049888900924293/117387290000354242242478453293056,-65818921124676081288586434901577/718173911042304093263232638124032,-17504752861222139751849206810273/55531375509146311310633118728192,14649609555556855783947125895503/85753133669755553209037369638912,5568244772200243418893459174713/21979263173228769099607296704512,4200029788309206216566429766083/163248282968370171531765703770112,958594560982315267324351401251/15238115520870347594590190567424,5866249943065863810320309015959/120997085920833975034182430097408,-184317331904959610807944722893841/2107076770684291605568441517342720,30493121432348587108320401197039/207800478781097672125771639947264,1916735419076947913376569417777/3507189391544484501135321726976,3027698756661366163018180830937/8407088162979660786601711108096,-7011527632143303094731617535407/60529405066214413631999000444928,-10482756091357278639441333963869/33079558545128512504336675241984,-965159233335537536425991905875/5566647651663035419385252544512,13376309246470951690545806114677/39833372356760948189488278929408,4564343773297586122039171452141/10297477139409446093810567217152,-15793527673123056144393926185453/71092178748939311748712818016256,-8221354054328067582757245587123/13186781381786743932667209711616,-2437813634252807701634925328085/12676187167602046311498103390208,-5140406134677040923539796454051/22330172410624900415243514019840,2458042923071639438442846292825/12015503788726137263023405924352,-3789094791399207327344993018767/14038373926412121248704065699840,6257144472287756017976805089781/13278017017308155570887822671872,-25783464512443489029188959979481/297301384767488892245853530488832,-3577775900968243460824476776259/13626429978652422701171427246080,-3225170165795349114600139599791/14417153276109919996714007134208,-6432237948330146725059274949787/27002349270170193071901044113408,-21980387199943016908202571578563/55244521481335885242251159273472,-1312598781402932334774240231143/24493605465252017998900454490112,771757724628389191007802759489/2954800964509709943996531343360,9252093020603125071725395118413/30102680717692945507721620226048,-3421517553588920851685218338275/6075780047758806310765337772032,-8036692178797023757605686441179/15264974475934793230310061375488,-1408547819020415528245761608033/8426585229363632014198038331392,-2888163320153565192792228189619/12270863185682834139122317131776,23666511371301369207549433700949/71908527417445158093525648146432,-17078200532190162108467175815925/70048221281549441475821474676736,-241158834857542647841603800665/8811099577965718237987857235968,1756547047062910353758907256195/5167806916396458831291787247616,3779762897251443412558726104807/25073404006336475953641061089280,2643594180021758600993361882637/4495308706526535123353414402048,7610427779577871392135493157/53353284770937469832515289088,-2356843676264837408634841390545/3364498954557387627470044266496,-2177432058634588749259083610997/24089555751094184357596028207104,-4665103526396426826200607265355/12298888463261224383733819244544,3201280806535036264196416067891/6391429882945713262398925701120,-161410608806950031197631073443/33438587528257717256800115359744,-141494208217268747366687361844545/1314435266742283698155930557349888,6021384781784930924475174971111/13267123924539589359086499528704,-1655063961840823909702848247403/20374235221048913858078802182144,8988879397989892264729788732879/22784308944490310123778545287168,17744240311776549859226518620683/49835241106515165441816149360640,3438262254087837165761280257337/5723522906357661900575063670784,15612377421943600808208043925621/60832267043382222623616446496768,-3256399937752352010456499227899/70978362085154930739311236087808,15792262983678394082376576142625/57254699334215041680874848387072,3811510477950869497280327681637/17567001597513162414547674333184,-50068557148236607658434945148185/96419638216912994360146228412416,-13800659133274797838172586886319/73121639948352222733029632638976,-621573735873454894736962112211/19738233838594968662252594397184,6030553951047768716924652303703/13263635089190270453443886841856,1105281663047271165454543376475/13560113502775768731358003724288,-4674699382006748209107495645715/41735046433547411963223591616512},{-327090633368046710171817382863/6548531642883594722928580624384,-1928047323478650571023557651423/7457588698564406325357473431552,-38234825510460434284189486821189/149803274986548849364593220255744,54747255947442246721669274435367/240798618823221205794130818498560,-216360402801461882898321272223/899028444768212163942896631808,4566015640403201130497595752895/20702069681064014520542914674688,-13909016036490644168562640661/25182753190501602439370637312,3141853344303619397229010260155/12408262319006954458531120545792,-18124542204482580065490405959935/39644614084027093832634873151488,-14618519850078241260678126794613/35353208641624775377318257360896,1482277437748954265127115941285/7336705625022140140154903330816,5565741942240903048049263617465/16832201040054002185857014956032,-19316195528464694586189049558627/55531375509146311310633118728192,13226066968572548015649946813/24033950019550323208810921984,-422274162265081296426255910799/1569947369516340649971949764608,-10482044688123605201385002249193/40812070742092542882941425942528,79009004507382614075163551839/5079371840290115864863396855808,-16351177517542272701364681898109/241994171841667950068364860194816,7726833739981561283835420578501/16461537270971028168503449354240,-9577766025658476406273604205847/34633413130182945354295273324544,3806812011997473757560783794885/14028757566177938004541286907904,-29412767781422284820792927730351/67256705303837286292813688864768,453393697925145878680721718811/15132351266553603407999750111232,11084083906989750041877375683/2362825610366322321738333945856,625915513667224567086643593015/2783323825831517709692626272256,9276379344889143309582885882551/19916686178380474094744139464704,73619022840303847331269059147/1287184642426180761726320902144,-55436741810449693628088929205/85447330227090518928741367808,17918628642906192932714862944827/52747125527146975730668838846464,-125315512742708073340229170281/198065424493781973617157865472,-6829019627452042968832442503237/44660344821249800830487028039680,27820155563622905575970201135261/360465113661784117890702177730560,15402786769052944015222277137987/112306991411296969989632525598720,-3009941685238754566062742757929/26556034034616311141775645343744,8441472729460510318066805448919/18581336547968055765365845655552,-16813803054441753609274960402473/54505719914609690804685708984320,19037085525017915915875680386053/57668613104439679986856028536832,627836482731444500857295588945/6750587317542548267975261028352,2335964884941219689184684227051/55244521481335885242251159273472,-26448241258577331568840314535005/48987210930504035997800908980224,-8863733085185177392277640807817/59096019290194198879930626867200,10918244757363458004678526681627/60205361435385891015443240452096,-397534564754024870790283002653/24303120191035225243061351088128,-39405774445668880447962453644741/122119795807478345842480491003904,3798321395834457122561332312559/6319938922022724010648528748544,-1025485779675473928142144807117/18406294778524251208683475697664,-295767610448657779564284015499/17977131854361289523381412036608,37253938373714922339479386493781/70048221281549441475821474676736,-1310013180864663971210484602345/5874066385310478825325238157312,1485082796732771939007534744163/2583903458198229415645893623808,-636137699796072280092797278787/6268351001584118988410265272320,26243836625394543169050790180455/143849878608849123947309260865536,-9672532141615173430492521223657/23048619021044986967646604886016,-1041222262092852620536555515131/10093496863672162882410132799488,9474362114689639953516487011965/24089555751094184357596028207104,-9870434576839988586954657496459/49195553853044897534935276978176,2615880358783174991751526209/6172738573485098349447086080,20777384628099139059579445752377/66877175056515434513600230719488,8171789141744145496993783981963/20538051042848182783686414958592,16147102055873167428354725952079/53068495698158357436345998114816,-1466999801344784960215450971011/20374235221048913858078802182144,1838616233596830228992213630619/22784308944490310123778545287168,-3049329233826455655772343254479/33223494071010110294544099573760,713752394277660764529756261667/45788183250861295204600509366272,1086316353360918510614733210983/2534677793474259275984018604032,-7961259621682300617364187955299/35489181042577465369655618043904,885863592338455966172965004989/2385612472258960070036452016128,285107274827231882339149210485/4391750399378290603636918583296,-13095012353813071510715969444249/48209819108456497180073114206208,10514190166624570382995689741163/146243279896704445466059265277952,6685455902168210414679764080331/19738233838594968662252594397184,510783025398787368764031926989/3315908772297567613360971710464,-7206245859511174300594948413631/108480908022206149850864029794304,-4782249048013131528346721170603/13911682144515803987741197205504}};
raysP = map(QQ^5, QQ^0, 0);
linealityP = map(QQ^5, QQ^0, 0);
P = convexHull(verticesP,raysP,linealityP);
ineqrhsPd = matrix {{6976038320585958943382221878240478210589188604952926241535934101033635557467978647586191509569440112961287429163974686923304073396176901521262558742381716915/3370290912466628758271228240016743168919963843150573004979778240431686489176891561326668587222289950851543643615658080165943671114216008446512945335015309312},{16745958713458550169401447685759404134467648205911570659676957726252573318604082200444404147278809153668975747298662157321090107148859076242924203760139100053/16955327604162665410680596397417175855135581165600204621562629612585170237768653926409776588009642272317714678811515380501091362255510515996025213775381528576},{16043851854520275628844742967085373752748426066917916464446039997040531665083915624378557600553180210516640158137928772421831916108238114013810511813787978617/21491554017015201569902660106372845256495438235005461169147075058154072782118349107983871451507662015275982038106705853811141258174096449266268579904532512768},{22557405702510878647143511315562707069285839647635805953204373858559619562037355060749913865279735159075007473540390661474775635483929589195269694187625773/16571368751348088036282113315595263603635209521119796445630589536595719933051115018812251892545020718604413364874169682809046219745819874037441848938594304},{1910488152950968644527514937303458473199752394140118544138165698860852523921339241506300245999527608091602563841650611519652523469455364925862797861269734132905/2510983827560896157847001991456602015022399838133335003406025180474542719888186739128727929968139273289867963173936319726979473843607115621425111933092177641472},{9011963002322290108024056024717348298262582583777989083971094251813067819247294441050853307140106791350042139518166530392044319631243326522898592966148307037/3661854316890830810287587222059289686992698719922999596761781266705751698451724565713985256582236341471107827925147860347254276518637822842222545631205916672},{74734186027152726203144460027021320085647508399848045703326917796417413875443594662263639590711794873248889273027316750047364678062526323327548995899540552165/72029416351606981316471979880440225973751708025878383671970961504875704236104628075547076665968760667194676284780515572780312367259009452291653990274360672256},{845906790760520320699082449155472409368878473149967134387871645472267233417927094833165329907562848772061848530709972551103086209707562794186505070754047237381/1182591783144144726277985614041538627177707873913106634214284044946116681899444458080587199967466206723678167706144898266792963394005958819325438634563475603456},{6432087965915836453102776946224918626712079298613355539614638281119164751067448579679575622738292231572541268524989446715011621248754155389454523680825457409/5398154375177147068108879694919065464953746539682953429542153678943210424948650100785666382867210543126617508043354214094039842912651221202981951077076697088},{507740918285844715134599210939257170650630595116213524560364266981824740277974479576705837371051068615160685272623688908686508361215427266026338380452585860247/275406721303048156448081973381667672591037981662363986320116540076984933489192107334262140700915179883879007560799323268968689913375410199408918763157042757632},{1870353715059533713643782516741784262706537185744179885891538233450729777237678921883497933210110280079100925705769861370215438872609738909013413384834688794323/2951530945289204079700865064876385298411309234813611960052193561673027687964313702465109197346113947750580358885280300415337915587591483016220459004173385465856},{188745636860474150305615254823426905157955124177571292513963725174281904008204284136114718595021991856796591919989276407561114600146168079277239340493265988140151/79563255737823484272347379101783670257352979381530418686961169382761239217061463627737499722960442412682607861600592219322884247991492659451038928814950997557248},{17525381848091802076670217568498072453287620559187747698860047326690075913594044890321477980823995037126845401597660267299295434612405281847944300612759993713/12889583701702291261094815855389160487194579817459710917391463367305131865249909361868375565955141976532537348484523008254729716332394101480572712352435142656},{188008438844738090291389313353531240781329361596575905299656672365711962623510743291130518226834204069513592993038943527884746185979382544667599665345540093777/234640515851618293442188897531952989098245850038144672904005222239779231330446472596303548891738958902486993383407947036299545157233166031581581596735779635200},{6404062185162024329797383433785140468310968954961682293448591991935885583775692341273832386651809346421020276776374849990942316463982763655836877711600345989/1464968787954230667120174383942915603756107513905867144677736995052467749462311484398861809485384376353363239520865984110880935734342015801831362058468196352},{285626545666505291142366525523534351208085308139864992878147404397244399272257283629092454323372965312031402413128255841819726640564238474143979184608757672665/89628645625603894326579565812142411229542171104763743477021198872873589297417909596786932880303090407222760018825019609162552866535728066677752431247871705088},{164992891996755579148147897016080155444001945476361466010755050770251014278935575151864913505242964840411379665686236288920784779130230041582301038880849563903/137090561471432185225621762058279505942310978387260505463434211611196663581585748582574096597116181416076444733447845236277488528513822015491482550648296177664},{204740616816318184327681570477514013166744690646012030492311421075617553591350090893800125272804689729897826279932368049584238443188852613160214484457811925737/319776765085959062675473945135754520825047791392191181560722746730121197077099080208614120044586441026100184762189318433465743107297263657468269224739280519168},{934165393075988289446557683838521079090233020949856183117369946237919653900263068241772888517352625501328579203486370873857119047889069614472216013631418433395/261145702541679391570041437845508207650430730469422379583123009013440796882484732104252064384710880232538838720570652728830745963896696437114862856971099308032},{80061866797245327493979013218044902015199165860990055598171091042912584247764019983574604383755630725000637856889826632018223925460950561155426451985820458323/18715064672288553931386935911462489511771826833082205520311954421916408748834245937399770661490936733097911138631474429176105215095494094757723646454384295936},{10660371538601814673110499847269074424594741083827505261513566348617208235510490012741996861854802457379717358548195360593190847828611652595628967186633549291/13759369254055392815444643285209240136270040825073905781564397077985140314446854626136714894349691421577836467793490697707171701807411081757040923644572729344},{2290881705913676560116635679783112001975607928847167352597991611132562482778182835674171837717559399263576651108745923960974457261684782081621541024196044097/1203812105304632932945285995628171660948566899739831771752457183886972874981279101134820391928670534916178823701309209805155566136707912735882153447219265536},{2785313680981728981075066361520661738306561553319206870548298985382111257835438494545193429083552562392202086779535708667800091794059574455104714719040778184569/3675942025058163447895051756015604265900025176680119451017492123485905130599818280402576590786732131949806908997931008389482668694414540531242230736082891177984},{793237949115618822013426129198727499205411780160485931960869562080415455604545364723788479390249521214395194234714075235481431032508319483980526191896498182157/1062829146873012813072300774642300933905570305331796051313817996534521534962003976800575672060200466467442062065991152909904877496466390099122750884929943371776},{5378728735547532314880164975132935426539561968471400847845975891434464440824024490801965010692656126506579701849886036743370493807551991760415218923428791284575/7257591819584835764830762573896950582046105187468299906609741823029929976034746224687216319165897448955047988776296178912508465234079708953083304381246819270656},{22405880898447523961242272119815093041868905487601276802454188437941853292181780857155130526716408583519646217884708618546812543635595794012259876925413133713/15829059745015116325259803507218605956305304724273202803512701770190314388294857790904265447854390618717561572576849590001609312693498529369190593397696495616},{37056184502010954806682774966514571012840295532809144569482106174623592033882648514263284188890631541395152633667947948079264463043584277870039462159257640205/14024190816053481667961097470496959389597672689995141222928659334202966628016059870448290546533040615867225137155119403904277387658205508840284593286132269056},{8123092004223330563837440717306231610410640061766179416412954731321184532911734465045592163056369165253780984313016963489309750196353684695767449283630556469/9259610765076694119678264612436293542527003556998800111560859588167813739927411296787834370168319497913504377331690743861131866068843736995877798763309301760},{6739766770576453415048122226792299824943652889508724681159050834445698108929903517688937974904139785132504106181990551757762585496562087685027369441312749529/6228199345399777648150176364478394047683888904009652437833514916641077581998837993007852994457574861401204804120955674324330795308290069133761395604938817536},{13089368897362909716713712683530705192553695118571150218253704703442526558411131418843021804488522313684613965548222633557823772592302308823864616867180731843229/8252986237614988801389023501722201315303595275464375834135951235858476186472479823625602617516378680396832003514495699340762488250066588392650849805836876775424},{862446935955739617864409864962353145172770004216581781225623020027321031511586765953890411330886240524972313600828413653083574239046636762563692925517154947881/931976524608072218650191831063348413240048893253069007273662564481511625568108081109676565419970552178288352470175178399667151360337197424882148559166099161088},{75676480660548510051583684686771233385869864815510392597955215006672552922478219464354180328358249281646285458229910283970059468214510757833402444381133923453/78083303865064186039384737560324407744249846006120293988706313925267028454604322536381333509434867820510188861038274685331088280110734735456698215783808892928},{222447946530250452788030947370033194676406681273858268356334980121214225242375081579039432427995598015613629681766292578643266699106366952493758199373614236333/79088280535657715816357026641670322981046860585373240036258727134931129976716110783526645093103782525128208286375500216132781943621117517298311546926543142912},{15082868964999316650013680716160024463203197906357120302230006324413324385295634803784072680584171140732893100642927592957501857577369331551810507883140751291/7421105310378207908147159352390331572571732514361491882014519798945739474488836859592957761221133593763214337897194049905094124066523555056598594312250327040},{840277739825139361601480299788721014073289394444072498536316453162461767867969200195272881140461283898524359348767364313183908835164773199087558005224265332427/1161221368346092480990732842849593828466847694629799484671936987168479840717289481460967902623388897279701070421676774781594187817690767520787009045842209275904},{2451680827646594238549038841495444417553405170292685350552770642747154545343999039546633212166861593380305115175519619450645881352634380489123130074780297975273/185143407774486633963349005028462363354357138279083696941814961441744726251940494217964363585319162242372537195878694406334888705568930697194938687230792695808},{5856972041150199377212352448293556185443191391505491226811889039539559437726823740012443873890136272259894380044572181903800080346038572835552087582043645289/6380909726192358256460656891066432433422896250888879148669795748903132329952230560845901015711663404673023697147135690506850724842148182731142347515293073408},{166209603512285564985755127988215917213955359506694732238369961303813872347645572517550998683202309965140218586610760401965315498285042224676880752216738170287/42235035457443964585036949602009301438732361506166822272275757796064473999230546942437151200672967996669486495731356785929558557316422629511945684770686500864},{2234604095059243355930315805950032534960851101897197971431920302003823782444018087589644499160571860082664606771989594406662822472089795987938096682956910667/3151850178215515994645907840267692857924458188937819964810094827862420662891515596430836667428929520693708852701453147944594901290665491817600850907424096256},{2548195438917164371527734034823022772882431213854400058374070346919255493712707055882407015639924738225900990252550024169556938072805810139791002592102697846245/3198968563480276298770019522835293440692919241044200539666375390029389520996847910454332220261962534648475011874959903465344211994237484370471185844994085748736},{3030440676145455096333153978396844516453022330172796590025219089180332608449188915018647628335240695141590395765810944196019946551305310372585134536805111603/1828192295062891763600734079892953074511739570656488307937321783991461844490788330069575468631626610097289100719559593070729053001079668966937538561008730112},{33114230741622556232852159504394620417016804048656389884455757197318540236496614277456683854057227353638188925367369150427928718464724993410535984676806090083/8220046915830441201916268752195364430101567695505824614546797773813052618595905615106982072927130920835588907286951009806549580676156939612266132453801852928},{22462389760856802963249465846770781765033059811028321018198061018996334791057730937158334479925931041663714180428196287670611073110627896873228693534701870111/29065618757438713172225288363833372071138136156341095792714508619143177447578269022411817610634555258617118160116605247038407227330277706834192363861527494656},{29673161844541048616365668982640087016333162621036189562447194473214999045818694160356894229367025387925570198478938048115874129671843677963371953363536706881/32232123953663906209354155936232307709189030481143880033564106324138034278643246430169071372072507148840665822778251355725879071506684479748658712386230812672},{1278147851975728114510954084700183124137927871517363964147574659377476980597750008674843887515034465113639894195109012503179437112051596044565566641527438387003/1663553674642744541896796638217668464281562173729492295425608567378668034224443827773988053056448123817768201082556337738975020565196108917636142144093977313280},{22579530929171691808716935474627698691597693850291966172035367780172707298369347918947429688725016669286476088732855688785555422287244316810121057305212146791/22259515806234191295085671132375159262786395698739651980313728624060329507723353315665549875343625067008268832220345020379514063404989743790781494113446395904},{54595805964057691830474154268195783371850518245755768433330812936135090393590499025878297011630068170246335699598924141272304521670476543092908830855888472147/91391046462600398251640427638621618596009279486630672617116085736566419354221453843985670643512657909983503740177690327410517801676542514719962305702634979328},{251403060579654577622490412141379504280998883148942723145727453359325203761285551351542142365161348104517965774438382335146445490939703252870337058045894829219/177959946095311615141951581130567639807915340593955187473827582711530174500658419503289539139215149164748337608153102691386558587469008913760918949180905357312},{1520677507303383833571887561185850862854310639968584762150698344564248684310711855127479097067915086246337454926700310802131273253311398535770793486601830696197/832326556488633531456149785034235652504652800636758106879948901728664287256946727294415014163094844724388269022380905906943419776629198368380009758247092224000},{19985373225426190518582176583136653196260318008501765059859694918072014512193514345102288311608140410482407357964023511366742369170686286930429628644644469825/31884849627040838168012807654303917207132835697391755885966180281503411805636487561404818168391808193488185806704493892693402001482803807482639514458377945088},{41734500361629514247938346991098521495709821150306218352027913835622196882044754527367890389530157876879378337606789499526589681323149563079711290713153195/43178779165668821877128104350336571130187731351824122445629124022533436352322541268834909060651510286448205222463925261192804583683920611155755088639688704},{71378582201632715191620757340618171033502207677657484250243606151784291244443733488734408488464989055962663039346908367228975593757111472360862574169909342721113/77193099388927649369226056290103756397336200309325919475612587198974933301572623957688830969090109256305844144750649182704929661904249832945562499245979805417472},{4083159473360579527599743763067590087611355888512366479031194555831382088778550915365487694252075713987240338537440496127165220805171167441303675782122933069965/1264999424118903606998167857705863849382043750819708175596420018599257372477285818439189031925807935713745208958464345503443153376716450115272166520989887234048},{5514238086429600089603666313090231347686926061539662653365375701227697808211368024861596916413645601104444847240222852223448293924879410560386561962426333323/7623160538686686506250964743466453364365668893234209645277779874487005312502991640665511957431118770602418539883814919981375581902575257939940708017535713280},{8746162207398162403475912222576082292458036581955156344454666777811710404732447416860430646072877096869125459341197903215214434370562125667205113688505156324771/322112258116799319306691280906251783346026728047891118396522832878595407031624976891241088913561053016676383909031889099791034349274721273047100549178098450432},{9209539452406413642974742491267554455231548261381327014663993474253870203645949888078809535574868315411360760386445482992733228970491328161581399015337021265/11806990138104453688684710106635701597625625101552206614691550167095428228391296728874095326324312923656953006837742160030143937072113533479547607375545368576},{8163885295634700648186955129976051617420120416127921989156965948254646693615019974292411077867321467544879186737168412640655419735736151927428474111895117157/7677449934194036426538209643601111345066820468992842963119933398803154860318127027631163547522548578237732765880448506041239768485894452413740461622235758592},{536493520950680894428454243849311206017009666208851169107838985947020603635228870255981416905660214057883667459747436384803891561172369112571865000313228388665/690261903617145989891445974401503711421448357937468525273558968317274301035117611906450704968704230493891381887140931285963173061833096568470743757589135228928},{76707834747195541034716456634786480293589008697337012754530186196515603743282339550267427730135834387463022098413593295498811182818654047722369729879710304097/73012445847420636179762833297542787064076815384576219975007023196208044688634980439832896116968675431437974171395233077429517376994728179952190313094105268224},{26761392080649222204217516599728542353203250700332937661247732886253766944292632998693982505382160673413777135119948303364029991693633221784453715586260455/41063156829768969743840736553971473795671938077752323147270632992173110942308873804195860774880836919211179168849132489615280938882581494016225037722320896},{26610200079187636211586469313522938843286631482518422847438762252507141752896695228955484929871606219154689656765530583094693431141291970021989064884617503311/32742204874850952793251095713733174915210325374763323284984977664737332363536141301842402352899292099980104110342766992174367792340604256664699097733299961856},{2957217234089470186971745944922504527786647444902887810601633107461999220844238935108525919015389509663678880483072853304076380312775328656203975422039114423/2481202169538554110092041797581160383216763120914021043751637597869400293239222610124858025361716436311182238010213714052550595818415860647334797133159071744},{8448046390598173967042858598438172266954025210393089016709292079063019306180019459046063058840774413640471690751866431738687550996552712164350910512630157929/11999619055581728867467930132374906995856694839035705252407946215587242300404826782832802878250677572145586006670648118093636881545357270617260119673838501888},{4267005893868287711912002338402409888798118099120769093310994418742845636858654797884710639946049704123898012822959385198026964936859216907298529494853405785/3258212292210236163547493017090586075111103435855961050448930688225786161533610551509655529166437650972035429510110134953595852167746831850843559146507206656},{330648540234623183044937410406913681197675702480356766142315913718153209668318720114790999111816230982647095784845040142807884424073965553421624705597565626839799/50127935079829908215631198370533158747830195518202580830785370328527237157528607605508059055538389999306171344754664520936119131320202886668184206744154235994112},{85519929519233062229957392053639234442415417669762231910827769965606166138499826767653634868809557209827770620634137212274211027873318237494937516839834037775403/100690789884799939661689239151458234132435534029697128541005688651216667745588752663853290212196943500296208991814745706412933491761842906587995535427948270059520},{82204586826722257509477549028320747169443176632830082175593461386412698475038534127151775007803515468007417236698756768137724930204786417882840038937629404639/15724271945659464920351316861038345530400895635177956639397599994813809229294586265841614629535441670833199935374172508911216478612132784821141858961196056576},{59105454400187807855686116571856146144450166600220693343338651552144280339093736157398942673989555485793372459529020832948156124965557960745975252602885139409/108791970171339405907481066916288148367139400463800282205178641116429175280645617637238159372732508533341029275688971679968213276285314130285935708991817515008},{3717100770403864541745511493159750330979882336762089616886424535974717719999335796176954858798857823032860415281685025341557411861323797209514986483958164909/3797409167718327215301584367316363151343163718970452891168499796579445571592855958107109506069888385787396267265910380582246672370890673303857657696789463040},{1793864414859148266215254194344722401128355125363331292635396513905095540614753387872209852355858071600876648634145457191391507064525672021907733892246960321695/3067334435123049494625652729712929596630371031993553177197673301111109698704673073711298145634379790716231809059754559251369882244013386667288698502879826149376},{440376544562298667273108474868200129417219321971499136698804225511070619024224344202363584306680916091502419692085885412943223760105656839568033433686071477887/544364890422410961130208875625052168360720789590003272261945384896510843673053148350722602931648085568135220028007753660891562743894314853218898807905887191040},{19487306393156687974653521745258534542026622724380555545297016360791386247059618560139747030612423773202438401560836913305691606515309036233237556412165421785/14980817972725489458622817630166876118283427269955252217611633924558315727908946553399125276467891539354983821511463722484346848565847025683141294949132140544},{18944213949840765207144361903013250513672343797995813898983962625646179735022655973035075361483276248879903890407770689408212828418518174008135898495661140128189/23517382239642966733397346127295682742685385770529262227151730798130841574357611042737809766532202782739181759421752652071825319974987831729345534910718413897728},{4600774367090931769872361319020358869390299083253536947839421329008396678707418002473904178959794784095024591372300184110869991500165493662968978277707448171633/1667003294768129582215104184005921702082280432144034504872766718723794680521053024679853456120023073313570908766743648782706583048133948052591494857888119128064},{42885299937549535096411294253182480272790657828169028574043063462727273216184829463910600488301627579730622265861198979890680393328766152633863201737209748161/39851710422112881068771387920026124267850316837962460807970040548235226920111649964499832401763268479734804294810200634617885157808125897768014805652742340608},{1894750630937248896948130677360100743216595361711926540852934787007893469436649039497686182183974672050892281905771039166592405133901707514400545825322182542655/3214769184816808841299022574273735100966705887541759714281513702482972605157186092295600250396406364339435759355265997844085483659951105807066703655574945398784},{11893271422634021204464726098229258383530226682385456691595591048202732920250172295979906287391580780742372682165986065018711302500628581420176689561786305243/5791578774016205430368546149036871934566502596398642477334515795967899655099282470695648829424265975395034775969309499021203209385170941032451642911065374720},{2026566507142819888818973155784584563283293792838972345254742894043278269030472785051949367116677058053409955890143542937070606453033851432798268347185403057/591185497566381250038342457346579758958109892562784637005638705076026741851584861014740844203446731455386277548473262779238439933868828714157131842004713472},{272590539635789488979616265886280677353806595312277769769524251579239360347761672948856108776820739514933691837132799590700350837177340912420604911574578650847/168752411153204378751665419345926258780915577424723353834007239222809984128578548879464823018157203227569498913916341263241585954217661357962159827654438551552},{37086407397771430492961130397694573008911740654555104425484690807050319388268338112000860402679132880841913754422213536589778365516036908636956313220416543022475/25048094956601378240257325239635897297907748745775965716413470430247495661340874583604408442703231110806918232464708411066318897044485652797307730062258653364224},{56441404396561961682578415881452581860109669901930696605642266348939978502440288768836692709245378534915050879053149592562326317588354825949043737157441564175/90538859701817691767706365397039420601540795712674223723641062829254774656461732787717747341521939214241083873995999173450503110055853777138023062025596305408},{119536598484646619159919732368076050854623330377672483736303329988276664742790798464363985255085136011438499376527515045089953655919400389122688196026940432759/129434261509778021117736904560503402153340398918623501025043029496620516855757229460932490326845769348787429594585809851177850677608305139797578434637329334272},{141245773863919833297030019956750472639316115386717454710813035777933286643952345787509284838548251851999351721224833664171184470375918536993746524421797841095/197116473934123607135229049035996300956006070283075433356753578054707147791461672017924176530648419752895828390435755579875521671662157724690266932190786879488},{942509371444640571080586829870999719360626582013321632953452367698400864539350832096759998284129000477132448539215214745961541695860076640720712230728163545/781156715892384181770597823980939851916765413807425015584775982528147563092098499527446654489686393396420282488761397477970445515211245204729981922552315904},{298536424389941043031779415493653434728223286453680254084142478149494489616048791324838809997772025215776748903101914289210487065682764453458797716331011138189775/294313120871814648544391208190701441116968084629802069915065359061290307731115941087961671435020649443094989462900254692526951758251685358157819015395247004844032},{233163795625336582681120870812531599183636922225008841115815888764109965921983901429465405574845594379021682331975810483096062163534441638252420942885238475/291560936746289638832832845970500505863772480627746623910903648669301803337716566134344460039175462102594871975576605367334488770439662765723785140772012032},{5329710288734488841969474374323963047471047893291820999509136619528445189150001004492866143817340379547980600961608061547732911696047441604679154690239565709765/5813311223469006279883697193679624546250370451155329200282923359440570609879990588007386231218852727743385626776682690835148154082215513421415334347439563866112},{1325090495111292482221605475583017263541289415369663288076748221230791428294089368552723864523962638168425334616475699759114461275236530150852144804995692034183/1663214124674820462975488526939691676216481208734038186467786186503756271321143392621989866949932991005452830108233041530503591883856412051148130774952081096704},{332431205136080227394386259284932166809965692178543578096648767449025839595963348623286929433367049708727156802693715461387870309154505490475276698116344243/443753989059028392494210911869758908972413030472730665406162456064558582469698978872349098012036126999201970689635957251551970358095537019920038288519331840},{520639466145208929273780422047682951261367796602652010986478968522200641834390904415155095802256438677491904655755968327630121653375050130010405211966710628419/924424825206103134120244719032072460936585288855563932887990651050711646086680862001685838578108803632619955285885710972541534254775719245572836737518063845376},{578344053687069451223708760090956728176202397570696891309508371736953450176801648934580982830371207102988068271696103724837885675597790061370233063312003369955/938080375456176084437015096758379622633795028707229531247205430989182014596986050553361148020421572207480113873732125450774936268914099849418107962063750955008},{16762576819367163142700022870879115446616543381566070995382086635086667165426378530347701666962594812372705605634316378688304030189619639625288670810161739164433/9516801945332987348891607099438934386907740016485559368360644401131559586972944221945146429920763297328368884382539757687454425050771268274776346647663239757824},{155949817737139176175337607895229512481538398995986453339192483821427867342484990805557185099541626576931113918832369647687965013842886404820460613292164582281/281951634084344447320556088878780313207365347534843178670785032463913367272648821718581980753001984675473885059389297426739246061138099524409039437962990845952},{228363512788548235412466500049145706060095795870200766048732488821783115756580995581252435668552449608549503899864034408347189675513901468948472647729174596943/5027631063609759225749975219945828166854675965211432437230461645917986414143206817519510051347952071532635592708657687862130143909853166675899853116425633792},{34899252260829465589434989398845060681492639514806327287617347551418365890803342417836942119808835326504920353321749185609382484638898762408576057085535439135321/13233772256156423105505380220230595414968591032948550774242899897654338723602566918806124802928425332463969409480963660455493372422239680732201625262692956635136},{52750765380644358611406286831991567481284577901391470100522074224122632723125598442451459053366013393021397923478023153981453620515885789823030697123143392777/66563922882800077103534355093163417726178075728454473529434728776333088008817970896966159389869946000834355484316290402665765087443742972173893724947320143872},{891630683439664127295495261390037869576250274105577374458210744935724823788852776454505831007701083708551465640697952614443732052547446603131712127643681233301/1218862828517026372804738373057039960108081429162621950956845694411714072568956049158786942124038709419016670905449578300169101950725442881682898374531494182912},{90089203920146310937551243190384632771759892453011445511856162715598750127572784727214916132951559316263700550659528562308848210160281206952599033800218260195/71912161343211762489301698108438750784612904569764052860551222475688132588958854257181949720229642265084177203691786861088035533802299732201291684740602003456},{245955573803568156799241818527691408440080371566221214646816628246407981010890920343900641043869399262945555704391385802691539949717986198759506393223941680781/181472094942281055863949737963539009188001907180258525610528225512937928715557169469273500561004583934839423805462004720055079070447692253802601544309570273280},{786877573097172141186999251495991810943516835525823279622631511109417241844139120782382503567657611275934281432828065110822001887978970398534666189889280061/353668623960343510580172943461629902292952741185694196471441709082808209346563178447796290850733705395791719910589791106894776659035120041450524982275735552},{11528165677719036761222978516067535729052507883553983620842532359993835135940230527792983609967350912748861937508387216753355207319443952720615272918264993247/11896128986181698113350238157822516667538804539356436331104828596301066713746620327823078821063336943089057723274442255323067217832493717486186554483177160704},{124109784638687269720651563116031737299773199854437080406596152832170848565255480148663945628907153142766939326361996405668470764280921957273947933765357324235/119336437388195016219032211179702448270611809300879248208150318170226286225998379718217506697694697444718868303926437064792334053445639578382121845800573075456},{179883973357213772577582208596458068960537015470026052454563980378441317892277677495159371948165376095247979826743119450921564981130469993005349619615632473801/35090090390162612648071570746664148743644274577708770781314975568295544204119111570543469369512033119062707846436806629663071173034796680906670963362566766592},{9295856555070354976708921722256450330736948932090580596091164620153342828612896055602288775428111417415366015386656462351600574168295516999416492861873022691/2246136692064169050776604964500787779506423141071655615040789820361097305475387415066613367685875235248982175807502369844202875877618756752317654220658442240},{79835725761026055413442414562656128582842189638948126826997359397477125200559179838291685885228289220027577123190168111966203681697164185840759048424500954941/11823403706184811572936237354655062752152148926702520966878086667173376531443003222332883443850256051541766760310600310401649956189164750720482392737518190592},{23550457440711780421389948456910061030591171135509103621344583510561368279514975545294371431831541335399653177127061643696952253339361056507615558347196001485/29649481142597554703676935773203585531793026845480452667576551866868513865043751825941837340125277986122658957611116905167100306419159505190009500014055784448},{866763712665993011529605743277886381821682608652236199286033402004095023216569512149494840059478191368419639038603037579857068637290583472715436833324947476693/1058505370815120253326431991372676236276512305068033121302404781532535753995187700260698120742682687395815454666199497723636845679465596637551408792417897283584},{1149482030956279119890771490764736429358925840584870868211545347806598132839795916905357012718005907675552811524473164427345869843036677346312345676310804343/1473014536141482796760297833538649716179689037462249733661655132456354183299283635174571028992279064864123073147925142695217738425025385608917127071033982976},{477491074593488164121830570595837175647598106092962574455518777650718244595665757928351108864657595489541239423141718819084014788567596089227114300707724840589/701873920033013756763125619146366034794356994400349750459933831139785227716093340480698420791809484900051456237990107314782194022920430878240284906921255763968},{4721372264011937623148699230145507441219494878596816775651926079560126922581933053978637400714017450343305074613194384212863913150223768949869438516058269065/4363101354081303480010190656801998826110591877147237733955981183924195628341589153752837805959651076611877706516474657053652209852516861761417120913467375616},{32121329347826685710624708910966475002611733702752401990149092570909377829180852857980642876661520151529590082028343672542103390328327830691116978210802078071/42875038982794504060652566414061869445558816776949535821482209766864505601109623519773399694212248941057681971334792751820370199330639392859485300504154603520},{1936124367850972563751232621229646017944001620230152572936484764134320515034523289854984760318217823214050037562178780853208083962598199211117821611787128867/2139353933718963182580519229990154373695202390863046836862072748310806561652548881066216942553494322370464778976166823442841048448751127995806381809120313344},{15933528924538629123037115673857636272875730856421554514174898721404607956171939958343275622368926304445007566567005294827168048036745345630325222776839658050865/27616516540626294473425101610426212206960816506148351034083798424056902068581354153584537637216394367644825367786796137492270715049265155006341813578672262610944},{1547728459069061496189827824687748277220218844148330921409276420889447029769432952807786501987676393088092898846311898838769163283931848974050303515886486839/1445412086201792164014372429000121560601514289177192978881122197332047991471346237705694381845398771704997092353220640673822621247366627808859794376723791872},{6814805130554345734345635642347969887851130245494693273110886993252284977995150630311080781946229597856642591696255122355185831398811094033777302127709847155427/7396091259090869463170503174480879992004568546532071299095106523799100563099250810295752106355517921492441377242359072658659762676641831591296256633804892930048},{1379294369587281181949746377319266038555404091507769190980092330707600771323640816604244719143233555424313947859101948911406188187001985150118545022597399029607/272608506209756892694053707237547876840896131582591284265844316090802717391215970000748385876349357837831567766533451776906821904312894944128886614910920818688},{610914706648337607434278772745115400243227353288847199640590743795757030100675534199426338624383639312578705649201980401188695782345233177726552134435379733601/560626500895714015316743399694106449638029044489478779193836238407813997288688493104195968236350566528258067880732371925049794307048349120943337920172673466368},{436901390764569280680158134795600030816828682164639458914938454896727283971927221553388533301060497594829481525160180822692689353389328506498236404633543318983/345023908067269269558359906086697461791761178890289576514606043632284541418844930852411594671300796331532098009742470152093532708148627275943304729663975194624},{28256285634944783748724186855617701618289087683980649315721660402607864809264680401266617694513406436860720157107375289352544292631949899244418816121221331065337/28618905637428223661751118243182109067487200098321228529290412132876510579356949651097812072166467551300135937868749799973174220820430340716684340396007756922880},{30533463907585579599885881126083084795240103421899570884727049376610785524559986839556613679523297864961479258636553969312376412957900852750669030912956250921/33752211740906352316735268681319653941027423182681896171415437345239433305315462168644532848781568931948749262157651714386205529890502952820006007522870689792},{68586206951383090851229530013790066236327678304725579518164352822125686074377906325633964204745701275462034566747943648184852359381748862034447704064544463241/39081560497251719423063992414554359816142850889967532689096218532123322653200829978248879816404026818527491061751819468338568378506251306815939003786790436864},{4977561129266752702737817892325887003115192464976315010390525600652524288060900950545687624617379146551769275261526847936964990297026911314117725383727188111/5265196519481837696159097251312856791356445104746745414962022505378951875252157260942659407201362432405075116855998889869290215166675066469090213748975599616},{2073915416029180323272538518431023418750134755930681202105653859962073396961442267795191506676340114515634067532148980675162692543370373897575367720156783498001/767788341538282304016015468501668955893374032317389855212256844966409505800161222061649244883911591848295165140152851185548165745900281877595078054032751722496},{117543347307370387866592760130567212113936672955157825119514558350387602866005905661764954819942570836966461444043258923538647945729984599498185836983571988295/105652006863126986925396089589354779544556612438992363276005050284863728402912215577188743030592813103482375106091855051519471222150497235890111869027549708288},{162368941976182018682553830584168090661254120702072742470218510307055437514627153331706203533349263634571651370264844026459345784103844204930542188106975382129/216747258073299659864168145267578039630778302893247546491098673196021851639055484077241389995492585332872956438786657604287060872098784034275566169854941593600},{7108665342475321951867989649812907178263618010702874801777551983330797093361184723138569679313379566512568568390944501262471765973680064224162592693562580623/12303673791599874413356240063945196316284968175415205616283234488055070304121377012407237019171721536885865899724316227395956530867585950774395346725774557184},{358018376807557441249770747919037616681224457808660282416384651699833258708240896093030614568186389280867560433992196210757626147786418317143650510584375414085/336282638232981737928332469860979038974660790685564994617292469649404659156874572108772556250695322563439253919279242210991645893797298162392247580545899298816},{1148293718004371869285141235240023919354526026358822486132368612082596456076829476872896825250025434887426606587877339339370530347495269487305924384323460238689/1547150368754608706165096483216670625206990085930291662646826667649348532374897057801962557732909833178975828923870425457014855328956093050201569989706444177408},{53732895333715129701680703433640731419125355446005383572409375832164620892406947178989391889788691282699644401344469550919854671177216089174788441315815844517/62740144802080485082801512662896402667117819818026706796818423935457499145107584607618019600381010388504052027995580624202570841501487787820862066116901470208},{274325373782701974326107925675259817232403418701636000598468163850817624316590369991196820647266997324312270994359872871157475682248731270548481643418516927817/464070993945486039929413796108004169093255540058095638483844740242403711071019277847576879419371701720161445471657618597188109154531801372640223797962034118656},{4556146060109058706938157688340023097980709707736623746278750408022302566261418546866914290351305948423366393864423666739901807784327505442426789910054192716721/5344247144486347246014853145481241252308902172459775498990873369188997978415999587883515955023827262938462375634904736118869547786543460517413197528898006417408},{592258286158736374868976383858911364021157204705319153518755973834472040578760297494627686719864347382377069202324978529593289812205488958755341994512076271531/336960839666957710885223923757183783461575325990265590388233728262983309789336495892764746659729416227460330226940815120431601135091135653953229700713077014528},{7890294108641458101663944561889247063318818778243565198708019960556981565318979054285423494607179618879843812539284976918324236981139850196059167343925633404589/3530636092059773335745952856596587313605334149604504567838352660380675276348814059499297582934557755287660468755972247040405795885174077625867322950132481130496},{378063180567472934089302383361257981001709277610889569758600648922135696983894908461906654402723794955764757061258709499082376784302665053332266364783698131923/220032430411782023409859158266878855987046865152329041951002123852392436139071909471898271020179947649551895539513585848030785883557898519720957406323708788736},{118430316161622891673349986995236489105099815171447871130388009191392146799988784367667152161969242276671030675619001924065161756778282126746474055245729579341/36952199525662777728981676989787568834873832408885049159687575014982292529277975947960280258265925462152897480509408952753917027215025705037627151727584083968},{20773148021948481468003537466828008753062902270644816576174109571460901799045495413479457952794682217651453914561611836916367365070549095463935099488763559327/28449851094070904060203079789617934057357232723256362747805341847018440736358456759932703371684371843471548796696821455660545223675280274372301787976947466240},{253055725245020979590908742452760567656859831292482552972981611692459757263888090495466661745895385509660004412806694021644276551912306052375442990722655415191/177671547659623092605734084698415305919150245845095735458969749457820275450005044459216818853703818444966333863343285175243395222430541337752000276036686184448},{320780929906794943017630681461274106180213664780746579655048442511538584607053277484947858646490800115175367679471554082519899674884501724696069564634281733603173/35482276122041097690021241476552715803502058204331797870708735293929202841606458950044663824424828792556331157744254103645535106209613523937064940629142169190400},{47298216010807371607030734698729982401684554379982409090697663856189529740846606697528355274151648340576784796044777822318060727050271333487659000530333639013/67144519308525308225515781161310114595116617507287868041423331974236793762003319419930658885294873912556502966708125704494266497522150137758540810954276339712},{2548074265558921093948552829877879218229747131735488596506002983255598308529002844166870382192417533515577344927666677924297282360219577236695606376198066110495/3422530062481979899344183045695255916131654176702795949687307712461720086599595910233859095795777014194831632256356746897921001494648765923412095910459742355456},{5731048943672390120346392739539046666491873076719407392299375377642094378470335793174855711190504113537801242826795837166071476094568034126949615146145403462307/7110420129445640114410694965445278997298947577964855717229084531250343749734256735007988122690182090586095973197870619295385925587181052751037440536278307176448},{326353218603647313207134450648784372584033550387861965262898827429744776172183742901295290258684622071539758653587911068723498043190697687653531166021793897918349/102920899500450187849831402940465087142479301521091561726879759021187741593801455542777226754672351029861207782161617688137082492082351985899111208413024145637376},{41310449720602787293961492315116094978377265781777095016812764730660446422999209888022290908610629977534770821235496988905692097014051507832298141894893925093615/47968251608954134185184533692000085121768288047191557691569185009637641573370146696750343651298176034031321298015282619229491716049545685302319476206413771964416},{11530455520785406275839826947023452043520862236486214959692201862706403960032694950943365061783215384946779548374941579929307112167488783466386475909799991549/19986701906217177761639562263553228908236486924012357870648816980245030182449600761691115064669133772032674831984512692039147982282225978012530647358746132480},{2943257189782902045281299578136769946400260525027142811421372445456188650698275076875428873604711232035780808846568653759802839437157268810672073146122695555781/3857511973175370839099517923864590679221029901063194471412214330316863197192320593985768808591330846509209411507082776280621640804913579425142345283686690193408},{69777216970272496772672094302383009380375018632420182987569164960909701163450129490834533359105164466770475076665713480665060567868502816110700315097237754359149/6278899005386636179282937602377707025154709146321568202986426873002182036666050993259236336339331554529481070492585676281306355292688207516042150037355809996800},{636399743259277787480484098975306974548572378398356736851981111502696635885344981515883410421156963993995241776880840152446850492996615127590402758898289942633/744880988884828640238991405800270910646154717111280593334477491240158752155095700753751445949491685882426490688026189418578417799089921424298372073127642398720},{368509642914645699859549786519596821239129117693889325772319661468668674591608787886714552152716169516693992289973743295248813233333625988154090109470034831727/577860119574870213328154922449219674108608125986641644095599970018725510828052555077560005694528284538937508178701108791747112466149092395695310171288165679104},{168475109963249077445690376287954182959901315973276969976974897555030646739844414128707257228507229217035330547940861758463032985645786212949425831978990124251/186285652776321078513714032787918141450771624282737649371283820309095805788469215994226300535034569292068451347114331821139133395484183719446497221623534321664},{6653960772874290390096956556414614073299654729029312003945227279909190186404091274812152660355937843610903387116475096575956450832966179315693868471959765429009/2758872722411697818568801722376079536690510040768869587304955112293709834474655112036993808548962287094806260245913907785416661774238696022082353432748598755328},{150511850541683303458834640116956783107846072779103323960739237782518471227431693201984104879899334030108389663142993866446078117346461025602920466938772414771/204199988293168923166946265965553205534718534049794784163698992926820195317038352106580893342894385781517268821830284271939586201597740499771249870954857234432},{74864399626494048907661997362605042875676302275973068960380162318015018679670940607372210469176858784658303369322431961899357863916808866046865439636138084373/6366223963151535413817697616890600447426711680611730008130057668979257384000703267992655748111717113257212403472636668556744839097031241819107967191746609152},{1556279231013324852868993743835906877266717975320738761945209553831860021501723755204094259129987197700387529394646398948997548743870861565385907968367275627283/1365249160463427044841137015754049691913378264955179629845280877411530398172122806002696993588851076330292712686651036542831027929937917267050197389116316844032},{416220512725100076481911903340258583781191489761683819736220321334482360474552691100264310169535709134103527258640809097023468308234280489317342130474802529427/398143190521579963958272916020665858522105194367735424468600003714472898253142496073511070896300602165863944987099702366548224972496833189025765707293149626368},{5380805738403170014475850695209427309037700683270193926027739496459157483143665905546324370969876140997852297369487353376242594713561651461022792801602626169/9393224840688013457056341239341332138744236753513003735518458196792789103372216215287293464526518150772618782424929346082983406548002396565802403965714825216},{52455545492481084544939450113537304408178447390297785725247154269784655492760605462823080872810018747538797215869764988657866178496581480425103082127566714117/35582454659721276562924701043910263636122377275698493597765457295791124718216529310363492344596329978993610695340673789154697856897873058669231214428378103808},{17019912111195117694774386579978437941455165579863682554895097310110410933809206001733612395154999768472415683634922771198852907814796631695265584959184799677/3482220246794676051373656858980560012450519327736807707156650030232887466242102529034795638159522622500099884555704390801468357032763952127784126806082191360},{551408739382181914619735696819376398093718565823961662704566993124087728248265651519701691188061588932388607275994870490738254784907579412584505082875583532979931/330215295045696312592595087314383722475348287430003403171787238342698293984276395608800153514902432214276976024799174688690044213787542833262486427237291445452800},{87838636464096179586292129835661585211573777855437719141657629227958383144625838141296813105062329120919030428860005035383676611926629777079588678017123496947/52605933371299922073876108942833127017699047778546937755813236197643245759987144229917339482379134137639847658074491786700586515862215517406900810015756517376},{174898691320193731683484205384883099874342628980417828724566170153065437744460843842905928227080456042739877708915091002982444775453606243821210543423720108247/110009468272839136841000942403282890271135980374321376766482469727623166954155517225091867549542510843037959626557660699063751506001195504639084344618617143296},{169541289781347105368733214296615772335533589719369726763264502686838337062842633900366359999222785618810940482312619095791194881095568535663150537803052829357/234433451781043717261893518452475317167834667637459943316893998367335674720695713376613753478381377079271228958096095047371330959363822824692671132007437172736},{39035482346607811344029186510114548464327417949688763202434998577867449243474590396874853483912880973233230847865190098584754897062915547181050368893404244953/49349213644839937435699288017449412058828556799868334146582010244121273177824837876409257859936177645246545527860072237391648867249495848110163551314881019904},{9627600438584174961319787083070113668504420953761377138423545806347682821031086043484061790176258467884440051617765862367369863013405862967652921778892135434181/1727992522982692369016817748010460138848273365995775981673172815795564482430106611500877492192444078288688419641104835733061696738313647347006115076807917043712},{46405618393394806947506770568530725595605133050486165549595595833741896024838252942195610259816494952251527376089202977815980231448799701629602120450075483757/40138122926360782813856275520288188928835355532414969019739897438946229977631710215066289402828712633090747025808373128003912429812821292659849258605345767424},{1136934106007667381168307014312026779237513122430952586159336149666136668855438814481154919249657620776375977643847640374946910551734161963603947497289820203/1591654747898776123124705485004229595294248110357956017121314693289816385241703130446066579506625414927885847517199597553847017639776347614075515405154123776},{482415623547291951893961275170340708688590943834979543821023663850843958836103669465331033689369479693886959138034859825688309746770561361754878791294320878639/462246359598132223247964391548701103254803723477228983796604979192140885293953605788362922888354239593364050513717114712817305183060212897322205799134477156352},{26806791995323504495659673670368218823511678539117714502775990178318368116112712014480387474967242121649888616443811848159961444202390257033859263251970518323/4690024288405253446052831197367306587037120970474491686044518369907264072688691164938580214189179094676602545643199386490879035436203158994176298074278723584},{539673310838310855764929486697854727047222804086365301652821808366409042592840532662570315580333546859003743915697936366460293237834494915770686606547734747/435066748861289686797798207359149471667250930310342818905998906749608086391490967079774354288600793427298166973677722837498667032543163745598384408019075072},{104177298809655891433420659314539932185181173422841098421658382689717106335407662104313117604636968173693084491944014352890946590077308042368553210773233612811/76182079263125321974423623388368317857349882515864801280257816171017402701544689117803728650198962913925614628297113249116293780728066895597954010922180673536},{23141093087014081740400915346504246280146960711584684334741459373640570734481279348625284760883812712547368490935940821542626892630449015028550733814616589745/31670542528045497868615890936766253355673548564196934537078281701061145094091465349787000312577136063392282733232711649384416689252755578710329028857437356032},{11885406889529175170422543320427106790360629454792829815380767876439339632288632703294162703552616226052620389516303924348697908537813222275951531790149273737/3477440731843439446309394052580148053498290007237740733039468719070111639208126460924197524064729857332222636873853789341496875172496519863214262981142511616},{170852868998765173407592066676070441293581118290538932626750687471670343565433244731862001703634466785054174497328462588583796535805440694262500918359237202203/90564411024722286446244477212499016583305303013105085059492146317842397121219345359534854031682674988286566418831732067685661615478298015918076696387861348352},{90192675660713852360431342472040856141228686916815493938914056775012659783213553413129646245560132191452544848359658738814959253099585710060149794350623888055/138868819608152216420581269165009697779022922334645605432171214136115722741439884619443208712363392678357570422002485407245288003065556853325952860134255362048},{2674700546253891384068369454663582601441586305370388987521985963877678630429873360171358909907764646659488814208902136153037708203121241603223426963505257221881/355954600888613002624874579127717894433385540412704894353144362788299208342714929215151864241726813250194506912636067963547146481803392194925223982033224597504},{6420858117087689664955019375298910810175158966844593940468588199571748083160903116387271117926336119615769333667998743502811706397457657731636588524356767902967/5581048187442771565427651340962530309179784020364733410623753542389248163181588903260866357353761850833739194913952784019968065618567053750118446039285534556160},{34281182994305092987089779973731922701236354827326294305467249745278202854176116868178300447942335765083716182873293688583518820489406158455646853913589568339/54080043769733373198753288884350373528039103976336880927747864348471720737681056831280536691966834159042501094601659138595223264532962553333983154573777305600},{8921557759937459839189257762631334119500236399930954529172977330556472102521402686444880509260477306087903583206109161155629653165691756754286291418658860213/7840077904785413995891244726707757007313645133650802730294804994263021560451810506814122559697012184978439644433097697227796955037638681299011516442639073280},{5371442290000365519981819746202224818739438736342621519895548958153756307131342362329882048963887643249801809481015323029990175311936978137454128679177992407355/4011859326879234727109616644895837308496979311236651942283006637427742847904323279919756656891665996394823768350348313608989348449676914896643284675344442851328},{382761501746808376122313433238496854841925230785851471310768660149193613738903620896284947471398651087740084441539830295137297429921779361383543059070608065445/381144618781291277986344420977378278685205424236178344648432724873593651859352688907610005149582963379758422600321360496590993464767027020660232718877244522496},{1190300924785314559221756887143320686539862835280781326785919299064215227130666805132786103363063706419923517652126651129737672154485451561546372093240741913/142558364616483780649129928812801483076874389647837211562675016751055723919196650663378737509682691233817592020305006787076138876861317245876086455224238080},{1921826744181376494367636865574807128527856364470162052368177332341980502739761299633861134216755266686927251573450818386122252303271149481361981192379769354013/738165874676354581091907155534242507813437014098744022283219711879651867087447417406046787571375567362459386495894963648759790862616660686438021626449734139904},{66183104911591956932198901127803424349819406262472178622540653928643180144882298691324845290710296029334288521571313093036386703477706540382926412460519846065/13695972731005301824522627299805761599856257028388693299543463169670598321108035416490221788886085155166993301680213432555173756519572485444009404010569138176},{161130853145081809449167478904873211584160589304504244737078073469673081880978729923126665078661665356789370848468628748102448108063966058418585455601059048367/274327889805701148048868162614466232045533280314925137702813202458196939809933542888269564738274113194129349137743180468480540317444368699765231049702102794240},{252866728973010255141362122523333111508622081392787630982910297163381047813734834460654884995368717288315490066466940195508205328922620921285858129903597580633/332761130665534629819815544577038297891461625793146125631721215438221024866631020368320416370260800603301683690565805040058474485103211269565957785989167448064},{16349012679205580278952833100912333603852718283637957750738194370137039772823159542749485516003491693390892130818519012593951681673606430612106157564734914351/15471771313460353698506453456451746275708217134502239061264086222207399646473594704080419261773258939418530768575535561509937084825350378947420286465546911744},{116159351411278165075164024920391756135045602999016369607293786853870176349722652632008309222092682338558136643948167927942563445195781970936192650636667681/225460782967157368319176281084937589398323184557737928017412938794072727028032468327418675538744918564685935521155894970632232781316593102981414553143738368},{2239282868347697749808639754643715577841154281649255871373041346457970133575517068807785479486805353623253328124904589794881135743819347014458323119549050503/3164549224412418429151967110880870832021778238086715282463036659511022623998320063893407102817509399947586603405892319455864913662661168114667388123951923200},{66488407475825787641406469011445936437719157957751459878372705169037556878006808951657154085994914281775600930648323991703990244688150276042625488280334893457/99972996917223009745221000115743778213569105818732120778147945051241069244224569852535196931763732926623066149754626056632642126648950117552043278229825388544},{69839469259355797309635233328837694745396999206195124446493030509455995546002877043130134687259797939855114185055337874668471166359830533184128597927729909969/123044359805024093293270299805625913350535278059869870129767176991850026050489123479620912458184440945835708931934538598878686066308445386211793223056768368640},{48414145473593491018609032619077852670310284906759031509047392701445947817451333337284716462013837682253430858941219992651415030804332406871042844865310983251/38684100516896192109718623745204664651277509897749846851661291187972606018399568273370325478964485690435802596309060187939559102287408749080534617270293364736},{56054226805565468979497925214268087955855084444286684431613687656794060797067757763788751492291426123816515740644862889671813629208874984693247778032070660435/77178220480776470175210828094833269159455235252635576475885706176676595861272336183674520089579566018158580399792861953546818806189617439469952978354865963008},{13932856624844341452217267262644251502198534324505972176286831383644759467212523821759570873214147970108262085730158500548018838983623851682132128968308021525433/16812820966432674961948881144747319614471864815616399712585813274700812737749379068406739860838024363791993481892266820188286504760412109654240233571331310878720},{47285031221656262399004174388950913106206657450618929913015422511895368764406420118938266767304641318646354883400305037844470977188298557040137798733505617631/53727146828288223897730788161307027133557036617577592221414605117638849997316152526402014341957398221780914975674146428030746433179360904669804968742947389440},{123411968035714173383091153221930791838966333844746742285998500571203812984280365763983187605808429226352766739322990639538175263892299536763414675104398364087/190309660950136477454566618855764278442779958267837457557310316084562539972558640190921071895912247311173761903099423739706510083040095689177916038370197241856},{937762687689896583218236682360836553813882279837378717361363490929761706362968715938600720476819829124574999932349066643356953839333054774541117697139129018753/1232806003363870462168246706985105243360278631586931360305451833258757752824961508959560837647416417986718398016266670062832071998729828265989268155621617696768},{249433887222434599103285751288820295321023469781749444552422908669362393817397715670150035350376922893596722411914094442998693994208316806089549047491394737345/235504325228142464107076248359465221465327876919764266691590372943183952841585608838569865603459122326339205170027529887692063146365924781484344005684558299136},{25191415438143994872293811364786932676469344398142377610395566727633306728025028496489636530469016529349055947837300005191928304670868675885388629878874811525/31326184420239546151266625182411860762156378727117999034891294808261968940847815591978464767466787192398811530744385328819713580183805610201574555216984932352},{35835556362960671000892012710349261438647619408350139980639077429082954885483975855382263444470504932736523700306136442107780347352413445253652921283236066169/10327674591892131734140415073961511473708288761069934303363101054214828734330623391938056147529116584012585627548356385720054940787401285342407693964301828096},{23674150494560079897982535899508781670595446935888427259869705184653199104407133399598740786164849696765521012397816560584321316232401985765450411981859775085099/30290776915015307383297665710415432450946207897846455452083013227253122254127140816672427517547997330526694426500600551819140177936358153478088115304629559885824},{219359521343298884471563600988989931647031698923320790284869708960287008800473490443691051146524356337536514942986067929209120905316672710259498759060326864779/128070273721003620727957272852339253368733708970571433314240178872043539615375899267582288791855443844072503080070975249134437664071156134198922755510953312256},{27102944477713594567523486727215103826901123647152304813997920227890266808802752493392366879719412348479041689049159816289380885628792506794963563288867163595/39294751172149525064664181923292290459416073995561041949057122522812063447000884735259496343313676170959930428471978290420399095452376856800989400634912407552},{2476173983611233801796753343068183172138498278939409979990410461931621850604160502166166254290642414010190008653692392959605295695052818114682174634785282897707/2088081457196334903075828029207030424896584518621176895455927208317675806820240895470947513385235252410953767259126587031534582715315734346260919843702295232512},{4666267008937334584916232304555034259760921028204068062154325522805730943819660369467873674491642674410100094993818528432596131789591808894657099584545017617083/4113429736144146565168995785322903262492610766230930448919456465400688074024285093896818317676141235070426306767834807497887897033640579424832747208877409304576},{694728943522400940301246317914960597639558741891752107907511917257765350129405389791193926296569193483595185108437890859565489740507703493731766420167173601/304099661819180998560644153836839937240180616690013755048470884258959483865977577754987987997510675348649132559966200311150320757456111619671988873363193856},{101014565472957932542490557493552153284194872035012033392545183194279048507882977692251962042763284374024203407199400720169445927434067535049932114302797140953/55524682415433236702322534476946212343924300353521220617998045458412157614794155285606875367519791481705262724033887182194117663056550641447461655949590134784},{957110703442229126458513448921160234220104582975745878086421565579260669041369776646607693589772086672653722175707911006406738018853708305446581249650086095425/853383514228767384678450124426384889390495951208531999423354378666979726705388368595751756301923537715932740418404564820716735044244505633643608523569495015424},{7787540004631874581269724736130777794534981471764350212008315055665823517272841614889055682518796053626170806196684309063639590238703353543733516117854375084533/9154292924182904944033026883232729423362816992818250566698099406664779533648867975679391171756297881705744488150237334781267757902497135331623438199669005484032},{67254077803161614791317161266604654535486678079158031535598144646243062100847241564844500288148251698037778786314352931987514270048627536459677381011428968373/85427373562533628632764350502864621568274397337084508801651134404276056159260409675829763166145141147574528703186369294474247826848271270297508289363075661824},{2998904275312392079899152207316891985759531582219961200408783608563489998929778448519234771088214115184052243328752246671014536902968877548376413739350771231555/1088108158387528952786393738425385623084888268222128617604292230449924983035118953418722034871140185111125056655668929658966855658124322541093399557377320026112},{88773306865434247437387517347060907679245048482555221486419407124536080303986210569936530750518566471698213341424268049892287182373018744082635058240780000999237/136693184346631104612300040807608922795479434139938285017404080065471455562393785334355057806526830861991753216213933441471361105884142485294942236309953936621568},{1849602735830054100903101884014002523747107039243752162939977317429597139052371472219893991933407756414927618740985031653101670107839229740200392891449502492367/860602437241356935927548589204455964777277604646792494142240212963503807384780351512292139979397928273022587377458887419102340673344801827037935946611021578240},{3904420297569804787258370846992301623031740958060830910603910398102827197314506846113563381665389266638636510626145990090556986715135920632513066868788525855711/3763405150781279161615268530367943108278594017807481431871961588470369128676133088935443562658115334809250972301686409684433767132172415704174529445416644116480},{1404858696314165550828042288224741465737506185403859177716771176630847660777565948064587283077023122147267369275054076340446786407586239705224918371992985202215/1981108664481881600020593828115787813991032076032834652898106782552145031707700890723301492366247637480844868850297471098754401578198096742740804931875023880192},{1843098447723904569890131226544350319984016847981752655896490693678592307966130453692125179917053065224531294851010832248662632424558730497005634638479963805/2427550387734211335251806504927180826140879745734198222034643014886583734561822358393390576686261088778647593536349187778488374702963215371528388060047736832},{14636442487856237386572646010940314179097269308930916222336687365415952066496145728238251214198428243432898737253527554332642055558605240942697692929547150409/19898351559058912809843578778039545710986890367115545661360181953909509177241965536648997024158431329870742774258868507756369588226602789533846925592961744896},{22879101203263081554100642104338630205620421025829369971815019129618631989835477644576406039217480346048333305294133075509130567609040341903386068958022693121/23310952335753331743955051534058673785595334888465564455959046689389511742520558426060986691369034868479660198919256760825329030170970869245389553329705582592},{248544137788780615341824278701680171444757296810008727082518655499751291714408298146050068904076312336666572229926041094242419517818051254013207152319599217595/234711129101811660927644902372931002379592052068356046590862370138274388039253417491113376335215512873202032577836018128186483309772552203587301570594448867328},{30413051846466606603126451333065743895123001222018520987494333751921384035709506440550147856333424893435441605252657944683645855300025435470874008702564386241/15549916072540834916847711638884282535719968317152213622966335072614850495096170424957575133887729795864302410509898762423232125170693232958660630333192404992},{2307867891687191219142651833505344755464970817590550867164947457069016006016467332167189546691577670501806273009882925358112950641752951452594989017847826939/3083360055542764717519598455608816201519703885702359743248675374241150195368266819735255294331431242989796758460818626604500819245236307930675947874610577408},{54753946300936381498003735497302625759922406890693367581919587813475751480015192396734633587759410289194711606050468126975077430636070571722396614798341195867/11424309329545926468265537312127300163760324892920433161539215284518860331545665184513036289122927102146613508826977189993727305752181489165586731592272838656},{6082769829996995245632008517736784052522291877881323687013995661461586188387361002314564508614579996921944766537262322076976382310887046465333435101592259241/5392678577812517645644996573157316992652777437719155116586300562158727765198558687068227145933298857936368247417215627993799122497477545593335958176232636416},{514877598493945152267562547810907116374006850984229215776150151492072134697402828699360025609469189809053109319950890056498700341402592061481034969652108087929/392761809566325729330239482309532140438752189257364698876917370737508937396986776117397743648334497644092865754368632015055587343078750058395301922053998772224},{76336373054221175021149470465339591292318131402535055217485317049312594814478377832856701165182624970507669946111528393282249117008289827285897516927021189925/71657572818351177036629642203272017825652352281912874060551818239461741031500101600848493542295347379956742884657270597154729908141332944199747932366802255872},{419938468078662701854032483791050516283249134873674043059772947844949321132952161657049352643655006186009930665780739503223701194767831949467973142690318343821/702706057763002725371618445727405257017399698152332657561164690234315297975340341091976250680303326778595604845091107471192711041046444411870031829665664991232},{1082354387125155839577275818373652777885099394221232111910855061799882405270826858783218621734345883252717058641916515658362343227972075121802659135749521272175/1420875121681103824477906971508679838482481584430389225806460021912710679937640497736278371535134358640995169568865521489562908690109097672541465637116658581504},{688313112885943557166745883235837330659909849935686163725064941725418239918233086909239011106177615779538931039678514258495709252287365956176200855728522819/757635390007741538265721552464276747407709911296859765698227882680242257500507062345047096138674354737787751468649916860763922174380344315353744453111119872},{35164246969197396999030601829878961219569429137238526862599809527139852063967271066004915000560638351874164741801680098835562213008654111899977634017430682495/33839699566765345057959913754277819323896303687471873953465937026553086189048357809066620516871453941722372933935422281831482281570413907053111551386162561024},{66475880839655870305008861702767474957777939941470985124699096148852462847743125170168323565141013412555770705610336770026137554077904136079886730254553773159/62441649391871067642752037395024462745126511618210680355334065106969649485551670749530467702112145642664735756763759962994513909889162159112227823966270521344},{7932380629453302363965461061954460450941489727298059677920838276366769190158254817969722524449064048492441582722867323226253172177458854446946565904859383356121/4359757276052358521896179837508682148370662572600015836199028474063373708428710701441031158747362815799876828066586584405724118034578331994861194966776770723840},{12014793191390756911105741253068792752216866239581772533457902156929882652098488634310276200592103932952785345050710952787789033994380573322416982801520603781431/7398945595065885338685477767021225100862194808226086017063625271842194148350849072847703805885524446595827566788781532160484462842296258568390839290889580314624},{413658167717610363858611715007633041437069033443431214387182019702859367582678308835081160183220297240744596449682612092961842037663083294794663917608928426751/460327356921271590243886722614657089032762357440145085044570540873258971366937825512463273130522571222016696668855342209868802226308053114447969325979105492992},{712864590162439043772097127796819048351884136913097286995537410800842296476802982458780694508638007151279865933069390908202538873448991691129346984016775650791/1121370662540968148189523696430761630724383287588719444939060465074105836208211684477967529349551173089414265813043468682326051422465215361773875451814687539200},{376818618995451910115373290288604589183868193662195574706909112159364030369130814036701298528103903116444678445898271105327929917222266786592650060711248229465/593085167560507485724048697520454919822792615271759149002604385146208276728303586602314573666749214756361731132360059976313793262573544013832268972991639977984},{220511987602345794201526425075074460097816996187738581605476884003050452647729156535728440469414457386367201904104550067913385578216192146717851170946558044911/64405235980315378266796747674244891807015872623020620712857170521843012161424521490511299094574473391358015790210910118536616385492331898192370560725527035904},{557829676408872843624450618291429899388695519668582951018859996117372509265386093205308604856460896153300272931380310354340556959409768840021256755356000402769/302995228416665747067426052777633715550393232224976508773874010407418896803178736225092696082317813343603797210778926244536008048370478311682305477919746555904},{236727144845138902440126843788240607837680174218736366613580805690205904386443521832207264821114625153079595281361204985477337442270206368028220576180914736013/296338496631567071012344235029581926903821582154116499842289473993587004706803007381358130617136175772836423413030223322694916373784644844549083672393751199744},{2930516948281208869819305875392417523674020093876461156873825735893657993298695758943916750127826884819714398946992283808300559757331128381935816427585057301/162496575961797066767329207319198040705426510299680352488367428938523778710679365145350323911172079859491495757019038469520009917649280231043802361487163392},{18205009400706119398204586033702690755948957192801175266370919798710185480226236604463272478158699334334274131803161034439958118672395696474867664765751341917/3172450380608187241588794535321981993084747126445182642273371220796510023894747384739110148845666236485620139831574113050308749397460938964562594835344654336},{8291463782721210282009290738232254433268329902500974597249874638196533823491982628660501531162222616505831361085469311519169562915993873671919089340291026855/9236048929783902375380952222636501936490816137315748909365720814972282325138354121921245272081797294857975729362843986183410006916608406981533800261572624384},{3713951347401975246448784773139441554672777362258328108577773846557719444040031930100707964259970093890243346527169338364417759893473714045985064272678133397/3250296774440412079560771576742395404527534170512272747035319901829010710051758631198536721240962868597474911328767565935787227519978925133074171974073188352},{2238189283046750391319733920738982453455243043865709570824220084468917381258846013390288125794578129592137540297518854856039095943826098026861436241053174892321/559187887296448248916626832969147602599172950552181082694568082256963919957341011906721493051801638390934146547327621555852200379679792694528069031716764778496},{750216550655259515847123665966678581404305688098401074879471662174720725421410481238331830049609210975939440639946137704621732259730599375350008637374765592295/225968460478686308934787254078874115252679461097300179254836104143151850357934891621712394258586375329117157804981935111892946931602802351694206782838048555008},{137548467861565094340238523846796089131021746456633733372403782723814937029926019434109344904968908377770105222966815398742497945367253434761056169678809299393/39542035516603799807995266870875484651274785724604396721291602581287641984847713745183520571437062331140310078501319075772402116149960576044820070483617644544},{638922670592944946717072750667406532391142488865717696685308841566483239051868081190898534233496155964571115406957584750543241805886623594991895282151966823/1149081513192721810263275242150392816205359419809803851918637883709258648400098501516898631640876343940022407500815565270589535865425445623081937383385464832},{21814370892970777688987590517848229493656328604914156846949449902557388922823710579721116841890527176262362459551805943478318843297530767004610387071148528750281/27987498208029066567030572400980294795098354724857727598600233333370530141063274008267531232355185780006019774198523507856179410477222347250118313879086624145408},{583069915968365406922354070336491749371842809931207407263781758309482219226435032237727781086031054375983275298471201245857654511168666138108407044418586943073/747953009198618897187466658383724255689767997287843428444163555414264516185570492658249530025550269971995363402575045422213316541944778220770205453794278375424},{63181742552784597730524455254804073137158516630007646282805892145772368957936002945222151248625475169584208789827054882063555244055176768135991707160358077807/44599144890334475333191040773799095949776210938499819244996932098284948982037304953465313704229970386469941710715179677929707522096906129332049069719796318208},{1096594126983778201921065699762854950944861794149805939339142381306101521383842470880959439560446099504524855483437690762876463617932266449464010752680693489/696849676841996346266053857008070478881474794946062328320560482039190507302279803325755726036837990587698226567839944804910942091519604599368525093426692096},{33191336343213087671291649858839132319397633692786195352699191827526809907895536807376654509239230508856170680807756423735795544388414386523382627988938596943/44418224270094834770996360849476455673770591619670921470498374101737026607367254723251315925310261406466726639522437539862276829981790421283138679828467679232},{535612806395913889244834129747375022649444171326964035216541392805312269040576100707302757068870996465373204499276425406192112047116101668998195391949039298249/910731646767575858166953010931674110737752864003280800174339439487789334134886877559519825805731698587828052494513905845090340440873279420586307407430757122048},{4974462569047141620135332144617363698909546599667098433507602053935389058813160612697871042301177935716826051829942727936146077514493105466468934090612002543/6816206803522106573166546152986812765321336091306503136164740398010909670611716306467152704541773394310918188208618578013084234137685820001786321069771063296},{8354750707056240649469189121401830886155887731559217380133523479195727465404720526442141791806555139607264030806350402566772694719036335888434767652821391711381/3272188180667868979761649097754970330598995110980537436008660231803211406683152792073306535532130944618666212557863780575425415679863629339539173772265108013056},{19906490811139324451811263875205480835704885711024021826865840615537057443297071741666143009409325606468354088355255023249646809002389348789887173130759218329/33093409748348536925948528830275828540159641531976024564940829599459061531699369710970684305755726296382452317019040148511845185275367835778960667837196140544},{24499059776081920067612460143944372263739050221447498536184437130360408688375984427238635010969575495924737874355362065875192566485374070209309522330729139/21009739738518348941850702487968795620489813245561281820392306333697037983704815915049011704460660296544499122640052459583752800476375032818420327391428608},{71683348691190361034385305651080030826443899235864609586093464369917370492370447625861447203673159771285775603612257111968652965525524395097306465829707080963/70633052979228482336003653591604775562133645737258325157281315754387591115524279120699227737491209627139426706632944959767475608528251833116609121527858724864},{9483590267587398440505408011543492945695618897770374484238812469190891587309577189744248986156999741245895113170001790362972863655613497615730235806652618257/16875902653937700518485695850050760902126079724905051255607425103501641131647280890822043373807678912807183844363523235087426279215374072832391346191087960064},{20705873828049544017479321067001464593575226647507153197855735044653217292698173061728675750626622950883808682626890883041292313037477460508321659591967074795819/14559634389809196307247423263882247848269939379850408250966770492074220720865403755418070774140545103663728631701428935805865104836271258648318513209609453305856},{445786228815597712486705758167444056954463463958868212736752000168889727282003559770333287400814738969034960339378376961255716921578263635359403556071536411179/490833535702313187193134058681458138667738341426088689281454073988512070042634399278700741750540789265615804717658337562081375570617034648632424491736331452416},{85855716733650735749285997987765453408217245340760844717298913003710951235094158252487119638294757417732168714082708618221762821111728674972365981839702984151/74941219187056484993790627713581520591762715114832126963435788556129760171227018214126349014711021988939796589131521424526706149676542378982722421979334311936},{240566315970249026052261560875142596843440915730624331648972507977036378604617877031934616852577792618827104621169626154542474817626215408426857965725973229281/198425673229670459011364283703446038068097063451014149012514232991035246409702240550967252753449614855082425976076165249791936218386280082572994194340736663552},{110004919873356757364204920851952292591305045422065686525616604690093265521569608992142937318083375848602299369439334870642538010473541700333175800326445373/5935788537393810918115049901758123887471910961801848931993502186083491134805152426248910652781175340787619177819834236660142427559089476867600682468769792},{149262387487978901077595074781252482967311741043286816674886384961013534376238223127808491135183930160097755104044993699089194365859342833406645984467714693105/186813330312520909074143569644092688476907326095041097539895856534255641455231015740896578069780675082002078244430803536439052179097925137569772992783896805376},{120537564149971400078137057363785511723245365885742299197019054841002221459217539209599405755814316587041786068457070949025027919944699127409004908518600194201/122557642937940009748114564081843871238440928506914788653871594131165447138306204583601127906124230261031687233656683947220336929123998370540566358404165009408},{216666377243938825651346071834754533257319726777968977199531407425434510534089759933496902964470757466656582967424164532406068683051468595506322357515858580789/316769032897067196051701616214191675055236940184306828687034687323554714523927195790051362104480513771020831465092394835476312287831423822021914842613583708160},{664780222470268826899484103402016224738799567907873922747647732968059099895985810551913710065851417938998493556614744573590040307947120999408681354614064987277/403676291690647862810416006957387085415030503507480104671661515173763125058757806660743706119016762596147301053008278743683617890493780131981574784344123244544},{1432539951217963409233610700124025678749792546920081492925323429406003418809218314489180856490027165557756339628194758860716927897076101069640711096474527179287/1211253338459439119636848227111320026004386833087846534764959364477411820247828217316879951013565307255937036368757728958677016999579142245938066594940401483776},{29116115052325343222031324106811720621261670298290389290781096780438507095032916449274938274802362358706426284332001419237059131610503423946608223447889776487/42854310741594383689716663062492518221606494359882532599402216527787673032938569322507348058928051591132791614781496424940441509111376928315082360575060082688},{155775595503282988533759707615454004009807741582647488752395959338998890709325497210425681335156625343939355790760181498698697298826064000188523546775331492703/122477358113295220490208013631579724933790660794518433369900551793060701985861561862993957446092957773939853037207940276862297272567450829208347835233155940352},{172425245164652877588758368087826628756364357011743173956492466680709900898807823284889514435537780243631757045535908099503917519270554323381081137553840062213/38660928235899401374654896423168752856122665544141329007552415382874801851350629483328090469115399640292883954135145786515586240863625866104232159134598823936},{308246547928966350962599450104967944959122032980776130080755604395689381538053337532135982849078880480381389738516853128692220976599066184920335582746963351595/27516482396329330776034835980704805479252731702054674410266422000371583159919430553761715175754752924091955559305932342149071250852612955577395427591751467008},{184306693248300410770331018086156922873086600622300482357301449938256281340172952501224116746323217801097623282831104741528419089340161573459311557975999294849/119291296169248761772216224120642239792870046084355595921130673357702699163856467228907142091202474543315074527738342727063384826793094538185859188244128727040},{98046368353757833511905914409189553307756292113539253928848690713236686183913763676688632623037357095537802297566531885288832663485163513697472569534688136025/103662874287749860893420298618424988510653660591448902980379822764534408846758014347278782382920865917281006551535496909314877532467484216199708279923307184128},{2390091681725655302166374864675764286934740780960436115400878426943329348976766764486091508865601103369019996949010667183279238415837139251229297017780536992557/1292216798589730162178787727770529833054706598593938419705081973719552551744132397448321420501603662464251889869292225517911539699322752076169854728354400305152},{1660007840858520098386343503276633487142598434073462618971815970339777608513719794012034278120670987719979349100095504421709636912643308632096647042507120821091/2256445459859433585441171761662186585992153472171061178072202770886339170617732107190603504611890620908104289938438816799255428117737739483749298358550286827520},{665629548303968013536573186299023240609506589404659549887538015238357434410298345602939839857246685614759121273542561472085776614314699166281008965872470510243/487438812977758916978507233992562851096000445656596364017161605805054403827316275464532240073011508046462449495265965378685549289791642722901832177277437214720},{101299269944098485255763041335841832743304676745533078318581529909980690504372161810722623651068022751626537059649969952594584688872175082473156923624721774785219/72693144009094602620603666635013958013633763771125832080967481052240087605901189203657727287455788007849992362554440010637015986124212750500458781604431670542336},{3162496757457862452433079400961284934386169709509032348536902119066257996065573265622830202617106937395553716580452344072780285102833936004995214186505415627/4471141963461591280849330419609239421384711642014670293369373582449964550845636569085033059384623582695326475906929907714601431147050083709098228083100483584},{91629632949677653676300441232085869748808607543722396971553297285354491815606435486781352753877131889101075767686489165694674686737253350937445107765733313371/39663615107341598845404400043526161421216729563565368545800522444158946225140304824210890353589478123737702137548472221468896969105026762762875387990763372544},{85745027592510450147976704278780328042629906942822856617558645179442012347868742292977697413471183016038571763871092460001920999782902642582661481601673167843/95810560835687467234703359673669199993693251716446107426633819980990475668395945053836785385934697961300155873334504235935789161292851861187563345553846697984},{2817661793989943481974217809407229655065877864993168098403982469701779393663713386446130863683567599394395312102168780140059005304372426753147495649224499767839/3907473748071851557039718689097713396731430369245040959613654500473326321555218169046351083832850326600338675689169497077185200428180675382413716445088580632576},{1014304625705579168336166411862479415063313924403760585122087518053151039405362438956623793017082596764638585124421754909392063479951563890985250963813750815167/1484595398538262097243183704048021133988241288378867881115395613532731913406568605823954332656127116029637850851130854311064912512898051099655706735514505707520},{79259991988826477728229485187413827104446296028463670906906068696965976322892785158813620096141029325367828524129628778232785989618205113041316695568767733491/17473028913082362776892431687654653267089737436569463668845600528274468345177092371884980578784491370036276589199587458609399097982622057089911023538572099584},{68368823712900206929326799175041697614891545235654942217453154834507094576783424099899654971335337477545677157667849027411366776599080841477233949757494628429/16773942846272228182323056500917890415089093427624060765451571532795725960658943312540100421651830444936055785896351756001575381933941003851531352048100442112},{11875817097990684597503460698577379218037472420736658042840830759223008094546386472536624952449590945390365567895090148494451229161139490864557045880524151171/16594779827437005466840822005535969278708418045922454913704474485121122295643608335792889099656421840691429549912090923363596076239615504720989905082650722304},{412474424750459810242576122904693540404696810310958918044532678575066973897379603466124402996522984788716851033569361607989265027184489428975461461484802922439/503185175296139666139515385174863208212172245713529356268667748278913383205702660998995073509538978376979618733981530560207425327807047997170915358501893046272},{12518278748723392517386055425453638396504614686155592197430281410707762941567612692477670808978215102215705593996861882207839648642419355921480340128506928328125/4008180463835994706856228328904252571303277295749690964561390234561527018015032819931897970234837265339852534131656061260008625822664739406936153086446601764864},{27408924270205156928616356017352227725940657338470370310284749386207173028774970547422767318244920871993832460858662440933597617147914488731596622908003762719/37914283089859827508905657476806515718825662638170763511924634623184373344893492629912651548908026815304472919610641138365065140849999464451643844438003810304},{789282449368898870131113392739960842077231054379741390044512719724492637885789800717540503280152856584090914444757703314913089526960062490107881664273150678929/1198446798407727436962718250080347047244802029579460377211286197900598941369260375931285293742368331670456163561218934687764924399099954295592596910568739176448},{20386192071019154699825665092777914199225577483378832724060228663939288858053972768940472502736278917634681363105279675567576880806131338449792912338193570973257/17836014930267882799035703382575842446307180641486825190482216370948917465510107399902029916518189506798032276630029733532056642278351573048592495742417744429056},{52297585430343043061512040406724535310797031283393778592412285912280417605440869524946051793823097128635105295730755497785022184745528888802681625239097705957/939096940924994026309816074927082035014768693583995152180192882820541712126716590419164302394703565123125603900187537702436187059856350574697908991686606848},{1765130632575284464676352314231681215269898545876855775951356273232797168219169371544430246441520305317894577118889654572936919448961931437521944995470059291/1437713219019079660803252207466942048152801605428425714152958005603366264863912219252904724466650154191446958070556616379203927601151537086545753125152620544},{10938672915721406966972062381955618561843818587400705658128736833198117386412846973363674163477517927776279489909844108177041143516995334257464961856391933215/11668864818455174962869622183502426373350236773774202441108941475199920404265851011085386404622758237517848374514590674461331512943958899684896472816443981824},{93930897356718697043774513581772842400777369832099283410081231477359784062813156721993489524989621909402818886471643834288777156630024941054665939750216216143/122171879499861334646693822968997866490818117628799716330294466641376569800537152357114415233839377609116930714050261794245047761714297708701449948495603564544},{1054576892090931081954155442308324648850161527867354457603154057913627003374414770580196148809054130595384349012557218729269183256326884993308205461512256566721/530975211558873410529292507075265467808256454377666144552759482163455637299785585216999121600681691362020812561161238833790958881142582879495105744160577028096},{827223487316801377204164085354534111896872331534665460356662021964993357047394753947923221061778236856146871408087445865577257866390552083639623389861989065089/32713686719904759124186836414884565782856495816425632195833340983530539098820511911552382573432294162810179808401704571957600689319334242312356071189096955904},{34668860844128793824439797374348771793165994938757531822544257334115992259277607088234327030098002765858771530701689838214818734616674019434488168861253266680227/27596807429493718080734772421877730935021301139522093372683460029685457444156687586110959283094176789001482187466764918012195883985442658595281210174921276129280},{228614174809480373307779387030458829706200747319831731391283348658920872274265271591785101171541598564881413707508018323832137933712527545900555072235726880477/208682167291479461902443930674847159309096245380926367117423016452910211864946922184189916349368909540577623831391717375724825332664519543432893565000304558080},{6505014793098584359131888343527052103149853858004365311882233347860423015038359009910553038756453659449370062654308703850966801979723395764591918381887065454205/8440650877415032293495654977395131263221989672723553195149744471436901650638631548791054078771050937209235932341529799056310318349241184526527518933563532640256},{93824903440278836657959593319936282150969502488792467953466487769901287884775089111227842412936395476451473824197050348784086096460635623868440746189492485151/65384609482779534675312100470369817100608695883607350720212248400275138416483591302013050785021545777714146477830205167947128061081492802158271354810417545216},{74667010234510803938624457690994707486384429840060019592886873962384457728506292454330314207290639845901447046434903567736300576023869957029024711403022785309/135064282586529519112016917292248676052540005380831817508989188265854078804480343723488695680940699923024414481763829524913256539064138285818941455812307976192},{202970758765708239535498382927463404082468557463990462266919808861506167987667373443141848593491450865586139343237078216831409094586591009478725716091616792101/114932420777526179089333163973759712078837530214819555281596720221043814710903138013003083599548145873363021163156451449088772048306500759667761254864366075904},{3817690942823940446297096218373306987483771655374216148139719241478785280981610758926128895245073214849885029632661090919488392203354411416198071225861417411/353446651287589872508264239287231804391799802724872093328998685652249421985419818521015097398455078502028561352454509010348436354678534657887509236636712960},{32334870956051064880043573102711753468025377715078382077148817578203614244897074331431000992334869695086576147426966728286283719513809796595786387516958383937/38915521334000406528852877149641856144151539574494410417155443969076834578446935209609494746272970779120209738718712760214811194821287806347777966269028892672},{1412160001194996757953921049868805993167881373621088290078019734469372579126096433174505907712955568721167438389363009242931024525516842661264975458030261299797/844125051705663798690352550599708038187076959663341050469514484356269529118088486448087647097415069400006069450369646919011831631294223451981148812288511180800},{3648568145308880199413159130989673685028082905161587980957171121082552499187943201093363541187291449541748202768015629137186866543568432839265306926086761059133/2195761787050557498545113923889392365533795419076706386070085236259228643334970375257516339460964573598945103574093113526888600155852139197104807924719167733760},{5845539914598326470748942580828088701197629468378915559527386907901348466402144477028565631493567858693658903526289880263384950978008726595400786516993651424891/4858144670176267094315106883011150364114830467367373037573211384276614714951295795233235887753130851498198726560104201665827073642175615363372394619530196811776},{3053243647372615378395984909679145081078569940565444225674331860828649580321223874312920424900944149072221268210784021686719544588836670476531853617216967688793/2711191830886369371198449784860164043331601495764159939823118464247759755112520424897721468178875974236592890103530564952841619064505087160449271327613536698368},{140910987485967916683735157007144780710437316535651092835209467300233707008044999220469116596139362915372858831042314204582096717999795401351598965971782333/152926081853527020015910865735637977941593923900031401643332852370448854978935377144756519696234100999962390978054702883476311663039136113049317806891335680},{3049712569201684986250129883670919681502543556961463409620448000872082762699503710115424112171846621867809370486922653703980916650032995051705363729607291139/3698528390130601965977972730117641606012896965201176652242498449927391599704726903441142611162612267530677811172951667263379768804375246945218103325673652224},{2803250239437750464764067801529954704917353274225578925641088043228322990417979565522123220664416155776373326403412313394485292119925613969914960198450708375/3502721757861657411631035134877572114042642880769138178235654383799321315057598029116594948039476488559158515539987242778463032020397804140106405299951239168},{33636483928174619010052278196174413966433347898624778233347805218347046957117333005248265073209475297086983220956169287157520232463849956840243467995794943977/30258038677979531373966580799269188117290155244521883709123497795503311488663877245127961978985002241487718289660508575355712730183096691110456128642584936448},{4987283021300783919415220581593710161411104405770254181504846042742925225179771320557805221199664613632370150635166680480220759847760302262873539892776994635/4481650338063683064713841462769611541448254810469628808155375746969573475581542450081511094859364751252044537960385019345727745958626866335734225701735235584},{2081836071067512722507154615396145586889914291871226983778784618866846485968667970786646262162676235960979922532639009016491141704756305200438220437899759437/1104806659417600158151347376617109495288545629701186032508471828819030265277093772675143340800156478168376116895908773049546343118193742599650374489477218304},{6166435299925581159494125752943768612608529408344523139030266116146391798920244232416530353873122452218730155764977048136696330217660794866086264047969305021/8977051474174186135609754937502101731187651643835638669573996776773271922924394081518185358305026083585225905837509971702906640576218404819133281919261736960},{1989470581891885189805359615291910370683771797794744122936202065226346418112539864360810153809106592379628690509187125501808718762770295144302619096764709536731/291640150111542387722878583894032159075435742285746473555183517189064292538950833843426593412218832579869517006557207822696915441915791436621818917909012938752},{1414677246860689074455427432786443706135707808186590791077371379910134780435255758066463138392925245876435886968602508535291882515593680537889963521855985546057/849416486752580253729272958213889660782823071671034099536895892146215644720712627027091548370768838873616740816821576088289490995149951432781479090355232047104},{167743147304856261547246825017174941225750700476252395975983045504384430501771861045133808272333675995443157396282900894697370956266185026080707227209242509751/77533581092988116160487510092131373533508239855124914790852699979920874168958876638713894796049672438320135851229209925182598957069223238065296963893645541376},{25150012053300036413042191602624624229024298929133494500277089853056921509791263002986229404413369373927094960802561349171070396387193135553082305641588911951/23783058343536895778947961888889271954594872479926392690272607809098004153708916736626337543531277022323237255051307207549705057598580266026139610483460669440},{428737244810441584241898389089882301306643567824476219956681297604414993051291375668442455673006163571563779800954358127371176175706458592429619040973009612133/104120992344306682311501392343959227778090924567604173830601093627782813387751635965856985407072632093977115179893847268154385332402654356310186060590717337600},{408521177414716820752436658460872934161771212420912996488656217921966396000785831126941802031782458765589357205356208563927195269384453563876329891525938586089/73536356767157920084847362804095241072971657872791739202485643592759232621859308224247499359294836182799819424370444839296814866442049839939772737525848735744},{11538354833005143020088679723843269026171292252383107508247255211512928923947284909623683894391863172729800608418473735477845809105886255598319749172036019/19431959614118387176751914653884815918286333711453809434950279993630347885320049639784046946203110198534312295612322630837362397760124694706649414993182720},{13906420273975876908591063135863508983788639079623845042132900182817097953664839619710544122940909799142261708134669302771485717989776547292034238961307136791/17618653165107272258930279128061376072898172313175082595503756500930040770018972073613096837755034988582527753543797236276856054048594157167328388456527167488},{432011112590004124057577434555775821922061057165182104584934386738330563386283513893140398568460362892452918734980148464264314965335820231997076032176279786523/559804792925117211680833030831567795072725429406737659578013348942620711455439286811278395113915986330662027780887785303807224208984665068167293333036513361920},{974726859585116289012798467487771387319388989157938340149234513618216670117834127385427919869860748208020977764758519880250098228054268158771217991560149752009/1487957644924568956819596986420510056265094945455290773790025643246346582650194427553623002863147383044326615940865842319959711050425642461196463555341491109888},{229002509413477797770581288789925610298522012620000948364336429968473812220200169034363644940975813773203698264202888575263900162042056684162873299107737489023/209995451882268015681221887134872120939811119148143952552410654293617644463655606292457672400282939326889824456977914219070453296737821612410642673774010302464},{185985592950549882710526984312645528897598063211397752626200773734832371463361588407793042523996236792587269157863442521499635354605562716224836707676575664295/61540043201087464043314420069193270659178217507055596520635876556153219245346135139847197229211500645502184852163838958549437073747322480751716062628952408064},{4212032841018832030863092672541918362512195211454091647503062765534167318812916258044550159892123990976379119837466388429838176377941365500365204330506232953/4030116163856375921234449171893359169768046890375727614099942292777465889731863909048341380279115362155442756380167957171409994988139320938828125332087767040},{100759895423985668514043465653399931748750762807963313731744286841018843272062512431287697790647256728615119095191934206180887389395345029143421057004678023527/117803811920696267307774791494489020874821778438998952691994303800590186079327044992399355878529814477420877456477545320611024107181695218092039432573177298944},{450344834882067103768460520491475689412536687126617305931790325671058520395564596572361513781053921997441089577437664887370400338944851415826177141275310669135/773351413026754276628874515870885378063210680310306996062387513240527134436617285334723855386989322928321896719550765362403856898460982631754775259980614336512},{4711528060861211328641777211458265539004796993714679193241441215885828286352356429583671948818699884889746243366245873126189529071121203746400541439719637981251/752954403573607248917117732603397318797894975390337063384819074155065349842979841263624259460811931312819281373601070730322560426850307844821434277944400805888},{1569790584886903271570831633862127455564596122242148695207155217998208875215244228024607494677437044743468520047520680418313707108365207686172359895995471919493/1568241575252710496349806169257405007596021253456273514854677150721395369210317785161497416571813796090918712074300800670580982898092705434605443222335842680832},{15341387471127765600640304799313994926872643063406500795300158609456915493629317292576245150565919694038127477640626837631238466659808640698014074279878274139/14990430279435696892764046418740697312686853105847749786904269217762180297711126301552529150654382815284258366905704079020646753517755291639178958507591335936},{6837552857248920551365653419712713958289664936301171073317989160254359716507366148814990945550956696555091104617716676952558388521444233019020605983637979641179/2030880936351892825131924519950924334549721295340810720675696620591691929732462416517708215261974896103382783017375901476053576064869327826636447665693854269440},{1057801328044097292165007670328096973942859485567855846783287812828851558498474929635647790220062939264392900432611023450857598742127206373379999188692189907/1165202321856606727939610874758313793964586060425630233272958476556614473423253804264838826796530894494619441952698774321731678544484783121717339171728654336},{7402772926571655117805055454151870044668122927632529218065103193513890792657355420043851951269603513754376193678567403075254433430631798500480520046847012998075/5007807918572511230781171926114435386621275182715159300441361242794646242910916387696886628293233101708190467939145916433432387166050815623290037612124806578176},{43209921384633326348301104679016428982415406809796502147472594169927154286771770357012326890840885643588824749132525467917416721151615215378809390941960336093/46984731347227171845604881679830031897397728308695035985635843842083345278600725565072684387206956450754135915308219470886610159700802457793760160100816781312},{18284440710226104450480428580927780845936992359122815324892051500016586847611615135854639331901974001308157398173601648813388931896544350848687155714608531747/35201307704696600075683511144764046883784374583907151206503265175156995119295824481128566863133318729958328591529283514438961035397332648885101365632578879488},{7009567771891752081441161130687710933098191416032906394052154804533975119094877031525792971484913034735212436322150439492682294288218873275382876630343293941/10222379200956624168032293200390642145167795884093264603361325857442194884796064408392107400676224641511865771392239360274048340408229565049579843025522130944},{415634201040356678020940354145228121042019899040343998759504717472679297949835740504618247101373310994472983752284254004863369101573058301703507930191229328491/558949876847575576626988758954657789813899802143735856799977395368329881636712306131427474066877067236386171232799219119704599115166974993997770010081124941824},{170448463251745997779313702934528152979469340384784705849764496692556390338385813993969527873864210416335422881584466707867559531952499483913666892851587951983/109794488778536285760521826517022103876893090652103740853376483176371312831957091201820009650365377162286427902118314875496671703076151135349570870293204828160},{667745579553967704299780420351717281822333436923168472089001357360931798886166934441247320950135808271679202992569379315508293744903896856496624954408668772795/798592355519964374310438703020589013191277663246077017825040219513090037371785991102893178604375274957783734030397221765860691718596071662760468669065097904128},{80605504553797799752946522425258812706362181105725607675302133504804135023138730614228489624334077531640938537111281937407291505082365780711189856528658156775/75936662538889605896858189294150001594384157680000083927318190746328719941472424200919895973928511455441415115661094425022744629100877172915035700681735405568},{4270328409875410386521065853719780508521471858705658196479092174443780257322806516579135743235117946599038385306081249935249903007017062797168765594137826171949/1305015802687331030528971323097559475468262273092503274140524195240091097418277896311513080310269584598601181126658625525192287364257896020359194864218686881792},{2449374649313631583285955236570266761692702395082626316698082877437353110815349788234202163392328828589505140997202703072209130601682641189876765366395290609/3110073122011170646577885852848413529318807540502753112983154223981369587194086425267213868078668885255348932595086632591132590153075129978995658326460596224},{60795570216325544076342983559268981828362712967786137435031875659900312392623101110357385742498523636004469779500934048319015150011607992141188764854955891217/70187495511661876851187274628092364269550065769695063195555195119933161273880307839355573026517386792610387091551397528838098236977143612572822259018097491968},{2329873094874202173915121252145057175789341316109162762579766946234287029448583676800672942797118568588879749298457872840126632015085326891412051042697343143/3237676979048931031001653173027946760510779445541088034006921731698383084780649506681464012685829998255942163840767570697128992185624741996568799111572094976},{6230038955585178329281083059942053392388992777073359691608427283676889780277080723684163808096162298761675123230502867024131294613508467312723990267520849001/4550681583464943414535191330193757017967039508860898889499010149723021090491458195911562680631553855132408423318019677560946548521839022306894808525656031232},{12940543951816471872721220875103785908665694936223991661207766940801384634894081070001249654074485128289613679808053025585853249194772523646455209155769478756895/22749932761509345756010670725429778197420113271490480096912998901779779804665587821060904577037049724897307342333654771190990395285919996834548776884730832355328},{3361117938048359214668667090606329974173521896897292152614757545062491630699761040989551562046197773116468429421376250377725060257379329577648089777490384238671/2953533254505294794588731080400117775979285565807255957285489610858371672927503046133878627126322251464700362896261270093177046179893269882343157663660276449280},{223411723667510489951369766441322959810964036899664723044304285130056027135353178342253969307003932970116359300658761153917705574149061820598497881914140692459/87583062870634102740675236011672943833680614505994399694753024476620363158210431770918439696211287512413790497508160525443040556745078284825644074521371607040},{49659267752349386272295057076347235710994725634577618019622043117762265770874428292220228976562931406515139833414913787101166637060306169409192483899239474373/55841858777494536722491127373870111396153697870290920890022249169297184976562247089909601316958703608552871475545031337856396988445098863652974611549030187008},{2519667937961441387844792299413270264233923018841105675417140663747341911901600694718445347820427298443463389726647982794216003047364302600554180640179246211077/471669591752756371281716217972759536783348666190440753502768966795769883812804240901482539676135214569334665560586523630992925853774042535696348692802114158592},{25621603066167728703573353686264205621482524518574695152352688852032198047216004526379270987887084325397739732757227481188536681690761741213323231659595203437/20261114053936725401645619552622565939238346403301885892048799264957177943136012474602924780789282499452523430614331674669890477829144678080122827393076297728},{349147199111239484286647922538769509572090397366340005325628691507806986558782090515631169496664125636178291716217802420600986493237836919547404985938034660611/92594980126539212062964313774863252687521571664621234681335953471087990899276041635608949088919733147072305657557741224045173795504201910130585484327085146112},{132672871099029901717413064939999893728757118486521850529455026553509390089304221434691960639774145980198638983048163217231106887212899054085590421957559035067/79428150490344880267883758731574442265014680072756324811260019370641781560388618630195067439330263312593087579365153178795348909517950349653816438082622193664},{309536265866976899663623994665428982641207690557576299672440295613003285208523240309329238237547585614753662067982737169901348652685876345114037269243477651741/334302720114552064945721533555864413095980053767495573211770693497216881016560427393885451626761500162502122441174811292085882708801701647120158924317415964672},{279791241460613629501024616106661357591552930995347344993137292675661773807277146307651828787747438967845042767653990239540667229666663259049669109205772591973/213897092592087594551762882570600547837467477368538455616234356914316832915540162686712630019120274060690636135130999904499855642880595393101509989394576048128},{1435586971764148234075221549793387491156638148799313628547489609931573798560650272550476291924570803254922734860494262533952822812371931677926591779262954279/1435451821113617464432797557761053808501063361395494096662200542908802877947707852128998074718275634680719249291094275660779384604616898375625898858517626880},{1881458010446882635145328542369820178013412913488898956271341833983789996572943496967029829564673017556757106587972658653374141303045302288468754365088665234065/2075666184322059317697020020996358497539821658555912208139233603074664804270579242622987853317985598384995012612306851564740908079882182693964362043961097846784},{1643854550837336110010717143332251793828009753578961806180526609767994512114836918114663354564369859785468263169496866368782030410518765866683270070415581617635/2822837325989266899752273691283042078653009046885309591131962003465057638338830296177861244381092613731484255935670737614430035316809001224609721240802224504832},{16247426827269435167609712050448181294675457781919534721335830745624327932556293973740344214329569855644455475693150228377567807030269208678727476390817700343235/20788358689521495476944181903390127722740983421579544164056696277693564311242993702234508499595547189691625782854111164189518150908223136942115952672193439072256},{16504682286714967325394117885215488865089120775599191901583357240181520992102638142128741684619198675298560156042000277846287614831931097743174634777492582189/7350825139986357148524958695800135213879894654708737210851742901078092303774880368580204997325288684314566182268269689722545302103636619335103317086822203392},{307941747432275723050634874794676549104748674307108773884870538467367959393869044490838235149689385650871587747647709600408534281465438707816323831790191098195/360722528758886885812189559733895533647703248434486336929950088501374240336911672768053626943675529458401294172656015299750334142233488196527009104945842487296},{68124429077157749493295222115330516309212786000090563707586280582970994422617603604915522997074778799733590845582734535816886106054932130686183774308770400955/88460673422740315753155078991818598444801260820690735633954744106244416603960842922526522691447019955287406543560358590135264482271625545033321310473313845248},{10400393546639464360121138556313419446609230523736479983378837884557825945080286155859907085490836322377797372595439903826156159888919768745034320840189356399097/12500769972894379839956523086283596473980158241032983554813388195786222775757237671822445865288444012368725368642898074521851996337495609893106238231580659154944},{1792681035056068864006397704314101354344450467213558929410027659558266273295067540749960576639566092117038855599821964838733809277863032491097891538433316499505/1848738410706603193581024830253533229973597572700415159314288816743716208799769707946565313613050249534751263876903756982703511051061746534362665976029996122112},{4585821043522187829276474592082390014199642961196474461697359944711762459295902891777452797229608320066937220690057555853010688219343492822701309954794546431/4911842850454001198102638728373326410366851155117463837216966710738157944759010008270558295272528933001288906266706881488141314805589010628644688168451309568},{7805434265191538243609034972658103870814153251674938775617859937350131946997238576569020486574766933143414737470667707791941573908565857521870411919371245891/215305349137688711723371766181648036369128088399511395060358721000907873017037682283353589996979893272648942218223588754663948500291594103993899184933568512},{4417760269020822182944098857278240899184349787086880552878550912998655170975175788841426241083489082752487790289197037760308270907029734312739746500728043057/5355965893682393977408436364864799115155636445716750909389985647472059310624228329301592475185350690664110774332592126808920574354201239872720638454752346112},{1638364847329456088774592527243014719175095845302871059997670008578426076767362812026809361473345441928836595333174068797282808670865244667648327076574839453/1491850336461567224134992265638108124961898729085533202526932727850909182871307425332749643483263251249383607655402603265204816944921538343810694392187977728},{90737030651248401624295499335631690468603383858031229485529465881119444204647617445735292470771373039029039969875382441669471633142627623105351856172678624787/70098471703739419742971288869644035391848216378135212626699461408880725076701714060061860421354363686723949043970303126562161863045265088132477162280819097600},{28902291541364213249047586906480090391744676201406055347784496689270869889430545902450214718861771679208391204648982427743928947167775225602571767273730614085/37309152565896367968875745696138268836631898982668651811721610930795421576544815636628063082871521466883201417060404673656637521422067947253720731827111460864},{32901975703769472821326955572435855280154001962711035237178579034582633326184682034137642920158524397783955808342678713308652028874365394544842170721878581357/29077059765140031623310693964465010194603116368725197326713365230344469039152236044958531838014326817923052610866574725326338148870686874029442574153890660352},{15730310029786541906263598670112389352200678342159253947175758867287553370691134840889688549266193413233880280403963194500221634701885077566093940061125105755547/10190503863105649163159080994998010924106284423127630011079393442855527354917316252276521379027715103867180567893680202375245756206109071101473988059979319345152},{34441153967547864913740457748188773486508300224015471392817137046979179807514805676270059805561985194761333831505466634786910317500290599690882225175865385973/21464551831120751100816324893328033571347330142427104714186694964131048216528218722491779335144707328528805379398450770969618791730830966852385466595498000384},{4607587087900698177498258353724889958411847152075327253885192941743964803287273845629207553932374889449247150909891226792110540600300830690328780330663192537/2553380996888339391366585267079325187669733594802655867576174431384684033675546918317139841248881189162867603817805526949655969300615144257884905133590970368},{7004727685410869441641881028835736922937600161642866933380662765989413930240084329240495316839986068773416344604847932428544034284235875465713911675367709344239/8350246522698163247738422004972022963499620621297667130205784105324331532785451819361490978693320873540509451528409396101464055506962313440157763437830846021632},{81395180098825905221240903900077792557920635192037471721273869507082897646315771679252300061514287232551513514773851785624480720898072497143214819483638736537/49882784358189250882239880125668617955721287268502590152831494439726073709111654352942039478772994363601254052465758836654012702852463778553952111541467742208},{144396775184917866679976430645765328999682133978600507148868668144686598317633029875057430570567118866547663261046987670502723072064390927705234746579375572269/100297580823749999695977764444786972735485213525524497151681515835148137983812875264059817578235442638004903222063111382209038433206775152353462807132985360384},{118215075886388361040089353422090758207005086602434410287695621714832965167995822421917246715543960792158574229481781386626582290003927122607976928184925690101/193064983333791601079064686792266512873668829018341625932136637454770607417440780565969000534819941130252404683555192801459858541498337605431486634915295920128},{558823289564358137855768457248878297192873452806852094496799445406410486866485539012517092376747159468583595206718328013614883381791558495203237975199343216483/663621349694431454930312723330006301702410920324596949941616125722750133370190092680791712383554222607829785955217963662129347017488800748226539806754582036480},{996850218844396608989673344140466548592417098603578211269085242061056571478455019992568532492830581174155427084460763715738162092146490574412400124761694686153/651871900454167463025232268253348211965138281075187417873759941924578594009951845158302427430389834513719581975647226895675857346734715999823969310977321598976},{1718752624472675576699039346877138444467173861685277621050599459862051295761694814968900257348637032923478642795574285284323090358849000714253821005533746335/1966271295076078349294806999521777249042106554334371360054667739368959653576388609197507173296844623549163917588028855570277650029287885399395795517982113792},{364778282472339757044170013872984564596939050892047333340568140097433706257433431270943597053563720569091625103011110240660620454048854825606475530524540867288799/157433394071276351300473495363204049764593406751070555497874996026614862373185764505836452894470526462903067014765420932785362302112690223989182634840736831373312},{756426644625233849603162500111374670629412728216144292104867281623509729556138641010507121059525506021225270381016038786649638708094937479533079168650746700959/422008979131723848655415858721593758457684582376250211550176446249100762609207745667947466574436940756388430257205988786706403351835321335195328022311082131456},{609811725389330922816037300784828316939513635686015515985655560345876863014807380940135887039998491160474721183994479823553133329286059987198855217732632493045/84708938798994678643537736820277237673451677046893769730371237818074988102224621436492657646747468260212083714953296189090667090750539580248467885487829286912},{5573956890367117227571609823370650194457210068691351261054724019864262524498037621505553479146855661437601933339371050332740565588121257934110486283891984490795/1651507018581543694208107160649529148047401268147332040886804628235332406192957997466364433200010564361100882536575254008929174202443593898763779450300924952576},{138744669250295956254171077271966641950323820618523045599245194045644243427249960932702283700768459393803721998278455433999884666073844741624835583052318159617/142613692581282141076158994518816457867583291194276436511121725801204372656295137188330171651991281844680773413557086687007761124436736506988979750076556509184},{1422506273934384407787871766157383181511931981834421999132241696787144445810747451070000034060998922332937487307231107035067187066978039820199936355142383232101/804546633755151285294326249068302586390737604632578790863747672837880879921069634881981516712526168301506896393909483260020461909653100034745944002958238679040},{9361189030960472773997050817348357925673179906886962252965111309888564367466954871018839964260274213813155208560616951580324116631434507592111918438017421457/8823048705898899047095053347656787740129040926954690706132991695630878479545393333351347340182890766150033479612159884626230279727963401224489494435925590016},{83359925621283819704047484227186045026785618560210319896838690074403391466792632437838468224798462065626562453322679165526486777983936655554612059220929014505/1955832735786632994918650676822595207871739556017175867149720928122653286550637626925572638469307694427597224259805724483941535017945178152118583460114726912},{735105549826726066308962162043138572729183316602043526291436606639635027636393424741768517476282371169589857640286780561270007208399006731166145598582247263/956185309052376257127713040017708835059462753580520486749337309367318122494009091849121778141273469237890637472212691580371224531875436251782414790209568768},{30973981997015609571572921190560046641248355717355543421620984417307318685495276172441196494045750224539793520811413699943080391511528935794997262400235422647/939390948788390529964579709109131207239873425913593117200005653943037045573976068242667239495267665443590723306495127774325685939699301699928803867776516096},{1226869431830213872675645354720603386630808492148491125655553108535348640637646046317018505948884148383118504241966062993594456923387235226637876295583526712435/2396133041401334203831578293589505968805049251609851384890743056762474739161883139990715254331813815533495080590825414121476581338172500682746750558601547350016},{74399351907005688161061190294414801939309184581933041585909646079879988662666907378966460330872521382027102253440510920668805785276089851282453012929436681871/76808984876793492146992350062304852709144633651575274899197386542213144078189706325906393370174461718388441236971653487630571781817459762467719574326781935616},{47235218132158257061843239211700488298539946069254049081825368601869933248212241654769165683813389700425561107682841299049305450190045963765452728976631001401/90468442228907693792263759479494044252533476541252620584816399980212910699795284841559590746649623474084494692032577480329840002858307622466909890849959575552},{4734367882199815579615698381733645560236593064120800616382678169984165452499570233292527998201587242391744326112043633408149004440935768510941847920494599601/4447478116155472662973294501931547874473760608820064171902517613321019654528873939410229846289535625305547003425623346941938924088267458846449839425019969536},{6733966796657532840303301751710401051779972695486822016963418338713716091998429286544555201365070139868546402670208008828258004833120419269786809675128944209/6752391500744370630702840936106814621098551454078745924517946690407003288833270793292935314437174323932804913602620759913755339813650923482448437246955618304},{35781264644512661329700358014778985783732358427396276786770872587675474040013951165322705706902805571765351039718874979543704007731209472196553179188878477597/24078600003541189220041815958609505963505180913092335149340033734724278233412238112629365749879389434977023523962484151216185286280919564375805997338530414592},{28674439033410913164841175404661893630276576910718658481537865244065045504904238024489684633580128906427152072064488357848946462309569327716951777933882061127/40853657401347095917621907822986713827566212715375028190441271933824447988799113394991134047583956179327510811537026448617764601603403613732150887750708494336},{128815614812198206097232474741641122360823756274037431037885201126936809807133622925085052036152724549079686264884179008458493941197985624774247434257313711015/103707087522192105432961480399359160642801817444602292055724857047877101418373980801569893287738691943839918077410325341188735291160594781920656699850758291456},{4551626905425325741785694154674150356205004049305482573024650951661916223702240369684631063923988503198372650363720624536261449345708186654224167874447957519/7143934155812933583229678227562595164802612858801535916260248677192234780440517117994498657229458595034196426040602044024697525540480523036301132424170962944},{461724243569012294032886516699537472892775038539621885332722367642009611438838350792061746881752528081579910725324976070519152974737474386444832628299458671863/774020893285924067874335399840395463512502257669329956436819378152012709816858570773174192713870834784182368533303068496284779699780631145451760608435232571392},{881356258944297390665528364627560398499520321709965532649766154638391979935008966425443257256346614643646081735995220052264061032700187196511819873039288058225/1057751600655378128598375671963986003687518748782163939868494990768510311051705252931717004754755499934644430670800543761659089885566007114670065016945116184576},{7107646643086150407551301752117244820033589563858684387626290888278559899505673634659433346666585853168940600139145125279838603372586504259392850848115999932355/9090233566496666516189545928873670847444093213546334407724832778321828852623837763168972639433357019951107555292094295961087607984139794052725294832816360521728},{6494153678784412051697182780720352423179758698176763209092272082772110751393577819200199132556213684403341524004223351307948490493467983024538215904641181523/729600230642444991504155263112993088501207114103935564243083496533886334331958713827214311259884911409460125964555638004210409952523088114020186277629919232},{3624606847771623601454024756726367145867980253206445299312542685757353864962222690051435257042089713778182384268529872665401370064783163802547155103588182663/5684593178290097665637201144462480491056730920197304960566060213359771261721592217876870823354622163352093267359873416611240348392420729805416361402277298176},{1440390458878153942395006144042817965789337223230007879258937718319268458963511562063680751140914575641150455692186184489196905485367881477706566832361222143675/251559100583847747602527619854084643432456848620964822874990067546897224003235880195568330781547079495785830477050734158639458629564492699062614479608643321856},{43019820627818738057110837451174337690447978019474193452146264064501686139749738413636314276157381321074694973785347933323602368119692926357149595697357436386805/1927359320817336780245929094858695119859033043741706418396416399463978419237704796141184188491295267605259930351901866446247900688038203418449538813700952031232},{3000440927347722417694427847646422244373829916822356337300437377366679085178273143193506122097003027380264771621366044302166644836170710654702713909521685419/784507633521241233035525603143824433407528715884836757725724230353793123468037809287800720778451184432766604440242115170887013821928878722472899167127601152},{42338125008278203084616475395131996631477436210171401574155831662933388824213787090580108103870193848576205406339074282275512350192349009930491799569915765507/29449785549155314788493084882046214171208010117395164027752644770806119134402442911921176929785557932015849225000009451943142519106016508641903150457014976512},{48568055408914833369202620777274471393275813059808569276318001728815402260583279429593222185117703449004755854316029165214164813325845562058413858630254005257/33745958135156302172735130019506696643384065953637292360862196273865497811388120053064206532781788030966219616409906512414993650858200194713377903762401656832},{1982729486384224078324287376181008352775171308986086520034325818181877439004920915050037006193823162380504689735114611322654705039114760827244626092623543317797/1835369488112603418197600402892751552614049587646212612177919693645155444530191203238636985752225196356616978302536518879241390223670570027793004443429658165248},{47778228904269559684306402563267681341521502264733044782294793620855231784269310186712846680570383707252423956572809287530140696410205033084833813105742156379/74341762746927652564877387878098557320242253578196883593745626993825826916188692767917348416869594692762141140610290070475666885356712786433227012311869292544},{119616613455747883587088833724167744931493920477461271502881910161548945859551996417344488184660737783966632973794447825730584059804384443374307541720141524913/20298013192210761543639898795222214143105473393232309598319305655140991138984006005314837294285018866239175034731009952583910176856675157422587758514626625536},{45458688259634425995021056281691478049447889844635472852846613589011934976151865841615223511149039374292665618284629894508328977696794329321420196734687495125/79698713725584670643124141852842425899249566473349779321448366884470447278648278386067654955645608749738649026022953920083176808913530396072950186729503981568},{132632489552285961950498525207148437912557797563669134185353052648795937663803350215014636109923650770814743610337031692212429981902797427781497766068695107233/41664091943209627386797918942228373422849674970365810159158355116672193190843322838359336344553536932331329042054026044207894779061641209587142775245321011200},{305174615867138815463335623421516457604654076397081340955911095856048724765588191553438789077260073059779990255609545351422653565513993911801361729967903201833/374920794374579836658401914551606223025973864948009412343991095667619024860705899912611774331296140371551554166548411350468446388486747498224433154977561575424},{32889976478170519870273099468774357933143901589139854934293912606357989659663038938813331757141620846417612676478911314494856864907580883866429411170003081923/26415588589566274424578080532399267392698527429756287100453527771792138612104528316575591285118476770120946078544660983661427185719968773876626171216443473920},{42751629853900754633676570490741295205770327254394862571820079752749476222175132712734000681334658001233839300424522374249586340883121454973797311716476190935/27203088018247421492998775102814679574351872796579521685989481379276187455177543296846211936944085715944841870797333718982332525904628682769701504308329054208},{93500596359726993679984058433298033079287588860497058975972632349384775551783136052417764399449335149978947370974436895252942446385027844205048095516758968205/57180109097823121482429582441269080373882195692880450401900819634946991098586173964019256813303242493395799971573459659620555169552790398212407681748595376128},{288260433661183609796946733297211986584761394837247815722116244293809747479613223743959215324057010767494602481907847054680541838598298014541473171215300453113/339331643115575463879034280999754387848353372468446988279733353389268311472090084861742722856887494470610887867014926750198801129308644875958043033370177306624},{887511469381393334852591572572977170744433077993682167359854920225518895207716583354225715869299112147864718333143228104876207964070290339130349974212330521/1058537211071981097206964910284309710738844396892979610297465944937598908914867673532613198308665901060842919843185201493136627298124280156903848859414822912},{21303214528550961992790166093031762677530427652222960426285828158003897699492085973325325446215153036445470478238452193919359445089898361709524014346264940599/17799073595592707190059310804365257115858298150439487625415167824883173016233778071260349632795506105206651355275502922531155107769913496019927653237578006528},{1875299415527033438095644969170314011561422752005284839099707901770394386453445066529688896054246047675332488014881598036249348325324445930388701920691299381/2186747718459956421712512483039228369005909319228702118847299442689989532551755270500036116649830925490036902489513545873663641158014621094085410017094139904},{55572011897189750465695171672524091720902119004509823815767069345207636791242255086960876290264804671635223861056519900544748410900074660573615376010608808069/39854443942312840065342842291506510277876330721622194823394924286820491894815768098732708397916692962042508710346469695912707389412765429946633458716239724544},{20287582755582424972899318993034737116352858613122729392834112558117782649919456116841699329648588260605339903863396590590729492446169558110245927256243952853/4482111112634784605924760476826871642029068979769845906008058804089573858776141825833762421201973404276541462538249844427798860438526907949053490931376848896},{285935716034697633663824232733898865433648008259199399473253327514438048856911477564389482226086900316968666756439937773509396856033429742279833542472539695373/476144853427114867646598018885511417336408875317491391251839048828713231142421099862554998334697003541555308320636855752569359149681270858605533897210030194688},{1893209866582052313032486954108946189139512195690170567818520560509357615244408701356912124624607545730542653948111162962451038385355917821654986260949210288815/1203457902807825739829918466788161629532672151998419864991980255621510807342366551776268406514435030691092734687970762988797479598186246898410025526116837490688},{145322899259872025392389255732967661061662222254395968101370771513255330198120782969215670591674917415284043032228672338953781712966792587440988568893444089355/134889889644501617082903149361018307366173143874722508872508217948222400864624060478400774065729469573859908793329620609495012241115131641632299483484529885184},{130411854257561553687232876748337349539246105589482388945555539286425184007112426820202914216905873230012752086882742183218162905885783752004271362830433031569/131178940606784269137369227369970043836439907468768087051199957961371298275219047175964392880288808623932839441345955587876435508539011860772102696567641735168},{5780700693167586135735381480593139724763620070241625624432430282918984111843610852953921467865974529426214056120932098609287316895836251297306587177755346511841/6128378364331742324665864035212071611330449343563491463057390718605253290885338756769182154699414288548438044132097951052359179710808695731718416734311516471296},{619636464371202752968071191824860996393388740411240225816963036220371828605895067377992497047302085128759918433571818888604841847910727367799588357050156599459/1206070765070709029292320787569819835247155707087067450855112360465704661057609126078259656841729964241939226904522236066448799225890019153104652611038641913856},{1048324555961302774424113513218483149578267293565747354099639334294998126349863465249015315195686613633205248346009923629436314134316653832708731401568846165/898716186861789479488862894454162164501040966417417195419041558348538900185634140085566090407995448768221536247188722336449966859053225419321878747939864576},{161806317355645624217255074389974368831295228831140642338016200535052619365418884543583728776994531832083956992830427872692309212513535261503794384441089819261/275361710627698232188544692431561057354559930669235313743120001244113780015917255340190667107942637588656063183333087193392006953495768014007043786231246749696},{85322933663448563812503313001519202419186625801884836026044507605895056646689098323555616798061310981638325350104578207952214617421044797038411056016886812517/85915050800515341446812571386447457838411670655544525200830974810726632800529299986177723955451628528761468133597761798673053188925500114479018267711060836352},{30325362043491254190388699747113957582767781719502380469859408259762499642956392466593588972704913424692770816706014812871738481028412989964479197439843094475/42860116088870742776647560443373372863891692249217174724315435091811619945292305535670691402971008253311629742881647455207005562199530793352256776479761760256},{3958512597067875845598407708059417344004309176475973899680871385045447107778127560174720164423131199130685496539259870730641738543912758122111866070167349217/2748581127019827294077169731645625074370434914596907262409570579240112059335456711452180433926652637390552985736590984614274979384059816140668464603898839040},{22310621042159460540268871117683178503449390449488150641138907101673246906968636260278962170681997758401505588142053104885595348238205407344299472653492288321/10859999204361253909354148181600791960787668328109983351453855952742222875903872929750457505340907367654206659792943560087434319445596402836824872546417508352},{537425922626290528028331233536490889267352347772454666946119288447222084601715100545385824066825964874750433877619208797239001251772027503360909602233542969649/692358259820443521888467153759987671858752712435642108651181158764115219576855343498751999102624041106774504850480737151584630892530177189632423107664092856320},{7879612638352894482035290000810164706295055835919398370071994937082732963078983318866946906232325170815166025083891485786656989047338843183420736169370263813/15458233803384335839136234469313859367536963731518485907182967147214856513243929396265764738963105077773056862801332056657895791336862684633126872861078716416},{86076131977616813434541134493935565927882045907835407374458776888015435920383883051989257098314078407853063314163853632885038951717885335733982400852787025537/84289865421257117529199363036129108506194660638205624069423494721455860880279512521917682264075465001769721721274778123953303550138711156514283231702232858624},{33833708473135301067075674910595217476601599122182123810995924399218662685037988952845167601897086075696804739491791204568604120092032600075005451938923095479/7476785256715534298233217998647069276431138640611558252777554353423351868071896576601226520087501256466888106132734399422248209727206829247414963294380752896},{63824960730919395133650015646856390577412403354143984861441508974073899761718467773763635939198758352945457812545765541373734436204799817188250622987995842941/38956643814842204783625044767877498377035405683330675260313776728358652966326389595371776872694956537303211494266374798248683529183770143037249634207019302912},{86950094672152855301029308331787957908230681320876119010288595370686698323563679543734630863002912942687217389883571310035935922236395170552893094453346894093/29757044840923275434893522854560350412577187105248565495647219129682363888303802984304174568947183997683555805627579025394623383382961094244060551502899445760},{97973431111738560988936606524051483357238431832591816121451458160172761473654938293140875194282516728448654134834653009856990331414894748043837683941877051/45012127858275142526335489579093883957253650759860607938259033474739050364136068228019500009603934570757532664590613643188118161405323171977624270640513024},{21476004404951916488826711242948022553233516627808753878036593986808514150053712873227577531334951279632597457244673342836598171650098925765467185839198143065965/24749049941705252338208606792687823246728697465102338918908831866357351572867800041274062365363897021658181343782340582393039609653870668207260913690036007862272},{7296027108651470660076174495452694009381416147448498555901711722839036664310900874538709302763958498189538139325656747982925703818054703885602141558192387938459/8784905690000902690855621364776130434702952624428145328568303660549134480923044544751331673575888892151635750032519422394603075846557396132117985151450764804096},{2601721472345313144406084675749337012773034877759531306398907453290937014841413279493287873990009643512308446052163253382377087485976532650294099633081715373699/2902809123629903195960985202026768968109275459378798056772654795325750046016608802256044626079660795134367643654370580096734954189455674445188103041993247031296},{1875387409008447175501438629152305184969572634574542367398615383019158444427117868710410753497821172933647516161573008761103273676878514318296593151658522070767/3047085909627518900151723975730060795458405474407955945847872108699474156289644361345785903465526430984723535406070446707989020064250911158151673159666221711360},{891386174458080478664935127107793671096695445380766027557168257300720232656027854236898924439275187156660982103070247518135591637577436902536975188607450895023/1701260079272280999813707544044895319032812431886911568500336223041743719643685742920731440754695306163115006489732032439168957943420156223861669755871303827456},{79030011599949732368031047152508213785303258516493028899577794346533240173171480589713794176678153604500827402385183267076339808890576808056582316262168335467/148273124805758269281307364918284932037233756014337476421719291697912047071305822786656222711131812613999608956518438872823149095374117039625772465717243805696},{55583224701346988893427121496463029547350706045664799796566719148903453297915541162226509219618372353194571775442607822020958160963724086809823550952024046059/53030887807189539082387317734900106175372705661290729599728206196045633794307614077687318073622779459892453108636886195571018413993080122022001275382169337856},{41529986687546545173303946643351812972826713538912304917159839025812827575984132486648451981205587359847558221539392902545102881182317523555992516623216642343/35331347188186295793561109540268115020464421847989508695213746057864820542244851179986441348772236453986455180070323177781369467529920878538634690810517389312},{16309558536398103380585392812599907544385926382576718532869002073216354049236867031957382356220087618743291779578273221964524321545065017659163748434341554249/10192842945912858635066586776331203686727738368337940078305834595956537855362411760394197126823457422136518437407636497139585036976494485704585578283796529152},{212139459693251495480230653459262135009838794348436055076283197244506454142341427876788440825820452901496440546086490900936415564418818471859911196434668219059/127484149021188355079674942528500922494179442884257476752127534137588034273727449653970837262399453488188489699139097679373357691507648409532229797224598470656},{2009371209291146940262955843457659000069227621611896457132590221581341688309714358505766271202665601100047817700346805905500224378358548929499297211872234299905/2105664013549223879016032395386978502880703947750221628329871387769342738713796151793753720120604066816763514405238200195891156559605784404358580190336507707392},{3233262105815326775806554206115367846325773349808769297010381037427113229577617763708073525297801204345963038424222451988532110742491650793372548100278111264259/214485449128492083115071763611095390889783063283809097506078579362359031904133368059686637180757293694937348774461469944824245860127084593527815250079714902016},{927013706428889352519773277418727463618289887033193397898298650899332078903265452923581484929348055491164205033125957418538794086812460174450976740645935665405/1133669267983383758817302651458860170373944882993299986721430062933196085193759196719217827621295525992696389325943282400732730331054684213280893757778445729792},{6252164881262309389101449782160897888736667896139559834255097315102793664887893403707102085086056292359564241293320177048126487006974073626316659545468469871/5921221112334423585663144622884965330964432767027269746101837485250860839939047856409681028989646818959544628635515760614450405907597139257074795066182074368},{15257433230658947024536330944867698444611718232902182938347225726217691686618167474767402605059413770537984976236932060891026478459060627244942302273529276199/25720736723000687162015441004956522516704364159493969523808177103827137219779629199483712245458782950321702539668285168123682145292087195826203903747338272768},{3059949365554665639022170945565418544378100626201241162089108694576192038325533324289132748297611611640644131319314776994101021710929057062006400240560427152881/5235409142827934930130979188650480788809988114424133824095557620034828279581649020759551827521225169683300353068214862651764595983208331662344515290745555386368},{5704760530639769179378292714351384725060237745403672203132912013371068673890888294732833357038537634011916386209336270585568108022846997360878098321474787505849/4973246877617584144649036908654461024226392567828187867515585204004129091044064267086745388624536995628859884823273377991615870413946894221825825874492065316864},{1426472031270753050859760835344127893242994201413873580039653303607669714401409911223331843333315063093995000267377512344976955527799988111234120217942344135/1838839720133304989435739234558721617840646159267427904313223740845880174628097494741546432637240708878118444553978073423769041435729191306826064137679798272},{4501171320507914226226323221422091226394120665614535758036646997032466284394120303768640502725909349411091363790422623486322802669508431748488770304554361951/5119859280849597565700211699002989413777592413505831850282544577249249600293041779586676311059890611570234730690860206315770668990041699495414993885576822784},{1878389934755535729381459000881436234055640683822750555029991492621308022608557434697738759512226935013876658652180562844417241391647231165333313150133758564219/1989672280804520683645691826100873225567760429889601277375637092182604999483321485164798535683319276457992893762400985448142153091880786595489915993087985844224},{276979432916908321023560495435706698565632321443454664375377664798566601670407509720303108141231685702841512933270978625140096629439946090525414226904623990233/382857825434701333562181615931830896394409748644645390464072043731126839854553513781920811791326424850096507825439312359681676696661851543963196436509244784640},{662079186781152399045911381903891786409680273655137162210078964855239393159126262551480033921524709487382382593978929648157681834812275431525640316391252421213/338864094843488353033149511174304824810413090522286733445383313330876234451163504471959802150986071273623980278062441319677252572847160043166901512870259851264},{362316933821673016939305505590464334216244832208936899666663440850204994591613595322584178320749941816501542955415708493478801651301245337010069254769957493729/566666140833707597014858492855701867725008021116297140843898408789704072914017900249897805453635613687624904589203345047308768780704268871485821353524941815808},{50679042805818432163802941372412787974874583833558879913800601175353297779608822296797184583837612617247862862108384318166635338400799632277566169460226661968713/17961407519073239038957512598838157880127980099524162400030143731597520815277796316851738477909977120790718334425233561771344570133555398214741963899108539236352},{926437372567263523526567213203980768150591609277286163695615955326539677363497342886746349656771868307717947694661592655914194729089058488059196407266338617641/1121260709423801176757583469122017263272774404410415753610495908429428484185582470842441123035069788364584115489228138956489683037118289265561641692628776910848},{920803654797698370606141962715439130726263262573950513572319793742595263369008174356108230096006718019659406224941531263615231532468963706995840693596965052355/486276235659598009772496791962331680761572753566410886488060190825548899644930217888213794227320085716583669608398257683562094360190614425012224519429251137536},{1343494544670838896557272001158248730325235445193964787610645053907586305145699092135145637844035359586614214147133765140432799798466515009066350919776071460389/1637781063330134219254470569044215818853271877103891130140732179365403531981191388142256053035070949803796976679615563227255963382518823034753160950366984994816},{65994492446765503749921925895653191206578395170146325640312370380457800840473153942567806646314178907203595024235790191736883211295586063747172788116847626253/18858938727635001155644528197623119227802412164825014438205819016735011955104400775061361255555134818289355456796677326447743259177138870602240492793256476672},{128383761392277743858600236876580399430707818547495408095783656407131672948241261645502005180468703797219450219779043104151950104867941349298085715686562576749/73958695041556405003158990063566458505076963718954884214071933206505429978713265494247745078086159902879823590668576558770275835256877400796991870892693782528},{403067025486969373148883634872964281815172902035781340091526965599916048016355581127112551809278502294805490719059562525608368513420263082953951182310715661035/169422995201239908080391995085372367223432257965291134239345014451012902231272774513956955294316252101742921052486532843148359855413817198936544392769588690944},{15151484095600815075023766989740304724197883337232960521893851321388581108774266811594081970227094949092481166339324111531226533400112880905479820797722388025/1138995306620997516953598511819942943808554909796074796758520633105950265514381967242611247542757453902892607503661864467774794220030781966005033871062597632},{1406997085196239971470891656875070173101894150006730846412398829743758936131714345058584186504905420739313342646237696326108108409069320139016380452233960365735/1512280700020923348713519843160196980126363188459747426269455275417573289750586822666977783856999940973648949048479564030625216559021317808603933207250996297728},{20214326516843308434887564978140752034659039533858271692816734942659098081806880060308687793404088238506279658461055945920926457183534802212002882798559095991/28742892962524406213336215772899733898866908066003957942377033322665713310155351772909645876875515146672285819699157871192071397457494381151840972162333671424},{11132564176672589643564493560245863220225404660910958661837188692896702101087172260243473400863677134287365457167755463734138742226631553854840177115978679711325/6245859472115283792809483511764956887386260902358659276137970698133012642810796451294104496095803348818340395891909080547490163449965248615514239220527268888576},{1782499038750248559750639695016615168744298093122310396503501642524047752219128032765896004043384498366651006979487056808470116807581895122088540280760573798835/1579648798069219423156774753186590499278771412681742611820839258849944063741114253286668525488354962340097854261874338694403918182953877823003018859255628300288},{221696512080344074403224559646658877250107541370460019249206156732517061306962379604490171046303754009218838369166697136354868045292594347251129454661229217247/138179873907962784321926500566130709382201677378126568509033770589552193144330496552952971406174139486522977127947985944491493124545069135060434752863620562944},{38161810463487363594551420308342704094033621341321340268620674442627871562161048522016523380853435146429819842771310892990148466177406696876540901387047315839135/23236032924937747000265600072914096149020454617912056163697905979195079651362360476665979734078212246045563974429527636104321361222879526058615261838684876439552},{13982245080984690342970156073763341375107055023187490609182760478496433078439433400750336637643987965302640352795310237069209958398196292178496860176817385747/7952084851746560739687891258545044573715412182866319363439896252984390859515156979093956477604471052654598217564013143810859705703004398806330705605589729280},{6520810298189796188055685503832367866897950867531961901970700732488053833779110898575559932746506201218896434077213711146641466821558171091369206225085934459887/3439676440295967646360739809849093977324317306322595680312965693739715470597810942058156086000482085366518589963374716865580819319788485925988837166369583661056},{7068553612570632714998690002036896498995412453668044562808693537955841083575279573960834780707607886286786168161107224066749266688595255224795488646493096009753/9951134735608500730950523539957289892594025680668388342495965091378410737653289019265573602263123035959137693304385627115627255196968067965971349626240070844416},{25507541617003970701513900546370251699916475333381145614976187924212297054137557323060038879112093751073581257113672142991339183859066315861194786248612464169/34295640136628698867586526459519679030262017060178165544805095665807098649021075021550343545326589795492015184461035307342467788186997966626140543059627081728},{677499393804113950705322996328497354353015332783791050995430342955985903341266902749483819967521393570363506832172976316931815692619887617510380288171054463821/160059348843472168174955129541502784959756835739103266003087999654827061372887026790499988805636180947428135883124657623804705295464521643522021089129512566784},{5285836649028406566992629955612100356454086933558700730707597394256251327333795527241515984815649500601702389472755932130879202046030620907734286954284324329/2450518022616027994826193989429809659652337846495950201390512558173514769946520734121298385195895793383214694040175049256990600988947082665210147067198963712},{346138803236547620055850234426629367034253286426979887974880621732459822920860334970826246818600459486509921837237586477044430379400956078098829562328618577665/541940555941798404887045950223767548928971110358056694702755019886904087710648852246719966824401055297851505609071642109446829277597464797463794246621356621824},{29155736379732371811904092873919628041078790304770273811979701235641151441864334599737521476529716658438428706951682915231996715444028341503652695066977289506661/14775775188165256392164161843708713480051825795181906532394209057134497481461130865418341366566497449364612810671127809748575394771097989027130923254385672716288},{438104143878256356530376204707310790006979107568809289266892199925258480540018082000601217253953830171774050545244131627805808426545893453154003076600051172003/740757409600109653211466485300237331832296775437631461033497969234972099923569716873560975434693727355998178571823342700109412144854106681091865137013401321472},{26237399386200295555619517546716105003165871463361212998385679086383499641950635352553938136325645487880393400088372534192109686281938642272524991114312132917/19572490226020469486630075206822072132115291017415754712106824597516194165456302399001979300004334185002536039393336323900245063105631890406496664253460643840},{3534405241427603236859956262465972177335636342021934288885043584374360019449274736766919365565697318968986662694991646650212626097835021429664986740445669618953/1799346972628996494632864923544698949496198240121756306876164823094794924597044502152291217922431230578938930568555410624230702871724511954208762052027768373248},{82371719142727307783328029788592481275662406166684841134993186752969019426522921860183992584004835572309551483175561350813059814512678144959709859821744202593/95042934743615604947516488915311054411695995121140323352413665572042718340640632730237654375867505347347967700083287874785238057635572925734990547161578471424},{104070736015422686109067192814686589183071711084425011990908069196280509144484688893566728828245176554242601799005394970136223838254910140956935990590835015625/150895482665423502631773402457026262350705185604261791709911792687311531640405269208569110748259449205944405914338104752614415233582853400673669125934218739712},{1594269244690578090527620922928924546549690146145206074842706638351548149539924243381778261341080208963681447821097157945636290640198559160299152445709145235/876584059002625742007940636369801386616380270757364958411403624392207043109414095495314342675378580745133956633933280542093929704842062427344553710081540096},{1826496718269897410813115184863322597239053424751830368681628454043392817605358347476056821588450904778121566750014140470785213702409150432106350516635319200809/3082467943741816366755475703146114403918171917245971598150717115818201420069963463774701429816418745362921913172139711868508677006316511457718471089827927293952},{6156693252826663001996798596524326510924258656519795772686741331666696121207295440583779517766830722718425672776680745569378780270423673973205801819748372391/3621076485952097797595843832468839073065137454531717159667636044377946951602679281453246186255612691026884390279361164804307752357762074537961743202044084224},{12729603081311268675153811485477047621815430473771460692931249135150654463617465597380212432468343009087675829409123501347265593008517043542661341478291224021957/11976892060039056442690513419568598987720198795211798156300945422516614189947997887596845173402842993195989793826502466653536511082601031269574958840668877750272},{59728342334418321656027194866525478377160859435845488888071372429211888255904137694500932872102751248142862613546604203434538124077827400010065317556725252161/20981297200028182201703781190148058234111611518563206362969774115630458455906114366647412896320872751906534675743280069096488941684302003271273140485890244608},{8790412873459251083582214948839150177268757847519769978617523022839778408165580614349431278628020670024637577599260740688912953976178464665442033049813795054759/2082433985974524968816145595554805049767762117862193027758246161224221649234546395222662175008650418482322562797596811994467444537037044463098859690338386182144},{385881623145670179266786313286134862158672362659561548993484752333511036048048638190946053580184076670779523773854512314810604426940597700706352044866894830583/511386334962461482401490327463179120172585251941377700099599075062571459655918096857258343811517800956457838698196782068476011416023882010870386612163611459584},{1593753570878771303699626344056480020051482875408732410453172048475176055630586768570569677854937720551504104629552321874232706184759875333606733639138917249/1942295269099577686241527114027636405743157895487492057788600242599930267064334929652874121240050969344815580076690161885636950215428864413240681970368774144},{65643990820606295708395853271112892397131309818663972852894035255430071925977140538267594185588033842013264496391204820908586613108803112944704227169235233837/21754100735018857991780694397770515001326263897869312844496811503789460613053983871166056281305101541446239096996379827615812305692954081114464859740834889728},{27017294842451204183567369357376574378710998365840808941780153500755900249026333455530392571347444498969584732771787726200522718402296182741571737407415789641/26413130211840341667584906938798188088458283501030098996226593240184189445942279706432670561040874209834923432712753456965916437874362545496603000488050294784},{124834746432637220932476222093589545739550234144343148720531117200165950779841054508459036411646906557213403770045141747250784889719934806274356985386301168733/164190792908040051550239611520897208354540104048363861182664198935655858376408550201934637776402890799817252499450774916101537026339579945987518713027119546368},{550191447311846226633953318910974438322669273282478043037669394613356164722195050022964341421149332923707858738379631067592551464636664834581523835697607739977/146365495103715897547873249030460766104724129393080258532789239722100583443788885164178453816076828980681558709749455539446236205731472545216350391211283972096},{3118997557714908381863247795822675597616126502042913307733784230748580790045780736910439182051907443890976370949382771118816927880567935178134610202978771449533/3215827381680509621742761119209722982552594701849951194582176164282858475825492253924299083605502555508774446974948055391732747865276942928946212131903446712320},{45219748952701454549822792098416310051905736861952508117882960319494234708677159273908906567688386818462144353299747296887686462868372352078341188289040021295/47414403127515301869981985998731739675530966561823788707304095790172473401461702045151483574164564294789074762657129609339549133257406607707349560955894759424},{30396874950262435077675967739377478836088575318210859170405264056900943596290720935846589985795323259833078992091960095485967741433166534000515734614222614877/32227750910008662638441846627473268496225485789724159405691119112029249582242079471676341670492691994667949501381397021798133058097261751203156217601135214592},{152724426137994670785950747697400802090536893426887357666856729086078910611993829884747027873106786981972639995309336195591633688940772356643667337143775401261/93742635317115471295872895380277763111560434596268604468707371693152434947771624259216152083749453830907049925851465152861732215796124125090815276440247861248},{38986953084834145311432062827276270318518986495189358098896116915522140330225997291972288078925374846068109702069747145294289213670908847560190906436637047913/50264441355403988835379101660884789466369094715185270742094118179689864295630912109181462811100606405611810820111796079870319110096346386411537717503660654592},{2006492603240092336387440582460106665208566154853558949821679366419688200593924684401281775398751057079874187873731447204091271669078289099254070504438699457323/2518551308288110300980775148372038797737932013070342288318914265648230929449621439038200042256895910071105956872697066779166443748192342256398013164517354110976},{6617622965821856279031699884639239449585181876721804631968582761073859847524774322591662504610274198739624444666926847619429423082785954750338381834083067751579/8751262533508173889043106102036331288194720876232017132613796792101131324132592739080483560084586747703952392977166357563699882637713363549142130419568268541952},{1471850971647541583367104465442057116441665673943525077905440631362094518344240092244748473041098040819259975367009810046960443952011224019134281546831531298693/801417362519164753137936572306166068988186236815136408465361900325157872511984442170417502985573196130761976314629751720323444451489722371839952646246703824896},{22673814458171484749261936695015778723141741981401692396700456253575207101921512715263648354098535543544574350541491981645722432121860293718240064097627550773/26690232777431391549876638686877118769460997777260129909516022660369638340158313131875947854123019281227963805988264539059523450485514682588526189247369052160},{2778487676002500432576745112749575041811080476769212749989820677404623016038551665927638322000656263273360648627813279244492854314743979721739643826557164234103/3340954675638428127467474489783313262515292611988945349392672296434657007330603709432194554869487268166190707742546683578612950035057565166292373127379669745664},{505716234368417748994482189071254529585612344526217387770587955996772032710611104184720299299090868083080027624961204890671512415204540027822920283697590387/32567408329067335428008849248786864318684917220374941756392980865132380167059651895088843306022274234294708746316655995201076261057860487294971955701088256},{1627534877573501596765921396206397148433483255593620372389628207619465284365747340427506509300073875977923199401854890636664807560522800668910528932313190337961/3131029504165467104345658473294754089229228308267564087890291819891811124278632754961151734399324498189114034760626613000943637885137589015573124024179133251584},{267892449136517215293482271538997331243494489620493593874791115663060019001903693876563814139413386838311485926224858553657255673382448410874906740630395588039/326412208598618849723144412885251311221341348430228733261173734687644892273674886781107900062667113564079151449756391811593413568722987279188403836367616868352},{23796868452070875957840549830405758151733905248285557021656451511753930099333227228313261988357794738840775715390784214540115453790182328395742672373841020511/29299168325096919613866827708044937606745881958834136983265760628683512727681776220027767493146269003032018686387879931204101775113738783903630627996568649728},{7060026201676262717231862337195827226354807534269152099808518048477323498147429642289210012458696211530682809173219755883563641947869207977847139679259476391/10449944402762772306835081947237939468287897247084964301599806751303732365595453421616627373615688242696964032652545895476465521908866049795841235452216475648},{16585712146685562362808658752662500772468975696749114305027960717773645573742669094132762998656398391996598382071129867630249467376565509701479989524835254438081/5631213159531691113315830666993658376182727566685790466628553532581804117577893375670722988620888952081003599928732547212504657982413075309105870788876830769152},{218885865639371636064879921521973765438310186793185266457337324570562070426854126255040816225229935123303289530911219208823870175496737429923585366965598098457/182266759119630135448969626663501884205862127249846168342650392521135200768238837002112618302637265670321977215333631688203830942428481191760896550178297741312},{7606435313221802049605889209263522790476400117624060925844204268108758295116758611465830624408941138352519646093719774671048215931310369878940960008450424649/67210272924567419794484822547892910179281515011460364546901466201677388703735095345708006757676663205501319205857460701747296566576432420934294812288352256},{11554345633964865171741062393323498560093279802241040818192080779501776222821277473764777316142433372573528560181860945976702525362437960118731413787208393797/4145726249549417249901714550184528367841485906169058409991931450291454429221210794621635367575976941257410316676257355483348905125858925279636492644203364352},{2089960558328880617939522526855478187327455882068620447010447310585935490022333861394675792841853992421632988564091931787786560802661677107618666228328263197321/3214379331252841063869303274089339744760582036045958856950291141789533146698462365144448422712284147182102724907644842181139952133369413306479599513230324531200},{5565870959032913498969117747886798283528365412030665118462078908274084722156291777257284795205474142449400770049056780918358220993579764304281769528535798315/7920757549738273605548238935434321848000408924499059975738425505512929744758689076254392163022190702179722453731009328101961548755626086915579213280909983744},{50936781184200010912003951369552123711177669473817294449468180438566463064492238759137642016995257746552339760783260988270958543009938445739603675517452491607/21988704232927489855516477718520684916601093591395228751621837152023358736342458157770986435911531569210720574595408613706981344873474904184189346615178297344},{71819169217287996952304521502047726427065917017183800014740401286804816720674887555935578978853954229269933858014740089690821602840343530619396530224626520326957/4326450963629978886244921582049089198471025588983562992382839503648681102096763087528169479962272587091604234943352644385568178431231135182754565240879157608448},{2964623044331826310464785054204430499754813253585167789133992254080606301022348720474539831480762748391684172261915807897059069317835140028144135568059732656513/5811409972178289727289109557732543887326192427922627563787847671842757880238917655316914456069199577253137509710311837421439388061931435751175527265883098447872},{136726969388677295766613631706786354901959655453422528360010464103468267785583734686824111989862540576082572293491674887967582609987186441171274227426942702431/228502982591811092259512164672285474092555957138824937203254151261267723603857515309235439707121013119089534099829961781233851059795116316397651312250223329280},{100602906209127141038724844743037345973980811932846513916169372368909529469409421269745374410386962033589828073452390217211578625561293797924904201485810819257935/84912221799480171136870230614572589258095920631425516017381550400047881208764998671726415564626540286500713146916874737285385661067651282036411311003911453671424},{30500566205914367073898295843511759358029592819483249822069622944449129112945530186793197308171649353047764109605025482138020376698415303750789396412288574915/31547666781814875215222068096888083494869572451353152853676012773667253757613612333171097738321611491079117459927414457934473147781096282662239643700335476736},{155879687828565083492813980725370991489576605082931762203537529831998020232173302558762454802722840172838782658574420607912178957707769615466225803917405120031/75950788739159179154423958122963119666024094615541105730873042822269289105181625629409320066079659012241096040915640274768762105198967317955115060846896611328},{1058685159807866501559355731050823976224173355467298419541387743653147705041272198038168041198442723883166755065483076146124335657266451570276773017907385773179/2073984544339698683837992731667859380752168494564667118267669088224036019863342807541784329630289534629030602159513083056735860726713631337232127421864560558080},{121289452997310254144016142052193329456805793782740796168930552448707956922561188278689934719163769315299714989223258062789820687124474037790370457640965009435/157963615882710329084715122734832734814370827286656738718231066285136357534495751191340248537030154543895006081075761924942584166404604362117205135593226371072},{27194630340981258583307458014359059962317106049011810792866054342403346918083338036922251670536488420103578729723415906182105883069527038737954827159793348385/22121483910896126542659360977183689357251384750950860289670560166286111439036423573125691913201572124961609881514299039688546217904905328232608954780087222272},{1687592810462032736150823337026035547157293256810575583220068240520129370725931652228465516034250433285850605956994398132634914114901482960695532083576135911591/2297061747530265009680940223288026960370317177018142944947050788329459536931629742774277061802122537245601586470150487919005327122491455869265210270571036672000},{30993923693963627324501184650570952270944945559811960651775988447249847279734914337042264091006160987845920574712397837955275493630886886286752707471609701571/8462499238027638276387816749561104235659079173937461183858150358847642750523362914940959399211178665987296004410572701022269937955627085081300194677892317184},{1180538644580668945114097997485672976549014838323426288404404961369792674699301728131604918650867981809319400275470371806687446372204132789944444602586109614507/1489803717319457776768073259815818046780822257593649576292487302984689628150463981913225121164239736319656070710906942895102029498272605185674977198676908179456},{256059443970574776168303234576681796811086113920466754312857464308546033561614672225416787203026116700535896116608101747079872393461905883597661277477723315/138654958159240250440241662489083743095550454578717086279909307337050257187587002283035243987377915377889230019552862545975373151032302708645053382392807424},{363135217669629234895632635077667474585337061367896062405096977252853403432374611197522346682895746932102064374961581060305271545680750841525484110680993914183/323209329778168881542711001041742615311563888278924137574546129052492489651131491353729357906087077959305902502536034860368912526199154418506198296470847750144},{540095084038812142112995765158824214670580161363806733970620359840094092328675981529758015871045488103369048144692437742045225224636909441026522967497796024639/484017679024652308595674888455776949076672781027222637779363868229575555043794935803872780625688953079240497490928111673600937426867084930149022417898648371200},{95076048344055085358777123849990566693596071239351971883175294633306277789438984059124330347824222694918956789529403796989364360201105577545609942031380529119/133449617998085676207834417180977621859176211461908097834330291960569795703203748803448289554397828432036823820237635542381560303198357061152384458258812764160},{147257188472732650737485846312763099705576722266691383731519840626721328288282834704859225198918344654184180887553029031615961394578233510700243840811054335765/63454871602718646240216424779760092760082692667227113397070856563834298212092276165368323764965813799620126544894150246886607608313610237323929607417702121472},{19601147721509452560380887487305266156647050153481055810360752972933949047538895112572314366666118030082197756162777405633551250086530730399616262696276519427/6889630673649759667845751404407621884561185877012788477944117868928326982504248291992615718161936324732425279154799081084255801694751346667887144923192885248},{2751588035138291151587831936485082433325564369480169886740508921344929635982314366278213291957271425935788991484789802247009508317458694434512737938170093833527/1195989043202677756792478581489858651994327327167717866883178025107707588067531830117302062763820571680843740541247477720858889795047402424174939306218064183296},{2861574617083203482224842379055569027140583917685275176123690342195675792105442618196806153458370699211210230511250553229371335842537925433838498046819669523/1188171121847237404955437189989871243158768686386598045456409173192310566165743979572842371762316287086307967802564536654825655974607128929875712079214673920},{1024909814435179028851213103752059517025844082623153234853990043609630188558608753676406511117099984217215851177217883193781305360190834442057412969508530324161/992376571673359932338263338074743339566487642704361103871296187503993582129691970260985670189530436612798337297192347046157123319665711055443064425771935203328},{49194657272811102594938920203843927572268837692401001235499647732136329771905661127950847381857566561228171816780428172918405945463356461121034574256940950853/38541781442385311354754772410580861569048529739477724722384849374677189182878126262677032017424493045447390583825389774077901291567516131921307198337954349056},{489518714132584419851086285364644813755280716261642493324569746016238464804948000096709631381280684221584215672185817126194627336542911071045687030645867747711/546055249775898346959017053453989691424991439650732981726583590180796726872668117328826718234427957831675718226386959825519390303318560522765991986841138495488},{118279716955604634301363096658272998128819950850465443608766902786763033752391101328611672442186072276328914404188968858378489225075999556658681666863155056285/203251719184941736336196351017713429960231448614916416096164682000240834406190944045848505483342424194275573917875115283156707384416068194542545648429842825216},{9034345078070890162320577027502891865118032976291594115887443762670948446533911667279525161618024458139640932900538381925624299477008993243655042862210469595/4772123292712598609369945724496160373670150595151606802926162173105204842265234063065425855435503670393056175552357880838147856376612933319539862478746812416},{14541949960035368677634581419117509691903307109321032060445250864342208775791677552872304528412858495497696047168547084425970848940588545041938960699278730601/24825964557775188831434884138129523272264641626702818961120084727803148823104069371597782313204564654939238861348064667415613755517310096077160260163094970368},{55139201419238672887105322788484914968442212429428863261435365198093026903902096843159721413041485277044743991131929440411380033685980007034597158656951334819/54796141732849725700562411261431726502867384904834116302501501322496399418146910161105494302860096873788497679627622236709023437138621335263960460957088481280},{2745994767179544176853541635575581151013220396266551454135921909288446212603599701564684261206498418317191513976901386465064505351020256343398642339777375793/1599391211372761367546428132966575679349999973464134988799452211414085570824509295941871541365561216536541050590044208230532672129006621796171449132185026560},{5380425628401288228753637141868702856642749145818520923106804859290607329516291345184388689370726292202439506424732903250305394065644174753966808196450300905/2191944806175829305969278655358464963591971894502414035649201557332113093979109432015399705521943643111960232483051220137268931595343663702900729850009485312},{70841134258244844586357705178336657417473348677968253742291008795098167187178512883108036455842044456868438805122268334427709747574889213059962461358537879283/55294903656272372803809027635066950732637390932853787761662054517560485589502584171973908642566115131091016539901666605655929193180772756461317010527814156288},{676446947332686376605694127806384667938625289665315676424751270972263977995369814016723587488757416354736761553137748539373204838104868084607121490063508383/837613438527202815040435562527370780761422778712349484560765392594607838294952091333747443213667838973586054308300901665280888423985498805153938936820662272},{139208691850141780335649227523240289507044371925765560604177651890803021750242822737555505166184346534866337315167595975668546618355627325927979264428719827171/133722099452940753046969471068788801688281442087300864155870869137569493591902378326081568146060215626431195439457982085541631584149240873539723039294860820480},{2742861748718762015514249807011190499711372754606872183527717745319119355386882380758730815119442145527610586241839219725802072899853105708465606756088620363157/3728990207318229676281494781187784670867651608175484265159614399937870113094436560226555289361158860424505353374443211448255074404837396534025112312354091040768},{13151451026812748711496832521567294414807272708837522022526620855404983444501531883022384731882032489516777629974485738902549285981044187196652137485308147/15938111339983957614864696449663706309866837352193743137979430679266167521357438579824112959487049210645645111441916192642749698490033917523239577984172032},{154634097124393524734548412225132616987886196526775034956894949612310613228634248396354279144752544793515352583795499075866532978291914717046913600448324703/134198895622884028728820444544948624052434451927878618239909349754587118277416665660247663297688550723992326537899025288107332055337627528849480994470756352},{153611725736986689215233412551979427285198772779434004423745491608212978207803111808656530643986890612611512681387014637413927440617018439052750395462565147127/95781067442334902779206323969822604956144979764486474201597213413675388140749454339063175784450562961560405274030607805645549842547174103551748912925554769920},{28680015715466950940868561691019835665740725164827846525440501970983061100839231521036086655512746677209279668311656258484139701911938230026745061714038358815/32120407792953430268723888181790017555547260094739779518056904673288788519361575267653268115826189108992577795474617718684532128402389373648332439062784770048},{6786732690645799100289501770803314860797183114109874337056175447628510415015935214616848353540451321003325039564862345088433102406680515404357813899999202115/9596901059441749673734773136126495875303356120850019547339120623873575748366614875826696601744563334255535021633600419114019949578739375840922598492377972736},{7942796038067895756127997149890578553731049559205843731411423267853919024853482973196338210818550115856729346401688136128886861521779141275982115393357362497/5678025596632513654125307436513685354767373581875727085394755760747039098107513696147043368817326015648442837166327422161406166243184103503457064750688501760},{2534894408727831959527789715329235572832626241495959274036510221165368900966681081958455380744256659802225831176701284678059912549162573385923227424919649037/726790000470289209311889053504182829387597251712516049039371725902910934329451802407055963610732985788318208885096542923355818926898945731559405857118617600},{532306530028537162448223593052205977365937396370541542034596536504104091106931601086205365508652312099578049272691527970501475755407216945743095523673170303547/636104892362333265993103278235909066959126289310371902727625150061755858082969486264095683840535694320435982421228545318140935258327890653981331512162250653696},{33731530213488138424467090703263496077724277858203747145608253756594694900767041084887940284030358713438662338808615162280179554805770025568334220087124105967/55667527368810315545880682226099497824263972552816584044085004184733553813486932114579341480709688951233390317452559182310919939650926932945043097301377089536},{349052408712563969631553953515210580367725108967134738483365278413791253969559294596254528766678235274574310759374613999111692126451699189887615066495462717845/432911009750442121562285663259328488501384998156330170989945340565867269613074973481733189382860616635247602171297984777430178600439635984234306301945316376576},{184335754166278216285295944654279351774146673253055961438867023250726773741238013851665215500433114243327864720218114195661158271594266999334854949864688021355/138056372250592776874698450044635619800718748331359434522209645385743285505362745451623838376524098010240052148351429240922559103376849284342508690420220821504},{165689372755583925315363234485701424820340430331349820734979307431049160245153125782847654405800873410085462518994989295308046971869304342727810671750097736971/167998762567107493566116881476974899553323085808828292181433429975220212705355732458697639612450395743966476088345146025801972490117597727682103985847706583040},{29287895730089614810273062133006032701283783548777824481016446688190766434070726988606894787472243136926395537522646223354210279758845956487676784755553929549/10019320971171544957894142630009185871103584182052748348882040434798925740433639297596000279764068785514057151335513490188224501803286857612838625730958983168},{576193309653282618177333787738238236057678737783361234300392868695517665110954440777241211928605689641664362340482807411445317824697865379429074661367397/400804836169170772960760168573666924879992008881279562570505709623408460827266087313031518559020959852412984907044935512138591655540874073850845495033856},{2779681026371962443183574444761845832536234814467039960977936188126023827400914235005328940352041414468538247855842691444032921969916656685546603615760055269569/752985321117237637406371120995838455074216865457691516073802312611378118119733018758021262383850898829349029249948152263766545037035359080872920020900804296704},{22636912574545800188431944830047439805384622379848885716187169103738571443824506049982113805957245534403506525669188623218331354934185039401081030665434672733/14651952904405546902274104309433271990548176172197402353920632085134310893546738594467285428993344178468160983267299473300794285104313171105478425887337611264},{33608280312435138876207851141139932970036377971353985380020127810188045196149466704304826928250119597523471504232189662392392194920211222771570075086588077843/22856012585801354097959837620719484026926579181622007679022395788653261556323860513718162157244466313497196690427026864971629323230421815997861555400265957376},{156364924678695697405248034838217463368917981282127780253318510063559304680342800727122779580045266165594702679029161895453696341147271107794297184904009590235/114677639546692403150044260408709472669641894312420493843100286685152932341308812644726954686433054645938367181117649706214674509308415962719082746951222951936},{128631815759381143655438837852132193395610821896933489408572419029712993200156645112840149277314023288603588808656592645690745426156348814744081166011237763/165633262743277588402940333662425038788929111561880944307327831006412900548692803565215000311545494205318027875473124932110886091197585664879069524622573568},{10871704841129478741470416499410620089136327192193191291028029418071859038802688499411811848943212353140477762913446329492709476673601039089814863656496841847/8488404158449794953816386863301097918215890385408676486064351967314185291761583338600879345401359981237121789853519425880294204262020328741853539788019204096},{2706944082703187960982124556748730767925141653064195369358812382630147216704411354133811735624724363980100481473895152142059678837980028351701051183338758410879/782010864044415449151940361310538037806945974263776559316443922881054532264859728571518924324573035424408496197933615443015796925267518875729986225428955660288},{13954400203836117823123982218379155608213518145846645593452688775087062735207399536577419639350619366645103051831290270394233799678123826156974699845588803881/12994298553931505017828046837688238092533699680789218514081013798866181386788311490023745276966163644217633450737726801503887978215854838590728181129328721920},{7395001145073140924915536815010393786421115206850480410021375833075298400949451038998600922848958082944872975413816605295718109177749664050128573631635500721/8633873033683867425742948495238953304911306455150716634775461883443314455082730166601852487052503156818462585961586497753703445780401286119587526697570795520},{650938574877460896330178233759084615683217863896588379867256526054408264584360699615867159392813638872734151671687180452763500892731566431022303545823700828331/175120159650072470200779298359558874429960294411731860177976575901199479971890229437008896438931275710808926414995271899626732885568788980820688906114312437760},{47124770488795428875054139051591076698328449286479255135502309323588249007734923905127003343247675763720385753287530215615860668388633218097819042783587371321/29379933029559217568987958100437970911554586840738158373818027228761512830660641511531425160011415170515214479841566538703045112002477292814786052557969031168},{15468884294023961918573893328839080200133914078347780328483782044657430740206067713240996986510168544481247048597506863833865993499762077802613967871453347129979/1457363362799337391977404125168464176991013499708113404856710434382308123775667715017960784007934342620936219190461976519435101757634253436064059894895489318912},{109910204629138625844226246344488447959173097596899371708025655095562299977238826057631129506275121290465902610931949397923299382598068385959476595794597894725/15147540976919987781494650513713461938510608641246621216534521430819481408394160200409847775009178357086087310784016574715914401568413512351263925802631168},{15717388778388631992596316822793942205007902759728547253072808100069496208027512488005069397752496937349701706448617942244373229800986259372842537461469210647/14571703735703619689336465762484994788469858103572809210127886578179295080192276874455309117510358073633872971713752104463245397819390685155458147941947015168},{24275018187613811842778221268391147547138945616302522115009012196630865349013131327913072863249448519444908197547974251303557713767927654793990515275004984065/9033585529619548924594768318582120138158405318584845656011354504648358270291843823469221665338676849319474393291353908414224882673940806514265938831935012864},{320489219870334558135629610031595770044075998266518445290305518220643348773862294777709222502087424218871339791733056186687006731147663346283306158224196698681/548913699204910335346231988444251260492297559485254490113752321022487045011262478720688446265011901778338351487389679823255769984441967413956214495459667869696},{1133904714538749381801303565071228118377249578526085624120217628435877678840819594195089874487852032461072712867117313910267323401887877738673603496454932029405/636353912399944488593475656084128423477575016011602731602305168008748443020666956651835996154430172782449593492717194706060863443807225942901680279861747253248},{117478963765743709973665961677389471610819506350027756174246982091764412466954385756444966357050926950409748808095279450075531043158172808082575618134585327291/12334330976559523069286260413003226303032026479733696156082355152239571859175699290524802178873131453461469082365182982419561213108649206582160259557906972672},{3208500467545129499427891072443690123281626776585369387340256934429943768321363813780802337858877219291881677036335362612894941282324226123249886421964712377/3104842933455216934631685888824151799655207541869371454462505933794139182717411283167646857920474085789715735168837435259050246393172447781611047560333492224},{10119509694784117414130577796274552258034643429128116963278464944407821737763494364214568064135087376255130092728482742360854048583535773848987527677953172745/5980487800493984004595590327910632499171185926843919466335667855603579222045418477260963661653405725121664016672432868697911129048925324546866639655212679168},{208231531265969816341731520972182692730889708545259709899507162460999686437152042687172833130827493923025508487888011916628296122200703600443968871563781258727/2326253508516315644175125973981512382297473439218413997760627452279112521418782680453256276018516790194043136138800752805826655194447814660518176366268514304},{2223539537432035382463768439284146490599220350278551949393747670300176102893532802395026529038965615112925603139428458905726700478671283758450693419318009573/2549540987237155071583222784678987956664744679144438569100537028334606931605780883577110746712983732288421029218573818823877929176255947081885658653606805504},{1386112785632208450748908938613054534753993993084194207979357297231534565663969032178379584850743052512348123357963753947624868281516538377646952804272317729611/417088488590291079550804331888876510182028449412659936444229577931256104456932497731225817687106923151789683163236793785039164205597775949430028263366760333312},{9600912849746422835572759972545327250685630418553505369140512334233936947297495467676050383640252981787575375303269700276984932450848188565102026908407772389/16525399875398094847493938912050088316714147518551477649938478526726842297990002340858992381732306534261124891849421416442483685435201396448301323540206977024},{85934232955254791947810468387922604642068894181526335529365427788072473518921584656543075532679453640956928012358648149632329888926281063790616615124473131/148270588549393796863996757826971542057745678619202669128310219560709590411488854813787397609150266473873988875928711583144482628532597391350504872059338752},{16434290615664040739304701402853850364766108425861241022587251528131486216051140905236891660783555351368292712334480255931161791817059718354024571638471660952497/28110962820538979204141710168654881708880857192699345641203052516906385964818678970785971683541831856556176123785991525089291633701206077046855225272437305245696},{8050511735998875420562828005545602603814626406519272476349627185696469498033296874596335849228999130088936139987648123073562267158195715859876175517205456498505/13266888966515052897892592830436039271538722076952733957026926910389549618304460101133952044820542471345884482285207407013951736173443644407083713369741409124352},{363723542737980162994510123546088567281779641088443099264211125322434197552080679711295970045329951558629082434427275154692037941240066218610405428394210573769/620244394535891732547495406054957829471500888652622775569347257951962105226844811011550157146845324466958030307811038652258915703622127287903208975884359827456},{142379466949911191541436048128087887408732615538758928085334280925387313091889931201992384803726250676462628820545986573513297661060441595357256244548445717209/139656467363693462587175331825457533007758125214543861465427195070070842712595964548432078383523127067450097488787841728038184429113525274000368599062145925120},{8361482947351324550792537033075682627190267615008310108998608901214535939555461132103967200515793365491005820172284231215308785537296635844196032307657553833985/456916282692699294713077424262775852660811691363162500884271534976110745623035470086695356872703596596973712770631349508347228845336137705336588796591161212928},{953864792776828641163802991860547382054748919180626463636917487041079054336970672817155582317957173727765889725222511847193628837360448284039597562750279859793/45310932020924951388854657427837306757857478392622463441753787182422205662404250725360483947092047123015526280050544379591510048512933569024283528008524890112},{75222159684082645972832726329694308576336029614758488908903395793708630232651800836239102082662817066547827482780194749469662236251947362566489525898600367009/10486087290797998478031174420218201111291168063050944424765657970410099501606875596380136579188648026289494459987814496286746088640103134845790865537640693760},{5880397956275141162790224049352889930259763729436651310677823597879748021908477704904934294981715037425126619423822259590687769585628197833015845206184575597469/2423305197616665915478857965004799106231423575238601975509704136170322127557261828896698484194407221708728195041681300759165536722667328382667142873948629237760},{21460685572248979149093224740878612920448702835839517218493045053573921869758839224976389464136943611409549730262718201685295718974886388900721767013438835435035/9526757046458548329334156846060790313980987865608511910899115606856938682680641560691603218846721287096886143879456557007145085028099687790868034889143153590272},{2983624545702404423604595629657652456315672222229790748256668339102172970854275290322289181952593617268639138372108584874925598643814706061513621647278283137/2342138723574986586678598062051002930389971221705148604014552542251589000754609082354148690499851978184422605803306386000433516033041998662795557202653347840},{429811289168301539274414420798542061781166030937316612490721825234187368223236921905154096477806455476842554115765307115717820840096205921801397201615704898027/684152643553089800554996744333656566725082550098143912146426487246243685987431074027932463206276392236888656438694214125296129458691245528898815272211944833024},{1249849324122538907578071171831859811033168777252435646784491882140099006358248476130930350384431410791720809142953847411926721006901211205539875340043794971/522462611532366393567980929412719396454236562782682866712214305706123839071626186853414118035358197521851355629848926237648651447068503517399207511769219072},{532547903614509181075678932411220456150451585850699947961624351176120895590208966464416959584031034613150986047594024947264686420212762997627319415054538502711627/304866566310982531886402571328585735392747641570857940987504061409688294105876996130477796862612587016106453057285082152068462301468596852658285544015132589293568},{1986188304292740010180316881460058401589712260549510559045328476885238351682152412254839489297542068613034008177018399585827621249185872254869642006138130566751/2040078484277014448922101999782907278775382873685289393509975075170898174186548065700300825627051065302282017433466925377421941711613337947289036768112376020992},{1242256874464807985229873371063009026832292269014505984665699575605234737263376031269856472756159607719986021572051650180689317475810913022458086508458105775/1657798512557329445622945511453718090714317847786302839071636618910414506526561338069697216652395202163614949916994018095375388180527421663744272279166517248},{7103161012933936746433116055149431825209289300089513836230787136674848017793179224900319381305208558652127089860520160527835361585357009194087139047364799969/7335620970996002880331852845482123560613173433943242847558366344057591709417560453551252369272855287908595811685893957073340923478801868733631672214087532544},{832532653614507266250719643724044071979859625811315815707651667421959909208828834747084711021598773256722381769213382635423551251719278531329075309549713093531/1065710417337343362362410610135920965710971476879978465942874768708927779440752514529899370872469916311344768591937529043099034453979310099473298297879122673664},{6062720658616094420430952868273685329629572889017252630404621756539987558786855839774522186155601278911750129530104951822121942406054498265559671454289047323969/9487508809671391648335513807116155664906558394501479668997523932408940922637729315949462841440691961719723909907741097543808841718764146179415153220541443211264},{265403136580366791834129607931893740275344223054658585763388256357198113006124202822251624020379076799581787383127992609195293833700680851432907737342764867/367019348007357660598929965343773161450611746870628352665460964201116486856507292628525410477196938699345615904609206727430187790591260846099218225698963456},{1822036502901575895178234770444384238211936030311046639789702404710476832856188058099267970537012950723661293452807773378126022627601353773523619601593020637/2595520414366464953829845675736623672849415484839431148806259819297418951512935258058031721147239255550173078963679680142471336077741025830384936090628784128},{35034170829933842389797523285550304695794988736711974031771480696485438821987923876884069996379227516308387411214518495764604651557189329153474570671777865251543/10182553487458898421352993065051742748610626483964263111564708461148092481064908815936997767459236926260335685314869922891466505347671373840973830274403906093056},{615119928908967320246055095298841302033203311379831463034946875866761163575853904191559298338897618124666850527777131452305020838887983663884517889305513288119/326493730307824622162770861005213680832343279008552745485684107652884696227109328488110661588788833653164031790676830951118583129475811736298033493870118961152},{98840749498045039318587431487318184096988333744443395673440266329240079840846749378249333905933769567678063530085997626695921017932236178754482971315116723945/114773790485636943989217672412166543794824620838613585152522419456921178801780456076697582673802280809796474329354159491975645192037653458938190006025319350272},{447545248721786240778809850939570803983407961463976208458200675503482313085823685702926351738681517855349006956951697242012670988966058716604455803588282891/744636229580929853708703499244606983742197423706744608257710900190630495737710718641824145682402255087440198195233270456816261891238699853048752354080325632},{170286675747912448798393805675769682953342429339383002227500349270333885772727257389808620335980280067041302583043639401660698927262763096335233995258876650401/159366254424366068741072228001875207208057426214104730808485799316399327169607282168887374962883611374578516994626679353455686975958382742799760368765620977664},{889628219520987174992621454989869862602602184889259959026806999754385823994579507332053938348043700441729031355716419312851136557233434092547056921153466611/1238808481237522586044321046724585102013436982843888365489730376437371505451820858589243262418095801373517803104302930114895523612792550888094286745312428032},{11080512794936320519485892288151672280535251922830046324787309546667364325787521074322062690268963442589665025757138662181387883181997179575816025646232528761/13141441081252469111075115720580819305454455089686696976282152755598470050086208708167525960263863787475821619156296341825055737868727581993257189111465246720},{1314659378697677002420600520912319394647119165587814111994403082927268460158093684449341357353345141923866626907962861272788567451754063574218640668551517138705/613235838743624544979903843779595744782808513534330069042474058580130350019422000451900660475563708294255315910636256467932903278237619821585524135564939362304},{39176555623670621421359755072486133635604589485720533601124885363033897466763302859032481687892735476378500838591489038208933004146842819444760342435518842503/33939379284418974724048331634879032329850872315512479849508611556018096565684151341037331228754379793367135792101173210290416046515258430248762960049569005568},{1584615998015255825804117175911192761597268189485172041894831959674507784908653028688014530600161067267507292143463160506638738645683788179892775789252029173/1537276658186783972198392536072563103491455217640959329397567396193033291798484716293560948050570904929666635551071600485491710934189205810445985329342054400},{126731027813869702157240805536331852008436589508891385748775355536867004992103127076308887297680745078016604245406817037790766111979850707668596468646633310817/124119585526229330078905055591634744520710573406967944776144977851336077339361786716999352080973618599885881762369915659760129674770537469807268928155046379520},{86127087905199059461216946426756279556292241313203105314460540142469562869175665254009693515774456049390290653834461041811901783211831790142598839485479170827/73962793594501713891049504165806293717776972080120226470472091950963031625475739584729634772984322396765781416154661278940967805954213739036165034195810254848},{3263622370974976674000997140884440857682619003084377576839673902797040025715824917851439172443729647867016713610904440946304671682364913254095345878221701555/3088711260295091499832694055852997522523209365803528940802510395920830219422604625991839260632645796899268144854779565437832510911871640637635441940672544768},{606599289804639768960895048204616198176345552035991221963472178981432572851087072681729129336910165373052273972608594514824325266489683537504015661545051929051/171684505516312083027762393887732772796963108342403093710023457330453144558444090793421418760621739268254159661968765130718015267304593365108357567274949804032},{4663762417857427779618958896524024434667083860588689406260718459779286860180294501459539061759217292159995717117505691088002496742124896217039708596860398415/3670203535975070140476792169569585191211461791926247924858557496512176930307853454677771997283575502788409247660283627638133563335909978631501867739821113344},{38929026772539546312548807771910001994643820571592279923720701440879662491788843630469403056482506358119899599842878140244770462685154559495285666985308464470875/22442224217472795322750498547259647781019891591156613005996683691700772325269034088807560142925143176348776050740363813089749458678584023135104127733321164652544},{137846160374652390410139631318260370981929725776823993704679982180882444407562028614488588434304360696126685408142466790569310068546468914441111989219643745675/91978985515815001565902310601399769592653059508593148376098449476027675527519277908974058416598414991814472447385451323773334624282657738313668234938310197248},{74582964607057080865434457393489199174400330725037204020491613475531298924646691914594992573858637790034538229483717913481737087576475573366337146553116489/117928619113618526584052443902748225528302433740267016850204392371117361541875566140754838900900160059742134096303146809439463818132906867815027039203753984},{27144075566238445514017975666901385364820004141575885443597976609940631920592040274315061726802110083091577549754984433647135171125268954794822979194624470311/7018923077370358142645089850578537204319816238751716127763671227303339865684513451549140880571765298294229312089461572691549355637261247780870226520993431552},{671029591826290376668446705568197683973873657716506939467414104433897788929587001915483114115736842314118426637435101113197871931490982946225738423607826519/116866232033144541241319855013655279422001278021626310088318395176006522747347708444569779344028130874833658683649920931686177159387884581084212693447999488},{32286822829849040074085989508497221324850420300872622070912393824406654343920399690451374183610485391790172380711046182869243102468378429139228066736103607192155/3196179795978470065052154800805746563762564413534300275827905775585449839967806222451486552183257501398138968913884166249039911158508537462422527304154530971648},{8416834272952065909115389296097213700971614904866530189757578431282395987720189392089156459421178961297328702961526694365816003955209630067531221533561623209/8540560328603969112174410301200570661580039621815235619643023553356883038734624052110202965297650672167696290793447076579385831343662443636812692521886941184},{60693399621471841797211719417577005974134070073206616405588304332606993508614373818253740300955007558850783714801723966237935576022628440970569518915749650487/34200074274711029177535874865243335632637177737994179352396865634892581481382174959864100556453435646143663461534204193214893857593293775826554840299318804480},{6984189607551319965355091874733489360612411455532675408052409221060165753012381197423436392470091641782657541374926316035159835349239555922872588714834313477/5613047045444886212996353250007054091790397176372376593317723507685403501892369961521287253098867028964016972453733447661570336640955002175429253055610093568},{369860369901722258058273601639997424072332577543879435209971417194871216825549986710481051744763270802135935690710144791003435009758899732061027726937631409/20895715500872468809168737001501067715929425846032365499738707049614708764167382484904540824880205685097552457923807384613495692216597808753531605775024128},{4225783040751499003813625881737108296726467329663426485123475904454341002827858515682703142471834196400135600079537794672844428197225590185102613586632399531/4356249858363926264725359709385115620633419053539575620629168817835580586813890982243727464002366887467530449819804387787283951954475686856974792072544387072},{20562648116438575221283232971574807749504442027252374751363134780838277325016912741968510922142010883594017321987031947143393759200859065541683724307761663265/20816398783204688422519500670627646767530800047483589195684603542552832109739094927756500755768955830517961557385185961462075971727421127946915420836743610368},{13409999707777720702080671162726332358774116355421992395241030157373534028491825099657776307017303060667753965409863499093364277400705176099217470412930945000481/17949731789366268639250938315942584924117431429044685029681402473860309545294067841370640084584083058651555874811922671122358791161533485490038299372825845170176},{22731206600614581419795695618252985011429946645954790980528566040294074332159568782887662049021368880661705747426324607824682291467019532229060885525649386022663/30457925515799937840441875642417653116758958770689125168887846106200340293079511404615356534921847204572046086120680300455295723830914134929662771838334279876608},{1043366173093636441393582392407896990071121774358176824543670429440004479884391925871484177085224402557307205817419414509566310716233035455454819935591795075/1041042013574846256314849890304420448646036908367198920435521769468395413517584457868062995462675220663360874218046900709446363936270627762688980832470695936},{23304745757710124608543633198519000503861150760032063650968749467316986328581918383050403036543865743590348365611829352127194745447267914922436282306261756377/652516230300791420909676830901516317265949369019299829426078957514566554343655512952484042003337564759084428232535861025455909887155856053369910528170786816},{358796503223820951692400538307878452981474886563498789237037412269860144546463198058754041489814540469265032839483044708447070413093222469917125054094604491613/234841933395892523961976331433122835232079402591211113794622048885815043133437123561673176592476621709772348526501342682688181941873069513335975483993296470016},{70065696695145194455457308077864742967580214076221035955441002242735085676952310198659045294409484746494490671380636094739649913636417870522755338727581981353/77031256400178706354541418907193968787946130701710648825867648680540209585808412010542701159718719695651139421152849588454670381744244127073366036143413395456},{3354041508165565397164443310488461363056556774326195829721863375231732816868565464114772333763806519660427200880471122453650134650434578839416712160383919797877/2027548446730216043448002808274786805937901502195070413036103153165978932081032365554599407369408586284268202775752425847972748709639239268338733207732453638144},{922956087795231004575442978790099798314368940637680354486049684023368590008169497220829580693196961969389631058224100302694282595058632420602763093929777230573/858246307138874932577349449013992854775934476644892444060743754121262769892858433031046686117188148079854059715351365903295578982642016338728705816296748482560},{120348628193048453044125254090673171107684711829487362911115252530394960235056409679931921276905416030514104058342767663356153389802504998380107116636869388719/195818616748854139881882411778503313767248353291037641997744280915014098208445281457286948971646260410687000015583326286239102869617384413927111709182226595840},{29341154914289947651089930048223031240454024518799584224372318106737235571360793421818211454498271039207494566318322459183656386526216349905448591030188296005/34606851036883259801162070379782723069158419388141327317241635174087435731466421296779087408793150681837840344883128387354302385518923847372775022167121199104},{62705243542884824198216377248957881642266301881945080829489866353028086320053953182110276101696732806865185981066537702904804472498352893282511009491977985377/42523500734419994283281422801969622080683995023613627210476704505175050303785310369191870453044295422361971060476917757742903163185712742195705349267888537600},{21856990035294016646682692048016159687011210986726483532897512675929899305764769648739555159996879428809226663041125437231106125024387715894836127554215163109/23460146236661900863527636197837456012592812746692285378977076162153681416021170155411073835149727092239231844193873492924094210623585716976473712306224627712},{13556980566456736884318302954696173766860738158954290130891726531543354118011648188518766918020511950797943788517359344472224061337392399427039850767532243509/4570184163726708946037699500504920630769604622482056012882346779151797044439703685559065976218298104666042902920203303901393539596403050646409034555331706880},{4138591075030657823059608447968045521910017886014234864586828510276568824289936765428438815343928805656218649704018526835746546748773927459868910452702293567/4134671298489967017564134383895762465218873093438698579235949401301083014626164986303816504142427810458753654857098798504592762892562171265356543126776315904},{35803812222184713391469033353690367589227365891514205950768738993054353188139470113926349187810326177522539485492919313119358654258574612414470854934723158981/16790733506384760939680249715147592112241878154643947176281553922144463543002824515486846923106117703836251363157383967391814632717588437645653144166341804032},{12211987700629854111944566003650170294254163279209538334298676343799992110233660056417580925035580898132681295366102703946895905358783324702859464510470883821/6362235184006664929631357245348898663630716603465642050024035530514330373863299455240056252293346511065194776148163790114492226606255488516358981394605539328},{12858802153680134913069680597288157261557834684260393635591125267030082994778826716703593444414017391060718593081247811248920732389710406516952632625531554111/2100823414975170471176297178252621996762308525725951268538191695680482664810843198994681057955833207081488004026247243670228780587997618983527538385410326528},{69661448880659115744278799796787734904761204472322738639260606000847486513501903709719067419917669878671032700227424737210431111707161945416296736520092032523/65427673618586746501940291558600058329941031279542740105536917960144210540249007381938517703573277433831369678203010445320410456874405354485683419457842577408},{68404572918226022361393610729633357944845662198674443064622814158275069067980374193669482819088182332781670475539474754005829496072597337562674786804483032209/88156490784941023051436310737346524338003992928766901016986170903780132324257207698779087172955291510650871555996265119421211514449148183154289764077523697664},{18171190978951914206607851430985823926413981623268229087842235013666925270027726122127421925649032937748600039853203277098246911198731324076568915811592315965/22399366653094328957043194631452190698113527059235619818972716752115371610888690931202252354107231448119454147392154213733736037615316560547379926051076440064},{288282763734325354175208874419681124416991671484319629038122969926491830513721049674783801901539970731403273552335812300579909212838894828791057993230585834337/564182215600904839637240369920291686095477912318898170990906284247110424656749208816737035836138402924467242960228420691561980649251410776281756180338194776064},{2588080205151259427863963893519262432446521992220609038955972856426419143612243674994840061203519143411187118778447701659509169785847268027409743668034412661/2014023633487441645233297735912358190144402960784432875628442398601839202654475012627059952269398698357241805714517170292597721940163354137875845656457773056},{2115930048612102343671018591701532521898736589543812654719874427007314973589205593164450098873954277062670170032414015646473540329133536405035079789982215355/2190737451604999244278826625784390392128788233615334051441075973527162761323522337281429966553308765337272570094934717603861766963163148682682669004092342272},{6384951018785456980041072549024526152114183369469807620331216564100469473237956360617943362006431172120352804455135634882004957348536519968416725426290327305/6998587102533628422388731623626619436408456431769689485695433548178050712943609510944039374940265955404802044551617095868094814206216589804577038816444940288},{577675109510900563522245303153016641565741173078476144544809610183782586728036974374019500588342231878135492075615964736699098330772790037482600183068806805/822780908728751421816245668575165130553195997771050439822482437734551449130536901915318177704968486682462857018276917160479352158546545028566342377035268096},{627083940067909631033840222498599842451225698730406844151966118663062394240284292656177807180744635595900000575559266582137203350047120968889607499838645447/866333005048416744696798127638324925336613095071613233518099131861765259906903685810438629767200829508010321097684533139861778525274504145605895474412257280},{336275127220870928597930581792153466125402630083322552430723276382200067838994037146641662902707345074017805786969627145160011400688412852194052563741364167067/349837984057848877594746708582153620811942229312656529449042144769698943552890687907873905273810023886860232672270135928134550394197242368085355381975811096576},{630876464468640171268858650824429871616794170797402773562390723860821525083156586010474619726838392765809101260241105312628327035894003335173317795421978225539/502391508922676090239800237887408033884306359098521463429633518735208417335489972887407599576142249363971794234542339185443113773531476421353388788158730076160},{15049346285534692739437093269429380063587675888333966792170629203080537704146721767298560447279204785447433199113036896614361235354864523873012559580899553133/20228194104838374690117532709211073534012719485258516803925193349999189742944497205158779483516239857819435437135289951989716987545614198620434907244254986240},{8876094290703944219076410457056328738795681381470607826368812283016906720869626919505638186894615754721611214868406037434246276774771791451386211043224471423/9517465639209332073388403394509076276598921631495938413043954834890178353376043376666527436703847928100053850525741985765735383765031755125881500961163182080},{746799287094284619143023424002166987963957345594378000810487061686131896768695347831270309119061805938229007908902692607463105160472599366191808291799757376277/1245594204242393970005949318644401921530898840749340646797296703497122308826083967894438612392129971642294747634075834905777148670713736098268740378630161432576},{19702067513787495059615557708982159306685978409402893589762408717931801479725776108454717993237515467161344715820894201365283825910885242776752537988367854813/19254271266580353963863091403988741877117840940871803725140352901827594745363047414693141308784089058267215339117625242658920056474804306235507763472091840512},{2166048987167098800302359639258532269289244618900922669801205528656542502195217199906809093810453152706953292263381937228640576823070966709495192165123867127/3488450440747251261011285972407925723008204979872597585517443980156122970166450456789634281041360116894506156646521391581339943008759149900870924112138076160},{44410480999291724964293231110632702197635128932041012395149729511996064958429133755879415512464860548850316108970137837470624930199491853613514673216601820097/53597809546146134641014561364021498301396550786470813908282358555871903071412673195559427198957022644889016406229358555077859858233597485033739268087043063808},{1932607428235549532027109088525082860352183224474217784295644676376191978095998024864522035851403850621890803527649038055747654302672554792013509362306166749863/597071350171122523500145935669121933220637869576221047265963365227634260882077529079672934525900913618593295114772877247656031669898339822046495281806619705344},{165996692384377246674078025782840688886796499849717673204482528679338258715761931527668972592663014087203707597176027480354951100148085330003821096007825824879/271066543292899853650746461828703211930947243335174282328560792453454620435398073691467109064469296575789566940417284083961608657247577086989491077762651258880},{11699668101320795466794991294700615643729228824240981236101254427223005663435832309402065637575106350749509696864407815993435001223235832635961436742395768339/5330986876808698698467999499525448713899611524638690441309205226018477415557688805890725493972539731800142486008629020658869164336677455275276797474867511296},{48147071644766272583370616608550242884246564888010428289270983446400763774960754960103942212418501918485088394808293308118584734738098285600789183422947246583/56278592671116323740704409105600455575912309150721326787100741082476882070767617195222310811737140540491813717896187620496153537988438971559009768297486876672},{112345983754448820307544623334331836686087475025201191876520435528464473584693382410448998247366243542458082195944364360440112105951693491796496341617840448043/131207290622233026730755589807036962458504471272123755989935147773412563297160629318082931205476798184351848769101982790207637763333766696344522461956144103424},{12043435827910447513610999707160434216250441947304705124037364070553489567522854501970270834183698992839016170206066699133205037596024815402655006425971425959943/14295764069780882321857780329586944870183798751794954769727814853080638819523589596400507584716516588657735382095190772299435612491692100907046828444994936766464},{139031130323653736215889996297957488704353817239617670481559324497530233529423246489462211021101580048632109547749406379160427605349210633355763355708016000370619/5867970495184137497548621074811262306277358854022578349821579131020081068262946984920235703318631520058742688250862103528434140338464372546843857391264674611200},{113505564883169379483646815962101461920616138337355563034420664241236439366342237179863520865232527164131629013764518769881386261530769369281486445798039645101/131801620721407407417398794520408573666207204354204176734774495160804662089378802965094368190760853913865559687201672160651460673945305784155602341969199628288},{54750449876406103394195783549324173222600954782100740618032030917588419051009611946957948331663348704234708070876063469201098150232885276534423039241058777957/63794351613096562774535060121645211446890779389717227917476235192477711884400518589399444227621900475099912069327340875936143764171125059703358378846229364736},{1509055058404691755394736254059448439234610725756233053758567214770331010262201429556792957545041842371826605792572251672641246165820349814605112715933331440165/561098040637905747328489948219139657731909541132214773743892039334370607352640100697755871131564302389003261293973099074043442403547138606613779671906911256576},{23908697191913543213623576388757523618370683394406915869356990143474159128829555913575892517610072344718868307125075408924599915455648475834187348919395382233513/13107330448762634192694291455443299523434896501248842001968355537845623133269448963692554533249811378421686701139560644385734373112430911213732105806345456844800},{628189542447582054776951356650869944032916663844327047312976977008632426349624987617690215601835256383467021717138067460212123571352458736623967052149952216113/823547479191589380711167948645376077522297140184046289615603079133736375589347002802713222740537397299827282888691819958303427118347925127405658564421016879104},{275141478318650857312883220570985891789950326631505310234598893025991599918694501592136564888850638118974361920287685535799214761770396645960362367048620821009/370830102299304267055201667982268614794853742734134375055799065280827702687198206781360349319366224061108507264732864594573442859324917107289280776906653827072},{1741867392957377442541403103078651351713164122385223713863521103523258117099434028172329575864047133423633538461972570959762376672331489772296605916693211686671/1545467952124544596737791605174600788537130788207692321133196212667027586402692576736382769718067000240153787887554907621096817729647721760104126961621462417408},{82772770800538632894093952332581942251850177527844601405392818548978758244868802138765195656303099350551857825889118692856162568663226998933646916055697074982717/5403778566730662816897947811351606207882097642650633601119426985440870852507318486459291398009128759866720713629544480274334381918119423749663371247991264378880},{1161377866018213587692913013546503330850784867046080175172244544079500686698907987467628054606018652146128273983561505844799783399271669837456901253255613781/1449589157264142770487541002649384350045908081065408583306316558870339646790031749437355912533001866168620436586121195369058116356679582016749015889940054016},{5594430761059547545740935933729695782752135483811849176068775898179995467689741746651378048571193211765883280864799299277372848267771739057676202697941980595/7333581453051799787023597250562360195475133367149049573948614376899352276670389918705569803999118996270612404266467129773547219118245675683589979134543003648},{10564075491933394127938844487746270810779734144645474337972107403618984116400760260050312800331528247125276324437538793721025101995614138212504337272809676941/16189646751543893492809103784163398722898955433934761152456448736379993925273206431462996517284100792716682036904149419399202441922722222231078462133372977152},{3845591525361784100623561702196756668268602369237302833330548458471271724893623030373511525365236522390480631921005105513637050824750419556514901281589656727/3493882357782127984962407639380473389626943563928169437198387049474635164509900809282763923889111718630177360825972719468713688556728357662352669788965175296},{1884283757492391448649370339671990383594631612331035452370315639672038758070788802224798533853698078858218644547925456765968268764214218276899788145954735917/3078587702982195627059283729511754565369475513778126006216081048122534883356524757840789558364650664253383837870157968515019023906491387475382142130708283392},{307461399486428666405421634192403800553078542332192672888489308813292063630790662170151791788177879990626990367421214621908198603793363714848616446189216380769/228483509783187119518897084617760514206445038764301385943396668469508592691727840415902025674423369331331285764541814124244672051433209556347026391083662704640},{5297574540310757898343148068305851405123545821930677214212945998211778281391045925657205195260293234330029268224490631344417717895566816583501852831353564315389/6371548962105659708288762738154077009204235426280056976129544883374423910176633917957570198876997939231766428656144337695148414177008560974121356783445134540800},{36614135505078652673655472895868826052283893181278799825014996354946845720729857097711127156390225558397019116181276297950579524571045940061794816733957937097/458418613207850881711175058696203444055991911892886112640340282464876878921527707236979856930847373730624123649883382945878154012390788322311438668037881856},{115335296985505935537305159867542305113663934800606588838067420031419705158090635654420032128337006081057233422112331302462581204820420298585909155792180450763/204187589802894122279256167874194788447452919871574927469931603688622293668882478043676771535681025462951044465602820869880792602193340694778295191163216330752},{1260983607972653693721876104988969846383641166860687980245350180681251319493717163605312353849766382052640085194507974552926262984707423260758445068355744819651/1470058382470140546614950030063043930695564854274532788535954572816012717437641456608173600036684614454211144529213180654242102396015852243934819051313655971840},{2823771240672073606875025219372837833370396645905593186900402738560625646852032831781421727310973515151071928401154472138596103568630729528576664762204571419/3120841602603625941355358202834770173130965636017317231908942944358578727053350579126504074356196841554003308880175882268339577575674010930491097108971520000},{114623654984339004772191450062221047683484882160868725546672890917079498526346849101351741836190914766156068993258934852563521699298805729583103522994222452763/205175222819032991775638885472456006917458554818265626729768895115126609879460266749442827240767908621772306022960736090755224662124592742052977926844050309120},{38788573888575157476099113233494044987189646279689492167193894565564371567298024235814348497340548510039976109119967884742361790314326130292670259668778719545/73673827234221749507728858365667412956091774783311028248841990198011404608362339652560067170533823978251920626881091177933709784839716085974010359210171170816},{2031236876842063463136819633615719333484053716821208159595970692055408556565511732614057998768322309541308478115236013820084346444184644181725598188987428337/1413544282548002766333077773178958445163099179430982732357351332558159762953812252431688228869458139941915030361480092456021785897787126938725124388018978816},{304210406956067224931863084339654080291494706074384491421806507388849037760901871498469753923094367560827549380007779398119596312233281912831966246673855013135/38292275544945370494391446113904110972397363542618755994333941072235251145844172516354126453463736452104545759171814170884116920308411501683870990286909865984},{22362592775366680688292225742672858281890536781101019002105170656343053032731296314631568327924259025164094438424072509060382192400011305075988036596453806287/2658614369420827878133684981950478276346778046429285500448805844521671709459692545132299209671868782184736053806047456764536682070506153864308379282523029504},{6069747509677136072785943356819242801573816630038314478073931366539685723806224482499723754803299849299140645305674958002742888174367696271764271086651529273/4399477358443279121365452724734177902537793909122444768252155333020948870860024055592806170272821907369876739784873037884170561721363008310530713429534048256},{956767782243086386452939533134384657472495501885224065589934028314280000118306456660195217552058451511775927519299543083878627415637818577534957962967289965/464247565920259648761325451695167146394656597392966774684145167632723840077947684399315740228598346847290210197158852003432918310245807049366328268958990336},{1641156775374175180414276287271755746013091384601064289123571187507697734282833529338013814577678837440144843518963923297478206744326651624092766982552701625/1136910994809826190293504784735964061208520867928832494328112038947790134356049670910396697309335178314803204337593676409513174122976221549037809022147756032},{69413012628403929491261788759371492890827953108407226314125074991528489390042016051238984884300130289661069947005123278699850868901355386386385577633121214931875/59698931995865163642156736410989184907740491132746090639169122699562634203442178493490308263286072165645620285445760908219160306384932607790838487549069563527168},{224510824725333060136310143834643417878169000579429594023436682519002892275191389501457102111383380713684876168274983524886796835628653287052869781018099053923/242211117870133300348706237106470268912678511092041747220705832415606203750572966523339147351296224780715398902542355108333890786307966267614435497009756504064},{9231935463660649280847144972639731107075940296296667209217127022850407176099403475147947201182296754301107925151479234839529652735154420750665296526666467587/9802900905881567564178716910001101293675892260223108703217169594968852931271865911881858579015482181689594641428247585814663093992423075950138326608640999424},{644455763732739011574633189251905506885461413090062250584859309173552452551873447015591795166344675465111155108759281140262101960174757650588103231199082674143/976378894457619808093300193818771737508567016508770588836263830091087018317742266948258727174083048390543286847942490716637959569138749075489951012522223992832},{631619787814317592754105052128991626721852698632978655822435108720994782662104372825974775205200425258027788822633138896655587925706867923587251958988612121713/164057476452260159533981528867825571605018973702911262462943421497067018356647187596882161907209827132068918006377969070308047326407609293167488695336100167680},{4327295116528114154074717715752210541978959879904757445405670118188427389554598654699238030123123535962560649315653915735845447236399154409537187073199228156951/2617531710812375987535738968769180455269556307586597674830197738760405025317775181290074544628590234187083321191474465720292373299639473018556982240356723589120},{202981713136784308565429099646194962791779765468312451012904584096240106943442187570768975207587494438358624171629159368677527469961805885515533161426754312245/165431379989579677191361919315714208418153198657028704540452128858482606158010506563606323767106749249377737919081798316081914690074694168920685679169259962368},{15882807410275909756981629848773471869390302000700949841978004854914479528328522985229044529830897263295238351367277286874821268915794214323966805452244464181/17577327955938238121372118522335521529838716571347399008237543154560477453579262903012526371847894954211213070771548792346598718604026345621694054314167762944},{799959514991040145113208287291363587588746138092361929286101733056617173001692117027845187940902251250541284081378748834920990553171455358486551658319535314793/834987058332392274035845590880620056636889470097018477238729876084146528448285132077355761416609279933281918753466795829312579819224461322318833426194799525888},{1665308465583694620038887355346450004703915637378710721206752785885720051902299238295578842959106088781816824626041662278596870660513813654551720000268815682055/1340553183858046057876563773165541955467319740226091243869922181057360543440625110309924790856737043619414332623162385899758795968084490490021325065136073867264},{1383997055441396557234257367583840443329358462314825895586882651013951599379792336720014208146473015561754037570470072243973976098286281635897257945389295364497/1101824918434121264814976240135406801436796072272814177891791660805208349931133526948701262986225283680229456946758214355853844861969961696209224363225882034176},{277207003029142613765170579505272766048759504859313198603930968247696462974104682461950266205645687072063316437356594330995352216125932411835168307308647501895/305094096241639996427653623980103659032051743835398753080570979573162892060298574427014214145362730826805024884764060637302830304334925434847567631778419048448},{68755512840633640882751047736737231199351646198654841111670476105717531481779420773512876476505791734654832723695497005830935706364850243222184746343913847495/78581210458523124635018397088507745767430153493286870712613461931852150830422483403105547907519128664203661842228558560394117242678424235992302716099988815872},{566778227898942525492962271444125768565011530931868115926173106565962352479315699992566443507368839139225984489030052957045917549957098540819010370282768626859/940764792737219965309750333948580827277361607833672700699653166343069044798286103014979418162299820404491181602962781473346905659431025409575494319667180208128},{11039245323435009433257696917214611896877371066243260154887742314351847674940440866669373174572862904589985811324003222859979323485463604311948505853807923143133/12805200955488746968386736420319628853608858379961146169127882406061840642543732580131859245840407477436370649175065194022573513075481842531306500237482145087488},{6102558843197855019169488745301147375680699119395807497372305321547361929640135734301710006762666925655554570785059491184946692616312174144888153078401559733/5999860266495493690171556054456453400568632747102414076076801441282820922174635102850323459750003671647884644090641587657620426327204795399064288189799727104},{9514910706718333387883094257474542952212968575920324903951003192343855578048968905075522479728001899241306259623260157232689111511631181076645498613002614451/1167584150587257980210244123657993384822144077976029296556234379356541583835174572667726459801387521622058973101912148532208350116214406097356231035540471808},{148788249785026161181280613730462578079960400206416681820891056453462539490665290627491811534622077420624491766150568765105162835840364780558584747759173687269/26771402826282218212719629399698157247248901766747044950538633125999323674604634145143851031637540076428243902418475616317670300689852473650198832355953606656},{3241313051348140281697590346205858077054359822493955227316290557700181154978852307756405935098749647322104096594212274711068189186749710594081949219355871641/3614834202245819714062958090688812753267157853395213596113077038418162520062505432254811357782916698510846648845656940358535472273454288196367071194399113216},{961552183856766427705293589758176295904275428637526835101926903154080966918303872763552425466406295545780836713810185439540411952241747946320940215708021744013/522908749287717074615543078505097736837188153843849047321336391001004844718374501466653768142942910455530280653406978118682258424334734933824361280968960507904},{8717410318335572246756973887009081908803079872500393864365808054232359823474242763478769268211575544951813260917728914204776474018831000141086143568175798189581/3842501389442937511968350729830280160026171871872281981060134424993319548405478344531585703061782661059525834725248789714691150423882878283094823987595354046464},{38437140246277751059390296374919669087981757547064529282920416170333629925330465101490756181432556106565747021575775674256359990490217899172375407914115312278711/27279972537777183246670413602397590808582961391601684407296840883498742407454125648502527087424100461251999098912627905263982514657635151526095278069126803750912},{9900543903186047990782876247276742431911177771505342994710597274730032408527260757907050274913901402659491389717556546656038668202065455574958096668627048121/3238511391634894499399304744800250159784922055331268223217269827139561380490769612819646918403530846180998222165526527969476075656687172838789794501645828096},{2403189405394857347546458588540738407977544695577429402616800101793220681482013779725118546838846172184380885660353272351257154880767196557117621944045050243/2174382300221426535973424053207593604118096099132364384480368167678899332503306283680054711035664379752610164743149761901463623053430705617796478935731535872},{1040269539720491494784577493667457244697807785225062661797999174433646492812025096966235741393556647068619695422027023815016862596987578272259587926158456778301/1030787501680163189306979145170246570068548800056597074188324185251122152240130636458324057813002061188243009570885539184577290170621739214129026190811533934592},{2240957509378794073176288167388835533096741708182063672634113235654333022943577708341924487324200797156405202343533722956382617551082582988878266051638128343/1642535159293186741289944069421212562276163144118273667977221008764230247824752553309615254498213723230472951073194813897945091469320612521881176085201158144},{3525398262189575791339543315158429425935816312942102224099735731597136753364865817506935060204441028712035087093013238676082266709842687519700084958868661161/1461401506277546537365665537390836814388419678619021207678648789899710294726686767793009735556911419145314888415662753646217741003120498853390411259059372032},{343134273306009583168063958983271075476153818433080859176001685227784492157148656706805539949370329164666262699945517127968457072387644330430619147421775741043/141504981927123595753966968639187423333103239870959344909399298059754349147934658990665037184312033258239847902699186190297380242015511066250060207936774340608},{941996297472222323030876493961339118376882031532051126122843339442338690207093294040272859073916548326894404595145364367432487031621102006065742454306880424947/174583546299952403907738486831999747554550246309614301449436549489137200512439053554819598172107950474190898310987684327928123810303188105361692467386524893184},{934804426086960227110333677900998068325151965845677738856497643435814494852195854913697815389009031685305420530847572062687774968201206303553401673813720587/616054929808175415809495102256446135337531055244470461901183219057092060444184057721543806681876211923407633609656433583210565808431658626832217777962483712},{5014410411300324650432725185534558803893817737803212443913280377246458127035360248690172116624714619040941747293440197874908907587954420889195621666926429721/5052130320196374227091233726669186004096667772034905164053833378816101730216219107137701279992978064099180704984795342048944440690216604325686083075549691904},{96888201914751089606871262398498373224763952851756247995025184856314917372164209277640937322238487519246361590297757339764064045276005251764108742407476214159/46664653438736149820111597874704358716145825014368300176217185605358251114450383877967843711539247685958043645007496610550854673707393106490675950767581429760},{9727258179686006544045807898947024542729549170148811786509012935636116772216470135164298108635951924756882635167102722378421855638996052312860347224355407772301/2946515717260002728054399492665718976941475707235139748585493387502158602347017392228688670385335901598300092978218552193700266054946438472631952221787436089344},{227109428541845077953877925419030406500793029162799177726479012587824011105299079621398979145972916439238891853950504147098377451817050486853805593694257502775/50053818822250671659683031033380611896985715013237675914025753064563986288109628510641435112442067643134034992480266244985222749279238451845045069546428301312},{30795774624186672612949474501386971325920183242820888869888243686050399264351934251754423567856872863108714450980645635101958729096039766587580408452863890535/7115649215436642847889019676438520103811016423012946763053652661012837046613471303624656098376336993031250346712054564579253919328334018035618023904192233472},{849923554903014309666835502882069501979848690772526728131963138002305480869485096547522922316977389586792817810239672389858072586842809738776176300528172761393/817900903654593469420561873250380771987096673760491572186345187542007038746496237825434151233125923085400985216834295470282202327037589723187750410828814221312},{585908024131290594374888855840646048358791225987610287745619987335229943684488519331520084467882765544568199947575752414691776858784771298804194255513711954029565/90401078332218371871042833255670305826822601426204478930655791980767834345371197152992163393702823984075198388145403492780752944737459049317196797216836328357888},{206820388355955146695595757528301229026503041320203799367359151403762777111117601957438671215568855157096161314093590208147608805944099246106553562311238975619111/20582053441384518651103195846416600075986641401156113493574463656677174861296018888739803521377744653749843695688750947055259540295688401353938495799772909142016},{379537923985803037837899305709472581183349674636876473989124718568053100859578533415554444566410440935744218338412326218515915518940409528809204714742234663563499/68448051855855835874713260058093356491040820458076777203585877953843155896770001454073548919137774310373348378554110077862022099207980976777314069646701316014080},{9670817678404179362164740530791482520089729692262845325559153970365028020905376731272651347496324118697658336680585567202514682959542142147645577166277455710981/2546860867550812395766723401142185989784947326658113730058195496482488661866442356464894039987069553952634508363457694979268572744261166097258817266335872974848},{4795853169652665317766710517194218977735382572007650495913815790452558213100989154939805277188325384572623193956862694920852709863309267458584367882330000581431/3163331510533231323459604303747657243719958282592389300087786682420135488630536094119244591991729595658825228681448135718943430173511763017424289538192732848128},{22750546931330554698072161566652026192498527075645704852516575085421778175365883778269577266090029042716311290984249373237360740088053724084810659611739165011505/14743025311495389872833999891903211150588636140330988862261029462671896415905164182623004891081969230137409231883196564490027645262251725560426336848463341289472},{34963218015285600259609892098029115540184946929712297822749818912616358771474695937152703829084454382543137059254101095765032750654279847568243805306353085833/12079987278709185131282437540539914152525599833552228795289170197595067051295129411177427918472317817395620653849846237819308446448570401863995892462323236864},{6029357845120319866854183379755779375730519650513639170670079644986949237585485730031575011208544619668194728186143848190297580347829832514463232399971166924363/537132099573007612253317670241856817355349355001699282245053274415745496754597953413471666522787806117654660415773389625778880055157424715832675144920966627328},{47003802003527066651669753380448459692472621888594930324661988055743222438653349457298205723763788050408391928289472195953712224193056961132523416230786936589/32581232246592318585921457827193994421059869587448792766680775726082288051171344621052167620593453149486593953624933878343654980574327399707617332489558687744},{676957134908050350901016994181694515355096694030976126904465054216306659537192709731390446688562577235201117461075041554865329582009878593466777526951255655/118067199138546853437116084338446125726034668598470164372847771813722878446515501217848825555896887835129393892781014312366257194846113112281838534445236224},{2792156647879339588452163353638396728181717075541123020049086728365636260700263765205614769298364944473510920269339931722994572880100452583035386169202050553/2936935898203473036819629239568042158379572237970447365238097308368868708728642696942953928432343792103449831844613357554430975706880717167290673174362980352},{5908132520692368448939438248577908747784879648476181449651297938423630966668441198408946117555589264486205337553605026780947693314360497735171249664308878301/6809011981376914007613558905096406565575290111481704268562582865437548566439332527456700549561685312593552474982744051420280320592123633390517192156369649664},{451544996618728012898788372608817944531478460887362757142590964885288242649911162643156002114423960011620374330060986459059307641271313971645933193030351691/487067288461969034407710429944098103674853440026847309078863881879236044169977752679758943672444462273508680638087264861890701519649786258705498882971795456},{149366289309795634833347753541234085511955142071211671556952113402458394495450877286458223946706853781532821013827620983530085989589446030215353707596712190623/20054530973196199607849094082413507267719381975975102970707588012856590301819157381701341770223152789238611959017383155946038961298110556743199699203847618560},{260753280393822846377567024586115049390228437923982014990338286620188265795337356980449999775107396266719586482479493527824678930679133853468279804327708596719/344276587904120688745690878911529937408230228885792442224361096958132115303996553871867920156328361980566096156918949536811234062955342757306190783875588816896},{4506766241990329493705055849812593918115869260797686048136972744661237856526090673767416818826443661445177359310597628409860305028457332143621906901728202049/1488587130963995966975431023849972307891022794314692621820967956063495019215136089411127369943455416562049960560704506271693156424532146561865952274241552384},{2451739241625173457343083657315615393430327576617288351314277541813611571574813616169786068776621616988061281657257076261454433254969903647886477145305016649893/1010529829999940228299770295614939509592535517577709347334127951380214587906182903671184640845455060495313650700841357248327032075610773660680507978429516218368},{295208272197872095806616864397627354902733737867885981839113963874215119532713888622850038687076566534054600094180525195495598374291699093364542622384016083253/439835549953115215300055197872605280761691137627162606790766259410463979299495624462099163730564325346855826740714098189545564777670200943854447822493123084288},{4836489356804756698456197301651892774624880842887793641014397268405171469843414584144907290295075952189638517624619813105951411035719665171656940513675282534171/7815049683247272073638739754704657761257370060838077780764942047337600726161204257552085496438669071694377996431189390454149797312211549118800986836574404608000},{428484982061397104696555239944328354762255576928959991796780774460866749216781066746439192014931270996604276751486956978645106109760474942308962779325806707595/826076754187289669140735029243969563267017076272664120594266168670396686464936563826638918344044528384739391717550129411647602566956497169606555473224300232704},{246716127018974821295896673012167984116847281760821906416627606028883673728407417154426344823651840691552058615982986479684657944507092284274570258991048780229/384368198665678968228109858050672522821098046467450177330858155567307885018607045705535752042446125599100616007334505827286592217412830524282791338929509367808},{1708717784988364788256671673110946284461446645637919552532162008288706160658703136283604485782359809956856475187381816656733230355331272371551941299738836556937/1284766311532190141715262635920869835222999863779265837041104618958967061866845154560585583253831615823990818464889924194820753924043664732916269990195182436352},{30966659736837408544732025666887216014787182430153254226224572469259364724219761827428884436566543010276834482229059300057858926211809462124401845302106402953/11819747213881828080836722709228316758455888514630154138914474290737315939961655601048695293209953196028458250471094405539810181817109028342290002070862823424},{1079541494012824457745915965969674448168702390282020518211391336605448267974012698762749700342570019929983770803645922924476680381747515277761193681221863107985/333853933603682632379344824532776565544403294034854393582436111457510713058394081269177861318344308351853337679129667548110424746557878241456348423302538067968},{2656674935130594900658939208607699024774100059355888010638503446795970684604722961625757468365356489767659629340881444621897992586424757930644334552109750124799/1873122263241696414441407020713964333627780797807765077383842369585953796075250553051620813181592626489851408181962094017591803589727033401476509925623732371456},{1632759775939042490769926698924241872673269214951127735750466371200333839070937020666234787260476048801871565922003007775627957082466859030893411221201437779759/1600291654547218834514930673657786507876981197532237596243929001428472961980497472903996437842979551930941780041288905200091425857287328722269404159526148505600},{669579061616897406436583470802046238211461385146036100749577989426857418973739191336007267301758804470209542963709298148666048895622835869348778108512728724237/882616459369263810977580073263708675435760781080285574450288128412396749956258722126410383940271470099966814042115963406266662725825689727202824813600399425536},{67423261222859288125378205279974825335433151664652317058407730078189108487753217268623261129914518862032740927370157373467201430734302542711820981269786895873/90777485528589224827398190559832954928119236896038108745983112162867706194052927753221642772563768498683217020668139762327689530161735239968273137459754172416},{610909126750443692255034081226544542584908902085864651475112865002285951662182506600161657204829837700280379993447277124218933560596773673885672831036583775701/953683733527876373261299897791824992159600245558907545887719379436047821408171612276009271680886051936102636884825225117184723986970900819373957626851875618816},{232982244693258958741735560711296271824685341703050393992497924460171403711865212078653209965519539361185591004368820778871459720713622196199417315275541583521/309078552371457063287190554407811086504722185398867915791560737731672105966507390366746711194963710305291345825445171962075709677742129759380443676618522624000},{930574148067290962791326054879082597780606169638969174289934620836253271633737252687007911369942734984244145507913284603926999723092460430519950852959948296127/256692097216490911471172424559254860539915672750215385247119836593160304330559636945705644450100487739310178157037365929517124675233308465875522662635671650304},{1075366529279485863426711409682000545149122632357327773743282400073383411778688586263424650068582143324252356462417862067750627761206328367713058270499738416685/1185140057042995588323457933656144008009783222829663427971405596259173192594554875049474545745845223638069965944763801770971655196059365261045408739292602171392},{142487892771894391471920941171101439435214177656229988425236017238122627270762936379674449506810478247971069764332119289599058542996664584856680933103869142645/176593343120291445722796979343674938170461311319926019668490910141118319511396718182690044590015217555825328155410372907446646877854257591984994204603794849792},{292444309621352343611190337351471486107720463951970094640199896246559281173813805706707201473733822095957536460960343079104928115581297944310001024054631881659/90409476421336621163308235728707740689464005233364288639793820490756886514449880860495176470651231020939239176009739469839880773787427588057530364038161104896},{823932152577202837243882739564282073767637658890436014342994048121417063315051565672718002054202881772411591370655610399079849070990075283839243922601496330201/204119133524171900773469159851403242107770190117570961124341239322167157804951692663152751115734489877060845014767065260225779411951991811748830370439056326656},{3829557924760333973335645197998532401910397756904908031281798089364438436261188266562937491063231163797281778001336973590188800515249631426109459892455149623161/1619255480374572949291955417389816332899259005454471789768078139865592063569628651718960211569766564202213873034459614346364701758987094181600219192290965454848},{8208745593158985218895089377966850028926759014195604943843213720864349674938257719605503078289442327882546402187895923441226729976387938473744933880595973261/5512749682999846347928800088459969621482105199249404627609417845277468362785766043782602551320158239526261147042676034226201079119213255419058444066936061952},{273401800733123069531605153636453851208210683722328876891274077659349118606810519263936924440412022265538440748628006984754058524708989953893073753855560537283/59997982975541381287533674971209197966982035828768930449187181892171252674785310902471958544341207579658057876661060355432790960048243469607717147379316555776},{4252869342620477088456311634466823100000091896263612038358742989956398390962316462496111986417373090653898901765121803982277161003285797531323472375191537467/958877618087413422454225359146391034261578898442861312427132844919117637084350514641268333545792824997763383165465557027465460748414426727846611884240273408},{77112791884592315640609384238951259865897401640509381198784637051615103160504305793863242782113206211384688041287536646539772714276092119743802109299017599267/29736083521723892321880055641401881775746616713617648530049550933367625321789433452917006910464991587422019425567653874456604768179105511276020350254688763904},{4250794170446860004946098963316034535643309719787010232050479936498562952684640713673741741023759357845299603094913974074556589799172777373775888123425265101/144484560494369340793305186447895513633298525525612939352491155800674406679058976078448003444592649292336754155252349600528795487433556527532061972353777664},{31092615999400484462319816736005351881965123105553713308726923146750014438766655888164621118867005560716039657756699400048615256928352531556944697448857979331/23967209354739380360575300287971747129325234893292977456852896392154388610211016211797478112802619967023712859157316565355020258870933862330289286533932384256},{30236225555856130232830742085948893789493655787932060530184233669346294316211618651086714253599343913397306583910964337915959514221021088963991575751558815343/26220428923614936458634837757823750269544060294420476579941008363058657296185836149161083606869963054163510009567187688120172641238732697507874829304525225984},{375136033362184599775429824568543461508888858155595869066774555960389680595988638202886821434460065716096618410155982446523904621167162606114736962601131155011/301002522844495367867334565223154838924431397217165185924374287483247158622979219982719271643302877513696688309964015351202911695188303818939627134090620174336},{33587335197553588646624515408650154739368512408803168443347770449186253605370392122507590888817903826244760374824300677260152329261386249366718050406565196697/41451258281037518932094889643242675883550413984225801028631142531429824814680880067068090723309917631370867532454837817470832875681772136969762982135517413376},{10968736577255825511196060546632968429977286319121634574169807458748245296844581906688390400073121048103698976015388173221740799811294388811970350180055623423079/7434356870860803773165658528344392753607946268757181829556924232625322250924426720791542441512039271998176600026982385729239678742250547148501007803137952579584},{10062706066885371655599212840516985599301250100586222920304093710420250347520855638727347444057243745046456837223820992752848167409832413688236220650770442667/8702280895642180298180118785230140849764439111422509655479133907839139218418799729982567837123520698071525470552664642916165599660919769359099560561756602368},{57809835985957218449226227012803350067284461164806582814609438430805111162964529255538847488375258576998254041586155179054421326055254982965926822381572337999/41290374045865001968307486764553056232445724910463926593725568084292906344206307631191424516767589728347501661703437690474071193476863053156569863946827202560},{10322364411836139556317165919590588108767921436368990144917051441255363053435553286741934114288390660105152993868470506333924160480477508634003218199460992563769/4759788093258310378837567203736298579547696153559404923925777490155870701723869568837885011968687652650279529382529455257432227729941286189202404728442017808384},{215202470136677010499796987558461572656474262386939778245575200410376385288245481097523115813374111023732637202463270466828994472704692349019096710952585704029/184077000802998082067958439994783891265399489580472982956816276401099946771433577479559032061866280592839380648511928252549569977354412930983199725347959472128},{68449503451437470981101314355594553078027433002752181497729961371408027221882434634766160984233018761159119637910190594153164941711261354429453447976226313129/49800709005678367800478717481309926449460968681127148381690635537011147967863346853261368528496771330343291114206996201617050698216049070852760950445277446144},{113727699510804437613737748569394885991128657812625763471068087415202347992566729860689557705416867207045705282043915896614089750590224991865091543148083609801/68351696088785059365583739168495845288301118165018345023996640614980243274326802280128034174183566395648856299500746902329897543817726846643447455182353334272},{295451451168384456128392770012546808611923788029017523016230164512067927125725433132224106876735376505750770919166761907154674799135983527212712390342064988063/223264908682161357882145800087812257637990692484301529479551868450956708114434706592284651381731120727029319347931793847912322132783127379090475814606526742528},{202367161797349636869002586189055910293874077835864191945773753019290589898987789960224995490434938726304972163090492837668940882256794718829983095681047087909/74007252583924668973924388726303486864472685646388364887707042706880208413864001002906832498009017979401347770040869457001196683136366626640298860158027563008},{10633589766399888218053207576965248564082505928414003580166229875314655706242416764628409700154301968832177566623838755923980050296441975936294177947697878023/306475431605226058258056082532650315008434636187302437362232289134785769176078901811614885662327971459636352201159278162202647728822884919499738526917853184},{24094464036900177820530118432534826219432565882206320906395588627322985731717289997732680011943314287657712510903834429568246388652490926568192865643376567/29748645027816060898898833639311246598132277889912897000322413912525848507372184995895097411559427846718183551401628165114684846080586941367429857013661696},{1124649860898939917064039272155353868389241276429022198302044170680629915769644652627735342246937733015451679677942729129131689408986808912873942391265585/1247930334869674372043119175711645099261219438504637535241099022593352689140643665985644691990462047670043268551655493577582584646508949479428171184996352},{15660693648750821253325980698960276105670206155425357987523685342905979954606707145794650015374847642734475777783450912922482501483352809890994545867167703/18013174816721410101395077033324871319425600197199612802291177435707121355776273629107615037384231159768354224995006711590393878742005041701867338820747264},{50711998142630391660389383178340612628613050551441799924366758282863357976310093351308393335516784894219168568736832032573956828586650414089549945976370810537/36701591316887921098170529556488737922599143495216021263293146900220427283095701254413636291282143277916608154109094553176371924157493465290890695251337412608},{36114688464050593188025947141897555182364911369853829924338259215439426292795971928714435801883724317501233896397202491961428579988436888161645541314821484745/32752807880148783493731543174530236595054593350718119414497472996315059545628355251907492804550345858041031891186690847186154287047821286723947749928444362752},{43517121296838881654126364711132963670889727638286613891624032109750478632341355896500768496606552635202626039045858902357251651994224746785954577611928186533/36685418689635744223706046171551720935674052923369243729761555685233879520293832846484920765191906902322146336074101617278520143343477942706757680100982915072},{22082256013930661122449119697250464790241041268492641083432135523345802157841646635033154191108448346977820742460713128246638377685955041137920072957061735257/14906971815122236394210769351686135916818677421517913086437974021143189167357158146586736446455668695835183788179809349790654877941166010153779395439466905600},{5221021142604074360171896735831375661357776519661473475681212814561285438899936848758900837874524514119989448967396811326339351791432016471226396483202765197/3259819965637783684175312243271443549174731088800604467815403544299624661127793265270181373660177071392684858832197763543206482847792598999521283417225297920},{19429578322168954776153612650191321976336652336809126715198403082912307255647987872172854799289769397496657626258914376962043817518451543075683363890069036051/27964004691268548772121412537892636189450109266585110226873702933905301826900734106783009331477117104365499311487010463005641582197974322192431320818748751872},{108276093058679115329707062089903905519409686173458769282522849408897790296512853128272032874653250172353336099215643259506580651396723146655586484954212899473/121700351737853335364282120135965857932377664196597955064998108927932891013678698216936290417710780172744371258474883263324951387611287008748892055885446119424},{134108827156410543376447098343184930270881948394260717822677558360450395418542574030514346049459993875540056895075424774201516924013397646109427753208775338821857/48068052549219629077367503449326589355606219738946601386219171885862186434468219231981047021033026370802906443941245735209227625763578517287200430739314901516288},{37210330665445386676240045198945165314224355657273819757271656759712969235293023760656797594600339579817381259632885759232815424196778765413103771326395728073/52345143358513110056151256276317984743719938005077809542657220536531279006559704877352107319906668956739902948410158413072549812941291750674858620567035052032},{854086184430194828291241884357493545786181183523535056719151990568814137641514999296610268503794343231454704245716685122485585089957891856757388279871485232053/432735972062910243381946806781215727991394276813598576374880436016963802821188246289370795133734014532502816720377780927790507114116653179983594608270239596544},{4998999919413295117993897555914686200258480845608226749902826312120543630287622128650993265380677663015192446118837486178859299258388797728017904841731269567/8425942725635983399906478396186067911810282977224619880747884252568620618337231693006144974338487360066381429916943534592059256524805511565117513689819250688},{2058692740950957849160219070494865954831322170224108323710617116788709971282245503527886881625697690335414781989246606739638954884376549683499372158012685875/3350269467633986986525321268666624695611175272645044535647866390210119014051601871551608423195485744468563413587732466798261007929723842762780553011907264512},{462539548793356113824757886672653606081799806810963493933133240172164590888544463990441483594030757260577818395593426240845063375779134637583402873351586843/416392826457722393091090689692045187343019336640875192854196489960857155923690850956871468190000711155329464975712323023661831751047325887136918740530102272},{12160683706973965706604805892873878972008235388575656472046228857547749187624206682156023431077002492892009796422599041773495664437885134390197259999779868591429/5549302793086370206170679625508245851590987628646286537384020118899629524125798912369106934566681463401883894148445559645051077145346902107407520947190779346944},{1912063586743401410572332208992731989451249737928066097373080856346246732925634191311231994466849722521808398918655072540558099004749814422873480585063275689793/1456432899553171903518914594032352563749215883743432564991758252876829961285142757094949836991056320507387618988666808830262178329379818984974953102932541702144},{118831146067485680872707759840090942530587570302860819719811150810920465277424190790741085228743153330462919437388571836454325539507091215476032951250809765313/195470537871573767443483419069953061547309422077342131118605448364407919844164165701637613597247515772108357756609215999066009690926526950591717284062366269440},{3578236809917960135099304997334942652668817018322712332709984568469656218878205549868371277357551344663306333104656748689563868978432031021887453282860915367/5139410147824535146650527294622099377954966584898293718080556352183542786521288614694142097126904911177013773547532234066135415004891332320386622009464848384},{3172761758592477352615502736499878267676719071034808069238368200936508300095458268848141368225150758960862798197858123280646068355918019572538739932841354013283/2991113622220258963358556024518975390721308564356325573069167962127113688200844498838868938330426031374653896308786092904775671877523839171353611235825203281920},{6408579628532230255551438681792553988625088480547899238883088265224720259821312548174932588220896303237054333805739255246531764077913793741822825818515641667/1879659423355235922159860750475526517397540395998560011717115186070794633391494150652668051383374677729482270411896206209706269279750022521065499039084052480},{34693286277128166073765062134736338818181179646058446356344949305424934479309279029985328612393168269409356397539925679704749552597006707460199979261751607275289/43645121172707728670305842700541663709767725653182378614072286374483936087586510235271767181207857200164802525608315170566781038315869679899475000242155869437952},{514147698750330734054310194300874206124813385114269139780209809742784862687136248190748131573851342081362461664694973856252444488969081839598893275462789006221/258439603545602220911167772752795130403437256639223425200316413035554059671183564346732629297532008865342837900036581917324188005104167122850414377090037055488},{7748371450909969342569318182251607750506797302540821042481311702118914985793451992917265910212356070366607069044034200101109586837879019091588436809355308655/3320487449224378644378523736441140127004558476238530901115332264075640671676442410500112697032573938759909333965335807036358984018070925748041347930030866432},{20662779175294662196302811386801498035869086327426372908441688015663557114405979633915471642203151463827074432183918705930637042890299286081563101038416244753/18659135985661612411067659866120261679702426801316650348872330766697503038779289863493402787426243437146482038693410600516709578764543972166760101823496847360},{40208736077746134653408298258076915499744208561721957247985256193106833798133940269579172702431959970608053555065399826388348312274944495023018371385315213168753/42830616988782824553188599821039577592312554895893646598194536474220465721430996122222338169807059981004610069633480734111897498036570054421034556780641486110720},{6454476068243441737134389888135265687112265127652156416154010361272255499987963108254923206227526221139660500318525648608091959416244245419807189222413750047529/5824801493024579665418460411191260767406368701679146533428903878735089417342305778900577404110339569977147823232866121298619714036501079088645034917512936947712},{702274409711721223256807300565769474052770212451164181974255231453246207237067930392197928393853024557338636506052160347866969275015028751750070408494730609781/636832997898577899698264548458323124476502994158277604298823411355619765438495870670259997941633708212335450004294971356103833635206287320937212492343502962688},{119381688296636192839397636916635473289768966561736952577398707115078020255555022621189296496131967099853641245766396565179719722240053480181904967553037725613949/100278108154515342956454461391556555731324336940384971624228122679843179439739373921431071692511837312563000523951984507869496974075360270435742518552347065450496},{793712391949725727476141782265534813007976332820493840695828697872353984896524784117632826782594457152303001494500772485118017077488004784813592165807046511781/1023182385397978817835800308329725943207302033726838507455423015940895915320968134972607110709692196614108775161014173491987793024437572024824764919401414656000},{92477680652524088112421848639375989786594435603918559260830667814494816259847342292069401327640459330523555098401195964251535804939802236639931955967464632763/102667968151891426289747275877927632192336415174294414405251268124234209270682525314923686898160319744414953880086470786124704111931112966812938511493422383104},{5454827372578355721883960280225697841420637305199538528812055414826058478115148986184699931082131421327128122581315915741525925525410997017953529720309892860521/4007306433616295406906545010540277529839748232180874223923093537936439893601656490042453025522054694884149794533812256256074120370127600250655416415918665760768},{3992529019341968962360438091224989891480686121729781232839925228691287117914486392757372950299867258832173221049701188254234123872145465100868174684331464705817/3288599930408461279219860123719224912966235962939652510306746132954541483003808094751020649932725838733005911931753688980468013378822844203970290576441810616320},{171304977047428351704240650026085166908959489920319956272157000834665918014441950654621187363537289576237879303892850650529436177963684986799084832988871058017/88719628574819172533287919266773380790676160189359065603323722437305899908443441819754471141764905574738142086925644737608504796481621167625616724737744961536},{136074924745887621204243015914705196317091482612696254444215583017829648076659466598220840791651581899408533534605709315081211617796277036926547973457681861/128889927470498069510073911706284698440966123612491437291063958152871679113471190951368190476784485608541060287272618805960790273062237795491781547206377472},{68638189235249365338973569285587429141191616198044173077817952237574979878481656799511368620259837381017248219941374403779356936481779558766414425776681350527/28754885450469940850576713691984526042093564680477385022486294769581238904059484070494464275073562910761041878519774376187809632294713410777968591071212470272},{163280471053667807643084698042292203959815273985229650506016311995999691499926646246608654360009806368215499610416767910356782737709349421674484354055726797651/222357157522707675569862498423340822698547517876114161728132344564570929434836182352918684341250954769926014460531163928078443117269340943187737792157042343936},{7028528194378719365174226370998373249858893465307509608486414717021380862947994547871216289857192565581233377498626449554387076942183571149055711175950242097/4637886162298244279463689830688060060164241749080261869095306535025327334717624132840855807985273126055302062906722610255569585484243364913743797362989465600},{1386036164175157002262775513034476461806750426194618113655195105703484763475133456768151109188343932399004104334614445344237443454891130205956747583872700097/664867124311920758861612944827391124432875406934652493423202009635621079311827367300755718293651622067373358934712542323862296263828727354735767661607649280},{580617813278272791486286683850015317705906303057659788375400147655219312393949248132740445399778883741397477351247033353202721392331292497268785274869299302073/236734706726860450953745375638592417130055679260342879588271936440958072693560783655224589387150354438271544674049981885231151738543132533058045509310062002176},{1260237569087531814389755964666719928999822043299066652622982782659893656979066694171431483151692827104736782628170597587772891718687323686412782169540194995171/1441893016748891781991058873274041054358781917347877629813474745675685942099303124422805428961181423731093406065500184198688985613837011844035597090169347375104},{1419876756873286821476438594612971402008743225553546207082393119886972353493040936229476039466095681205178520388641281879624798412762658377949202907425117857/99733268393028319060760330891737850056687535674914808268587573450708848702663961017587128590074800153194907017437800795086965159210083595096757544457076736},{13043179190770611208941263969610430845354952627578646406860204878121241432297586553234609150233418630826103517286852829180506923264004868175018789982146722443/1992319662345475321250911715219720034635591281331622693834791752038658311896890473305682008966221828488321665926107731123214645633818302337570152229891997696},{8366760297262033667173484284124420549425051052833721987448332850713557022050993744787398487365986008761860539716292438922949145618413122719145173769019646987/10946659076974251710349388049625074450042242525244231227549790869155194040425392033482289365573158715269101125717080180595135285486470627246422880519815430144},{2756375954606422772810915456244904714270297090251369882370115067350427215093420417497386495035094325029866487717564650728735842152347059726632043293197425831/249020565503843192054436097039842809137697947216180767179356075062290437163025593888954729199163737803189666621469736489862327579940567563431505430619095040},{99532020252651016016221056883361260686145873223981582029096770783724911222759907510865451293748935409174355125670967430208376163908177391867606878742086720041/58581222468311254994682386909513002751476063498325935840944761951734102334846068549325669432324131290367837463713269036038943860918164513766184193649362862080},{20102851316922787974164123176914176798066755234294997328920514234491444854886435782862314265869641334883492897675725742823319621623807352229457778370819898943/7262465646748797939660906873502700820613594167130965213663708358884111606419644376704230515378134992289636258921452714914557857278188689790236119386967506944},{1463971545560305352171953268264732628086829116628565922103537680230258006498140532842193908406785953070985888340623456939945499119101178335778321895445343473835/146211908426758468643970307770992439615864003004366678267608956553988945503292705693947248374124896942064318634330501460422530433998101188765318535502580678656},{16927564520530088636182953840134779374864708945462105416314767794873640666099021321000185779743905709937486264216799770210031074455910346917276423820123179632425/215839772290722890293282377429410931764439613312758522899363623709306395894473934738199535593500271879179178341497077012928788019662264620788597872714164731904},{6550596747027444718257610718401819392417781734758200225744207061641584288432928821903784407750664471350730893559491922408822353539333375657112981794725422361795/2098941695867599757842164616316127438111487494384633302069896629811904799847318162596826776701977438525150419145140795299618952618478666708044608259895565746176},{556413962780631257865686165834202232792002693985095441775236899896227330389737789320500065802529394675128065056812875733224397813294122827716446259280700188353/387551058014371171736699501466150100158024748673293350354593255762301601743782779625782221736813439845990135751986986394987008648796086535211135617343256788992},{79922269148046806154380005566538619208325953411553381743709887481065871566417155961869658629822931037751062175828722634220762051819108930576980090864161527309/46147106765542529950592778858556240099284056463184795852461653271442952178578599011251570410995666150545392833261830457941325052642631833486014258510452228096},{22887036823155406455291777845599138190420548060930611786610129767816933981552024768553731695000592613616825797901267535854243789373623608593578455712235912575/12246381122328116360152328197269858964381472883301172657207211490746675057932331795787133713658180181012264153870944506382253297411059203503217382278167527424},{26725273443084412669402289192668033337397297583656135357548405117675488326947690846293713522235070792457150611508658885614144663949826486533357113236978722535/14970818110202749090739377544159589595652383010796520569319158271963869841432527005717490389817085195879051243726033429863738604049886172371304159740786900992},{6844296227523539985217404218777731165378721815150920950422624423958918596544602973500605902091941625683613097024405418586970897705886735253210810428634834297/4092912280114580853629317445864911162103592634043231195724696057849351350506336797286140380681213199189686383256029474631127331963025611019563921715721404416},{3298188074379428907618759132027821768634092002360520742505459754142211089089637887608670139128619055447323999977443695222353762082612495100622569609128074141/602280595865902138561204217209310392622424542778529879897173908449059043695919334048929592233749494401173006412525032514717265721787249450478828589166886912},{12352651449849911089761210949434421580199504914191669955065306485804319389877613017613881270380445641932128022386699625843230412384384177893363862031382110373/17654776525420158882264025239076108445157211689747279620073988194615374831358649326359305584106117477487169511233251859990801165936949826410417544601270222848},{1094162495900984123641523237171050745456892175732394450475373869046581503727713830835150638578682307394851396330600300845574830479800204159336111954982668601/1298381889103871041100757950281040687626539364786207506608117757949130603063000532957496908322935713512713245717220372679210969819536499931685009220595351552},{148575586391868177503969392249921975999507196700498024454005572996758406590355304384793937787940207609562614202340195473611854889697660766013782472722768308353/114257712468182636148619004923818729638458111465821037311837882035833107842385273900637665376177885369811797955130950241480700553102721898675387467706855849984},{74298171044739403039999028595871655205780871884093022571968879920184292684022162411531994806174479263066786480460475428575601665443883793937623161327872273072795/38592258980352039043570530985780273301060747306824623259526338914938521743930549517734961022453938070294287268506745615441679004001497111734136106755918856716288},{4772013494657967943492838270131813157315084740759752667278847252127163856374387245032236554487074507736867103523605956501331563644150467921631120135826625045/5241902696212881070838525205423131856202631660600898537167530794776218616293916219970071369095087202159477693041314293521148901311953989761520421534722686976},{1996581268751058425238446044578252965378087175248115629499142140183566183922951755991114048499586966427351148386669735714997638451361260459795029653583007486785/19966811389648967450062039272342548363506246332454864605471147621988386371749423612408669904799407202320056290606433903797002217133691010531493252593134600192},{193152668296513150971451660052253601836390715667066276787089823831484933499997876616885009524216770732188989686730434615428432412144758648617939213919886365/51382136540718758640390619925695089721280310303167100299471858342378400369031518761178579943756341117674234636735716748019807480387926632089677698835677184},{367886678263110297451007684443722537000435384358671587607054328988311659335968336942261625261269229661752709281672104439998485659532435982556039781237923782183/486095340845025758886782402513291760222267917575766902631536889533035249310025357055506554811778302278002738411440551700391293107590918889958527088594497044480},{1363370827382590885416338306496596200622465273691102554695001801088109303392381991809958634033591849040315099544414722190960018049449241111456591747708455383905/1767731854520650028590145860251759199177928173226167608241189430095747822576123327003607067717866336300961151553170805550693667745407298471810284170976054738944},{747695424467630416809475109578525206089428044421632496491173298728067290616207339151142251404996859278539255352992854425749299893672630039945491722602901192769/931769655330831935302348513934926578397752870118865132913257931821950934986991139291155184550417580954531395988235164573936309231629249964982087737214819631104},{468825664488038432428349405934210283876128773701622667156050556227399785095643773047447844816314472313293009652101139523327006294138833875986310546772741912943/206128830117924915777332070581977374126388104096774071475430901992484823834775255461394236534083366776325575112990997932108128130095226416025449013972217888768},{14310104532182558537986876788773236723140687563608140145566061653532515918762916195055951177819857213375566701991360120052545998687231028693196544106173184963/6342939403035071197307771870450493916295201358806148358733285671404243271778853195444165232068193062211125272161895984143248929740185137505512335294710415360},{163381030115009155867064900782376125411166219246685774757335579481857069181350936136761987753183003584826543189242424259556383344474500583380241404073553486431/52345723088482671936161590971184434582684671764709831720085409567126924753773905622094549287655393343688284116599960203235943762989912562773432025617010786304},{9351199785382865581551917743158442551902671559564594495720886697399011881386327053471162094087823996933841137311112007629394622874225817335028702153241524863/6847161905649790812637108538246775445292690987095521372750705556542402627516690265099535095500523444095584935497121066058598259823700100319187391186743590912},{2159557374787528393533901799915879021043304844169180098446961245061652746205297145672989682680876536135336288681682242442168663036541792134352036384235754693/1650181598580324936481470809496231829678769403822582268468145767828642831348759947391159594913773864221065793160825182094809342408408277795047684476873211904},{27910771823233680326355203019401634087943526779039917965856563185272544807304054917989900258891843659212715390629870710442078130437537036100858441910678576211/43758426980203346913021132990654910453375051913548070653356716232913607030028183724460064178620010246283398361299892221107656157298016221411033470816273563648},{1308894620433580645669015734771456070279446843286649465081017777528795794247569421354027503113632826011922046141781954765596264359100663895073897743815261509/1410709303323388494651252508297981538723615803889071781556793329813909100084673201855824649349650049789803145153272088267293897157401048026169695495895646208},{4048217044458965443780338298915581860933465702447480720944823179279544993504941303086410509855623449282551556474195869751887313822097885247796730831554465653/5146597618924250974614944687476969860946783149752785434142979422494526294248835849511750037890761528664295989444118652405176306668563138371424671262788550656},{104856453915985807039307542548887123299301764124920274289745880013161616264553225147697541981281489307227279068692330165957202706581271989103112903630606590341/12681551499236395061209604320486884963431132030453165161455169884099843746999561885715214371128600768300606158926521820339221100817225007366916671844980359168},{50617223167797986171799604297319332266844000495769483560343071634878406691175306294823489032904346047254927321997596317259004369580872414802593108962112210069/39183326250732375458179572295196881290781950336126094987184629240953093224624090968061262649367162837127454430895436523121791422844589194776599583301826510848},{28957300778575281759938145045970798970388488594590540478292013170872842403591251594401707052181312500904289692333169617770883661673995881268499481190855641047/22488902355601404672797869293648818992075703315912838408624558590099604483717114966048252477055439479950486145707157845503507391702120488123452842182344441856},{72686701467983029727103920454248884667399123179267695624662580853400157716234606251346044350914822681921056990971049797541058018146429633310918489728777806957/83331897441319550659001181058950728007234441245590877356432068326890943458518855015954937064475328320859160872447776549292718017888334109747333866567111802880},{126589980552010327928898315227336810511759246001290688949114445659633421709248045923287836429520470163325459193010582443147863558349701954726231489314330566809/52987624615940121592265468802784003867569388257869641495262550825894336736001970717700595937173977802537380103267382762450526705240076300978463279628433752064},{53134302048662752309294665265532401373541775746957343398001113413147469823993739504760905336628847663898742287575009272163560357947629530986927668985814279247/55036893929809398916548163279512771872227543882674054346989245740300706458423865956134627869902494025022472458013381885099261538179778832986487032024504729600},{8944178979638556684024696161817745712957938031423214170593801476687811777716991701607590518129262877189085875153426659447956850016754478396142206075419829541/1791773745723968847054925103018893835762968718263730461195431655935041504864703277654308339947378125310414221841191801817622031979924998035489057515625250816},{20984028744254922812910132105222026342429357264027181669947913881086939314593081860213011048322086311434236531221972516421705295636475881918432285551563889363/4117327552052305658330948156417490358381547325396740847822924330608001218238261416176658883220836287379192197130668081025860946020503832916710933393309696000},{224709258793077545696460479749269079683137589392524278254650030534012754332530668828591685026272676741801447616503844732857687353194188458410364561432842639147/2549235927908214458677572001426340591353419283581687084206611750767118436489260765236821481288966613319644544915550947302744680674853205495647739195426865152},{4680605397312743677627301713157712599950100233797228815377700267042860514153065592165483100843142130405475741801543263942211408319598963991187281318436369487775/1201582408868568794303869248995596380662007344778043187493552752614212376487859052598115103559224381156325196553393906509885579768070098751409329821277820551168},{815265568858553876972809524722434769986623726188014496602377360040461211235942235740544809400533419073235079893345144229251761729051168406385943304721594392245/497525609654880131643849931795872323314561298732050320845180201449604752754551894013086291580821993034744366605241387726115594718783423185064173763270654558208},{52819865864419925492280478232151081956116550620854373228403134087344695171561636105106322172473987501145780186784134896510730287556499490225808600175328738357/28690750843637536596257618177237608540018059313043371454323073315260943436443778056833899417555364051639589673022401053783216889550945081298937604012414861312},{1213980208950859956335150492834435762283674102015193531109203323993073717597431317001666923149610078556187378876358642598619532514752321701742945427391534473/1740505490193884729810857984892240127021789267063899618679883457266533908522522032305790867533226055364304011273083542819769746801111854415376809488602365952},{18034217692223647871268109644149076817727022501951175497145641617315772488668489548038740322793522372589130641241616978979340544392394284547088465217541252977/23575785128185278864846344063217369517486150699029124665650206654616481223964677830789089781190630733452391059985154682385400107354415809337070277667495870464},{26691908075962115497145002842432348831078695112881824170083233001837104130815692077518372816833153327260850321099559757351755338657263541865107799242964174711/35754891143205295551150182460366191489121473029459793861052083847878828928673552607304511957885347195902331152677427991779909589654524766263480851316026114048},{5456453898289570109426787633413846181971450106781297427709721714174511361974457327734061666498078651648437895087587639535699198372475231219141243117124839133/7471777090948783872598185019632108844376857147891645479912138931257634265769781981185377433780553501167598545567777653250038930461670231454784679287671750656},{173497201562597440763388842302102479822684831783993757367916326935848709282963134575506590218166543279235834146174687748093298913577099959032991269459426471577/199122210843994361319432962517001143267014085663588602226971186381001829999698760597134681736521640458833481606081573130116389741639457644327348357205550170112},{15957016771669429726424831476998141417028195172362174155323859807299241361532200068378833964539320480770527839714566509759091435465444792791409944600427557573/14112497409750041353300935035740903650479903256018068960947762862353479859211924192411576400350001166346747786100060938577459273807871121721894099993441075200},{261596731971066870712321508360238677072202605952618532775515140561643535750882767955675166915406175559983815260568563815292905726774122540218997109094946426933/257636648486607844473601639703166673807642346458206103283072007339213084802425468942012836213980743650499326088175844205275506865730254979980183805527549018112},{1512545052457188725668288318883384197273660595553591446255455657585552262486422636893931971356311535184241345634723596119155505554120401257548388124018035101/1783701388408279008009851849931117605157847487633586957399324104945791588339276999080142511876216040221852011105968609903080330865475367315980320765164650496},{206661908159666369949421890480965205880306015014259073995320237075195211978510244220693289081872609683003682230312034607421590323463705368534060095683778987251/167510071061021107217661131573737858657184020945407795574279160517230149591320476142914829272508034155139215267743598244661381001681088382753210594399939985408},{1760824886843413086579410674083022216598451564349782241775604949753684992139395036121013159922602037422518658150820162180376234591238212819771571583553732909/1870232970991059676269313668747311383460156439337128616735724343222308408202727353723907968247170730223276477402334128295037050889293114474010293771750604800},{516461359515494823417777516667622194072305325483510782347857199071538032684565018062040373113167229355057416583829765019576882470325715012242324100497909371/542192648624495036852865316149322831499439714113693225060124992052870739695554801736467466048270078252050954396042634886878579724430065528586931787765121024},{19981391259762321488419650942748118013721637458793128554104548296300187934250546130120039030903258109784507401492541997814208231666610650021769478239870231537/21938926334895563870562197576891545843670167975334894995424852547632732780700554726022408821127066930304345254321780120707701223140996734036318898998741041152},{6057510254840505112774237460960550912923829130447529794199221091163601791403462970059392447403228599719805070983890030918605416257196653714154615317434347037/9620938066724558322468255742055271261112564709844296930300897141039258764017522392525346529861519421779842932597678932829256134543927790201854750705716297728},{23731193830246769845315957341546018498497102051672520417234279745710500775667837691769108147501859893772149550827775516319803131285689683814652929364842916557/23068029654251462445994074159233280097462561437664258078993852767010470100125002776966570374359452042439237759225649107217799712572562549524609531197863755776},{79105798450504577522698741996894259455803749015248388008541274302701945144565385439807394018943783595965235366992339576411222846153889449690438945546873403003/109344761174664801275691626116234825380117153782036515401437845138074997307277731586672154586244495208789334648794647545258542052962686802027716433133728432128},{18075749779449023269563434194797641895177641750956769065705020269152437358647378306581613816319130229162498949037626418979198212314032322264690871236239318141/10473356786478866559444277953979422287593810135000222520803832442868612169617030510689903819994355525418601308510230892005482846167005715885697396751579742208},{2122661613834652939783051183749455093774978702327363434297837638798085558318963821906376758022688991219552399426957007284986670879443812503149409841295625479/1923991096122859249381700565624186823646896733036077854555915384295274450315468933666045920537499084241411804567945305624225908859739378174127707443530039296},{19146343618852531777580075219042597448084187852120460431288628707921260769183077849990505647046366701479512389339150856578639395932682475157602234091756365547/9421476255759393075540597664403977273877820352967287413702287100732054255978574129814606395994520617960877168193457858862973979331746647229159427303677100032},{19398774168629736838005353437184239878836884025801095066896283143397006537719531652998708288509173051240024354738754652833779294938726787702470405025641318717/14236123044676131422794891194789749631757862806365933070870261480280441829309766921298619682148349121704439159863804226919542425467893993054497393840647831552},{5067739655019008865110468064356168859915506834321672406426284284760237746388952942653387895874180978275231500093648814940129318122982963145001664061736742803/3676335614361043697996473609525292605829889534832438246734807172650919235383396749193827781322235548922190731401229653194165859814095526382459301662080106496},{47187806933813039372577446881153085952466808799509570870274057604285321816648181426745957669039583157047987545940850768602829321091947852927126326843546435235/23882547361502182543445569767966017943277671557875642240918192867999307001558477357409769337193766191001191481634320643206559622752456581816712741621130592256},{40331429452722660509749300142989488688424364793468460520325410553550446326859805819521626057080037129663618743896286533337625600151215250363968849738936145777/43630548831388679585929135496460322414198753972698309209080124617433823096285418156627839321303339868967376388558968541048634150336031083269842788674009300992},{8807904342630491066797895423685290576959823836926196613703391258035132183808899709490607807112533932552306415562769212143040911011679275207850410042731344149/8637964881795079716176184236366897501238856308532809801305073658561720211512237057918113296474712494372928933661023280081552772468245332782252953429264564224},{3675855505759167812898574608977412771951725874089616523863274749528164962014546415427636300478594695115516681016593621980304936244324678339683790796048321857/3819672089211233692270747222027394501051429836344630288297328229180473003173853220446300894694056347684389395129579351970694786387974428984293802712324112384},{543163875070749276958163499143222002012504472630016576274923285169207332723633513740231396083285560933287344057596337672360780857816723424744168286657747613/387221584669520860460778772084315011841709765969002953145877842311104014144064123564894766097993647740022622133845239648823877404291662604569267283299401728},{1832355175844750126694860137832799397642650375199040477926193788675703257200048757196999304157048924514625101055770384794498562746003046176785567354434023123959/2826344448074528479233679242982594198470271445826904705956697220580899931297409748898866337906759042731357721684210019484535449375736403400541508846139608137728},{178293559916764635119003334073214697222738616869109527467815710762752951242202418965585475583881784819784991521540848717679340889650889095816150449012161243/303098652767064382862365053529399406487938127615974317283973364359937540281549192334057080611471561621414087896894382272841998563285474045038589007330541568},{3795863992526800167285262168804302481756827026421569947175168410922037598668368435671814103864041876569683496829675402482458558800004087228805890209606813739/5842082092539110603116699497924315960589312809034886299991960426194556685444990644102686344858457639674399276045428400149323060158176526372428898141986619392},{107909254214645803393333925938129600709653289646461161229499338999175438222696923851359166443816585693849296899535804869755271495395652504073933440137046343/176844280837338620944463238641406546593858846193922020987284837914958286780829208969391077187844102902105370349368701462357404142678726469782611043112976384},{147674627942117465250811312189332161411386312952128537910637186631846877353920747889690264030509966233296230443156790475198730373887405141307261612109342529829/190423137948540603488695540784243699210518249856095941633872676334144119775797504626672199941635048643880774668018848157406337060988147012716223730608446111744},{1212385986141650502292649258925982761547436140055967667460092251027792463642438818108027578313606270364136908608984999754662640819146563761728968560174565411277/1462630266544474702884689069524602252943539654676077479636346073949730779285316266475249792906672899020227503716600567966479336160841042889928885001395455918080},{43446178735918907132861185338450782055413356335904080391692294669106314224287563083104665268172671291358038660154010847259687464511922106904854524428917860697/51801873230434716283147821217676976941873819110313725850403557363442818046536337294569486131664719117955206820665904033876994227198078583917880326794399186944},{146837609798270634831405907684745770882576568996199833167132844816610849103634907744788324429857930363841229598073148438400761609624830147761751087564826359231/188161255573830496133701824870985422427705872063597335038666594810728572519575816889787018186814023402079278335965264935909795630255786592989615361768279769088},{306364479873245451745512769944571324700448574381918272900779386596206846990511384990732002713948193865906572273822047104729978984606331222630036826869884984751/100816593275917632593191761180962970201914855664006935505749609179353544077114972053694248974299053423281739331620644196065672136611628589542851333521261199360},{2062612033382079228687557810533794250009080612285213044283053926653586777662748951641276326085363008919311463209037503233839121487995450207765703026323370330007/868109225240527077919012127635780009929518982761603753542089020767948202751506591818991898659465234262785916972175844079367270467022635891465640968026056556544},{51455211944418867577494513721195712264534617283509090031295604028564393835069008731108805530477960290112110350947297917927777942744772316676107144486712505992267/21747039308032159198479917956841751946362813766191306742643709756254910676001261143237803233921540191186945919908346372306982673503733211617987199504487941144576},{39991082523770859667039143522880675687519586714116913434427823323740016823915135682241760004384494247071508637521597517532753878968022223703861770198374461157/48020963562826381089755231176474219699524039010301311867678052525262343980300380098124829630395263250176350082019546389478166626630595062392757692288097845248},{90101038718440364977827586252017496305858389673987311847132611718172777946206096198436176423107907987267323427940811499125039373229484747447138828258292182107/10299231686203112839042573993676984038770135704058239779436570131251145663387136942119915956798882193454734666764984263601045935299123894427696385203986497536},{253145508423345927037520507529208492328865604838000351162410699640925974468238199916955164311944850427714751705133681580013470895968492699833907491863149411531/37565853440034362151261565774599070213710484619394992011944277875536572131570567537563024146644072710564981803197298481661960456426724806434145873839354019840},{362961292881682918162356916716198954113289887338304839678417724779255339529659988584069554062543935884899043560289195458503614126077181395319266710121391992935/187598801566854847400420971732736018499831275047669003937071713997808171150165858119394440602746832087721858041437347434050703798288412939876159832406230040576},{38076863320117688761821609995681066731373320516617157355283181196502699077971230952118305053008785654017060063180131987531324913726465814319675389142746091768943/46548793506786063858955637436332870783501594118431750658651137230918438016641222739038213673240413909026501510169241719699871371148093368750101846905265235951616},{15580627348914853179926186745076575842277006572568273901387305315604272175948059836958802077777026884593235266609602632058696960083544592296782063244644282176495/18406210686839803132649486248546604532270316667441819020073452301138671747551520391181753040988348332013240953695415525736386277485195392228497823408853073526784},{677388345999627256214520815620584948442547211718580163626143935098251414143567459809641265322333216442965557363601666475816301460672429343457316623196507792675/981752623020810368072874703070258473604152375596712616677968005145162330526122603246029423110131167636547472065591477631637196432697915704883629364593119199232},{60693892491879228917212113608102244320575886985267625947292051987458698577833499836882280293409840934733879533643282112623068927782548294800696866397140226217/36074251541790463319970743474687912559754700169064081251403403753041651967143912262189933373118995238042877386626827110546536054889615983395995911424302907392},{13759723914831422079528203885258223827892294800653464336613431621643559367661176899825273167256203221590478048298507914521898161431012753081876596623677725187/14945994388086437769011430253313151823153117071911524247832770141525236720006719281798140852427279718095930038909221239606327478344888618053984529628910321664},{10905800490144425209837114446998007186658199330212945564109374684052801216977261873230653901123954978778999394127934351354055809290541755412012946319463790135953/13881917763101452581652850477878168344360015922703849798547477097452988020733079563314932382078174711424790948732467300957875020097756615365811600126954734354432},{251041172819713443725592358020181002712383523801153756101207323768030520378055617374887064331394104040212460576237023659934646953753614830323722519978381622101/313257290890330672067964502514173746959182607336614967615104099837992868982247230625296608236404015605045992598710155085824053517636544689261676144531257425920},{866699789166199204183056244785177937332849833260743560246691740338290064989520365289512960975226481456655113871868999865800937919903070372147836763756223545191/1100952466785995313290532703288719610559663904879889110883685107333984534392890743381928168466833784261929911409190357021090040753852316773815154850646796009472},{8616057628919110827886809146571179246660244410955179794438197389412399039199869653144725438119741835237608657724548890710341185544658135362260902091639061157/8397931496703540954350424561781171857308456031916638783749395818568961146605538042048774610378510831249486634440097031732935692670098248145223715658227777536},{8005989841192647662384724526472935921999912481385122091719517754815361471250908963668283328444648638039319732193673295864208861906567803280861179788548561385/2083135739487663144298266860141528609621630539960450965875125879320287382205014773680856190023780760663981215935623622242832914038224582795113755813330550784},{8780035208598627346081392029321974472025533288674489813558806632110012785404914656034246520904993376354439696890281687996743924375068214138068065166647434251/1447700522082308551679740165385465569472335600841463219487118096102337747881800470667756901734608665117723171168638720718069208095533213624939101906958024704},{389824438142029576641592817582370254214086967712534526529630750739400539116733844742667237490631625879431759306564667077311409500870529457054272432686530704849/461425465616439995460556702054613238437702794302690470544428321408819947679478476065652667334946816112922252906365072623730989529711533262572679320457397141504},{12945664885331024757274240500873822942173779068279000020581977107595459266325914949161902751202401215464525163990683804318645556393779551446888957336779636927/4406369224961884829940933946443467272524185355857394243394662792729811525521688920597341182452701172981815637992464308024691608848942570364757007155962314752},{41900213176306201984771235839907720265028923974300323691923410821849803978317404850818784950691320347133301733063981771144629019001483453147764259157960667883/9629519205880504361088633019008317182474443986950834905972874447315690231051461414350461246657099174176517959191455628972888020059590489335746192921409355776},{77795560562898884897511486887850372943918398335456357160146603136752382010810150565610045286577094724668908808799446787010897167051915006934104203012122058615/9944078543453215179724422933269932842528161633907844893664822715824811245386507694600611950840206545365268820619766476596353279319789273330326147909793349632},{6228758846638142515416586212109204276995541302733370125732564965622417160329831569503271533672274692980427649961781890512085327609380778679967795543859755656165/8131980661987623614252160604132676139420783804720754914251623431599707229895858588293004431427742211644292094236151816898096824671063697314766839313952898285568},{8477567247024780794889673500965501953324656552922409383523306677684442691779624463268644947418461190073245301954943499096865875696373167211237765692884121425511/11492710677339503479257836378573748448638779031456404761233180250238846610722412309461980726895825446099426800632426708196256581064707603187300104682016451067904},{14068777057825506150127060867284515703567136782273079140775492741100515457379832436265880554192609225790004128070133271145146598080087792840460492381291054749/4397473285811159232746081728298671997583143670615888316935504754704380134893560505715955782425363518904877378599239062040717517182873760105218051466105716736},{23950955760356012157042046077349615127684162572606654513179704204407705790705188244714639250903248299168966144677756493536993691038765793346066203648278614203/5986180962833356805675063993789555624438563959159224654972541164173962478633147623475429149763075671060656245181798490681108040053249623896902030052527439872},{431519916472467952227785792430096743012302309411806102929536458165434491396317230789154348817440906286754151072166977837284053441856818909140998684364999408779/180949457701975225365984493528002233958365886950414532476302852574730554742241118605110910165832768937398337462534472366966636075273433934890345531185618223104},{2166943917484864684989377022061172148145881170538315232849684028079658870072993332390857586712956973339702787281627419045916668916248002040010153317660086073/3096786052878271873957993526469911341107678090354001369027279042033297369365331988333279505769540237429682789323860520831095141292904445265437472377926582272},{1750511384650950242712842145014416038019114337759678384846632860769599848902030873095455867986942355299975740674629292257745424100448285341321875136965800756461/1942168380794794785897581108685388808750352513028194003929954387260072706509235159422240377738707207079591785995692939088586486904922535323144756956052517289984},{60370304917433533006379272621483171701696153313233508397977784737193975803191954552830982077828328658105434710519266535855271917566054662734591551559772818930281/53827866097303613010539985132386170195828897097181595259717832220116958985505635758922310582795551860255961157907149084224285613911218177518033142874875933753344},{103803783386077643626822380451448876952481785245019739559915001361351267148081991603752532503581884846427637877227605112343150551236526108024692819888722229847/77525194730702539116762051129547340768774700553638110926143710206558746274702365954636718008437629851468011658375865384598707582723456627918471934923453235200},{13895121544545023912987649337549344730870453414077834842793469630517229461520860000081279233765207182053094657994811112783395706568564153181462351911768239760121/6808645285369471153769083751580423592608631351177142153304963045004234925606509680571649510100066651327737204285908992534170647166257122453765794243024942792704},{339799488950137131000963384913403545326660537273945575638380240610850645260786195093689521379313093125956901207365537323003626630213988316814688041200096643851/154109530986063931244353409044661380406644960969767694473486541137183580323133020161320285922359808181616380116474070441739108004253382745694250112278622896128},{1086095757828201861007698208456322472772598433835005375121754699204248918988269239590274509059727882938775806550225857691740332992376201259435831774091366983771/599604523662333563169681030895659615621255492016713693264804587660756260616142329514964277309750584726391381714770027582109545872262210298367151310213379784704},{251973624347300799350669585589975181173950836019510097319915462505258533688666318070332498354372545777448105983040098249711117489551040829029515678560793337933/66645818209236109295585004513074510238489935536549864463000081542537590270376003776594571483018732119108983414808372511833908960965399703457378083013413306368},{134525883013323482151202830008499462178013034320483218623509851741860872293948920566921508179945561556583566604116714798702589472252318038186763760220150361443/33576590869635379066503080605416307024595913277131402847773304942868588407387488277831382425564440721195691779269257493093291267085283966592088567852419776512},{154719725952932949131638616612235461316318749628710782187000533706116414761307806831319309540263218988648541000193991531471256573169209596364280815716990432783/160349142225761939687505173414185146291539599676618754335595713370691712392887566440132086766595424926325428176012614843304203443291446834246038717897159737344},{85438293684376127406985339935800785825937911446045182846874930225506822454656415146869393365429410860916148924699426923215316206368924123636205148850277218461/92376431497254551162109308790891449487884542189171012592589604119456943810007103470431577135418893118108590390159259717128781744823102289458042789617622581248},{26989877917082349148786784423560094504502900492266797078196677508633058675885236187672727935181906552875856951577032538910357926044186200004830108445542026551/38019471143094979077124825065985803504496385466474878012279443022259683572016217393888476204054310907141454655473077728778588375578616402393508568595660013568},{88645784638532015929845099215775704943584007843292583627849052511512982150333629990711735918425812096832453206557984742355421295053635354786256600497885985713/108995048273098360794275171420523301353219402667825175970513930354139629726193876637403275770400282325924895983575972334079282330173586819839793352792648712192},{112793444939316507579550059956143690124321939949788002211817382401579258601644856369334719241100330669234020022150579412949093931793563829893293977108136465665/165976147700434378708758866745346975497995162373044518257192862013583244768294070541254878312215237872824008277905531354547150466969587512956633157614105001984},{335302892069012160325984347649201884333987049935540152777523278148966375135454618656439821280537330498657398432886338232374424293909391862489864937159920443361/172201164641992710519027650247116561502823643729702850067360540723593002248409916050596101846998822689054126987293146923424075921705928525303465638512352559104},{1486300419112159463830192173424672979421064105910920022348662337228718220127315197206139969132152086170552713396932915575604378695442786096537509164496415363557/740575831033931640018428893177170329411234086835841122784337853023772273831612034454227317625196498103727142020234751060071679570403502053330471213367713333248},{390351795874784135552043930337153393325884548117760241663811840060744112786833115307796084955335500890539860377536527295834845996488713523770215127083431572163/518151740931468650413523916507125877332940244014769912274150291817928999975308467982498454372349103921248096445655060890246130043719181780427103345942500212736},{61553484171713547618863038203167090382395863537434702577344623393835452259031306457306873503146235346647148939709257553208172297711758531720217094638026190719/69362331717530854940345987077346895423804070432615040267716634821267008124396902290217323169866534652216333360695686281607686693147509637335211384131597369344},{2539534234329152275911519029973377907204538783973730616358241620824852403150693411120900784423602651434592190914727239056052697276076906544681073260410355393/3333248342532543366025324977937043495705190453737599594680648706269083494933770397048152377073218802609727962851057206229078686128766875436937529546453811200},{480828263917302861723572947325928506512320527684139792848795334159963470552534150417568190610632298835601639033300689508972435616986347211530419855890575694591/606911566847026707191575196050588450477111603366126278319666299179382213287363675940993528360701895784551334629593127371848512791866236689183633824070007521280},{15191754756768617048052197201808783058916292318166819550098798026973251716735518148936595985460081451646722364601087127690068526519537674610261261173945882269/20071767039143915051936282167310162031154179996255837384854645844200493808646672889702324414418653667792177607086001190718142359431025710012025467032850923520},{64639003031596327427398785223879724103238672240463765285551561122117761862859812759479504681192893078095832143065043348411780946237725020537906896778288894066571/54168088888799283404693612473868634676036828416528697842687804527063251441197470834709337028107055344213739851993192680603913972034087644314099881068022768599040},{21890860532781191333028666231267134974097752624942136008700712498314379889753887614999544198878822118159374514053197429815501464629471257444409245241864324968289/11381605993884558690153138856280048691092857359989188582252180578541764131204292054017889537037625671501061888392930118664837776861921553738387492839686710231040},{8585569606760588846041447649598016034024076986616386520406357514416818922790591239565707046493846948842386431887333889692425472844869110282131813597445447013149/232050591862539556339495850576504016371974089400350626420926387582000600150698670371664578417046846078351873240560439073919004857192983973461957271111115735040},{111661981027672259264763453519626080725341679182065971207338371086085476948993585163457415503401412629642132817588975323181167818483022437591341643480926896139729/71878295856628494247914305139472345053068352497827470007981488167930713864766193926530578678841727289366484296876404160784854469962908512515666490062572942786560},{5123211527507272648125285579071495014929512386988514643792539874173972254289772236154339550920552174247437868832999244186107518905231875090558069249303311167961/5397154558039829433266044439814868436334249267042452488314888864310588386590474883234697936919039553360265234999480366728600513681675387738436011787979142987776},{11883076745208642231048064510107332908845950690696233086779205555725811224817012330787232018750697858015635067361765813898677367366142534937576454390144814168049/18862039835843719842512364898805521316572459054163933738939275349272905804144350274953718752529515066092328327299494551080397467460309032746587950686420848345088},{23601327277590887512228884877993518699526649500804713291360648442107422970520739745042674970332898153796967521158050612961501185939667672787671330778132992080515/29511664432309210214643273005566899299989702915534348151167961587256474499199551449229824241767976806143314374801684913825838757410851373038995634696446025924608},{11701411929843476065739624838546881436529674873512116576765610236297963697037440023489411108239760457295859420344085738017868475894410455625224964583053859347495/14456186147139378322643362585160746549945837968238906334438335085136398831783571518811513264735629825449594821092800341558891072289413154755284601121678624292864},{6517983494205696647599138959964336551581680243703891923709457215597218051286743762317452331502487384790815232661856979773739063350896433434971113849868833151719/1813928683505930552652520512111338523974396511841430154513989682668410967087006010924323178298941947625645394988043263199197612796205625312452466535575276486656},{3314894622312040863091802250210069338825070003466660760498042212138066433947695507251124679081511966782055931504710494716954053667029305114514416068939872113/1588376463940114370920004287445729756782790336246764493619111426653163770595712284793698519508382437140780579141944691105610648527902100874755077372672737280},{43081699153199264909807493240744207518948281613802731206106849138343271040915634406172908783416469403901277255020304357753088934121678540138070926193183285335/20876125468481859627835562771649462829110329602668593806952981605864134489084587261752007482604753955769681215516871705679672345496302489103540827280065429504},{32991549569759341869701593864831684540294645427588725561754089169074987129646620509817805605345393972926801533074324966890409588103263578177021612118840518689/16489916867514396222515504235171345240394581370189338978141189186598652191051156913937030201758046685766949804268139185171770782362664334297375715343202254848},{423681185872086338088367007142780983722351996671077893037138945636184802661685998514821973513685400245912901466943265696093978291690973109830361469468228372597/166135074935827474882710014190933829923574533980127656661842951786504204928964347615477539057028783472542170688708710183747093830981372448291891364742152847360},{118154489699339607981104781206682990441972886018875401776634316761542505605577526677653897332282221938454219895503311738796181120130592284940984739053762324421/34107757139505513258945494674538684180139610609883535652139323965169457751118023190955754720660588865367072014945882445730229013171902973392796979822457782272},{6847159234975234593779028149873290514026006569845935993143511778682732608578825567194022599288883948628420651377149747875637027427995746466017540264639575017/2476514236020760682144626059327896671640525892027682238279627699394676071441931934400683312857156323175136652450829843575237140334119348973434236142803222528},{2755514023501875313600183581796472876956247950970820598980680551370091051077540999948880804318766215636253745080723497477424455803755882303796554661275691901/1832877797160575830309874613285891547506320854630040232773470357520029741672443434628541799558559718535220206904825457806712146249524296063293874258212552704},{27394118924161250450366944142604411078213240063535655128570537029135052137855642685234646634848624758998607062860985905065927198067920910426093973919872029093/34863542012465519989661509789279894166487562180472854044348782220810541309446918112802332712404629455013812165550426128327037338154835237244978002068718485504},{1133044676003279414332511403087969399935340277431725334821568851522686699271712967525590507801580464143486796223764068544552456409046679490741865159388539862101/462759378900048634222553733285603551166876502117090987374228914791881638733283178823726321896966316764874301722555372338881180906813832265831012941430750445568},{8196472079720637143323062347285704531613052353185367303325681090857159818408189138851168494409225567978003488249940731212896990385570992225445352505842103857/425468115368131434060905913322781804784099508363885680857646936812931941328597883688022699502772290649768906475715306099387286093722753815431970323303497728},{46779283039962169430490585299487732278131119142468547919087525361051519393970941465821940818921839227550412510342399615142524053226182820336935130177868869567/66966355820147338730696526107854671607184471998635864041408957531552554715572215743263797969262392844950255545628142097488746036324251628703940988549805899776},{32516567829628820147954646110055982646056649874126546606995567935596467230909532068890398698501723787100204575547411782495961481306469766199059353826789563253/41746549991502230872522924606261139149620187377367139819642787562426641630839910147259088946437836534604119904055944161121872825010191073925857334000892248064},{156098177161756140214866173259463719003920694455174724033253267243410928494521306632303662609944430296054720761249462706981904610438093490369333174693689768215/38376912387481525575877373205299374533564672697224122399307222392414839350546326291299569089700042031510161547185906345307243828133725299645998378332417163264},{56035262301471445071009290353314428105986327625284056435593304499020629945219917352174994640480387754782366729032187208431473617056115564010961799323313214405/16397388995134521938965274977607946303040843942055059354204648694886950797987521280971740923483716682420363888536793543932811410685035664308290492702421680128},{1548069775704790330242190040996115330588767768114735967872600811210299353429141523190991956862436973096386472849297946640467446727075191526115338996642097827965/1931439738645564466779567661285174585444355193727066909692833968079358060213288175301304871477017259969580592186204582368154603679672886115154847469754149502976},{1149019931457418067208457882653910288641421471297987926959069412406841743839465769564463972198646553924496771005775886941184732555959236395878072533109810465457/1431724301899457857095254829222805608164602043394470327336938979480284061585776624262540695056066006567955794029831183817557351019677595680512163143690171187200},{276240762583047511926696420022047645366384188542599081037585965343235554759669987914311025343495047272926397288548229666702488385393630504082046921572613234875/345068407833004740600345180688798823190985109095377974598080536390404685078287567576522352705019658582614392003016470356103195582567954934601361206365067935744},{44871632923048699406557302158984111805088733183879338047769232144626471959493987615698908849735624093772057229121749307784489257163111009015944489474181332511/9556424064446680847486821595491331018390691330649854879042991042100962343898934991374898144447743537287580477419618106001281485362278863926483029875765542912},{341114516042130625314041417056574997713428987727742309263244155889918033869544506875756166914274982598308753344284631972999215857583047651818292463640441630717/131081959129788333982915074796840339499225659738825750440668741472820780958038659640267004423388077219089765416189889825731108780660329092312800444890015072256},{550010050114999023502962253278359775801154767001809283888595203888809242418084505384916449256689941917948865853715725286020480236514221760156796143037995863037/906734710254831672838818930274442646109313198314170161262685453329029793351053087989852669963758084841211647553258488499051319607799676790434327889383593607168},{44920837644682491470624357817502981723348515010012881228709791632333486224547639182582459840016918692874268905881714558036122752774853939659638068789990330161/40365929632800830413597562277960716954410868842950545098441212687993269838334737831665573411916856059151502087370161988124753969618201315012625799865564135424},{4571383765439406217598186814500890778476095934285573577942658209433510786700951708411802633419722151211927872961109737500982039821181882005033524083354278595823/3292684129558244076346871356582459824342037608901218770114929440194276150535830076553416719281105515001618831885496867362737526820701128262900943978410419421184},{6365359818700180510785620112980990449254697881347879564968935616010214728489060037972057726414208045482400647299530325881047898401873388011270842230218949503121/10467492539350751819740647109593586151288667479058917847831328416860488590564804032962214370560897206465072871595249164149697117077222858458077348422729284452352},{4771941987496963921120250719254095425670756071907372693816173180092753005010281758085517513916685863103720696226134669047655276390208400851091076248481555585695/7458567976766994120796460961883975037522923822191314132808545833126091760541229647549506504497562203769954588263753188072150088261693728218038856963281626595328},{97481522046193665012972172871215020980136053893809195363625774711750646388900194805163619660841998377121656126625553119231966858453411538395596068581626916623/137642425727255237224106691366113566248390868010864326918452339307547916946564957842579337156787082717206974320862613634100418519548988562670082162618894647296},{14256471871941712420587465549952464514703652034753009471393585159040891141579242633130679130497631365383699820351111144113627860544246770188370401114862961221/20715974773265493034619979203520797081629204089423311346389340579270598252860770169492584057707709669721646332223985858434712088964201698276487297418051715072},{1511443892341567388862805016561969955229577601217066250878397963160922815080021631677668212079082486189837231234773047769556963102915758723826931341003117160853/283984754084015724262378429281905322586927972381066748195097881204830500547641901854044311363802555903493522998496227385477261905604509028407523540097103298560},{14389955350893698898133138740040598327541188659060238535170927949757059942044983988566770647855131581127046869148358780800630580765138112734084213845335534010827/10121469464561986509037016191699148976227676119312706871594838942359401134008896465713310827867164492914768096187730159728652433260903155057567763828253325787136},{3099368602881777088259577359550980964402960615504417162350079364220361207209697967999064435720977090536643972544735575682409929295864654279139574540820161065393/2310504923749110663324883269461811376408286855956542815891798451761441690159350758422419236347538195348650683572861863940075790017224889827054404070846018617344},{7588722257832725311257573661872750741252847669584882295633510469986197653253967286491069862838348657059293279607759658163717610890091568555491950258713580715/5484458205595705653374097847146242514923628156757999906655490748198808351626574920930606266657844329281613028926020114842885657535270507357384318135121739776},{44096603198967032486038454141402109388957254054711906745840765173041908979424571110328947951116780627139130638106333251547542368011926402505716518889979181129/23575487526572302401616199854843339290034214889475278412140019608105266022798447239919016446979366771384483134706329377553613207196949615197427690726599163904},{35726726173385907176668766940100064442857610745576052506075605549616720936770570760224252297379814845214027785195477050179086122991955723226029290174541727543/31028550650485136079230753369337435696283541871780215392201965439881643792443142865025219314557695768615157324604777852277008989576608003601981592204547719168},{4222594318976150342688823082515538396052278008189477980678414975146674749524229109024603548842973279252628183714927754426815063361205323822514725391978142215891/3266531601556191343531394264761718181148708040993379670049666548660286700375753028613774695404306963069679645264683408443117770722915875478799776606365151133696},{9021671422169186850504447818467916321376399835855871674381000695859604566481952074110665554743747656435878275563027025551964783598273557009795971309268753571685/6689882728335042341119374362501125189798455295398917514008425490295255642751470654715022069248660600745375362354109246641741820548704739195450427318306808528896},{1324883044860883752138243071737721792111414459351455001387386061723291393155759216219705126944822265215718108529361341989608298186906219801078565911033990441663/268679987420617424096658455790204562206217079001247483330811991665054789608659385055976044456569342077531675668489664099211770092894108669932576620141556531200},{17508610716370208438850105953344519367858549017510739681688891514840650645810188181011025298341780986636634977570402761220524059724429196114158233159779281482011/8925393131230064760844925625020739890267631377162787608981203891656334705908816081525500240318984302303306680705510100706560737551969677604083730976189930536960},{11892505574169913296905094586133761756051236791802025221190868456053320795898709923784146874071239980811625051266181286884564913576744934298130393377936870989/8822208765604256907404456403668908491297139200094812415750723175520652066785658965672595598541197775610614293671602427953200213714480215109310562692417716224},{240796207265674565044050384539553433965609239337001352326815641806214900017008731207356309491670101836749654006570896180053357276892852862573328109397406024755/126503996284412521853259507816269345919191210687527012328597652358992208740911974559605845562258636718257976905219903283501309561609475554753755749690906771456},{628341535283092812068760077095586511805007545635024998706666706504290730932636936879911672657613257482788390688047054702173618957427595764897771789368979680901/245302998943698246804647082327919502230434523934113774028626226946085872319166350036019495161942664610139107099415024540455016311266260754736399947838747836416},{78791359552188469169056789619938813699328115060581905502233675866622706989744231122437922832530740416681249939516068794172807333983086881270454082461584817783/110092687661052727628296614042842408380114845628268332011222114351792823382718868230674600986527245803986420532596522096392820117180539927635367937419864178688},{177974491215781877599795522224108333096736577820315034722475029384090418744954461986399881792773969471277346944072028850029348481394376448152121094704595243139/251679240226766315310751972481650292824694126780251736925007258194658560750439524877354277822250892556326866677794188652974522325403741938868903792507457896448},{1461992568944360288916139253907890899336745664803245608699748666747061937663054591815675536145212372721187191449115636684149190120268154501313672605600546938797/2604570021871665269427844403459737439514054883779315678334511167769396962190413513063129840551032307612325560289324505995699678578120617904635994196896638828544},{4258879172880286136862746501930151803469727853211149640992491875982307760313359351251271511818815624721593770772062360754076635289414478606319878024884162481/5720428352562019155968917336869335115360520106567169949349512099629371140135251157271863781994091210830518105905319443998116443402991420948459984862379507712},{16609547568247466495837109569940404339523558448785371472635551187931419241390992272220915871183335893818679334606525065524972551551983016106922376597248533357/5221861082060482102419590445567907932109537714456552010810286155629728409758661266941361598350654622591916498416693744649777617475545054268702925254004047872},{1145324006686745452633088358569756979269879403779128260329920472668060784602272862506375335054384449079464483242647798686432410832086697468136091053855828899/121737929079351992913470003573852895411689543585469976961055335738950271770508441143168489203705394149769266913744144333850328120630537064089920849215225856},{43272699611061529524655951311035708546644942411034576625064386688852284745590471379394829162704766066209211994344510974173054294270981819242018927163576538993/4439682992279324737176345364730794811488728582046572709539779649931133374593435460665461220332463474329680424549637118529955009934955964069853446634641293312},{3814213830461673204937146032236268644711547895386526775622163648423141271803647478939531147307084343144049644845467182161370412517709223685337103457896151379/377708421818960719019289880475166588558806634916501070001560761649880872328755537487313270388054664363880776707036563279118358625624770740539984987394408448},{19138299073569802231458493239856033201624104303148729521761876746491107296822573446842581044287381840904579164918781133059171073004249932448198107692350320189729/10479691320459703588743034647487637885648757719344451967953034990029209820343639011412789490685075382935516583178531330800895758058551169975712311052004543692800},{172644503626012092299197343744927324808337462668840471774796734095067687368493868475381530432639655399104046378228196074043563687777779064245615321801774840130921/90128546884037845544595954551485384230422951207440153415246045287693490650785061904624490876142083010743386448093626236477188784390654954919078643888374367125504},{7143252158843809866640604787261781322355093120981295693609744153197535909230447883185222247086702971747015693607312564499816027435599446160993069738413519450935/5931278545229784727086651915442169483623844808952844133210382351851633709345604498025080793504467262163629400980475340436235065541226224886730380598077022011392},{4344555302463006430770386948697303601134548637157204075366428861611644647404746905232308695143618982275120806789930880972315502761026270547210576993481804452469/3917223088612128034300698023881727914429577214574140289085603940349814139389578116625747716757336206681932345908270380325936579710261350474428588720306056593408},{1014923832850109446422960095285997007460148261862650035679172536517895196768101892148386777740567068909288379497397638101376238119613997083004261602898575559929/234618388432882047174638578186160858363815036624844619753929742775772168676106232867148768760959544263230663730128293022337208822293088863698396360067927506944},{18822177263928872545956804862797004134899334027753376387610971126678532472390178422119965345176721878138635684291037919624630121945658206050061160526966343729/6385896753457149365902593153224984298415051677380727707847862376631896937597559010179278877275701265999097847379574046671286984473778449453324947293024026624},{98093146914342431259190649626313760656901586040626272850951138923389005327257532439023032287768258977733166058228638815418083851118199374352241755658352483365/115661061591597873720943999836079017565618097646419131621582700944485944189224167922422850578453809279962123466220738445973766874326470145184688305261512753152},{4461088270834934874413711329578983191761347839085632937593509682840609674896685876290005072423480225080624131209888592702130891486724897045797469795908597883/5434290734440919391708301862098053373329264047614991980196536342617265953355663053835711755917092329115128096303947729757927574305729226920318932239157035008},{20043462976235399128975457075978387609157508801602359643623983723355200170873274775472046247920955441407116554369956714679077450410041161948740025195280374523/33558327347373828177651227530309374103293179507127780475265310555716902868318980417348940052346023454702994896491622520519449566004599472367860038844469477376},{1454037758255863837392025075011703407169733565603913102646799056097217938779905802540355012934224559856149813027900928501786718875767351909458400994544082017/2088079275916462999023663587000062963234870775722956695113487164826180524442810953597519037104288870420153223096146105195509511468385190668953941223097237504},{17333124146018232228719608557369123800086388657399256391775923515451590304410874636161940118700826088645255373038422849873188939738953081275065938560251552371/3431084565530592451024702765682811176737947439612206030474407337539339792976768885884603321637240162622614392898682362841787900554125828606356368284364832768},{6224368017051054436400414641483337054696306190785871959522562508884554118614770630719748319370635270632624303244401355798146587577258296020971165497365704097/2056668339280273191677299475541531642925625058983457881475395722694903528810847839868736280595742786190015951609563750568981531474876813700438859900247867392},{2391117432163091330788523202570889178882694724022368147205241713895381828848749840717353801014358822930784679215842371308074642335653052546066347975034604818541/301373899238536877465912642669553752413224159618750940545699756598329667238823257874904760390306041592166464797495083319575843373606718971038788312993816903680},{1579243802563384174597303321766790806571750981579148047684713027433388432701972996881571303687292658892489544615724814001421154719261558622005701015074453363891/1198422018227742211864907289240895483533216344987935178229921099587613619277094015260399620792385084021943822087279772049148191355349984950742577523468114329600},{1593035451086947608533688086088614016638844814608925334681490648104946804737744204680378289238360806513658237640508284105862165767070479414618549388018878319497/885306370323749214478060209738546120194394073120877288827963025439847789324110214148198974060882516063032493336384709081968934274429573069953125360518509363200},{2379907237617596118228509937469413560396679243581817036396794523877139384455847169550664605485723006275478242854815858407869853834449689509440267225106831299841/8675545224103970372072802345883669071841428736179339515428242408357655913355386000215797589073007751671431161347011964157914017041492872208505380095179882496},{1635318638186395649148490363689285426397025836683264345444116147002257375062175688124233231793659085995627196811064739491563012756364544356940707156405075271191/616546910898992524539469035586585394751490194345388072406480131700102401224020983259141346639681247672241087774054086190761701383900868556051395726573733675008},{323171011982534758996561284437404643212693346715371415130098429258100410664955795023646423083573268801584873629229365456218914625048961728052788110653575457859/76438915063046507790985262632303965303863706076618562848175975245098551881281992207507593924397048172358630225019985647584373302893692338581764686487949934592},{692287176581167545285888640772338071323982430125161938636633579470870332221776482538314452952706095003802443945397558691105212394049900720081284170658779723/791601471422508599753442813648896567197977399732735073460819893163368633609771730533047273582151493438490131904031735385119927510475263939182286991823732736},{646181366916154281909492113817036468647504463963402923659315467651091315334093991717760616001284432979567017323843362869282956761975630127609532488306854896941/738116255859118599332046275996621919873266246549403370417106081557121446140972608615859974157108968392948229726899624377603498065908097988014123758101903441920},{3131002436771730480809566885054218176609094678482427180978871295863671639038438910814924005816084968074092961956531003747591147744684416384496205402066289119457/3545045510977252585464530812922022677692563059678653621290120138950071969328056639797098761699311535488025444245556179248022077588082507717166093355975075430400},{157892368728987038424743183753288704718980025812810419261670698411630330044538816884158038640390746041461123581915783438622176058410933738012937072966116854639/156044125522150620489922415548136075983815175827711168181638261625906621180976983329530613138449265160950286042078632638552632515917737414154752720407598465024},{205699390884176429109325705757389492970794422549200741146003768585998245667079651021056316807770369362786708385588241541300935580289006080666513680511658852793/239371192911782901227264511748782096775229414764000296013423506801041413699044283180187059036123416865117956184204961128003005164090061633717738590428114976768},{10434271882114816257080323266725927784996317530968770408994292471232981119991813091665696743101781346879363038031455730103951444435325322405117581050243839627/9600536503669155538662381982968910523009639290297416105459548445228319484964570362411224749784178865474440268797678646334822014938416010133977930742183755776},{1002867936516214968114094737790178202533242885903099380027399628652290293078000336341822329369652917346393468065295292738459221751557795360380342863933039311/1094282074496727959490119613483566685334034338591292852902367524745776369521593777610655394514484858636017362527288740776708089427476700428847770266524188672},{2620154430514118487809195562332231193452726020901664732943673605077095730577397574971437078328505245537045113090581755023724345262893898762478834593016494927691/1719985864110288933090208596820102387160826904977184468368380343148848251126378480045669378160088472770959686373012653215132973288703735989264051036589467369472},{112959023230694513403829976668093986287351620681246676320287703002478063675337760483341611849360911097230244582747680271456919332463541827167690775757458088189/90709877702940377418320529921713076661999212316738632951803619308129187338657232522207095647148757494317299430488616033400928064767994063103737185779950026752},{22386154094821938171598169484510090323719015757268812633230730444911768044836550628633238302226440619902058920200248190884975166825517780043029699079560810081/25545145444789047804716297393541720734936281652443034952726503976960265443566807959182357064079060448844280191200508445196035128892880318234365885234448171008},{152728661986989402710674770979112664571531789718941708586782988501568447977810932580529597097107614780698059438434602709090150823606401274394137753844924930097/161548968285998687518024223733533331221512253074488171278132057428630049980376251764170651992046265968845555967800850647878281878983780396447093250214472974336},{27918740987641820360998766816439595548220034927699879832578731253793882028021156019393392670047164648346477616294594430849910349164920857508906779907967467815/27483905867064824921695899744015407797139297063725753967798790878905556719072322520344693313627548859467628726499603884766669677371490843497934310058517069824},{24146546505791314029053210847370829074861524272611779332849258815948568188693679388216986776539486320583287110296952822007188695422182290912359789023882479665141/4206396484847373785041025510341911049770667981745612347639703914107697933779537912184358959215327379951331047343426806052275844481669106685590425387469499269120},{107199283880203165683114402876819526901035365316547168043334617541074927976119007339293879371758593633911517801235257359956217964628342456915294749636795931431/17340336704893755597331960496188185132722925993758918178173805746881179723784055121850811162767974670225875683341793022005033097539008707476980984794767687680},{89721250939076974804982601865249399708133650241791552129393546124460594034315705877478462233887820213535513945477553103787145951907102640941758259982183350415667/63697344768165394365247826630055423700815950857849383586265635724306945189367562251259059253310536078984440221179859736816851395169642387407427018918782941790208},{25390694025650948951289639049588488139191711117333012787608197930112424555431546229200530226845499000236895078334153336404875201720978053926716879996124543164431/13068375147377908232572598370077085960354343055999018646560184759346717940865921586018534871396449729668672628998355608369652799162511252617824770830938555809792},{7447500740331807441798733416053823285531469839960650911682810836391848776234075722480094647714025231193903629491465989559002205951064638621860940334635442638191/1177640543594425952224020294350981851197911353361285007125858301223484588070204725157057087865514526913677878994979251118078054594793821165564550575389170204672},{2315138071362888286187553296234301286461387844537565791813436665977178034816725420257294480084163632454708476315420686857039209917982949999478219665839494428443/1459931344568102423221308479849659901657053561772564725132736566252500291015827030611651050268056851015799254735867882431360452859642592001779438283600858972160},{185419895537620748928199932709060606717394192690046880559014335208430269597733897315623546965549164444797936626304495309091570920769559758725288317158217256005/82386667729762843248864470730594269309667403695699735080441500671142530529622707080918598083186492601921014006235373048698834146333094688442889803713271037952},{2247239168608336671413500209208779115968690113673198781542969020131350146741006605203267251649938728571103143426663731447544257218625288435360123193758379448515/1002563016440900285244563859401917410455170707030501325386050373073583458057375229604176617145615058742671382197922562190914807614901507463243315704030325374976},{58343632713886189904094356683271686045806277106385898061023584740590071320871933959346334990141666829302390899791887086102158650284796576514411404864227153319/37752492641342199604712415661115909741798159052954668147685023957220649231444999826629988773664442513932246481700812914350475745396029870805066729754514685952},{26123957907446060246408912628795398731759013680661330714055278712430321938070753721346697753997108320227280842762643970531794098936054922177659633614469637135/18052549342592463211945002433689055414274089890026917992850958532115173796097662240211808490943739568907466460565565783847633221530229346018547698502033473536},{680069657237976014516160625042940173019801040686631352313917623667476680765483507814612623835074222460744092568344383474162196118342756760516723236465598545751/146159905979074015108752696572053938907614071046230666435034817534146691553831892332608348434380522091524725787636169700068433758498205230566656439765652471808},{36302604703003981771373934845311823814836235329567826874960139862018569350885832185269283693578318841075481945751223168138565444785004529590439876068410236692839/4490068032996895345533081138143260387498627000705644511726328798187119030342253894090247547540662153983954531521330523839110738391997893920923773435037154803712},{22195775340019645524800323236933959877808433459466353893025998672533119106729626003334763484296160685402583264396704718972843289535538476701047756842669452111283/13101871198199916981663075079232276925972158804997088444193568756202605716652289176141967518819984394554125332126517355733819145032267141764395872494631688077312},{41673655919115116072043983460876410455399281812129411892193704952043241136889290715988298302857807803767656164476506472608215568428709308159949074187604971444749/23976351473726231394655900286178365665621102906039320091563469058507173424229921529088610865510084503372617493375630388134159678179193886603432524275960239882240},{4315158934285601911440991315167102848250371608571329329937787387494316700343498907176900296002313135172891020921484309843452481104412172794132096970766479299279/2898062522205863304375655161064180018807503069825071670227803059402013859819528367520924221865114774255288904735936418690541276934872832334377315436996962811904},{26928601268460386644913880647064677331200528877396196447079409969653995570767427447029765467639607623257875060841821453159552780824013968272995872717918658453847/18551248734177534920839130478149259345199815387288132787734392296464117908007417995729262187881045512385846228101079015449688532859099261922867235433806408187904},{1881741638941037079984643769914905399539603205834380550037584610289765129331564902991614291551836729157615084969487779961424596878207628281649334345696684338507/1399878087514381708850090858788814065299446409730461371671656969975876472173955129719964678231517147059749588079867623808814963163661219625780387081251580805120},{40809676569028246985361157876370768964339420294538526161555887093118377930707387636446451070538886323403009105229731588713083714520085580850755839314158382933789/25090696144946597319415352454567522749469132405305760923334935249270844737454570917305540396095497961685137390333686652511863061974583604381220746647225591922688},{522960955522963107731506039659179942205206989407668245570075477992884642055163212307530815702636118008766656718857645104443868203035943663455898950161375826063/360536040127288722671475784521902203890662083146047019505486702944173208959291375256427406258511156779675313176947441827704163926470746366264762090598162235392},{236289384711023866051044230424297653884198222491112295952301718256954724567534356148753698068585725159466273126132506344902656998525684991053514908127006912461031/147295085135304663559869834142214030952723053071371936441910358997600310984215481503338097424662217722330271114321852786764193237855666891394828420603008939720704},{213580196473996088700308451876082904760614059042829050100159352549928289755591637595844321582196801338416529718895615300765400507044753530379805446029553302513/204417100882467485108161438316799895452195340295560412617807990388388099794357651695139108428805686224748124538440268228508620059335421439972963433996651331584},{141700671710383448188739435211844549726621073492423550769822087205829907847816339224662831393120828129191601879564217944614567448266220251083756833416536769133/132035308732189492901595164126654378474265966414054544471872983170917510698709052813443512171845997513732502383307138257922149720181806493708393946324951629824},{255905979595478458343756215931541667817390598957217174533627539069844066305751761122552851839783162973554359690577383064325072029104341364496565065421886313829/180949407429398155173808226218187735000955622058213711256670990976112028372189298034544089259642192740447237002623050242431959038535808324722279831669164736512},{18108939448949227967868966017557880875786308379327218749723182176744522110034408631916112716343369582579808579838912058031819205907111096372347690959598470651/13220654006022562539180508153058639045594282157117264257608862331294766971490041888538564299212729797185179306061254675333630228562805040924450422804367540224},{3510053867810965274403959264316019952251024843699655973829225360892834386869200947123369025398681537545252531704949713547159779769214450943273080653960511341/2012013418362418228122695925505287004848153785743346223021327903182224385429750582766258779156268992550222379488222238620617822530664016750454058549245902848},{136246160245518966033415099438040020892891275446989025732159497536520385538723326499314101257993752863752920560123566813398050579731930472021918263044432818787/66639532558652577089149747560000746477084411267612558913392995441085111602205280664448439936963999737501206542212503729268278643401553156680297657609314369536},{1200933372048413949376282448279358429340978645364748532022731600055097582390784788551939250521092313117613135628571899326469675742178592736882082081440507285049/763786653614789544243526129466345594868827570738633357665527338075608010154084376985491944786457400481027128899319780118079487672431540022776100584162353741824},{116301833286186718834831451315960491981299462209147236228593034971758265101414495390689357542589821106468551139035150015357148678126953996498209362785544790893/208078536178638909702782265819577009230953621778592151532225984395282269462677723544942298351054037305379180873393012341766548901790404880611001619412894613504},{966573560307850606413293798704918512070763943811923825315549127571689815133439194296743844483318144075253976679946164728031517335570044000626066223398905152115/211592355453841117910582823995189238872209587756828713921358027146939453952640665305463521081966399312036051341064656649377621485766350078698836446470949830656},{92680778963129380462790073728390908248689060229617729305679624639660259499765486942388823560442664714744633358072017382274197688117281597301632967234388362323/129066531575970299875172432717853287066494628452719648082362785535090171103431543830809769088226294418693698670401259424218802286699805978890069914986555113472},{2439093242478012961934725127626976068277746026265095322765543440839123440295109550692685192264351938341643713976679662050568987225405678183850871878910549069/2687506160888582052129854465856555703987174278851475040521512079077985467216502255236444773822296341494886033635809216954546169911210747159588183498351968256},{218276906427675380772765525768650785958897793474024474693427411521282366168298413292055552919424635764249166222168501986369877615888807446090402549022134174415/308481188859234670027975639636396103358916446272460570763841070519779750742189406331668472510268442809620377723763427096437026747807352118312919747240084373504},{41907767971480445195788833393385945672702178895339019234480110153899187851709710848810300434592044512893776525990953578491759797743880816726107953247162