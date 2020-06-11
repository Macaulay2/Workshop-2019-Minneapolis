
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
        (cellComplex,Ring,List)
    Headline
        create an cell complex
    Usage
        cellComplex(R,maxCells)
    Inputs
        R : Ring
            that specifies the base ring for interpreting the labels
        naxCells : List
            that specifies the maximal cells for the complex
    Outputs
        : CellComplex
            a cell complex
    Description
        Text
            This function constructs a cell complex from a list of cells
        Example
            c = newCell {}
            C = cellComplex(QQ,{c});
            R = QQ[x,y];
            C = cellComplex(R,{c});
    SeeAlso
        (newCell,List)
        (newCell,List,Thing)
///

doc ///
    Key
        newCell
        (newCell,List,Thing)
        (newCell,List)
    Headline
        creates a new cell
    Usage
        newCell(boundary,label)
        newCell(boundary)
    Inputs
        boundary : List
            that gives the boundary of the new cell either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was created
    Description
        Text
            This function creates a new cell, to be added to a cell complex,
            if given a list of cells without any orientation information, it
            attempts to infer the orientation.
        Text
            If given an empty list for the boundary, the function creates a
            0-cell (a vertex).
        Text
            If not given a label, and the labels on the boundary are monomials
            or monomial ideals, from the same ring, then the label is the lcm
            of the labels of the boundary
        Example
            R = QQ[x,y]
            a = newCell({},x);
            b = newCell({},y);
            c1 = newCell {(a,1),(a,-1)};
            c2 = newCell {a,a};
            c3 = newCell {a,b};
            C = cellComplex(R,{c1,c2,c3});
    Caveat
        This function does not check that there is a valid map from the boundary
        of a n-cell to the given boundary. It only checks that the boundary
        forms a cycle in homology.
    SeeAlso
        cellComplex
        newSimplexCell
///

doc ///
    Key
        newSimplexCell
        (newSimplexCell,List,Thing)
        (newSimplexCell,List)
    Headline
        create a new cell
    Usage
        newSimplexCell(boundary,label)
        newSimplexCell(boundary)
    Inputs
        boundary : List
            that gives the boundary of the new cell either as a list of pairs
            of cells and their orientation, or a list of cells.
        label : Thing
            that gives a label to associate to the cell, otherwise attempt to
            infer it based on the labels on the boundary
    Outputs
        : Cell
            that was created
    Description
        Text
            This function will only create simplices, and it will verify that
            the new cell is a simplex, as such does not have the caveat of
            @TO newCell@. Otherwise it has the same behavior. This is
            particularly useful in constructing \Delta-complexes.
    SeeAlso
        cellComplex
        newCell
///

doc ///
    Key 
    	(ring,CellComplex)
    Headline 
    	return the base ring of a cell complex
    Usage 
    	ring C
    Inputs 
    	C : CellComplex 
    Outputs
    	: Ring 
	    the base ring of the cell complex C
    Description 
    	Text
	    This returns the base ring associated to a cell complex C, which is used to interpret labels. 
    SeeAlso 
    	cellComplex
///

doc /// 
    Key 
    	cellLabel 
	(cellLabel,Cell)
    Headline 
    	return the label of a cell 
    Usage 
    	cellLabel C
    Inputs
    	C : Cell 
    Outputs
    	: Thing 
	    the label of the cell
    SeeAlso
    	Cell
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
        : ZZ
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
        : ZZ
            the dimension of the cell
    SeeAlso
        (dim,CellComplex)
///

doc ///
    Key
    	cells 
    	(cells,CellComplex)
    Headline 
    	return the cells of a cell complex as a hashtable whose keys are cell dimensions
    Usage
    	cells(C)
    Inputs 
    	C : CellComplex  
	    the CellComplex whose cells are to be returned 
    Outputs 
    	: HashTable 
	    the cells of C 
    SeeAlso 
    	(cells,ZZ,CellComplex)
    	
///

doc /// 
    Key
	(cells,ZZ,CellComplex)
    Headline 
    	return the cells of a cell complex 
    Usage 
    	cells(ZZ,C)
    Inputs 
    	r : ZZ
	    the dimension of the cells to be returned 
    	C : CellComplex 
	    the CellComplex whose r-cells are to be returned
    Outputs
    	: List 
	    the r-cells of C
///

doc ///
    Key 
    	(boundary,Cell)
    Headline 
    	returns the boundary cells along with relative orientations
    Usage 
    	boundary(C)
    Inputs 
    	C : Cell 
    Outputs 
    	: List 
	    two-element sequences where the first element is the boundary cell and the second element is the orientation of the boundary cell relative to C    
    Description 
    	Text 
	    Given a cell C, this command returns a list whose elements are two-element sequences. The first element of each tuple is a boundary cell of C and the second element is the orientation of that boundary cell relative to C. 
    SeeAlso
    	(boundaryCells,Cell)
///

doc /// 
    Key 
    	boundaryCells
	(boundaryCells,Cell)
    Headline 
    	returns the boundary cells of the given cell
    Usage 
    	boundaryCells(C)
    Inputs
    	C : Cell
	    the cell whose boundary is to be returned 
    Outputs 
    	: List
	    the cells in the boundary of C 
    SeeAlso 
    	(boundary,ZZ,CellComplex)
///

doc /// 
    Key  
	(boundary,ZZ,CellComplex) 
    Headline 
    	compute the boundary map of a cell complex from r-faces to (r-1)-faces 
    Usage 
    	boundary(ZZ,CellComplex) 
    Inputs 
    	r : ZZ
	    the source cell-dimension 
	C : CellComplex 
	    the CellComplex 
    Outputs 
    	: Matrix 
	    the boundary map from r-faces to (r-1)-faces of C
    SeeAlso 
        (chainComplex,CellComplex)
        (boundary,SimplicialComplex)
        (chainComplex,SimplicialComplex) 
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
            the cell complex for which to compute the chain complex
    Outputs
        : ChainComplex
            the dimension of the complex
    Description
        Text
            This constructs the cellular chain complex for a cell complex,
            taking into account the labels on the cells
        Example
            R = QQ[x]
            a = newSimplexCell({},x);
            b1 = newCell {a,a};
            b2 = newCell {a,a};
            C = cellComplex(R,{b1,b2});
            chainComplex C
    SeeAlso
        (boundary,ZZ,CellComplex)
	(boundary,SimplicialComplex)
        (chainComplex,SimplicialComplex)
///

doc ///
    Key
        (homology,CellComplex)
    Headline
        compute the homology modules of a cell complex
    Usage
        homology(C)
    Inputs
        C : CellComplex
            a cell complex with labels in ring(C)
    Outputs
        : GradedModule
            the graded module, the homology of C with coefficients in ring(C)
    SeeAlso
        (homology,ZZ,CellComplex)
///

doc ///
    Key
        (homology,ZZ,CellComplex)
    Headline
        compute the homology modules of a cell complex
    Usage
        homology(r,C)
    Inputs
        r : ZZ
	    an integer
	C : CellComplex
            a cell complex with labels in ring(C)
    Outputs
        : Module
            the r-th homology module of C with coefficients in ring(C)
    SeeAlso
        (homology,ZZ,CellComplex)
///

doc ///
    Key 
    	isMinimal
    	(isMinimal,CellComplex)
    Headline 
    	check if a labeled cell complex supports a minimal resolution 
    Usage 
    	isMinimal C
    Inputs
    	C : CellComplex 
	    a cell complex 
    Outputs
    	: Boolean
	    true if the cellular resolution supported on C is minimal
	    false otherwise
    Description
    	Text 
	    This determines whether the cell complex C supports a minimal free resolution of the monomialIdeal I generated by the labels of the vertices of C. 
	    Note: we assume the cell complex is free. 
	Example 
	    R = QQ[x,y,z];
	    v1 = newCell({},x^2*y);
	    v2 = newCell({},y*z);
	    v3 = newCell({},z^3);
	    e12 = newCell({v1,v2});
	    e13 = newCell({v1,v3});
	    e23 = newCell({v2,v3});
	    f123 = newCell({e12,e13,e23});
	    C = cellComplex(R,{e12,e23});
	    isMinimal C
	    D = cellComplex(R,{f123});
	    isMinimal D
///

doc /// 
    Key
        isSimplex
        (isSimplex,Cell)
    Headline 
    	check if a cell is a simplex 
    Usage 
    	isSimplex C 
    Inputs 
    	C : Cell 
	    a cell
    Outputs 
    	: Boolean 
	    true if the cell C is a simplex 
	    false otherwise 
    Description 
    	Text
	    This determines whether a cell C is a simplex. 
	Example
	    v1 = newCell {};
	    v2 = newCell {};
	    e1 = newCell {v1,v2};
	    isSimplex e1
	    e2 = newCell {v1,v1};
	    isSimplex e2
///

doc /// 
    Key 
    	(facePoset,CellComplex)
    Headline
    	generates the face poset of a cell complex 
    Usage 
        facePoset C
    Inputs 
    	C : CellComplex 
	        a cell complex
    Outputs 
    	: Poset
	        the face poset of C 
    Description 
    	Text
	    The face poset of a cell complex is the poset of cells with partial ordering given by inclusion
	Example
	    R = QQ;
	    v1 = newCell {};
	    v2 = newCell {};
	    v3 = newCell {};
	    v4 = newCell {};
	    e12 = newCell({v1,v2});
	    e23 = newCell({v2,v3});
	    e34 = newCell({v3,v4});
	    e41 = newCell({v4,v1});
	    f = newCell({e12,e23,e34,e41});
	    C = cellComplex(R,{f});
	    facePoset C
///

doc /// 
    Key 
    	(skeleton,ZZ,CellComplex)
    Headline 
    	returns the $r$-skeleton of a cell complex 
    Usage 
    	skeleton(ZZ,CellComplex)
    Inputs 
    	r : ZZ 
    	C : CellComplex
	    a cell complex 
    Outputs 
    	: CellComplex 
	    the $r$-skeleton of the cell complex $C$
    Description 
    	Text 
	    The $r$-skeleton of a cell complex is the union of its cells whose dimension is at most $r$. 
	Example 
	    R = QQ[x];
	    P = hypercube 3;
	    C = cellComplex(R,P);
	    dim C
	    cells C
	    S = skeleton(2,C);
	    dim S
	    cells S
///

doc /// 
    Key 
    	(cellComplex,Ring,Polyhedron)
    Headline 
    	creates cell complex from given polyhedron 
    Usage 
    	cellComplex(Ring,Polyhedron) 
    Inputs 
    	R : Ring 
	    that specifies the base ring 
	P : Polyhedron 
    Outputs 
    	: CellComplex 
	    whose cells are the faces of the given polyhedron P
    Description 
    	Text 
	    Given a polyhedron, this command returns the cell complex whose cells correspond to the faces of the polyhedron. 
	Example 
	    R = QQ;
	    P = convexHull matrix {{1,1,-1,-1},{1,-1,1,-1}};
	    faces P
	    C = cellComplex(R,P);
	    cells C
///

doc /// 
    Key 
        relabelCellComplex 
        (relabelCellComplex,CellComplex,HashTable)
    Headline 
    	relabels a cell complex 
    Usage 
    	relabelCellComplex(CellComplex,HashTable) 
    Inputs 
    	C : CellComplex 
	H : HashTable 
	    whose keys are cells of C and whose corresponding value is the cell's new label
    Outputs 
    	: CellComplex 
	    whose cells are relabeled by the values in the hashtable H 
    Description 
    	Text 
	    Given a cell complex C and a hashtable, whose key-value pairs are a cell from C and a new label for that cell, this command relabels C accordingly. Labels for cells not provided in the hashtable are inferred, unless the option "inferLabel" is turned off. 
	Example 
	    R = QQ[a,b,c];
	    P1 = convexHull matrix {{0,1,0},{0,0,1}};
	    P2 = convexHull matrix {{1,0,1},{0,1,1}};
	    P = polyhedralComplex {P1,P2};
	    C = cellComplex(R,P);
	    verts = cells(0,C);
	    v1 = verts#0;
	    v2 = verts#1;
	    v3 = verts#2;
	    v4 = verts#3;
	    T = new HashTable from {v1 => a^2*b, v2 => b*c^2, v3 => b^2, v4 => a*c};
	    relabeledC = relabelCellComplex(C,T);
	    for c in cells(0,relabeledC) list cellLabel(c)
	    for c in cells(1,relabeledC) list cellLabel(c)
///