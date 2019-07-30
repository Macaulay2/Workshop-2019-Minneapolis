restart
load "nops.m2"

R = QQ[x,y]
I = ideal(x^2-y, y^2)


socleMonomials = I -> (
	R:= ring I;
	inI := monomialIdeal leadTerm gens gb I;
	irreducibleDecomposition inI / gens / 
						entries / first / 
						product / exponents / first / 
						(i -> (i / (j -> j-1))) / (e -> R_e)
)

SM = socleMonomials I
G = gens gb I

rewriteBackwards = (g, f) -> (
	m := (sort terms f)#0;
	if g%m != 0 then error("Least term of f does not divide g");
	g//m * (m - f)
)

writeBack = m -> (
	f := m;
	while true do (
		smallest := (sort terms f)#0;
		))