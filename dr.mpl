# Dixon Resultant Package
# Author: Manfred Minimair
# E-mail: manfred@minimair.org 
# Released under GPL Version 3  
# Copyright (C) Manfred Minimair, 2015";

DR := module()

  export DixonResultant, DixonPolynomial, DixonMatrix, RankSubMatrix, DixonExtract, 
	 #Sub-modules
	 Info, Extract, EDF, DRes, GenPoly, Samples;

  local MonomialCoeff, PolynomialToMatrix;

  description 
  "Dixon Resultant Computation.", 
  "Date: 12/18/15, Version 2.0",
  "Author: Manfred Minimair",
  "Additional contributors: Arthur Cherba, Hoon Hong",
  "E-mail: manfred@minimair.org",
  "Released under GPL Version 3",
  "Copyright (C) Manfred Minimair, 2015";


  DixonResultant := proc(pols::list(polynom), vars::list(symbol))
    description
    "Compute the Dixon resultant.",
    "pols::list(polynom) ... list of (d+1) polynomials",
    "vars::list(symbol) ... list d variables of pols"
    "Return the Dixon resultant of pols.";
    return DRes:-RankMaxMinorPivotRow(pols, vars);
  end proc;


  DixonExtract := proc(M::Matrix)
    description
    "Extract the Dixon resultant from a Dixon matrix."
    "M ... Dixon matrix",
    "Return message, expression. If message is empty then the expression is the Dixon resultant of pols Otherwise message is an error message and the expression is only a factor of the Dixon resultant.";
    return Extract:-MaxMinorElim(M);
  end proc;



  RankSubMatrix := proc(M::Matrix, {params::list:=[]})
    description
    "Extract a maximal-rank minor of a parametric matrix.",
    "M ... parametric matrix", 
    "params ... parameters of M",
    "Return a maximal-rank minor of the parametric matrix M with high likelihood. If params is empty, the parameters will be automatically determined.";
    return Extract:-RankSubMatrix(M, parameters=params);
  end proc;

  DixonPolynomial := proc(pols::list(polynom), vars::list(symbol), auxVars::list(symbol):=[])
    local nVars, nPols, new_vars, CM, i, j, dv, dixp;
    description
    "Compute the Dixon polynomial.",
    "pols ... list of (d+1) polynomials (different order just changes sign)",
    "vars ... list d variables (different order produce completely different polynomials!)",
    "auxVars ... list of d new variables. If it is empty, the new variables will be automatically generated.",
    "Returns dixPol,newVarList, vars,",
    "where dixPol is the expanded Dixon polynomial, vars is copied from the input parameter vars and newVarList is the list of d new variables used by the Dixon polynomial.";
    #Implementation:
    #Expands cancelation matrix using maple LinearAlgebra, det function.
    #Examples:
    #  dixPol, newVarList, vars := DixonPolynomial(myPolSystem, [x, y]);
    #  (inewVarList = [x_, y_], vars = [x, y])
    #  dixPol, newVarList, vars := DixonPolynomial(myPolSystem, [x, y],[alpha, beta]);
    #  (newVarList = [alpha, beta], vars = [x, y])   
    # Author: A. Cherba, 2003
    # Modified: M. Minimair, 5/18/2015

    nPols := nops(pols);
    nVars := nops(vars);

    if not(nPols = nVars + 1) then
       ERROR(`Number of equations should be number of variables plus one, but number of equations is %1 and number of vairables is %2`, nPols, nVars);
    end if;
   
    if auxVars=[] then
	new_vars := [seq(cat(vars[i],`_`),i=1..nVars)];
    else
      new_vars := auxVars;
    end if;
    #print(new_vars);

    # Construct Cancelation Matrix
    CM := Matrix(nPols,nPols);

    for i from 1 to nPols do
      for j from 1 to nPols do
	if (i = 1) then
	  CM[i,j] := pols[j];
	else
	  CM[i,j] := subs(vars[i-1]=new_vars[i-1], CM[i-1,j]);
	end if;
      end do;
    end do;

    # Now, divide the determinant by the product of (vars[i-1] - new_vars[i-1]).
    # This is accomplished by subtracting ith row from the (i+1)^th, and
    # dividing the result by (vars[i-1] - new_vars[i-1]).

    for i from nPols by -1 to 2 do
      dv := (vars[i-1] - new_vars[i-1]);
      for j from 1 to nPols do
	CM[i, j] := simplify((CM[i, j] - CM[i-1, j])/dv);
      end do;
    end do;

    #print("CM", CM);
    #dixp := expand(LinearAlgebra:-Determinant(CM), method=fracfree):  # Dixon Polynomial
    dixp := expand(LinearAlgebra:-Determinant(CM)):  # Dixon Polynomial, let Maple choose the most appropriate method

    return dixp, new_vars, vars;
  end proc:


  MonomialCoeff := proc(f, vs, m)
    local v, i, n, zero, x, k;
    description
    "Coefficient of a monomial in a polynomial.",
    "f ... polynomial", 
    "vs ... list of variables of f or a single variable",
    "m ... monomial in variables vs",
    "Return the coefficient of m in f.";
    if type(vs, list) then
      v := vs;
    else
      v := [vs];
    end if;
    k := [];
    for x in v do
      k := [op(k), degree(m, x)];
    end do;
    n := nops(v);
    zero := [seq(0, i=1..n)];
    return coeftayl(f, v=zero, k);
  end proc:


  PolynomialToMatrix := proc(pol::polynom, rowVars::list(name), colVars::list(name), roworder, colorder)
    local rowPols, nRows, nCols, rowMon, colMon, rowMonS, colMonS, 
      Arow, Acol, Trow, Tcol,
      M, i, j, allVars, p, mon;
    description
    "Convert a polynomial to a coefficient matrix.",
    "pol ... polynomial to be converted",
    "rowVars ... polynomial variables for the row indices of the matrix",
    "colVars ... polynomial variables for the collumn indices of the matrix",
    "roworder ... polynomial ordering for the row monomials, compatible with Groebner:-MonomialOrder",
    "colorder ... polynomial ordering for the column monomials, compatible with Groebner:-MonomialOrder",
    "Returns matrix, rowMons, colMons,",
    "where matrix is the coefficient matrix, rowMons and colMons are the monomials indexing the rows and respectively columns of the matrix.";
    p := pol;
    rowPols := [coeffs(collect(p, rowVars, distributed), rowVars, 'rowMon')];
    nRows := nops(rowPols);
    rowMon := [rowMon];

    coeffs(collect(p, colVars, distributed), colVars, 'colMon');
    colMon := [colMon];
    nCols := nops(colMon);

    Arow := Ore_algebra:-poly_algebra(op(rowVars));
    Trow:=Groebner:-MonomialOrder(Arow, roworder(op(rowVars)));
    Acol := Ore_algebra:-poly_algebra(op(colVars));
    Tcol:=Groebner:-MonomialOrder(Acol, colorder(op(colVars))):

    colMonS := sort(colMon, (t1,t2)->Groebner:-TestOrder(t1,t2, Tcol));
    rowMonS := sort(rowMon, (t1,t2)->Groebner:-TestOrder(t1,t2, Trow));

    allVars := [op(rowVars), op(colVars)];
    p := collect(p, allVars, distributed);
    M := Matrix(nRows, nCols);
    for i from 1 to nRows do
      for j from 1 to nCols do
	mon := rowMonS[i] * colMonS[j];
	M[i, j] := MonomialCoeff(p, allVars, mon);
      end do;
    end do;

    return M, rowMonS, colMonS;
  end:


  DixonMatrix := proc(dixonPoly, rowVars, colVars, roworder:=tdeg, colorder:=tdeg)
    local dixMat, rowMons, colMons, outRowMons, outColMons, outRowVars, outColVars;
    description
    "Compute the Dixon matrix of a Dixon polynomial.",
    "dixonPoly ... Dixon polynomial",
    "rowVars ... variables for monomials to index the Dixon matrix rows",
    "colVars ... variables for monomials to index the Dixon matrix columns",
    "roworder ... ordering of the row monomials indexing the Dixon matrix, compatible with Groebner:-MonomialOrder",
    "colorder ... ordering of the column monomials indexing the Dixon matrix, compatible with Groebner:-MonomialOrder",
    "Return matrix, rowMons, rowVars, colMons, colVars,",
    "where matrix is the Dixon matrix, rowMons and colMons are the monomials indexing the rows and respectively columns of the matrix";

    dixMat, rowMons, colMons := PolynomialToMatrix(dixonPoly, rowVars, colVars, roworder, colorder);
    if type(infolevel[DR], numeric) and infolevel[DR]>=2 then
      printf("  %a: Dixon matrix of size %d x %d.\n", procname, LinearAlgebra:-RowDimension(dixMat), LinearAlgebra:-ColumnDimension(dixMat));
    end if;

    return dixMat, rowMons, rowVars, colMons, colVars;
  end proc:


    ###########################################################################################################################
    Info := module()
      export System, OpenFile, SaveStart, SaveEnd, RunMeasure, DixonMatrix, Extraction, Measure, Rank;
      description 
      "Various functions to gather information about the computation.",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";


      OpenFile := proc(fname::string, {header::truefalse:=true, discard::truefalse:=true})::filedesc;
	local fh;
        description
        "OPen file to write computational data.",
        "fname ... name of the file",
        "header ... if true, print a csv header into the file",
        "discard ... if true, delete the file if it exists before opening a new file",
        "Return the file descriptor of the file.";
	if discard then
          if FileTools:-Exists(fname) then
            #printf("%s found\n", fname);
            fclose(fname);
	    FileTools:-Remove(fname);
          end if;
	end if;
	fh := fopen(fname, APPEND, TEXT);
	if header then
	  fprintf(fh, "\"name\",\"numVars\",\"degrees\",\"numParams\",\"paramDegrees\",\"extractor1\",\"extractor2\",\"rows\",\"columns\",\"rank\",\"memory\",\"time\"");
	  fflush(fh);
	end if;
	return fh;
      end proc;


      SaveStart := proc(fh::filedesc, t::table)
        description
        "Save initial data.",
        "fh ... file descriptor of output file",
        "t ... table wit data";
	fprintf(fh, "\n\"%s\",%d,\"%s\",%d,\"%s\"",
          t["system"]["name"],
	  t["system"]["numVars"],
	  convert(t["system"]["degrees"], string),
	  t["system"]["numParams"],
	  convert(t["system"]["paramDegrees"], string));
	fflush(fh);
      end proc;


      SaveEnd := proc(fh::filedesc, t::table)
	local extr2;
        description
        "Save remaining data.",
        "fh ... file descriptor of output file",
        "t ... table wit data";
	if nops(t["DixonExtraction"]["extractor"])>1 then
	  extr2 := t["DixonExtraction"]["extractor"][2];
	else
	  extr2 := "None";
	end if;

	fprintf(fh, ",\"%s\",\"%s\",%d,%d,%d,%.3f,%.6f",
	  t["DixonExtraction"]["extractor"][1], 
          extr2,
	  t["DixonMatrix"]["rows"],
	  t["DixonMatrix"]["columns"],
	  t["DixonMatrix"]["rank"],
	  t["DixonExtraction"]["memory"],
	  t["DixonExtraction"]["time"]);
	fflush(fh);
      end proc;


      System := proc(sysName::string, numVars::integer, degrees, numParams::integer:=undefined, paramDegrees:=undefined)
	local measures;
	description
	"Polynomial system information.",
        "numVars ... number of variables",
	"degrees ... list of degrees or degree vectors of the system",
        "numParams ... number of parameters",
	"paramDegrees ... list of degrees or degree vectors of the parameters of the system",
	"Return a table with the inormation."; 
	measures["system"] := table();
        measures["system"]["name"] := sysName;
	measures["system"]["numVars"] := numVars;
	measures["system"]["degrees"] := degrees;
	if numParams=undefined then
	  measures["system"]["numParams"] := 0;
        else
	  measures["system"]["numParams"] := numParams;
	end if;
	if paramdDgrees=undefined then
	  measures["system"]["paramDegrees"] := 0;
        else
	  measures["system"]["paramDegrees"] := paramDegrees;
	end if;
	return measures;
      end proc;


      RunMeasure := proc(fh::filedesc, fun::{procedure, function}, sysName::string, f::list, deg::{list, integer}, vars::list, paramDeg::integer, params::list)
	local t, dr;
        description
        "Run a function and store performance measures in a file.",
        "fh ... file descriptor of csv file to append peformance measures, according to SaveStart and SaveEnd",
        "sysName ... name of the polynomial system",
        "fun ... function to run, usually Dixon resultant computation",
        "f ... polynomial system",
        "deg ... degree list or degree in the variables of the system",
        "vars ... list of variables of the system",
        "paramDeg ... total degree of the parameters of the system",
        "params ... list of parameters of the system",
        "Return a table with the performace measure, as prduced by DRes:-ComposeDixon.";
	printf("%s\n", convert(fun, string)); 
	t := DR:-Info:-System(sysName, nops(vars), deg, nops(params), paramDeg): 
	DR:-Info:-SaveStart(fh, t):
	dr := fun(f, vars, t):
	DR:-Info:-SaveEnd(fh, t):
	printf("%.6f\n", t["DixonExtraction"]["time"]); 
	unassign('dr'):
	gc();
	return eval(t);
      end proc:


      DixonMatrix := proc(measures::{undefined, name, table}:=undefined)
	description
	"Gather and store information about the Dixon  matrix.",
	"measures ... table where the information about the Dixon matrix will be stored: number of rows and columns and rank",
	"Return a function with input matrix M that returns M.";
	return 
	  proc(M::Matrix)
	    description
	    "M ... Dixon matrix",
	    "Return the matrix M.";
	    if measures<>undefined then
	      measures["DixonMatrix"]["rows"] := LinearAlgebra:-RowDimension(M);
	      measures["DixonMatrix"]["columns"] := LinearAlgebra:-ColumnDimension(M);
	      measures["DixonMatrix"]["rank"] := Rank(M);
	    end if;
	    return M;
	  end proc;
      end proc;


      Extraction := proc(extractorList::list, measures::{undefined, name, table}:=undefined)
	description
	"Store the extraction time and the name of the function used for Dixon resultant extraction.",
	"extractorList ... list of functions that are composed to extract the Dixon resultant from the Dixon matrix.",
	"measures ... table where the information about the Dixon matrix will be stored: number of rows and columns and rank",
	"Return a function that returns its own input."; 
	return
	  proc()
	    description
	    "Return the sequence of inputs to this function.";
	    local e, extractorName;
	    if measures<>undefined then
	      measures["DixonExtraction"] := table();
	      measures["DixonExtraction"]["time"] := 0.0;
	      measures["DixonExtraction"]["memory"] := 0.0;
	      measures["DixonExtraction"]["extractor"] := [];
	      for e in extractorList do
		extractorName := convert(e, string);
		measures["DixonExtraction"]["time"] := measures["DixonExtraction"]["time"] + measures[extractorName]["time"];
		measures["DixonExtraction"]["memory"] := measures["DixonExtraction"]["memory"] + measures[extractorName]["memory"];
		measures["DixonExtraction"]["extractor"] := [op(measures["DixonExtraction"]["extractor"]), extractorName];
		measures["DixonExtraction"][extractorName] := table(measures[extractorName]);
		measures[extractorName] := evaln(measures[extractorName]);
	      end do;
	    end if;
	    return _rest;
	  end proc;
      end proc;

      Measure := proc(fun, measures::{undefined, name, table}:=undefined)
	description
	"fun ... procedure to be timed and whose byte usage is to be collected",
	"Measure the execution time of a function.",
	"measures ... output parameter, a table. If measures is a name then a table will be assigned to it. If measures is undefined then no timing will be produced.",
	"other paramters ... will be passed as inputs to fun",
	"Return a function that returns the output of fun and assigns the timing to the table measures with key convert(fun, string).";
	
	return
	  proc()
	    local result, usageResults, usedTime, usedBytes, funName;
	    if measures=undefined then
	      result := fun(_passed);
	    else
	      usageResults := CodeTools:-Usage(fun(_passed), output = ['cputime', 'bytesused', 'output'], quiet = true);
	      usedTime := usageResults[1];
	      usedBytes := usageResults[2];
	      result := usageResults[3..nops([usageResults])];
	      funName := convert(fun, string);
	      measures[funName]["time"] := usedTime/60.0;
	      measures[funName]["memory"] := usedBytes/(1024.0*1024.0);
	    end if;
            #print([result]);
	    return result;
	  end proc;
      end proc;


      Rank := proc(M::Matrix, tries::integer:=10, params::list(name):=[])
	local rows, cols, maxseq, count, varLst, rk, ultimateRank, Ms;
	description
	"Determine the rank of a parametric matrix.",
	"M ... matrix with parametric entries", 
	"tries ... number of times the matrix parameters should be instantiated with random numbers and instantiated matrix should have the same rank",
	"params ... optional list of parameters of matrix, if absent then automatically determined",
	"Return the rank (with high probability) of the parametric matrix M";
	(rows, cols) := LinearAlgebra:-Dimension(M);

	maxseq := tries;

	if params=[] then
	  varLst := indets(M);
	else
	  varLst := params;
	end if;

	randomize();
	ultimateRank := -1;
	count := 0;
	while count<=maxseq do
	  Ms := Extract:-InstantiateMatrix(M, varLst);
	  #printf("Ms; %a\n", Ms);
	  rk := LinearAlgebra:-Rank(Ms);
	  #printf("rk: %a\n", rk);
	  if rk <> ultimateRank then
	    ultimateRank := rk;
	    count := 0;
	  else
	    count := count + 1;
	  end if;
	end do;

	return ultimateRank;
      end proc:


    end module;

    ###########################################################################################################################
    Extract := module()
      export RankSubMatrix, MapleDet, MaxMinorElim, MaxMinorElimFull, InstantiateMatrix;
      local notFullRankMsg, IsNonZeroRow, PermutationFromMatrix, 
        ColNumersFactors, ColDenomsFactors, RowNumersFactors, RowDenomsFactors;
      description 
      "Various functions used in extracting the Dixon resultant from the Dixon matrix.",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";


      RankSubMatrix := proc(M::Matrix, {parameters::list:=[]})
	local rows, cols, varLst, substLst, i, j, Ms, Um, Pm, perm, indepRows, indepCols;
	description
	"Extract a maximal-rank minor of a parametric matrix.",
	"M ... parametric matrix", 
	"parameters ... parameters of M",
	"Return a maximal-rank minor of the parametric matrix M with high likelihood. If params is empty, the parameters will be automatically determined.";
	# Author: M. Minimair, 5/18/2015
	# Based on function RSC by A. Chtcherba, 2003
	
	(rows, cols) := LinearAlgebra:-Dimension(M);

	if parameters = [] then
	  varLst := indets(M);
	else
	  varLst :=  [op(parameters)];
	end if;

        #printf("varLst: %a\n", varLst);
        #printf("Rank of M: %d\n", Info:-Rank(M));
	Ms := InstantiateMatrix(M, varLst);
        #printf("Rank of Ms: %d\n", LinearAlgebra:-Rank(Ms));	

	Pm, Um := LinearAlgebra:-LUDecomposition(Ms, output=['P', 'U']);

	#printf("Rank of Um: %d\n", LinearAlgebra:-Rank(Um));

	perm := PermutationFromMatrix(Pm);
	
	indepRows := [];
	indepCols := [];
	for i from 1 to rows do
	  j := IsNonZeroRow(Um, i, i);
	  if j>0 then
	    indepRows := [op(indepRows), perm(i)];
	    indepCols := [op(indepCols), j];
	  else
	    #Matrix is upper triangular
	    break;
	  end if;
	end do;
	
	#return indepRows, indepCols;
	return LinearAlgebra:-SubMatrix(M, indepRows, indepCols);
      end proc:



      notFullRankMsg := "Dixon sub-matrix does not have full rank. Rerun the computation.\n";


      MapleDet := proc(A::Matrix)
        local n, m, det;
	description
	"Compute the determinant of a matrix with Maple's determinant function.",
	"A ... given matrix",
	"Return lowerRank::string, the determinant of A. non-empty message iff the determinant is 0 and A is not empty.";
        n, m := LinearAlgebra:-Dimension(A);
        if n=0 or m=0 then
          return "", 0;
        else
	  det := LinearAlgebra:-Determinant(A);
          if det=0 then
            return notFullRankMsg, 0;
          else
            return "", det;
          end if;
        end if;
      end proc;


      ColNumersFactors := proc(B, d, r, n, c, m)
	local jj, g, ii, updated;
	description
	"Extract common factors from numerators of enries in matrix columns."
	"B ... matrix to work on",
	"d ... partial determinant",
	"r ... first row to work on",
	"n ... number of rows of B",
	"c ... first column to work on",
	"m ... last column to work on",
	"Return d multiplied with the common factors in factored form.";
	updated := d;
	for jj from c to m do
	  g := 0;
	  for ii from 1 to n do
	    #printf("ii: %d, jj: %d\n", ii, jj);
	    #if ii=2 and jj=2 then
	    #  print(B[ii, jj]);
	    #end if;
	    g := gcd(g, numer(B[ii, jj]));
	  end do;
	  if g=0 then
	    next;
	  end if;
	  g := factor(g);
	  for ii from r to n do
	    B[ii, jj] := normal(B[ii, jj]/g);
	  end do;
	  updated := updated * g;
	end do;
	return updated;
      end proc;


      ColDenomsFactors := proc(B, d, r, n, c, m)
	local jj, g, ii, updated;
	description
	"Extract common factors from denominators of enries in matrix columns."
	"B ... matrix to work on",
	"d ... partial determinant",
	"r ... first row to work on",
	"n ... number of rows of B",
	"c ... first column to work on",
	"m ... last column to work on",
	"Return d multiplied with the common factors in factored form.";
	updated := d;
	for jj from c to m do
	  g := 0;
	  for ii from 1 to n do
	    #printf("ii: %d, jj: %d\n", ii, jj);
	    #if ii=2 and jj=2 then
	    #  print(B[ii, jj]);
	    #end if;
	    g := gcd(g, denom(B[ii, jj]));
	  end do;
	  if g=0 then
	    next;
	  end if;
	  g := factor(g);
	  for ii from r to n do
	    B[ii, jj] := normal(B[ii, jj]*g);
	  end do;
	  updated := updated / g;
	end do;
	return updated;
      end proc;


      RowNumersFactors := proc(B, d, r, n, c, m)
	local ii, g, jj, updated;
	description
	"Extract common factors from numerators of enries in  a matrix row."
	"B ... matrix to work on",
	"d ... partial determinant",
	"r ... row to work on",
	"n ... last row to work on",
	"c ... first column to work on",
	"m ... last column to work on",
	"Return d multiplied with the common factors in factored form.";
	updated := d;
	for ii from r to n do
	  g := 0;
	  for jj from c to m do
	    g := gcd(g, numer(B[ii, jj]));
	  end do;
	  if g<>0 then
	    g := factor(g);
	    for jj from c to m do
	      B[ii, jj] := normal(B[ii, jj]/g);
	    end do;
	    updated := updated * g;
	  end if;
	end do;
	return updated;
      end proc;


      RowDenomsFactors := proc(B, d, r, n, c, m)
	local ii, h, jj, updated;
	description
	"Extract common factors from denominators of enries in  a matrix row."
	"B ... matrix to work on",
	"d ... partial determinant",
	"r ... row to work on",
	"n ... last row to work on",
	"c ... first column to work on",
	"m ... last column to work on",
	"Return d divided by the common factors in factored form.";
	updated := d;
	for ii from r to n do
	  h := 0;
	  for jj from c to m do
	    h := gcd(h, denom(B[ii, jj]));
	  end do;
	  if h<>0 then
	    h := factor(h);
	    for jj from c to m do
	      B[ii, jj] := normal(B[ii, jj]*h);
	    end do;
	    updated := updated / h;
	  end if;
	end do;
	return updated;
      end proc;


      PermutationFromMatrix := proc(P::Matrix)
	local size, perm, i, j;
	description
	"Generate a permutation function from a permutation matrix.",
	"P ... permutation matrix",
	"Returns a procedure that implements the permutation given by P";

	size := LinearAlgebra:-RowDimension(P);

	for i from 1 to size do
	  for j from 1 to size do
	    if P[i,j]<>0 then
	      perm(j):=i;
	      break;
	    end if;
	  end do;
	end do;

	return eval(perm);
      end proc:


      InstantiateMatrix := proc(M::Matrix, varLst)
	local substLst, p, v, i, j, Ms, rows, cols;
	description
	"Randomly instantiate a parametric matrix."
	"M ... matrix with parametric entries",
	"varLst ... list of variables in M",
	"Return M with parameters varLst randomly specialized.";
	(rows,cols) := LinearAlgebra:-Dimension(M);
	
	# Create Substitution List
	substLst := []:
	#p := 137;
	for v in varLst do 
	  #p := nextprime(p);
	  p := rand(1..10000)();
	  #p := rand(1..100)();
	  substLst := [op(substLst),v=p]: 
	end do:

	# Create Numeric Matrix
	Ms := Matrix(rows,cols);
	for i from 1 to rows do
	  for j from 1 to cols do
	    Ms[i,j] := subs(substLst,M[i,j]);
	  end do:
	end do:

	#return eval(Ms);
	return Ms;
      end proc:


      IsNonZeroRow := proc(M::Matrix, rowIndex, startColIndex)
	local startIndex, j, cols, nonZeroIndex;
	description
	"Find first non-zero entry in matrix row.",
	"M", 
	"rowIndex ... index of row",
	"startColIndex ... optional column index",
	"Return col index of first non-zero entry of M in row rowIndex; 0 iff row is zero",
	"If column index startColIndex is present then start checking row from index startColIndex.";
	if _npassed>2 then
	  startIndex := startColIndex;
	else
	  startIndex := 1;
	end if;

	nonZeroIndex := 0;
	cols := LinearAlgebra:-ColumnDimension(M);
	for j from startIndex to cols do
	  if M[rowIndex, j]<>0 then
	    nonZeroIndex := j;
	    break;
	  end if;
	end do;
	return nonZeroIndex;
      end proc:



      MaxMinorElimFull := proc(A::Matrix)
	description
	"Compute a maximal-rank minor of a matrix via Gaussian elimination.",
	"A ... given matrix",
        "Named optional arguments: like MaxMinorElim, except tellFullrank which is not allowed.",
        "Return lowerRank::string, output from MaxMinorElim",
	"If lowerRank is empty then return a maximal minor of full rank. lowerRank is not empty iff A is not square of full rank and is not empty.";
        return MaxMinorElim(A, detFactor, tellFullrank=true, _options);
      end proc;


      MaxMinorElim := proc(A::Matrix,
                                     {tellFullrank::truefalse:=false, moreRows::truefalse:=true,
                                      preColNumers::truefalse:=false, preColDenoms::truefalse:=false, 
                                      preRowNumers::truefalse:=false, preRowDenoms::truefalse:=false, 
                                      colNumers::truefalse:=false, colDenoms::truefalse:=false, 
				      pivotRowNumers::truefalse:=false, pivotRowDenoms::truefalse:=false,
				      restRowNumers::truefalse:=false, restRowDenoms::truefalse:=false})
	local  n, m, tmp_n, rankCount, r, c, d, dr, B, i, j, t, det;
	description
	"Compute a maximal-rank minor of a matrix via Gaussian elimination.",
	"A ... given matrix",
        "Named optional arguments:",
        "return_fullrank ... indicate whether a non-zero minor of full rank was found",
        "moreRows ... transpose A such that it has at least as many rows as columns",
        "*row* ... extract common factors of rows",
        "*col* ... extract common factors of columns",
        "*Numers ... extract common factors from numerators of matrix entries",
        "*Denoms ... extract common factors from denominators of matrix entries",
        "pre* ... extract common factors before the Gaussian elimination starts",
        "pivot* ... extract common factors from the pivot row during Gaussian elimination",
        "rest* ... extract common factors from rows other than the pivot row during Gaussian elimination",
	"Return lowerRank::string, a maximal-rank minor of A.", 
        "If moreRows=false and rest* is activated then the output can be a multiple of a maximal-rank minor, because the output may contain factors of rows that are not needed for the determinant computation.";
	"If tellFullrank is true and lowerRank is empty then return a maximal minor of full rank. lowerRank is not empty iff A is not square of full rank and is not empty.",
        "If tellFullrank is false then lowerRank is empty and any maximal minor is returned.";
	#Original version: Hoon Hong (procedure drsm).
	#Modified: M. Minimair, 5/19/2015
	#          Various EDF options and other added.

        #print([_options]);

	(n, m) := LinearAlgebra:-Dimension(A);
        rankCount := 0;

        if n=0 or m=0 then
          return "", 0;
        end if;

	# Make sure that there are at least as many rows as columns.
        # It may result in fewer factors in the output if rest* is activiated.
	if n<m and moreRows then
	  B := LinearAlgebra:-Transpose(A);
          tmp_n := n;
          n := m;
          m := tmp_n;
	else
	  B := copy(A);
	end if;

	r := 1;
	d := 1;

	#Special case of Matrix([0])
	if m=1 and n=1 and B[1,1]=0 then
          if tellFullrank then
            return notFullRankMsg, 0;
          else
            return "", 0;
          end if;
	end if;

	#Other cases

        # Remove factors of the numerators of the matrix entries along the columns
        if preColNumers then
	    d := ColNumersFactors(B, d, 1, n, 1, m);
        end if;

        # Remove factors of the denominators of the matrix entries along the columns
        if preColDenoms then
	    d := ColDenomsFactors(B, d, 1, n, 1, m);
        end if;

        # Remove factors of the numerators of the matrix entries along the rows
        if preRowNumers then
	    d := RowNumersFactors(B, d, 1, n, 1, m);
        end if;

        # Remove factors of the denominators of the matrix entries along the rows
        if preRowDenoms then
	    d := RowDenomsFactors(B, d, 1, n, 1, m);
        end if;

	for c from 1 to m while r <= n do

	  if type(infolevel[DR], numeric) and infolevel[DR]>=2 then
	    printf("  %a: Working on column %d................\n", procname, c);
	  end if;

	  # Column numerator EDF
	  if colNumers then
	    d := ColNumersFactors(B, d, r, n, c, m);
	  end if;

	  # Column denominator EDF
	  if colDenoms then
	    d := ColDenomsFactors(B, d, r, n, c, m);
	  end if;

	  #--- Find a pivot ----
	  for i from r to n while B[i,c] = 0 do end do;
	  if i > n then
	    #printf(`    no pivot\n`);
	    next;
          else
            rankCount := rankCount + 1;
	  end if;
	  for j from i+1 to n do
	    if B[j,c] = 0 then next end if;
	    if length(B[j,c]) < length(B[i,c]) then i := j end if
	  end do;

	  #--- Swap rows ----
	  if i <> r then
	    for j from c to m do
	      t := B[i,j]; B[i,j] := B[r,j]; B[r,j] := t
	    end do
	  end if;
          # r is the row and c is the column index of the pivot
	 
	  # Pivot row numerator EDF
	  if pivotRowNumers then
	    d := RowNumersFactors(B, d, r, r, c, m);
	  end if;

	  # Pivot row denominator EDF
	  if pivotRowDenoms then
	    d := RowDenomsFactors(B, d, r, r, c, m);
	  end if;

	  # Rest row numerator EDF
	  if restRowNumers then
	    d := RowNumersFactors(B, d, r+1, n, c, m);
	  end if;

	  # Rest row denominator EDF
	  if restRowDenoms then
	    d := RowDenomsFactors(B, d, r+1, n, c, m);
	  end if;

	  #--- Update determinant ---
	  B[r,c] := factor(B[r,c]);
	  d := d * B[r,c];
	  #print(B[r,c]);
	  #print(`-----------------`);
	  #print(d);

	  #--- Eliminate ---
	  for i from r+1 to n do
	    if B[i,c] = 0 then next end if;
	    t := normal(B[i,c]/B[r,c]);
	    for j from c+1 to m do
	      B[i,j] := normal(B[i,j]-t*B[r,j]);
	    end do;
	    B[i,c] := 0
	  end do;
	  r := r+1
	end do;

	det := normal(d);
        if tellFullrank then
          if rankCount=n and rankCount=m then
            return "", det;
          else
            return notFullRankMsg, det;
          end if;
        else
          return "", det;
        end if;
      end;

    end module;

    ###########################################################################################################################
    EDF := module()
      export ColNumers, ColDenoms, Col, PivotRowNumers, PivotRowDenoms, PivotRow, ColPivotRow, Row, ColRow,
             RowFull, ColRowFull, PivotRowFull, ColPivotRowFull;
      description 
      "Versions of maximal minor with EDF computation (Early Detection of Factors).",
      "EDF from: Lewis, Rober H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
      "All functions work accordingly:",
      "Compute a maximal-rank minor of a matrix via Gaussian elimination.",
      "A ... given matrix",
      "Return a maximal-rank minor of A.",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";


      ColNumers      := proc(A) return Extract:-MaxMinorElim(A, colNumers=true) end proc;
      ColDenoms      := proc(A) return Extract:-MaxMinorElim(A, colDenoms=true) end proc;
      Col            := proc(A) return Extract:-MaxMinorElim(A, colNumers=true, colDenoms=true) end proc;
      PivotRowNumers := proc(A) return Extract:-MaxMinorElim(A, pivotRowNumers=true) end proc;
      PivotRowDenoms := proc(A) return Extract:-MaxMinorElim(A, pivotRowDenoms=true) end proc;
      PivotRow       := proc(A) return Extract:-MaxMinorElim(A, pivotRowNumers=true, pivotRowDenoms=true) end proc;
      ColPivotRow    := proc(A) return Extract:-MaxMinorElim(A, colNumers=true, colDenoms=true, pivotRowNumers=true, pivotRowDenoms=true) end proc;
      Row            := proc(A) return Extract:-MaxMinorElim(A, pivotRowNumers=true, pivotRowDenoms=true,  
                                                       restRowNumers=true, restRowDenoms=true) end proc;
      ColRow         := proc(A) return Extract:-MaxMinorElim(A, preColNumers=true, preRowNumers=true,
                                                       colNumers=true, colDenoms=true, pivotRowNumers=true, pivotRowDenoms=true,
                                                       restRowNumers=true, restRowDenoms=true) end proc;

      PivotRowFull       := proc(A) return Extract:-MaxMinorElimFull(A, pivotRowNumers=true, pivotRowDenoms=true) end proc;
      ColPivotRowFull    := proc(A) return Extract:-MaxMinorElimFull(A, colNumers=true, colDenoms=true, pivotRowNumers=true, pivotRowDenoms=true) end proc;
      RowFull            := proc(A) return Extract:-MaxMinorElimFull(A, pivotRowNumers=true, pivotRowDenoms=true,  
                                                       restRowNumers=true, restRowDenoms=true) end proc;
      ColRowFull         := proc(A) return Extract:-MaxMinorElimFull(A, preColNumers=true, preRowNumers=true,
                                                       colNumers=true, colDenoms=true, pivotRowNumers=true, pivotRowDenoms=true,
                                                       restRowNumers=true, restRowDenoms=true) end proc;

    end module;


    ###########################################################################################################################
    DRes := module()
      export ComposeDixon, MaxMinor, RankMaxMinor, RankMapleDet, RankFactorsMapleDet,
        MaxMinorColNumers, MaxMinorColDenoms, MaxMinorCol, MaxMinorPivotRowNumers,  MaxMinorPivotRowDenoms, MaxMinorPivotRow, MaxMinorColPivotRow, MaxMinorColRow,
        RankRow, RankColRow, RankMaxMinorPivotRow, RankMaxMinorColPivotRow, GetFirst, DetectFail;
      description 
      "Versions of Dixon resultant computation.",
      "All functions work accordingly:",
      "  pols ... list of (d+1) polynomials",
      "  vars ... list d variables of pols"
      "  measures ... output parameter, a table. If measures is a name then a table will be assigned to it. If measures is undefined then no timing will be produced.",
      "  Return the Dixon resultant of pols.",
      "  If randomized maximal-rank submatrix extraction is used, it is possible, however extremely unlikely, that Dixon resultant computation fails.",
      "  In this case, the computation terminates with an error message indicating to run the computation again or not to use randomized maximal-rank submatrix extraction.",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";


      GetFirst := proc()
	description
	"Return the first argument passed to the function."; 
	return _passed[1] 
      end proc;


      DetectFail := proc(msg::string)
        description
        "Detect failure case.",
        "msg ... error message",
        "Return rest arguments or, if msg is a non-empty string with the error message msg.";
        if msg<>"" then
          error msg;
        end if;
        return _rest;
      end proc;


      ComposeDixon := proc(extractor:=Extract:-MaxMinorElim, {measures_::{undefined, name, table}:=undefined})
	local getFirst, infoDixonMatrix, extractorList, extractionTime, DixonProc;
	description
	"Compose a function that computes the Dixon resultant.",
	"extractor ... function, with Dixon matrix input, to extract the Dixon resultant, or list of functions that will be composed",
	"measures ... output parameter, a table. If measures is a name then a table will be assigned to it. If measures is undefined then no timing will be produced.",
	"Return the timed Dixon resultant function.";


	extractorList := [extractor, op([_rest])];

	DixonProc := `@`(DetectFail,
                         Info:-Extraction(extractorList, measures_), seq(Info:-Measure(e,measures_), e in extractorList), 
			 Info:-DixonMatrix(measures_), GetFirst, Info:-Measure(DixonMatrix,measures_), 
			 Info:-Measure(DixonPolynomial,measures_));
	return DixonProc;
      end proc;


      MaxMinor := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction.";
        return ComposeDixon(Extract:-MaxMinorElim, measures_=measures)(pols, vars);
      end proc;


      RankMaxMinor := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Maximal-rank submatrix extraction and MaxMinor function.";
	return ComposeDixon(Extract:-MaxMinorElimFull, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      RankMapleDet := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Maximal-rank submatrix extraction and Maple determinant function.";
	return ComposeDixon(Extract:-MapleDet, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      RankFactorsMapleDet := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Maximal-rank submatrix extraction and Maple determinant function with PreEDF.";
	return ComposeDixon(Extract:-MapleDet, Extract:-DetectFactors, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      MaxMinorColNumers := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColNumersFactors detection.";
	return ComposeDixon(EDF:-ColNumers, measures_=measures)(pols, vars);
      end proc;


      MaxMinorColDenoms := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColDenomsFactors detection.";
	return ComposeDixon(EDF:-ColDenoms, measures_=measures)(pols, vars);
      end proc;


      MaxMinorCol := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColDenomsFactors detection.";
	return ComposeDixon(EDF:-Col, measures_=measures)(pols, vars);
      end proc;


      MaxMinorPivotRowNumers := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with PivotRowNumersFactors detection.";
	return ComposeDixon(EDF:-PivotRowNumers, measures_=measures)(pols, vars);
      end proc;


      MaxMinorPivotRowDenoms := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with PivotRowDenomsFactors detection.";
	return ComposeDixon(EDF:-PivotRowDenoms, measures_=measures)(pols, vars);
      end proc;


      MaxMinorPivotRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with PivotRowNumersFactors and PivotRowDenomsFactors detection.";
	return ComposeDixon(EDF:-PivotRow, measures_=measures)(pols, vars);
      end proc;


      MaxMinorColPivotRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColNumersFactors, ColDenomsFactors, PivotRowNumersFactors and PivotRowDenomsFactors detection.";
	return ComposeDixon(EDF:-ColPivotRow, measures_=measures)(pols, vars);
      end proc;


      MaxMinorColRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColNumersFactors, ColDenomsFactors, PivotRowNumersFactors and PivotRowDenomsFactors detection.";
	return ComposeDixon(EDF:-ColRow, measures_=measures)(pols, vars);
      end proc;


      RankRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with Pivot/RestRowNumersFactors and Pivo/ResttRowDenomsFactors detection on RankSubMatrix";
	return ComposeDixon(EDF:-RowFull, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      RankColRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Gaussian elimination based minor extraction with ColNumersFactors, ColDenomsFactors, Pivot/RestRowNumersFactors and Pivot/RestRowDenomsFactors detection on RankSubMatrix.";
	return ComposeDixon(EDF:-ColRowFull, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      RankMaxMinorPivotRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Maximal-rank submatrix extraction and MaxMinor function with PivotRowFactors detection.";
	return ComposeDixon(EDF:-PivotRowFull, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


      RankMaxMinorColPivotRow := proc(pols::list(polynom), vars::list(symbol), measures::{undefined, name, table}:=undefined)
	description
	"Maximal-rank submatrix extraction and MaxMinor function with Col/PivotRowFactors detection.";
	return ComposeDixon(EDF:-ColPivotRowFull, Extract:-RankSubMatrix, measures_=measures)(pols, vars);
      end proc;


    end module;


    ###########################################################################################################################
    GenPoly := module()
      export RandParamTotalDeg, RandParamMultiDeg;
      local NextTuple, FirstTuple, IsTuple, ToList, PowerBasis, GenPolynom;
      description 
      "Generate polynomials.",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";

      RandParamTotalDeg := proc(deg::{integer, list(integer)}, vars::list, paramDeg::integer, params::list, {sparse::truefalse:=true})
	local cfun;
	description
	"Random parametric total-degree integer polynomial.",
	"deg ... total degree or list of total degrees",
	"vars ... list of variables",
	"paramDeg ... total degree of the coefficient polynomials in the parameters",
	"params ... list of parameters",
        "sparse ... if true generate sparse polynomials and sparse parametric coefficients; otherwise dense."
	"Return a random dense integer polynomial of total degree deg, or a list of such polynomials if deg is a list with corresponding degrees. The degrees of the parameters are of total degree paramDeg. If parmDeg=0 or params=[] then create random ineger coefficients.";
	if type(deg, integer) then
	  #return randpoly([op(vars), op(params)], dense, degree = deg);
	  if paramDeg=0 or params=[] then
	    cfun := proc() rand(-99..99)() end proc;
	  else
            if sparse then
  	      cfun := proc() randpoly(params, degree=paramDeg) end proc;
            else
  	      cfun := proc() randpoly(params, degree=paramDeg, dense) end proc;
            end if;
	  end if;
          if sparse then
	    return randpoly(vars, coeffs=cfun, degree=deg);
          else
	    return randpoly(vars, coeffs=cfun, degree=deg, dense);
          end if;
	else
	  return map(d->RandParamTotalDeg(d, vars, paramDeg, params, sparse), deg);
	end if;
      end proc;


      NextTuple := proc(N::list, K::list)
	local n, zero, J, i;
	description
	"Next tuple in tuple iteration.",
	"N, K ... lists of non-negative integers",
	"Return next tuple after K <= N in lexicographic order. Return [0,..,0] if K = []. Return [] if there is no tuple <= N, i.e. if K = N";
	if nops(N)=0 then
	  return K;
	end if;

	n := nops(N);
	zero := [seq(0,i=1..n)];

	if K = [] then
	  return zero;
	end if;

	J := [seq(K[i], i=1..n)];
	#print(J);
	for i from n to 1 by -1 do
	  #print(i);
	  if J[i]<N[i] then
	    J[i] := J[i]+1;
	    break;
	  else
	    J[i] := 0;
	  end if;
	end do;

	if J = zero then
	  J := [];
	end if;

	return J;
      end proc:

      FirstTuple := proc(N::list)
	description
	"First tuple.",
	"N ... maximum tuple",
	"Return first tuple for NextTuple iteration.";
	return NextTuple(N, []);
      end proc:

      IsTuple := proc(L)
	description
	"Check if tuple.",
	"L ... tuple for tuple iteration",
	"Return true iff L is a valid tuple for NextTuple iteration.";
	return evalb(L <> []);
      end proc:

      ToList := proc(x)::list;
	local L;
	description
	"Wrap in list if not list.",
	"x ... to be wrapped in list",
	"Return [x] if x is not a list, otherwise return x.";
	if type(x, list) then 
	  L := x;
	else
	  L := [x];
	end if;
	return L;
      end proc:

      PowerBasis := proc(Ns,
      Js, Vs)
	local b, i, n, N, J, V;
	description
	"Power basis computation"
	"Ns ... maximal multi-degree of the basis. Not really needed, however, kept for symmetry.",
	"Js ... multi-degree of the basis",
	"Vs ... list of basis variables",
	"Return power basis of degree verctor Js in variables Vs.";    
	N := ToList(Ns);
	J := ToList(Js);
	V := ToList(Vs);
	n := nops(N);
	b := 1;
	for i from 1 to n do
	  b := b * V[i]^J[i];
	end do;
	return b;
      end proc:

      GenPolynom := proc(degs, cfun, vars)::polynom;
	local N, vs, i, J, p;
	description
	"Generate a polynomial.",
	"degs is a degree or list of degrees",
	"cfun is a function list(integer)->coefficient that maps an exponent tuple to a coefficient vars is a variable, variable-prefix or list of variables",
	"Return a polynomial of max degrees degs with coefficients generated by cfun, if vars is not a list then it is treated as a prefix from which the variables are formed, if degs is just one number then it is repeated according to the number of variables";
	if type(degs, list) and not type(vars, list) then
	  N := degs;
	  if nops(N) >= 2 then
	    vs := [seq(cat(vars, i), i=1..nops(N))];
	  else
	    vs := [vars];
	  end if;
	elif not type(degs, list) and type(vars, list) then
	  vs := vars;
	  N := [seq(degs, i=1..nops(vs))];
	elif not type(degs, list) and not type(vars, list) then
	  vs := [vars];
	  N := [degs];
	else
	  N := degs;
	  vs := vars;
	end if;

	J := FirstTuple(N);
	p := 0;
	while IsTuple(J) do
	  p := cfun(J) * PowerBasis(N, J, vs) + p;
	  J := NextTuple(N, J);
	end do;
	return p;
      end proc:


      RandParamMultiDeg := proc(degs, vars, paramDeg, params)::polynom;
	local cfun, J;
	description
	"Random dense parametric multi-degree integer polynomial.",
	"degs ... multi-degree list or list of multi-degree lists",
	"vars ... list of variables",
	"paramDeg ... total degree of the coefficient polynomials in the parameters",
	"params ... list of parameters",
	"Return a random dense integer polynomial of multi-degree degs, or a list of such polynomials if degs is a list with corresponding multi-degrees. The degrees of the parameters are of total degree paramDeg. If parmDeg=0 or params=[] then create random ineger coefficients.";
	if paramDeg=0 or params=[] then
	  cfun := J -> rand(-99..99)();
	else
	  cfun := J -> randpoly(params, dense, degree = paramDeg);
	end if;
	if type(degs, list) then
	  return GenPolynom(degs, cfun, vars);
	else
	  return map(d->RandParamMultiDeg(d, vars, paramDeg, params), deg);
	end if;
      end proc:


    end module;


    ###########################################################################################################################
    Samples := module()
      export Get, repo, Available;
      local AddSystem;
      description 
      "Collection of sample polynomial systems.",
      "To iterate over repo use: for sys in op(op(repo)) do ... Get(lhs(sys)) ... end do",
      "or iterator over the list Available()",
      "Author: Manfred Minimair",
      "E-mail: manfred@minimair.org",
      "Released under GPL Version 3",
      "Copyright (C) Manfred Minimair, 2015";


      # Table of polynomial systems
      repo := table();


      AddSystem := proc(sysName::string, polys::list(polynom), degrees::list, variables::list, paramDegree::integer, parameters::list, citation::{string, list(string)})
	description
	"Add a sample polynomial system."
	"sysName ... name of the system",
	"polys ... list of the polynomials",
	"degrees ... total or multi-degree list; nops(polys)=nops(degrees)",
	"variables ... list of variables of the polynomials",
        "paramDegree ... total degrees of the paramters",
	"parameters ... list of parameters of the polynomials",
	"citation ... citation string of the system";
	repo[sysName]["system"] := polys;
	repo[sysName]["degrees"] := degrees;
	repo[sysName]["variables"] := variables;
	repo[sysName]["paramDegree"] := paramDegree;
	repo[sysName]["parameters"] := parameters;
	repo[sysName]["citation"] := citation;
	return;
      end proc;
	

      Get := proc(sysName::string)
	description
	"Retrieve a polynomial system for use in Dixon resultant computation",
	"Return name, system, degrees, variables, param degree, parameters";
	return sysName, repo[sysName]["system"], repo[sysName]["degrees"], repo[sysName]["variables"][1..nops(repo[sysName]["variables"])], repo[sysName]["paramDegree"], repo[sysName]["parameters"];
      end proc;

      
      Available := proc({include::{undefined, list(string), string}:=undefined, exclude::list(string):=[]})::list;
        local sel, inLst, sys, sysName;
        description
        "List names of available sample systems.",
        "paramters are sequence of arbitrary strings; list of names of systems to be considered; if undefined, then consider all in repo",
        "include ... list of names of systems or single name to be considered; if undefined, then consider all in repo",
        "exclude ... list of names of systems to be excluded be default",
        "Return a list of system names available.";
        inLst := [_rest];
        if inLst=[] and type(include, undefined) then
          inLst := undefined;
        elif type(include, list) then
          inLst := [op(inLst), op(include)];
        else
          inLst := [op(inLst), include];
        end if;
        sel := [];
        for sys in op(op(repo)) do
          #printf("sys: %a\n", sys);
          sysName := lhs(sys);
          #printf("sysName: %s\n", sysName);
          if (type(inLst, undefined) or sysName in inLst) and not (sysName in exclude) then
            sel := [op(sel), sysName];
            #printf("sel: %a\n", sel);
         end if;
       end do;
       return sort(sel);
     end proc;


#     AddSystem(
#	"Bricard",
#	[sa^2 + ca^2 - 1,
#	 (dx - cx)^2 + (dy - cy)^2 - s4^2,
#	 sg^2 + cg^2 - 1,
#	 (ix - hx)^2 + (iy - hy)^2 - s6^2,
#	 sb^2 + cb^2 - 1,
#	 (fx - gx)^2 + (fy - gy)^2 - s5^2],
#	[2, 2, 2, 2, 2, 2, 2],
#	[ca, sa, cb, sb, cg, sg],
#        2,
#	[b, e, s1, s2, s3, s4, s5, s6, s7, s8, s9],
#	"Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
#     );


      AddSystem(
        "Bricard",
        [a1*t1^2*t2^2 + b1*t1^2 + 2*c1*t1*t2 + d1*t2^2 + e1,
          a2*t2^2*t3^2 + b2*t2^2 + 2*c2*t2*t3 + d2*t3^2 + e2,
          a3*t1^2*t3^2 + b3*t1^2 + 2*c3*t1*t3 + d3*t3^2 + e3],
        [4, 2, 2],
        [t1, t2],
        3,
        [t3, a1, b1, c1, d1, e1, a2, b2, c2, d2, e2, a3, b3, c3, d3, e3],
        "Lewis, Robert H., https://home.bway.net/lewis/dixon"
      );


      AddSystem(
        "SB L1 M5 K1",
        [2*mm1^2-2*mm2,
         6*mm2^2+2*mm4-8*mm1*mm3,
         20*mm3^2-2*mm6+12*mm1*mm5-30*mm2*mm4,
         70*mm4^2+2*mm8-12*mm1*mm7+56*mm2*mm6-112*mm3*mm5,
         252*mm5^2-2*mm10+20*mm1*mm9-90*mm2*mm8+240*mm3*mm7-420*mm4*mm6-mm3+mm1*mm2,

         -315+14496*mm1+23912*mm3-9310*mm4+8*mm7-196*mm6+1904*mm5-30184*mm2,
         2*mm8-728*mm6+9408*mm5-51632*mm4+141120*mm3-185152*mm2+91392*mm1-2205,
         4*mm9-17052*mm6+247380*mm5-1445010*mm4+4105160*mm3-5529048*mm2+2784096*mm1-72765,
         mm10-43407*mm6+670320*mm5-4070200*mm4+11869200*mm3-16288944*mm2+8326080*mm1-231525],
         [2, 2, 2, 2, 2, 2, 1, 1, 1],
         [mm1, mm2, mm3, mm4, mm5, mm6, mm7, mm8],
         1,
         [mm9, mm10],
         ["Little, John B., Solving the Selesnick–Burrus filter design equations using computational algebra and algebraic geometry, Advances in Applied Mathematics 31 (2003) 463–500", "Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"]
       );


      AddSystem(
        "Cyclic 6-Root",
        [x1          + x2          + x3          + x4       + x5       + x6,
         x1*x2       + x2*x3       + x3*x4       + x4*x5    + x5*x6    + x6*x1,
         x1*x2*x3    + x2*x3*x4    + x3*x4*x5    + x4*x5*x6 + x5*x6*x1 + x6*x1*x2,
         x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x1 + x5*x6*x1*x2 + x6*x1*x2*x3,
         x1*x2*x3*x4*x5 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x1 + x4*x5*x6*x1*x2 + x5*x6*x1*x2*x3 + x6*x1*x2*x3*x4,
         x1*x2*x3*x4*x5*x6 - 1],
         [1, 2, 3, 4, 5],
         [x1, x2, x3, x4, x5],
         6,
         [x6],
	 "Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
       );


      AddSystem(
        "Cyclic 7-Root",
        [x1          + x2          + x3          + x4       + x5       + x6 + x7,
         x1*x2       + x2*x3       + x3*x4       + x4*x5    + x5*x6    + x6*x7 + x7*x1,
         x1*x2*x3    + x2*x3*x4    + x3*x4*x5    + x4*x5*x6 + x5*x6*x7 + x6*x7*x1 + x7*x1*x2,
         x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x7 + x5*x6*x7*x1 + x6*x7*x1*x2 + x7*x1*x2*x3,
         x1*x2*x3*x4*x5    + x2*x3*x4*x5*x6    + x3*x4*x5*x6*x7    + x4*x5*x6*x7*x1    + x5*x6*x7*x1*x2    + x7*x1*x2*x3*x4,
         x1*x2*x3*x4*x5*x6 + x2*x3*x4*x5*x6*x7 + x3*x4*x5*x6*x7*x1 + x4*x5*x6*x7*x1*x2 + x5*x6*x7*x1*x2*x3 + x6*x7*x1*x2*x3*x4 + x7*x1*x2*x3*x4*x5,
         x1*x2*x3*x4*x5*x6*x7 - 1],
         [1, 2, 3, 4, 5, 6],
         [x1, x2, x3, x4, x5, x6],
         7,
         [x7],
	 "Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
       );


      AddSystem(
	"KK5",
	[-2*x2 + x1^2 - mu,
	 -2*x3 + x2^2 - mu,
	 -2*x4 + x3^2 - mu,
	 -2*x5 + x4^2 - mu,
	 -2*x1 + x5^2 - mu,
	 x1*x2*x3*x4*x5 + 1],
	[2, 2, 2, 2, 2],
	[x1, x2, x3, x4, x5],
        1,
	[mu],
	"Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
       );


       AddSystem(
	"KK6",
	[-2*x2 + x1^2 - mu,
	 -2*x3 + x2^2 - mu,
	 -2*x4 + x3^2 - mu,
	 -2*x5 + x4^2 - mu,
	 -2*x6 + x5^2 - mu,
	 -2*x1 + x6^2 - mu,
	 x1*x2*x3*x4*x5*x6 + 1],
	[2, 2, 2, 2, 2, 2],
	[x1, x2, x3, x4, x5, x6],
        1,
	[mu],
	"Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
      );


       AddSystem(
	"KK7",
	[-2*x2 + x1^2 - mu,
	 -2*x3 + x2^2 - mu,
	 -2*x4 + x3^2 - mu,
	 -2*x5 + x4^2 - mu,
	 -2*x6 + x5^2 - mu,
	 -2*x7 + x6^2 - mu,
	 -2*x1 + x7^2 - mu,
	 x1*x2*x3*x4*x5*x6*x7 + 1],
	[2, 2, 2, 2, 2, 2, 2],
	[x1, x2, x3, x4, x5, x6, x7],
        1,
	[mu],
	"Lewis, Robert H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
      );

       AddSystem(
	"KK8",
	[-2*x2 + x1^2 - mu,
	 -2*x3 + x2^2 - mu,
	 -2*x4 + x3^2 - mu,
	 -2*x5 + x4^2 - mu,
	 -2*x6 + x5^2 - mu,
	 -2*x7 + x6^2 - mu,
	 -2*x8 + x7^2 - mu,
	 -2*x1 + x8^2 - mu,
	 x1*x2*x3*x4*x5*x6*x7*x8 + 1],
	[2, 2, 2, 2, 2, 2, 2, 2],
	[x1, x2, x3, x4, x5, x6, x7, x8],
        1,
	[mu],
	"Lewis, Rober H., Heuristics to accelerate the Dixon resultant, Mathematics and Computers in Simulation, 77, 4 (2008), 400-407"
      );


      AddSystem(
	"ParamElim",
	[vu^2 - 2*x*vu + v^2 - 2*y*v + x^2 + y^2 - 1,
	 v^2 - vu^3,
	 -3*vu^2*v + 3*y*vu^2 - 2*v*vu + 2*x*v,
	 6*w^2*vu^2*v - 3*vu^2*w - 2*w*v + 1],
	[[2, 2, 0], [3, 2, 0], [2, 1, 0], [2, 1, 2]],
	[vu, v, w],
        2,
	[x, y],
	"Lewis, Robert H., Comparing acceleration techniques for the Dixon and Macaulay resultants, Mathematics and Computers in Simulation, 80, 6 (2010), 1146-1152"
      );


      AddSystem(
        "Enneper",
        [3*vss - vss^3 + 3*vss*vtt^2 - 3*x,
         -3*vtt - 3*vss^2*vtt + vtt^3 - 3*y,
         vss^2 - vtt^2 - z],
        [3, 3, 2],
        [vss, vtt],
        1,
        [x, y, z],
        "Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );



      AddSystem(
	"bicubic",
	[-x + 3*tt*(tt-1)^2 + (vss-1)^3 + 3*vss,
	-y + 3*vss*(vss-1)^2 + tt^3 + 3*tt,
	-z - 3*vss*(vss^2-5*vss+5)*tt^3 - 3*(vss^3+6*vss^2-9*vss+1)*tt^2
	      + tt*(6*vss^3+9*vss^2-18*vss+3) - 3*vss*(vss-1)],
	[3, 6, 5],
	[vss,tt],
	1,
	[x, y, z],
	"Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );


      AddSystem(
	"cubic",
	[vtt^3 - 3*vss*vtt - vss^3 - vss + x,
	vtt*vss^2 - 3*vtt + 1 - y,
	2*vtt^3 - 5*vss*vtt + vtt - vss^3 - z],
	[4, 2, 3],
	[vss, vtt],
	1,
	[x,y,z],
	"Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );

      AddSystem(
	"sec11_2_2",
	subs([
	X = 2*vtt^3 + 4*vtt^2 + 2*vtt + 4*vss*vtt + vss^2*vtt + 2 + 3*vss + vss^2,
	Y = -2*vss*vtt^2 - 2*vtt - vss*vtt + 2 + vss - 2*vss^2 - vss^3,
	Z = 2*vtt^2 - 3*vss*vtt^2 - 2*vtt - 3*vss*vtt - 2*vss^2*vtt - 2*vss - 3*vss^2 - vss^3,
	W = vtt^3 + vtt^2 - vtt + vss^2*vtt - 1 - vss + vss^2 + vss^3],
	[x*W - X, y*W - Y, z*W - Z]),
	[3, 3, 3, 3],
	[vss, vtt],
	1,
	[x, y, z],
	"Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );

      AddSystem(
	"sphere",
	[x*(vrr^2+1) - (vrr^2-1),
	y*(vrr^2+1)*(vrr^2+vss^2) - 4*vrr^2*vss,
	z*(vrr^2+1)*(vrr^2+vss^2) - 2*vrr*(vrr^2-vss^2)],
	[2, 2, 2],
	[vrr, vss],
	1,
	[x, y, z],
	"Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );


      AddSystem(
	"stophoid",
	subs(S = C*T, [
	x-a*S,
	y-a*T*(1+S),
	C^2+S^2-1
	]),
	[2, 3, 4],
	[C, T],
	1,
	[x, y],
	"Kapur, Deepak, and Manfred Minimair. Multivariate resultants in Bernstein basis. Automated Deduction in Geometry. Springer Berlin Heidelberg, 2011. 60-85."
      );


      # Random sparse systems

      AddSystem(
        "sparse_3v2_1p2",
        [87*p1-56*p2+(-62*p1+97*p2-73)*x2+(-4*p1-83*p2-10)*x1^2+(62*p1-82*p2+80)*x1*x2+(-44*p1+71*p2-17)*x1^2*x2+(-75*p1-10*p2-7)*x1*x2^2, -23*p1+87*p2+44+(29*p1+98*p2-23)*x2+(10*p1-61*p2-8)*x1^2+(-29*p1+95*p2+11)*x1*x2+(-49*p1-47*p2+40)*x2^2+(-81*p1+91*p2+68)*x1^3, (p1+55*p2-28)*x1+(16*p1+30*p2-27)*x1^2+(-15*p1-59*p2-96)*x2^2+(72*p1-87*p2+47)*x1^3+(-90*p1+43*p2+92)*x1^2*x2+(-91*p1-88*p2-48)*x1*x2^2],
        [3, 3, 3],
        [x1, x2],
        1,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_4v2_1p2",
        [71*p1+16*p2+83+(9*p1-60*p2-83)*x1+(98*p1-48*p2-19)*x1^3+(62*p1+37*p2+5)*x1^2*x2+(96*p1-17*p2+25)*x1*x2^2+(91*p1+98)*x2^3, (-47*p1-39*p2-53)*x1+(-72*p1-97*p2+33)*x1^2+(10*p1+7*p2-89)*x1^3+(65*p1+12*p2-25)*x1^2*x2+(-96*p1+50*p2-60)*x2^3+(-42*p1+7*p2-89)*x1^4, (69*p1+80*p2+28)*x1^2+(-42*p1-33*p2+21)*x2^2+(-35*p1+97*p2+30)*x1^3+(-64*p1+89*p2-16)*x1^2*x2+(59*p1-69*p2-46)*x1^4+(-33*p1+87*p2-34)*x1*x2^3],
        [4, 4, 4],
        [x1, x2],
        1,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_5v2_1p2",
        [(18*p1+52*p2+36)*x2+(91*p1-22*p2+51)*x1^2+(-27*p1+50*p2+60)*x1*x2+(-91*p1-47*p2-97)*x1^3*x2+(-2*p1-31*p2+25)*x2^4+(31*p1-27*p2+65)*x1^4*x2, (73*p1+95*p2+68)*x2+(-29*p1+5*p2-26)*x1*x2^2+(-51*p1+88*p2+97)*x1*x2^3+(-67*p1+58*p2+29)*x1^5+(37*p1+5*p2-36)*x1^2*x2^3+(-57*p1+85*p2+80)*x1*x2^4, (65*p1-12*p2+78)*x1*x2^2+(5*p1-63*p2-5)*x1^3*x2+(36*p1-8*p2+30)*x1^2*x2^2+(-3*p1-56*p2-91)*x1*x2^3+(-70*p1+42*p2+9)*x1^3*x2^2+(-21*p1-27*p2-79)*x2^5],
        [5, 5, 5],
        [x1, x2],
        1,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_6v2_1p2",
        [(-31*p1+45*p2+49)*x1*x2^3+(-58*p1+49*p2+1)*x2^4+(-95*p1+86*p2-97)*x1^4*x2+(-14*p1+83*p2-96)*x1*x2^4+(-8*p1-54*p2+62)*x1^4*x2^2+(96*p1-51*p2+89)*x2^6, -2*p1+86*p2+57+(-35*p1+57*p2+28)*x1*x2+(63*p1+21*p2-71)*x2^3+(-66*p1-34*p2+72)*x1^3*x2^2+(-40*p1-68*p2-15)*x2^5+(-32*p1+17*p2+87)*x1^6, (-50*p1-22*p2+48)*x1^3+(-82*p1+46*p2-49)*x1^2*x2+(-66*p1+18*p2+16)*x1^4+(-22*p1+38*p2-48)*x1^3*x2+(-87*p2-13)*x1^4*x2+(-87*p1-63*p2-31)*x2^6],
        [6, 6, 6],
        [x1, x2],
        1,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_1v2_1p3",
        [(-7*p1+22*p2-55*p3-94)*x1+(87*p1-56*p2-62)*x2+97*p1-73*p2-4*p3-83, (-10*p1+62*p2-82*p3+80)*x1+(-44*p1+71*p2-17*p3-75)*x2-10*p1-7*p2-40*p3+42, (-50*p1+23*p2+75*p3-92)*x1+(6*p1+74*p2+72*p3+37)*x2-23*p1+87*p2+44*p3+29],
        [1, 1, 1],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_2v2_1p3",
        [-16*p1-33*p2-54*p3-65+(30*p1+37*p3+35)*x1+(14*p1+77*p2+53*p3+82)*x2+(-58*p1+49*p2+7*p3+68)*x1^2+(72*p1-54*p2-10*p3-48)*x1*x2+(-30*p1+32*p2+68*p3-47)*x2^2, -91*p1-58*p2-2*p3+35+(-10*p1-43*p2+46*p3-69)*x1+(68*p1-37*p2+42*p3-40)*x2+(32*p1-27*p2-22)*x1^2+(-61*p1+28*p2+32*p3-69)*x1*x2+(-60*p1+58*p2-49*p3+22)*x2^2, 82*p1-62*p2+50*p3+46+(88*p1+99*p2-12*p3+27)*x1+(65*p1-32*p2-83*p3-86)*x2+(-42*p1-54*p2-6)*x1^2+(-11*p1-6*p2-6*p3-32)*x1*x2+(-88*p1+58*p2+81*p3+25)*x2^2],
        [2, 2, 2],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_3v2_1p3",
        [-70*p1+6*p2-40*p3+44+(25*p1+98*p2+92*p3+71)*x2+(91*p1+16*p2+10*p3+37)*x1^2+(-45*p1-30*p2+66*p3-17)*x1*x2+(43*p1+54*p2+4*p3+41)*x1^3+(89*p1+74*p2+59*p3-66)*x1*x2^2, 87*p1-25*p2+70*p3+26+(46*p1-25*p2-50*p3-91)*x2+(57*p1+53*p2+87*p3+68)*x1*x2+(-17*p1-17*p2-71*p3+44)*x2^2+(71*p1+5*p2+10*p3+89)*x1^3+(91*p1+3*p2-21*p3+21)*x1^2*x2, -51*p1-79*p2+55*p3+58+(-60*p1+65*p2-50*p3+45)*x1+(60*p1+64*p2+69*p3-9)*x2+(-33*p1-85*p2+85*p3-65)*x2^2+(-48*p1-17*p2-62*p3-92)*x1^2*x2+(-57*p1-63*p2+59*p3+28)*x2^3],
        [3, 3, 3],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_4v2_1p3",
        [(37*p1-30*p2+59*p3-32)*x2+(12*p1+42*p2-9*p3+40)*x1^3+(38*p1+43*p2+6*p3-80)*x1*x2^2+(63*p1+11*p2+9*p3-49)*x2^3+(-67*p1+43*p2-25*p3-37)*x1^4+(91*p1-53*p2-13*p3+38)*x1*x2^3, -90*p1+17*p2+26*p3+52+(-52*p1-28*p2-35*p3+15)*x1+(-11*p1+72*p2-45*p3-52)*x2^2+(-9*p1+66*p2+70)*x1^3+(7*p1-32*p2-5)*x1^3*x2+(5*p1+37*p2-40*p3-82)*x1*x2^3, 72*p1+7*p2+36*p3-94+(17*p1-55*p2+81*p3+29)*x2+(-90*p1-85*p2-97*p3+11)*x1^2+(-97*p1+22*p2-59*p3-60)*x2^3+(-11*p1+2*p2+59*p3-51)*x1^3*x2+(14*p1+53*p2+46*p3-41)*x1*x2^3],
        [4, 4, 4],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_5v2_1p3",
        [(92*p2-47*p3-58)*x1^2+(-72*p1+43*p2-94*p3-74)*x1*x2^2+(-48*p1+20*p2+83*p3-25)*x1^4+(55*p1-53*p2+29*p3+98)*x1*x2^3+(49*p1-72*p2+43*p3+9)*x1^5+(-50*p1-5*p2-97*p3+33)*x1^4*x2, -12*p1-12*p2-82*p3-80+(-p1+66*p2-50*p3-22)*x2+(p1+95*p2-68*p3+18)*x2^3+(-72*p1+21*p2-57*p3+71)*x1^3*x2+(54*p1+56*p2+52*p3-52)*x2^4+(50*p1+75*p2-16*p3-99)*x1*x2^4, (76*p1+9*p2+53*p3+55)*x2+(-63*p1+8*p2-54*p3+2)*x1^2+(-57*p1+83*p2+70*p3+18)*x1*x2+(34*p1-74*p2-22*p3-51)*x1^3+(-93*p1+77*p2+12*p3-71)*x1*x2^3+(-50*p1-75*p2-84*p3+74)*x1^4*x2],
        [5, 5, 5],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_6v2_1p3",
        [(54*p1-41*p2+68*p3-76)*x2^2+(62*p1+97*p2-65*p3+34)*x1^2*x2+(67*p1-88*p2-64*p3+45)*x1^4+(77*p1-93*p2-58*p3+90)*x1^2*x2^3+(95*p1+10*p2+43*p3-16)*x1^3*x2^3+(-71*p1+38*p2-82*p3+18)*x1^2*x2^4, (9*p1-13*p2-14*p3+49)*x2+(-45*p1-23*p2+97*p3+45)*x1^2*x2+(-43*p1-87*p2+73*p3+11)*x1*x2^2+(57*p1+5*p2+40*p3+19)*x1^5+(-79*p1-76*p2-22*p3+78)*x1^5*x2+(-91*p1-41*p2+23*p3-69)*x1^4*x2^2, (65*p1+14*p2+61*p3+56)*x1*x2^2+(-91*p1+63*p2+68*p3+43)*x2^3+(-28*p1+22*p2-84*p3+63)*x1^4*x2+(19*p1-57*p2+95*p3-75)*x1^2*x2^3+(2*p1+9*p2+58*p3+27)*x1^6+(-87*p1+12*p2-95*p3+56)*x1^4*x2^2],
        [6, 6, 6],
        [x1, x2],
        1,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_1v2_2p2",
        [(-94*p1^2+87*p1*p2-56*p2^2+22*p1-55*p2-7)*x1+(-73*p1^2-4*p1*p2-83*p2^2-62*p1+97*p2)*x2+80*p1^2-44*p1*p2+71*p2^2+62*p1-82*p2-10, (-7*p1^2-40*p1*p2+42*p2^2-75*p1-10*p2-17)*x1+(-92*p1^2+6*p1*p2+74*p2^2+23*p1+75*p2-50)*x2+87*p1^2+44*p1*p2+29*p2^2+37*p1-23*p2+72, (-61*p1^2-8*p1*p2-29*p2^2-23*p1+10*p2+98)*x1+(-47*p1^2+40*p1*p2-81*p2^2+11*p1-49*p2+95)*x2+31*p1^2-51*p1*p2+77*p2^2+68*p1-10*p2+91],
        [1, 1, 1],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_2v2_2p2",
        [40*p1^2-15*p1*p2+5*p2^2+30*p1+51*p2-59+(59*p1^2-82*p1*p2+75*p2^2+32*p1+17*p2+60)*x1+(63*p1^2+16*p1*p2-86*p2^2-80*p1+19*p2-8)*x2+(-43*p1^2-79*p1*p2-93*p2^2+22*p1+49*p2+82)*x1^2+(94*p1^2-43*p1*p2+43*p2^2-93*p1+99*p2-54)*x1*x2+(37*p1^2-44*p1*p2+33*p2^2-51*p1+55*p2+84)*x2^2, -23*p1^2-45*p1*p2+86*p2^2-78*p1-58*p2-85+(-86*p1^2-53*p1*p2-45*p2^2-15*p1-87*p2+21)*x1+(39*p1^2-77*p1*p2-58*p2^2-39*p1+58*p2+84)*x2+(55*p1^2-81*p1*p2+12*p2^2+39*p1-65*p2+35)*x1^2+(62*p1^2-16*p1*p2-44*p2^2-33*p1-84*p2-73)*x1*x2+(65*p1^2+80*p1*p2-65*p2^2-53*p1+40*p2-36)*x2^2, 7*p1^2+59*p1*p2+37*p2^2-3*p1-42*p2-33+(-74*p1^2-25*p1*p2-13*p2^2+61*p1-24*p2+27)*x1+(-40*p1^2+87*p1*p2-4*p2^2-15*p1+91*p2+34)*x2+(-64*p1^2+17*p1*p2-27*p2^2+31*p1+79*p2-77)*x1^2+(-59*p1^2-88*p1*p2+7*p2^2+20*p1-17*p2-7)*x1*x2+(28*p1^2+49*p1*p2+p2^2+46*p1-63*p2+4)*x2^2],
        [2, 2, 2],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_3v2_2p2",
        [99*p1^2-35*p1*p2+41*p2^2-6*p1-96*p2-73+(-36*p1^2+63*p1*p2+81*p2^2-52*p1-61*p2+94)*x1+(95*p1^2-63*p1*p2+28*p2^2-48*p1-90*p2-45)*x2+(27*p1^2+7*p1*p2+50*p2^2+26*p1-37*p2-34)*x1^2+(10*p1^2+13*p1*p2-14*p2^2-57*p1-79*p2+66)*x2^2+(-75*p1^2+16*p1*p2+78*p2^2+95*p1+83*p2+32)*x1^2*x2, -24*p1^2+41*p1*p2+64*p2^2-39*p1-27*p2+98+(32*p1^2-20*p1*p2+37*p2^2+24*p1+66*p2+48)*x2+(12*p1^2-40*p1*p2-98*p2^2+76*p1-p2-93)*x1^2+(-83*p1^2-46*p1*p2+12*p2^2-55*p1-28*p2+34)*x1*x2+(26*p1^2-69*p1*p2-58*p2^2+26*p1-8*p2+34)*x1*x2^2+(-43*p1^2-8*p1*p2+67*p2^2+66*p1-67*p2+67)*x2^3, (92*p1^2-76*p1*p2-p2^2+52*p1+22*p2+13)*x1+(90*p1^2-7*p1*p2-85*p2^2+98*p1+69*p2+48)*x2+(-62*p1^2-47*p2^2+29*p1+42*p2+69)*x1*x2+(99*p1^2+16*p1*p2-83*p2^2+23*p1-52*p2-99)*x1^2*x2+(-81*p1^2-51*p1*p2-35*p2^2-74*p1+95*p2+32)*x1*x2^2+(-61*p1^2-99*p1*p2+59*p2^2+10*p1+52*p2-31)*x2^3],
        [3, 3, 3],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_4v2_2p2",
        [(-62*p1^2+71*p1*p2+74*p2^2-41*p1+12*p2-13)*x1^2+(-77*p1^2+32*p1*p2+54*p2^2+17*p1+19*p2-36)*x1*x2+(93*p1^2+55*p1*p2-20*p2^2+80*p1-49*p2+79)*x1^2*x2+(2*p1^2-82*p1*p2-43*p2^2-17*p1-30*p2+27)*x2^3+(-41*p1^2+85*p1*p2+73*p2^2-5*p1+56*p2-56)*x1^2*x2^2+(83*p1^2-63*p1*p2+64*p2^2+18*p1+81*p2-13)*x1*x2^3, -41*p1^2+76*p1*p2-17*p2^2-59*p1-68*p2-39+(-73*p1^2+18*p1*p2+10*p2^2+70*p1-46*p2-56)*x2+(-15*p1^2+83*p1*p2+31*p2^2-60*p1+58*p2+26)*x1^3+(86*p1^2-84*p1*p2+38*p2^2-96*p1-27*p2+11)*x1^4+(-47*p1^2-80*p1*p2+27*p2^2-96*p1-72*p2+4)*x1^3*x2+(-60*p1^2+86*p1*p2-27*p2^2-11*p1-38*p2+71)*x1*x2^3, (38*p1^2+40*p1*p2-68*p2^2+90*p1-13*p2+48)*x1+(-9*p1^2+8*p1*p2-89*p2^2+94*p1+30*p2-58)*x2+(21*p1^2-85*p1*p2+9*p2^2+14*p1-40*p2+80)*x1*x2+(72*p1^2-29*p1*p2-12*p2^2+46*p1-86*p2-84)*x1*x2^2+(24*p1^2+44*p1*p2-18*p2^2+69*p1-17*p2-1)*x2^3+(29*p1^2-49*p1*p2-9*p2^2+96*p1-53*p2+62)*x2^4],
        [4, 4, 4],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_5v2_2p2",
        [(53*p1^2-68*p1*p2+79*p2^2-63*p1-68*p2+66)*x2+(-50*p1^2+80*p1*p2-11*p2^2+6*p1-34*p2-73)*x1^2+(18*p1^2-20*p1*p2-16*p2^2+35*p1+7*p2+7)*x1^3+(-8*p1^2-47*p1*p2-3*p2^2-11*p1-67*p2+13)*x1^2*x2^2+(20*p1^2-9*p1*p2+83*p2^2+92*p1+40*p2-75)*x1^5+(10*p1^2-72*p1*p2-6*p2^2+3*p1+16*p2-43)*x1^4*x2, -81*p1^2-25*p1*p2+p2^2+94*p1-50+(-81*p1^2-13*p1*p2+53*p2^2+63*p1+29*p2+75)*x2+(-48*p1^2-42*p1*p2+22*p2^2+97*p1+34*p2-48)*x1^2+(-58*p1^2-63*p1*p2+97*p2^2-48*p1-14*p2-87)*x1*x2+(-8*p1^2-99*p1*p2+29*p2^2+25*p1+78*p2+88)*x2^2+(28*p1^2-88*p1*p2-82*p2^2+40*p1+30*p2+88)*x2^3, (-17*p1^2+23*p1*p2-47*p2^2+63*p1-19*p2-89)*x1^3*x2+(-50*p1^2+75*p1*p2-95*p2^2-64*p1+47*p2-65)*x2^4+(-80*p1^2+48*p1*p2-9*p2^2-87*p1+40*p2-79)*x1^4*x2+(78*p1^2+15*p1*p2-81*p2^2-21*p1-71*p2-71)*x1^3*x2^2+(-46*p1^2-33*p1*p2-70*p2^2+93*p1-36*p2-86)*x1*x2^4+(-10*p1^2-29*p1*p2+p2^2+36*p1-57*p2+37)*x2^5],
        [5, 5, 5],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_6v2_2p2",
        [(31*p1^2+77*p1*p2-10*p2^2-50*p1-61*p2-44)*x1+(50*p1^2-55*p1*p2+11*p2^2-69*p1+33*p2+47)*x1^2+(-37*p1^2-53*p1*p2+77*p2^2+12*p1+26*p2+17)*x1^2*x2+(-33*p1^2-55*p1*p2-89*p2^2-65*p1+88*p2-44)*x1^3*x2+(-58*p1^2+14*p1*p2+43*p2^2-6*p1-70*p2+35)*x1^5+(13*p1^2+35*p1*p2+57*p2^2+46*p1-58*p2+88)*x1^2*x2^3, (-90*p1^2-24*p1*p2-97*p2^2-61*p1-76*p2+98)*x1^2+(41*p1^2+37*p1*p2-68*p2^2-90*p1-5*p2+48)*x1*x2^2+(-66*p1^2-59*p1*p2+49*p2^2-43*p1-9*p2+65)*x1^3*x2+(25*p1^2-51*p1*p2+30*p2^2-15*p1+5*p2-73)*x2^5+(53*p1^2-11*p1*p2-73*p2^2+23*p1-39*p2+64)*x1^5*x2+(78*p1^2+30*p1*p2-80*p2^2-87*p1+41*p2+18)*x2^6, (-7*p1^2-56*p1*p2+84*p2^2-35*p1+32*p2-46)*x1^2+(37*p1^2-47*p1*p2+34*p2^2-12*p1+81*p2+18)*x1*x2^2+(-89*p1^2+88*p1*p2-90*p2^2+30*p1+92*p2-93)*x1^3*x2+(-67*p1^2-49*p1*p2-99*p2^2-55*p1+83*p2+71)*x1*x2^3+(-67*p1^2+p1*p2-98*p2^2-43*p1+42*p2+95)*x1^4*x2+(-88*p1^2+15*p1*p2+76*p2^2+52*p1+68*p2+31)*x1^4*x2^2],
        [6, 6, 6],
        [x1, x2],
        2,
        [p1, p2],
        "random"
      );


      AddSystem(
        "sparse_1v2_2p3",
        [(-62*p1^2+97*p2^2-73*p2*p3-56*p2+87)*x1+(-44*p1*p2+71*p1*p3-17*p2*p3+62*p1-82*p2+80*p3)*x2+37*p1^2-23*p1*p2+87*p1*p3+74*p2+72*p3+6, (-47*p1^2+40*p2^2-81*p2*p3+91*p3^2+11*p1-49*p3)*x1+(p1^2+p1*p2+55*p1*p3-28*p2*p3+77*p1+95*p3)*x2+72*p1*p2-87*p1*p3+47*p2^2-90*p3^2-96*p2-59, (-10*p1^2-82*p1*p2+71*p2^2+5*p1+13*p2-28)*x1+(-19*p1^2+62*p1*p2+37*p2*p3+5*p3^2+98*p1-48*p2)*x2-13*p1^2+44*p1*p3-2*p2^2+71*p2*p3-34*p1-60],
        [1, 1, 1],
        [x1, x2],
        2,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_2v2_2p3",
        [43*p1^2-39*p1*p2-81*p1*p3-22*p2^2-p1-51+(-39*p1*p2+7*p1*p3+26*p2*p3-47*p3^2+47*p2+74)*x1+(61*p1^2-29*p1*p2+89*p1*p3+34*p3^2-64*p3-92)*x2+(88*p1^2-47*p1*p3-p2^2+60*p2*p3-76*p3^2-72)*x1^2+(43*p1^2-87*p2^2+19*p2*p3+24*p3^2+17*p1+83)*x1*x2+(-40*p1^2-95*p1*p3-93*p2*p3+35*p3^2+70*p2+3*p3)*x2^2, -87*p1*p2-22*p2^2-18*p2*p3+81*p3^2+45*p1-43*p3+(87*p1^2+50*p1*p3+43*p3^2+78*p1-68*p2+85)*x1+(-35*p2^2-95*p2*p3-76*p1-23*p2+34*p3+3)*x2+(-79*p1*p2-p1*p3-37*p2^2+24*p2*p3+35*p2+65*p3)*x1^2+(59*p1^2+32*p1*p3-99*p2^2-96*p2+43*p3-79)*x1*x2+(-22*p1^2-31*p1*p2+30*p2^2-18*p3^2-33*p3-71)*x2^2, 82*p1*p2-10*p2^2+78*p2*p3+83*p3^2-74*p1-32*p3+(-26*p1^2+77*p1*p3+52*p2*p3+96*p3^2+79*p1-53)*x1+(-40*p1^2-31*p1*p3-62*p2*p3-55*p3^2+35*p2+36*p3)*x2+(-89*p1^2+16*p1*p3+49*p3^2+34*p2+22*p3-30)*x1^2+(p1^2+35*p1*p3+56*p2*p3-89*p3^2-91*p1+64*p2)*x1*x2+(-58*p1^2-45*p1*p3+98*p2^2-77*p2*p3+28*p3^2+60*p1)*x2^2],
        [2, 2, 2],
        [x1, x2],
        2,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_3v2_2p3",
        [-20*p1^2-16*p2*p3-16*p3^2-20*p1+88*p3-72+(-69*p1*p2+8*p2^2-30*p1+43*p2+45*p3-1)*x1^2+(-87*p1*p2+42*p2^2+60*p1+35*p2-38*p3-83)*x1^3+(-61*p1*p2-83*p1*p3-43*p2^2-91*p1-66*p2+29)*x1^2*x2+(89*p1^2+98*p2^2+75*p2*p3+77*p3^2-57*p1+8*p3)*x1*x2^2+(-11*p1^2-21*p1*p2+29*p1-42*p2-53*p3-65)*x2^3, (-16*p1*p2-58*p1*p3+25*p3^2-56*p2+29*p3-83)*x1+(-31*p1^2-29*p1*p2-11*p2^2-20*p2*p3+21*p3^2-17*p2)*x1^2+(16*p1^2+59*p1*p2-97*p3^2+91*p2+16*p3+91)*x1*x2+(52*p1^2+97*p1*p3+14*p2^2+87*p1+99*p2+48)*x2^2+(-80*p1^2-49*p1*p2-17*p1*p3-59*p1+86*p2+22*p3)*x1^3+(95*p1^2+3*p1*p2+89*p1*p3-42*p2^2-p2*p3+56*p2)*x1^2*x2, -67*p1*p3-36*p2^2-74*p3^2-96*p1-99*p3-87+(61*p1^2-51*p1*p2+67*p2^2+78*p3^2+86*p3-6)*x1^2+(-66*p1*p2-79*p1*p3-72*p3^2-75*p1-71*p2-68)*x1*x2+(56*p1^2+44*p1*p3-27*p2*p3+81*p1-29*p3+16)*x2^2+(91*p1^2-59*p2*p3+11*p3^2+64*p2-28*p3-96)*x1*x2^2+(-57*p1^2+13*p1*p2+9*p1*p3-41*p2^2-26*p3^2+81*p1)*x2^3],
        [3, 3, 3],
        [x1, x2],
        2,
        [p1, p2, p3],
        "random"
      );


      AddSystem(
        "sparse_4v2_2p3",
        [(-72*p1^2-86*p1*p3-74*p2*p3+19*p3^2+45*p2-p3)*x2+(43*p1^2-32*p1*p2+71*p1*p3-74*p2^2+96*p1)*x1*x2+(4*p1*p3-94*p2^2-25*p3^2-93*p1-21*p3-61)*x1*x2^2+(59*p1^2-18*p1*p3+60*p2*p3+82*p1-36*p2-71*p3)*x1^4+(-77*p1*p3-50*p2^2-69*p1-44*p2-25*p3-45)*x1^3*x2+(-74*p1^2-12*p1*p2-32*p2^2-18*p2*p3-65*p3^2+6*p2)*x2^4, -93*p2*p3+72*p3^2-98*p1+93*p2+58*p3+7+(82*p1*p2-45*p2^2+35*p3^2-70*p1+9*p2-65)*x1^2+(43*p1^2-50*p1*p2+38*p1*p3-p3^2+24*p3+8)*x1*x2+(-45*p1^2+16*p2^2-98*p3^2+43*p2-36*p3+34)*x2^2+(61*p1*p3+42*p2^2-44*p3^2+25*p1-48*p2+49)*x1^4+(-2*p1^2-46*p1*p2-98*p2^2-23*p3^2+61*p2-38)*x1*x2^3, (-62*p1*p3+15*p2*p3-44*p3^2-58*p1-97*p2-21)*x1^2+(-65*p1^2+55*p1*p2-80*p1*p3+45*p1+24*p2-5)*x1*x2+(54*p1^2+17*p2*p3-21*p3^2-32*p1+62*p2+71*p3)*x1^3+(-54*p1^2+99*p1*p3-35*p2^2-27*p2*p3+39*p3+38)*x1^2*x2+(-90*p1*p3+19*p2*p3-7*p3^2+80*p1+74*p2-59)*x1*x2^2+(24*p1*p3+19*p2*p3-6*p3^2+62*p1+74*p2+88)*x1^2*x2^2],
        [4, 4, 4],
        [x1, x2],
        2,
        [p1, p2, p3],
        "random"
      );



      end module;


end module;

