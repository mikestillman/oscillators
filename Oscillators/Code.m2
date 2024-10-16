
--------------------------------------------------------------------
----- Internal Helper Functions
--------------------------------------------------------------------

hasNoLeaf = method()
hasNoLeaf Graph := G -> all(vertices G, i -> degree(G,i) > 1)


--------------------------------------------------------------------
----- Creation of main objects (OscRing, oscQuadrics)
--------------------------------------------------------------------

oscRing = method(Options => {
        Start => 1,
        Symbols => (getSymbol "x", getSymbol "y"),
        CoefficientRing => QQ
        })

oscRing(ZZ, List) := Ring => opts -> (n, params) -> (
    -- n: #angle functions, indexed 1..n by default
    -- params: a list of parameter symbols (or ring variables that can be used as variables in a new ring)
    -- output is a ring R
    -- R_0..R_(n-1): are the cosine functions
    -- R_n..R_(2n-1): are the sine functions
    -- R_(2n).. are the parameters
    kk := opts.CoefficientRing;
    s := opts#Symbols#1;
    c := opts#Symbols#0;
    start := opts#Start;
    R := if #params > 0
      then kk[c_start..c_(start+n-1), s_start..s_(start+n-1),params,MonomialOrder=>{2*n,#params}]
      else kk[c_start..c_(start+n-1), s_start..s_(start+n-1)];
    R.numOscillators = n;
    R
    )
oscRing ZZ := Ring => opts -> n -> oscRing(n, {}, opts)
oscRing Graph := Ring => opts -> G -> oscRing(# vertices G , opts)

-- Helper function to return the ideal on cos and sin for n variables.
trig = method()
trig Ring := (R) -> (
    if not R.?numOscillators then error "expected ring generated by 'oscRing'";
    n := R.numOscillators;
    ideal for i from 0 to n-1 list (R_i)^2+(R_(i+n))^2-1
    )

---------------------------------------------------
-- Jacobian Constructions
--------------------------------------------------

oscJacobian = method(Options => {Reduced=>false})
oscJacobian Ideal := opts -> P -> (
    -- P is a system created with oscSystem, or oscSystemReduced
    -- each polynomial should be <= linear in R_i, R_(n+i), all i=0..n-1
    R := ring P;
    P1 := ideal drop(P_*,1);
    if not R.?numOscillators then error "expected ring generated by 'oscRing'";
    n := R.numOscillators;
    matrix for i from 0 to n-1 list (
        first entries (R_i * diff(R_(n+i), gens P1) - R_(n+i) * diff(R_i, gens P1))
        )
    )

oscJacobian(Graph, Ring) := opts -> (G, R) -> (
    -- let w(1), ..., w(n-1) be w(i) = theta(i) - theta(0).  These are the angles that are represented in R.
    -- let w(0) = theta_0 + ... + theta_(n-1).
    -- We want the Jacobian of the system after this change of coordinates.
    -- i.e. dw(i)/dt = 0 if i == 0
    --      dw(i)/dt = d theta(i) / dt - d theta(0) / dt
    --  in terms of w(i), i=1..n-1, this is
    --      P_i := dw(i)/dt = sum(j in N_i, j != 0) sin(w(j)-w(i)) - sum(j in N_i, j = 0) sin(w(i)) - sum(j in N_0) sin(w(j))
    --           (note: the second to last term is either 0, or one term).
    --      P_0 := 0;
    -- Jacobian: 
    --  d P_i / d w(j) = 0                                                (case if i or j == 0)
    --  d P_i / d w(j) = cos(w(j)-w(i)) - (if (i,j) = edge) cos(w(j))     (case if i, j >= 0 and i != j)
    --  d P_i / d w(i) = - sum(j in N_i, j != 0) cos(w(j)-w(i)) 
    --                   - (if (i,0) is an edge) 2 * cos(w(i))            (case: if i > 0)
    if not opts.Reduced then return fullJacobian(G,R);
    if not R.?numOscillators then error "expected ring generated by 'oscRing'";
    n := R.numOscillators;
    sine := (i) -> if i === 0 then 0_R else R_(n+i-1);
    cosine := (i) -> if i === 0 then 1_R else R_(i-1);
    N0 := neighbors(G,0);
    matrix for i from 1 to n list for ell from 1 to n list (
        Ni := neighbors(G,i);
        if i == ell then (
            -- do the diagonal element
            part1 := sum(sort toList Ni, j -> - cosine j * cosine i - sine j * sine i);
            part2 := if member(i, N0) then - cosine i else 0_R;
            part1 + part2
            )
        else (
            part1 = if member(ell, Ni) then cosine ell * cosine i + sine ell * sine i else 0_R;
            part2 = if member(ell, N0) then - cosine ell else 0_R;
            part1 + part2
            )
        )
    )    

oscJacobian(Graph) := opts -> G -> oscJacobian(G, oscRing(G), opts)


fullJacobian = method()
fullJacobian(Graph, Ring) = (G,R) -> (
    -- G is a grpah
    -- R is a ring created with oscRing(# vertices G, params)
    if not R.?numOscillators then error "expected ring generated by 'oscRing'";
    n := R.numOscillators;
    if # vertices G =!= n+1 then error("expected "|(n+1)|" oscillators");
    sine := (i) -> if i === 0 then 0_R else R_(n+i-1);
    cosine := (i) -> if i === 0 then 1_R else R_(i-1);
    --sine := (i) -> R_(n+i);
    --cosine := (i) -> R_i;
    matrix for i from 0 to n list for j from 0 to n list (
        if i =!= j then (
            if member(j, neighbors(G,i)) then (
                -- entry is cos(theta_j - theta_i)
                cosine j * cosine i + sine j * sine i
            ) else 0_R
        ) else 
          - sum for j in toList neighbors(G, i) list (cosine j * cosine i + sine j * sine i)
        )
    )

---------------------------------------------------
-- Ideal without weights, i.e. weights are all 1 --
---------------------------------------------------
oscSystemReduced(Graph, Ring) := (G, R) -> (
    n := # vertices G - 1;  -- R should have == 2n variables.
    if toList(0..n) =!= sort vertices G then
      error "got incompatible graph";
    sine := (i) -> if i === 0 then 0_R else R_(n+i-1);
    cosine := (i) -> if i === 0 then 1_R else R_(i-1);
    I := ideal for i from 0 to n list (
        N := toList neighbors(G,i);
        sum for j in N list (
            e := sort {i,j};
            (sine j * cosine i - sine i * cosine j)
            )
        );
    I + trig R
    )

oscQuadrics = method()
oscQuadrics(Graph, Ring) := (G, R) -> (
    -- R should have == 2n variables.
    n := # vertices G;
    if 2*n =!= numgens R then error("expected ring with "|2*n|" variables");
    if toList(0..(n-1)) =!= sort vertices G then
      error "got incompatible graph";
    sine := (i) -> R_(n+i);
    cosine := (i) -> R_i;
    I := ideal for i from 0 to n-1 list (
        N := toList neighbors(G,i);
        sum for j in N list (
            (- sine j * cosine i + sine i * cosine j)
            )
        );
    I
    )

oscQuadrics(Graph) := G -> oscQuadrics(G, oscRing(G))
    



----------------------------------------------------------------------------------
-- Finding solutions: beware, these might give duplicates and/or miss solutions --
----------------------------------------------------------------------------------
isReal List := (L) -> all(L, p1 -> imaginaryPart p1 == 0)

findRealSolutions = method()
findRealSolutions Ideal := (I) -> (
    -- I is an ideal over the complexes, which is zero-dimensional (might be higher dimensional though...!
    sols := solveSystem I_*;
    bads := positions(sols, p -> p.cache.SolutionStatus != Regular);
    if #bads > 0 then (<< "warning: some solutions are not regular: " << bads << endl);
    coords := matrix(sols/coordinates);
    coords = clean(.00000001, coords);
    select(entries coords, p -> all(p, p1 -> imaginaryPart p1 == 0))
    )
findRealSolutions(Graph) := (G) -> (
    n := # vertices G;
    RC := oscRing(n-1,{}, CoefficientRing => CC);
    IC := trim oscSystemReduced(G,RC);
    JC := oscJacobian(G, RC);
    pts := findRealSolutions IC;
    for p in pts list {p, sub(JC, matrix{p}), identifyStability(JC, p, Tolerance=>1e-10)}
    )

allUniquePrincipalMinors = method(Options => {Modulo=>null})
allUniquePrincipalMinors Matrix := opts -> (M) -> (
    I := if opts.?Modulo and opts.Modulo =!= null then opts.Modulo else 0_(ring M);
    S := drop(subsets numcols M, 1); -- drop the empty subset
    -- I is an ideal in a polynomial ring over e.g. QQ
    -- M is a symmetric n by n matrix over the sam ering as I.
    if ring M =!= ring I then error "expected ideal and matrix over the same ring";
    if M - transpose M != 0 then error "expected a symmetric matrix";
    M1 := M % I;
    unique for s in S list (
        (det submatrix(M1, s, s)) % I
        )
    )

identifyStability = method(Options => {Tolerance => 1e-15})
-- returns one of the symbols Stable, Unstable, Semistable
identifyStability(Matrix, List) := Symbol => opts -> (Jac,pt) -> (
    J0 := substitute(Jac, matrix{pt});
    ev := eigenvalues J0;
    identifyStability(opts, ev)
    )
identifyStability BasicList := Symbol => opts -> eigenvals -> (
    signs := for e in eigenvals list
        if abs(realPart e) < opts.Tolerance then 0 
        else if realPart e < 0 then -1 
        else 1;
    nzero := # select(signs, s -> s === 0);
    npos := # select(signs, s -> s === 1);
    nneg := # select(signs, s -> s === -1);
    assert(nzero + npos + nneg == # eigenvals);
    if npos > 0 then Unstable
    else if nzero == 1 then Stable
    else Semistable
    )
isStableSolution = method()
isStableSolution(Matrix, List) := (J,pt) -> (
    J0 := substitute(J, matrix{pt});
    ev := eigenvalues J0;
    all(ev, e -> realPart e < 0.0)
    )

-----------------------------
-- twisted angle points -----
-----------------------------
twisted = method()
twisted(ZZ, ZZ, Ring) := (p, n, R) -> (
    a := cos(p*pi/n);
    c := for i from 1 to n-1 list cos(1_R * i*p*2*pi/n);
    s := for i from 1 to n-1 list sin(1_R * i*p*2*pi/n);
    c | s
    )
twisted(ZZ, ZZ, ZZ) := (p, n, prec) -> twisted(p, n, RR_prec)

----------------------------------------------------------------------------
-- Useful display functions and methods to see results in terms of angles --
----------------------------------------------------------------------------
outputSolutions = method()
-- L is a double list of real numbers
outputSolutions List := (L) -> demark("\n", L/(p->demark("  ", p/toString)))

getAngles = method(Options => {Radians => true})
getAngles(ZZ,List) := opts -> (nangles, sols) -> (
    -- assume 'sols' is a list of
    -- (cos theta_1, cos theta_1, ..., cos theta_nangles, sin theta_1, ..., sin theta_nangles, ...)
    for p in sols list (
        for i from 0 to nangles-1 list (
            --a := atan2(realPart p#(i), realPart p#(i+nangles));
            a := atan2(realPart p#(i+nangles), realPart p#(i));
            b := if a < 0 then a + 2*pi else a;
            if opts.Radians then b else b * 180/pi
            )
        )
    )

getLinearlyStableSolutions = method()
getLinearlyStableSolutions Graph := G -> (
    <<  "---- doing graph " << G << " --------------" << endl;
    n := # vertices G;
    RC := oscRing(n-1,{}, CoefficientRing => CC);
    IC := trim oscSystemReduced(G,RC);
    --<< see IC << endl;
    JC := oscJacobian(G, RC);
    elapsedTime realsols := findRealSolutions IC;
    stablesols := for pt in realsols list (
        result := identifyStability(JC, pt, Tolerance=>1e-10);
        if result == Stable then pt else continue
        );
    if #stablesols > 1 then (
        << "-- extra stable solutions! --" << endl;
        << netList stablesols << endl;
        << netList (getAngles(n-1, stablesols, Radians=>false)) << endl;
        );
    stablesols
    )
