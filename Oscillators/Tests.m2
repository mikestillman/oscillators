
TEST ///
  -- A simple example to make sure this is working
-*
  restart
  needsPackage "Oscillators"
*-
  G = graph ({2, 0, 1, 3}, {{0, 2}, {1, 2}, {2, 3}, {0, 3}, {1, 3}})

  S = oscRing 3
  assert(numgens S == 6) -- basically: we set theta_0 to 0.

  IG = oscSystem(G,S)
  Ians = ideal(y_2+y_3,
      -x_2*y_1-x_3*y_1+x_1*y_2+x_1*y_3,
      x_2*y_1-x_1*y_2-x_3*y_2+x_2*y_3-y_2,
      x_3*y_1+x_3*y_2-x_1*y_3-x_2*y_3-y_3,
      x_1^2+y_1^2-1,
      x_2^2+y_2^2-1,
      x_3^2+y_3^2-1
      )
  assert(IG == Ians)
  assert(gens IG == gens Ians) -- checking that we did not remove any of the equations.
  
  assert(dim IG == 1)
  
  -- now let's make sure the Jacobian(s) of this system are OK.
  Jac = oscJacobian(G,S)
  assert(det Jac == 0)
  assert(rank Jac == 3)
  assert((numrows Jac, numcols Jac) == (4,4))
///

TEST ///
  -- The square
-*
  restart
  needsPackage "Oscillators"
*-
  G = graph ({2, 0, 1, 3}, {{0, 2}, {1, 2}, {2, 3}, {0, 3}, {1, 3}})
  G = graph( {0,1,2,3}, {{0,1}, {1,2}, {2,3}, {0,3}})

  S = oscRing 3
  assert(numgens S == 6) -- basically: we set theta_0 to 0.

  IG = oscSystem(G,S)
  Ians = ideal(y_1+y_3,
      -x_2*y_1+x_1*y_2-y_1,
      x_2*y_1-x_1*y_2-x_3*y_2+x_2*y_3,
      x_3*y_2-x_2*y_3-y_3,
      x_1^2+y_1^2-1,
      x_2^2+y_2^2-1,x_3^2+y_3^2-1
      )
  assert(IG == Ians)
  assert(gens IG == gens Ians) -- checking that we did not remove any of the equations.
  
  assert(dim IG == 1)
  
  -- now let's make sure the Jacobian(s) of this system are OK.
  Jac = oscJacobian(G,S)
  assert(det Jac == 0)
  assert(rank Jac == 3)
  assert((numrows Jac, numcols Jac) == (4,4))
  
  numericalIrreducibleDecomposition IG
///



TEST /// 
-*
  restart
  needsPackage "Oscillators"
*-
  -- The triangle
  G = graph({0,1,2},{{0,1},{1,2},{0,2}})
  R = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,R);
  I = ideal gens gb trim I
      C = decompose I
      standardI = ideal(y_1, y_2, x_1^2-1, x_2^2-1)
      I1 = I : standardI
      I2 = I : I1
      assert(I == intersect(I1, I2))
      assert(I2 == standardI)
      decompose I1
      assert(I == radical I)
      assert(dim I == 0)
      assert(#C == 5) -- 4 standard plus 1 more (which has 2 solutions)
    J = oscJacobian(G,R)
    Jacs = for i in C list J % i
    EVs = for j in Jacs list eigenvalues(lift(j, QQ)) -- self-synch solution is stable, all others are unstable.
    cleanedEVs = EVs/(evs -> evs/(x -> clean(1e-10, x)))
    stable = positions(cleanedEVs, x -> Stable == identifyStability x )
    semistable = positions(cleanedEVs, x -> Semistable == identifyStability x)    
    unstable = positions(cleanedEVs, x -> Unstable == identifyStability x)
    assert(#semistable == 0)
    assert(#unstable == 4)
    assert(#stable == 1)

    unstablesols = select(findRealSolutions I, p -> Unstable === identifyStability(J,p))
    tally ((findRealSolutions I)/(p -> identifyStability(J, p)))
    getAngles(2, unstablesols, Radians=>false)
    
    -- note: the 2 nonstandard solutiona are twisted solutions:
    -- theta_0 = 0, theta_1 = 2 pi/3, theta_2 = 4 pi/3
    -- theta_0 = 0, theta_1 = 4 pi/3, theta_2 = 8 pi/3
    -- upshot: ideal is radical zero-dim, 4 standard solutions, 2 twisted solutions (unstable), SS
///

///
-*
  restart
  needsPackage "Oscillators"
*-
  -- The square
  -- TODO: get the same info as in email...
  G = graph({0,1,2,3},{{0,1},{1,2},{2,3},{0,3}})
  R = oscRing(# vertices G - 1,{}, CoefficientRing => QQ)
  I = oscSystem(G,R);
  --I = ideal gens gb trim I
  C = decompose I
      assert(#C == 5)
      assert(C/dim == {1,1,1,0,0})
      assert(C/degree == {2,2,4,1,1})
      assert(I == intersect C) -- therefore a radical ideal
      
      -- the 3 codim 1 components all intersect at 2 (twisted) points:
      C01 = trim ideal gens gb (C#0 + C#1)
      C02 = trim ideal gens gb (C#0 + C#2)
      C12 = trim ideal gens gb (C#1 + C#2)
      C01 == C02
      C01 == C12
      decompose C01
      -- the 2 isolated components: are both standard points.
      standardI = ideal(y_1..y_3) + ideal(x_1^2-1, x_2^2-1, x_3^2-1)
      assert(C#3 : standardI == 1)
      assert(C#4 : standardI == 1)
    -- The isolated (smooth) points:
    J = oscJacobian(G,R)
    J4 = J % C_4
    eigenvalues(lift(J4, QQ)) -- all positive, theta_0 = theta_1 = 0, theta_2 = theta_3 = pi
    J3 = J % C_3
    eigenvalues(lift(J3, QQ)) -- all negative: all angles 0.
    -- now for the positive dimensional components:
    J0 = J % C_0
      eigenvalues lift(sub(J0, {x_3 => 1}),QQ) -- one positive, one negative, one zero.
        -- UNLESS: x_3 = 0, in which we get two solutions whose Jacobian is ZERO:
        -- (theta_1, theta_2, theta_3) = (pi, -pi/2  , pi/2), (pi, pi/2 , - pi/2) these have Jacobian == 0
    J1 = J % C_1
      eigenvalues lift(sub(J1, {x_3 => 1}),QQ) -- one positive, one negative, one zero.
      -- same points as on J0...
    J2 = J % C_2
      eigenvalues lift(sub(J2, {x_3 => 1}),QQ) -- one positive, one negative, one zero.
      C2 = decompose trim(C_2 + ideal(x_3))
      assert(#C2 == 2)
      assert(C2/degree == {1,1})
      assert(C2/dim == {0,0})
      -- gives two points with all non-negative eigenvalues:
      J2a = J % C2_0 -- ZERO, point is (pi, -pi/2, pi/2)
      J2b = J % C2_1  
    -- upshot:
    -- 1 stable solution angles = (0,0,0).
    -- 1 unstable solution: (0,pi,pi)
    -- 3 families of 1-d solutions, each of which consists of unstable solutions, together
    --   with 2 points which have zero jacobian (it is the same 2 points for each component).
    --   The intersection of these three components consists of two points:
    --   (pi, -pi/2, pi/2)
    --   (pi, pi/2, -pi/2)
    -- of the 3 1D components:
    --   (pi, -theta, theta).
    --   (pi, theta + pi, theta)
    --   (, , ) -- theta_2 = theta_3 + pi.
    
    
///

TEST ///
  -- analyze stability of all solutions given a graph G
  G = graph({0,1,2,3},{{0,1},{1,2},{2,3},{0,3}})
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- by looking at each one, all points but one are unstable.
///

TEST ///
-*
  restart
  needsPackage "Oscillators"
  
  needsPackage "Visualize" 
  openPort "8080"
  visualize Gs_3 -- 
  closePort()
*-    
  needsPackage "NautyGraphs"
  Gstrs = flatten for i from 4 to 9 list generateGraphs(5,i, OnlyConnected=>true)
  Gs = Gstrs/stringToGraph
  assert(#Gs == 20)
  Gs = select(Gs, G -> all(vertices G, v -> degree(G,v) > 1))
  assert(#Gs == 10)  

  -- Gs_0 -- pentagon  
  G = Gs_0
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim I == 0, 30 points, radical.
  -- 3 stable points.  All the rest are not stable
  -- C_-1: 4 real solutions (actually the 4 complex solutions are all real)
  realsols = findRealSolutions C_-1
  for x in realsols list identifyStability(J,x)
  -- 2 of these points are stable.  
  -- total of 3 stable points.
  
  -- Gs_1
  G = Gs_1
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  netList allUniquePrincipalMinors(-J, Modulo=>C_0)
  I0 = trim(C_0 + ideal(y_4^2-1, x_1*x_4 + x_2*x_4 + y_1*y_4 + y_2*y_4 + x_4))
  netList for i in decompose I0 list allUniquePrincipalMinors(-J, Modulo=>i)  
  -- Gs_1: one component of dim 2, 8 components of dim 0, radical.  All points except one are unstable.
  -- all points on C_0 are UNSTABLE.
  -- the stable point is one of the 0 dim components, so a smooth point.

  -- Gs_2
  G = Gs_2
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- Gs_2: 26 components, all of dim 0, radical.
  --  36 points total. all unstable except SS point.

  -- Gs_3
  G = Gs_3
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- Gs_3: 19 components, all dim 0, radical.
  --  total # points = 30.
  --  exactly one is stable.  The last ideal consists of 10 points, 6 of which are real.
  --  one can look at stability at these real points.
  --  one cal=n also eyeball that all such points are unstable:
  --   x3 >= 0
  --   x2 + x3 + x4 >= 0
  --   -x2 - 3x3 - x4 >= 0
  --   therefore: -x3 >= 0
  --   therefore x3 == 0. therefore x2+x4 == 0
  --   but C_-1 + ideal(x3, x2+x4) == 1, so no stable points.

  -- Gs_4
  G = Gs_4
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- radical, one comp: dim 2.  rest: 12 components, dim 0 degree == 16 total.
  -- the comp of dim 2: has only UNSTABLE points.
  -- (because: look at sum of 2 princ minors of -J)

  -- Gs_5
  G = Gs_5
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  C/degree
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim=0, radical, 19 components, 36 points.
  -- one stable, all others UNSTABLE
  -- there are 2 components which one needs to check:
  -- both consist of unstable points (using princ minors, see my notes).     

  -- Gs_6
  G = Gs_6
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  C/degree
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim=1, not radical, 3 components of dim 1, 10 of dim 0.
  -- all points are unstable, except SS.


  -- Gs_7
  G = Gs_7
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  C/degree
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim 1, not radical
  -- 3 comps of dim 1.
  -- 12 of dim 0 (16 points here)
  -- exactly one stable point, all others UNSTABLE.

  -- Gs_8
  G = Gs_8
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  C/degree
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim 0, radical, 19 components, 22 points.
  -- princ minors all simple, exacxtly one stable point.  

  -- Gs_9
  G = Gs_9
  S = oscRing(G, CoefficientRing => QQ)
  I = oscSystem(G,S);
  netList I_*
  dim I
  I == radical I
  C = decompose I;
  C = C/trim;
  netList C
  C/dim
  C/degree
  J = oscJacobian(G,S)
  netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
  -- dim 1, radical
  -- 2 comps of dim 1.
  -- 16 comps of dim 0, (16 smooth points).
  -- exactly one point stable.
///

TEST ///
-*
  restart
  needsPackage "Oscillators"
  
  needsPackage "Visualize" 
  openPort "8080"
  visualize Gs_3 -- 
  closePort()
*-    
  needsPackage "NautyGraphs"
  Gstrs = flatten for i from 4 to 9 list generateGraphs(5,i, OnlyConnected=>true)
  Gs = Gstrs/stringToGraph
  assert(#Gs == 20)
  Gs = select(Gs, G -> all(vertices G, v -> degree(G,v) > 1))
  assert(#Gs == 10)  
  
  -- now let's compute these numerically:
  -- note: this may completely miss the points on
  --  higher dimensional components.  I'm not sure,
  --  but for singular isolated points, I think it
  --  will find all of those, maybe with higher multiplicity.
  n = 5
  RC = oscRing(n-1,{}, CoefficientRing => CC)
  R = oscRing(n-1,{}, CoefficientRing => QQ)
  realsolsGs = for G in Gs list (
      IC = oscSystem(G,RC);
      JC = oscJacobian(G, RC);
      elapsedTime realsols = findRealSolutions IC;
      realsols = realsols/(pt -> prepend(identifyStability(JC, pt), pt));
      stablesols = select(realsols, x -> x#0 === Stable);
      (stablesols, realsols)
      )
///

TEST ///
  restart
  needsPackage "Oscillators"
  needsPackage "NautyGraphs"
  n = 5;
  R = oscRing (n, Start => 0);

  Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true);
  #Gstrs
  #Gstrs == 20
  G5s = Gstrs/stringToGraph;
  G5s = select(G5s, hasNoLeaf); #G5s 
  assert(#G5s == 10)

  I0 = oscQuadrics(G5s_0, R)
  ans0 = ideal(x_2*y_0+x_3*y_0-x_0*y_2-x_0*y_3,x_3*y_1+x_4*y_1-x_1*y_3-x_1*y_4,-x_2*y_0+x_0*y_2+x_4*y_2-x_2*y_4,-x_3*y_0-x_3*y_1+x_0*y_3+x_1*y_3,-x_4*y_1-x_4*y_2+x_1*y_4+x_2*y_4)
  assert(gens I0 ==  gens ans0)

  A = diff(transpose matrix{{y_0,y_1,y_2,y_3,y_4}}, gens I0)
  det A
  det A_{0,1,2,3}^{0,1,2,3}
  
  I0 = oscQuadrics(G5s_4, R)
  M = genericMatrix(R, R_0, n, 2)
  Seg = minors(2, M)
  I1 = I0 : Seg
  trim(I1 + Seg)
  (codim oo, degree oo)
  res coker A
  
  for G in G5s list G => getLinearlyStableSolutions G
  pt = first getLinearlyStableSolutions G5s_0
  pt1 = {1.0} | pt_{0..3} | {0.0} | pt_{4..7}
  RC = CC (monoid R)
  I = oscQuadrics(G5s_0, RC)
  sub(I, matrix {pt1})
  comps = decompose I
  netList comps
  comps/(i -> sub(i, matrix{pt1}))  
  sub(comps_1, x_0 => 1, y_0 => 0)
  dim comps_1
  oscSystem(G5s_0, R)
  
  for g in G5s list (
      I := oscQuadrics(g, R);
      radical I == I
      )

  for g in G5s list (
      I := oscQuadrics(g, R);
      dim I
      )
  for g in G5s list hasNoLeaf g
  decompose(oscQuadrics(G5s_1, R))

  for g in G5s list (
      I := trim(oscQuadrics(g, R) + ideal(x_0-1,y_0));
      radical I == I
      )
  for g in G5s list (
      I := trim(oscQuadrics(g, R) + ideal(x_0-1,y_0));
      decompose I
      )
  netList oo
  for g in G5s list (
      I := trim(oscQuadrics(g, R) + ideal(x_0-1,y_0));
      dim I
      )
  for g in G5s list (
      I := ideal gens gb trim(oscQuadrics(g, R) + ideal(x_0-1,y_0));
      I
      )
  
  G5s_4 -- square with extra edge to 5th vertex. 
  I = oscQuadrics(G5s_4, R) 
  decompose I 
  intersect oo == I
  (decompose I)/codim
  codim I  
  mingens I
  Itrig = ideal apply(n, i->x_i^2+y_i^2-1);  
  J = I + Itrig    
  for c in decompose I list decompose (c+Itrig)
  
  Itrig = ideal apply(n, i->x_i^2+y_i^2-1);  
  for g in G5s list (
      I := oscQuadrics(g, R);
      J := trim(I + Itrig + ideal(y_0, x_0-1));
      dim J
      )
  for g in G5s list (
      I := oscQuadrics(g, R);
      J := trim(I + Itrig + ideal(y_0, x_0-1));
      ans := (J == radical J);
      print ans;
      ans
      )

  I = oscQuadrics(G5s_1, R)
  codim I
  decompose I
  Itrig = ideal apply(n, i->x_i^2+y_i^2-1);  
  for c in decompose I list dim trim(c + Itrig)
  (decompose I)_1 + Itrig
  for c in decompose I list decompose trim(c + Itrig)


  I = (G5s_2, R)
  compsI = decompose I
  for c in compsI list gens gb trim(c + Itrig)
  for c in compsI list dim trim(c + Itrig)
  for c in compsI list degree trim(c + Itrig)
  for c in compsI list ((decompose trim(c + Itrig + ideal(y_0, x_0-1)))/degree)
  netList oo

  J = trim(I + Itrig + ideal(x_0-1,y_0))
    
  R1 = oscRing(4, {}, CoefficientRing => ZZ/32003)
  for g in G5s list (
      I := oscSystem(g, R1);
      elapsedTime radical I == I
      )

  id1 = oscSystem(G5s_6, R1)
  rad1 = radical id1
  decompose id1
  decompose rad1
  isSubset(id1, rad1)
  rad1 : id1
  id1 : rad1
  intersect decompose id1 == id1
restart
CkSegreTrig = (n)->(R=ZZ/32003[x_0..x_(n-1),y_0..y_(n-1)];
                                Itrig = ideal apply(n, i->x_i^2+y_i^2-1);
                                SegreI =minors(2,genericMatrix(R,n,2));
                                I=Itrig+SegreI;
                                radical(I)==I)

I1 = ideal I_*
intersect decompose I1 == I1
CkSegreTrig 5        
///


TEST ///
  --Let's do all connected graphs with 4 vertices, which do not have a vertex with degree 1
  restart
    needsPackage "Oscillators"
    needsPackage "NautyGraphs"
    needsPackage "Visualize"

    n = 4
    -- There are 6 connected graphs:
    Gstrs = flatten for i from 3 to 6 list generateGraphs(n,i, OnlyConnected=>true)
    Gs = Gstrs/stringToGraph
    assert(#Gs == 6)
    -- the first 2 are trees:
    assert(Gs/isTree == splice{2:true, 4:false})
    Gs = select(Gs, G -> all(vertices G, i -> degree(G,i) > 1))
    assert(#Gs == 3)
    -- there are (at least) two ways to view graphs in Macaulay2.  The first 
    -- just pops up a jpeg picture using graphviz (which you must have
    -- installed on your system):
    displayGraph Gs_0
    Gs/displayGraph

    -- The next lines are for viewing graphs in your browser.
    -- This generally gives a nicer picture than graphviz, and you can change the look,
    -- and output tikz code too.  However, you need to end the session in the browser
    -- in order to move on to the next graph.
    -- SS: self synchro (only one stable real solution)
    -- POS: has positive dimensional components
    -- NONRAD: is not reduced
    -- NONREG: some solutions are not regular (from homotopy continuation)
    openPort "8080"
    visualize Gs_0 -- square NONREG SS POS
    visualize Gs_1 -- square with one diagonal NONREG SS NONRAD POS
    visualize Gs_2 -- complete graph NONREG SS NONRAD POS
    closePort()

    RC = oscRing(n-1,{}, CoefficientRing => CC)
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    realsolsGs = for G in Gs list (
        IC = oscSystem(G,RC);
        JC = oscJacobian(G, RC);
        elapsedTime realsols = findRealSolutions IC;
        realsols = realsols/(pt -> prepend(identifyStability(JC, pt), pt));
        stablesols = select(realsols, x -> x#0 === Stable);
        (stablesols, realsols)
        )

    netList realsolsGs_0_1 -- one stable, 4 semistable (Jacobian == 0 at these), all rest unstable
    netList realsolsGs_1_1 -- 6 standard solutions, 6 others.
    netList realsolsGs_2_1 -- only standard solutions, yet some of the unstable ones must be singular
                             -- since some solutions appear more than once.

    realsolsGsQQ = for G in Gs list (
        I = oscSystem(G,R);
        JC = oscJacobian(G, RC);
        J = oscJacobian(G,R);
        C = decompose I;
        C = select(C, i -> dim i == 0);
        C = C/(i -> sub(i, RC));
        realsols = C/findRealSolutions//flatten;
        stablesols = select(realsols, p -> isStableSolution(JC, p));
        {stablesols, realsols}
        )

    netList realsolsGsQQ_0_1 -- this one has 2 isolated solutions, both with sin theta_i = 0. (all i)
    netList realsolsGsQQ_1_1 -- this one has 4 isolated solutions, all with sin theta_i = 0
    netList realsolsGsQQ_2_1 -- this one has 5 isolated solutions, all with sin theta_i = 0

    Is = for G in Gs list (
        I = oscSystem(G,R);
        J = oscJacobian(G,R);
        C = decompose I;
        (I, C, J)
        )

    Is/first/dim
    Cs = Is/first/primaryDecomposition
    select(Cs_0, i -> dim i == 0 and i == radical i)
    select(Cs_0, i -> dim i == 0 and i != radical i)

    select(Cs_1, i -> dim i == 0 and i == radical i)
    select(Cs_1, i -> dim i == 0 and i != radical i)

    select(Cs_2, i -> dim i == 0 and i == radical i)
    select(Cs_2, i -> dim i == 0 and i != radical i)

    Is/(x -> x#1)/length

    -- let's consider the higher dimensional components, and radical-ness of ideals
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    for i in Is/first list  radical i == i
    Is/first/dim -- {0,0,0,1,1,1}
    -- Notes:
    -- (1) the first 3 graphs have edges sticking out.
    -- (2) the first 3 graphs have radical ideals, of dim 0.
    -- We now describe the features of each example.
    --
    ---------------------
    -- Graph Gs_0
    -- square
    G = graph ({0, 1, 2, 3}, {{0, 1}, {1,2}, {2, 3}, {0,3}})
    I = oscSystem(G,R)
    C = decompose I
    trim first Is_0 -- 8 solutions: all sin=0, cos=+-1.
    assert(I == radical I) -- true
    assert(dim I == 0) -- dim 0.
    C = decompose I
    J = oscJacobian(G,R)
    Jacs = for i in C list J % i
    for j in Jacs list eigenvalues(lift(j, QQ))
      -- upshot: this is a tree, ideal is a reduced set of the 8 "standard" points.
      --   self-synch pt: yes
      --   other 7 are unstable.
    ---------------------
    -- Graph Gs_1
    -- square with one diagonal
    -- sort vertices Gs_1, (edges Gs_1)/toList/sort//sort
    G = graph ({0, 1, 2, 3}, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}})
    I = oscSystem(G,R);
    C = decompose I
    PD = primaryDecomposition I
    gbI = ideal gens gb trim I
    J = oscJacobian(G,R)
    Jacs = for i in C list J % i
    for j in Jacs list eigenvalues(lift(j, QQ))
      -- upshot: this is a tree, ideal is a reduced set of the 8 "standard" points.
      --   self-synch pt: yes
      --   other 7 are unstable.
    ---------------------
    -- Graph Gs_2
    -- complete graph on 4 vertices
    -- sort vertices Gs_2, (edges Gs_2)/toList/sort//sort
    G = graph ({0, 1, 2, 3}, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}})
    I = oscSystem(G,R);
    C = decompose I
    C/dim
    #C
    PD = primaryDecomposition I
    #PD
    PD/dim
    PD/(i -> i == radical i)
    J = oscJacobian(G,R)
    assert(#C == 8) -- 3 of dim 1, 5 standard points
    for i in C list allUniquePrincipalMinors(-J, Modulo => i)
    C0 = select(C, i -> dim i == 0)
    C1 = select(C, i -> dim i >  0)
    Jacs0 = for i in C0 list J % i
    for j in Jacs0 list eigenvalues(lift(j, QQ)) -- one stable, the other 4 are unstable.
    trim(C1_0 + C1_1) -- a standard point
    trim(C1_0 + C1_2) -- a different standard point
    trim(C1_1 + C1_2) -- a different standard point
    assert(1 == trim(C1_0 + C1_1 + C1_2)) -- don't intersect
    -- geometry: there are 8 components, all reduced.
    --  have 3 curve solutions C0, C1, C2: each consisting solely of unstable points
    --  have 5 standard isolated solutions, one of which is singular.
    --  4 are unstable, 1 is stable.
    -- the pairwise intersection of 2 of the curves is a single (unstable standard point) (all different).
    -- the three curves do not have a common intersection.
    -- there are 3 embedded points (the intersection points of the curves), all standard.
    -- all 8 standard points occur as associated primes.
///

TEST ///
-- Test: oscJacobian and reduced version return the same eigenvalues, except for
--  one extra eigenvalue of 0.0 for the non-reduced form.
-- XXX
  restart
    needsPackage "Oscillators"
    needsPackage "NautyGraphs"

    n = 4
    -- There are 6 connected graphs:
    Gstrs = flatten for i from 3 to 6 list generateGraphs(n,i, OnlyConnected=>true)
    Gs = Gstrs/stringToGraph

    G = Gs_4
    RC = oscRing(3, {}, CoefficientRing => CC)
    J1 = oscJacobian(G,RC,Reduced=>false)
    J2 = oscJacobian(G,RC)

    PC = oscSystem(G,RC);
    IC = PC + trig RC;
    elapsedTime realsols = findRealSolutions IC;
    netList realsols    
    evs = for p in realsols list map(CC, RC, p)
    evsJ1 = for phi in evs list eigenvalues phi J1 -- I think J1 is not correct..
    evsJ2 = for phi in evs list eigenvalues phi J2
///

TEST ///
  --Let's do all connected graphs with 5 vertices -- YYY
  restart
    needsPackage "Oscillators"
    needsPackage "NautyGraphs"
    needsPackage "Visualize"

    n = 5
    Gstrs = flatten for i from 4 to 10 list generateGraphs(n,i, OnlyConnected=>true)
    Gs = Gstrs/stringToGraph
    assert(#Gs == 21)
    openPort "8080"
    visualize Gs_3
    closePort()

    RC = oscRing(n-1,{}, CoefficientRing => CC)
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    realsolsGs = for G in Gs list (
        PC = oscSystem(G,RC);
        JC = oscJacobian(G, RC);
        IC = PC + trig RC;
        elapsedTime realsols = findRealSolutions IC;
        stablesols = select(realsols, p -> isStableSolution(JC,p));
        (stablesols, realsols)
        )
    netList realsolsGs
    realsolsGs/(x -> (#x#0, #x#1))//netList 

    realsolsGsQQ = for G in drop(Gs,-1) list (
        << "doing " << G << endl;
        P = oscSystem(G,R);
        JC = oscJacobian(G, RC);
        J = oscJacobian(G,R);
        I = P + trig R;
        C = elapsedTime decompose I;
        C = select(C, i -> dim i == 0);
        C = C/(i -> sub(i, RC));
        elapsedTime realsols = C/findRealSolutions//flatten;
        stablesols = select(realsols, p -> isStableSolution(JC, p));
        << "stablesols = " << stablesols << endl << " #real sols = " << #realsols << endl;
        {stablesols, realsols}
        )
    netList realsolsGsQQ
    realsolsGsQQ/(x -> (#x#0, #x#1))//netList 

    -- the only one with more than one stable (isolated) solution is the pentagon ring.
    G = Gs_7
    openPort "8081"
    visualize G
    closePort()
    
    dims = for G in drop(Gs,-1) list (
        << "doing " << G << endl;
        P = oscSystem(G,R);
        I = P + trig R;
        C = elapsedTime decompose I;
        G => tally (C/dim)
        )
    positions(dims, x -> x#1#?1 or x#1#?2)
    openPort "8081"
    visualize dims#4#0
    visualize dims#8#0
    visualize dims#9#0
    visualize dims#12#0
    visualize dims#13#0
    visualize dims#15#0
    visualize dims#16#0
    visualize dims#17#0
    visualize dims#19#0
    
    visualize dims#5#0
    visualize dims#6#0    
    visualize dims#7#0    
    visualize dims#10#0    
    visualize dims#11#0    
    visualize dims#14#0    
    closePort()
    
///

TEST ///
-*
  restart
    needsPackage "Oscillators"
*-
    needsPackage "NautyGraphs"
    needsPackage "Visualize"

    n = 6
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true)
    Gcomplete = stringToGraph first generateGraphs(n,binomial(n,2), OnlyConnected=>true)
    assert(#Gstrs == 111)
    Gs = Gstrs/stringToGraph;
    
    G = stringToGraph Gstrs_50
    G = stringToGraph Gstrs_51
    G = stringToGraph Gstrs_30 -- has non-regular solutions.
    openPort "8081"
    visualize G
    closePort()

    getSols = (G) -> (
        << "---- doing graph " << G << endl;
        n := # vertices G;
        RC = oscRing(n-1,{}, CoefficientRing => CC);
        PC = oscSystem(G,RC);
        JC = oscJacobian(G, RC);
        IC = PC + trig RC;
        elapsedTime realsols = findRealSolutions IC;
        stablesols = select(realsols, p -> isStableSolution(JC,p));
        << netList stablesols << endl;
        (stablesols, realsols)
        )
    getSols Gs_0 -- one sol
    getSols Gs_1
    for i from 2 to 10 list
          getSols Gs_i
    for i from 11 to 20 list
          getSols Gs_i
    -- 21 takes too long?
    for i from 22 to 30 list
          getSols Gs_i
    for i from 31 to 40 list
          getSols Gs_i
    for i from 41 to 50 list
          getSols Gs_i
    for i from 51 to 60 list
          getSols Gs_i
    for i from 61 to 70 list
          getSols Gs_i
    for i from 71 to 80 list
          getSols Gs_i
    for i from 81 to 90 list
          getSols Gs_i
    for i from 91 to 100 list
          getSols Gs_i
    for i from 91 to 100 list
          getSols Gs_i
    for i from 101 to 110 list
          getSols Gs_i
    getSols Gs_21 -- ok
    openPort "8081"
    visualize Gs_12 -- 3 stable solutions
    visualize Gs_18 -- 3 stable solutions
    visualize Gs_33 -- 3 stable solutions
    closePort()
 
    getSols Gs_33
    
    G = Gs_21
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    P = oscSystem(G,R);
    J = oscJacobian(G,R);
    I = P + trig R;
    dim I
    degree I
    C0 = select(decompose I, i -> dim i == 0)
    C = primaryDecomposition I
    for i in C0 list lift(J % i, QQ)
    for i in C0 list lift((vars ring i)%i, QQ)
    oo/eigenvalues
    Isat = saturate(I, product(5, i -> y_(i+1)))
    I0 = I : Isat
    I == intersect(I0, Isat)
    solveSystem Isat_*
    ideal gens gb I    
    C = elapsedTime decompose I;
    C = minprimes(I, Verbosity=>2);
    C0 = select(C, i -> dim i == 0);
    
    Rlex = newRing(R, MonomialOrder=>Lex)
    Ilex = sub(I, Rlex)
    gens gb Ilex;
    
    use ring I
    I1 = trim saturate(I, x_2+1)
    I2 = trim eliminate({x_2,y_2,y_4}, I1)
    J1 = I : I1
    I : J1 == I1 -- true
    degree I1
///

TEST ///
-- MES XXX
-*
    restart
    needsPackage "Oscillators"
*-
    needsPackage "NautyGraphs"
    needsPackage "Visualize"
    n = 7
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true)
    Gs = Gstrs/stringToGraph;
    #Gs
    assert(#Gs == 852)
    positions(Gs, isTree) -- the first 11 are trees, that is it...
    
    getSols Gs_11

    -- example: Gs_11 (one triangle, with edges sticking out)
    G = Gs_11
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    P = oscSystem(G,R);
    J = oscJacobian(G,R);
    I = P + trig R;
    dim I
    C = decompose I
    for i in C list ((vars R) % i)
    I == intersect C -- true, so reduced, zero dimensional
    realroots = flatten for i in C list findRealSolutions i -- 96 real points.  Interesting, since degree is 96...
    RC = oscRing(n-1,{}, CoefficientRing => CC)
    JC = oscJacobian(G,RC);
    realroots/(p -> isStableSolution(JC,p))
    -- upshot: 96 roots, all real, exactly one is stable.
    realroots/(p -> (map(CC,RC,p)) JC)
    oo/eigenvalues

    -- example: Gs_500 (planar symmetric triangulation of a triangle)
    G = Gs_500
    R = oscRing(n-1,{}, CoefficientRing => ZZ/32003)
    P = oscSystem(G,R);
    J = oscJacobian(G,R);
    I = P + trig R; -- appears to have dim 0, degree 232 (these are OK in char 32003).
    getSols G -- only 
    
    -- let's do one with each number of edges?
    G = stringToGraph first generateGraphs(n,13, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,14, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,15, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,16, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,17, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,18, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,19, OnlyConnected=>true)    
    getSols G -- SS

    G = stringToGraph first generateGraphs(n,20, OnlyConnected=>true)    
    getSols G -- SS

    #Gs

    -- The main line:
    results7 = for i from 0 to #Gs - 1 list (
        << "----doing graph " << i << " -----------------" << endl;
        getSols Gs_i
        )
    positions(results7, a -> # first a > 1)
    -- {19, 26, 29, 30, 31, 39, 43, 72, 81, 82, 88, 98, 99, 101, 104, 147, 187, 215, 313}
    for a in results7 list (#a#0, #a#1)
    tally oo
    
    dim I
    C = decompose I
    for i in C list ((vars R) % i)
    I == intersect C -- true, so reduced, zero dimensional
    realroots = flatten for i in C list findRealSolutions i -- 96 real points.  Interesting, since degree is 96...
    RC = oscRing(n-1,{}, CoefficientRing => CC)
    JC = oscJacobian(G,RC);
    realroots/(p -> isStableSolution(JC,p))
    -- upshot: 96 roots, all real, exactly one is stable.
    realroots/(p -> (map(CC,RC,p)) JC)
    oo/eigenvalues
///

///
-*
    restart
    needsPackage "Oscillators"
*-
    needsPackage "NautyGraphs"
    needsPackage "Visualize"
    n=7
    G = graph ({0, 4, 1, 5, 2, 6, 3}, {{4, 0}, {0, 5}, {4, 1}, {4, 6}, {5, 2}, {5, 6}, {6, 3}})
    
    R = oscRing(n-1,{}, CoefficientRing => QQ)
    P = oscSystem(G,R);
    J = oscJacobian(G,R);
    I = P + trig R;
    C = decompose I;
    C/dim -- {1,1,1,0,0}
    C/degree -- {2,2,4,1,1}
    C/(i -> i == radical i)

    #C
    eliminate({x_1,y_1,x_2,y_2,x_3,y_3},C_10)
    positions(C, i -> (y_3+y_6) % i == 0)
    positions(C, i -> (y_3-y_6) % i == 0)
    
    getSols = (G) -> (
        << "---- doing graph " << G << endl;
        n := # vertices G;
        RC = oscRing(n-1,{}, CoefficientRing => CC);
        PC = oscSystem(G,RC);
        JC = oscJacobian(G, RC);
        IC = PC + trig RC;
        elapsedTime realsols = findRealSolutions IC;
        stablesols = select(realsols, p -> isStableSolution(JC,p));
        << netList stablesols << endl;
        (stablesols, realsols)
        )
    getSols G
    
    use ring I
    (gens I) % ideal(x_1, x_2, x_3+1, x_4, x_5, x_6 + 1, y_1+1, y_2-1, y_3, y_4+1, y_5-1,y_6)
    J % ideal(x_1, x_2, x_3+1, x_4, x_5, x_6 + 1, y_1+1, y_2-1, y_3, y_4+1, y_5-1,y_6)
///

TEST ///
    n = 8
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true)
    #Gstrs
    assert(#Gstrs == 11116)
///

TEST ///
    n = 9
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true);
    #Gstrs
    assert(#Gstrs == 261079)
///

TEST ///
    n = 10
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true);
    #Gstrs
    assert(#Gstrs == 11716570)
///

TEST ///
    -- this one might take a while...
    -- maybe we need to save these to disk, not try to load them in.
    -- maybe we need to remove graphs without a cycle.
    -- this one basically made me need to reboot my computer! (I think...)
    n = 11
    Gstrs = flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true);
    #Gstrs
    assert(#Gstrs == 1006700564)
///