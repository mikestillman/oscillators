-*
doc ///
  Key
  Headline
  Usage
  Inputs
  Outputs
  Description
    Text
    Example
  Caveat
  SeeAlso
///
*-

undocumented {
  (oscRing, Graph, List),
  Symbols,
  trig,
  (trig, Ring),
  Stable,
  Semistable,
  Unstable,
  Modulo,
  Radians,
  numOscillators
}

doc ///
Key
  Oscillators
Headline
  generation and analysis of oscillator steady states for small graphs
Description
  Text
    @SUBSECTION "Computations from the paper [Harrington, Schenck, Stillman, Algebraic aspects of homogeneous Kuramoto oscillators]"@
  Text
    @UL {
        TO "Generation of all SCT (simple, connected, 2-connected) graphs on small numbers of vertices",
        TO "B",
        TO "C",
        TO "D",
        TO "E",
        TO "F1"
        }@
  Text
    This package supports computations with Kuramoto oscillators, including computations for the paper [HSS], Harrington, Schenck, Stillman, @arXiv("2312.16069", "Algebraic aspects of homogeneous Kuramoto oscillators")@.
  Text
    We show a possible workflow using this package.  We use @TO "NautyGraphs::NautyGraphs"@
    to generate graphs of small size.  We use @TO "Visualize::Visualize"@ to look at these graphs.
  Example
    needsPackage "Oscillators"
    needsPackage "NautyGraphs"
    needsPackage "Visualize"
    vert4edges4 = generateGraphs(4, 4)
    Gs = vert4edges4/stringToGraph
  Text
    Now, for each graph we consider the oscillator dynamical system.
    
    Let's do an example, using a graph with 5 vertices and 7 edges
  Example
    Gstrs = generateGraphs(5, OnlyConnected => true, MinDegree => 2);
    Gs = Gstrs/stringToGraph
    Gcycle5 = Gs_3
    K5 = Gs_-1
  Pre
    openPort "8083"
    visualize G -- important: click on End session in browser window before continuing
    closePort()
  Example
    G = Gcycle5
    R = oscRing(G, CoefficientRing => CC, Reduced => true)
    I = oscSystem(G, R);
    netList I_*
    Jac = oscJacobian(G,R)
    realsols = findRealSolutions I;
    netList realsols
    assert(# realsols == 30)
    stablesols = select(realsols, p -> Stable === identifyStability(Jac,p))
    tally realsols/(p -> identifyStability(Jac, p))
    getAngles(4, stablesols, Radians=>false)
    sols = solveSystem I_*;
    sols = sols/coordinates;
    sols = matrix sols -- every row is a solution
    sols = clean(1e-6, sols) -- set to 0 numbers very close to 0
    assert(numrows sols == 30)
    tally (entries sols)/(p -> identifyStability(Jac, p))
    sols
///

///
Key
  Oscillators
Headline
  generation and analysis of oscillator steady states for small graphs
Description
  Text
    We show a possible workflow using this package.  We use @TO "NautyGraphs::NautyGraphs"@
    to generate graphs of small size.  We use @TO "Visualize::Visualize"@ to look at these graphs.
  Example
    needsPackage "Oscillators"
    needsPackage "NautyGraphs"
    needsPackage "Visualize"
    vert4edges4 = generateGraphs(4, 4)
    Gs = vert4edges4/stringToGraph
  Text
    A slightly larger example.
  Example
    vert7edges10 = generateGraphs(7, 10);
    Gs = vert7edges10/stringToGraph
  Text
    Now, for each graph we consider the oscillator dynamical system.
    
    Let's do an example, using a graph with 5 vertices and 7 edges
  Example
    Gstrs = generateGraphs(5,7);
    Gs = Gstrs/stringToGraph
    G = Gs_1
  Pre
    openPort "8083"
    visualize G -- important: click on End session in browser window before continuing
    closePort()
  Example
    R = oscRing(G, CoefficientRing => CC, Reduced => true)
    I = oscSystem(G,R);
    Jac = oscJacobian(G,R)
    realsols = findRealSolutions I
    assert(# realsols == 24)
    stablesols = select(realsols, p -> Stable === identifyStability(Jac,p))
    tally realsols/(p -> identifyStability(Jac, p))
    getAngles(4, stablesols, Radians=>false)
    sols = solveSystem I_*;
    sols = sols/coordinates
    assert(#sols == 36)
    tally sols/(p -> identifyStability(Jac, p))
    netList sols
///

doc ///
  Key
    oscRing
    (oscRing, Graph)
    (oscRing, ZZ)
    [oscRing, Reduced]
    [oscRing, CoefficientRing]
    [oscRing, Symbols]
  Headline
    create a polynomial ring for a given graph or number of oscillators
  Usage
    S = oscRing G
    S = oscRing n
  Inputs
    n:ZZ
      The number of oscillators
    G:Graph
      The number of oscillators is the number of vertices in the Graph
    CoefficientRing => Ring
      the coefficient ring to use, for numerical work, @TO "CC"@ is a good choice
    Symbols => Sequence
      a sequence of two symbols.  The first refers to the cosine of the given angles
        and the second refers to the sine of the given angles
    Reduced => Boolean
      if true, then the angles are reduced to the first angle being 0. That is, the number of oscillators will be one less than the number of vertices of $G$
  Outputs
    :Ring
  Description
    Text
      This function returns a polynomial ring in $2n$ variables
      representing the cosine and the sine of $n$ angles. Inputting a number assumes that the vertices are labeled $0, \ldots, n$.
    Example
      oscRing(3, CoefficientRing => CC)
      oscRing(3, CoefficientRing => CC, Symbols=>{c, s})
    Text
      If one gives a graph with Reduced => true, then the number of angles considered
      is one less than the number of vertices of the graph.
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing G
      oscRing(G, CoefficientRing => CC, Reduced => true)
  SeeAlso
    oscSystem
    oscJacobian
    displayGraph
///

doc ///
  Key
    oscSystem
    (oscSystem, Graph, Ring)
    (oscSystem, Graph)
    [oscSystem, Reduced]
  Headline
    the ideal of the reduced equilibrium points of a dynamical system of oscillators
  Usage
    oscSystem(G,R)
    oscSystem(G)
  Inputs
    G:Graph
      an undirected, connected graph $G$
    R:Ring
      created with @TO oscRing@
    Reduced => Boolean
      if true, then the angles are reduced to the first angle being 0. That is, the number of oscillators will be one less than the number of vertices of $G$
  Outputs
    :Ideal
      in the ring R
  Description
    Text
      $R$ should be a ring created with @TO "oscRing"@.  The dynamical system
      involved is the oscillator system associated to $G$: one angle per vertex.
      If $a_{ij} = 1$ if $(i,j)$ is an edge of the undirected graph $G$, 
      and is zero otherwise, then the system is $d\theta_i/dt = \sum_j a_{ij} \sin(\theta_j - \theta_i)$
      where we consider only reduced equilibrium solutions $\theta_0 = 0$.
      
      This function returns the ideal of equilibrium points, where angles $(0, \theta_1, ..., \theta_{n-1})$
      are represented via cosines and sines of the angles.
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing(G, CoefficientRing => CC)
      R = oo
      I = oscSystem(G,R)
      netList I_*
    Text
      We can find approximations to the 26 complex solutions to this system.
      If the system has positive dimension (not the case here), the idea is that this set of points should
      contain at least one on each component.
    Example
      solveSystem I_*
      #oo
    Text
      We can find approximations to the 6 real solutions to this system.  
    Example
      findRealSolutions I
      #oo
    Text
      The angles of these solutions (in degrees, not radians, and the 3 refers to the
      numbner of oscillators).
    Example
      netList getAngles(3, findRealSolutions I, Radians=>false)
  SeeAlso
      oscRing
      oscJacobian
      findRealSolutions
      getAngles
      solveSystem
///

doc ///
  Key
    oscJacobian
    (oscJacobian,Graph,Ring)
    (oscJacobian, Graph)
    (oscJacobian, Ideal)
    [oscJacobian, Reduced]
  Headline
    create the Jacobian for the oscillator system associated to a graph
  Usage
    oscJacobian(G,R)
    oscJacobian(I)
    oscJacobian(G)
  Inputs
    G:Graph
      an undirected, connected graph $G$ 
    S:Ring
      created with @TO oscRing@
    I:Ideal
      an ideal, intended to be ideal created by oscSystem
    Reduced => Boolean
      if true, then the angles are reduced to the first angle being 0. That is, the number of oscillators will be one less than the number of vertices of $G$
  Outputs
    :Matrix
      the $n \times n$ Jacobian matrix of the system, as a matrix of polynomials
      involving the cosines and sines of angles $\theta_1, \ldots, \theta_{n-1}$, where we
      set $\theta_0 = 0$ if $\texttt{Reduced => true}$
  Description
    Text
      The matrix is a symmetric $n \times n$  matrix (with determinant zero).
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing(G, CoefficientRing => CC)
      S = oo
      I = oscSystem(G,S)
      Jac = oscJacobian(G,S)
      assert(det Jac == 0)
      assert(Jac - transpose Jac == 0)
    Text
      We can find the eigenvalues of the Jacobian at approximate points, and see if they are
      stable (all eigenvalues negative, except for the one required 0), unstable (a positive 
      eigenvalue), or semistable (no positive eigenvalues, up to a certain tolerance).
    Example
      realsols = findRealSolutions I
      jacs = for pt in realsols list sub(Jac, matrix{pt})
      jacs/eigenvalues
      jacs/eigenvalues/identifyStability
  SeeAlso
    oscRing
    oscSystem
    findRealSolutions
    identifyStability
    eigenvalues
///

doc ///
  Key
    findRealSolutions
    (findRealSolutions, Ideal)
    (findRealSolutions, Graph)
  Headline
    find real solutions, at least one per component for well-conditioned systems
  Usage
    findRealSolutions I
    findRealSolutions G
  Inputs
    I:Ideal
      an ideal in a polynomial ring over QQ or RR or CC.
    G:Graph
      an undirected, connected graph $G$
  Outputs
    :List
      of all the real solutions that were found
  Description
    Text
      In this package the main use is to zero in on equilibrium points for oscillator 
      systems associated to graphs.
      
      We use this in the following example.
    Example
      G = graph({0,1,2,3,4}, {{0,1},{1,2},{2,3},{0,3},{0,4}})
      S = oscRing(G, CoefficientRing => QQ, Reduced => true)
      I = oscSystem(G,S)
      dim I
      Jac = oscJacobian(G,S)
      assert(det Jac == 0)
      assert(Jac - transpose Jac == 0)
      realsols = findRealSolutions I      
      netList getAngles(4, realsols, Radians=>false)
      C = decompose I
      for i in C list Jac % i
      M = Jac % C_1
      M1 = matrix{{-1,0,0,0,1},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,-1}}
      use S
      eigenvalues lift(sub(M - M1, x_3=>1), QQ)
      eigenvalues lift(sub(M, x_3=>1), QQ)
      eigenvalues lift(sub(M, x_3=>1/100000), QQ)
      eigenvalues lift(sub(M, x_3=>100), QQ)
  Caveat
    If the system is positive dimensional, then of course it won't find all roots. There might be
    real roots on a component, but none are found. If a component (even a point) is singular, then
    it might have problems, depending on the situation. In the latter case, there
    is a warning message emitted (about non-regular solutions).
  SeeAlso
    oscRing
    oscSystem
    oscJacobian
    identifyStability
///

doc ///
  Key
    oscQuadrics
    (oscQuadrics, Graph, Ring)
    (oscQuadrics, Graph)
    [oscQuadrics, Reduced]
  Headline
    find the homogeneous quadrics in the homogeneous Kuramoto ideal
  Usage
    oscQuadrics(G,R)
    oscQuadrics(G)
  Inputs
    G:Graph
      an undirected, connected graph $G$ on vertices $0, \ldots, n-1$, where $n$
      is the number of vertices of $G$
    R:Ring
      created with @TO oscRing@
    Reduced => Boolean
      if true, then the angles are reduced to the first angle being 0. That is, the number of oscillators will be one less than the number of vertices of $G$
  Outputs
    :Ideal
      in the ring R
  Description
    Text
      The ideal of quadrics of the oscillator system associated to the graph $G$.
      This is the ideal generated by the quadrics of the system, where the system is
      $d\theta_i/dt = \sum_j a_{ij} \sin(\theta_j - \theta_i)$
      where we consider only reduced equilibrium solutions $\theta_0 = 0$.
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing(G, CoefficientRing => CC)
      S = oo
      I = oscQuadrics(G,S)
      netList I_*
///

doc ///
Key
  standardSols
  (standardSols, Graph, Ring)
  (standardSols, Graph)
  [standardSols, Reduced]
Headline
  find the "standard solutions" for the oscillator system associated to a graph
Usage
  standardSols(G,R)
  standardSols(G)
Inputs
  G:Graph
    an undirected, connected graph $G$
  R:Ring
    created with @TO oscRing@
  Reduced => Boolean
    if true, then the angles are reduced to the first angle being 0. That is, the number of oscillators will be one less than the number of vertices of $G$
Outputs
  :Ideal
    a list of the "standard solutions" that you get where all the angles are the same. This is the same as a certain Segre variety.
Description
  Text
    The oscillator ideal associated to a graph constructed by @TO oscQuadrics@ always contains the "standard solutions" as a minimal prime. These standard solutions are the Segre variety $\mathbb{P}^1\times \mathbb{P}^{n-1}$, where $n$ is the number of vertices of the graph. These are the solutions where all the angles are the same.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    R = oscRing(G);
    I = oscQuadrics(G, R);
    any(decompose I, P -> P == standardSols(G, R))
SeeAlso
  oscRing
  oscSystem
  oscJacobian
///

doc ///
Key 
  vertexSpanningPolynomial
  (vertexSpanningPolynomial, Graph)
  (vertexSpanningPolynomial, Graph, Ring)
Headline
  computes the vertex spanning polynomial
Usage
  vertexSpanningPolynomial(G)
  vertexSpanningPolynomial(G, R)
Inputs
  G:Graph
    an undirected, connected graph $G$
  R:Ring
    created with @TO oscRing@
Outputs
  :PolynomialRing
    the vertex spanning polynomial of the graph $G$ inside $R$.
Description
  Text
    Let $S$ be the set of all spanning trees of the graph $G$. For each spanning tree $T$, let $d_i$ be the degree of $v_i$ in $T$. The vertex spanning polynomial of a graph $G$ is defined as $\sum_{T\in S} \prod_{v_i \in T} x_i^{d_i-1}$. The factorization of this polynomial is related to the number of components of the oscillator ideal of the graph computed via @TO oscQuadrics@.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    vertexSpanningPolynomial G
SeeAlso
  oscQuadrics
///

doc ///
Key
  getAngles
  (getAngles, ZZ, List)
  [getAngles, Radians]
Headline
  Compute angles from a list of solutions
Usage
  getAngles(n, S)
Inputs
  n: ZZ
    An integer representing the number of angles.
  sols: List
    A list of solutions, where each solution is a list of cosine and sine values.
  Radians => Boolean
    An optional boolean flag (default is true). If true, the angles are returned in radians; otherwise, they are returned in degrees.
Outputs
  : List
    A list of angles computed from the given solutions.
Description
  Text
    The function getAngles computes the angles from a list of solutions. Each solution is assumed to be a list of cosine and sine values for the angles. The function returns the angles in radians or degrees based on the Radians option.
  Example
    getAngles(2, {{1, 0, 0, 1}, {0, 1, 1, 0}}, Radians => false)
SeeAlso
  getLinearlyStableSolutions
  findRealSolutions
///

doc ///
Key
  identifyStability
  (identifyStability, Matrix, List)
  (identifyStability, BasicList)
  [identifyStability, Tolerance]
Headline
  Identify the stability of a list of eigenvalues, or of potential solutions to the oscillator system
Usage
  identifyStability(Jac, sols)
  identifyStability(eigenvals)
Inputs
  Jac: Matrix
    A matrix representing the Jacobian of the oscillator system.
  sol: List
    A potential solution to the oscillator system.
  eigenvals: BasicList
    An eigenvalue
  Tolerance => RR
    An optional real number representing the tolerance for the stability check (default is 1e-10).
Outputs
  : BasicList
    A list of stability values for the given solutions or eigenvalues.
Description
  Text
    The function identifyStability computes the stability of a list of solutions or eigenvalues. For solutions, the function computes the eigenvalues of the Jacobian at each solution and determines the stability based on the sign of the real part of the eigenvalues. For eigenvalues, the function determines the stability based on the sign of the real part of the eigenvalues.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    R = oscRing(G, Reduced => true);
    I = oscSystem(G, R);
    Jac = oscJacobian(G, R);
    realsols = findRealSolutions I;
    jacs = for pt in realsols list sub(Jac, matrix{pt});
    eigenvals = jacs/eigenvalues;
    eigenvals / identifyStability
    realsols/(pt -> prepend(identifyStability(Jac, pt), pt))
SeeAlso
  oscJacobian
  findRealSolutions
///

doc ///
Key
  hasNoLeaf
  (hasNoLeaf, Graph)
Headline
  Check if a graph has no leaf vertices
Usage 
  hasNoLeaf(G)  
Inputs
  G: Graph
    An undirected graph
Outputs
  : Boolean
    True if the graph has no leaf vertices, false otherwise
Description
  Text
    The function hasNoLeaf checks if a given graph has no leaf vertices. A leaf vertex is a vertex with degree 1, meaning it is connected to only one other vertex in the graph.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    hasNoLeaf(G)
///

doc ///
Key
  getLinearlyStableSolutions
  (getLinearlyStableSolutions, Graph)
Headline
  Compute linearly stable solutions for the Kuramoto oscillator system associated to a graph
Usage
  getLinearlyStableSolutions(G)
Inputs
  G: Graph
    An undirected, connected graph
Outputs
  : List
    A list of linearly stable solutions for the Kuramoto oscillator system associated to the graph
Description
  Text
    The function getLinearlyStableSolutions computes the linearly stable solutions for the Kuramoto oscillator system associated to a given graph. The Kuramoto oscillator system is a system of coupled phase oscillators, where the dynamics of each oscillator is given by the Kuramoto model. The linear stability of a solution is determined by the eigenvalues of the Jacobian matrix of the system evaluated at the solution.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    getLinearlyStableSolutions(G)
SeeAlso
  findRealSolutions
  identifyStability
///

doc ///
Key
  allUniquePrincipalMinors
  (allUniquePrincipalMinors, Matrix)
  [allUniquePrincipalMinors, Modulo]
Headline
  Compute all unique principal minors of a given matrix
Usage
  allUniquePrincipalMinors(M)
Inputs
  M: Matrix
    A square matrix
  Modulo => Ideal
    An optional ideal to compute the principal minors modulo the given ideal
Outputs 
  : List
    A list of all unique principal minors of the given matrix
Description
  Text
    The function allUniquePrincipalMinors computes all unique principal minors of a given square matrix. A principal minor of a matrix is the determinant of a submatrix obtained by deleting the same set of rows and columns from the matrix. The function returns a list of all unique principal minors of the given matrix.
  Example
    G = graph({0,1,2,3},{{0,1},{1,2},{2,3},{0,3}})
    S = oscRing(G, CoefficientRing => QQ, Reduced => true)
    I = oscSystem(G,S);
    C = decompose I;
    J = oscJacobian(G,S)
    netList for i in C list allUniquePrincipalMinors(-J, Modulo=>i)
    -- by looking at each one, all points but one are unstable.
SeeAlso
  oscJacobian
  identifyStability
///

doc ///
Key
  isStableSolution
  (isStableSolution, Matrix, List)
Headline
  Check if a given solution is stable for the Kuramoto oscillator system
Usage
  isStableSolution(Jac, sol)
Inputs
  Jac: Matrix
    A matrix representing the Jacobian of the Kuramoto oscillator system
  sol: List
    A potential solution to the Kuramoto oscillator system
Outputs
  : Boolean
    True if the solution is stable, false otherwise
Description
  Text
    The function isStableSolution checks if a given solution is stable for the Kuramoto oscillator system. The stability of a solution is determined by the eigenvalues of the Jacobian matrix of the system evaluated at the solution. If all eigenvalues have negative real parts, the solution is considered stable.
  Example
    G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}});
    R = oscRing(G, Reduced => true);
    I = oscSystem(G, R);
    Jac = oscJacobian(G, R);
    realsols = findRealSolutions I;
    select(realsols, S -> isStableSolution(Jac, S))
SeeAlso
  oscJacobian
  findRealSolutions
  identifyStability
///


doc ///
  Key
    "Generation of all SCT (simple, connected, 2-connected) graphs on small numbers of vertices"
  Headline
    generating all SCT graphs on n vertices
  Description
    Text
      Using the NautyGraphs package, we generate the list of
      isomorphism classes of the SCT (simple, connected, 2-connected)
      graphs with a fixed number of vertices.
    Example
      needsPackage "Oscillators"
      needsPackage "NautyGraphs"
      Gstrs = generateGraphs(5, OnlyConnected => true, MinDegree => 2);
      #Gstrs == 11
      Gstrs = generateGraphs(6, OnlyConnected => true, MinDegree => 2);
      #Gstrs == 61
      Gstrs = generateGraphs(7, OnlyConnected => true, MinDegree => 2);
      #Gstrs == 507
      Gstrs = generateGraphs(8, OnlyConnected => true, MinDegree => 2);
      #Gstrs == 7442
      Gstrs = generateGraphs(9, OnlyConnected => true, MinDegree => 2);
      #Gstrs == 197772
    Text
      Here is a simple table with all of these numbers.
    Example
      netList for n from 5 to 9 list {n, #generateGraphs(n, OnlyConnected => true, MinDegree => 2)}
  SeeAlso
///

-- XXX

doc ///
  Key
    "Checking the codimension and irreducible decomposion of the IG ideal"
  Headline
    generating all SCT graphs on n vertices
  Description
    Text
      We first construct the ideal $I_G$ for a specific graph $G$ on 5 vertices.
      We use the 5-cycle as the specific example.
    Example
      needsPackage "Oscillators"
      needsPackage "NautyGraphs"
      Gstrs = generateGraphs(5, OnlyConnected => true, MinDegree => 2);
      Gs = Gstrs/stringToGraph
      G = Gs_3
      R = oscRing(5, Reduced => false);
      IG = oscQuadrics(G, R)
      netList IG_*
      codim IG
    Text
      Each $I_G$ on $n$ vertices has $n-1$ minimal generators.
      This particular $I_G$ has the same codimension $n-1=4$, so is a complete intersection.
      This ideal decomposes as an intersection of 2 prime ideals.
    Example
      comps = decompose IG;
      netList comps_0_*, netList comps_1_*
      comps/codim
      comps/degree
      comps/isPrime
    Text
      Let's see how many graphs are not complete intersections, i.e. have codimension
      $\le n-2$.
    Example
      # select(Gs, G -> (
          IG = oscQuadrics(G, R);
          codim IG <= #vertices G - 2
          ))
    Example
      for G in Gs list (
          IG = oscQuadrics(G, R);
          elapsedTime comps := decompose IG;
          {comps/codim, comps/degree}
          );
      netList oo
    Example
      n = 6
      Gstrs = generateGraphs(n, OnlyConnected => true, MinDegree => 2);
      Gs = Gstrs/stringToGraph;
      R = oscRing(n, Reduced => false);
      # select(Gs, G -> (
          IG = oscQuadrics(G, R);
          codim IG <= #vertices G - 2
          ))
    Example
      allcomps = for G in Gs list (
          IG = oscQuadrics(G, R);
          elapsedTime comps := decompose IG;
          {comps/codim, comps/degree}
          );
      netList ({{"codimensions", "degrees"}} | allcomps)
    Pre
      n = 7
      Gstrs = generateGraphs(n, OnlyConnected => true, MinDegree => 2);
      Gs = Gstrs/stringToGraph;
      R = oscRing(n, Reduced => false);
      # select(Gs, G -> (
          IG = oscQuadrics(G, R);
          elapsedTime codim IG <= #vertices G - 2
          ))
    Text
      The next one takes some time.
    Pre
      n = 8
      Gstrs = generateGraphs(n, OnlyConnected => true, MinDegree => 2);
      Gs = Gstrs/stringToGraph;
      R = oscRing(n, Reduced => false);
      # select(Gs, G -> (
          IG = oscQuadrics(G, R);
          codim IG <= #vertices G - 2
          ))
  
  SeeAlso
///
