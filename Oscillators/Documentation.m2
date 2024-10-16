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

doc ///
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
    R = oscRing(G, CoefficientRing => CC)
    I = oscSystemReduced(G,R);
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
  Headline
    create a polynomial ring for a given graph or number of oscillators
  Usage
    S = oscRing G
    S = oscRing n
  Inputs
    n:ZZ
      The number of oscillators
    G:Graph
      The number of oscillators will be one less than the number of vertices of $G$
    CoefficientRing => Ring
      the coefficient ring to use, for numerical work, @TO "CC"@ is a good choice
    Symbols => Sequence
      a sequence of two symbols.  The first refers to the cosine of the given angles
        and the second refers to the sine of the given angles
  Outputs
    :Ring
  Description
    Text
      This function returns a polynomial ring in $2n$ variables
      representing the cosine and the sine of $n$ angles.
    Example
      oscRing(3, CoefficientRing => CC)
      oscRing(3, CoefficientRing => CC, Start=>0)
      oscRing(3, CoefficientRing => CC, Start=>0, Symbols=>{c, s})
    Text
      If one gives a graph, then the number of angles considered
      is one less than the number of vertices of the graph.
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing G
      oscRing(G, CoefficientRing => CC)
  SeeAlso
    oscSystemReduced
    oscJacobian
    displayGraph
///

doc ///
  Key
    oscSystemReduced
    (oscSystemReduced, Graph, Ring)
  Headline
    the ideal of the reduced equilibrium points of a dynamical system of oscillators
  Usage
    oscSystemReduce
  Inputs
    G:Graph
    S:Ring
  Outputs
    :Ideal
      in the ring S
  Description
    Text
      $S$ should be a ring created with @TO "oscRing"@.  The dynamical system
      involved is the oscillator system associated to $G$: one angle per vertex.
      If $a_{ij} = 1$ if $(i,j)$ is an edge of the undirected graph $G$, 
      and is zero otherwise, then the system is $d\theta_i/dt = \sum_j a_{ij} \sin(\theta_j - \theta_i)$
      where we consider only reduced equilibrium solutions $\theta_0 = 0$.
      
      This function returns the ideal of equilibrium points, where angles $(0, \theta_1, ..., \theta_{n-1})$
      are represented via cosines and sines of the angles.
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing(G, CoefficientRing => CC)
      S = oo
      I = oscSystemReduced(G,S)
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
  Headline
    create the Jacobian for the oscillator system associated to a graph
  Usage
    J = oscJacobian(G,S)
  Inputs
    G:Graph
      an undirected, connected graph $G$ on vertices $0, \ldots, n-1$, where $n$
      is the number of vertices of $G$
    S:Ring
      created with @TO oscRing@
  Outputs
    :Matrix
      the $n \times n$ Jacobian matrix of the system, as a matrix of polynomials
      involving the cosines and sines of angles $\theta_1, \ldots, \theta_{n-1}$, where we
      set $\theta_0 = 0$.
  Description
    Text
      The matrix is a symmetric $n \times n$  matrix (with determinant zero).
    Example
      G = graph({0,1,2,3}, {{0,1},{1,2},{2,3},{0,3}})
      oscRing(G, CoefficientRing => CC)
      S = oo
      I = oscSystemReduced(G,S)
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
    oscSystemReduced
    findRealSolutions
    identifyStability
    eigenvalues
///

doc ///
  Key
    findRealSolutions
    (findRealSolutions, Ideal)
  Headline
    find real solutions, at least one per component for well-conditioned systems
  Usage
    findRealSolutions I
  Inputs
    I:Ideal
      an ideal in a polynomial ring over QQ or RR or CC.
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
      oscRing(G, CoefficientRing => QQ)
      S = oo
      I = oscSystemReduced(G,S)
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
    oscSystemReduced
    oscJacobian
    identifyStability
///

doc ///
  Key
    oscQuadrics
    (oscQuadrics, Graph, Ring)
  Headline
    the quadrics of the oscillator system associated to a graph
  Usage
    oscQuadrics(G,R)
  Inputs
    G:Graph
      an undirected, connected graph $G$ on vertices $0, \ldots, n-1$, where $n$
      is the number of vertices of $G$
    R:Ring
      created with @TO oscRing@
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