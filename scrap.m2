
needsPackage("NautyGraphs")
needsPackage "Visualize"


generateOscGraphs = method(Options => {noLeaves => true})
generateOscGraphs(ZZ) := List => opts -> n -> (
    Gstrs := flatten for i from n-1 to binomial(n,2)-1 list generateGraphs(n,i, OnlyConnected=>true);
    if opts.noLeaves then select(Gstrs / stringToGraph, hasNoLeaf) else Gstrs / stringToGraph
)


listStats = method()
listStats(Graph) := G -> (
    R = oscRing(G);
    {G, minimalPrimes oscQuadrics(G,R), factor vertexSpanningPolynomial G}
    )
listStats(List) := L -> L / listStats

checkConj = method()
checkConj(Graph) := G -> (
    R = oscRing(G);
    I = saturate(oscQuadrics(G,R), standardSols(G,R));
    assert(# decompose I  >= # factor vertexSpanningPolynomial G)
)
checkConj(List) := L -> L / checkConj

rationalCircleParametrization = method()
rationalCircleParametrization(Ideal) := Ideal => I -> (
    R := ring I;
    G := R.Graph;
    vertexList := if R.Reduced then drop(vertices G, 1) else vertices G;
    tlist := for v in vertexList list t_v;
    S := QQ[tlist,z];
    f := for v in vertexList list (select(R.generators, k -> toString k == toString (R.Symbols_0)_v))_0 => z^2-t_v^2;
    g := for v in vertexList list (select(R.generators, k -> toString k == toString (R.Symbols_1)_v))_0 => t_v;
    phi := map(S, R, f | g );
    phi(I)
)

matrixTreeGeneralization = method()
matrixTreeGeneralization(Graph, Ring) := (G,R) -> (
    vertexList := if R.Reduced then drop(vertices G, 1) else vertices G;
    varList := for v in vertexList list l_v;
    S := R[varList];
    cosine := (i) -> (select(R.generators, k -> toString k == toString (R.Symbols_0)_i))_0;
    VG := matrix for i in vertices G list for j in vertices G list (
        if i === j then (
            sum for n in toList neighbors(G,i) list cosine n
        ) else if member(j, neighbors(G,i)) then - cosine i else 0
    );
    det(VG + diagonalMatrix(vars S))
)
matrixTreeGeneralization(Graph) := G -> matrixTreeGeneralization(G, oscRing(G))

texResults = method()
texResults(List, String) := (Gstats, fileName) -> (
    results = openOut(baseDirectory | "/examples/" | fileName);
    results << "\\documentclass[12pt]{amsart}" << endl;
    results << "\\usepackage{amsmath, amsfonts, amsthm, amssymb}" << endl;
    results << "\\usepackage[top=1in, bottom=1in, left=1in, right=1in]{geometry}" << endl;
    results << "\\usepackage{breqn}" << endl;
    results << "\\usepackage{tkz-graph}" << endl << endl;

    results << "\\begin{document}" << endl;
    for G in Gstats do (
        results << "\\begin{center}" << endl;
        results << "\\scalebox{0.75}{" << endl;
        results << "\\begin{tikzpicture}" << endl;
        results << "\\SetGraphUnit{2}" << endl;
        results << "\\Vertices{circle}{";
        vertexlist = vertices G_0;
        for i from 0 to #vertexlist-2 do results << vertexlist_i << ",";
        results << vertexlist_(#vertexlist-1) << "}" << endl;
        for e in edges G_0 do results << "\\Edge(" << (toList e)_0 << ")(" << (toList e)_1 << ")" << endl;
        results << "\\end{tikzpicture}" << endl;
        results << "}" << endl;
        results << "\\end{center}" << endl;

        results << "\\texttt{decompose} $I_G : I_\\Sigma^\\infty$: \\\\" << endl ;
        for i in G_1 do (
            results << "\\begin{dmath*}" << endl;
            latex := tex i;
            strippedlatex := substring(1,#latex-2,latex);
            results << strippedlatex << endl;
            results << "\\end{dmath*}" << endl; 
        );
        results << "$F_G$: \\\\" << endl ;
        results << "\\begin{dmath*}" << endl;
        latex := tex G_2;
        strippedlatex := substring(1,#latex-2,latex);
        results << strippedlatex << endl;
        results << "\\end{dmath*}" << endl; 
        results << endl << "\\hspace{1cm}" << endl;
    );
    results << "\\end{document}" << endl;
    closeOut(results);
)
texResults(List) := Gstats -> texResults(List, "results.txt")

texGraphs = method(Options => {numOscComponents => false})
texGraphs(List, String) := opts -> (Gstats, fileName) -> (
    results = openOut(baseDirectory | "/examples/" | fileName);
    results << "\\documentclass[12pt]{amsart}" << endl;
    results << "\\usepackage{amsmath, amsfonts, amsthm, amssymb}" << endl;
    results << "\\usepackage[top=1in, bottom=1in, left=1in, right=1in]{geometry}" << endl;
    results << "\\usepackage{breqn}" << endl;
    results << "\\usepackage{tkz-graph}" << endl << endl;

    results << "\\begin{document}" << endl;
    for G in Gstats do (
        results << "\\begin{center}" << endl;
        results << "\\scalebox{0.75}{" << endl;
        results << "\\begin{tikzpicture}" << endl;
        results << "\\SetGraphUnit{2}" << endl;
        results << "\\Vertices{circle}{";
        vertexlist = vertices G_0;
        for i from 0 to #vertexlist-2 do results << vertexlist_i << ",";
        results << vertexlist_(#vertexlist-1) << "}" << endl;
        for e in edges G_0 do results << "\\Edge(" << (toList e)_0 << ")(" << (toList e)_1 << ")" << endl;
        results << "\\end{tikzpicture}" << endl;
        results << "}" << endl;
        results << "\\end{center}" << endl;
        if opts.numOscComponents then (
            results << "Number of Oscillator Components: " << #G_1 << endl;
        );
    );
    results << "\\end{document}" << endl;
    closeOut(results);
)

updateStats = method()
updateStats(Graph) := G -> (
    R = oscRing(G);
    I = saturate(oscQuadrics(G,R), standardSols(G,R));
    G.cache.oscIdeal = decompose I;
    G.cache.vertexSpanningPolynomial = factor vertexSpanningPolynomial G;
)

