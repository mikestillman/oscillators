-- This is a temporary file until I figure out what to do with this code. 



---------------------------
-- Versions with weights --
-- need to check these   --
---------------------------
oscSystem = method()
oscSystem(Graph, Ring, HashTable) := (G, R, wts) -> (
    -- 
    n := # vertices G;  -- R should have 2n + #params variables.
    if toList(0..n-1) =!= sort vertices G then
      error "got incompatible graph";
    sine := (i) -> R_(n+i);
    cosine := (i) -> R_(i);
    ideal for i from 0 to n list (
        N := toList neighbors(G,i);
        sum for j in N list (
            e := sort {i,j};
            wt := if wts#?e then wts#e else 1_ZZ;
            --<< "(e,wt) = " << e << " " << wt << endl;
            wt*(sine j * cosine i - sine i * cosine j)
            )
        )
    )
oscSystemReduced = method()
oscSystemReduced(Graph, Ring, HashTable) := (G, R, wts) -> (
    -- 
    n := # vertices G - 1;  -- R should have 2n + #params variables.
    if toList(0..n) =!= sort vertices G then
      error "got incompatible graph";
    sine := (i) -> if i === 0 then 0_R else R_(n+i-1);
    cosine := (i) -> if i === 0 then 1_R else R_(i-1);
    ideal for i from 0 to n list (
        N := toList neighbors(G,i);
        sum for j in N list (
            e := sort {i,j};
            wt := if wts#?e then wts#e else 1_ZZ;
            --<< "(e,wt) = " << e << " " << wt << endl;
            wt*(sine j * cosine i - sine i * cosine j)
            )
        )
    )
oscSystemReduced(Graph, Ring) := (G, R) -> (
    n := # vertices G - 1;  -- R should have == 2n variables.
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

ringOscillatorGraph = method()
ringOscillatorGraph(ZZ,ZZ) := (n,k) -> (
    graph unique flatten for i from 0 to n-1 list (
        for j from i-k to i+k list (
            if j == i then continue;
            j1 := if j < 0 then j + n else if j >= n then j-n else j;
            sort {i,j1}
            )
        )
    )

egRingOsc = method(Options => options oscRing)
egRingOsc(ZZ,ZZ) := opts -> (n,k) -> (
    G := ringOscillatorGraph(n,k);
    R := oscRing(n-1,{},opts);
    P := oscSystemReduced(G,R,hashTable{});
    J := oscJacobian P;
    (J, P + trig R)
    )

egRingOscWeighted = method(Options => options oscRing
        )
egRingOscWeighted(ZZ,ZZ,Symbol) := opts -> (n,k,symb) -> (
    G := ringOscillatorGraph(n,k);
    R := oscRing(n-1,{symb},opts);
    a := R_(2*(n-1));
    wts := flatten for k0 from 2 to k list
      for i from 0 to n-1 list sort {i,if i+k0 >= n then i+k0-n else i+k0} => a^(k0-1);
    P := oscSystemReduced(G,R,hashTable wts);
    J := oscJacobian P;
    (J, P + trig R)
    )