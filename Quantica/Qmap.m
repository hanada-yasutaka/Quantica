(* Mathematica Package *)


Begin["Qmap`"]
(* Exported symbols added here with SymbolName::usage *)

(*infomation*)
Qmap::usage-"context of the quantum mapping"  
Evolve::usage="Evelve[vec]: |vec'>=U|vec>"
SymmetricEvolve::usage="SymmetricEvelve[vec]: |vec'>=U|vec>"
AbsorbedEvolve::usage="AbsorbedEvelve[vec]: |vec'>=PU|vec>"
AbsorbedEvolve::AbsorberError = "request matrix forms, e.g. {{vec1,vec2,...}}"
Unitary::usage="return unitary matrix <q'|U|q>"
SymmetricUnitary::usage="return symetric unitary matrix <q'|U|q>"
AbsorbedUnitary::usage="return absorbed unitary matrix <q'|PU|q>"

(*instantiation*)
Qmap`Evolve             :=  Qmap`Private`Evolve
Qmap`SymmetricEvolve    :=  Qmap`Private`SymmetricEvolve
Qmap`AbsorbedEvolve     :=  Qmap`Private`AbsorbedEvolve
Qmap`Unitary            :=  Qmap`Private`Unitary
Qmap`SymmetricUnitary   :=  Qmap`Private`SymmetricUnitary
Qmap`AbsorbedUnitary    :=  Qmap`Private`AbsorbedUnitary

End[]

Begin["Qmap`Private`"]


Evolve[vec_] := Module[{q,p,funcT,funcV,pd,qvec,pvec},
(*	SetSharedVariable[q,p,funcT,funcV,pd];*)	
    q = Quantica`X[[1]];
    p = Quantica`X[[2]];
    {funcT,funcV} = Quantica`GetSystem[];
    pd=Domain[[2]];
    pvec = SetPrecision[Fourier[Exp[ -I*funcV[q]*Tau/Hbar]*vec],Quantica`MP`dps];
    If[  pd[[1]]*pd[[2]] < 0, 
        qvec=InverseFourier[Exp[-I*RotateLeft[ funcT[p],Dim/2 ]*Tau/Hbar]*pvec],
        qvec=InverseFourier[Exp[-I*funcT[p]*Tau/Hbar]*pvec]
    ];
    Return[ SetPrecision[qvec,Quantica`MP`dps] ]
]

SymmetricEvolve[vec_] := Module[{q,p,funcT,funcV,pd,qvec,pvec},
    q = Quantica`X[[1]];
    p = Quantica`X[[2]];
    {funcT,funcV} = Quantica`GetSystem[];
    pd=Domain[[2]];
    pvec = SetPrecision[Fourier[ Exp[-I*funcV[q]*Tau/Hbar/2 ]*vec],Quantica`MP`dps];
    If[  pd[[1]]*pd[[2]] < 0, 
        qvec=InverseFourier[ Exp[-I*RotateLeft[ funcT[p], Dim/2 ]*Tau/Hbar]*pvec],
        qvec=InverseFourier[ Exp[-I*funcT[p]*Tau/Hbar]*pvec]
    ];
    qvec = SetPrecision[qvec,Quantica`MP`dps];
    pvec = SetPrecision[Fourier[ Exp[ -I*funcV[q]*Tau/Hbar/2]*qvec],Quantica`MP`dps];
    qvec = SetPrecision[InverseFourier[ pvec ],Quantica`MP`dps];
    Return[qvec];
]
AbsorbedEvolve[vec_,ab_,gamma_,isSymmetric_:True] := Module[{i,qvec},
    If[Not[MatrixQ[ab]],Message[AbsorbedEvolve::AbsorberError];Abort[]];
    If[isSymmetric,
        qvec = SymmetricEvolve[vec];
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec ],
        qvec = Evolve[vec];
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec ]];
    Return[qvec]
]

(*Unitary[] := Transpose @ Map[ Evolve[State`Unit[#]] &, Range[Quantica`Dim] ]*)
Unitary[] := Transpose @ Map[ Evolve[State`Unit[#]] &, Range[Quantica`Dim] ]
SymmetricUnitary[] := Transpose @ Map[ SymmetricEvolve[State`Unit[#]] &, Range[Quantica`Dim] ]
AbsorbedUnitary[ab_,gamma_,isSymmetric_:True] := Transpose @ Map[ AbsorbedEvolve[State`Unit[#], ab, gamma, isSymmetric] &, Range[Quantica`Dim] ]
(*
SymmetricUnitary[] := Module[{i,mat,eye},
    eye = N[IdentityMatrix[Dim],Quantica`MP`dps];
    mat = Table[0,{Dim,Dim}];
    For[i=0,i<Dim,i++;
        mat[[i]] = SymmetricEvolve[eye[[i]]];
    ];
    Return[Transpose[mat]]
]
AbsorbedUnitary[ab_,gamma_,isSymmetric_:True]:= Module[{i,mat,eye},
    eye = N[IdentityMatrix[Dim], Quantica`MP`dps];
    mat = Table[0,{Dim,Dim}];
    For[i=0,i<Dim,i++;
        mat[[i]] = AbsorbedEvolve[eye[[i]],ab,gamma,isSymmetric]
    ];
    Return[Transpose[mat]]
]
*)
End[] (* End Private Context *)
