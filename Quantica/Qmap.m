(* Mathematica Package *)

BeginPackage["Quantica`Qmap`"]
(* Exported symbols added here with SymbolName::usage *)
Quantica`Qmap::usage-"context of the quantum mapping"  
Evolve::usage="Evelve[vec]: |vec'>=U|vec>"
SymmetricEvolve::usage="SymmetricEvelve[vec]: |vec'>=U|vec>"
AbsorbedEvolve::usage="AbsorbedEvelve[vec]: |vec'>=PU|vec>"
AbsorbedEvolve::AbsorberError = "request matrix forms, e.g. {{vec1,vec2,...}}"
Unitary::usage="return unitary matrix <q'|U|q>"
SymmetricUnitary::usage="return symetric unitary matrix <q'|U|q>"
AbsorbedUnitary::usage="return absorbed unitary matrix <q'|PU|q>"

Begin["`Private`"] (* Begin Private Context *) 
dim = Quantica`Dim
hbar = Quantica`Hbar
tau = Quantica`Tau
q = Quantica`X[[1]]
p = Quantica`X[[2]]
dps = Quantica`MP`dps
domain = Quantica`Domain
{funcT,funcV} = Quantica`GetSystem[]

Evolve[vec_] := Module[{pd,qvec,pvec},
    pd=domain[[2]];
    pvec = SetPrecision[Fourier[Exp[ -I*funcV[q]*tau/hbar]*vec],dps];
    If[  pd[[1]]*pd[[2]] < 0, 
        qvec=InverseFourier[Exp[-I*RotateLeft[ funcT[p],dim/2 ]*tau/hbar]*pvec],
        qvec=InverseFourier[Exp[-I*funcT[p]*tau/hbar]*pvec]
    ];
    Return[ SetPrecision[qvec,dps] ]
]
SymmetricEvolve[vec_] := Module[{pd,qvec,pvec},
    pd=domain[[2]];
    pvec = SetPrecision[Fourier[ Exp[-I*funcV[q]*tau/hbar/2 ]*vec],dps];
    If[  pd[[1]]*pd[[2]] < 0, 
        qvec=InverseFourier[ Exp[-I*RotateLeft[ funcT[p], dim/2 ]*tau/hbar]*pvec],
        qvec=InverseFourier[ Exp[-I*funcT[p]*tau/hbar]*pvec]
    ];
    qvec = SetPrecision[qvec,dps];
    pvec = SetPrecision[Fourier[ Exp[ -I*funcV[q]*tau/hbar/2]*qvec],dps];
    qvec = SetPrecision[InverseFourier[ pvec ],dps];
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
Unitary[] := Module[{i,mat,eye},
    eye = N[IdentityMatrix[dim],dps];
    mat = Table[0,{dim,dim}];
    For[i=0,i<dim,i++;
        mat[[i]] = Evolve[eye[[i]]];
    ];
    Return[Transpose[mat]]
]
SymmetricUnitary[] := Module[{i,mat,eye},
    eye = N[IdentityMatrix[dim],dps];
    mat = Table[0,{dim,dim}];
    For[i=0,i<dim,i++;
        mat[[i]] = SymmetricEvolve[eye[[i]]];
    ];
    Return[Transpose[mat]]
]
AbsorbedUnitary[ab_,gamma_,isSymmetric_:True]:= Module[{i,mat,eye},
    eye = N[IdentityMatrix[dim], dps];
    mat = Table[0,{dim,dim}];
    For[i=0,i<dim,i++;
        mat[[i]] = AbsorbedEvolve[eye[[i]],ab,gamma,isSymmetric]
    ];
    Return[Transpose[mat]]
]

End[] (* End Private Context *)
EndPackage[]