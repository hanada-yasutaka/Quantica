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
SetSystem::usage="sss"
GetSystem::usage="aaa"

(*instantiation*)
Qmap`SetSystem			:=	Qmap`Private`SetSystem
Qmap`GetSystem			:=	Qmap`Private`GetSystem
Qmap`Evolve             :=  Qmap`Private`Evolve
Qmap`SymmetricEvolve    :=  Qmap`Private`SymmetricEvolve
Qmap`AbsorbedEvolve     :=  Qmap`Private`AbsorbedEvolve
Qmap`Unitary            :=  Qmap`Private`Unitary
Qmap`SymmetricUnitary   :=  Qmap`Private`SymmetricUnitary
Qmap`AbsorbedUnitary    :=  Qmap`Private`AbsorbedUnitary

End[]

Begin["Qmap`Private`"]

(*GetSystem[]:= {FuncT[#] &, FuncV[#] &};*)
GetSystem[]:= {FuncT, FuncV};
SetSystem[T_, V_] := Module[{},
	FunctionPrecision[#] &/@{T,V};
	FuncT[x_]=T[x];
	FuncV[x_]=V[x];
];

Evolve[vec_] := Module[{q,p,funcT,funcV,pd,qvec,pvec,input},
    q = Quantica`X[[1]];
    p = Quantica`X[[2]];
    {funcT,funcV} = GetSystem[];
    pd=Domain[[2]];
    pvec = FFT@ (Exp[-I*funcV[q]*Tau/Hbar]*vec);
    If[ pd[[1]]*pd[[2]] < 0,
		input = Exp[-I * RotateLeft[ funcT[p], Dim/2] * Tau / Hbar]*pvec,
        input = Exp[ -I * funcT[p]*Tau/Hbar]*pvec;
    ];
    
    qvec=IFFT[input];
    Return[ qvec ]
]

SymmetricEvolve[vec_] := Module[{q,p,funcT,funcV,pd,qvec,pvec, input},
    q = Quantica`X[[1]];
    p = Quantica`X[[2]];
    {funcT,funcV} = GetSystem[];
    pd=Domain[[2]];
    pvec = FFT @ (Exp[-I*funcV[q]*Tau/Hbar/2 ]*vec);

    If[  pd[[1]]*pd[[2]] < 0, 
        input = Exp[-I*RotateLeft[ funcT[p], Dim/2 ]*Tau/Hbar]*pvec,
        input = Exp[-I*funcT[p]*Tau/Hbar]*pvec
    ];
    qvec = IFFT[input];    
    pvec = FFT@(Exp[-I*funcV[q]*Tau/Hbar/2 ]*qvec);
	qvec = IFFT[pvec];    
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


Unitary[] := Transpose @ Map[ Evolve[State`Unit[#]] &, Range[Quantica`Dim] ]
SymmetricUnitary[] := Transpose @ Map[ SymmetricEvolve[State`Unit[#]] &, Range[Quantica`Dim] ]
AbsorbedUnitary[ab_,gamma_,isSymmetric_:True] := Transpose @ Map[ AbsorbedEvolve[State`Unit[#], ab, gamma, isSymmetric] &, Range[Quantica`Dim] ]

End[] (* End Private Context *)
