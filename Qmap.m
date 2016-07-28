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
Qmap`SymplecticUnitary  :=  Qmap`Private`SymplecticUnitary
Qmap`SymplecticEvolve	:=	Qmap`Private`SymplecticEvolve
Qmap`SymmetricAbsorbedEvolve:=Qmap`Private`SymmetricAbsorbedEvolve
Qmap`SymmetricAbsorbedUnitary:=Qmap`Private`SymmetricAbsorbedUnitary
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

SymmetricEvolve[vec_,z_:1] := Module[{c,q,p,funcT,funcV,pd,qvec,pvec, input},
    q = Quantica`X[[1]];
    p = Quantica`X[[2]];
    c = z*{1/2, 1};
    {funcT,funcV} = GetSystem[];
    pd=Domain[[2]];
    pvec = FFT @ (Exp[-I*funcV[q]*Tau/Hbar*c[[1]] ]*vec);

    qvec = If[  pd[[1]]*pd[[2]] < 0, 
        IFFT @ (Exp[-I*RotateLeft[ funcT[p], Dim/2 ]*Tau/Hbar*c[[2]]]*pvec),
        IFFT @ (input = Exp[-I*funcT[p]*Tau/Hbar]*pvec)
    ];
    pvec = FFT@(Exp[-I*funcV[q]*Tau/Hbar*c[[1]] ]*qvec);
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
SymmetricAbsorbedEvolve[vec_,ab_,gamma_,isSymmetric_:True] := Module[{i,qvec},
    If[Not[MatrixQ[ab]],Message[AbsorbedEvolve::AbsorberError];Abort[]];
    qvec = vec;
    If[isSymmetric,
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec; ];
        qvec = SymmetricEvolve[qvec];
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec ]
        ,
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec ];        
        qvec = Evolve[qvec];
        For[i=0,i<Length[ab],i++; qvec -= gamma*ab[[i]]*Conjugate[ab[[i]]].qvec ]
    ];
    Return[qvec]
]


Unitary[] := Transpose @ Map[ Evolve[State`Unit[#]] &, Range[Quantica`Dim] ]
SymmetricUnitary[] := Transpose @ Map[ SymmetricEvolve[State`Unit[#]] &, Range[Quantica`Dim] ]
AbsorbedUnitary[ab_,gamma_,isSymmetric_:True] := Transpose @ Map[ AbsorbedEvolve[State`Unit[#], ab, gamma, isSymmetric] &, Range[Quantica`Dim] ]
SymmetricAbsorbedUnitary[ab_,gamma_,isSymmetric_:True] := Transpose @ Map[ SymmetricAbsorbedEvolve[State`Unit[#], ab, gamma, isSymmetric] &, Range[Quantica`Dim] ]


Symplectic4[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/3);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = SymmetricEvolve[vec, c[[1]] ];
	qvec = SymmetricEvolve[qvec, c[[2]] ]; 	     
	qvec = SymmetricEvolve[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic6[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/5);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic4[vec, c[[1]] ];
	qvec = Symplectic4[qvec, c[[2]] ]; 	     
	qvec = Symplectic4[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic8[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/7);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic6[vec, c[[1]] ];
	qvec = Symplectic6[qvec, c[[2]] ]; 	     
	qvec = Symplectic6[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic10[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/9);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic8[vec, c[[1]] ];
	qvec = Symplectic8[qvec, c[[2]] ]; 	     
	qvec = Symplectic8[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic12[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/11);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic10[vec, c[[1]] ];
	qvec = Symplectic10[qvec, c[[2]] ]; 	     
	qvec = Symplectic10[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic14[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/13);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic12[vec, c[[1]] ];
	qvec = Symplectic12[qvec, c[[2]] ]; 	     
	qvec = Symplectic12[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic16[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/15);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic14[vec, c[[1]] ];
	qvec = Symplectic14[qvec, c[[2]] ]; 	     
	qvec = Symplectic14[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic18[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/17);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic16[vec, c[[1]] ];
	qvec = Symplectic16[qvec, c[[2]] ]; 	     
	qvec = Symplectic16[qvec, c[[1]] ];
    Return[qvec]
]
Symplectic20[vec_,z_:1] := Module[{qvec,beta,c},
    beta = 2^(1/19);
    c = z*{1/(2-beta),-beta/(2-beta)};
	qvec = Symplectic18[vec, c[[1]] ];
	qvec = Symplectic18[qvec, c[[2]] ]; 	     
	qvec = Symplectic18[qvec, c[[1]] ];
    Return[qvec]
]



SymplecticEvolve[order_:1] := Module[{res},
	res = Which[order==1, Evolve,
		order==2,SymmetricEvolve,
		order==4,Symplectic4,
		order==6,Symplectic6,
		order==8,Symplectic8,
		order==10,Symplectic10,
		order==12,Symplectic12,
		order==14,Symplectic14,
		order==16,Symplectic16						
	];
	Return[res]
]
SymplecticUnitary[order_:1]:=Module[{evolve},
	evolve=SymplecticEvolve[order];
	Transpose @ Map[ evolve[State`Unit[#]] &, Range[Quantica`Dim] ]
]

End[] (* End Private Context *)
