 (* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/14 *)


Begin["PositionBase`"]
(*infomation*)
Quantica`PositionBase::usage="演算子の行列表現を座標基底で評価するcontextです"
(*
MatrixV::usage="return Potential term <q'|V(q)|q>の位置表現を返します"
MatrixT::usage="return Kinetic term <q'|T(p)|q>の位置表現を返します"
*)
Hamiltonian::usage="return Hamilton matrix <q'|H(q,p)|q>"
Func2Matrix::usage="Func2Matrix[f,x] xはQ or P でなければならない"
Func2Matrix::RepError="usage: Func2Matrix[f, Q (or P)], x =`1`"

(*instantiation*)
(*
PositionBase`MatrixV     := PositionBase`Private`MatrixV
PositionBase`MatrixT     := PositionBase`Private`MatrixT
*)
PositionBase`Hamiltonian := PositionBase`Private`Hamiltonian
PositionBase`Func2Matrix	:= PositionBase`Private`Func2Matrix
PositionBase`ParallelFunc2Matrix	:= PositionBase`Private`ParallelFunc2Matrix
PositionBase`AbsorbedHamiltonian :=PositionBase`Private`AbsorbedHamiltonian
End[]

Begin["PositionBase`Private`"]
Func2Matrix[f_, x_]:=
	Which[
	 SymbolName[x]=="Q",MatrixQFunction[f],
	 SymbolName[x]=="P",MatrixPFunction[f],
	 True,Message[PositionBase`Func2Matrix::RepError,x]
	 ]
ParallelFunc2Matrix[f_, x_]:=
	Which[
	  SymbolName[x]=="Q",MatrixQFunction[f],
	 	SymbolName[x]=="P",ParallelMatrixPFunction[f],
	 	True,Message[PositionBase`Func2Matrix::RepError,x]
	]


MatrixQFunction[F_]:=Module[{mat,q},
    q = Quantica`X[[1]];
    mat= N[ IdentityMatrix[Dim]*F[q],Quantica`MP`dps];
    Return[mat]
]
MatrixPFunction[F_]:=Module[{mat,eye,pd,cond,qbase,pbase,i,p},
    p = Quantica`X[[2]];
    eye = N[ IdentityMatrix[Dim], Quantica`MP`dps];
    mat = Table[0,{Dim,Dim}];
    pd = Domain[[2]];
    cond=pd[[1]]*pd[[2]];
    For[ i=0, i<Dim, i++;
        pbase = SetPrecision[ Fourier[ eye[[i]] ], Quantica`MP`dps];
        If[ cond < 0,
            qbase = InverseFourier[ pbase*RotateLeft[ F[p],Dim/2] ],
            qbase = InverseFourier[ pbase*F[p] ]
        ];
        mat[[i]] = SetPrecision[qbase, Quantica`MP`dps];
    ];
    Return[ Transpose[mat] ]
]
ParallelMatrixPFunction[F_]:=Module[{mat,eye,pd,cond,qbase,pbase,i,p},
    p = Quantica`X[[2]];
    eye = N[ IdentityMatrix[Dim], Quantica`MP`dps];
    (*mat = Table[0,{Dim,Dim}];*)
    pd = Domain[[2]];
    cond=pd[[1]]*pd[[2]];
		SetSharedVariable[eye, cond, Dim, Quantica`MP`dps];
    mat = ParallelTable[
			pbase = SetPrecision[ Fourier[ eye[[i]] ], Quantica`MP`dps];
			If[ cond < 0,
					qbase = InverseFourier[ pbase*RotateLeft[ F[p],Dim/2] ],
					qbase = InverseFourier[ pbase*F[p] ]
			];
			SetPrecision[qbase, Quantica`MP`dps]
    ,{i,1,Dim}];
		UnsetShared[eye, cond, Dim, Quantica`MP`dps];
    Return[ Transpose[mat] ]
]


(*to do
ParallelMatrixPFunction[F_]:=Module[{mat,eye,pd,cond,qbase,pbase,i,p},
    p = Quantica`X[[2]];
    eye = N[ IdentityMatrix[Dim], Quantica`MP`dps];
    mat = Table[0,{Dim,Dim}];
    pd = Domain[[2]];
    cond=pd[[1]]*pd[[2]];
    For[ i=0, i<Dim, i++;
        pbase = SetPrecision[ Fourier[ eye[[i]] ], Quantica`MP`dps];
        If[ cond < 0,
            qbase = InverseFourier[ pbase*RotateLeft[ F[p],Dim/2] ],
            qbase = InverseFourier[ pbase*F[p] ]
        ];
        mat[[i]] = SetPrecision[qbase, Quantica`MP`dps];
    ];
    Return[ Transpose[mat] ]
]
*)


AbsorbedHamiltonian[Ham_,ab_,gamma_,isSymmetric_:True] := Module[{i,res,eye,j,mat},
    If[Not[MatrixQ[ab]],Message[AbsorbedEvolve::AbsorberError];Abort[]];
	eye = N[ IdentityMatrix[Dim], Quantica`MP`dps];
	mat = ConstantArray[0,{Dim,Dim}];
    For[j=0,j<Dim,j++;
		res = ConstantArray[0,Dim];
    	For[i=0,i<Length[ab],i++;
			res += I*gamma*ab[[i]]*Conjugate[ab[[i]]]*Ham.eye[[j]];
    	];
    	mat[[j]] = res;
    ];
    res = Ham - Transpose@mat;

    Return[res]
]
(*
MatrixV[] := Module[{mat,q,funcV,funcT},
    q = Quantica`X[[1]];
    {funcT,funcV} = Quantica`GetSystem[];
    mat= N[ IdentityMatrix[Dim]*funcV[q],Quantica`MP`dps];
    Return[mat]
    ]

MatrixT[] := Module[{pd,cond,i,eye,mat,pbase,qbase,p,funcT,funcV},
    p = Quantica`X[[2]];
    {funcT,funcV} = Quantica`GetSystem[];
    eye = N[ IdentityMatrix[Dim], Quantica`MP`dps];
    mat = Table[0,{Dim,Dim}];
    pd = Domain[[2]];
    cond=pd[[1]]*pd[[2]];
    For[ i=0, i<Dim, i++;
        pbase = SetPrecision[ Fourier[ eye[[i]] ], Quantica`MP`dps];
        If[ cond < 0,
            qbase = InverseFourier[ pbase*RotateLeft[ funcT[p],Dim/2] ],
            qbase = InverseFourier[ pbase*funcT[p] ]
        ];
        mat[[i]] = SetPrecision[qbase, Quantica`MP`dps];
    ];
    Return[ Transpose[mat] ]
    ]
*)
Hamiltonian[T_,V_] := Func2Matrix[T,Symbol["P"]] + Func2Matrix[V,Symbol["Q"]]

End[]
