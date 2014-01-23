 (* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/14 *)


Begin["FourierBase`"]
(*infomation*)
Quantica`FourierBase::usage="演算子の行列表現を座標基底で評価するcontextです"
MatrixV::usage="return Potential term <q'|V(q)|q>の位置表現を返します"
MatrixT::usage="return Kinetic term <q'|T(p)|q>の位置表現を返します"
Hamiltonian::usage="return Hamilton matrix <q'|H(q,p)|q>"

(*instantiation*)
FourierBase`MatrixV     := FourierBase`Private`MatrixV
FourierBase`MatrixT     := FourierBase`Private`MatrixT
FourierBase`Hamiltonian := FourierBase`Private`Hamiltonian
End[]

Begin["FourierBase`Private`"]


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

Hamiltonian[] := MatrixT[] + MatrixV[] 

End[]
