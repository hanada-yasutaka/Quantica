 (* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/14 *)

BeginPackage["Quantica`FourierBase`"]
Quantica`FourierBase::usage="演算子の行列表現を座標基底で評価するcontextです"
MatrixV::usage="return Potential term <q'|V(q)|q>の位置表現を返します"
MatrixT::usage="return Kinetic term <q'|T(p)|q>の位置表現を返します"
Hamiltonian::usage="return Hamilton matrix <q'|H(q,p)|q>"



Begin["`Private`"]
dim = Quantica`Dim
hbar = Quantica`Hbar
q = Quantica`X[[1]]
p = Quantica`X[[2]]
dps = Quantica`MP`dps
domain = Quantica`Domain

{funcT,funcV} = Quantica`GetSystem[]

MatrixV[] := Module[{mat},
    mat=IdentityMatrix[dim]*funcV[q];
    Return[mat]
    ]

MatrixT[] := Module[{pd,cond,i,eye,mat,pbase,qbase},
    eye = N[ IdentityMatrix[dim], dps];
    mat = Table[0,{dim,dim}];
    pd = domain[[2]];
    cond=pd[[1]]*pd[[2]];
    For[ i=0, i<dim, i++;
        pbase = SetPrecision[ Fourier[ eye[[i]] ], dps];
        If[ cond < 0,
            qbase = InverseFourier[ pbase*RotateLeft[ funcT[p],dim/2] ],
            qbase = InverseFourier[ pbase*funcT[p] ]
        ];
        mat[[i]] = SetPrecision[qbase, dps];
    ];
    Return[ Transpose[mat] ]
    ]

Hamiltonian[] := MatrixT[] + MatrixV[] 

End[]
EndPackage[]