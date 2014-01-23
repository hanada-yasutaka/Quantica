 (* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/19 *)



Begin["Util`"]
(*infomation*)
Util::usage="context of the utility(I/O) operation"
SaveState::usage="SaveState[x,vector,fname]"
SaveEigenvalue::usage="SaveEigenvalue[eigenvalues,fname]"
SaveEigen::usage="SaveEigen[eigenvalues,eigenvectors,fname->eigen, verbose->True]\n"
SaveSplitting::usage"SaveSplitting[energy] save energy splitting"

(*instantiation*)
Util`SaveState      := Util`Private`SaveState
Util`SaveEigenvalue := Util`Private`SaveEigenvalue
Util`SaveEigen      := Util`Private`SaveEigen
Util`SaveSplitting  := Util`Private`SaveSplitting

End[]
Begin["Util`Private`"]

(*
Domain = Quantica`Domain
Planck = Quantica`Planck 
*)
WriteScaleInfo[str_] := Module[{d},
    d = Domain;
    Export[str,"# date: "<>DateString[] <>" \n"];
    Export[str,"# QMIN: "<> ToString[ N[d[[1]][[1]],10] ] <>" \n"];
    Export[str,"# QMAX: "<> ToString[ N[d[[1]][[2]],10] ] <>" \n"];
    Export[str,"# PMIN: "<> ToString[ N[d[[2]][[1]],10] ] <>" \n"];
    Export[str,"# PMAX: "<> ToString[ N[d[[2]][[2]],10] ] <>" \n"];
    Export[str,"# PLANCK: "<> ToString[N[Planck,10]] <>" \n"];
    Export[str,"# DIM: "<> ToString[Dim] <> " \n"];
]
Options[SaveState] = {Verbose->False}
SaveState[x_,evec_,fname_,OptionsPattern[]] := Module[{str},
    str=OpenWrite[fname];
    WriteScaleInfo[str];
    Export[str,"# x, |<x|vec>|^2, Re[<x|vec>], Im[<x|vec>] \n"];        
    Export[str,Transpose[{x,Abs[evec*Conjugate[evec]],Re[evec],Im[evec]}]];
    Close[str];
    If[OptionValue[Verbose],Print["outdate: "<>fname]]
]
Options[SaveEigenvalue] = {Verbose->False}
SaveEigenvalue[evals_,fname_,OptionsPattern[]] := Module[{str,ene,abs},
    str=OpenWrite[fname];
    WriteScaleInfo[str];
    ene = I*Hbar*Log[evals]/Tau;
    abs = Abs[Conjugate[evals]*evals];
    Export[str,"# index, Re[val], Im[val], Re[q-ene], Im[q-ene], |conj(val)val|^2 \n"];
    Export[str,Transpose[{Range[0,Length[evals]-1], Re[evals], Im[evals], Re[ene], Im[ene], abs }] ];
    Close[str];
    If[OptionValue[Verbose],Print["outdata: "<>fname]]
    ]

Options[SaveEigen] = {Head->"eigen",Verbose->False}
SaveEigen[evals_,evecs_,OptionsPattern[]] := Module[{fname,i},
    fname=OptionValue[Head] <>"_evals.dat";
    SaveEigenvalue[evals,fname, Verbose->OptionValue[Verbose]];
    For[i=0,i<Length[evals],i++;
        fname=OptionValue[Head] <>"_qrep_"<> ToString[i-1] <> ".dat";
        SaveState[X[[1]], evecs[[i]], fname,Verbose->OptionValue[Verbose]];
    ]
]
MirrorParity := State`MirrorParity
TranslationParity := State`TranslationParity
Options[SaveSplitting] = {Mirror->False, Translation->False, Fname->"splitting.dat"}
SaveSplitting[energy_,vecs_, OptionsPattern[]] := Module[{i, fop, splitting, parity0=" ",parity1=" ",trans0=" ",trans1=" "},
    fop=OpenWrite[OptionValue[Fname]];
    Export[fop,"# doublet index, splitting (low ene parity), (high ene parity)\n"];
    For[i=0,i<Dim,i++;
        If[ Mod[i,2]==0,
            splitting = Abs [ energy[[i-1]] - energy[[i]] ];
            If[OptionValue[Mirror], 
                If[MirrorParity[ vecs[[i-1]] ]>0, parity0="+", parity0="-"]
            ];
            If[OptionValue[Mirror],
                If[MirrorParity[ vecs[[i]]   ]>0, parity1="+", parity1="-"]
            ];
            If[OptionValue[Translation],
                If[TranslationParity[ vecs[[i-1]] ]>0, trans0="+", trans0="-"]
            ];
            If[OptionValue[Translation],
                If[TranslationParity[ vecs[[i]]   ]>0, trans1="+", trans1="-"]
            ];
            Export[fop, {{i/2, splitting, {parity0,trans0}, {parity1,trans1} }}];
            Export[fop,"\n"];
        ];
    ];
    Close[fop];
]

End[]