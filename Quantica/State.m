(* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/14 *)

BeginPackage["Quantica`State`"]
Quantica`State::usage="context of the quantum states functions and operation"

(*information*)
Zeros::usage="Zero[] return zero list"
Unit::usage="Unit[n] :return　(n番目の成分が1の) unit vector\n
n:integer\n
Note:listのindexは1から始まる事に注意!"
CS::usage="CS[qc,pc]: return coherent state centered at (qc,pc)\n
qc,pc: number\n
Note:周期境界条件は課していません．｡
"
MirrorTr::usage="MirrorTr[vec]: mirror transformation. MirrorTr[vec(q)]:->vec(-q)"
MirrorParity::usage="MirrorParity[vec]: Re[<vec(-q)|vec(q)>]"
TranslationTr::usage="TranslationTr[vec,shift]:translation transformation. TranslationTr[vec(q)]:-> vec(q+shift)"
TranslationParity::usage="TranslationParity[vec,shift:dim/2]:Re[<vec(q+shift)|vec(q)>]"
Q2P::usage="Q2P[vec(q)]: Q2P[<q|vec>] :-> <p|vec>"
P2Q::usage="P2Q[vec(q)]: P2Q[<p|vec>] :-> <q|vec>"
Abs2::usage="Abs2[vec]:=|<x|vec>|^2"

Begin["`Private`"]
dim=Quantica`Dim
hbar=Quantica`Hbar
q = Quantica`X[[1]]
p = Quantica`X[[2]]
dps = Quantica`MP`dps
domain=Quantica`Domain

Zeros[] := Table[0, {dim}]
Unit[n_]:= Module[{x},
    x = Zeros[];
    x[[n]] = N[1,dps];
    Return[x]
]
CS[qc_,pc_]:= Module[{vec,norm},
    vec = Exp[-(q-qc)^2/(2*hbar)+I*pc*(q-qc)/hbar];
    norm=Abs[Inner[Times,Conjugate[vec],vec,Plus]];
    Return[vec/Sqrt[norm]]
]

Q2P[vec_]:=Module[{pd,pvec},
    pd=domain[[2]];
    pvec = SetPrecision[Fourier[vec], dps];
    If[pd[[1]]*pd[[2]]<0, pvec = RotateLeft[ pvec , Length[vec]/2] ];
    Return[pvec]
]
P2Q[vec_]:=Module[{pd,vec1},
    pd=domain[[2]];
    If[pd[[1]]*pd[[2]]<0, vec1 = RotateLeft[ vec , Length[vec]/2] ];    
    vec1 = SetPrecision[InverseFourier[vec1], dps];
    Return[vec1]
]
InnerProduct=Quantica`InnerProduct
MirrorTr[vec_] := Reverse[vec]
MirrorParity[vec_] := Module[{vec1,vec2,inner},
    vec1 = Normalize[Append[vec,vec[[1]]]];
    vec2 = MirrorTr[vec1];
    inner = InnerProduct[vec1,vec2];
    Return[Re[inner]]
]

TranslationTr[vec_,shift_] := RotateLeft[vec, shift]
TranslationParity[vec_, n_:-1] := Module[{vec1,inner,shift},
    If[n==-1,shift=dim/2,shift=n];
    vec1 = TranslationTr[vec,shift];
    inner = InnerProduct[vec1,vec];    
    Return[Re[inner]]
]
Abs2[vec_]:=Abs[Conjugate[vec]*vec]
End[]
EndPackage[]