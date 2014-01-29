(* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/14 *)

Begin["State`"]
State::usage="context of the quantum states functions and operation"
(*information*)
State`Zeros::usage="Zero[] return zero list"
State`Unit::usage="Unit[n] :return　(n番目の成分が1の) unit vector\n
n:integer\n
Note:listのindexは1から始まる事に注意!"
State`CS::usage="CS[qc,pc]: return coherent state centered at (qc,pc)\n
qc,pc: number\n
Note:周期境界条件は課していません．｡
"

MirrorTr::usage="MirrorTr[vec]: mirror transformation. MirrorTr[vec(q)]:->vec(-q)"
MirrorParity::usage="MirrorParity[vec]: Re[<vec(-q)|vec(q)>]"
TranslationTr::usage="TranslationTr[vec,shift]:translation transformation. TranslationTr[vec(q)]:-> vec(q+shift)"
TranslationParity::usage="TranslationParity[vec,shift:Dim/2]:Re[<vec(q+shift)|vec(q)>]"
Q2P::usage="Q2P[vec(q)]: Q2P[<q|vec>] :-> <p|vec>"
P2Q::usage="P2Q[vec(q)]: P2Q[<p|vec>] :-> <q|vec>"
Abs2::usage="Abs2[vec]:=|<x|vec>|^2"
Linear::usage="linear..."
FFT::usage="FFT[vec], Forward Fourier Transformation"
IFFT::usage="IFFT[vec],Inverse Fourier Transformation"

(*instantiation*)
State`Zeros             :=  State`Private`Zeros
State`Unit              :=  State`Private`Unit
State`CS                :=  State`Private`CS
State`MirrorTr          :=  State`Private`MirrorTr
State`MirrorParity      :=  State`Private`MirrorParity
State`TranslationTr     :=  State`Private`TranslationTr
State`TranslationParity :=  State`Private`TranslationParity
State`Q2P               :=  State`Private`Q2P
State`P2Q               :=  State`Private`P2Q
State`Abs2              :=  State`Private`Abs2
State`Linear			:=	State`Private`Linear
State`FFT				:=	State`Private`FFT
State`IFFT				:=	State`Private`IFFT

End[]
Begin["State`Private`"]
testPrec::precisionError = "input value is `1`, usege infinite precision or > dps "
testPrec[x_]:= If[ Precision[x] < Quantica`MP`dps, Message[testPrec::precisionError,x];Abort[],True]

Zeros[] := PadLeft[{}, Dim]
Unit[n_Integer]:= Module[{x},
    x = Zeros[];
    x[[n]] = N[1,Quantica`MP`dps];
    Return[x]
]
CS[qc_,pc_]:= Module[{vec,norm,q},
    q = Quantica`X[[1]];
    vec = Exp[-(q-qc)^2/(2*Hbar)+I*pc*(q-qc)/Hbar];
    norm=Abs[Inner[Times,Conjugate[vec],vec,Plus]]; (*replace Normalize*)
    Return[vec/Sqrt[norm]]
]
Linear[pc_, k_, omega_]:= Module[{q,h, pre,vec,x,norm},
	Map[testPrec, {pc,k,omega}];
	q = Quantica`X[[1]];
	h = Quantica`Planck;
	pre = -k/(16*Pi*Pi*Sin[Pi*omega]);
	x = pre*Sin[2*Pi*(q - omega/2)] + q*pc;
	vec = Exp[I*2*Pi/h * x];
    norm=Abs[Inner[Times,Conjugate[vec],vec,Plus]]; (*replace Normalize*)
	Return[vec/Sqrt[norm]];
]
FFT[vec_]  := SetPrecision[Fourier        [ SetPrecision @@ {vec, Quantica`MP`dps}], Quantica`MP`dps] 
IFFT[vec_] := SetPrecision[InverseFourier [ SetPrecision @@ {vec, Quantica`MP`dps}], Quantica`MP`dps] 

Q2P[vec_]:=Module[{pd,pvec},
    pd=Domain[[2]];
    pvec = FFT[vec];
    If[pd[[1]]*pd[[2]]<0, pvec = RotateLeft[ pvec , Length[pvec]/2]];
    Return[pvec];
]
P2Q[vec_]:=Module[{pd, vec1},
    pd=Domain[[2]];	
    If[pd[[1]]*pd[[2]]<0, vec1 = RotateLeft[ vec , Length[ vec ]/2], vec1=vec ];
    (*
	invec = SetPrecision[vec, Quantica`MP`dps];        
    outvec = SetPrecision[InverseFourier[invec], Quantica`MP`dps];
    *)
    Return[IFFT[vec1]]
]
InnerProduct:=Quantica`InnerProduct
MirrorTr[vec_] := Reverse[vec]
MirrorParity[vec_] := Module[{vec1,vec2,inner},
    vec1 = Normalize[Append[vec,vec[[1]]]];
    vec2 = MirrorTr[vec1];
    inner = InnerProduct[vec1,vec2];
    Return[Re[inner]]
]

TranslationTr[vec_,shift_] := RotateLeft[vec, shift]
TranslationParity[vec_, n_:-1] := Module[{vec1,inner,shift},
    If[n==-1,shift=Dim/2,shift=n];
    vec1 = TranslationTr[vec,shift];
    inner = InnerProduct[vec1,vec];    
    Return[Re[inner]]
]
Abs2[vec_]:=Abs[Conjugate[vec]*vec]

End[]
