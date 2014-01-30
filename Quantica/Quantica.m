(* Mathematica Package *)

(* Created by the Wolfram Workbench 2013/12/11 *)

BeginPackage["Quantica`"]
Quantica::usage=
    "A Mathematica Package for calculations of quantum systems in multiple precision." 

Begin["Quantica`MP`"]
Quantica`MP::usage="context of the multiple precision settings"
dps::usage="the decimal precision (defalt 20)"
dps=20
End[]
(*
Get["Quantica`BCH`"]
*)
Get["Quantica`Systems`"]

(* infomation *)
QuanticaSetting::usage = "QuanticaSetting[dim, domain, tau=1 ],\n
dim: real integer\n
domain: 2x2 list such as {{qmin,qmax}, {pmin,pmax}}. 各変数の精度は\[Infinity]でなければならない\n,
tau:(optional) 精度は\[Infinity]でなければならない\n"
QuanticaSetting::dimError = "require real integer value, but you input `1`"
QuanticaSetting::domainError = "require 2x2 list such that {{qmin,qmax},{pmin,pmax}}"
QuanticaSetting::QdomainError = "require {qmin<qmax), but you input `1`"
QuanticaSetting::PdomainError = "require {pmin<pmax}, but you input `1`"
QuanticaSetting::precisionError = "don't use floting point value(1.0 or 2.0).\n use infinite precision representation such as 2, Sqrt[2] or 1/10."

SetSystem::usage="SetSystem[T,V]: Kinetic function and Potential function"
SetSystem::precisionError = "Precision of `1` is MachinePrecision"
SetSystem::defError = "Function `1` is not defined."
GetSystem::usage="return [T,V]: Kinetic and Potential term"
Eigen::usage="{evalues,evectors}=Eigen[mat,sort:True]\n return eigenvalue and eigenvectors"
ParallelEigen::usage=""

QuasiEnergy::usage="return quasi-energies"
SortIndex::usage="SortIndex[basis, evecs] <evecs_m|basos_n>が最大の値を返すmをn=0から順に求め，mのリストし返します"
SortEigen::usage="SortEigen[evals, evecs, index] indexの順番で並び替えされたevalsとevecsを返します"
InnerProduct::usage="InnerProduct[vec1,vec2]: return <vec1|vec2>"
InnerProducts::usage="InnerProducts[vec1,basis]:return {<basis[1]|vec2>,<basis[2]|vec2>,...}"
Overlap::usage="InnerProduct[vec1,vec2]: return |<vec1|vec2>|^2"


Dim::usage="Hilbert Space dimension (matrix dimension)"
Domain::usage="Domain of the phase space {{q interval}, {p interval}}"
Tau::usage="discretization of time(dt): (optional, defalt Tau=1)"
X::usage="coordinate of the position and momentum X = [q,p]"
Area::usage="Area of the phase space"
Planck::usage="effective Planck's constant"
Hbar::usage="effective Planck's constant divided by 2pi"

(*infomation of Pakages*)
Quantica`Help::usage="Helpを開きます"
Quantica`Help[] := NotebookOpen[ToFileName[{$UserBaseDirectory, "Applications", "Quantica","Documentation", "English", "Guides"}, "Quantica.nb"]]

(* constant *)
Dim=Dim
Domain=Domain
Tau=Tau
tau=tau
X=X
Area=Area  
Planck=Planck
Hbar=Hbar

Protect["tau"]
Begin["`Private`"]
(*Privateに入れないとNames["Quantica`アスタリスク"]に(局所)変数が表示されてしまうのだ*)

(*QuanticaSetting[dim_, domain_, tau_:1, Verbose_:False] :=  Module[ {},*)
Options[QuanticaSetting] = {tau->1, Verbose->False}
(*todo:
    optionの局所変数化*)
QuanticaSetting[dim_, domain_, OptionsPattern[]] :=  Module[ {},
    Which[
        Not[IntegerQ[dim]], 
                Message[QuanticaSetting::dimError,dim],
        TrueQ[dim<0],
                Message[QuanticaSetting::dimError,dim],
        Dimensions[domain]!={2,2},
                Message[QuanticaSetting::domainError,domain],
        Length[Select[Flatten[domain],Precision[#] == Infinity &]]!=4,
                Message[QuanticaSetting::precisionError,domain],
        Precision[OptionValue[tau]] != Infinity,
                Message[QuanticaSetting::precisionError,domain],
        domain[[1]][[1]]>=domain[[1]][[2]],
                Message[QuanticaSetting::QdomainError,domain[[1]]],
        domain[[2]][[1]]>=domain[[2]][[2]],
                Message[QuanticaSetting::PdomainError,domain[[2]]],
        True,
            Unprotect["Quantica`*"];
            Dim=dim; Domain=domain;X=setX[];
            Area=setArea[];Planck=setPlanck[];Tau=OptionValue[tau];Hbar=setHbar[];
            (*Quantica`MP`dps=OptionValue[Dps];*)            
            Print[{"Dim:", Dim,"Domain:",Domain,"Tau:",Tau,"dps:",Quantica`MP`dps}];
            If[Verbose==True,
                Print[{"Dim:", Dim,"Domain:",Domain,
                "Tau:",Tau,"dps:",Quantica`MP`dps}]
            ];
           
            Get["Quantica`Util`"];
            Get["Quantica`State`"];
            
    ];
    Protect["Quantica`*"];
]
GetSystem[]:= {FuncT, FuncV};
SetSystem[T_, V_] := Module[{},
    Which[
    Precision[ T[ X[[2]] ] ]==MachinePrecision,
        Message[SetSystem::precisionError, ToString[T]],
    Precision[ V[ X[[1]] ] ]==MachinePrecision,
        Message[SetSystem::precisionError, ToString[V]],
    Not[ ListQ[T[ X[[2]] ] ] ],
        Message[SetSystem::defError, ToString[T]],
    Not[ ListQ[V[ X[[1]] ] ] ],
        Message[SetSystem::defError, ToString[V]],
    True,
        Clear[FuncT,FuncV];        
        FuncT[x_]=T[x];FuncV[x_]=V[x];
        (*Get[Quantica`HarmonicBasep"]*)
        Get["Quantica`FourierBase`"];
        Get["Quantica`Qmap`"];        
    ];
]
Options[Eigen] = {Sort->True, Vector->True}
Eigen[mat_,OptionsPattern[]] := Module[{index,evals,evecs},
    If[ OptionValue[Vector],
        {evals, evecs} = Eigensystem[mat],
        evals=Eigenvalues[mat]
    ];
    If[OptionValue[Sort],
        index = Ordering[Re[evals]];
        evals = evals[[index]];
        If[OptionValue[Vector],
            evecs = evecs[[index]]
        ];
    ];
    If[OptionValue[Vector],
        Return[{evals,evecs}],
        Return[evals]
    ];
]
Options[ParallelEigen] = {Sort->True, Vector->True}
ParallelEigen[mats_,OptionsPattern[]]:=Module[{i,res,index},
	If[ OptionValue[Vector],
		res = ParallelMap[ Eigensystem, mats];	
		index = Ordering[ Re[res[[#]][[1]]] ] & /@ Range[Length[mats]];
		For[i=1,i<=Length[mats],i++,
			res[[i]][[1]] = res[[i]][[1]][[ index[[i]] ]];	
			res[[i]][[2]] = res[[i]][[2]][[ index[[i]] ]];			
		],
		res = ParallelMap[Eigenvalues, mats];
		index = Ordering[ Re[res[[#]] ] ] & /@ Range[Length[mats]];
		For[i=1,i<=Length[mats],i++,
			res[[i]] = res[[i]][[ index[[i]] ]]	
		];
	];
	Return[res];
]
	
QuasiEnergy[evals_]:= Hbar*I*Log[evals]/Tau

SortIndex[ref_,evecs_]:= Ordering [ Overlap[evecs, ref[[#]] ] ][[-1]] & /@ Range[Dim]

SortEigen[evals_,evecs_,index_]:=Module[{vals,vecs},
    vals=evals[[index]];
    vecs=evecs[[index]];
    Return[{vals,vecs}]
]

(* Operation *)
InnerProduct[vec1_,vec2_] := Inner[Times,Conjugate[vec1],vec2,Plus]
InnerProducts[vec_, basis_] := Table[InnerProduct[basis[[i]], vec],{i,1,Length[basis]}]
Overlap[vec1_,vec2_] := Module[{ovl},
    ovl = InnerProduct[vec1,vec2];
    Return[ Abs[ ovl*Conjugate[ovl] ] ]
]


(*Private functions*)
setX[] := Module[{q,p,dq,dp,i}, 
    dq = (Domain[[1]][[2]] - Domain[[1]][[1]])/Dim;
    dp = (Domain[[2]][[2]] - Domain[[2]][[1]])/Dim;
    q = N[Table[Domain[[1]][[1]]+i*dq,{i,0,Dim-1}], Quantica`MP`dps];
    p = N[Table[Domain[[2]][[1]]+i*dp,{i,0,Dim-1}], Quantica`MP`dps];
    Return[{q,p}]
]
setArea[]  := N[(Domain[[1]][[2]] - Domain[[1]][[1]])*(Domain[[2]][[2]] - Domain[[2]][[1]]), Quantica`MP`dps]
setPlanck[]:= N[Area/Dim,Quantica`MP`dps]
setHbar[]  := N[Planck/(2*Pi),Quantica`MP`dps]


End[](*end private*)
EndPackage[]

Column[{
"Quantica`"<>
"A Mathematica Package for calculations of quantum systems in arbitary precision. \n" <>
"by Yasutaka Hanada.\n",
"注意事項：\n"<>
"任意精度計算の為に，機械精度の浮動小数点表現(0.5や1.2等)を用いないで下さい．\n" <>
"有理数表現(1/2, 2/10)もしくは無理数表現(Sqrt[3])等を用いて下さい．\n" <>
"数値を代入する場合は精度Precision[x]が"<>ToString[Precision[1/2]]<>"となる様にして下さい．\n"<>
"\nMATHEMATICA "<>$Version,
"\nQuantica version 0.3.3 [beta]",
"\nTODAY IS "<>DateString[]
}]  
