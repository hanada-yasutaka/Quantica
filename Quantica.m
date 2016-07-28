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
(*
SetSystem::usage="SetSystem[T,V]: Kinetic function and Potential function"
SetSystem::precisionError = "Precision of `1` is MachinePrecision"
SetSystem::defError = "Function `1` is not defined."
GetSystem::usage="return [T,V]: Kinetic and Potential term"
*)
Eigen::usage="{evalues,evectors}=Eigen[mat,sort:True]\n return eigenvalue and eigenvectors"
ParallelEigen::usage=""

QuasiEnergy::usage="return quasi-energies"
SortIndex::usage="SortIndex[basis, evecs] <evecs_m|basos_n>が最大の値を返すmをn=0から順に求め，mのリストし返します"
SortEigen::usage="SortEigen[evals, evecs, index] indexの順番で並び替えされたevalsとevecsを返します"


Dim::usage="Hilbert Space dimension (matrix dimension)"
Domain::usage="Domain of the phase space {{q interval}, {p interval}}"
Tau::usage="discretization of time(dt): (optional, defalt Tau=1)"
X::usage="coordinate of the position and momentum X = [q,p]"
PSArea::usage="Area of the phase space"
Planck::usage="effective Planck's constant"
Hbar::usage="effective Planck's constant divided by 2pi"

(*infomation of Pakages*)
Quantica`Help::usage="Helpを開きます"
Quantica`Help[] := NotebookOpen[ToFileName[{$UserBaseDirectory, "Applications", "Quantica","Documentation", "English", "Guides"}, "Quantica.nb"]]

Linspace[xmin_,xmax_,n_Integer,endpoint_:False]:=Module[{dx,res},
	If[endpoint,
		dx=(xmax-xmin)/(n-1);
		res=Range[xmin,xmax,dx],
		dx=(xmax-xmin)/n;
		res=Range[xmin, xmax-dx,dx]
	];
	Return[res];
]
FFT[vec_]  := SetPrecision[Fourier        [ SetPrecision @@ {vec, Quantica`MP`dps},FourierParameters->{1,-1}], Quantica`MP`dps] 
IFFT[vec_] := SetPrecision[InverseFourier [ SetPrecision @@ {vec, Quantica`MP`dps},FourierParameters->{1,-1}], Quantica`MP`dps] 

(* constant *)
Dim=Dim
Domain=Domain
Tau=Tau
X=X
PSArea=PSArea
Planck=Planck
Hbar=Hbar
Protect[Dim, Domain,Tau, X, PSArea, Planck, Hbar,tau]
Protect[Q,P]

FunctionPrecision::PrecisionError= "関数の精度 (`1`) はです．目標精度(Quantica'MP'dps)より低いため評価を中断しました"    
FunctionPrecision[f_]:=Module[{x},
	If[ Precision@f[x] < Quantica`MP`dps,
		 Message[FunctionPrecision::PrecisionError,Precision@f[x]]; Abort[];
	];
]
NumberPrecision::PrecisionError= "代入した値の精度は (`1`) です．目標精度(Quantica'MP'dps)より低いため評価を中断しました"
NumberPrecision[x_]:=If[ Precision@x < Quantica`MP`dps, Message[NumberPrecision::PrecisionError, Precision@x]; Abort[]]

Begin["`Private`"]
(*Privateに入れないとNames["Quantica`アスタリスク"]に(局所)変数が表示されてしまうのだ*)

(*QuanticaSetting[dim_, domain_, tau_:1, Verbose_:False] :=  Module[ {},*)


(*todo:
    optionの局所変数化*)
Options[QuanticaSetting] = {tau->1, Verbose->True,a->1}    
QuanticaSetting[dim_, domain_, OptionsPattern[]] :=  Module[ {},
    Which[
        Not[IntegerQ[dim]], 
                Message[QuanticaSetting::dimError,dim],
        TrueQ[dim<0],
                Message[QuanticaSetting::dimError,dim],
        Dimensions[domain]!={2,2},
                Message[QuanticaSetting::domainError,domain],
		(*                
        Length[Select[Flatten[domain],Precision[#] == Infinity &]]!=4,
                Message[QuanticaSetting::precisionError,domain],
        Precision[OptionValue[tau]] != Infinity,
                Message[QuanticaSetting::precisionError,domain],
        *)
        domain[[1]][[1]]>=domain[[1]][[2]],
                Message[QuanticaSetting::QdomainError,domain[[1]]],
        domain[[2]][[1]]>=domain[[2]][[2]],
                Message[QuanticaSetting::PdomainError,domain[[2]]],
        True,
        	NumberPrecision[#] &/@ {dim, domain, OptionValue[tau]};
			Unprotect[Dim, Domain,Tau, X, PSArea, Planck, Hbar];        
            Dim=dim; Domain=domain;X=setX[];
            PSArea=setArea[];Planck=setPlanck[];Tau=OptionValue[tau];Hbar=setHbar[];
            (*Quantica`MP`dps=OptionValue[Dps];*)            
            If[OptionValue[Verbose]==True,
                Print[{"Dim:", Dim,"Domain:",Domain,
                "Tau:",Tau,"dps:",Quantica`MP`dps}]
            ];
           
            Get["Quantica`Util`"];
            Get["Quantica`State`"];
        	Get["Quantica`PositionBase`"];
        	Get["Quantica`HarmonicBase`"];
        	Get["Quantica`Qmap`"];        	
    ];
	Protect[Dim, Domain,Tau, X, PSArea, Planck, Hbar];        	
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
Options[ParallelEigen] = {Sort->True, Vector->True,KernelNum->1}
ParallelEigen[mats_,OptionsPattern[]]:=Module[{i,res,index,num,MAP},
	num=OptionValue[KernelNum];
	If[num>1, LaunchKernels[num]];
	If[Length@Kernels[] > 1, MAP=ParallelMap, MAP=Map];
	If[ OptionValue[Vector],
		res = MAP[ Eigensystem, mats];	
		index = Ordering[ Re[res[[#]][[1]]] ] & /@ Range[Length[mats]];
		For[i=1,i<=Length[mats],i++,
			res[[i]][[1]] = res[[i]][[1]][[ index[[i]] ]];	
			res[[i]][[2]] = res[[i]][[2]][[ index[[i]] ]];			
		],
		res = MAP[Eigenvalues, mats];
		index = Ordering[ Re[res[[#]] ] ] & /@ Range[Length[mats]];
		For[i=1,i<=Length[mats],i++,
			res[[i]] = res[[i]][[ index[[i]] ]]	
		];
	];
	Return[res];
]

QuasiEnergy[evals_]:= Hbar*I*Log[evals]/Tau

(*duplex problem occur using this method
SortIndex1[ref_,evecs_]:= Ordering [ State`Overlap[evecs, ref[[#]] ] ][[-1]] & /@ Range[Length[ref]]
*)

SortIndex[ref_,evecs_]:=Module[{i,j,index, maxindex},
  index = {};
  For[i=1,i<=Length@evecs,i++,
    maxindex = Ordering[State`Overlap[evecs,ref[[i]] ] ] [[-1]];
    j=1;
    While[MemberQ[index,maxindex],
      maxindex = Ordering[State`Overlap[evecs, ref[[i]] ] ] [[-j]];
      j+=1;
    ];
    index = Append[index, maxindex];
  ];
  Return[index];
]

                                        
SortEigen[evals_,evecs_,index_]:=Module[{vals,vecs},
    vals=evals[[index]];
    vecs=evecs[[index]];
    Return[{vals,vecs}]
]

(* Operation *)



(*Private functions*)
setX[] := Module[{q,p,dq,dp}, 
    dq = (Domain[[1]][[2]] - Domain[[1]][[1]])/Dim;
    dp = (Domain[[2]][[2]] - Domain[[2]][[1]])/Dim;
    q = N[Linspace[Domain[[1]][[1]], Domain[[1]][[2]], Dim,False],Quantica`MP`dps]; 
    p = N[Linspace[Domain[[2]][[1]], Domain[[2]][[2]], Dim,False],Quantica`MP`dps];
    Return[{q,p}]
]
setArea[]  := N[(Domain[[1]][[2]] - Domain[[1]][[1]])*(Domain[[2]][[2]] - Domain[[2]][[1]]), Quantica`MP`dps]
setPlanck[]:= N[PSArea/Dim,Quantica`MP`dps]
setHbar[]  := N[Planck/(2*Pi),Quantica`MP`dps]


End[](*end private*)
EndPackage[]

Column[{
"Quantica`"<>
"A Mathematica Package for calculations of quantum systems in arbitary precision. \n" <>
"by Yasutaka Hanada.\n",
"注意事項：\n"<>
"任意精度計算の為に，機械精度の浮動小数点表現(0.5や1.2等)を用いないで下さい．\n" <>
"代わりに，有理数表現(1/2, 2/10)もしくは無理数表現(Sqrt[3])等を用いて下さい．\n" <>
"数値を代入する場合は精度Precision[x]が"<>ToString[Precision[1/2]]<>"となる様にして下さい．\n"<>
"\nMATHEMATICA "<>$Version,
"\nQuantica version 0.6.3 \[Beta]",
"\nTODAY IS "<>DateString[]
}]  
