(* Mathematica Package *)



Begin["BCH`"]
(* Exported symbols added here with SymbolName::usage *)
Quantica`BCH::usage="context of the BCH seriese"
Translated2::usage="Translated2[n] return n-th BCH serise for 2 word (x,y)"
Translated3::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,w)"
Translated3S::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,x)"
BCH2::usage="BCH2[n,A,B] return n-th order BCH matrix"
BCH3::usage="BCH3[n,A,B,C] return n-th order BCH matrix"
BCH2List::usage""
(*instantiation*)
BCH`Translated2     :=  BCH`Private`Translated2
BCH`Translated3     :=  BCH`Private`Translated3
BCH`Translated3S    :=  BCH`Private`Translated3S
BCH`BCH2            :=  BCH`Private`BCH2
BCH`BCH2List        :=  BCH`Private`BCH2List
BCH`BCH3            :=  BCH`Private`BCH3
BCH`BCH3List        :=  BCH`Private`BCH3List

End[]

Begin["BCH`Private`"]

 
p2[n_] := p2[n] = ( F = Table[1/(j-i)!,{i,n+1},{j,n+1}];
          G = Table[1/(j-i)! Product[s[k],{k,i,j-1}],{i,n+1},{j,n+1}]; 
          qthpower = IdentityMatrix[n+1]; 
          FGm1 = F.G - qthpower;
          Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]])
  
Translated2[n_] := (temp = Expand[Product[s[k]^2, {k,n}] p2[n]];
             Sum[term = Apply[List, temp[[i]]]; term[[1]] Apply[StringJoin, Take[term,-n] /. {s[i_]^2->"x",s[i_]^3->"y"}], {i,Length[temp]}])

p3[n_] := p3[n] = (
           F = Table[1/(j-i)! Product[s[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}];
           G = Table[1/(j-i)! Product[s[k,"y"],{k,i,j-1}],{i,n+1},{j,n+1}];
           H = Table[1/(j-i)! Product[s[k,"w"],{k,i,j-1}],{i,n+1},{j,n+1}]; 
           qthpower = IdentityMatrix[n+1]; 
           FGm1 = F.G.H - qthpower;
           Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]]
           )
Translated3[n_] := (temp = p3[n]; Sum[term = Apply[List, temp[[i]]]; term[[1]]* Apply[StringJoin, Take[term,-n] /. s[j_,k_]->k], {i,Length[temp]}])
ps3[n_] := ps3[n] = (
           F = Table[1/(j-i)! Product[s[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}];
           G = Table[1/(j-i)! Product[s[k,"y"],{k,i,j-1}],{i,n+1},{j,n+1}];
           H = Table[1/(j-i)! Product[s[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}]; 
           qthpower = IdentityMatrix[n+1]; 
           FGm1 = F.G.H - qthpower;
           Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]]
           )
Translated3S[n_] := (temp = ps3[n]; Sum[term = Apply[List, temp[[i]]]; term[[1]]* Apply[StringJoin, Take[term,-n] /. s[j_,k_]->k], {i,Length[temp]}])

Ravel[n_Integer, bch_] :=  bch[[n]][[1]]*Fold[Dot, "Iden", StringSplit @@ {bch[[n]][[2]], ""}]

Decode2Matrix[A_,B_,bch_]:=Module[{i,decode, list,mat},
	decode = Table[Ravel[i,bch],{i,1,Length[bch]}];
	SetSharedVariable[decode];
	(*list = ParallelTable[ メモリ食いすぎ*)
	mat = ParallelSum[ (*bottle neck*)
		ReplaceAll[decode[[i]],
				{ "x"->A,
			  	  "y"->B,
			      "Iden"->IdentityMatrix[Length[A]] 
			   }], 
		{i,1,Length[decode]}];
	UnsetShared[decode];
	Return[mat];
]
Decode3Matrix[A_,B_,C_,bch_]:=Module[{i,decode, list,mat},
	decode = Table[Ravel[i,bch],{i,1,Length[bch]}];
	SetSharedVariable[decode];
	list = ParallelTable[
		ReplaceAll[decode[[i]],
				{ "x"->A,
			  	  "y"->B,
			  	  "w"->C,
			      "Iden"->IdentityMatrix[Length[A]] 
			   }], 
		{i,1,Length[decode]}];
	UnsetShared[decode];		
	mat = Apply[Plus, list];
	Return[mat];
] 
Options[BCH2] = {Verbose->True,KernelNum->1}
BCH2[n_Integer,A_,B_,OptionsPattern[]]:=Module[{mat, i,res,num},
	num = OptionValue[KernelNum];
	Which[Length[Kernels[]]==0 && $ProcessorCount > num > 0, LaunchKernels[num],
		  Length[Kernels[]]==0,LaunchKernels[]
	];
	res = Table[
		If[ OptionValue[Verbose], Print[{"2 word BCH series:",i}] ];				
		Decode2Matrix[A,B,Translated2[i]],
		{i,2,n}];
	mat = Fold[Plus, A + B, res];
	Return[mat];
]
Options[BCH2List] = {Verbose->True}
BCH2List[n_Integer,A_,B_,OptionsPattern[]] := Module[{res, list,i},
	res = Table[
			If[ OptionValue[Verbose], Print[{"2 word BCH series:",i}] ];				
			Decode2Matrix[A,B,Translated2[i]],
			{i,2,n}];
	list = FoldList[Plus, A+B,res];
	Return[list];
]

Options[BCH3] = {Verbose->False,KernelNum->1}
BCH3[n_,A_,B_,C_,opts:OptionsPattern[]] := Module[{mat},
    If[A==C,mat = SymmetricBCH3[n,A,B,opts], mat = GeneralBCH3[n,A,B,C,opts] ];
    Return[mat];
]
Options[BCH3List] = {Verbose->False,KernelNum->1}
BCH3List[n_,A_,B_,C_,opts:OptionsPattern[]] := Module[{mat},
    If[A==C,mat = SymmetricBCH3List[n,A,B,opts], mat = GeneralBCH3List[n,A,B,C,opts] ];
    Return[mat];
]
	
Options[SymmetricBCH3] = {Verbose->False,KernelNum->1}
SymmetricBCH3[n_,A_,B_,opt:OptionsPattern[]]:=Module[{mat,num,i,res},
	num = OptionValue[KernelNum];
	Which[Length[Kernels[]]==0 && $ProcessorCount > num > 0, LaunchKernels[num],
		  Length[Kernels[]]==0,LaunchKernels[]
	];
	res = Table[
		If[ OptionValue[Verbose], Print[{"2 word Symmetric BCH series:",i}] ];				
		Decode2Matrix[A,B,Translated3S[i]],
		{i,3,n,2}];
	mat = Fold[Plus, 2*A + B, res];
	Return[mat];	
]
Options[SymmetricBCH3List] = {Verbose->False,KernelNum->1}
SymmetricBCH3List[n_,A_,B_,OptionsPattern[]]:=Module[{mat,num,i,res},
	num = OptionValue[KernelNum];
	Which[Length[Kernels[]]==0 && $ProcessorCount > num > 0, LaunchKernels[num],
		  Length[Kernels[]]==0,LaunchKernels[]
	];
	res = Table[
		If[ OptionValue[Verbose], Print[{"2 word Symmetric BCH series:",i}] ];				
		Decode2Matrix[A,B,Translated3S[i]],
		{i,3,n,2}];
	mat = FoldList[Plus, 2*A + B, res];
	Return[mat];	
]

Options[GeneralBCH3] = {Verbose->False, KernelNum->1}
GeneralBCH3[n_,A_,B_,C_,OptionsPattern[]] := Module[{mat,res,num,i},
	num = OptionValue[KernelNum];
	Which[Length[Kernels[]]==0 && $ProcessorCount > num > 0, LaunchKernels[num],
		  Length[Kernels[]]==0,LaunchKernels[]
	];
	res = Table[
		If[ OptionValue[Verbose], Print[{"3 word BCH series:",i}] ];				
		Decode3Matrix[A,B,C,Translated3[i]],
		{i,2,n}];
	mat = Fold[Plus, A + B + C, res];
	Return[mat];	
	
]
Options[GeneralBCH3List] = {Verbose->False, KernelNum->1}
GeneralBCH3List[n_,A_,B_,C_,OptionsPattern[]] := Module[{mat,res,num,i},
	num = OptionValue[KernelNum];
	Which[Length[Kernels[]]==0 && $ProcessorCount > num > 0, LaunchKernels[num],
		  Length[Kernels[]]==0,LaunchKernels[]
	];
	res = Table[
		If[ OptionValue[Verbose], Print[{"3 word BCH series:",i}] ];				
		Decode3Matrix[A,B,C,Translated3[i]],
		{i,2,n}];
	mat = Fold[Plus, A + B + C, res];
	Return[mat];	
	
]


End[] 