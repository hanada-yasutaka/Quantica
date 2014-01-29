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
(*
q = Quantica`q
p = Quantica`p
*)

OrderFunc[f_] := Length[CoefficientList[f[x], x]] - 1

RevalToOperator[word_, opT_,opV_] := Module[{oplist,res, Psi,q,n},
	Psi = Quantica`\[Psi];
	q = Quantica`q;
	oplist=StringSplit[word[[2]], ""] /. {"x" -> opT, "y" -> opV};
	res = word[[1]]*Fold[#2[#1] &, Psi[q], Reverse[oplist]];
	Return[res];
]
BCH2Operation[n_,T_,V_]:=BCH2Operation[n,T,V]=Module[
	{bch,len,Psi,HBar,Tau, order, opT, opV,res,psilist, replacelist},
	Clear[Quantica`q];Clear[Quantica`p];Clear[Quantica`\[HBar]];Clear[Quantica`\[Psi]];
	Clear[Quantica`S];Clear[Quantica`\[Tau]];
	Psi = Quantica`\[Psi];		
	HBar = Quantica`\[HBar];
	Tau = Quantica`\[Tau];

	order=OrderFunc[T];
	
	opT[f_, x_:q] := S*(T[p] /. p -> HBar/I) * Fold[D, f, Table[x, {i, 1, order}]];
	opV[f_, x_:q] := S*V[x]*f;
	bch = Translated2[n];	
	res = Sum[RevalToOperator[bch[[i]],opT,opV], {i,1,Length[bch] } ]/S;

	len = 2*StringLength @ bch[[1]][[2]];
	psilist = FoldList[D, Psi[q], Table[q, {i, 1, len}]];
	replacelist = Table[psilist[[i]] -> (I/HBar*p)^(i - 1), {i, 1, Length[psilist]}];
	Return[{res//Expand,replacelist}]
]
BCH2HamiltonianTerm[1,T_,V_] = T[p] + V[q] 
BCH2HamiltonianTerm[n_, T_,V_]:=BCH2HamiltonianTerm[n,T,V]=Module[{bch, op, list, HBar,Tau},
	bch = Translated2[n];
	{op,list} = BCH2Operation[n,T,V] // Expand;
(*	HBar = Quantica`\[HBar];*)
(*	Tau = Quantica`\[Tau];*)	
	Return[op /. list ]
]
BCH2Hamiltonian[n_,T_,V_]:=BCH2Hamiltonian[n,T,V]=Fold[Plus, T[p] + V[q], Table[BCH2HamiltonianTerm[i,T,V], {i,2,n}]] /. S->(-I*\[Tau]/\[HBar])

BCH3SOperation[n_,T_,V_]:=BCH3SOperation[n,T,V]=
	Module[{Psi,HBar,Tau, order, opT, opV,res,len,psilist, replacelist,bch},
	Clear[Quantica`q];Clear[Quantica`p];Clear[Quantica`\[HBar]];Clear[Quantica`\[Psi]];
	Clear[Quantica`S];Clear[Quantica`\[Tau]];
	Psi = Quantica`\[Psi];		
	HBar = Quantica`\[HBar];
	Tau = Quantica`\[Tau];
	order=OrderFunc[T];
	opT[f_, x_:q] := S*(T[p] /. p -> HBar/I) * Fold[D, f, Table[x, {i, 1, order}]];
	opV[f_, x_:q] := S*V[x]*f;
	bch = Translated3S[n];
	If[Length[bch]!=0,
		res = Sum[RevalToOperator[bch[[i]],opT,opV], {i,1,Length[bch] } ]/S;
		len = 2*StringLength @ bch[[1]][[2]];
		psilist = FoldList[D, Psi[q], Table[q, {i, 1, len}]];
		replacelist = Table[psilist[[i]] -> (I/HBar*p)^(i - 1), {i, 1, Length[psilist]}];
		Return[{res//Expand,replacelist}],
		Return[{0,0}]
	]
]
BCH3SHamiltonianTerm[1,T_,V_] = T[p] + V[q]; 
BCH3SHamiltonianTerm[n_, T_,V_]:=BCH3SHamiltonianTerm[n,T,V]=Module[{bch,op, res, list, HBar,Tau},
	bch = Translated3S[n];
	If[Length[bch]!=0,
		{op,list} = BCH3SOperation[bch,T,V] // Expand;
		HBar = Quantica`\[HBar];
		Tau = Quantica`\[Tau];
		res = op /. list,
		res=0
	];
	Return[res]
]
BCH3SHamiltonian[n_,T_,V_,order_:"VTV"]:=Module[{H0},
	If[ order=="VTV", H0 = T[p] + 2*V[q], H0 = 2*T[p] + V[q] ];
	BCH3SHamiltonian[n,T,V]= Fold[Plus, H0,Table[BCH3SHamiltonianTerm[i,T,V], {i,3,n,2}]]
] 
(*/.S->(I*\[Tau]/Quantica`\[HBar])*)

ClassicalTerm[H_]:=Coefficient[H, Quantica`\[HBar], 0]



End[]

Begin["BCH`Private`"]

 
p2[n_] := p2[n] = Module[{F,G,i,j,k,qthpower,FGm1,q},( F = Table[1/(j-i)!,{i,n+1},{j,n+1}];
          G = Table[1/(j-i)! Product[ss[k],{k,i,j-1}],{i,n+1},{j,n+1}]; 
          qthpower = IdentityMatrix[n+1]; 
          FGm1 = F.G - qthpower;
          Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]])
]
  
Translated2[n_] := Module[{tmp,temp,term,k,i},(temp = Expand[Product[ss[k]^2, {k,n}] p2[n]];
             Sum[term = Apply[List, temp[[i]]]; term[[1]] Apply[StringJoin, Take[term,-n] /. {ss[i_]^2->"x",ss[i_]^3->"y"}], {i,Length[temp]}])
]

p3[n_] := p3[n] = Module[{F,G,H,i,j,k,qthpower,FGm1,q},(
           F = Table[1/(j-i)! Product[sss[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}];
           G = Table[1/(j-i)! Product[sss[k,"y"],{k,i,j-1}],{i,n+1},{j,n+1}];
           H = Table[1/(j-i)! Product[sss[k,"w"],{k,i,j-1}],{i,n+1},{j,n+1}]; 
           qthpower = IdentityMatrix[n+1]; 
           FGm1 = F.G.H - qthpower;
           Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]]
           )
]
Translated3[n_] := Module[{temp, term,i},
	(temp = p3[n]; Sum[term = Apply[List, temp[[i]]]; term[[1]]* Apply[StringJoin, Take[term,-n] /. sss[j_,k_]->k], {i,Length[temp]}])
]
ps3[n_] := ps3[n] = Module[{F,G,H,i,j,k,qthpower,FGm1,q},(
           F = Table[1/(j-i)! Product[ssss[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}];
           G = Table[1/(j-i)! Product[ssss[k,"y"],{k,i,j-1}],{i,n+1},{j,n+1}];
           H = Table[1/(j-i)! Product[ssss[k,"x"],{k,i,j-1}],{i,n+1},{j,n+1}]; 
           qthpower = IdentityMatrix[n+1]; 
           FGm1 = F.G.H - qthpower;
           Expand[ -Sum[qthpower=qthpower.FGm1; (-1)^q / q qthpower,{q,n}][[1,n+1]]]
           )
]
Translated3S[n_] := Module[{temp, term,i},
	(temp = ps3[n]; Sum[term = Apply[List, temp[[i]]]; term[[1]]* Apply[StringJoin, Take[term,-n] /. ssss[j_,k_]->k], {i,Length[temp]}])
]

RavelToMatrix[n_Integer, bch_] :=bch[[n]][[1]]*Fold[Dot, "Iden", StringSplit @@ {bch[[n]][[2]], ""}]

Decode2Matrix[A_,B_,bch_]:=Decode2Matrix[A,B,bch]=Module[{i,decode, list,mat},
	decode = Table[RavelToMatrix[i,bch],{i,1,Length[bch]}];
	SetSharedVariable[decode];
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
Decode3Matrix[A_,B_,C_,bch_]:=Decode3Matrix[A,B,C,bch]=Module[{i,decode, list,mat},
	decode = Table[RavelToMatrix[i,bch],{i,1,Length[bch]}];
	SetSharedVariable[decode];
	mat = ParallelSum[ (*bottle neck*)	
		ReplaceAll[decode[[i]],
				{ "x"->A,
			  	  "y"->B,
			  	  "w"->C,
			      "Iden"->IdentityMatrix[Length[A]] 
			   }], 
		{i,1,Length[decode]}];
	UnsetShared[decode];
	Return[mat];
] 
Options[BCH2] = {Verbose->True,KernelNum->1}
BCH2[n_Integer,A_,B_,OptionsPattern[]]:=BCH2[n,A,B]=Module[{mat, i,res,num},
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
BCH2List[n_Integer,A_,B_,OptionsPattern[]] :=BCH2List[n,A,B]=Module[{res, list,i},
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
SymmetricBCH3[n_,A_,B_,opt:OptionsPattern[]]:=SymmetricBCH3[n,A,B]=Module[{mat,num,i,res},
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
SymmetricBCH3List[n_,A_,B_,OptionsPattern[]]:=SymmetricBCH3List[n,A,B]=Module[{mat,num,i,res},
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
GeneralBCH3[n_,A_,B_,C_,OptionsPattern[]] := GeneralBCH3[n,A,B,C]=Module[{mat,res,num,i},
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
GeneralBCH3List[n_,A_,B_,C_,OptionsPattern[]] := GeneralBCH3List[n,A,B,C]=Module[{mat,res,num,i},
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