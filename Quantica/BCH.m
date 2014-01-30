(* Mathematica Package *)


BeginPackage["Quantica`BCH`"]
Begin["BCH`"]

(* Exported symbols added here with SymbolName::usage *)
Quantica`BCH::usage="context of the BCH seriese"
Translated2::usage="Translated2[n] return n-th BCH serise for 2 word (x,y)"
Translated3::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,w)"
Translated3S::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,x)"
BCH2::usage="BCH2[n,A,B] return n-th order BCH matrix"
BCH3::usage="BCH3[n,A,B,C] return n-th order BCH matrix"
BCH2List::usage"BCH[n,A,B] n次までのBCH展開による行列表現を与える"

(*instantiation*)
BCH`Translated2     :=  BCH`Private`Translated2
BCH`Translated3     :=  BCH`Private`Translated3
BCH`Translated3S    :=  BCH`Private`Translated3S
BCH`BCH2            :=  BCH`Private`BCH2
BCH`BCH2List        :=  BCH`Private`BCH2List
BCH`BCH3            :=  BCH`Private`BCH3
BCH`BCH3List        :=  BCH`Private`BCH3List

BCH`BCH2Operation			:=	BCH`Private`BCH2Operation
BCH`BCH2HamiltonianTerm		:=	BCH`Private`BCH2HamiltonianTerm
BCH`BCH2Hamiltonian			:=	BCH`Private`BCH2Hamiltonian
BCH`BCH3SOperation			:=	BCH`Private`BCH3SOperation
BCH`BCH3SHamiltonianTerm	:=	BCH`Private`BCH3SHamiltonianTerm
BCH`BCH3SHamiltonian		:=	BCH`Private`BCH3SHamiltonian

BCH`ClassicalTerm	:=	BCH`Private`ClassicalTerm

End[]

Begin["BCH`Private`"]

Psi[x_] := Global`\[Psi][x];
HBar = Global`\[HBar];
Tau = Global`\[Tau];
S = Global`S;
q = Global`q;
p = Global`p;
Protect[Global`\[Psi],Global`\[HBar], Global`S, Global`q, Global`p];
 
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


OrderFunc[f_] := Length[CoefficientList[f[x], x]] - 1
RevalToOperator[word_, opT_,opV_] := Module[{oplist,res},
	oplist=StringSplit[word[[2]], ""] /. {"x" -> opT, "y" -> opV};
	res = word[[1]]*Fold[#2[#1] &, Psi[q], Reverse[oplist]];
	Return[res];
]
BCH2Operation[n_Integer,T_,V_]:=BCH2Operation[n,T,V]=Module[{bch,len, order, opT, opV,res,psilist, replacelist},	
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

Options[BCH2HamiltonianTerm] = {S->Symbol["S"]}
BCH2HamiltonianTerm[1,T_,V_,OptionsPattern[]] = T[p] + V[q] 
BCH2HamiltonianTerm[n_Integer, T_,V_,OptionsPattern[]]:=Module[{op, list},
	{op,list} = BCH2Operation[n,T,V] // Expand;
	Return[op /. list /.S->OptionValue[S] ]
]
Options[BCH2Hamiltonian] = {S->Symbol["S"]}
BCH2Hamiltonian[n_Integer,T_,V_,opts:OptionsPattern[]]:=Fold[Plus, BCH2HamiltonianTerm[1,T,V,opts], Table[BCH2HamiltonianTerm[i,T,V,opts], {i,2,n}]] /. S->OptionValue[S]

BCH3SOperation[n_Integer,T_,V_]:=BCH3SOperation[n,T,V]=
	Module[{order, opT, opV,res,len,psilist, replacelist,bch},
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

Options[BCH3SHamiltonianTerm] = {TVT->False,S->Symbol["S"]}
BCH3SHamiltonianTerm[1,T_,V_,OptionsPattern[]] = If[OptionValue[TVT], 2*T[p] + V[q], T[p] + 2*V[q] ]; 
BCH3SHamiltonianTerm[n_Integer, T_,V_,OptionsPattern[]] := Module[{op, res, list,factor},
	If[ ( OptionValue[TVT] && Mod[n+1,4]==0 ),factor=-1,factor=1]; 
	If[Mod[n,2]!=0,
		{op,list} = BCH3SOperation[n,T,V] // Expand;
		res = op /. list /. S->OptionValue[S],
		res=0
	];
	Return[factor*res]
]
Options[BCH3SHamiltonian] = {TVT->False, S->Symbol["S"]}
BCH3SHamiltonian[n_Integer,T_,V_,opts:OptionsPattern[]]:=Module[{},
	BCH3SHamiltonian[n,T,V] = Fold[Plus, BCH3SHamiltonianTerm[1,T,V,opts], Table[BCH3SHamiltonianTerm[i,T,V,opts], {i,3,n,2}]] /.S->OptionValue[S]
] 


ClassicalTerm[H_]:=Coefficient[H, HBar, 0]

End[]
EndPackage[] 