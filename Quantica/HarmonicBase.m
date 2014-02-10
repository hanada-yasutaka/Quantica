(* Mathematica Package *)

Begin["HarmonicBase`"]
(* Exported symbols added here with SymbolName::usage *)
HarmonicBase::usage="調和振動子の基底で演算子のの行列表現を評価するcontextです"
HarmonicBase`Q::usage="return Potential term <n'|q|n>の位置表現を返します"
HarmonicBase`P::usage="return Kinetic term <n'|p|n>の位置表現を返します"

HarmonicBase`Q	:=	HarmonicBase`Private`QMatrix  
HarmonicBase`P	:=	HarmonicBase`Private`PMatrix
HarmonicBase`MatrixRep	:=	HarmonicBase`Private`MatrixRep
End[]

Begin["HarmonicBase`Private`"] (* Begin Private Context *)
hbar=Quantica`Hbar
dim =Quantica`Dim
dps =Quantica`MP`dps 
QMatrix[] := Module[{m},
	m=N[Sqrt[Range[dim]],dps];
	Return[N[Sqrt[hbar/2],dps]*(DiagonalMatrix[m,-1,dim] + DiagonalMatrix[m,1,dim])]
]

PMatrix[] := Module[{m},
	m=N[Sqrt[Range[dim]],dps];
	Return[N[Sqrt[hbar/2],dps]*(DiagonalMatrix[m,-1,dim] - DiagonalMatrix[m,1,dim])*I]
]
MatrixRep[f_,x_]:=Module[{coef,list,i,res},
	coef = CoefficientList[f[x], x];
	list = Table[coef[[i]]*Fold[Dot, "Iden", ConstantArray[x, i - 1]], {i, 1, Length[coef]}];
	res = Fold[Plus, 0, list] /."Iden"->IdentityMatrix[dim];
	Return[res];
]

End[] (* End Private Context *)

