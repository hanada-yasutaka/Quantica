(* Mathematica Package *)

BeginPackage["Quantica`BCH`"]
(* Exported symbols added here with SymbolName::usage *)
Quantica`BCH::usage="context of the BCH seriese"
Translated2::usage="Translated2[n] return n-th BCH serise for 2 word (x,y)"
Translated3::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,w)"
Translated3S::usage="Translated2[n] return n-th BCH serise for 3 word (x,y,x)"
BCH2::usage="BCH2[n,A,B] return n-th order BCH matrix"
BCH3::usage="BCH3[n,A,B,C] return n-th order BCH matrix"

Begin["`Private`"] (* Begin Private Context *)
dim = Quantica`Dim
dps = Quantica`MP`dps
 
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

Options[BCH2] = {Verbose->False}
BCH2[n_,A_,B_,OptionsPattern[]] := Module[{i,j,l,word,term,bch,mat},
    mat = A + B;
    For [l=1,l<n,l++;
        bch = Translated2[l];
        If[OptionValue[Verbose],Print[{"2 word BCH series:",l}]];
        For [i=0,i<Length[bch],++i;
            word=StringSplit[bch[[i]][[2]], ""];
            term = DiagonalMatrix[Table[bch[[i]][[1]],{dim}]];
            For[j=Length[word]+1, j>1, j--;
                Which[ word[[j]]=="x", term=Dot[A,term],
                       word[[j]]=="y", term=Dot[B,term]
                ];
            ];
        mat = mat + term;
        ];
    ];
    Return[mat];
]
Options[BCH3] = {Verbose->False}
BCH3[n_,A_,B_,C_,OptionsPattern[]] := Module[{mat},
    mat = A + B + C;
    If[A==C,mat = SymmetricBCH3[n,A,B,OptionValue[Verbose]], mat = GeneralBCH3[n,A,B,C,OptionValue[Verbose]] ];
    Return[mat];
]
SymmetricBCH3[n_,A_,B_,verbose_:False]:=Module[{i,j,l,mat,word,term,bch},
    mat = 2*A + B;
    For [l=1,l<n,l+=1;
        bch = Translated3S[l];
        If[verbose,Print[{"2 word Symmetric BCH series:",l}]];
        For [i=0,i<Length[bch],++i;
            word=StringSplit[bch[[i]][[2]], ""];
            term = DiagonalMatrix[Table[bch[[i]][[1]],{dim}]];
            For[j=Length[word]+1, j>1, j--;
                Which[ word[[j]]=="x", term=Dot[A,term],
                       word[[j]]=="y", term=Dot[B,term]
                ];
            ];
        mat = mat + term;
        ];
    ];
    Return[mat];
]

GeneralBCH3[n_,A_,B_,C_,verbose_:False] := Module[{i,j,l,mat,word,term,bch},
    mat = A + B + C;
    For [l=1,l<n,l+=1;
        bch = Translated3[l];
        If[verbose, Print[{"3 word BCH series:",l}]];
        For [i=0,i<Length[bch],++i;
            word=StringSplit[bch[[i]][[2]], ""];
            term = DiagonalMatrix[Table[bch[[i]][[1]],{dim}]];
            For[j=Length[word]+1, j>1, j--;
                Which[ word[[j]]=="x", term=Dot[A,term],
                       word[[j]]=="y", term=Dot[B,term],
                       word[[j]]=="w", term=Dot[C,term]
                ];
            ];
        mat = mat + term;
        ];
    ];
    Return[mat];
]

End[] (* End Private Context *)

EndPackage[]