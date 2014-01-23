(* Mathematica package *)


Begin["Systems`"]
Harmonic::usage="usage Haromic[omega]"
Henon::usage="usage Haromic[a]"
DoubleWell::usage="usage Haromic[a]"
Pendulum::usage="usage Pendulum[k]"
Standard::usage="usage Standard[k]\n
Return functions {T,V}, where T[x]=x^2, V[x]=k*Cos(2*Pi*x)/(4*Pi^2)"
LinearRotor::usage="usage LinearStandrd[k,omega]\n
Return functions {T,V}, where T[x]=omega*x, V[x]=k*Cos(2*Pi*x)/(4*Pi^2)"
LogCoshRotor::usage="usage LogCoshRotor[k,omega,beta]\n"


(*instantation*)
Systems`Harmonic		= Systems`Private`Harmonic
Systems`Henon			= Systems`Private`Henon
Systems`DoubleWell		= Systems`Private`DoubleWell
Systems`Pendulum		= Systems`Private`Standard
Systems`Standard 		= Systems`Private`Standard
Systems`LinearRotor 	= Systems`Private`LinearRotor
Systems`LogCoshRotor 	= Systems`Private`LogCoshRotor
End[]

Begin["Systems`Private`"]
testPrec::precisionError = "input value is `1`, usege infinite precision or > dps "
testPrec[x_]:= If[ Precision[x] < Quantica`MP`dps, Message[testPrec::precisionError,x];Abort[],True]

Harmonic[omega_:1]:= Module[{T,V},
	Map[testPrec,{omega}];
	T[x_]:=x^2/2;
	V[x_]:=omega*x^2/2;
	Return[{T,V}];
]
Henon[a_:1]:= Module[{T,V},
	Map[testPrec,{a}];
	T[x_]:=x^2/2;
	V[x_]:= a*(2*x^2 + x^3/3);
	Return[{T,V}];
]
DoubleWell[a_]:= Module[{T,V},
	Map[testPrec,{a}];
	T[x_]:=x^2/2;
	V[x_]:=(x^2 - a^2)^2;
	Return[{T,V}];
]
	

Standard[k_:(Pi/2)^2*(35/100)^2] := Module[{T,V},
	Map[testPrec, {k}];
	T[x_]:=x^2/2;
	V[x_]:=k*Cos[2*Pi*x]/(4*Pi^2);
	Return[{T,V}];
]


LinearRotor[k_:2,omega_:Sqrt[2]]:= Module[{T,V},
	Map[testPrec, {k,omega}];
	T[x_]:=omega*x;
	V[x_]:=k*Cos[2*Pi*x]/(4*Pi^2);
	Return[{T,V}];
]
LogCoshRotor[k_:2,omega_:Sqrt[2],beta_:10]:= Module[{T,V},
	Map[testPrec, {k,omega.beta}];
	T[x_]:= Log[Cosh[beta*x]]/(2*beta) + omega*x  + x/2;
	V[x_]:=k*Cos[2*Pi*x]/(4*Pi^2);
	Return[{T,V}];	
]


End[]