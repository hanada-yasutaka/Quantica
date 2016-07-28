(* Mathematica package *)

BeginPackage["Quantica`Systems`"]
Begin["Systems`"]
Harmonic::usage="usage Haromic[omega]"
Henon::usage="usage Haromic[omega]"
DoubleWell::usage="usage DoubleWell[a]"
Pendulum::usage="usage Pendulum[k]"
Standard::usage="usage Standard[k]"
LinearRotor::usage="usage LinearStandrd[k,omega]"
LogCoshRotor::usage="usage LogCoshRotor[k,omega,beta]"
JeremyNormal::usage"JeremyNormal[a,b]"

(*instantation*)
Systems`Harmonic		= Systems`Private`Harmonic
Systems`Henon			= Systems`Private`Henon
Systems`DoubleWell		= Systems`Private`DoubleWell
Systems`Pendulum		= Systems`Private`Pendulum
Systems`Standard 		= Systems`Private`Standard
Systems`LinearRotor 	= Systems`Private`LinearRotor
Systems`LogCoshRotor 	= Systems`Private`LogCoshRotor
Systems`JeremyNormal 	= Systems`Private`JeremyNormal
Systems`JeremyNormal2 	= Systems`Private`JeremyNormal2


(*a1->1,a2->-55/100,b->5/100,phi->0}*)


End[]

Begin["Systems`Private`"]
testPrec::precisionError = "input value is `1`, usege infinite precision or > dps "
testPrec[x_]:= If[ Precision[x] < Quantica`MP`dps, Message[testPrec::precisionError,x];Abort[],True]

Harmonic[omega_]:= Module[{T,V},
	Map[testPrec,{omega}];
	T=(#^2/2) &;
	V=(omega*#^2/2) &;
	Return[{T,V}];
]

Henon[eps_]:= Module[{T,V},
	Map[testPrec,{eps}];
	T=(#^2/2) &;
	V=eps*(2*#^2 + #^3/3) &;
	Return[{T,V}];
]

DoubleWell[a_]:= Module[{T,V},
	Map[testPrec,{a}];
	T= (#^2/2) &;
	V=(#^2 - a^2)^2 &;
	Return[{T,V}];
]

Pendulum[k_]:=Module[{T,V},
	Map[testPrec, {k}];
	T = (#^2/2) &;
	V = (k*Cos[#]) &;
	Return[{T,V}]
]

Standard[k_:(Pi/2)^2*(35/100)^2] := Module[{T,V},
	Map[testPrec, {k}];
	T=(#^2/2) &;
	V:=(k*Cos[2*Pi*#]/(4*Pi^2)) &;
	Return[{T,V}];
]

LinearRotor[k_:2,omega_:Sqrt[2]]:= Module[{T,V},
	Map[testPrec, {k,omega}];
	T=(omega*#) &;
	V=(k*Cos[2*Pi*#]/(4*Pi^2)) &;
	Return[{T,V}];
]
LogCoshRotor[k_:2,omega_:Sqrt[2],beta_:10]:= Module[{T,V},
	Map[testPrec, {k,omega.beta}];
	T= (Log[Cosh[beta*#]]/(2*beta) + omega*#  + #/2) &;
	V= (k*Cos[2*Pi*#]/(4*Pi^2)) &;
	Return[{T,V}];	
]
(*a1->1,a2->-55/100,b->5/100,phi->0}*)
Options[JeremyNormal] = {a1->1, a2->-55/100,b->5/100, phi->0}
JeremyNormal[OptionsPattern[]]:=Module[{A1,A2,B,Phi,T,V},
	A1 = OptionValue[a1];
	A2 = OptionValue[a2];
	B = OptionValue[b];
	Phi = OptionValue[phi];
	Map[testPrec, {A1,A2,B,Phi}];				
	T=((A1/2)*(Cos[#1]^2 + Cos[#2]^2) + A2*(Cos[#1]^2+Cos[#2]^2)^2 )&;
	V=(Abs[B]*( (Cos[#2]^4 + Cos[#1]^4 - 6*Cos[#2]^2*Cos[#1]^2)*Cos[Phi] -4*(Cos[#2]^3*Cos[#1] - Cos[#2]*Cos[#1]^3)*Sin[Phi] ))&; 
	Return[{T,V}];
]

Options[JeremyNormal2] = {r1 -> 4, r2 -> 11, I1->5/100, I2 -> 45/100, \[Phi]1 -> 7/10*Pi, \[Phi]2 -> 15/10*Pi, b1 -> 3/1000,b2 -> 5/10000}
JeremyNormal2[q_,p_,OptionsPattern[]]:=Module[{H0,J,Theta,V,H,R1,R2,J1,J2,B1,B2,phi1,phi2},

	R1 = OptionValue[r1];
	R2 = OptionValue[r2];
	J1 = OptionValue[I1];
	J2 = OptionValue[I2];
	B1 = OptionValue[b1];
	B2 = OptionValue[b2];
	phi1 = OptionValue[\[Phi]1];
	phi2 = OptionValue[\[Phi]2];

	J[x_, y_] = (Cos[x]^2 + Cos[y]^2)/2;
	Theta[x_, y_] = ArcTan[Cos[x]^2/Cos[y]^2];	

	H0[J_] :=  J^3/3 - (J1 + J2)*J^2/2 + J1*J2*J;
	V[J_,\[Theta]_]:= B1*(J - J2)*(2*J)^(R1/2)*Cos[R1*\[Theta] + phi1] + B2*(J - J1)*(2*J)^(R2/2)*Cos[R2*\[Theta] + phi2];
	H = H0[J[q,p]] + V[J[q,p],Theta[q,p]]
]

End[]
EndPackage[] 