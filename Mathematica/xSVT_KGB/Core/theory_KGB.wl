(* ::Package:: *)

(* ::Title:: *)
(*Theory KGB*)


(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*Scalar Field*)


DefTensor[scalarcov[], M4, PrintAs -> "\[Phi]"]
DefTensorPerturbation[pertscalarcov[LI[order]], scalarcov[], M4, PrintAs -> "\[Delta]\[Phi]"]


DefTensor[Xcov[], M4, PrintAs -> "X"]
DefTensorPerturbation[pertXcov[LI[order]], Xcov[], M4, PrintAs -> "\[Delta]X"]


DefTensorSVT[scalar[], M1, PrintAs -> "\[Phi]", BackgroundQ->True]
DefTensorSVT[pertscalar[LI[order]], {M1, M3}, PrintAs -> "\[Delta]\[Phi]", ScalarPertQ->True]
DefTensorSVT[pertV[LI[order]], {M1, M3}, PrintAs -> "\!\(\*SubscriptBox[\(v\), \(X\)]\)", ScalarPertQ->True]
DefTensorSVT[pertP[LI[order]], {M1, M3}, PrintAs -> "\[Pi]", ScalarPertQ->True]


DefTensorSVT[X[], M1, BackgroundQ->True]
DefTensorSVT[pertX[LI[order]], {M1, M3}, PrintAs -> "\[Delta]X", ScalarPertQ->True]


(* ::Subsection::Closed:: *)
(*Alphas*)


DefTensorSVT[currentS[], M1, PrintAs -> "\!\(\*SubscriptBox[\(J\), \(\[Phi]\)]\)", BackgroundQ->True]
DefTensorSVT[shiftS[], M1, PrintAs -> "\!\(\*SubscriptBox[\(S\), \(\[Phi]\)]\)", BackgroundQ->True]


DefTensorSVT[densityS[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(S\)]\)", BackgroundQ->True]
DefTensorSVT[pressureS[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[ScriptCapitalP]\), \(S\)]\)", BackgroundQ->True]


DefTensorSVT[alphaK[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(K\)]\)", BackgroundQ->True]
DefTensorSVT[alphaB[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(B\)]\)", BackgroundQ->True]


DefTensorSVT[alphaKK[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(KK\)]\)", BackgroundQ->True]
DefTensorSVT[alphaBB[], M1, PrintAs -> "\!\(\*SubscriptBox[\(\[Alpha]\), \(BB\)]\)", BackgroundQ->True]


DefTensorSVT[wS[], M1, PrintAs -> "\!\(\*SubscriptBox[\(w\), \(S\)]\)", BackgroundQ->True]
DefTensorSVT[kinD[], M1, PrintAs -> "D", BackgroundQ->True]
DefTensorSVT[cs2[], M1, PrintAs -> "\!\(\*SuperscriptBox[SubscriptBox[\(c\), \(S\)], \(2\)]\)", BackgroundQ->True]
DefTensorSVT[Qs[], M1, PrintAs -> "\!\(\*SubscriptBox[\(Q\), \(S\)]\)", BackgroundQ->True]
DefTensorSVT[kB[], M1, PrintAs -> "\!\(\*SubscriptBox[\(k\), \(B\)]\)", BackgroundQ->True]


(* ::Subsection::Closed:: *)
(*Scalar functions*)


DefScalarFunction[G2fun, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(2\)]\)"]
DefScalarFunction[G3fun, PrintAs -> "\!\(\*SubscriptBox[\(G\), \(3\)]\)"]


(* ::Subsection::Closed:: *)
(*Sources*)


DefTensor[source11[LI[order]], {M1, M3}, PrintAs -> "\!\(\*SubscriptBox[\(S\), \(11\)]\)"]


(* ::Section::Closed:: *)
(*Expansion Rules*)


(****   Scalar Field   ****)

scalarrules = Flatten[{
	scalarcov[] :> scalar[],
	pertscalarcov[LI[order_]] :> pertscalar[LI[order]],
	Xcov[]:>X[],
	pertXcov[LI[order_]]:>pertX[LI[order]]}];
If[fill,
	$SVTDecompositionRules[[1]] = Union[$SVTDecompositionRules[[1]], scalarrules];
]

Print[Column[{"Scalar field Decomposition", ScreenDollarIndices[scalarrules]}]]
Clear[scalarrules]


(****   X to Scalar Field   ****)
XcovToScalarcov[expr_]:= expr //.MakeRule[{Xcov[], Evaluate[-1/2 metricg[\[Mu], \[Nu]] CD[-\[Mu]]@scalarcov[] CD[-\[Nu]]@scalarcov[]]}];
ScalarcovToXcov[expr_]:= expr //.Gfun_[scalarcov[],x_]:>Gfun[scalarcov[],Xcov[]];


(* ::Section::Closed:: *)
(*To alphas*)


(* ::Subsection:: *)
(*Definitions*)


EQalphaK = - hubbleC[]^2/scale[]^2 alphaK[]*Mpl^2 +
	(6*hubbleC[]*primescalar[]^3*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^4 +
	(primescalar[]^2*Derivative[0, 1][G2fun][scalar[], X[]])/scale[]^2 +
	(3*hubbleC[]*primescalar[]^5*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^6  +
	(primescalar[]^4*Derivative[0, 2][G2fun][scalar[], X[]])/scale[]^4  -
	(2*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]])/scale[]^2  -
	(primescalar[]^4*Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^4  // ToCanonical;


EQalphaB = - alphaB[] hubbleC[]/scale[] Mpl^2 +
	(primescalar[]^3*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^3// ToCanonical;


EQalphaKK = - hubbleC[]^2*Mpl^2*alphaKK[]/scale[]^2  -
	(primescalar[]^4*Derivative[0, 2][G2fun][scalar[], X[]])/(4*scale[]^4) -
	(hubbleC[]*primescalar[]^7*Derivative[0, 3][G3fun][scalar[], X[]])/(4*scale[]^8)  -
	(primescalar[]^6*Derivative[0, 3][G2fun][scalar[], X[]])/(12*scale[]^6)  +
	(primescalar[]^4*Derivative[1, 1][G3fun][scalar[], X[]])/(3*scale[]^4)  +
	(primescalar[]^6*Derivative[1, 2][G3fun][scalar[], X[]])/(12*scale[]^6)// ToCanonical;


EQalphaBB = - hubbleC[]*Mpl^2*alphaBB[]/scale[] -
	(primescalar[]^5*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^5// ToCanonical;


EQdensityS = -densityS[] -
	G2fun[scalar[], X[]]+
	(3*hubbleC[]*primescalar[]^3*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^4  +
	(primescalar[]^2*Derivative[0, 1][G2fun][scalar[], X[]])/scale[]^2  -
	(primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]])/scale[]^2// ToCanonical;


EQpressureS = - pressureS[] +
	G2fun[scalar[], X[]] -
	(pprimescalar[]*primescalar[]^2*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^4 +
	(hubbleC[]*primescalar[]^3*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^4  -
	(primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]])/scale[]^2// ToCanonical;


(* ::Subsection:: *)
(*Functions*)


(****   NoG's   ****)
NoG2[expr_]:=expr//.Derivative[__][G2fun][__]:>0//.G2fun[__]:>0
NoG3[expr_]:=expr//.Derivative[__][G3fun][__]:>0//.G3fun[__]:>0


PertScalarToPertP[expr_] := Module[{tmp}, tmp = expr;
	tmp = tmp //.pertscalar[LI[2]] :> primescalar[] pertP[LI[2]] + pprimescalar[] pertP[LI[1]]^2 // Expand;
	tmp = tmp //.pertscalar[LI[1]] :> primescalar[] pertP[LI[1]] // Expand;
	tmp
]


PertPToPertScalar[expr_] := Module[{tmp}, tmp = expr;
	tmp = tmp //.pertP[LI[2]] :> pertscalar[LI[2]]/primescalar[] - pprimescalar[]/primescalar[]^3 pertscalar[LI[1]]^2 // Expand;
	tmp = tmp //.pertP[LI[1]] :> pertscalar[LI[1]]/primescalar[] // Expand;
	tmp
]


subSV1 = -primescalar[]/scale[] PD[-i]@pertV[LI[2]]+
	2/primescalar[] timevec[a] PD[-a]@pertscalar[LI[1]] PD[-i]@pertscalar[LI[1]]-
	2 pertpsi[LI[1]] PD[-i]@pertscalar[LI[1]];
subSV2 = PD[-j]@subSV1 // SVTExpand // Symmetrize // SVTExpand;
subSV4 = PD[-j]@PD[-k]@PD[-l]@subSV1 // SVTExpand // Symmetrize // SVTExpand;


PertScalarToPertV[expr_] := Module[{tmp}, tmp = expr;
	tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@PD[-j]@PD[-k]@PD[-l]@pertscalar[LI[2]]], Evaluate[subSV4]}] // Expand;
	tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@PD[-j]@pertscalar[LI[2]]], Evaluate[subSV2]}] // Expand;
	(*tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@pertscalar[LI[2]]], Evaluate[subSV1]}] // Expand;*)
	tmp = tmp //.pertscalar[LI[1]] :> -primescalar[]/scale[] pertV[LI[1]] // Expand;
	tmp
]


subPV1 = -PD[-i]@pertV[LI[2]]/scale[]+
	2 timevec[a] PD[-a]@pertP[LI[1]] PD[-i]@pertP[LI[1]]-
	2 pertpsi[LI[1]] PD[-i]@pertP[LI[1]];
subPV2 = PD[-j]@subPV1 // SVTExpand // Symmetrize // SVTExpand;
subPV4 = PD[-j]@PD[-k]@PD[-l]@subPV1 // SVTExpand // Symmetrize // SVTExpand;


PertPToPertV[expr_] := Module[{tmp}, tmp = expr;
	tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@PD[-j]@PD[-k]@PD[-l]@pertscalar[LI[2]]], Evaluate[subPV4]}] // Expand;
	tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@PD[-j]@pertscalar[LI[2]]], Evaluate[subPV2]}] // Expand;
	(*tmp = xSVTUtilities`FirstS[tmp] //.MakeRule[{Evaluate[PD[-i]@pertscalar[LI[2]]], Evaluate[subPV1]}] // Expand;*)
	tmp = tmp //.pertP[LI[1]] :> -pertV[LI[1]]/scale[] // Expand;
	tmp
]


(* ::Section::Closed:: *)
(*Tex Corrections*)


Tex[pertV] ^= "v_\\textrm{X}";


Tex[densityS] ^= "\\mathcal{E}_\\textrm{S}";
Tex[pressureS] ^= "\\mathcal{P}_\\textrm{S}";


Tex[alphaK] ^= "\\alpha_\\textrm{K}";
Tex[alphaB] ^= "\\alpha_\\textrm{B}";


Tex[alphaKK] ^= "\\alpha_\\textrm{KK}";
Tex[alphaBB] ^= "\\alpha_\\textrm{BB}";


Tex[cs2] ^= "c_\\textrm{S}^2";
Tex[Qs] ^= "Q_\\textrm{S}";
