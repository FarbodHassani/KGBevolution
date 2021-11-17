(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[1]] - 3*hubbleC[]*pertdensity[LI[1]]*pressure[] + 
 density[]*timevec[h$2306318]*PD[-h$2306318][pertdensity[LI[1]]] - 
 3*density[]*timevec[h$2306318]*PD[-h$2306318][pertphi[LI[1]]] - 
 3*pressure[]*timevec[h$2306318]*PD[-h$2306318][pertphi[LI[1]]] - 
 density[]*metric\[Delta][p$2306319, p$2306320]*
  PD[-p$2306320][PD[-p$2306319][pertB[LI[1]]]] - 
 metric\[Delta][p$2306319, p$2306320]*pressure[]*
  PD[-p$2306320][PD[-p$2306319][pertB[LI[1]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$2306319, p$2306320]*
   PD[-p$2306320][PD[-p$2306319][pertshear[LI[1]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$2306319, p$2306320]*pressure[]*
   PD[-p$2306320][PD[-p$2306319][pertshear[LI[1]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$2306319, p$2306320]*
   PD[-p$2306320][PD[-p$2306319][pertvelocity[LI[1]]]])/scale[] + 
 (metric\[Delta][p$2306319, p$2306320]*pressure[]*
   PD[-p$2306320][PD[-p$2306319][pertvelocity[LI[1]]]])/scale[] + 
 density[]*metric\[Delta][p$2306319, p$2306320]*timevec[h$2306318]*
  PD[-p$2306320][PD[-p$2306319][PD[-h$2306318][pertE[LI[1]]]]] + 
 metric\[Delta][p$2306319, p$2306320]*pressure[]*timevec[h$2306318]*
  PD[-p$2306320][PD[-p$2306319][PD[-h$2306318][pertE[LI[1]]]]]
