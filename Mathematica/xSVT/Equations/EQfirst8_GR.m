(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[1]] - 3*hubbleC[]*pertdensity[LI[1]]*pressure[] + 
 density[]*timevec[h$70046]*PD[-h$70046][pertdensity[LI[1]]] - 
 3*density[]*timevec[h$70046]*PD[-h$70046][pertphi[LI[1]]] - 
 3*pressure[]*timevec[h$70046]*PD[-h$70046][pertphi[LI[1]]] - 
 density[]*metric\[Delta][p$70047, p$70048]*
  PD[-p$70048][PD[-p$70047][pertB[LI[1]]]] - metric\[Delta][p$70047, p$70048]*
  pressure[]*PD[-p$70048][PD[-p$70047][pertB[LI[1]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$70047, p$70048]*
   PD[-p$70048][PD[-p$70047][pertshear[LI[1]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$70047, p$70048]*pressure[]*
   PD[-p$70048][PD[-p$70047][pertshear[LI[1]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$70047, p$70048]*
   PD[-p$70048][PD[-p$70047][pertvelocity[LI[1]]]])/scale[] + 
 (metric\[Delta][p$70047, p$70048]*pressure[]*
   PD[-p$70048][PD[-p$70047][pertvelocity[LI[1]]]])/scale[] + 
 density[]*metric\[Delta][p$70047, p$70048]*timevec[h$70046]*
  PD[-p$70048][PD[-p$70047][PD[-h$70046][pertE[LI[1]]]]] + 
 metric\[Delta][p$70047, p$70048]*pressure[]*timevec[h$70046]*
  PD[-p$70048][PD[-p$70047][PD[-h$70046][pertE[LI[1]]]]]
