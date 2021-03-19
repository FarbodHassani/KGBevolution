(* Created with the Wolfram Language : www.wolfram.com *)
-3*density[]*hubbleC[]*pertdensity[LI[1]] - 3*hubbleC[]*pertpressure[LI[1]] - 
 pertdensity[LI[1]]*primedensity[] - density[]*timevec[h$69685]*
  PD[-h$69685][pertdensity[LI[1]]] + 3*density[]*timevec[h$69685]*
  PD[-h$69685][pertphi[LI[1]]] + 3*pressure[]*timevec[h$69685]*
  PD[-h$69685][pertphi[LI[1]]] + density[]*metric\[Delta][p$69686, p$69687]*
  PD[-p$69687][PD[-p$69686][pertB[LI[1]]]] + metric\[Delta][p$69686, p$69687]*
  pressure[]*PD[-p$69687][PD[-p$69686][pertB[LI[1]]]] + 
 (density[]*hubbleC[]*metric\[Delta][p$69686, p$69687]*
   PD[-p$69687][PD[-p$69686][pertshear[LI[1]]]])/scale[]^2 + 
 (hubbleC[]*metric\[Delta][p$69686, p$69687]*pressure[]*
   PD[-p$69687][PD[-p$69686][pertshear[LI[1]]]])/scale[]^2 - 
 (density[]*metric\[Delta][p$69686, p$69687]*
   PD[-p$69687][PD[-p$69686][pertvelocity[LI[1]]]])/scale[] - 
 (metric\[Delta][p$69686, p$69687]*pressure[]*
   PD[-p$69687][PD[-p$69686][pertvelocity[LI[1]]]])/scale[] - 
 density[]*metric\[Delta][p$69686, p$69687]*timevec[h$69685]*
  PD[-p$69687][PD[-p$69686][PD[-h$69685][pertE[LI[1]]]]] - 
 metric\[Delta][p$69686, p$69687]*pressure[]*timevec[h$69685]*
  PD[-p$69687][PD[-p$69686][PD[-h$69685][pertE[LI[1]]]]]
