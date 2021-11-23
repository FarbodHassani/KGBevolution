(* Created with the Wolfram Language : www.wolfram.com *)
-3*density[]*hubbleC[]*pertdensity[LI[1]] - 3*hubbleC[]*pertpressure[LI[1]] - 
 pertdensity[LI[1]]*primedensity[] + 2*density[]*timevec[h$2305276]*
  PD[-h$2305276][pertpsi[LI[1]]] + 2*pressure[]*timevec[h$2305276]*
  PD[-h$2305276][pertpsi[LI[1]]] - density[]*timevec[h$2305330]*
  PD[-h$2305330][pertdensity[LI[1]]] + 3*density[]*timevec[h$2305824]*
  PD[-h$2305824][pertphi[LI[1]]] + 3*pressure[]*timevec[h$2305824]*
  PD[-h$2305824][pertphi[LI[1]]] - 2*density[]*timevec[h$2305824]*
  PD[-h$2305824][pertpsi[LI[1]]] - 2*pressure[]*timevec[h$2305824]*
  PD[-h$2305824][pertpsi[LI[1]]] + 
 density[]*metric\[Delta][p$2305277, p$2305278]*
  PD[-p$2305278][PD[-p$2305277][pertB[LI[1]]]] + 
 metric\[Delta][p$2305277, p$2305278]*pressure[]*
  PD[-p$2305278][PD[-p$2305277][pertB[LI[1]]]] - 
 (density[]*metric\[Delta][p$2305277, p$2305278]*
   PD[-p$2305278][PD[-p$2305277][pertvelocity[LI[1]]]])/scale[] - 
 (metric\[Delta][p$2305277, p$2305278]*pressure[]*
   PD[-p$2305278][PD[-p$2305277][pertvelocity[LI[1]]]])/scale[] + 
 (density[]*hubbleC[]*metric\[Delta][p$2305470, p$2305471]*
   PD[-p$2305471][PD[-p$2305470][pertshear[LI[1]]]])/scale[]^2 + 
 (hubbleC[]*metric\[Delta][p$2305470, p$2305471]*pressure[]*
   PD[-p$2305471][PD[-p$2305470][pertshear[LI[1]]]])/scale[]^2 - 
 density[]*metric\[Delta][p$2305825, p$2305826]*timevec[h$2305824]*
  PD[-p$2305826][PD[-p$2305825][PD[-h$2305824][pertE[LI[1]]]]] - 
 metric\[Delta][p$2305825, p$2305826]*pressure[]*timevec[h$2305824]*
  PD[-p$2305826][PD[-p$2305825][PD[-h$2305824][pertE[LI[1]]]]]
