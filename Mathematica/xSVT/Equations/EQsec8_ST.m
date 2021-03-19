(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[2]] - 3*hubbleC[]*pertdensity[LI[2]]*pressure[] - 
 source8[LI[2]] + density[]*timevec[h$424021]*
  PD[-h$424021][pertdensity[LI[2]]] - 3*density[]*timevec[h$424021]*
  PD[-h$424021][pertphi[LI[2]]] - 3*pressure[]*timevec[h$424021]*
  PD[-h$424021][pertphi[LI[2]]] - 
 density[]*metric\[Delta][p$424017, p$424018]*
  PD[-p$424018][PD[-p$424017][pertB[LI[2]]]] - 
 metric\[Delta][p$424017, p$424018]*pressure[]*
  PD[-p$424018][PD[-p$424017][pertB[LI[2]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$424017, p$424018]*
   PD[-p$424018][PD[-p$424017][pertshear[LI[2]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$424017, p$424018]*pressure[]*
   PD[-p$424018][PD[-p$424017][pertshear[LI[2]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$424017, p$424018]*
   PD[-p$424018][PD[-p$424017][pertvelocity[LI[2]]]])/scale[] + 
 (metric\[Delta][p$424017, p$424018]*pressure[]*
   PD[-p$424018][PD[-p$424017][pertvelocity[LI[2]]]])/scale[] + 
 density[]*metric\[Delta][p$424017, p$424018]*timevec[h$424021]*
  PD[-p$424018][PD[-p$424017][PD[-h$424021][pertE[LI[2]]]]] + 
 metric\[Delta][p$424017, p$424018]*pressure[]*timevec[h$424021]*
  PD[-p$424018][PD[-p$424017][PD[-h$424021][pertE[LI[2]]]]]
