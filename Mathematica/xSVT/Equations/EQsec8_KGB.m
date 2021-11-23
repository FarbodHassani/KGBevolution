(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[2]] - 3*hubbleC[]*pertdensity[LI[2]]*pressure[] - 
 source8[LI[2]] + density[]*timevec[h$269063]*
  PD[-h$269063][pertdensity[LI[2]]] - 3*density[]*timevec[h$269063]*
  PD[-h$269063][pertphi[LI[2]]] - 3*pressure[]*timevec[h$269063]*
  PD[-h$269063][pertphi[LI[2]]] - 
 density[]*metric\[Delta][p$269059, p$269060]*
  PD[-p$269060][PD[-p$269059][pertB[LI[2]]]] - 
 metric\[Delta][p$269059, p$269060]*pressure[]*
  PD[-p$269060][PD[-p$269059][pertB[LI[2]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$269059, p$269060]*
   PD[-p$269060][PD[-p$269059][pertshear[LI[2]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$269059, p$269060]*pressure[]*
   PD[-p$269060][PD[-p$269059][pertshear[LI[2]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$269059, p$269060]*
   PD[-p$269060][PD[-p$269059][pertvelocity[LI[2]]]])/scale[] + 
 (metric\[Delta][p$269059, p$269060]*pressure[]*
   PD[-p$269060][PD[-p$269059][pertvelocity[LI[2]]]])/scale[] + 
 density[]*metric\[Delta][p$269059, p$269060]*timevec[h$269063]*
  PD[-p$269060][PD[-p$269059][PD[-h$269063][pertE[LI[2]]]]] + 
 metric\[Delta][p$269059, p$269060]*pressure[]*timevec[h$269063]*
  PD[-p$269060][PD[-p$269059][PD[-h$269063][pertE[LI[2]]]]]
