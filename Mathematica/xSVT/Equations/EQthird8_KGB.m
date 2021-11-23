(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[3]] - 3*hubbleC[]*pertdensity[LI[3]]*pressure[] + 
 (hubbleC[]*pertshearcov[LI[3], p$3144549, -p$3144549])/scale[]^2 - 
 source8[LI[3]] + density[]*timevec[h$3144555]*
  PD[-h$3144555][pertdensity[LI[3]]] - 3*density[]*timevec[h$3144555]*
  PD[-h$3144555][pertphi[LI[3]]] - 3*pressure[]*timevec[h$3144555]*
  PD[-h$3144555][pertphi[LI[3]]] - 
 density[]*metric\[Delta][p$3144549, p$3144550]*
  PD[-p$3144550][PD[-p$3144549][pertB[LI[3]]]] - 
 metric\[Delta][p$3144549, p$3144550]*pressure[]*
  PD[-p$3144550][PD[-p$3144549][pertB[LI[3]]]] + 
 (density[]*metric\[Delta][p$3144549, p$3144550]*
   PD[-p$3144550][PD[-p$3144549][pertvelocity[LI[3]]]])/scale[] + 
 (metric\[Delta][p$3144549, p$3144550]*pressure[]*
   PD[-p$3144550][PD[-p$3144549][pertvelocity[LI[3]]]])/scale[] + 
 density[]*metric\[Delta][p$3144549, p$3144550]*timevec[h$3144555]*
  PD[-p$3144550][PD[-p$3144549][PD[-h$3144555][pertE[LI[3]]]]] + 
 metric\[Delta][p$3144549, p$3144550]*pressure[]*timevec[h$3144555]*
  PD[-p$3144550][PD[-p$3144549][PD[-h$3144555][pertE[LI[3]]]]]
