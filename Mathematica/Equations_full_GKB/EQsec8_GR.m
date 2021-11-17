(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[2]] - 3*hubbleC[]*pertdensity[LI[2]]*pressure[] - 
 source8[LI[2]] + density[]*timevec[h$48880]*
  PD[-h$48880][pertdensity[LI[2]]] - 3*density[]*timevec[h$48880]*
  PD[-h$48880][pertphi[LI[2]]] - 3*pressure[]*timevec[h$48880]*
  PD[-h$48880][pertphi[LI[2]]] - density[]*metric\[Delta][p$48876, p$48877]*
  PD[-p$48877][PD[-p$48876][pertB[LI[2]]]] - metric\[Delta][p$48876, p$48877]*
  pressure[]*PD[-p$48877][PD[-p$48876][pertB[LI[2]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$48876, p$48877]*
   PD[-p$48877][PD[-p$48876][pertshear[LI[2]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$48876, p$48877]*pressure[]*
   PD[-p$48877][PD[-p$48876][pertshear[LI[2]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$48876, p$48877]*
   PD[-p$48877][PD[-p$48876][pertvelocity[LI[2]]]])/scale[] + 
 (metric\[Delta][p$48876, p$48877]*pressure[]*
   PD[-p$48877][PD[-p$48876][pertvelocity[LI[2]]]])/scale[] + 
 density[]*metric\[Delta][p$48876, p$48877]*timevec[h$48880]*
  PD[-p$48877][PD[-p$48876][PD[-h$48880][pertE[LI[2]]]]] + 
 metric\[Delta][p$48876, p$48877]*pressure[]*timevec[h$48880]*
  PD[-p$48877][PD[-p$48876][PD[-h$48880][pertE[LI[2]]]]]
