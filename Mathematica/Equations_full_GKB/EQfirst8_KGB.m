(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[1]] - 3*hubbleC[]*pertdensity[LI[1]]*pressure[] + 
 density[]*timevec[h$2781740]*PD[-h$2781740][pertdensity[LI[1]]] - 
 3*density[]*timevec[h$2781740]*PD[-h$2781740][pertphi[LI[1]]] - 
 3*pressure[]*timevec[h$2781740]*PD[-h$2781740][pertphi[LI[1]]] - 
 density[]*metric\[Delta][p$2781741, p$2781742]*
  PD[-p$2781742][PD[-p$2781741][pertB[LI[1]]]] - 
 metric\[Delta][p$2781741, p$2781742]*pressure[]*
  PD[-p$2781742][PD[-p$2781741][pertB[LI[1]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$2781741, p$2781742]*
   PD[-p$2781742][PD[-p$2781741][pertshear[LI[1]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$2781741, p$2781742]*pressure[]*
   PD[-p$2781742][PD[-p$2781741][pertshear[LI[1]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$2781741, p$2781742]*
   PD[-p$2781742][PD[-p$2781741][pertvelocity[LI[1]]]])/scale[] + 
 (metric\[Delta][p$2781741, p$2781742]*pressure[]*
   PD[-p$2781742][PD[-p$2781741][pertvelocity[LI[1]]]])/scale[] + 
 density[]*metric\[Delta][p$2781741, p$2781742]*timevec[h$2781740]*
  PD[-p$2781742][PD[-p$2781741][PD[-h$2781740][pertE[LI[1]]]]] + 
 metric\[Delta][p$2781741, p$2781742]*pressure[]*timevec[h$2781740]*
  PD[-p$2781742][PD[-p$2781741][PD[-h$2781740][pertE[LI[1]]]]]
