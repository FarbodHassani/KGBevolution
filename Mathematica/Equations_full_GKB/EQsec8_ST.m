(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[2]] - 3*hubbleC[]*pertdensity[LI[2]]*pressure[] - 
 source8[LI[2]] + density[]*timevec[h$102804]*
  PD[-h$102804][pertdensity[LI[2]]] - 3*density[]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[2]]] - 3*pressure[]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[2]]] - 
 density[]*metric\[Delta][p$102800, p$102801]*
  PD[-p$102801][PD[-p$102800][pertB[LI[2]]]] - 
 metric\[Delta][p$102800, p$102801]*pressure[]*
  PD[-p$102801][PD[-p$102800][pertB[LI[2]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102801][PD[-p$102800][pertshear[LI[2]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$102800, p$102801]*pressure[]*
   PD[-p$102801][PD[-p$102800][pertshear[LI[2]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[2]]]])/scale[] + 
 (metric\[Delta][p$102800, p$102801]*pressure[]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[2]]]])/scale[] + 
 density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
  PD[-p$102801][PD[-p$102800][PD[-h$102804][pertE[LI[2]]]]] + 
 metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
  PD[-p$102801][PD[-p$102800][PD[-h$102804][pertE[LI[2]]]]]
