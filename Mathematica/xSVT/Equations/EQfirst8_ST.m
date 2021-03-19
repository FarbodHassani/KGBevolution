(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]*pertpressure[LI[1]] - 3*hubbleC[]*pertdensity[LI[1]]*pressure[] + 
 density[]*timevec[h$2917844]*PD[-h$2917844][pertdensity[LI[1]]] - 
 3*density[]*timevec[h$2917844]*PD[-h$2917844][pertphi[LI[1]]] - 
 3*pressure[]*timevec[h$2917844]*PD[-h$2917844][pertphi[LI[1]]] - 
 density[]*metric\[Delta][p$2917845, p$2917846]*
  PD[-p$2917846][PD[-p$2917845][pertB[LI[1]]]] - 
 metric\[Delta][p$2917845, p$2917846]*pressure[]*
  PD[-p$2917846][PD[-p$2917845][pertB[LI[1]]]] - 
 (density[]*hubbleC[]*metric\[Delta][p$2917845, p$2917846]*
   PD[-p$2917846][PD[-p$2917845][pertshear[LI[1]]]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][p$2917845, p$2917846]*pressure[]*
   PD[-p$2917846][PD[-p$2917845][pertshear[LI[1]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$2917845, p$2917846]*
   PD[-p$2917846][PD[-p$2917845][pertvelocity[LI[1]]]])/scale[] + 
 (metric\[Delta][p$2917845, p$2917846]*pressure[]*
   PD[-p$2917846][PD[-p$2917845][pertvelocity[LI[1]]]])/scale[] + 
 density[]*metric\[Delta][p$2917845, p$2917846]*timevec[h$2917844]*
  PD[-p$2917846][PD[-p$2917845][PD[-h$2917844][pertE[LI[1]]]]] + 
 metric\[Delta][p$2917845, p$2917846]*pressure[]*timevec[h$2917844]*
  PD[-p$2917846][PD[-p$2917845][PD[-h$2917844][pertE[LI[1]]]]]
