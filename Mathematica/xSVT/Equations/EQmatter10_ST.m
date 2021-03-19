(* Created with the Wolfram Language : www.wolfram.com *)
-3*density[]*hubbleC[]*pertdensity[LI[1]] - 3*hubbleC[]*pertpressure[LI[1]] - 
 pertdensity[LI[1]]*primedensity[] + 2*density[]*timevec[h$2916807]*
  PD[-h$2916807][pertpsi[LI[1]]] + 2*pressure[]*timevec[h$2916807]*
  PD[-h$2916807][pertpsi[LI[1]]] - density[]*timevec[h$2916860]*
  PD[-h$2916860][pertdensity[LI[1]]] + 3*density[]*timevec[h$2917351]*
  PD[-h$2917351][pertphi[LI[1]]] + 3*pressure[]*timevec[h$2917351]*
  PD[-h$2917351][pertphi[LI[1]]] - 2*density[]*timevec[h$2917351]*
  PD[-h$2917351][pertpsi[LI[1]]] - 2*pressure[]*timevec[h$2917351]*
  PD[-h$2917351][pertpsi[LI[1]]] + 
 density[]*metric\[Delta][p$2916808, p$2916809]*
  PD[-p$2916809][PD[-p$2916808][pertB[LI[1]]]] + 
 metric\[Delta][p$2916808, p$2916809]*pressure[]*
  PD[-p$2916809][PD[-p$2916808][pertB[LI[1]]]] - 
 (density[]*metric\[Delta][p$2916808, p$2916809]*
   PD[-p$2916809][PD[-p$2916808][pertvelocity[LI[1]]]])/scale[] - 
 (metric\[Delta][p$2916808, p$2916809]*pressure[]*
   PD[-p$2916809][PD[-p$2916808][pertvelocity[LI[1]]]])/scale[] + 
 (density[]*hubbleC[]*metric\[Delta][p$2916998, p$2916999]*
   PD[-p$2916999][PD[-p$2916998][pertshear[LI[1]]]])/scale[]^2 + 
 (hubbleC[]*metric\[Delta][p$2916998, p$2916999]*pressure[]*
   PD[-p$2916999][PD[-p$2916998][pertshear[LI[1]]]])/scale[]^2 - 
 density[]*metric\[Delta][p$2917352, p$2917353]*timevec[h$2917351]*
  PD[-p$2917353][PD[-p$2917352][PD[-h$2917351][pertE[LI[1]]]]] - 
 metric\[Delta][p$2917352, p$2917353]*pressure[]*timevec[h$2917351]*
  PD[-p$2917353][PD[-p$2917352][PD[-h$2917351][pertE[LI[1]]]]]
