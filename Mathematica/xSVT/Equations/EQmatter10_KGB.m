(* Created with the Wolfram Language : www.wolfram.com *)
-3*density[]*hubbleC[]*pertdensity[LI[1]] - 3*hubbleC[]*pertpressure[LI[1]] - 
 pertdensity[LI[1]]*primedensity[] + 2*density[]*timevec[h$2780304]*
  PD[-h$2780304][pertpsi[LI[1]]] + 2*pressure[]*timevec[h$2780304]*
  PD[-h$2780304][pertpsi[LI[1]]] - density[]*timevec[h$2780358]*
  PD[-h$2780358][pertdensity[LI[1]]] + 3*density[]*timevec[h$2780852]*
  PD[-h$2780852][pertphi[LI[1]]] + 3*pressure[]*timevec[h$2780852]*
  PD[-h$2780852][pertphi[LI[1]]] - 2*density[]*timevec[h$2780852]*
  PD[-h$2780852][pertpsi[LI[1]]] - 2*pressure[]*timevec[h$2780852]*
  PD[-h$2780852][pertpsi[LI[1]]] + 
 density[]*metric\[Delta][p$2780305, p$2780306]*
  PD[-p$2780306][PD[-p$2780305][pertB[LI[1]]]] + 
 metric\[Delta][p$2780305, p$2780306]*pressure[]*
  PD[-p$2780306][PD[-p$2780305][pertB[LI[1]]]] - 
 (density[]*metric\[Delta][p$2780305, p$2780306]*
   PD[-p$2780306][PD[-p$2780305][pertvelocity[LI[1]]]])/scale[] - 
 (metric\[Delta][p$2780305, p$2780306]*pressure[]*
   PD[-p$2780306][PD[-p$2780305][pertvelocity[LI[1]]]])/scale[] + 
 (density[]*hubbleC[]*metric\[Delta][p$2780498, p$2780499]*
   PD[-p$2780499][PD[-p$2780498][pertshear[LI[1]]]])/scale[]^2 + 
 (hubbleC[]*metric\[Delta][p$2780498, p$2780499]*pressure[]*
   PD[-p$2780499][PD[-p$2780498][pertshear[LI[1]]]])/scale[]^2 - 
 density[]*metric\[Delta][p$2780853, p$2780854]*timevec[h$2780852]*
  PD[-p$2780854][PD[-p$2780853][PD[-h$2780852][pertE[LI[1]]]]] - 
 metric\[Delta][p$2780853, p$2780854]*pressure[]*timevec[h$2780852]*
  PD[-p$2780854][PD[-p$2780853][PD[-h$2780852][pertE[LI[1]]]]]
