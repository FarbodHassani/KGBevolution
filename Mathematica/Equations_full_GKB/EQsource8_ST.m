(* Created with the Wolfram Language : www.wolfram.com *)
(2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
   pertshearten[LI[1], -p$102801, -p$102803])/scale[]^2 + 
 (2*density[]*hubbleC[]*metric\[Delta][-p$102800, -p$102801]*
   pertvelocityvec[LI[1], p$102800]*pertvelocityvec[LI[1], p$102801])/
  scale[]^2 + (2*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
   pertshearten[LI[1], -p$102801, -p$102803]*pressure[])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][-p$102800, -p$102801]*
   pertvelocityvec[LI[1], p$102800]*pertvelocityvec[LI[1], p$102801]*
   pressure[])/scale[]^2 - (2*metric\[Delta][-p$102800, -p$102801]*
   pertvelocityvec[LI[1], p$102800]*pertvelocityvec[LI[1], p$102801]*
   primepressure[])/scale[]^2 + 
 (2*pertS[LI[1], -p$102800]*pertvelocityvec[LI[1], p$102800]*primepressure[])/
  scale[] - source8[LI[2]] + density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
  timevec[h$102804]*PD[-h$102804][perth[LI[1], -p$102801, -p$102803]] + 
 metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  perth[LI[1], -p$102800, -p$102802]*pressure[]*timevec[h$102804]*
  PD[-h$102804][perth[LI[1], -p$102801, -p$102803]] - 
 (density[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*pertshearten[LI[1], -p$102800, -p$102802]*timevec[h$102804]*
   PD[-h$102804][perth[LI[1], -p$102801, -p$102803]])/scale[]^2 - 
 (metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pertshearten[LI[1], -p$102800, -p$102802]*pressure[]*timevec[h$102804]*
   PD[-h$102804][perth[LI[1], -p$102801, -p$102803]])/scale[]^2 + 
 6*density[]*pertdensity[LI[1]]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[1]]] + 12*density[]*pertphi[LI[1]]*
  timevec[h$102804]*PD[-h$102804][pertphi[LI[1]]] + 
 6*pertpressure[LI[1]]*timevec[h$102804]*PD[-h$102804][pertphi[LI[1]]] + 
 12*pertphi[LI[1]]*pressure[]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[1]]] + 
 (4*density[]*pertvelocityvec[LI[1], p$102800]*timevec[h$102804]*
   PD[-h$102804][pertS[LI[1], -p$102800]])/scale[] + 
 (4*pertvelocityvec[LI[1], p$102800]*pressure[]*timevec[h$102804]*
   PD[-h$102804][pertS[LI[1], -p$102800]])/scale[] - 
 2*density[]*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*
  timevec[h$102804]*PD[-h$102804][pertS[LI[1], -p$102801]] - 
 2*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*pressure[]*
  timevec[h$102804]*PD[-h$102804][pertS[LI[1], -p$102801]] + 
 (2*density[]*pertS[LI[1], -p$102800]*timevec[h$102804]*
   PD[-h$102804][pertvelocityvec[LI[1], p$102800]])/scale[] + 
 (2*pertS[LI[1], -p$102800]*pressure[]*timevec[h$102804]*
   PD[-h$102804][pertvelocityvec[LI[1], p$102800]])/scale[] - 
 (4*density[]*metric\[Delta][-p$102800, -p$102801]*
   pertvelocityvec[LI[1], p$102800]*timevec[h$102804]*
   PD[-h$102804][pertvelocityvec[LI[1], p$102801]])/scale[]^2 - 
 (4*metric\[Delta][-p$102800, -p$102801]*pertvelocityvec[LI[1], p$102800]*
   pressure[]*timevec[h$102804]*PD[-h$102804][pertvelocityvec[LI[1], 
     p$102801]])/scale[]^2 + (2*pertvelocityvec[LI[1], p$102800]*
   primepressure[]*PD[-p$102800][pertB[LI[1]]])/scale[] - 
 2*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
  PD[-h$102804][pertS[LI[1], -p$102801]]*PD[-p$102800][pertB[LI[1]]] - 
 2*metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
  PD[-h$102804][pertS[LI[1], -p$102801]]*PD[-p$102800][pertB[LI[1]]] + 
 (2*density[]*timevec[h$102804]*PD[-h$102804][pertvelocityvec[LI[1], 
     p$102800]]*PD[-p$102800][pertB[LI[1]]])/scale[] + 
 (2*pressure[]*timevec[h$102804]*PD[-h$102804][pertvelocityvec[LI[1], 
     p$102800]]*PD[-p$102800][pertB[LI[1]]])/scale[] - 
 (2*density[]*pertvelocityvec[LI[1], p$102800]*
   PD[-p$102800][pertdensity[LI[1]]])/scale[] + 
 (6*density[]*pertvelocityvec[LI[1], p$102800]*PD[-p$102800][pertphi[LI[1]]])/
  scale[] + (6*pertvelocityvec[LI[1], p$102800]*pressure[]*
   PD[-p$102800][pertphi[LI[1]]])/scale[] - 
 (2*pertvelocityvec[LI[1], p$102800]*PD[-p$102800][pertpressure[LI[1]]])/
  scale[] - (4*density[]*pertvelocityvec[LI[1], p$102800]*
   PD[-p$102800][pertpsi[LI[1]]])/scale[] - 
 (4*pertvelocityvec[LI[1], p$102800]*pressure[]*
   PD[-p$102800][pertpsi[LI[1]]])/scale[] + 
 (4*density[]*hubbleC[]*pertvelocityvec[LI[1], p$102800]*
   PD[-p$102800][pertvelocity[LI[1]]])/scale[]^2 + 
 (4*hubbleC[]*pertvelocityvec[LI[1], p$102800]*pressure[]*
   PD[-p$102800][pertvelocity[LI[1]]])/scale[]^2 - 
 (4*pertvelocityvec[LI[1], p$102800]*primepressure[]*
   PD[-p$102800][pertvelocity[LI[1]]])/scale[]^2 - 
 (4*density[]*timevec[h$102804]*PD[-h$102804][pertvelocityvec[LI[1], 
     p$102800]]*PD[-p$102800][pertvelocity[LI[1]]])/scale[]^2 - 
 (4*pressure[]*timevec[h$102804]*PD[-h$102804][pertvelocityvec[LI[1], 
     p$102800]]*PD[-p$102800][pertvelocity[LI[1]]])/scale[]^2 - 
 (4*density[]*pertvelocityvec[LI[1], p$102800]*timevec[h$102804]*
   PD[-p$102800][PD[-h$102804][pertvelocity[LI[1]]]])/scale[]^2 - 
 (4*pertvelocityvec[LI[1], p$102800]*pressure[]*timevec[h$102804]*
   PD[-p$102800][PD[-h$102804][pertvelocity[LI[1]]]])/scale[]^2 - 
 2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
  pertS[LI[1], -p$102800]*PD[-p$102801][pertB[LI[1]]] - 
 2*hubbleC[]*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*
  pressure[]*PD[-p$102801][pertB[LI[1]]] - 
 2*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*primepressure[]*
  PD[-p$102801][pertB[LI[1]]] + 2*density[]*metric\[Delta][p$102800, 
   p$102801]*PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][pertdensity[LI[1]]] + 
 (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*pertshearten[LI[1], -p$102802, 
    -p$102803]*PD[-p$102801][pertF[LI[1], -p$102800]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*pertshearten[LI[1], -p$102802, 
    -p$102803]*pressure[]*PD[-p$102801][pertF[LI[1], -p$102800]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102802]*
  metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
  PD[-h$102804][perth[LI[1], -p$102802, -p$102803]]*
  PD[-p$102801][pertF[LI[1], -p$102800]] + 
 2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
  pressure[]*timevec[h$102804]*PD[-h$102804][perth[LI[1], -p$102802, 
    -p$102803]]*PD[-p$102801][pertF[LI[1], -p$102800]] - 
 6*density[]*metric\[Delta][p$102800, p$102801]*PD[-p$102800][pertB[LI[1]]]*
  PD[-p$102801][pertphi[LI[1]]] - 6*metric\[Delta][p$102800, p$102801]*
  pressure[]*PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][pertphi[LI[1]]] + 
 2*metric\[Delta][p$102800, p$102801]*PD[-p$102800][pertB[LI[1]]]*
  PD[-p$102801][pertpressure[LI[1]]] + 
 4*density[]*metric\[Delta][p$102800, p$102801]*PD[-p$102800][pertB[LI[1]]]*
  PD[-p$102801][pertpsi[LI[1]]] + 4*metric\[Delta][p$102800, p$102801]*
  pressure[]*PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][pertpsi[LI[1]]] + 
 (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*perth[LI[1], -p$102802, -p$102803]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*perth[LI[1], -p$102802, -p$102803]*
   pressure[]*PD[-p$102801][pertshearvec[LI[1], -p$102800]])/scale[]^2 - 
 (2*density[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
   PD[-h$102804][perth[LI[1], -p$102802, -p$102803]]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
   pressure[]*timevec[h$102804]*PD[-h$102804][perth[LI[1], -p$102802, 
     -p$102803]]*PD[-p$102801][pertshearvec[LI[1], -p$102800]])/scale[]^2 + 
 (2*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*
   primepressure[]*PD[-p$102801][pertvelocity[LI[1]]])/scale[] + 
 (4*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
   PD[-h$102804][pertS[LI[1], -p$102800]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] + (4*metric\[Delta][p$102800, p$102801]*pressure[]*
   timevec[h$102804]*PD[-h$102804][pertS[LI[1], -p$102800]]*
   PD[-p$102801][pertvelocity[LI[1]]])/scale[] + 
 (2*metric\[Delta][p$102800, p$102801]*primepressure[]*
   PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/scale[] - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102800][pertdensity[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] + (6*density[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102800][pertphi[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] + (6*metric\[Delta][p$102800, p$102801]*pressure[]*
   PD[-p$102800][pertphi[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] - (2*metric\[Delta][p$102800, p$102801]*
   PD[-p$102800][pertpressure[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] - (4*density[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102800][pertpsi[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] - (4*metric\[Delta][p$102800, p$102801]*pressure[]*
   PD[-p$102800][pertpsi[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[] + (2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   PD[-p$102800][pertvelocity[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[]^2 + (2*hubbleC[]*metric\[Delta][p$102800, p$102801]*pressure[]*
   PD[-p$102800][pertvelocity[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[]^2 - (2*metric\[Delta][p$102800, p$102801]*primepressure[]*
   PD[-p$102800][pertvelocity[LI[1]]]*PD[-p$102801][pertvelocity[LI[1]]])/
  scale[]^2 + (2*density[]*metric\[Delta][p$102800, p$102801]*
   pertS[LI[1], -p$102800]*timevec[h$102804]*
   PD[-p$102801][PD[-h$102804][pertvelocity[LI[1]]]])/scale[] + 
 (2*metric\[Delta][p$102800, p$102801]*pertS[LI[1], -p$102800]*pressure[]*
   timevec[h$102804]*PD[-p$102801][PD[-h$102804][pertvelocity[LI[1]]]])/
  scale[] + (2*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
   PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][
    PD[-h$102804][pertvelocity[LI[1]]]])/scale[] + 
 (2*metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
   PD[-p$102800][pertB[LI[1]]]*PD[-p$102801][
    PD[-h$102804][pertvelocity[LI[1]]]])/scale[] - 
 (4*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
   PD[-p$102800][pertvelocity[LI[1]]]*PD[-p$102801][
    PD[-h$102804][pertvelocity[LI[1]]]])/scale[]^2 - 
 (4*metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
   PD[-p$102800][pertvelocity[LI[1]]]*PD[-p$102801][
    PD[-h$102804][pertvelocity[LI[1]]]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102801]*pertdensity[LI[1]]*
  PD[-p$102801][PD[-p$102800][pertB[LI[1]]]] + 
 2*metric\[Delta][p$102800, p$102801]*pertpressure[LI[1]]*
  PD[-p$102801][PD[-p$102800][pertB[LI[1]]]] + 
 2*density[]*metric\[Delta][p$102800, p$102801]*pertpsi[LI[1]]*
  PD[-p$102801][PD[-p$102800][pertB[LI[1]]]] + 
 2*metric\[Delta][p$102800, p$102801]*pertpsi[LI[1]]*pressure[]*
  PD[-p$102801][PD[-p$102800][pertB[LI[1]]]] - 
 4*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[1]]]*PD[-p$102801][PD[-p$102800][pertE[LI[1]]]] - 
 4*metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
  PD[-h$102804][pertphi[LI[1]]]*PD[-p$102801][PD[-p$102800][pertE[LI[1]]]] + 
 (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*pertphi[LI[1]]*
   PD[-p$102801][PD[-p$102800][pertshear[LI[1]]]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102801]*pertphi[LI[1]]*pressure[]*
   PD[-p$102801][PD[-p$102800][pertshear[LI[1]]]])/scale[]^2 - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*timevec[h$102804]*
   PD[-h$102804][pertphi[LI[1]]]*PD[-p$102801][
    PD[-p$102800][pertshear[LI[1]]]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102801]*pressure[]*timevec[h$102804]*
   PD[-h$102804][pertphi[LI[1]]]*PD[-p$102801][
    PD[-p$102800][pertshear[LI[1]]]])/scale[]^2 - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*pertdensity[LI[1]]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[1]]]])/scale[] - 
 (2*metric\[Delta][p$102800, p$102801]*pertpressure[LI[1]]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[1]]]])/scale[] - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*pertpsi[LI[1]]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[1]]]])/scale[] - 
 (2*metric\[Delta][p$102800, p$102801]*pertpsi[LI[1]]*pressure[]*
   PD[-p$102801][PD[-p$102800][pertvelocity[LI[1]]]])/scale[] - 
 2*density[]*metric\[Delta][p$102800, p$102801]*pertdensity[LI[1]]*
  timevec[h$102804]*PD[-p$102801][PD[-p$102800][
    PD[-h$102804][pertE[LI[1]]]]] - 
 4*density[]*metric\[Delta][p$102800, p$102801]*pertphi[LI[1]]*
  timevec[h$102804]*PD[-p$102801][PD[-p$102800][
    PD[-h$102804][pertE[LI[1]]]]] - 2*metric\[Delta][p$102800, p$102801]*
  pertpressure[LI[1]]*timevec[h$102804]*
  PD[-p$102801][PD[-p$102800][PD[-h$102804][pertE[LI[1]]]]] - 
 4*metric\[Delta][p$102800, p$102801]*pertphi[LI[1]]*pressure[]*
  timevec[h$102804]*PD[-p$102801][PD[-p$102800][
    PD[-h$102804][pertE[LI[1]]]]] + 
 2*density[]*metric\[Delta][p$102800, p$102802]*
  metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
  PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102802][PD[-h$102804][pertF[LI[1], -p$102803]]] + 
 2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
  pressure[]*timevec[h$102804]*PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102802][PD[-h$102804][pertF[LI[1], -p$102803]]] - 
 (2*density[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102802][PD[-h$102804][pertF[LI[1], -p$102803]]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102802][PD[-h$102804][pertF[LI[1], -p$102803]]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*timevec[h$102804]*
  PD[-h$102804][perth[LI[1], -p$102801, -p$102803]]*
  PD[-p$102802][PD[-p$102800][pertE[LI[1]]]] + 
 2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  pressure[]*timevec[h$102804]*PD[-h$102804][perth[LI[1], -p$102801, 
    -p$102803]]*PD[-p$102802][PD[-p$102800][pertE[LI[1]]]] - 
 (density[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*timevec[h$102804]*PD[-h$102804][perth[LI[1], -p$102801, 
     -p$102803]]*PD[-p$102802][PD[-p$102800][pertshear[LI[1]]]])/
  (2*scale[]^2) - (metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pressure[]*timevec[h$102804]*
   PD[-h$102804][perth[LI[1], -p$102801, -p$102803]]*
   PD[-p$102802][PD[-p$102800][pertshear[LI[1]]]])/(2*scale[]^2) - 
 (2*density[]*metric\[Delta][p$102801, p$102802]*
   pertvelocityvec[LI[1], p$102800]*PD[-p$102802][
    PD[-p$102801][PD[-p$102800][pertE[LI[1]]]]])/scale[] - 
 (2*metric\[Delta][p$102801, p$102802]*pertvelocityvec[LI[1], p$102800]*
   pressure[]*PD[-p$102802][PD[-p$102801][PD[-p$102800][pertE[LI[1]]]]])/
  scale[] + (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102803]*
   metric\[Delta][p$102801, p$102802]*PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][pertshearvec[LI[1], -p$102802]])/scale[]^2 + 
 (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][pertshearvec[LI[1], -p$102802]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102803]*
   metric\[Delta][p$102801, p$102802]*pressure[]*
   PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][pertshearvec[LI[1], -p$102802]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*pressure[]*
   PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][pertshearvec[LI[1], -p$102802]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
  timevec[h$102804]*PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102801]]] + 
 2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  perth[LI[1], -p$102800, -p$102802]*pressure[]*timevec[h$102804]*
  PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102801]]] - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pertshearten[LI[1], -p$102800, 
    -p$102802]*timevec[h$102804]*PD[-p$102803][
    PD[-h$102804][pertF[LI[1], -p$102801]]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pertshearten[LI[1], -p$102800, -p$102802]*pressure[]*timevec[h$102804]*
   PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102801]]])/scale[]^2 + 
 4*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*timevec[h$102804]*
  PD[-p$102802][PD[-p$102800][pertE[LI[1]]]]*
  PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102801]]] + 
 4*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  pressure[]*timevec[h$102804]*PD[-p$102802][PD[-p$102800][pertE[LI[1]]]]*
  PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102801]]] + 
 2*density[]*metric\[Delta][p$102800, p$102802]*
  metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
  PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102802]]] + 
 2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
  pressure[]*timevec[h$102804]*PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102802]]] - 
 (2*density[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102802]]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102803][PD[-h$102804][pertF[LI[1], -p$102802]]])/scale[]^2 + 
 (4*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pertshearten[LI[1], -p$102800, 
    -p$102802]*PD[-p$102803][PD[-p$102801][pertE[LI[1]]]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pertshearten[LI[1], -p$102800, 
    -p$102802]*pressure[]*PD[-p$102803][PD[-p$102801][pertE[LI[1]]]])/
  scale[]^2 + (density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
   PD[-p$102803][PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 + 
 (hubbleC[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*perth[LI[1], -p$102800, -p$102802]*pressure[]*
   PD[-p$102803][PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 - 
 (density[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*timevec[h$102804]*PD[-p$102802][
    PD[-h$102804][pertF[LI[1], -p$102800]]]*
   PD[-p$102803][PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 - 
 (metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102802][
    PD[-h$102804][pertF[LI[1], -p$102800]]]*
   PD[-p$102803][PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 + 
 (2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*PD[-p$102802][
    PD[-p$102800][pertE[LI[1]]]]*PD[-p$102803][
    PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pressure[]*
   PD[-p$102802][PD[-p$102800][pertE[LI[1]]]]*
   PD[-p$102803][PD[-p$102801][pertshear[LI[1]]]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
  timevec[h$102804]*PD[-p$102803][PD[-p$102801][
    PD[-h$102804][pertE[LI[1]]]]] + 2*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*perth[LI[1], -p$102800, -p$102802]*
  pressure[]*timevec[h$102804]*PD[-p$102803][
   PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]] - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pertshearten[LI[1], -p$102800, 
    -p$102802]*timevec[h$102804]*PD[-p$102803][
    PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 - 
 (2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pertshearten[LI[1], -p$102800, -p$102802]*pressure[]*timevec[h$102804]*
   PD[-p$102803][PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 + 
 4*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*timevec[h$102804]*
  PD[-p$102802][PD[-p$102800][pertE[LI[1]]]]*
  PD[-p$102803][PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]] + 
 4*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  pressure[]*timevec[h$102804]*PD[-p$102802][PD[-p$102800][pertE[LI[1]]]]*
  PD[-p$102803][PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]] - 
 (density[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*timevec[h$102804]*PD[-p$102802][
    PD[-p$102800][pertshear[LI[1]]]]*PD[-p$102803][
    PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 - 
 (metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102802][
    PD[-p$102800][pertshear[LI[1]]]]*PD[-p$102803][
    PD[-p$102801][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 + 
 (8*density[]*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*PD[-p$102801][
    pertshearvec[LI[1], -p$102800]]*PD[-p$102803][
    PD[-p$102802][pertE[LI[1]]]])/scale[]^2 + 
 (8*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*pressure[]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102803][PD[-p$102802][pertE[LI[1]]]])/scale[]^2 + 
 (2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][PD[-p$102802][pertshear[LI[1]]]])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*pressure[]*
   PD[-p$102801][pertF[LI[1], -p$102800]]*
   PD[-p$102803][PD[-p$102802][pertshear[LI[1]]]])/scale[]^2 - 
 (2*density[]*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*PD[-p$102801][
    PD[-p$102800][pertE[LI[1]]]]*PD[-p$102803][
    PD[-p$102802][pertshear[LI[1]]]])/scale[]^2 - 
 (2*hubbleC[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*pressure[]*
   PD[-p$102801][PD[-p$102800][pertE[LI[1]]]]*
   PD[-p$102803][PD[-p$102802][pertshear[LI[1]]]])/scale[]^2 + 
 4*density[]*metric\[Delta][p$102800, p$102802]*
  metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
  PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102803][PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]] + 
 4*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
  pressure[]*timevec[h$102804]*PD[-p$102801][pertF[LI[1], -p$102800]]*
  PD[-p$102803][PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]] - 
 (4*density[]*metric\[Delta][p$102800, p$102802]*
   metric\[Delta][p$102801, p$102803]*timevec[h$102804]*
   PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102803][PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 - 
 (4*metric\[Delta][p$102800, p$102802]*metric\[Delta][p$102801, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102801][pertshearvec[LI[1], -p$102800]]*
   PD[-p$102803][PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 + 
 (density[]*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, 
    p$102803]*timevec[h$102804]*PD[-p$102801][
    PD[-p$102800][pertshear[LI[1]]]]*PD[-p$102803][
    PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 + 
 (metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pressure[]*timevec[h$102804]*PD[-p$102801][
    PD[-p$102800][pertshear[LI[1]]]]*PD[-p$102803][
    PD[-p$102802][PD[-h$102804][pertE[LI[1]]]]])/scale[]^2 + 
 2*density[]*metric\[Delta][p$102800, p$102801]*
  metric\[Delta][p$102802, p$102803]*PD[-p$102800][pertB[LI[1]]]*
  PD[-p$102803][PD[-p$102802][PD[-p$102801][pertE[LI[1]]]]] + 
 2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
  pressure[]*PD[-p$102800][pertB[LI[1]]]*
  PD[-p$102803][PD[-p$102802][PD[-p$102801][pertE[LI[1]]]]] - 
 (2*density[]*metric\[Delta][p$102800, p$102801]*
   metric\[Delta][p$102802, p$102803]*PD[-p$102800][pertvelocity[LI[1]]]*
   PD[-p$102803][PD[-p$102802][PD[-p$102801][pertE[LI[1]]]]])/scale[] - 
 (2*metric\[Delta][p$102800, p$102801]*metric\[Delta][p$102802, p$102803]*
   pressure[]*PD[-p$102800][pertvelocity[LI[1]]]*
   PD[-p$102803][PD[-p$102802][PD[-p$102801][pertE[LI[1]]]]])/scale[]
