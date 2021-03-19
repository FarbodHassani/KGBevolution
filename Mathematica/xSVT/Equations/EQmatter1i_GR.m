(* Created with the Wolfram Language : www.wolfram.com *)
-4*density[]*hubbleC[]*pertS[LI[1], -i] - 4*hubbleC[]*pertS[LI[1], -i]*
  pressure[] - pertS[LI[1], -i]*primedensity[] - 
 pertS[LI[1], -i]*primepressure[] + 
 (3*density[]*hubbleC[]*metric\[Delta][-i, -p$73537]*
   pertvelocityvec[LI[1], p$73537])/scale[] + 
 (3*hubbleC[]*metric\[Delta][-i, -p$73537]*pertvelocityvec[LI[1], p$73537]*
   pressure[])/scale[] + (metric\[Delta][-i, -p$73537]*
   pertvelocityvec[LI[1], p$73537]*primedensity[])/scale[] + 
 (metric\[Delta][-i, -p$73537]*pertvelocityvec[LI[1], p$73537]*
   primepressure[])/scale[] - density[]*timevec[h$73538]*
  PD[-h$73538][pertS[LI[1], -i]] - pressure[]*timevec[h$73538]*
  PD[-h$73538][pertS[LI[1], -i]] + 
 (density[]*metric\[Delta][-i, -p$73537]*timevec[h$73538]*
   PD[-h$73538][pertvelocityvec[LI[1], p$73537]])/scale[] + 
 (metric\[Delta][-i, -p$73537]*pressure[]*timevec[h$73538]*
   PD[-h$73538][pertvelocityvec[LI[1], p$73537]])/scale[] + 
 (density[]*timevec[h$73538]*PD[-h$73538][PD[-i][pertvelocity[LI[1]]]])/
  scale[] + (pressure[]*timevec[h$73538]*
   PD[-h$73538][PD[-i][pertvelocity[LI[1]]]])/scale[] + 
 PD[-i][pertpressure[LI[1]]] + density[]*PD[-i][pertpsi[LI[1]]] + 
 pressure[]*PD[-i][pertpsi[LI[1]]] + 
 (3*density[]*hubbleC[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (3*hubbleC[]*pressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primedensity[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primepressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (density[]*metric\[Delta][p$73537, p$73539]*
   PD[-p$73539][PD[-p$73537][pertshearvec[LI[1], -i]]])/scale[]^2 + 
 (metric\[Delta][p$73537, p$73539]*pressure[]*
   PD[-p$73539][PD[-p$73537][pertshearvec[LI[1], -i]]])/scale[]^2
