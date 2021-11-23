(* Created with the Wolfram Language : www.wolfram.com *)
-4*density[]*hubbleC[]*pertS[LI[1], -i] - 4*hubbleC[]*pertS[LI[1], -i]*
  pressure[] - pertS[LI[1], -i]*primedensity[] - 
 pertS[LI[1], -i]*primepressure[] + 
 (3*density[]*hubbleC[]*metric\[Delta][-i, -p$2786950]*
   pertvelocityvec[LI[1], p$2786950])/scale[] + 
 (3*hubbleC[]*metric\[Delta][-i, -p$2786950]*pertvelocityvec[LI[1], 
    p$2786950]*pressure[])/scale[] + 
 (metric\[Delta][-i, -p$2786950]*pertvelocityvec[LI[1], p$2786950]*
   primedensity[])/scale[] + (metric\[Delta][-i, -p$2786950]*
   pertvelocityvec[LI[1], p$2786950]*primepressure[])/scale[] + 
 (density[]*metric\[Delta][-i, -p$2786950]*timevec[h$2786951]*
   PD[-h$2786951][pertvelocityvec[LI[1], p$2786950]])/scale[] + 
 (metric\[Delta][-i, -p$2786950]*pressure[]*timevec[h$2786951]*
   PD[-h$2786951][pertvelocityvec[LI[1], p$2786950]])/scale[] - 
 density[]*timevec[h$2786951]*PD[-h$2786951][PD[-i][pertB[LI[1]]]] - 
 pressure[]*timevec[h$2786951]*PD[-h$2786951][PD[-i][pertB[LI[1]]]] + 
 (density[]*timevec[h$2786951]*PD[-h$2786951][PD[-i][pertvelocity[LI[1]]]])/
  scale[] + (pressure[]*timevec[h$2786951]*
   PD[-h$2786951][PD[-i][pertvelocity[LI[1]]]])/scale[] - 
 density[]*timevec[h$2787186]*PD[-h$2787186][pertS[LI[1], -i]] - 
 pressure[]*timevec[h$2787186]*PD[-h$2787186][pertS[LI[1], -i]] + 
 density[]*timevec[h$2787186]*PD[-h$2787186][PD[-i][pertB[LI[1]]]] + 
 pressure[]*timevec[h$2787186]*PD[-h$2787186][PD[-i][pertB[LI[1]]]] + 
 PD[-i][pertpressure[LI[1]]] + density[]*PD[-i][pertpsi[LI[1]]] + 
 pressure[]*PD[-i][pertpsi[LI[1]]] + 
 (3*density[]*hubbleC[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (3*hubbleC[]*pressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primedensity[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primepressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (density[]*metric\[Delta][p$2787435, p$2787436]*
   PD[-p$2787436][PD[-p$2787435][pertshearvec[LI[1], -i]]])/scale[]^2 + 
 (metric\[Delta][p$2787435, p$2787436]*pressure[]*
   PD[-p$2787436][PD[-p$2787435][pertshearvec[LI[1], -i]]])/scale[]^2
