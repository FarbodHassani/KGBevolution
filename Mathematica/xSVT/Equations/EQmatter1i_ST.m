(* Created with the Wolfram Language : www.wolfram.com *)
-4*density[]*hubbleC[]*pertS[LI[1], -i] - 4*hubbleC[]*pertS[LI[1], -i]*
  pressure[] - pertS[LI[1], -i]*primedensity[] - 
 pertS[LI[1], -i]*primepressure[] + 
 (3*density[]*hubbleC[]*metric\[Delta][-i, -p$2310953]*
   pertvelocityvec[LI[1], p$2310953])/scale[] + 
 (3*hubbleC[]*metric\[Delta][-i, -p$2310953]*pertvelocityvec[LI[1], 
    p$2310953]*pressure[])/scale[] + 
 (metric\[Delta][-i, -p$2310953]*pertvelocityvec[LI[1], p$2310953]*
   primedensity[])/scale[] + (metric\[Delta][-i, -p$2310953]*
   pertvelocityvec[LI[1], p$2310953]*primepressure[])/scale[] + 
 (density[]*metric\[Delta][-i, -p$2310953]*timevec[h$2310954]*
   PD[-h$2310954][pertvelocityvec[LI[1], p$2310953]])/scale[] + 
 (metric\[Delta][-i, -p$2310953]*pressure[]*timevec[h$2310954]*
   PD[-h$2310954][pertvelocityvec[LI[1], p$2310953]])/scale[] - 
 density[]*timevec[h$2310954]*PD[-h$2310954][PD[-i][pertB[LI[1]]]] - 
 pressure[]*timevec[h$2310954]*PD[-h$2310954][PD[-i][pertB[LI[1]]]] + 
 (density[]*timevec[h$2310954]*PD[-h$2310954][PD[-i][pertvelocity[LI[1]]]])/
  scale[] + (pressure[]*timevec[h$2310954]*
   PD[-h$2310954][PD[-i][pertvelocity[LI[1]]]])/scale[] - 
 density[]*timevec[h$2311189]*PD[-h$2311189][pertS[LI[1], -i]] - 
 pressure[]*timevec[h$2311189]*PD[-h$2311189][pertS[LI[1], -i]] + 
 density[]*timevec[h$2311189]*PD[-h$2311189][PD[-i][pertB[LI[1]]]] + 
 pressure[]*timevec[h$2311189]*PD[-h$2311189][PD[-i][pertB[LI[1]]]] + 
 PD[-i][pertpressure[LI[1]]] + density[]*PD[-i][pertpsi[LI[1]]] + 
 pressure[]*PD[-i][pertpsi[LI[1]]] + 
 (3*density[]*hubbleC[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (3*hubbleC[]*pressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primedensity[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (primepressure[]*PD[-i][pertvelocity[LI[1]]])/scale[] + 
 (density[]*metric\[Delta][p$2311438, p$2311439]*
   PD[-p$2311439][PD[-p$2311438][pertshearvec[LI[1], -i]]])/scale[]^2 + 
 (metric\[Delta][p$2311438, p$2311439]*pressure[]*
   PD[-p$2311439][PD[-p$2311438][pertshearvec[LI[1], -i]]])/scale[]^2
