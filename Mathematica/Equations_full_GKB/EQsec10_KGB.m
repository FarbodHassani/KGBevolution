(* Created with the Wolfram Language : www.wolfram.com *)
-source10[LI[2], -i] + density[]*hubbleC[]*metric\[Delta][p$545469, p$545470]*
  PD[-p$545470][PD[-p$545469][pertS[LI[2], -i]]] + 
 hubbleC[]*metric\[Delta][p$545469, p$545470]*pressure[]*
  PD[-p$545470][PD[-p$545469][pertS[LI[2], -i]]] + 
 metric\[Delta][p$545469, p$545470]*primepressure[]*
  PD[-p$545470][PD[-p$545469][pertS[LI[2], -i]]] + 
 density[]*metric\[Delta][p$545469, p$545470]*timevec[h$545471]*
  PD[-p$545470][PD[-p$545469][PD[-h$545471][pertS[LI[2], -i]]]] + 
 metric\[Delta][p$545469, p$545470]*pressure[]*timevec[h$545471]*
  PD[-p$545470][PD[-p$545469][PD[-h$545471][pertS[LI[2], -i]]]] - 
 (metric\[Delta][-i, -p$545469]*metric\[Delta][p$545470, p$545472]*
   primepressure[]*PD[-p$545472][PD[-p$545470][pertvelocityvec[LI[2], 
      p$545469]]])/scale[] - (density[]*metric\[Delta][-i, -p$545469]*
   metric\[Delta][p$545470, p$545472]*timevec[h$545471]*
   PD[-p$545472][PD[-p$545470][PD[-h$545471][pertvelocityvec[LI[2], 
       p$545469]]]])/scale[] - (metric\[Delta][-i, -p$545469]*
   metric\[Delta][p$545470, p$545472]*pressure[]*timevec[h$545471]*
   PD[-p$545472][PD[-p$545470][PD[-h$545471][pertvelocityvec[LI[2], 
       p$545469]]]])/scale[] - (density[]*metric\[Delta][p$545469, p$545470]*
   metric\[Delta][p$545472, p$545473]*PD[-p$545473][
    PD[-p$545472][PD[-p$545470][PD[-p$545469][pertshearvec[LI[2], -i]]]]])/
  scale[]^2 - (metric\[Delta][p$545469, p$545470]*
   metric\[Delta][p$545472, p$545473]*pressure[]*
   PD[-p$545473][PD[-p$545472][PD[-p$545470][PD[-p$545469][
       pertshearvec[LI[2], -i]]]]])/scale[]^2
