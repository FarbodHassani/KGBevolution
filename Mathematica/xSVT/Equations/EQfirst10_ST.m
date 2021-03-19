(* Created with the Wolfram Language : www.wolfram.com *)
density[]*hubbleC[]*metric\[Delta][p$2925312, p$2925313]*
  PD[-p$2925313][PD[-p$2925312][pertS[LI[1], -i]]] + 
 hubbleC[]*metric\[Delta][p$2925312, p$2925313]*pressure[]*
  PD[-p$2925313][PD[-p$2925312][pertS[LI[1], -i]]] + 
 metric\[Delta][p$2925312, p$2925313]*primepressure[]*
  PD[-p$2925313][PD[-p$2925312][pertS[LI[1], -i]]] + 
 density[]*metric\[Delta][p$2925312, p$2925313]*timevec[h$2925314]*
  PD[-p$2925313][PD[-p$2925312][PD[-h$2925314][pertS[LI[1], -i]]]] + 
 metric\[Delta][p$2925312, p$2925313]*pressure[]*timevec[h$2925314]*
  PD[-p$2925313][PD[-p$2925312][PD[-h$2925314][pertS[LI[1], -i]]]] - 
 (metric\[Delta][-i, -p$2925312]*metric\[Delta][p$2925313, p$2925315]*
   primepressure[]*PD[-p$2925315][PD[-p$2925313][pertvelocityvec[LI[1], 
      p$2925312]]])/scale[] - (density[]*metric\[Delta][-i, -p$2925312]*
   metric\[Delta][p$2925313, p$2925315]*timevec[h$2925314]*
   PD[-p$2925315][PD[-p$2925313][PD[-h$2925314][pertvelocityvec[LI[1], 
       p$2925312]]]])/scale[] - (metric\[Delta][-i, -p$2925312]*
   metric\[Delta][p$2925313, p$2925315]*pressure[]*timevec[h$2925314]*
   PD[-p$2925315][PD[-p$2925313][PD[-h$2925314][pertvelocityvec[LI[1], 
       p$2925312]]]])/scale[] - 
 (density[]*metric\[Delta][p$2925312, p$2925313]*metric\[Delta][p$2925315, 
    p$2925316]*PD[-p$2925316][PD[-p$2925315][PD[-p$2925313][
      PD[-p$2925312][pertshearvec[LI[1], -i]]]]])/scale[]^2 - 
 (metric\[Delta][p$2925312, p$2925313]*metric\[Delta][p$2925315, p$2925316]*
   pressure[]*PD[-p$2925316][PD[-p$2925315][PD[-p$2925313][
      PD[-p$2925312][pertshearvec[LI[1], -i]]]]])/scale[]^2
