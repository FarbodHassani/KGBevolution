(* Created with the Wolfram Language : www.wolfram.com *)
density[]*hubbleC[]*metric\[Delta][p$75734, p$75735]*
  PD[-p$75735][PD[-p$75734][pertS[LI[1], -i]]] + 
 hubbleC[]*metric\[Delta][p$75734, p$75735]*pressure[]*
  PD[-p$75735][PD[-p$75734][pertS[LI[1], -i]]] + 
 metric\[Delta][p$75734, p$75735]*primepressure[]*
  PD[-p$75735][PD[-p$75734][pertS[LI[1], -i]]] + 
 density[]*metric\[Delta][p$75734, p$75735]*timevec[h$75736]*
  PD[-p$75735][PD[-p$75734][PD[-h$75736][pertS[LI[1], -i]]]] + 
 metric\[Delta][p$75734, p$75735]*pressure[]*timevec[h$75736]*
  PD[-p$75735][PD[-p$75734][PD[-h$75736][pertS[LI[1], -i]]]] - 
 (metric\[Delta][-i, -p$75734]*metric\[Delta][p$75735, p$75737]*
   primepressure[]*PD[-p$75737][PD[-p$75735][pertvelocityvec[LI[1], 
      p$75734]]])/scale[] - (density[]*metric\[Delta][-i, -p$75734]*
   metric\[Delta][p$75735, p$75737]*timevec[h$75736]*
   PD[-p$75737][PD[-p$75735][PD[-h$75736][pertvelocityvec[LI[1], p$75734]]]])/
  scale[] - (metric\[Delta][-i, -p$75734]*metric\[Delta][p$75735, p$75737]*
   pressure[]*timevec[h$75736]*PD[-p$75737][
    PD[-p$75735][PD[-h$75736][pertvelocityvec[LI[1], p$75734]]]])/scale[] - 
 (density[]*metric\[Delta][p$75734, p$75735]*metric\[Delta][p$75737, p$75738]*
   PD[-p$75738][PD[-p$75737][PD[-p$75735][PD[-p$75734][
       pertshearvec[LI[1], -i]]]]])/scale[]^2 - 
 (metric\[Delta][p$75734, p$75735]*metric\[Delta][p$75737, p$75738]*
   pressure[]*PD[-p$75738][PD[-p$75737][PD[-p$75735][
      PD[-p$75734][pertshearvec[LI[1], -i]]]]])/scale[]^2
