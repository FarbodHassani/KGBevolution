(* Created with the Wolfram Language : www.wolfram.com *)
density[]*hubbleC[]*metric\[Delta][p$2313797, p$2313798]*
  PD[-p$2313798][PD[-p$2313797][pertS[LI[1], -i]]] + 
 hubbleC[]*metric\[Delta][p$2313797, p$2313798]*pressure[]*
  PD[-p$2313798][PD[-p$2313797][pertS[LI[1], -i]]] + 
 metric\[Delta][p$2313797, p$2313798]*primepressure[]*
  PD[-p$2313798][PD[-p$2313797][pertS[LI[1], -i]]] + 
 density[]*metric\[Delta][p$2313797, p$2313798]*timevec[h$2313799]*
  PD[-p$2313798][PD[-p$2313797][PD[-h$2313799][pertS[LI[1], -i]]]] + 
 metric\[Delta][p$2313797, p$2313798]*pressure[]*timevec[h$2313799]*
  PD[-p$2313798][PD[-p$2313797][PD[-h$2313799][pertS[LI[1], -i]]]] - 
 (metric\[Delta][-i, -p$2313797]*metric\[Delta][p$2313798, p$2313800]*
   primepressure[]*PD[-p$2313800][PD[-p$2313798][pertvelocityvec[LI[1], 
      p$2313797]]])/scale[] - (density[]*metric\[Delta][-i, -p$2313797]*
   metric\[Delta][p$2313798, p$2313800]*timevec[h$2313799]*
   PD[-p$2313800][PD[-p$2313798][PD[-h$2313799][pertvelocityvec[LI[1], 
       p$2313797]]]])/scale[] - (metric\[Delta][-i, -p$2313797]*
   metric\[Delta][p$2313798, p$2313800]*pressure[]*timevec[h$2313799]*
   PD[-p$2313800][PD[-p$2313798][PD[-h$2313799][pertvelocityvec[LI[1], 
       p$2313797]]]])/scale[] - 
 (density[]*metric\[Delta][p$2313797, p$2313798]*metric\[Delta][p$2313800, 
    p$2313801]*PD[-p$2313801][PD[-p$2313800][PD[-p$2313798][
      PD[-p$2313797][pertshearvec[LI[1], -i]]]]])/scale[]^2 - 
 (metric\[Delta][p$2313797, p$2313798]*metric\[Delta][p$2313800, p$2313801]*
   pressure[]*PD[-p$2313801][PD[-p$2313800][PD[-p$2313798][
      PD[-p$2313797][pertshearvec[LI[1], -i]]]]])/scale[]^2
