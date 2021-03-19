(* Created with the Wolfram Language : www.wolfram.com *)
-source10[LI[2], -i] + density[]*hubbleC[]*metric\[Delta][p$374789, p$374790]*
  PD[-p$374790][PD[-p$374789][pertS[LI[2], -i]]] + 
 hubbleC[]*metric\[Delta][p$374789, p$374790]*pressure[]*
  PD[-p$374790][PD[-p$374789][pertS[LI[2], -i]]] + 
 metric\[Delta][p$374789, p$374790]*primepressure[]*
  PD[-p$374790][PD[-p$374789][pertS[LI[2], -i]]] + 
 density[]*metric\[Delta][p$374789, p$374790]*timevec[h$374791]*
  PD[-p$374790][PD[-p$374789][PD[-h$374791][pertS[LI[2], -i]]]] + 
 metric\[Delta][p$374789, p$374790]*pressure[]*timevec[h$374791]*
  PD[-p$374790][PD[-p$374789][PD[-h$374791][pertS[LI[2], -i]]]] - 
 (metric\[Delta][-i, -p$374789]*metric\[Delta][p$374790, p$374792]*
   primepressure[]*PD[-p$374792][PD[-p$374790][pertvelocityvec[LI[2], 
      p$374789]]])/scale[] - (density[]*metric\[Delta][-i, -p$374789]*
   metric\[Delta][p$374790, p$374792]*timevec[h$374791]*
   PD[-p$374792][PD[-p$374790][PD[-h$374791][pertvelocityvec[LI[2], 
       p$374789]]]])/scale[] - (metric\[Delta][-i, -p$374789]*
   metric\[Delta][p$374790, p$374792]*pressure[]*timevec[h$374791]*
   PD[-p$374792][PD[-p$374790][PD[-h$374791][pertvelocityvec[LI[2], 
       p$374789]]]])/scale[] - (density[]*metric\[Delta][p$374789, p$374790]*
   metric\[Delta][p$374792, p$374793]*PD[-p$374793][
    PD[-p$374792][PD[-p$374790][PD[-p$374789][pertshearvec[LI[2], -i]]]]])/
  scale[]^2 - (metric\[Delta][p$374789, p$374790]*
   metric\[Delta][p$374792, p$374793]*pressure[]*
   PD[-p$374793][PD[-p$374792][PD[-p$374790][PD[-p$374789][
       pertshearvec[LI[2], -i]]]]])/scale[]^2
