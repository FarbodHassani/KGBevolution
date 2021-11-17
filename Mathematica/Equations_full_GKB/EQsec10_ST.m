(* Created with the Wolfram Language : www.wolfram.com *)
-source10[LI[2], -i] + density[]*hubbleC[]*metric\[Delta][p$378961, p$378962]*
  PD[-p$378962][PD[-p$378961][pertS[LI[2], -i]]] + 
 hubbleC[]*metric\[Delta][p$378961, p$378962]*pressure[]*
  PD[-p$378962][PD[-p$378961][pertS[LI[2], -i]]] + 
 metric\[Delta][p$378961, p$378962]*primepressure[]*
  PD[-p$378962][PD[-p$378961][pertS[LI[2], -i]]] + 
 density[]*metric\[Delta][p$378961, p$378962]*timevec[h$378963]*
  PD[-p$378962][PD[-p$378961][PD[-h$378963][pertS[LI[2], -i]]]] + 
 metric\[Delta][p$378961, p$378962]*pressure[]*timevec[h$378963]*
  PD[-p$378962][PD[-p$378961][PD[-h$378963][pertS[LI[2], -i]]]] - 
 (metric\[Delta][-i, -p$378961]*metric\[Delta][p$378962, p$378964]*
   primepressure[]*PD[-p$378964][PD[-p$378962][pertvelocityvec[LI[2], 
      p$378961]]])/scale[] - (density[]*metric\[Delta][-i, -p$378961]*
   metric\[Delta][p$378962, p$378964]*timevec[h$378963]*
   PD[-p$378964][PD[-p$378962][PD[-h$378963][pertvelocityvec[LI[2], 
       p$378961]]]])/scale[] - (metric\[Delta][-i, -p$378961]*
   metric\[Delta][p$378962, p$378964]*pressure[]*timevec[h$378963]*
   PD[-p$378964][PD[-p$378962][PD[-h$378963][pertvelocityvec[LI[2], 
       p$378961]]]])/scale[] - (density[]*metric\[Delta][p$378961, p$378962]*
   metric\[Delta][p$378964, p$378965]*PD[-p$378965][
    PD[-p$378964][PD[-p$378962][PD[-p$378961][pertshearvec[LI[2], -i]]]]])/
  scale[]^2 - (metric\[Delta][p$378961, p$378962]*
   metric\[Delta][p$378964, p$378965]*pressure[]*
   PD[-p$378965][PD[-p$378964][PD[-p$378962][PD[-p$378961][
       pertshearvec[LI[2], -i]]]]])/scale[]^2
