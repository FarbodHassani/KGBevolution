(* Created with the Wolfram Language : www.wolfram.com *)
-source9[LI[2]] - metric\[Delta][p$5070528, p$5070529]*
  PD[-p$5070529][PD[-p$5070528][pertpressure[LI[2]]]] - 
 density[]*metric\[Delta][p$5070528, p$5070529]*
  PD[-p$5070529][PD[-p$5070528][pertpsi[LI[2]]]] - 
 metric\[Delta][p$5070528, p$5070529]*pressure[]*
  PD[-p$5070529][PD[-p$5070528][pertpsi[LI[2]]]] - 
 (metric\[Delta][p$5070528, p$5070529]*primepressure[]*
   PD[-p$5070529][PD[-p$5070528][pertvelocity[LI[2]]]])/scale[] - 
 (density[]*metric\[Delta][p$5070528, p$5070529]*timevec[h$5070530]*
   PD[-p$5070529][PD[-p$5070528][PD[-h$5070530][pertvelocity[LI[2]]]]])/
  scale[] - (metric\[Delta][p$5070528, p$5070529]*pressure[]*
   timevec[h$5070530]*PD[-p$5070529][PD[-p$5070528][
     PD[-h$5070530][pertvelocity[LI[2]]]]])/scale[]
