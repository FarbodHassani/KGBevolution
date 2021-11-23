(* Created with the Wolfram Language : www.wolfram.com *)
-(metric\[Delta][p$2789488, p$2789489]*PD[-p$2789489][
    PD[-p$2789488][pertpressure[LI[1]]]]) - 
 density[]*metric\[Delta][p$2789488, p$2789489]*
  PD[-p$2789489][PD[-p$2789488][pertpsi[LI[1]]]] - 
 metric\[Delta][p$2789488, p$2789489]*pressure[]*
  PD[-p$2789489][PD[-p$2789488][pertpsi[LI[1]]]] - 
 (metric\[Delta][p$2789488, p$2789489]*primepressure[]*
   PD[-p$2789489][PD[-p$2789488][pertvelocity[LI[1]]]])/scale[] - 
 (density[]*metric\[Delta][p$2789488, p$2789489]*timevec[h$2789490]*
   PD[-p$2789489][PD[-p$2789488][PD[-h$2789490][pertvelocity[LI[1]]]]])/
  scale[] - (metric\[Delta][p$2789488, p$2789489]*pressure[]*
   timevec[h$2789490]*PD[-p$2789489][PD[-p$2789488][
     PD[-h$2789490][pertvelocity[LI[1]]]]])/scale[]
