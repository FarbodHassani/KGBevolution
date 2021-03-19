(* Created with the Wolfram Language : www.wolfram.com *)
-(metric\[Delta][p$2924447, p$2924448]*PD[-p$2924448][
    PD[-p$2924447][pertpressure[LI[1]]]]) - 
 density[]*metric\[Delta][p$2924447, p$2924448]*
  PD[-p$2924448][PD[-p$2924447][pertpsi[LI[1]]]] - 
 metric\[Delta][p$2924447, p$2924448]*pressure[]*
  PD[-p$2924448][PD[-p$2924447][pertpsi[LI[1]]]] - 
 (metric\[Delta][p$2924447, p$2924448]*primepressure[]*
   PD[-p$2924448][PD[-p$2924447][pertvelocity[LI[1]]]])/scale[] - 
 (density[]*metric\[Delta][p$2924447, p$2924448]*timevec[h$2924449]*
   PD[-p$2924448][PD[-p$2924447][PD[-h$2924449][pertvelocity[LI[1]]]]])/
  scale[] - (metric\[Delta][p$2924447, p$2924448]*pressure[]*
   timevec[h$2924449]*PD[-p$2924448][PD[-p$2924447][
     PD[-h$2924449][pertvelocity[LI[1]]]]])/scale[]
