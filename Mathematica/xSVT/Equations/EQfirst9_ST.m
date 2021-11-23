(* Created with the Wolfram Language : www.wolfram.com *)
-(metric\[Delta][p$2312931, p$2312932]*PD[-p$2312932][
    PD[-p$2312931][pertpressure[LI[1]]]]) - 
 density[]*metric\[Delta][p$2312931, p$2312932]*
  PD[-p$2312932][PD[-p$2312931][pertpsi[LI[1]]]] - 
 metric\[Delta][p$2312931, p$2312932]*pressure[]*
  PD[-p$2312932][PD[-p$2312931][pertpsi[LI[1]]]] - 
 (metric\[Delta][p$2312931, p$2312932]*primepressure[]*
   PD[-p$2312932][PD[-p$2312931][pertvelocity[LI[1]]]])/scale[] - 
 (density[]*metric\[Delta][p$2312931, p$2312932]*timevec[h$2312933]*
   PD[-p$2312932][PD[-p$2312931][PD[-h$2312933][pertvelocity[LI[1]]]]])/
  scale[] - (metric\[Delta][p$2312931, p$2312932]*pressure[]*
   timevec[h$2312933]*PD[-p$2312932][PD[-p$2312931][
     PD[-h$2312933][pertvelocity[LI[1]]]]])/scale[]
