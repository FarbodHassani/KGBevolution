(* Created with the Wolfram Language : www.wolfram.com *)
-(metric\[Delta][p$57881, p$57882]*scale[]^2*
    PD[-p$57882][PD[-p$57881][pertpressure[LI[1]]]])/(2*Mpl^2) + 
 3*hubbleC[]^2*metric\[Delta][p$57881, p$57882]*
  PD[-p$57882][PD[-p$57881][pertpsi[LI[1]]]] - 
 (density[]*metric\[Delta][p$57881, p$57882]*scale[]^2*
   PD[-p$57882][PD[-p$57881][pertpsi[LI[1]]]])/Mpl^2 - 
 (metric\[Delta][p$57881, p$57882]*pressure[]*scale[]^2*
   PD[-p$57882][PD[-p$57881][pertpsi[LI[1]]]])/Mpl^2 + 
 2*hubbleC[]*metric\[Delta][p$57881, p$57882]*timevec[h$57883]*
  PD[-p$57882][PD[-p$57881][PD[-h$57883][pertphi[LI[1]]]]] + 
 hubbleC[]*metric\[Delta][p$57881, p$57882]*timevec[h$57883]*
  PD[-p$57882][PD[-p$57881][PD[-h$57883][pertpsi[LI[1]]]]] + 
 metric\[Delta][p$57881, p$57882]*timevec[h$57883]*timevec[h$57884]*
  PD[-p$57882][PD[-p$57881][PD[-h$57884][PD[-h$57883][pertphi[LI[1]]]]]]
