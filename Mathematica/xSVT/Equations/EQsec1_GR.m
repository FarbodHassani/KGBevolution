(* Created with the Wolfram Language : www.wolfram.com *)
3*hubbleC[]^2*pertpsi[LI[2]] + (density[]*pertdensity[LI[2]]*scale[]^2)/
  (2*Mpl^2) - source1[LI[2]] + 3*hubbleC[]*timevec[h$154225]*
  PD[-h$154225][pertphi[LI[2]]] + 
 hubbleC[]*metric\[Delta][p$154221, p$154222]*
  PD[-p$154222][PD[-p$154221][pertB[LI[2]]]] - 
 metric\[Delta][p$154221, p$154222]*PD[-p$154222][
   PD[-p$154221][pertphi[LI[2]]]] - 
 hubbleC[]*metric\[Delta][p$154221, p$154222]*timevec[h$154225]*
  PD[-p$154222][PD[-p$154221][PD[-h$154225][pertE[LI[2]]]]]
