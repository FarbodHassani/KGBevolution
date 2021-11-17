(* Created with the Wolfram Language : www.wolfram.com *)
-source2[LI[2]] - (2*density[]*pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281787][pertdensity[LI[1]]])/Mpl^2 + 
 (4*density[]*pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281787][pertphi[LI[1]]])/Mpl^2 + 
 (4*pertvelocityvec[LI[1], p$281787]*pressure[]*scale[]*
   PD[-p$281787][pertphi[LI[1]]])/Mpl^2 - 
 (2*pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281787][pertpressure[LI[1]]])/Mpl^2 - 
 (2*density[]*pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281787][pertpsi[LI[1]]])/Mpl^2 - 
 (2*pertvelocityvec[LI[1], p$281787]*pressure[]*scale[]*
   PD[-p$281787][pertpsi[LI[1]]])/Mpl^2 + 
 (2*density[]*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*
   scale[]^2*PD[-p$281788][pertdensity[LI[1]]])/Mpl^2 - 
 (4*density[]*metric\[Delta][p$281787, p$281788]*scale[]^2*
   PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][pertphi[LI[1]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281788]*pressure[]*scale[]^2*
   PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][pertphi[LI[1]]])/Mpl^2 + 
 (2*metric\[Delta][p$281787, p$281788]*scale[]^2*PD[-p$281787][pertB[LI[1]]]*
   PD[-p$281788][pertpressure[LI[1]]])/Mpl^2 + 
 12*hubbleC[]^2*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*
  PD[-p$281788][pertpsi[LI[1]]] - 
 (4*density[]*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*
   scale[]^2*PD[-p$281788][pertpsi[LI[1]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*pressure[]*
   scale[]^2*PD[-p$281788][pertpsi[LI[1]]])/Mpl^2 - 
 12*hubbleC[]^2*metric\[Delta][p$281787, p$281788]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][pertpsi[LI[1]]] + 
 (6*density[]*metric\[Delta][p$281787, p$281788]*scale[]^2*
   PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][pertpsi[LI[1]]])/Mpl^2 + 
 (6*metric\[Delta][p$281787, p$281788]*pressure[]*scale[]^2*
   PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][pertpsi[LI[1]]])/Mpl^2 + 
 8*hubbleC[]*metric\[Delta][p$281787, p$281788]*PD[-p$281787][pertpsi[LI[1]]]*
  PD[-p$281788][pertpsi[LI[1]]] + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*timevec[h$281791]*
  PD[-h$281791][perth[LI[1], -p$281789, -p$281790]]*
  PD[-p$281788][pertS[LI[1], -p$281787]] + metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*timevec[h$281791]*timevec[h$281792]*
  PD[-h$281792][PD[-h$281791][perth[LI[1], -p$281789, -p$281790]]]*
  PD[-p$281788][pertS[LI[1], -p$281787]] - 
 (2*density[]*metric\[Delta][p$281787, p$281788]*scale[]*
   PD[-p$281787][pertdensity[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/
  Mpl^2 + (4*density[]*metric\[Delta][p$281787, p$281788]*scale[]*
   PD[-p$281787][pertphi[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/Mpl^2 + 
 (4*metric\[Delta][p$281787, p$281788]*pressure[]*scale[]*
   PD[-p$281787][pertphi[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*scale[]*
   PD[-p$281787][pertpressure[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/
  Mpl^2 - (2*density[]*metric\[Delta][p$281787, p$281788]*scale[]*
   PD[-p$281787][pertpsi[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*pressure[]*scale[]*
   PD[-p$281787][pertpsi[LI[1]]]*PD[-p$281788][pertvelocity[LI[1]]])/Mpl^2 - 
 (2*density[]*metric\[Delta][p$281788, p$281789]*
   perth[LI[1], -p$281787, -p$281789]*scale[]*
   PD[-p$281788][pertvelocityvec[LI[1], p$281787]])/Mpl^2 - 
 (2*metric\[Delta][p$281788, p$281789]*perth[LI[1], -p$281787, -p$281789]*
   pressure[]*scale[]*PD[-p$281788][pertvelocityvec[LI[1], p$281787]])/
  Mpl^2 + 8*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  pertS[LI[1], -p$281787]*timevec[h$281791]*
  PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]] - 
 8*hubbleC[]*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]] - 
 16*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-p$281787][pertphi[LI[1]]]*PD[-p$281788][
   PD[-h$281791][pertphi[LI[1]]]] + 4*metric\[Delta][p$281787, p$281788]*
  timevec[h$281791]*PD[-p$281787][pertpsi[LI[1]]]*
  PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]] + 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*
  timevec[h$281791]*PD[-p$281788][PD[-h$281791][pertpsi[LI[1]]]] - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][PD[-h$281791][pertpsi[LI[1]]]] + 
 4*metric\[Delta][p$281787, p$281788]*pertS[LI[1], -p$281787]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281788][
   PD[-h$281792][PD[-h$281791][pertphi[LI[1]]]]] - 
 4*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*timevec[h$281792]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281788][
   PD[-h$281792][PD[-h$281791][pertphi[LI[1]]]]] - 
 12*hubbleC[]^2*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*
  PD[-p$281788][PD[-p$281787][pertB[LI[1]]]] - 
 (4*density[]*metric\[Delta][p$281787, p$281788]*pertphi[LI[1]]*scale[]^2*
   PD[-p$281788][PD[-p$281787][pertB[LI[1]]]])/Mpl^2 + 
 (2*metric\[Delta][p$281787, p$281788]*pertpressure[LI[1]]*scale[]^2*
   PD[-p$281788][PD[-p$281787][pertB[LI[1]]]])/Mpl^2 + 
 (6*density[]*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*scale[]^2*
   PD[-p$281788][PD[-p$281787][pertB[LI[1]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281788]*pertphi[LI[1]]*pressure[]*scale[]^2*
   PD[-p$281788][PD[-p$281787][pertB[LI[1]]]])/Mpl^2 + 
 (6*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*pressure[]*scale[]^2*
   PD[-p$281788][PD[-p$281787][pertB[LI[1]]]])/Mpl^2 - 
 8*hubbleC[]*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-h$281791][pertphi[LI[1]]]*PD[-p$281788][PD[-p$281787][pertB[LI[1]]]] - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-h$281791][pertpsi[LI[1]]]*PD[-p$281788][PD[-p$281787][pertB[LI[1]]]] - 
 4*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*timevec[h$281792]*
  PD[-h$281792][PD[-h$281791][pertphi[LI[1]]]]*
  PD[-p$281788][PD[-p$281787][pertB[LI[1]]]] - 
 8*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-h$281791][pertphi[LI[1]]]*PD[-p$281788][
   PD[-p$281787][pertphi[LI[1]]]] + 
 8*hubbleC[]*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*
  PD[-p$281788][PD[-p$281787][pertpsi[LI[1]]]] + 
 4*metric\[Delta][p$281787, p$281788]*timevec[h$281791]*
  PD[-h$281791][pertphi[LI[1]]]*PD[-p$281788][
   PD[-p$281787][pertpsi[LI[1]]]] - 
 (2*density[]*metric\[Delta][p$281787, p$281788]*pertdensity[LI[1]]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 + 
 (4*density[]*metric\[Delta][p$281787, p$281788]*pertphi[LI[1]]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*pertpressure[LI[1]]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 - 
 (2*density[]*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 + 
 (4*metric\[Delta][p$281787, p$281788]*pertphi[LI[1]]*pressure[]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*pertpsi[LI[1]]*pressure[]*scale[]*
   PD[-p$281788][PD[-p$281787][pertvelocity[LI[1]]]])/Mpl^2 - 
 8*metric\[Delta][p$281787, p$281788]*pertphi[LI[1]]*timevec[h$281791]*
  PD[-p$281788][PD[-p$281787][PD[-h$281791][pertphi[LI[1]]]]] - 
 (2*density[]*metric\[Delta][p$281788, p$281789]*scale[]*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281789][pertvelocityvec[LI[1], p$281787]])/Mpl^2 - 
 (2*metric\[Delta][p$281788, p$281789]*pressure[]*scale[]*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281789][pertvelocityvec[LI[1], p$281787]])/Mpl^2 - 
 (2*density[]*metric\[Delta][p$281787, p$281789]*scale[]*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281789][pertvelocityvec[LI[1], p$281788]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281789]*pressure[]*scale[]*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281789][pertvelocityvec[LI[1], p$281788]])/Mpl^2 + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*timevec[h$281791]*
  PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281790]]] + 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281789][PD[-h$281791][pertS[LI[1], -p$281790]]] + 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281789][PD[-h$281792][PD[-h$281791][pertF[LI[1], -p$281790]]]] - 
 2*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*
  PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]*
  PD[-p$281789][PD[-p$281787][pertB[LI[1]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-h$281792][
   PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]]*
  PD[-p$281789][PD[-p$281787][pertB[LI[1]]]] - 
 (4*density[]*metric\[Delta][p$281788, p$281789]*scale[]*
   PD[-p$281788][pertvelocityvec[LI[1], p$281787]]*
   PD[-p$281789][PD[-p$281787][pertE[LI[1]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281788, p$281789]*pressure[]*scale[]*
   PD[-p$281788][pertvelocityvec[LI[1], p$281787]]*
   PD[-p$281789][PD[-p$281787][pertE[LI[1]]]])/Mpl^2 + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]*
  PD[-p$281789][PD[-p$281787][pertphi[LI[1]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]*
  PD[-p$281789][PD[-p$281787][pertpsi[LI[1]]]] - 
 (2*density[]*metric\[Delta][p$281788, p$281789]*
   pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281789][PD[-p$281788][pertF[LI[1], -p$281787]]])/Mpl^2 - 
 (2*metric\[Delta][p$281788, p$281789]*pertvelocityvec[LI[1], p$281787]*
   pressure[]*scale[]*PD[-p$281789][PD[-p$281788][pertF[LI[1], -p$281787]]])/
  Mpl^2 - (4*density[]*metric\[Delta][p$281788, p$281789]*
   pertvelocityvec[LI[1], p$281787]*scale[]*
   PD[-p$281789][PD[-p$281788][PD[-p$281787][pertE[LI[1]]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281788, p$281789]*pertvelocityvec[LI[1], p$281787]*
   pressure[]*scale[]*PD[-p$281789][PD[-p$281788][
     PD[-p$281787][pertE[LI[1]]]]])/Mpl^2 + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281790]*
  metric\[Delta][p$281788, p$281789]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][pertS[LI[1], -p$281789]] - 
 2*hubbleC[]*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][pertS[LI[1], -p$281789]] - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281788]]] + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*timevec[h$281791]*
  PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281789]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-h$281791][pertS[LI[1], -p$281788]]] + 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-h$281791][pertS[LI[1], -p$281789]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281789][
   PD[-p$281787][pertB[LI[1]]]]*PD[-p$281790][
   PD[-h$281792][PD[-h$281791][pertF[LI[1], -p$281788]]]] + 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-h$281792][PD[-h$281791][pertF[LI[1], -p$281789]]]] + 
 (2*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*perth[LI[1], -p$281787, -p$281789]*
   scale[]^2*PD[-p$281790][PD[-p$281788][pertB[LI[1]]]])/Mpl^2 + 
 (2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   perth[LI[1], -p$281787, -p$281789]*pressure[]*scale[]^2*
   PD[-p$281790][PD[-p$281788][pertB[LI[1]]]])/Mpl^2 + 
 (4*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]^2*
   PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
   PD[-p$281790][PD[-p$281788][pertE[LI[1]]]])/Mpl^2 + 
 (4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]^2*PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
   PD[-p$281790][PD[-p$281788][pertE[LI[1]]]])/Mpl^2 + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281787]]]*
  PD[-p$281790][PD[-p$281788][pertphi[LI[1]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281787]]]*
  PD[-p$281790][PD[-p$281788][pertpsi[LI[1]]]] - 
 (2*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*perth[LI[1], -p$281787, -p$281789]*
   scale[]*PD[-p$281790][PD[-p$281788][pertvelocity[LI[1]]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   perth[LI[1], -p$281787, -p$281789]*pressure[]*scale[]*
   PD[-p$281790][PD[-p$281788][pertvelocity[LI[1]]]])/Mpl^2 - 
 (4*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]*
   PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]*
   PD[-p$281790][PD[-p$281788][pertvelocity[LI[1]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]*PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]*
   PD[-p$281790][PD[-p$281788][pertvelocity[LI[1]]]])/Mpl^2 + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertB[LI[1]]]]] - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-p$281787][pertphi[LI[1]]]]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-p$281787][pertpsi[LI[1]]]]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  perth[LI[1], -p$281787, -p$281789]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]]] - 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281789][
   PD[-p$281787][pertB[LI[1]]]]*PD[-p$281790][
   PD[-p$281788][PD[-h$281792][PD[-h$281791][pertE[LI[1]]]]]] + 
 (4*density[]*metric\[Delta][p$281787, p$281789]*
   metric\[Delta][p$281788, p$281790]*scale[]^2*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281790][PD[-p$281789][pertB[LI[1]]]])/Mpl^2 + 
 (4*metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
   pressure[]*scale[]^2*PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281790][PD[-p$281789][pertB[LI[1]]]])/Mpl^2 - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*PD[-p$281788][
   PD[-p$281787][pertB[LI[1]]]]*PD[-p$281790][PD[-p$281789][pertB[LI[1]]]] - 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281788][PD[-h$281791][pertphi[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281787]]] + 
 (2*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]^2*PD[-p$281787][pertB[LI[1]]]*
   PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281788]]])/Mpl^2 + 
 (2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]^2*PD[-p$281787][pertB[LI[1]]]*
   PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281788]]])/Mpl^2 - 
 (2*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]*
   PD[-p$281787][pertvelocity[LI[1]]]*PD[-p$281790][
    PD[-p$281789][pertF[LI[1], -p$281788]]])/Mpl^2 - 
 (2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]*PD[-p$281787][pertvelocity[LI[1]]]*
   PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281788]]])/Mpl^2 + 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  PD[-p$281788][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][pertphi[LI[1]]]] - 
 2*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*pertS[LI[1], -p$281787]*
  PD[-p$281790][PD[-p$281789][pertS[LI[1], -p$281788]]] + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*PD[-p$281787][pertB[LI[1]]]*
  PD[-p$281790][PD[-p$281789][pertS[LI[1], -p$281788]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  PD[-p$281787][pertphi[LI[1]]]*PD[-p$281790][
   PD[-p$281789][pertS[LI[1], -p$281788]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  PD[-p$281787][pertpsi[LI[1]]]*PD[-p$281790][
   PD[-p$281789][pertS[LI[1], -p$281788]]] - 
 (4*density[]*metric\[Delta][p$281787, p$281789]*
   metric\[Delta][p$281788, p$281790]*scale[]*
   PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281790][PD[-p$281789][pertvelocity[LI[1]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
   pressure[]*scale[]*PD[-p$281788][pertF[LI[1], -p$281787]]*
   PD[-p$281790][PD[-p$281789][pertvelocity[LI[1]]]])/Mpl^2 - 
 2*metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertB[LI[1]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281788][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertB[LI[1]]]]] + 
 4*hubbleC[]*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*timevec[h$281791]*
  PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertE[LI[1]]]]] + 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*
  PD[-p$281788][PD[-p$281787][pertB[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertE[LI[1]]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281788][PD[-p$281787][pertphi[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertE[LI[1]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281788][PD[-p$281787][pertpsi[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertE[LI[1]]]]] + 
 2*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*pertS[LI[1], -p$281787]*
  timevec[h$281791]*PD[-p$281790][PD[-p$281789][
    PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 2*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281790][
   PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281787][pertphi[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281788]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281787][pertpsi[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 4*metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertphi[LI[1]]]]] + 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281788][PD[-p$281787][pertE[LI[1]]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertphi[LI[1]]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  pertS[LI[1], -p$281787]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertS[LI[1], -p$281788]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281787][pertB[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281791][pertS[LI[1], -p$281788]]]] + 
 2*metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281788][pertS[LI[1], -p$281787]]*
  PD[-p$281790][PD[-p$281789][PD[-h$281792][PD[-h$281791][pertE[LI[1]]]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*timevec[h$281792]*PD[-p$281788][
   PD[-p$281787][pertB[LI[1]]]]*PD[-p$281790][
   PD[-p$281789][PD[-h$281792][PD[-h$281791][pertE[LI[1]]]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  pertS[LI[1], -p$281787]*timevec[h$281791]*timevec[h$281792]*
  PD[-p$281790][PD[-p$281789][PD[-h$281792][PD[-h$281791][
      pertF[LI[1], -p$281788]]]]] - metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*timevec[h$281791]*timevec[h$281792]*
  PD[-p$281787][pertB[LI[1]]]*PD[-p$281790][
   PD[-p$281789][PD[-h$281792][PD[-h$281791][pertF[LI[1], -p$281788]]]]] + 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*pertS[LI[1], -p$281787]*
  PD[-p$281790][PD[-p$281789][PD[-p$281788][pertB[LI[1]]]]] - 
 4*hubbleC[]*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*PD[-p$281787][pertB[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-p$281788][pertB[LI[1]]]]] + 
 (4*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]^2*PD[-p$281787][pertB[LI[1]]]*
   PD[-p$281790][PD[-p$281789][PD[-p$281788][pertE[LI[1]]]]])/Mpl^2 + 
 (4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]^2*PD[-p$281787][pertB[LI[1]]]*
   PD[-p$281790][PD[-p$281789][PD[-p$281788][pertE[LI[1]]]]])/Mpl^2 - 
 (4*density[]*metric\[Delta][p$281787, p$281788]*
   metric\[Delta][p$281789, p$281790]*scale[]*
   PD[-p$281787][pertvelocity[LI[1]]]*PD[-p$281790][
    PD[-p$281789][PD[-p$281788][pertE[LI[1]]]]])/Mpl^2 - 
 (4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   pressure[]*scale[]*PD[-p$281787][pertvelocity[LI[1]]]*
   PD[-p$281790][PD[-p$281789][PD[-p$281788][pertE[LI[1]]]]])/Mpl^2 - 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  pertS[LI[1], -p$281787]*PD[-p$281790][PD[-p$281789][
    PD[-p$281788][pertphi[LI[1]]]]] + 4*metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*PD[-p$281787][pertB[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-p$281788][pertphi[LI[1]]]]] + 
 4*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  timevec[h$281791]*PD[-p$281787][pertphi[LI[1]]]*
  PD[-p$281790][PD[-p$281789][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281790][PD[-h$281791][perth[LI[1], -p$281788, -p$281794]]]*
  PD[-p$281793][perth[LI[1], -p$281787, -p$281789]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertF[LI[1], -p$281794]]]]*
  PD[-p$281793][perth[LI[1], -p$281787, -p$281789]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281788][PD[-h$281791][perth[LI[1], -p$281790, -p$281794]]]*
  PD[-p$281793][PD[-p$281789][pertF[LI[1], -p$281787]]] - 
 (3*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
   PD[-p$281793][perth[LI[1], -p$281787, -p$281789]]*
   PD[-p$281794][PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]])/2 - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][perth[LI[1], -p$281787, -p$281789]]*
  PD[-p$281794][PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]]*
  PD[-p$281794][PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281788]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-h$281791][perth[LI[1], -p$281787, -p$281789]]]*
  PD[-p$281794][PD[-p$281790][PD[-p$281788][pertE[LI[1]]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][perth[LI[1], -p$281787, -p$281789]]*
  PD[-p$281794][PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281793][PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]]*
  PD[-p$281794][PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]*
  PD[-p$281794][PD[-p$281793][perth[LI[1], -p$281787, -p$281789]]] - 
 (metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
   metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
   PD[-h$281791][perth[LI[1], -p$281787, -p$281789]]*
   PD[-p$281794][PD[-p$281793][perth[LI[1], -p$281788, -p$281790]]])/2 - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281789][PD[-h$281791][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281793][perth[LI[1], -p$281788, -p$281790]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*PD[-p$281790][
   PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281793][pertS[LI[1], -p$281788]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281793][PD[-h$281791][pertF[LI[1], -p$281788]]]] + 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281788][PD[-p$281787][pertE[LI[1]]]]]*
  PD[-p$281794][PD[-p$281793][PD[-h$281791][pertF[LI[1], -p$281790]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*perth[LI[1], -p$281787, -p$281789]*
  timevec[h$281791]*PD[-p$281794][PD[-p$281793][
    PD[-h$281791][perth[LI[1], -p$281788, -p$281790]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]*
  PD[-p$281794][PD[-p$281793][PD[-h$281791][perth[LI[1], -p$281788, 
      -p$281790]]]] - 2*metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*metric\[Delta][p$281793, p$281794]*
  timevec[h$281791]*PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281794][PD[-p$281793][PD[-h$281791][perth[LI[1], -p$281789, 
      -p$281790]]]] - metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*metric\[Delta][p$281793, p$281794]*
  PD[-p$281790][PD[-p$281789][pertS[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281788][pertE[LI[1]]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281790][PD[-p$281789][pertF[LI[1], -p$281787]]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281788][PD[-h$281791][pertE[LI[1]]]]]] - 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  metric\[Delta][p$281793, p$281794]*PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281789][pertS[LI[1], -p$281790]]]] - 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281789][PD[-h$281791][
      pertF[LI[1], -p$281790]]]]] - metric\[Delta][p$281787, p$281788]*
  metric\[Delta][p$281789, p$281790]*metric\[Delta][p$281793, p$281794]*
  perth[LI[1], -p$281787, -p$281789]*PD[-p$281794][
   PD[-p$281793][PD[-p$281790][pertS[LI[1], -p$281788]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*PD[-p$281789][
   PD[-p$281787][pertE[LI[1]]]]*PD[-p$281794][
   PD[-p$281793][PD[-p$281790][pertS[LI[1], -p$281788]]]] - 
 metric\[Delta][p$281787, p$281789]*metric\[Delta][p$281788, p$281790]*
  metric\[Delta][p$281793, p$281794]*PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281790][pertS[LI[1], -p$281789]]]] + 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281788][PD[-p$281787][pertE[LI[1]]]]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281790][PD[-h$281791][pertE[LI[1]]]]]] - 
 metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*perth[LI[1], -p$281787, -p$281789]*
  timevec[h$281791]*PD[-p$281794][PD[-p$281793][
    PD[-p$281790][PD[-h$281791][pertF[LI[1], -p$281788]]]]] - 
 2*metric\[Delta][p$281787, p$281788]*metric\[Delta][p$281789, p$281790]*
  metric\[Delta][p$281793, p$281794]*timevec[h$281791]*
  PD[-p$281789][PD[-p$281787][pertE[LI[1]]]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281790][PD[-h$281791][
      pertF[LI[1], -p$281788]]]]] - metric\[Delta][p$281787, p$281789]*
  metric\[Delta][p$281788, p$281790]*metric\[Delta][p$281793, p$281794]*
  timevec[h$281791]*PD[-p$281788][pertF[LI[1], -p$281787]]*
  PD[-p$281794][PD[-p$281793][PD[-p$281790][PD[-h$281791][
      pertF[LI[1], -p$281789]]]]]
