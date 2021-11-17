(* Created with the Wolfram Language : www.wolfram.com *)
-2*hubbleC[]^3*metric\[Delta][-i, -j]*pertP[LI[1]] + 
 alphaB[]*hubbleC[]^3*metric\[Delta][-i, -j]*pertP[LI[1]] + 
 2*hubbleC[]^2*metric\[Delta][-i, -j]*pertpsi[LI[1]] - 
 alphaB[]*hubbleC[]^2*metric\[Delta][-i, -j]*pertpsi[LI[1]] - 
 (density[]*pertshearten[LI[1], -i, -j])/(2*Mpl^2) + 
 metric\[Delta][-i, -j]*pertP[LI[1]]*pprimehubbleC[] - 
 (pertshearten[LI[1], -i, -j]*pressure[])/(2*Mpl^2) + 
 (hubbleC[]^2*metric\[Delta][-i, -j]*pertP[LI[1]]*primealphaB[])/2 - 
 (hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primealphaB[])/2 + 
 alphaB[]*hubbleC[]*metric\[Delta][-i, -j]*pertP[LI[1]]*primehubbleC[] + 
 metric\[Delta][-i, -j]*pertpsi[LI[1]]*primehubbleC[] - 
 (alphaB[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primehubbleC[])/2 + 
 (density[]*hubbleC[]*metric\[Delta][-i, -j]*pertP[LI[1]]*scale[]^2)/
  (2*Mpl^2) - (metric\[Delta][-i, -j]*pertpressure[LI[1]]*scale[]^2)/
  (2*Mpl^2) - (density[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*scale[]^2)/
  (2*Mpl^2) + (hubbleC[]*metric\[Delta][-i, -j]*pertP[LI[1]]*pressure[]*
   scale[]^2)/(2*Mpl^2) - (metric\[Delta][-i, -j]*pertpsi[LI[1]]*pressure[]*
   scale[]^2)/(2*Mpl^2) + (metric\[Delta][-i, -j]*pertP[LI[1]]*
   primepressure[]*scale[]^2)/(2*Mpl^2) + 
 (hubbleC[]*timevec[h$2752640]*PD[-h$2752640][perth[LI[1], -i, -j]])/2 - 
 hubbleC[]^2*metric\[Delta][-i, -j]*timevec[h$2752640]*
  PD[-h$2752640][pertP[LI[1]]] + 
 (3*alphaB[]*hubbleC[]^2*metric\[Delta][-i, -j]*timevec[h$2752640]*
   PD[-h$2752640][pertP[LI[1]]])/2 + 
 (hubbleC[]*metric\[Delta][-i, -j]*primealphaB[]*timevec[h$2752640]*
   PD[-h$2752640][pertP[LI[1]]])/2 + metric\[Delta][-i, -j]*primehubbleC[]*
  timevec[h$2752640]*PD[-h$2752640][pertP[LI[1]]] + 
 (alphaB[]*metric\[Delta][-i, -j]*primehubbleC[]*timevec[h$2752640]*
   PD[-h$2752640][pertP[LI[1]]])/2 + 
 (density[]*metric\[Delta][-i, -j]*scale[]^2*timevec[h$2752640]*
   PD[-h$2752640][pertP[LI[1]]])/(2*Mpl^2) + 
 (metric\[Delta][-i, -j]*pressure[]*scale[]^2*timevec[h$2752640]*
   PD[-h$2752640][pertP[LI[1]]])/(2*Mpl^2) + 
 2*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$2752640]*
  PD[-h$2752640][pertphi[LI[1]]] + hubbleC[]*metric\[Delta][-i, -j]*
  timevec[h$2752640]*PD[-h$2752640][pertpsi[LI[1]]] - 
 (alphaB[]*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$2752640]*
   PD[-h$2752640][pertpsi[LI[1]]])/2 + 
 (hubbleC[]*timevec[h$2752640]*PD[-h$2752640][PD[-i][pertF[LI[1], -j]]])/2 + 
 (timevec[h$2752640]*PD[-h$2752640][PD[-i][pertS[LI[1], -j]]])/4 + 
 (hubbleC[]*timevec[h$2752640]*PD[-h$2752640][PD[-j][pertF[LI[1], -i]]])/2 + 
 (timevec[h$2752640]*PD[-h$2752640][PD[-j][pertS[LI[1], -i]]])/4 - 
 (timevec[h$2752640]*PD[-h$2752640][PD[-j][PD[-i][pertB[LI[1]]]]])/2 + 
 hubbleC[]*timevec[h$2752640]*PD[-h$2752640][PD[-j][PD[-i][pertE[LI[1]]]]] + 
 (timevec[h$2752640]*timevec[h$2752641]*PD[-h$2752641][
    PD[-h$2752640][perth[LI[1], -i, -j]]])/4 + 
 (alphaB[]*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$2752640]*
   timevec[h$2752641]*PD[-h$2752641][PD[-h$2752640][pertP[LI[1]]]])/2 + 
 metric\[Delta][-i, -j]*timevec[h$2752640]*timevec[h$2752641]*
  PD[-h$2752641][PD[-h$2752640][pertphi[LI[1]]]] + 
 (timevec[h$2752640]*timevec[h$2752641]*PD[-h$2752641][
    PD[-h$2752640][PD[-i][pertF[LI[1], -j]]]])/4 + 
 (timevec[h$2752640]*timevec[h$2752641]*PD[-h$2752641][
    PD[-h$2752640][PD[-j][pertF[LI[1], -i]]]])/4 + 
 (timevec[h$2752640]*timevec[h$2752641]*PD[-h$2752641][
    PD[-h$2752640][PD[-j][PD[-i][pertE[LI[1]]]]]])/2 + 
 (hubbleC[]*PD[-i][pertS[LI[1], -j]])/2 - 
 (density[]*PD[-i][pertshearvec[LI[1], -j]])/(2*Mpl^2) - 
 (pressure[]*PD[-i][pertshearvec[LI[1], -j]])/(2*Mpl^2) + 
 (hubbleC[]*PD[-j][pertS[LI[1], -i]])/2 - 
 (density[]*PD[-j][pertshearvec[LI[1], -i]])/(2*Mpl^2) - 
 (pressure[]*PD[-j][pertshearvec[LI[1], -i]])/(2*Mpl^2) - 
 hubbleC[]*PD[-j][PD[-i][pertB[LI[1]]]] + PD[-j][PD[-i][pertphi[LI[1]]]]/2 - 
 PD[-j][PD[-i][pertpsi[LI[1]]]]/2 - 
 (density[]*PD[-j][PD[-i][pertshear[LI[1]]]])/(4*Mpl^2) - 
 (pressure[]*PD[-j][PD[-i][pertshear[LI[1]]]])/(4*Mpl^2) + 
 hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
  PD[-p$2752643][PD[-p$2752642][pertB[LI[1]]]] - 
 (metric\[Delta][p$2752642, p$2752643]*PD[-p$2752643][
    PD[-p$2752642][perth[LI[1], -i, -j]]])/4 - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
   PD[-p$2752643][PD[-p$2752642][pertphi[LI[1]]]])/2 + 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
   PD[-p$2752643][PD[-p$2752642][pertpsi[LI[1]]]])/2 + 
 (density[]*metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
   PD[-p$2752643][PD[-p$2752642][pertshear[LI[1]]]])/(4*Mpl^2) + 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*pressure[]*
   PD[-p$2752643][PD[-p$2752642][pertshear[LI[1]]]])/(4*Mpl^2) + 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
   timevec[h$2752640]*PD[-p$2752643][PD[-p$2752642][
     PD[-h$2752640][pertB[LI[1]]]]])/2 - hubbleC[]*metric\[Delta][-i, -j]*
  metric\[Delta][p$2752642, p$2752643]*timevec[h$2752640]*
  PD[-p$2752643][PD[-p$2752642][PD[-h$2752640][pertE[LI[1]]]]] - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2752642, p$2752643]*
   timevec[h$2752640]*timevec[h$2752641]*PD[-p$2752643][
    PD[-p$2752642][PD[-h$2752641][PD[-h$2752640][pertE[LI[1]]]]]])/2
