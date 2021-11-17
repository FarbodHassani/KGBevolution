(* Created with the Wolfram Language : www.wolfram.com *)
Mpl^2*hubbleC[]^2*pertS[LI[1], -i] + 2*Mpl^2*pertS[LI[1], -i]*
  primehubbleC[] + density[]*metric\[Delta][-i, -p$2738779]*
  pertvelocityvec[LI[1], p$2738779]*scale[] + metric\[Delta][-i, -p$2738779]*
  pertvelocityvec[LI[1], p$2738779]*pressure[]*scale[] - 
 density[]*pertS[LI[1], -i]*scale[]^2 + G2fun[scalar[], X[]]*pertS[LI[1], -i]*
  scale[]^2 + 2*Mpl^2*timevec[h$2738665]*
  PD[-h$2738665][PD[-i][pertphi[LI[1]]]] - 
 Mpl^2*hubbleC[]^2*PD[-i][pertB[LI[1]]] - 2*Mpl^2*primehubbleC[]*
  PD[-i][pertB[LI[1]]] - G2fun[scalar[], X[]]*scale[]^2*
  PD[-i][pertB[LI[1]]] - pressure[]*scale[]^2*PD[-i][pertB[LI[1]]] + 
 2*Mpl^2*hubbleC[]*PD[-i][pertpsi[LI[1]]] + 
 density[]*scale[]*PD[-i][pertvelocity[LI[1]]] + 
 pressure[]*scale[]*PD[-i][pertvelocity[LI[1]]] + 
 (Mpl^2*metric\[Delta][p$2738666, p$2738667]*
   PD[-p$2738667][PD[-p$2738666][pertS[LI[1], -i]]])/2 + 
 (Mpl^2*metric\[Delta][p$2738666, p$2738667]*timevec[h$2738665]*
   PD[-p$2738667][PD[-p$2738666][PD[-h$2738665][pertF[LI[1], -i]]]])/2 - 
 primescalar[]*PD[-i][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], 
   X[]] - (pertS[LI[1], -i]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (hubbleC[]*pertS[LI[1], -i]*primescalar[]^3*Derivative[0, 1][G3fun][
    scalar[], X[]])/scale[]^2 + (primescalar[]^2*timevec[h$2738945]*
   PD[-h$2738945][PD[-i][pertscalar[LI[1]]]]*Derivative[0, 1][G3fun][
    scalar[], X[]])/scale[]^2 + (pprimescalar[]*primescalar[]^2*
   PD[-i][pertB[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*primescalar[]^3*PD[-i][pertB[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (primescalar[]^3*PD[-i][pertpsi[LI[1]]]*Derivative[0, 1][G3fun][scalar[], 
    X[]])/scale[]^2 - (3*hubbleC[]*primescalar[]^2*PD[-i][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 pertS[LI[1], -i]*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]] + 
 primescalar[]^2*PD[-i][pertB[LI[1]]]*Derivative[1, 0][G3fun][scalar[], 
   X[]] + 2*primescalar[]*PD[-i][pertscalar[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]]
