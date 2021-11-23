(* Created with the Wolfram Language : www.wolfram.com *)
-(Mpl^2*hubbleC[]^2*perth[LI[1], -i, -j]) + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertphi[LI[1]] + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertpsi[LI[1]] - 
 density[]*pertshearten[LI[1], -i, -j] - pertshearten[LI[1], -i, -j]*
  pressure[] - 2*Mpl^2*perth[LI[1], -i, -j]*primehubbleC[] + 
 4*Mpl^2*metric\[Delta][-i, -j]*pertphi[LI[1]]*primehubbleC[] + 
 4*Mpl^2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primehubbleC[] - 
 G2fun[scalar[], X[]]*perth[LI[1], -i, -j]*scale[]^2 + 
 2*G2fun[scalar[], X[]]*metric\[Delta][-i, -j]*pertphi[LI[1]]*scale[]^2 - 
 metric\[Delta][-i, -j]*pertpressure[LI[1]]*scale[]^2 - 
 perth[LI[1], -i, -j]*pressure[]*scale[]^2 + 2*metric\[Delta][-i, -j]*
  pertphi[LI[1]]*pressure[]*scale[]^2 + Mpl^2*hubbleC[]*timevec[h$2746989]*
  PD[-h$2746989][perth[LI[1], -i, -j]] - 5*Mpl^2*hubbleC[]*
  metric\[Delta][-i, -j]*timevec[h$2746989]*PD[-h$2746989][pertphi[LI[1]]] - 
 Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$2746989]*
  PD[-h$2746989][pertpsi[LI[1]]] + Mpl^2*hubbleC[]*timevec[h$2746989]*
  PD[-h$2746989][PD[-i][pertF[LI[1], -j]]] + 
 (Mpl^2*timevec[h$2746989]*PD[-h$2746989][PD[-i][pertS[LI[1], -j]]])/2 + 
 Mpl^2*hubbleC[]*timevec[h$2746989]*PD[-h$2746989][
   PD[-j][pertF[LI[1], -i]]] + 
 (Mpl^2*timevec[h$2746989]*PD[-h$2746989][PD[-j][pertS[LI[1], -i]]])/2 - 
 Mpl^2*timevec[h$2746989]*PD[-h$2746989][PD[-j][PD[-i][pertB[LI[1]]]]] + 
 2*Mpl^2*hubbleC[]*timevec[h$2746989]*PD[-h$2746989][
   PD[-j][PD[-i][pertE[LI[1]]]]] + 
 (Mpl^2*timevec[h$2746989]*timevec[h$2746990]*
   PD[-h$2746990][PD[-h$2746989][perth[LI[1], -i, -j]]])/2 - 
 Mpl^2*metric\[Delta][-i, -j]*timevec[h$2746989]*timevec[h$2746990]*
  PD[-h$2746990][PD[-h$2746989][pertphi[LI[1]]]] + 
 (Mpl^2*timevec[h$2746989]*timevec[h$2746990]*
   PD[-h$2746990][PD[-h$2746989][PD[-i][pertF[LI[1], -j]]]])/2 + 
 (Mpl^2*timevec[h$2746989]*timevec[h$2746990]*
   PD[-h$2746990][PD[-h$2746989][PD[-j][pertF[LI[1], -i]]]])/2 + 
 Mpl^2*timevec[h$2746989]*timevec[h$2746990]*
  PD[-h$2746990][PD[-h$2746989][PD[-j][PD[-i][pertE[LI[1]]]]]] + 
 9*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$2747307]*
  PD[-h$2747307][pertphi[LI[1]]] + 3*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*
  timevec[h$2747307]*PD[-h$2747307][pertpsi[LI[1]]] + 
 3*Mpl^2*metric\[Delta][-i, -j]*timevec[h$2747307]*timevec[h$2747308]*
  PD[-h$2747308][PD[-h$2747307][pertphi[LI[1]]]] - 
 Mpl^2*hubbleC[]^2*PD[-i][pertF[LI[1], -j]] - 
 2*Mpl^2*primehubbleC[]*PD[-i][pertF[LI[1], -j]] - 
 G2fun[scalar[], X[]]*scale[]^2*PD[-i][pertF[LI[1], -j]] - 
 pressure[]*scale[]^2*PD[-i][pertF[LI[1], -j]] + 
 Mpl^2*hubbleC[]*PD[-i][pertS[LI[1], -j]] - 
 density[]*PD[-i][pertshearvec[LI[1], -j]] - 
 pressure[]*PD[-i][pertshearvec[LI[1], -j]] - 
 Mpl^2*hubbleC[]^2*PD[-j][pertF[LI[1], -i]] - 
 2*Mpl^2*primehubbleC[]*PD[-j][pertF[LI[1], -i]] - 
 G2fun[scalar[], X[]]*scale[]^2*PD[-j][pertF[LI[1], -i]] - 
 pressure[]*scale[]^2*PD[-j][pertF[LI[1], -i]] + 
 Mpl^2*hubbleC[]*PD[-j][pertS[LI[1], -i]] - 
 density[]*PD[-j][pertshearvec[LI[1], -i]] - 
 pressure[]*PD[-j][pertshearvec[LI[1], -i]] - 
 2*Mpl^2*hubbleC[]*PD[-j][PD[-i][pertB[LI[1]]]] - 
 2*Mpl^2*hubbleC[]^2*PD[-j][PD[-i][pertE[LI[1]]]] - 
 4*Mpl^2*primehubbleC[]*PD[-j][PD[-i][pertE[LI[1]]]] - 
 2*G2fun[scalar[], X[]]*scale[]^2*PD[-j][PD[-i][pertE[LI[1]]]] - 
 2*pressure[]*scale[]^2*PD[-j][PD[-i][pertE[LI[1]]]] + 
 Mpl^2*PD[-j][PD[-i][pertphi[LI[1]]]] - 
 Mpl^2*PD[-j][PD[-i][pertpsi[LI[1]]]] - 
 (density[]*PD[-j][PD[-i][pertshear[LI[1]]]])/2 - 
 (pressure[]*PD[-j][PD[-i][pertshear[LI[1]]]])/2 - 
 Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$2746991, p$2746992]*
  PD[-p$2746992][PD[-p$2746991][pertB[LI[1]]]] - 
 (Mpl^2*metric\[Delta][p$2746991, p$2746992]*
   PD[-p$2746992][PD[-p$2746991][perth[LI[1], -i, -j]]])/2 + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$2746991, p$2746992]*
  PD[-p$2746992][PD[-p$2746991][pertphi[LI[1]]]] + 
 Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$2746991, p$2746992]*
  timevec[h$2746989]*PD[-p$2746992][PD[-p$2746991][
    PD[-h$2746989][pertE[LI[1]]]]] + 3*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*
  metric\[Delta][p$2747309, p$2747310]*PD[-p$2747310][
   PD[-p$2747309][pertB[LI[1]]]] - 2*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$2747309, p$2747310]*PD[-p$2747310][
   PD[-p$2747309][pertphi[LI[1]]]] + Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$2747309, p$2747310]*PD[-p$2747310][
   PD[-p$2747309][pertpsi[LI[1]]]] + Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$2747309, p$2747310]*timevec[h$2747307]*
  PD[-p$2747310][PD[-p$2747309][PD[-h$2747307][pertB[LI[1]]]]] - 
 3*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$2747309, 
   p$2747310]*timevec[h$2747307]*PD[-p$2747310][
   PD[-p$2747309][PD[-h$2747307][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$2747309, p$2747310]*
  timevec[h$2747307]*timevec[h$2747308]*PD[-p$2747310][
   PD[-p$2747309][PD[-h$2747308][PD[-h$2747307][pertE[LI[1]]]]]] + 
 (density[]*metric\[Delta][-i, -j]*metric\[Delta][p$2747502, p$2747503]*
   PD[-p$2747503][PD[-p$2747502][pertshear[LI[1]]]])/2 + 
 (metric\[Delta][-i, -j]*metric\[Delta][p$2747502, p$2747503]*pressure[]*
   PD[-p$2747503][PD[-p$2747502][pertshear[LI[1]]]])/2 + 
 metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^2*
  Derivative[0, 1][G2fun][scalar[], X[]] - metric\[Delta][-i, -j]*
  primescalar[]*timevec[h$2747668]*PD[-h$2747668][pertscalar[LI[1]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] + 
 (perth[LI[1], -i, -j]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*pertphi[LI[1]]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*perth[LI[1], -i, -j]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*primescalar[]^3*timevec[h$2747812]*
   PD[-h$2747812][pertpsi[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]*
   timevec[h$2747812]*PD[-h$2747812][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (3*hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$2747812]*
   PD[-h$2747812][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$2747812]*
   timevec[h$2747813]*PD[-h$2747813][PD[-h$2747812][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (pprimescalar[]*primescalar[]^2*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*primescalar[]^3*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (pprimescalar[]*primescalar[]^2*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*primescalar[]^3*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*pprimescalar[]*primescalar[]^2*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*hubbleC[]*primescalar[]^3*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^4*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^5*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]^3*timevec[h$2747668]*
   PD[-h$2747668][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^4*
   timevec[h$2747668]*PD[-h$2747668][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 metric\[Delta][-i, -j]*pertscalar[LI[1]]*scale[]^2*
  Derivative[1, 0][G2fun][scalar[], X[]] + perth[LI[1], -i, -j]*
  primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]] - 
 2*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^2*
  Derivative[1, 0][G3fun][scalar[], X[]] - 2*metric\[Delta][-i, -j]*
  pertpsi[LI[1]]*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]] + 
 2*metric\[Delta][-i, -j]*primescalar[]*timevec[h$2747866]*
  PD[-h$2747866][pertscalar[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 primescalar[]^2*PD[-i][pertF[LI[1], -j]]*Derivative[1, 0][G3fun][scalar[], 
   X[]] + primescalar[]^2*PD[-j][pertF[LI[1], -i]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 
 2*primescalar[]^2*PD[-j][PD[-i][pertE[LI[1]]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 
 (metric\[Delta][-i, -j]*pertscalar[LI[1]]*pprimescalar[]*primescalar[]^2*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^3*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (metric\[Delta][-i, -j]*primescalar[]^3*timevec[h$2747668]*
   PD[-h$2747668][pertscalar[LI[1]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^2*
  Derivative[2, 0][G3fun][scalar[], X[]]
