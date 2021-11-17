(* Created with the Wolfram Language : www.wolfram.com *)
-(Mpl^2*hubbleC[]^2*perth[LI[2], -i, -j]) + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertphi[LI[2]] + 
 4*Mpl^2*hubbleC[]^2*perth[LI[1], -i, -j]*pertpsi[LI[1]] - 
 8*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]] - 
 8*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2 + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*pertpsi[LI[2]] + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*metric\[Delta][p$347314, 
   p$347315]*pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315] - 
 density[]*pertshearten[LI[2], -i, -j] - 
 2*density[]*metric\[Delta][-i, -p$347403]*metric\[Delta][-j, -p$347404]*
  pertvelocityvec[LI[1], p$347403]*pertvelocityvec[LI[1], p$347404] - 
 pertshearten[LI[2], -i, -j]*pressure[] - 2*metric\[Delta][-i, -p$347403]*
  metric\[Delta][-j, -p$347404]*pertvelocityvec[LI[1], p$347403]*
  pertvelocityvec[LI[1], p$347404]*pressure[] - 
 2*Mpl^2*perth[LI[2], -i, -j]*primehubbleC[] + 
 4*Mpl^2*metric\[Delta][-i, -j]*pertphi[LI[2]]*primehubbleC[] + 
 8*Mpl^2*perth[LI[1], -i, -j]*pertpsi[LI[1]]*primehubbleC[] - 
 16*Mpl^2*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*
  primehubbleC[] - 16*Mpl^2*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*
  primehubbleC[] + 4*Mpl^2*metric\[Delta][-i, -j]*pertpsi[LI[2]]*
  primehubbleC[] + 4*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$347314, p$347315]*pertS[LI[1], -p$347314]*
  pertS[LI[1], -p$347315]*primehubbleC[] + 
 2*density[]*metric\[Delta][-j, -p$347481]*pertS[LI[1], -i]*
  pertvelocityvec[LI[1], p$347481]*scale[] + 
 2*density[]*metric\[Delta][-i, -p$347481]*pertS[LI[1], -j]*
  pertvelocityvec[LI[1], p$347481]*scale[] + 2*metric\[Delta][-j, -p$347481]*
  pertS[LI[1], -i]*pertvelocityvec[LI[1], p$347481]*pressure[]*scale[] + 
 2*metric\[Delta][-i, -p$347481]*pertS[LI[1], -j]*
  pertvelocityvec[LI[1], p$347481]*pressure[]*scale[] - 
 G2fun[scalar[], X[]]*perth[LI[2], -i, -j]*scale[]^2 + 
 2*G2fun[scalar[], X[]]*metric\[Delta][-i, -j]*pertphi[LI[2]]*scale[]^2 - 
 2*perth[LI[1], -i, -j]*pertpressure[LI[1]]*scale[]^2 + 
 4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpressure[LI[1]]*scale[]^2 - 
 metric\[Delta][-i, -j]*pertpressure[LI[2]]*scale[]^2 - 
 2*density[]*pertS[LI[1], -i]*pertS[LI[1], -j]*scale[]^2 - 
 perth[LI[2], -i, -j]*pressure[]*scale[]^2 + 2*metric\[Delta][-i, -j]*
  pertphi[LI[2]]*pressure[]*scale[]^2 - 2*pertS[LI[1], -i]*pertS[LI[1], -j]*
  pressure[]*scale[]^2 + Mpl^2*hubbleC[]*timevec[h$344999]*
  PD[-h$344999][perth[LI[2], -i, -j]] + 
 (Mpl^2*timevec[h$344999]*timevec[h$345000]*
   PD[-h$345000][PD[-h$344999][perth[LI[2], -i, -j]]])/2 + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*timevec[h$345217]*
  PD[-h$345217][pertphi[LI[2]]] + 2*Mpl^2*metric\[Delta][-i, -j]*
  timevec[h$345217]*timevec[h$345218]*PD[-h$345218][
   PD[-h$345217][pertphi[LI[2]]]] - 4*Mpl^2*hubbleC[]*pertpsi[LI[1]]*
  timevec[h$345683]*PD[-h$345683][perth[LI[1], -i, -j]] + 
 4*Mpl^2*hubbleC[]*perth[LI[1], -i, -j]*timevec[h$345683]*
  PD[-h$345683][pertpsi[LI[1]]] - Mpl^2*timevec[h$345683]*timevec[h$345684]*
  PD[-h$345683][pertpsi[LI[1]]]*PD[-h$345684][perth[LI[1], -i, -j]] - 
 2*Mpl^2*pertpsi[LI[1]]*timevec[h$345683]*timevec[h$345684]*
  PD[-h$345684][PD[-h$345683][perth[LI[1], -i, -j]]] - 
 16*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*timevec[h$346074]*
  PD[-h$346074][pertphi[LI[1]]] - 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*
  pertphi[LI[1]]*timevec[h$346074]*PD[-h$346074][pertpsi[LI[1]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*timevec[h$346074]*timevec[h$346075]*
  PD[-h$346074][pertphi[LI[1]]]*PD[-h$346075][pertpsi[LI[1]]] - 
 8*Mpl^2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*timevec[h$346074]*
  timevec[h$346075]*PD[-h$346075][PD[-h$346074][pertphi[LI[1]]]] - 
 16*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*timevec[h$346451]*
  PD[-h$346451][pertpsi[LI[1]]] + 2*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*
  timevec[h$346694]*PD[-h$346694][pertpsi[LI[2]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
  pertS[LI[1], -p$347314]*timevec[h$347316]*
  PD[-h$347316][pertS[LI[1], -p$347315]] + 
 2*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
  metric\[Delta][p$348241, p$348242]*perth[LI[1], -p$348239, -p$348241]*
  timevec[h$348243]*PD[-h$348243][perth[LI[1], -p$348240, -p$348242]] - 
 Mpl^2*metric\[Delta][p$348239, p$348240]*timevec[h$348243]*timevec[h$348244]*
  PD[-h$348243][perth[LI[1], -i, -p$348239]]*
  PD[-h$348244][perth[LI[1], -j, -p$348240]] + 
 (3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
   metric\[Delta][p$348241, p$348242]*timevec[h$348243]*timevec[h$348244]*
   PD[-h$348243][perth[LI[1], -p$348239, -p$348241]]*
   PD[-h$348244][perth[LI[1], -p$348240, -p$348242]])/4 + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
  metric\[Delta][p$348241, p$348242]*perth[LI[1], -p$348239, -p$348241]*
  timevec[h$348243]*timevec[h$348244]*PD[-h$348244][
   PD[-h$348243][perth[LI[1], -p$348240, -p$348242]]] + 
 12*Mpl^2*hubbleC[]*perth[LI[1], -i, -j]*timevec[h$348533]*
  PD[-h$348533][pertphi[LI[1]]] + Mpl^2*timevec[h$348533]*timevec[h$348534]*
  PD[-h$348533][pertphi[LI[1]]]*PD[-h$348534][perth[LI[1], -i, -j]] + 
 6*Mpl^2*perth[LI[1], -i, -j]*timevec[h$348533]*timevec[h$348534]*
  PD[-h$348534][PD[-h$348533][pertphi[LI[1]]]] - 
 4*Mpl^2*hubbleC[]*pertpsi[LI[1]]*timevec[h$349419]*
  PD[-h$349419][PD[-i][pertF[LI[1], -j]]] - 4*Mpl^2*hubbleC[]*pertpsi[LI[1]]*
  timevec[h$349419]*PD[-h$349419][PD[-j][pertF[LI[1], -i]]] - 
 2*Mpl^2*pertpsi[LI[1]]*timevec[h$349419]*timevec[h$349420]*
  PD[-h$349420][PD[-h$349419][PD[-i][pertF[LI[1], -j]]]] - 
 2*Mpl^2*pertpsi[LI[1]]*timevec[h$349419]*timevec[h$349420]*
  PD[-h$349420][PD[-h$349419][PD[-j][pertF[LI[1], -i]]]] - 
 Mpl^2*timevec[h$349419]*timevec[h$349420]*PD[-h$349419][pertpsi[LI[1]]]*
  PD[-h$349420][PD[-i][pertF[LI[1], -j]]] - Mpl^2*timevec[h$349419]*
  timevec[h$349420]*PD[-h$349419][pertpsi[LI[1]]]*
  PD[-h$349420][PD[-j][pertF[LI[1], -i]]] + Mpl^2*hubbleC[]*timevec[h$349763]*
  PD[-h$349763][PD[-i][pertF[LI[2], -j]]] + Mpl^2*hubbleC[]*timevec[h$349763]*
  PD[-h$349763][PD[-j][pertF[LI[2], -i]]] + 
 (Mpl^2*timevec[h$349763]*timevec[h$349764]*
   PD[-h$349764][PD[-h$349763][PD[-i][pertF[LI[2], -j]]]])/2 + 
 (Mpl^2*timevec[h$349763]*timevec[h$349764]*
   PD[-h$349764][PD[-h$349763][PD[-j][pertF[LI[2], -i]]]])/2 - 
 Mpl^2*metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  timevec[h$350243]*PD[-h$350243][PD[-i][perth[LI[1], -j, -p$350242]]] - 
 Mpl^2*metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  timevec[h$350243]*PD[-h$350243][PD[-j][perth[LI[1], -i, -p$350242]]] + 
 2*Mpl^2*pertS[LI[1], -j]*timevec[h$350467]*
  PD[-h$350467][PD[-i][pertphi[LI[1]]]] + 2*Mpl^2*pertS[LI[1], -i]*
  timevec[h$350467]*PD[-h$350467][PD[-j][pertphi[LI[1]]]] - 
 2*Mpl^2*pertpsi[LI[1]]*timevec[h$350660]*
  PD[-h$350660][PD[-i][pertS[LI[1], -j]]] - 2*Mpl^2*pertpsi[LI[1]]*
  timevec[h$350660]*PD[-h$350660][PD[-j][pertS[LI[1], -i]]] + 
 (Mpl^2*timevec[h$350737]*PD[-h$350737][PD[-i][pertS[LI[2], -j]]])/2 + 
 (Mpl^2*timevec[h$350737]*PD[-h$350737][PD[-j][pertS[LI[2], -i]]])/2 + 
 4*Mpl^2*pertpsi[LI[1]]*timevec[h$350989]*
  PD[-h$350989][PD[-j][PD[-i][pertB[LI[1]]]]] - 
 Mpl^2*timevec[h$351102]*PD[-h$351102][PD[-j][PD[-i][pertB[LI[2]]]]] - 
 8*Mpl^2*hubbleC[]*pertpsi[LI[1]]*timevec[h$351798]*
  PD[-h$351798][PD[-j][PD[-i][pertE[LI[1]]]]] - 
 4*Mpl^2*pertpsi[LI[1]]*timevec[h$351798]*timevec[h$351799]*
  PD[-h$351799][PD[-h$351798][PD[-j][PD[-i][pertE[LI[1]]]]]] - 
 2*Mpl^2*timevec[h$351798]*timevec[h$351799]*PD[-h$351798][pertpsi[LI[1]]]*
  PD[-h$351799][PD[-j][PD[-i][pertE[LI[1]]]]] + 
 2*Mpl^2*hubbleC[]*timevec[h$352089]*PD[-h$352089][
   PD[-j][PD[-i][pertE[LI[2]]]]] + Mpl^2*timevec[h$352089]*timevec[h$352090]*
  PD[-h$352090][PD[-h$352089][PD[-j][PD[-i][pertE[LI[2]]]]]] - 
 2*Mpl^2*metric\[Delta][p$352837, p$352838]*pertS[LI[1], -p$352837]*
  timevec[h$352839]*PD[-h$352839][PD[-j][PD[-i][pertF[LI[1], -p$352838]]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*timevec[h$354171]*timevec[h$354172]*
  PD[-h$354171][PD[-j][pertF[LI[1], -p$354169]]]*
  PD[-h$354172][perth[LI[1], -i, -p$354170]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*timevec[h$354171]*timevec[h$354172]*
  PD[-h$354171][PD[-i][pertF[LI[1], -p$354169]]]*
  PD[-h$354172][perth[LI[1], -j, -p$354170]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*timevec[h$354338]*timevec[h$354339]*
  PD[-h$354338][pertphi[LI[1]]]*PD[-h$354339][pertphi[LI[1]]] + 
 Mpl^2*timevec[h$354932]*timevec[h$354933]*PD[-h$354932][pertphi[LI[1]]]*
  PD[-h$354933][PD[-i][pertF[LI[1], -j]]] + Mpl^2*timevec[h$354932]*
  timevec[h$354933]*PD[-h$354932][pertphi[LI[1]]]*
  PD[-h$354933][PD[-j][pertF[LI[1], -i]]] - 
 Mpl^2*metric\[Delta][p$355733, p$355734]*timevec[h$355735]*timevec[h$355736]*
  PD[-h$355735][PD[-i][pertF[LI[1], -p$355733]]]*
  PD[-h$355736][PD[-j][pertF[LI[1], -p$355734]]] + 
 2*Mpl^2*timevec[h$356253]*timevec[h$356254]*PD[-h$356253][pertphi[LI[1]]]*
  PD[-h$356254][PD[-j][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*timevec[h$356527]*PD[-h$356527][PD[-j][pertphi[LI[1]]]]*
  PD[-i][pertB[LI[1]]] + 4*Mpl^2*hubbleC[]^2*pertpsi[LI[1]]*
  PD[-i][pertF[LI[1], -j]] + 8*Mpl^2*pertpsi[LI[1]]*primehubbleC[]*
  PD[-i][pertF[LI[1], -j]] - 2*pertpressure[LI[1]]*scale[]^2*
  PD[-i][pertF[LI[1], -j]] + 4*Mpl^2*hubbleC[]*timevec[h$349419]*
  PD[-h$349419][pertpsi[LI[1]]]*PD[-i][pertF[LI[1], -j]] + 
 12*Mpl^2*hubbleC[]*timevec[h$354932]*PD[-h$354932][pertphi[LI[1]]]*
  PD[-i][pertF[LI[1], -j]] + 6*Mpl^2*timevec[h$354932]*timevec[h$354933]*
  PD[-h$354933][PD[-h$354932][pertphi[LI[1]]]]*PD[-i][pertF[LI[1], -j]] - 
 Mpl^2*hubbleC[]^2*PD[-i][pertF[LI[2], -j]] - 
 2*Mpl^2*primehubbleC[]*PD[-i][pertF[LI[2], -j]] - 
 G2fun[scalar[], X[]]*scale[]^2*PD[-i][pertF[LI[2], -j]] - 
 pressure[]*scale[]^2*PD[-i][pertF[LI[2], -j]] - 
 2*Mpl^2*hubbleC[]*metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  PD[-i][perth[LI[1], -j, -p$350242]] - 
 Mpl^2*metric\[Delta][p$350241, p$350242]*timevec[h$350243]*
  PD[-h$350243][pertS[LI[1], -p$350241]]*
  PD[-i][perth[LI[1], -j, -p$350242]] + 4*Mpl^2*hubbleC[]*pertS[LI[1], -j]*
  PD[-i][pertphi[LI[1]]] + 2*Mpl^2*timevec[h$350467]*
  PD[-h$350467][pertS[LI[1], -j]]*PD[-i][pertphi[LI[1]]] - 
 2*Mpl^2*timevec[h$356527]*PD[-h$356527][PD[-j][pertB[LI[1]]]]*
  PD[-i][pertphi[LI[1]]] - 4*Mpl^2*hubbleC[]*pertpsi[LI[1]]*
  PD[-i][pertS[LI[1], -j]] - Mpl^2*timevec[h$350467]*
  PD[-h$350467][pertphi[LI[1]]]*PD[-i][pertS[LI[1], -j]] - 
 Mpl^2*timevec[h$350660]*PD[-h$350660][pertpsi[LI[1]]]*
  PD[-i][pertS[LI[1], -j]] + Mpl^2*hubbleC[]*PD[-i][pertS[LI[2], -j]] - 
 density[]*PD[-i][pertshearvec[LI[2], -j]] - 
 pressure[]*PD[-i][pertshearvec[LI[2], -j]] - 
 2*density[]*metric\[Delta][-j, -p$356728]*pertvelocityvec[LI[1], p$356728]*
  PD[-i][pertvelocity[LI[1]]] - 2*metric\[Delta][-j, -p$356728]*
  pertvelocityvec[LI[1], p$356728]*pressure[]*PD[-i][pertvelocity[LI[1]]] + 
 2*density[]*pertS[LI[1], -j]*scale[]*PD[-i][pertvelocity[LI[1]]] + 
 2*pertS[LI[1], -j]*pressure[]*scale[]*PD[-i][pertvelocity[LI[1]]] - 
 2*Mpl^2*timevec[h$356527]*PD[-h$356527][PD[-i][pertphi[LI[1]]]]*
  PD[-j][pertB[LI[1]]] - 4*Mpl^2*hubbleC[]*PD[-i][pertphi[LI[1]]]*
  PD[-j][pertB[LI[1]]] + 4*Mpl^2*hubbleC[]^2*pertpsi[LI[1]]*
  PD[-j][pertF[LI[1], -i]] + 8*Mpl^2*pertpsi[LI[1]]*primehubbleC[]*
  PD[-j][pertF[LI[1], -i]] - 2*pertpressure[LI[1]]*scale[]^2*
  PD[-j][pertF[LI[1], -i]] + 4*Mpl^2*hubbleC[]*timevec[h$349419]*
  PD[-h$349419][pertpsi[LI[1]]]*PD[-j][pertF[LI[1], -i]] + 
 12*Mpl^2*hubbleC[]*timevec[h$354932]*PD[-h$354932][pertphi[LI[1]]]*
  PD[-j][pertF[LI[1], -i]] + 6*Mpl^2*timevec[h$354932]*timevec[h$354933]*
  PD[-h$354933][PD[-h$354932][pertphi[LI[1]]]]*PD[-j][pertF[LI[1], -i]] - 
 Mpl^2*hubbleC[]^2*PD[-j][pertF[LI[2], -i]] - 
 2*Mpl^2*primehubbleC[]*PD[-j][pertF[LI[2], -i]] - 
 G2fun[scalar[], X[]]*scale[]^2*PD[-j][pertF[LI[2], -i]] - 
 pressure[]*scale[]^2*PD[-j][pertF[LI[2], -i]] - 
 2*Mpl^2*hubbleC[]*metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  PD[-j][perth[LI[1], -i, -p$350242]] - 
 Mpl^2*metric\[Delta][p$350241, p$350242]*timevec[h$350243]*
  PD[-h$350243][pertS[LI[1], -p$350241]]*
  PD[-j][perth[LI[1], -i, -p$350242]] + 
 (Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
   PD[-i][perth[LI[1], -p$348239, -p$348241]]*
   PD[-j][perth[LI[1], -p$348240, -p$348242]])/2 + 
 4*Mpl^2*hubbleC[]*pertS[LI[1], -i]*PD[-j][pertphi[LI[1]]] + 
 2*Mpl^2*timevec[h$350467]*PD[-h$350467][pertS[LI[1], -i]]*
  PD[-j][pertphi[LI[1]]] - 2*Mpl^2*timevec[h$356527]*
  PD[-h$356527][PD[-i][pertB[LI[1]]]]*PD[-j][pertphi[LI[1]]] - 
 4*Mpl^2*hubbleC[]*PD[-i][pertB[LI[1]]]*PD[-j][pertphi[LI[1]]] + 
 6*Mpl^2*PD[-i][pertphi[LI[1]]]*PD[-j][pertphi[LI[1]]] - 
 2*Mpl^2*PD[-i][pertpsi[LI[1]]]*PD[-j][pertphi[LI[1]]] - 
 2*Mpl^2*PD[-i][pertphi[LI[1]]]*PD[-j][pertpsi[LI[1]]] + 
 2*Mpl^2*PD[-i][pertpsi[LI[1]]]*PD[-j][pertpsi[LI[1]]] - 
 4*Mpl^2*hubbleC[]*pertpsi[LI[1]]*PD[-j][pertS[LI[1], -i]] - 
 Mpl^2*timevec[h$350467]*PD[-h$350467][pertphi[LI[1]]]*
  PD[-j][pertS[LI[1], -i]] - Mpl^2*timevec[h$350660]*
  PD[-h$350660][pertpsi[LI[1]]]*PD[-j][pertS[LI[1], -i]] - 
 Mpl^2*metric\[Delta][p$347314, p$347315]*PD[-i][pertS[LI[1], -p$347314]]*
  PD[-j][pertS[LI[1], -p$347315]] + Mpl^2*hubbleC[]*
  PD[-j][pertS[LI[2], -i]] - density[]*PD[-j][pertshearvec[LI[2], -i]] - 
 pressure[]*PD[-j][pertshearvec[LI[2], -i]] - 
 2*density[]*metric\[Delta][-i, -p$356728]*pertvelocityvec[LI[1], p$356728]*
  PD[-j][pertvelocity[LI[1]]] - 2*metric\[Delta][-i, -p$356728]*
  pertvelocityvec[LI[1], p$356728]*pressure[]*PD[-j][pertvelocity[LI[1]]] + 
 2*density[]*pertS[LI[1], -i]*scale[]*PD[-j][pertvelocity[LI[1]]] + 
 2*pertS[LI[1], -i]*pressure[]*scale[]*PD[-j][pertvelocity[LI[1]]] - 
 2*density[]*PD[-i][pertvelocity[LI[1]]]*PD[-j][pertvelocity[LI[1]]] - 
 2*pressure[]*PD[-i][pertvelocity[LI[1]]]*PD[-j][pertvelocity[LI[1]]] + 
 8*Mpl^2*hubbleC[]*pertpsi[LI[1]]*PD[-j][PD[-i][pertB[LI[1]]]] + 
 2*Mpl^2*timevec[h$350989]*PD[-h$350989][pertpsi[LI[1]]]*
  PD[-j][PD[-i][pertB[LI[1]]]] + 2*Mpl^2*timevec[h$356527]*
  PD[-h$356527][pertphi[LI[1]]]*PD[-j][PD[-i][pertB[LI[1]]]] - 
 2*Mpl^2*hubbleC[]*PD[-j][PD[-i][pertB[LI[2]]]] + 
 8*Mpl^2*hubbleC[]^2*pertpsi[LI[1]]*PD[-j][PD[-i][pertE[LI[1]]]] + 
 16*Mpl^2*pertpsi[LI[1]]*primehubbleC[]*PD[-j][PD[-i][pertE[LI[1]]]] - 
 4*pertpressure[LI[1]]*scale[]^2*PD[-j][PD[-i][pertE[LI[1]]]] + 
 8*Mpl^2*hubbleC[]*timevec[h$351798]*PD[-h$351798][pertpsi[LI[1]]]*
  PD[-j][PD[-i][pertE[LI[1]]]] + 24*Mpl^2*hubbleC[]*timevec[h$356253]*
  PD[-h$356253][pertphi[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]] + 
 12*Mpl^2*timevec[h$356253]*timevec[h$356254]*
  PD[-h$356254][PD[-h$356253][pertphi[LI[1]]]]*PD[-j][PD[-i][pertE[LI[1]]]] - 
 2*Mpl^2*hubbleC[]^2*PD[-j][PD[-i][pertE[LI[2]]]] - 
 4*Mpl^2*primehubbleC[]*PD[-j][PD[-i][pertE[LI[2]]]] - 
 2*G2fun[scalar[], X[]]*scale[]^2*PD[-j][PD[-i][pertE[LI[2]]]] - 
 2*pressure[]*scale[]^2*PD[-j][PD[-i][pertE[LI[2]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$352837, p$352838]*pertS[LI[1], -p$352837]*
  PD[-j][PD[-i][pertF[LI[1], -p$352838]]] - 
 2*Mpl^2*metric\[Delta][p$352837, p$352838]*timevec[h$352839]*
  PD[-h$352839][pertS[LI[1], -p$352837]]*
  PD[-j][PD[-i][pertF[LI[1], -p$352838]]] + 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  perth[LI[1], -p$348239, -p$348241]*
  PD[-j][PD[-i][perth[LI[1], -p$348240, -p$348242]]] + 
 4*Mpl^2*pertphi[LI[1]]*PD[-j][PD[-i][pertphi[LI[1]]]] + 
 Mpl^2*PD[-j][PD[-i][pertphi[LI[2]]]] + 4*Mpl^2*pertpsi[LI[1]]*
  PD[-j][PD[-i][pertpsi[LI[1]]]] - Mpl^2*PD[-j][PD[-i][pertpsi[LI[2]]]] - 
 2*Mpl^2*metric\[Delta][p$347314, p$347315]*pertS[LI[1], -p$347314]*
  PD[-j][PD[-i][pertS[LI[1], -p$347315]]] - 
 (density[]*PD[-j][PD[-i][pertshear[LI[2]]]])/2 - 
 (pressure[]*PD[-j][PD[-i][pertshear[LI[2]]]])/2 - 
 (Mpl^2*metric\[Delta][p$345001, p$345002]*
   PD[-p$345002][PD[-p$345001][perth[LI[2], -i, -j]]])/2 - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$345219, p$345220]*
  PD[-p$345220][PD[-p$345219][pertphi[LI[2]]]] + 
 Mpl^2*metric\[Delta][p$345685, p$345686]*PD[-i][perth[LI[1], -j, -p$345686]]*
  PD[-p$345685][pertpsi[LI[1]]] + Mpl^2*metric\[Delta][p$345685, p$345686]*
  PD[-j][perth[LI[1], -i, -p$345686]]*PD[-p$345685][pertpsi[LI[1]]] - 
 Mpl^2*metric\[Delta][p$345685, p$345686]*PD[-p$345685][pertpsi[LI[1]]]*
  PD[-p$345686][perth[LI[1], -i, -j]] + 
 2*Mpl^2*metric\[Delta][p$345685, p$345686]*perth[LI[1], -i, -j]*
  PD[-p$345686][PD[-p$345685][pertpsi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$345685, p$345686]*
  metric\[Delta][p$345687, p$345688]*perth[LI[1], -p$345685, -p$345687]*
  PD[-p$345688][PD[-p$345686][pertpsi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$346452, p$346453]*
  PD[-p$346452][pertpsi[LI[1]]]*PD[-p$346453][pertpsi[LI[1]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$346452, p$346453]*
  pertpsi[LI[1]]*PD[-p$346453][PD[-p$346452][pertpsi[LI[1]]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$346695, p$346696]*
  PD[-p$346696][PD[-p$346695][pertpsi[LI[2]]]] - 
 Mpl^2*metric\[Delta][p$347314, p$347315]*PD[-p$347314][pertS[LI[1], -i]]*
  PD[-p$347315][pertS[LI[1], -j]] + Mpl^2*metric\[Delta][p$347314, p$347315]*
  pertS[LI[1], -p$347314]*PD[-p$347315][PD[-i][pertS[LI[1], -j]]] + 
 Mpl^2*metric\[Delta][p$347314, p$347315]*pertS[LI[1], -p$347314]*
  PD[-p$347315][PD[-j][pertS[LI[1], -i]]] - 
 (Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347318]*
   metric\[Delta][p$347315, p$347317]*PD[-p$347315][pertS[LI[1], -p$347314]]*
   PD[-p$347318][pertS[LI[1], -p$347317]])/2 + 
 (3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347317]*
   metric\[Delta][p$347315, p$347318]*PD[-p$347315][pertS[LI[1], -p$347314]]*
   PD[-p$347318][pertS[LI[1], -p$347317]])/2 + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
  metric\[Delta][p$347317, p$347318]*pertS[LI[1], -p$347314]*
  PD[-p$347318][PD[-p$347317][pertS[LI[1], -p$347315]]] - 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  PD[-p$348240][perth[LI[1], -j, -p$348242]]*
  PD[-p$348241][perth[LI[1], -i, -p$348239]] + 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  PD[-p$348241][perth[LI[1], -i, -p$348239]]*
  PD[-p$348242][perth[LI[1], -j, -p$348240]] - 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  perth[LI[1], -p$348239, -p$348241]*PD[-p$348242][
   PD[-i][perth[LI[1], -j, -p$348240]]] - 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  perth[LI[1], -p$348239, -p$348241]*PD[-p$348242][
   PD[-j][perth[LI[1], -i, -p$348240]]] + 
 Mpl^2*metric\[Delta][p$348239, p$348240]*metric\[Delta][p$348241, p$348242]*
  perth[LI[1], -p$348239, -p$348241]*PD[-p$348242][
   PD[-p$348240][perth[LI[1], -i, -j]]] + 
 (Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
   metric\[Delta][p$348241, p$348242]*metric\[Delta][p$348245, p$348246]*
   PD[-p$348242][perth[LI[1], -p$348240, -p$348246]]*
   PD[-p$348245][perth[LI[1], -p$348239, -p$348241]])/2 - 
 (3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
   metric\[Delta][p$348241, p$348242]*metric\[Delta][p$348245, p$348246]*
   PD[-p$348245][perth[LI[1], -p$348239, -p$348241]]*
   PD[-p$348246][perth[LI[1], -p$348240, -p$348242]])/4 - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$348239, p$348240]*
  metric\[Delta][p$348241, p$348242]*metric\[Delta][p$348245, p$348246]*
  perth[LI[1], -p$348239, -p$348241]*PD[-p$348246][
   PD[-p$348245][perth[LI[1], -p$348240, -p$348242]]] + 
 Mpl^2*metric\[Delta][p$348535, p$348536]*PD[-i][perth[LI[1], -j, -p$348536]]*
  PD[-p$348535][pertphi[LI[1]]] + Mpl^2*metric\[Delta][p$348535, p$348536]*
  PD[-j][perth[LI[1], -i, -p$348536]]*PD[-p$348535][pertphi[LI[1]]] - 
 3*Mpl^2*metric\[Delta][p$348535, p$348536]*PD[-p$348535][pertphi[LI[1]]]*
  PD[-p$348536][perth[LI[1], -i, -j]] + 
 2*Mpl^2*metric\[Delta][p$348535, p$348536]*perth[LI[1], -j, -p$348535]*
  PD[-p$348536][PD[-i][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$348535, p$348536]*perth[LI[1], -i, -p$348535]*
  PD[-p$348536][PD[-j][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$348535, p$348536]*pertphi[LI[1]]*
  PD[-p$348536][PD[-p$348535][perth[LI[1], -i, -j]]] - 
 4*Mpl^2*metric\[Delta][p$348535, p$348536]*perth[LI[1], -i, -j]*
  PD[-p$348536][PD[-p$348535][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$349421, p$349422]*
  PD[-j][PD[-i][pertF[LI[1], -p$349422]]]*PD[-p$349421][pertpsi[LI[1]]] + 
 2*Mpl^2*metric\[Delta][p$349421, p$349422]*PD[-i][pertF[LI[1], -j]]*
  PD[-p$349422][PD[-p$349421][pertpsi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$349421, p$349422]*PD[-j][pertF[LI[1], -i]]*
  PD[-p$349422][PD[-p$349421][pertpsi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$349421, p$349422]*
  metric\[Delta][p$349423, p$349424]*PD[-p$349421][pertpsi[LI[1]]]*
  PD[-p$349424][PD[-p$349423][pertF[LI[1], -p$349422]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$349421, p$349423]*
  metric\[Delta][p$349422, p$349424]*PD[-p$349422][pertF[LI[1], -p$349421]]*
  PD[-p$349424][PD[-p$349423][pertpsi[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$350241, p$350242]*timevec[h$350243]*
  PD[-h$350243][perth[LI[1], -j, -p$350242]]*
  PD[-p$350241][pertS[LI[1], -i]] - Mpl^2*metric\[Delta][p$350241, p$350242]*
  timevec[h$350243]*PD[-h$350243][perth[LI[1], -i, -p$350242]]*
  PD[-p$350241][pertS[LI[1], -j]] + 2*Mpl^2*hubbleC[]*
  metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  PD[-p$350242][perth[LI[1], -i, -j]] + 
 Mpl^2*metric\[Delta][p$350241, p$350242]*timevec[h$350243]*
  PD[-h$350243][pertS[LI[1], -p$350241]]*
  PD[-p$350242][perth[LI[1], -i, -j]] + 4*Mpl^2*hubbleC[]*
  metric\[Delta][-i, -j]*metric\[Delta][p$350241, p$350244]*
  metric\[Delta][p$350242, p$350245]*perth[LI[1], -p$350244, -p$350245]*
  PD[-p$350242][pertS[LI[1], -p$350241]] + Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$350241, p$350244]*metric\[Delta][p$350242, p$350245]*
  timevec[h$350243]*PD[-h$350243][perth[LI[1], -p$350244, -p$350245]]*
  PD[-p$350242][pertS[LI[1], -p$350241]] + 
 2*Mpl^2*metric\[Delta][p$350241, p$350242]*pertS[LI[1], -p$350241]*
  timevec[h$350243]*PD[-p$350242][PD[-h$350243][perth[LI[1], -i, -j]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$350241, p$350242]*
  metric\[Delta][p$350244, p$350245]*perth[LI[1], -p$350241, -p$350244]*
  timevec[h$350243]*PD[-p$350245][PD[-h$350243][pertS[LI[1], -p$350242]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$350468, p$350469]*
  pertS[LI[1], -p$350468]*timevec[h$350467]*
  PD[-p$350469][PD[-h$350467][pertphi[LI[1]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$350661, p$350662]*
  pertS[LI[1], -p$350661]*PD[-p$350662][pertpsi[LI[1]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$350990, p$350991]*
  PD[-p$350990][pertB[LI[1]]]*PD[-p$350991][pertpsi[LI[1]]] - 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$350990, p$350991]*
  pertpsi[LI[1]]*PD[-p$350991][PD[-p$350990][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$350990, p$350991]*
  timevec[h$350989]*PD[-h$350989][pertpsi[LI[1]]]*
  PD[-p$350991][PD[-p$350990][pertB[LI[1]]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$350990, p$350991]*
  pertpsi[LI[1]]*timevec[h$350989]*PD[-p$350991][
   PD[-p$350990][PD[-h$350989][pertB[LI[1]]]]] + 
 2*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$351103, p$351104]*
  PD[-p$351104][PD[-p$351103][pertB[LI[2]]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$351103, p$351104]*
  timevec[h$351102]*PD[-p$351104][PD[-p$351103][
    PD[-h$351102][pertB[LI[2]]]]] + 
 2*Mpl^2*metric\[Delta][p$351800, p$351801]*PD[-p$351800][pertpsi[LI[1]]]*
  PD[-p$351801][PD[-j][PD[-i][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][p$351800, p$351801]*PD[-j][PD[-i][pertE[LI[1]]]]*
  PD[-p$351801][PD[-p$351800][pertpsi[LI[1]]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$351800, p$351801]*
  pertpsi[LI[1]]*timevec[h$351798]*PD[-p$351801][
   PD[-p$351800][PD[-h$351798][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$351800, p$351801]*
  timevec[h$351798]*timevec[h$351799]*PD[-h$351798][pertpsi[LI[1]]]*
  PD[-p$351801][PD[-p$351800][PD[-h$351799][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$351800, p$351801]*
  pertpsi[LI[1]]*timevec[h$351798]*timevec[h$351799]*
  PD[-p$351801][PD[-p$351800][PD[-h$351799][PD[-h$351798][pertE[LI[1]]]]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$351800, p$351801]*
  metric\[Delta][p$351802, p$351803]*PD[-p$351802][
   PD[-p$351800][pertE[LI[1]]]]*PD[-p$351803][
   PD[-p$351801][pertpsi[LI[1]]]] - 2*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$351800, p$351801]*metric\[Delta][p$351802, p$351803]*
  PD[-p$351800][pertpsi[LI[1]]]*PD[-p$351803][
   PD[-p$351802][PD[-p$351801][pertE[LI[1]]]]] - 
 2*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$352091, p$352092]*
  timevec[h$352089]*PD[-p$352092][PD[-p$352091][
    PD[-h$352089][pertE[LI[2]]]]] - Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$352091, p$352092]*timevec[h$352089]*timevec[h$352090]*
  PD[-p$352092][PD[-p$352091][PD[-h$352090][PD[-h$352089][pertE[LI[2]]]]]] - 
 Mpl^2*metric\[Delta][p$352837, p$352838]*timevec[h$352839]*
  PD[-h$352839][PD[-j][pertF[LI[1], -p$352838]]]*
  PD[-p$352837][pertS[LI[1], -i]] - Mpl^2*metric\[Delta][p$352837, p$352838]*
  timevec[h$352839]*PD[-h$352839][PD[-i][pertF[LI[1], -p$352838]]]*
  PD[-p$352837][pertS[LI[1], -j]] - Mpl^2*metric\[Delta][p$352837, p$352838]*
  timevec[h$352839]*PD[-p$352837][pertS[LI[1], -j]]*
  PD[-p$352838][PD[-h$352839][pertF[LI[1], -i]]] - 
 Mpl^2*metric\[Delta][p$352837, p$352838]*timevec[h$352839]*
  PD[-p$352837][pertS[LI[1], -i]]*PD[-p$352838][
   PD[-h$352839][pertF[LI[1], -j]]] + 
 Mpl^2*metric\[Delta][p$352837, p$352838]*pertS[LI[1], -p$352837]*
  timevec[h$352839]*PD[-p$352838][PD[-h$352839][PD[-i][pertF[LI[1], -j]]]] + 
 Mpl^2*metric\[Delta][p$352837, p$352838]*pertS[LI[1], -p$352837]*
  timevec[h$352839]*PD[-p$352838][PD[-h$352839][PD[-j][pertF[LI[1], -i]]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352840]*
  metric\[Delta][p$352838, p$352841]*timevec[h$352839]*
  PD[-p$352838][pertS[LI[1], -p$352837]]*
  PD[-p$352840][PD[-h$352839][pertF[LI[1], -p$352841]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352840]*
  metric\[Delta][p$352838, p$352841]*timevec[h$352839]*
  PD[-p$352838][pertF[LI[1], -p$352837]]*
  PD[-p$352840][PD[-h$352839][pertS[LI[1], -p$352841]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352841]*
  metric\[Delta][p$352838, p$352840]*PD[-p$352838][pertF[LI[1], -p$352837]]*
  PD[-p$352841][pertS[LI[1], -p$352840]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352840]*
  metric\[Delta][p$352838, p$352841]*PD[-p$352838][pertF[LI[1], -p$352837]]*
  PD[-p$352841][pertS[LI[1], -p$352840]] + Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$352837, p$352840]*metric\[Delta][p$352838, p$352841]*
  timevec[h$352839]*PD[-p$352838][pertS[LI[1], -p$352837]]*
  PD[-p$352841][PD[-h$352839][pertF[LI[1], -p$352840]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352840]*
  metric\[Delta][p$352838, p$352841]*timevec[h$352839]*
  PD[-p$352838][pertF[LI[1], -p$352837]]*
  PD[-p$352841][PD[-h$352839][pertS[LI[1], -p$352840]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352838]*
  metric\[Delta][p$352840, p$352841]*pertS[LI[1], -p$352837]*
  PD[-p$352841][PD[-p$352840][pertF[LI[1], -p$352838]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352838]*
  metric\[Delta][p$352840, p$352841]*timevec[h$352839]*
  PD[-h$352839][pertS[LI[1], -p$352837]]*
  PD[-p$352841][PD[-p$352840][pertF[LI[1], -p$352838]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$352837, p$352838]*
  metric\[Delta][p$352840, p$352841]*pertS[LI[1], -p$352837]*
  timevec[h$352839]*PD[-p$352841][PD[-p$352840][
    PD[-h$352839][pertF[LI[1], -p$352838]]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*timevec[h$354171]*timevec[h$354172]*
  PD[-h$354172][perth[LI[1], -j, -p$354170]]*
  PD[-p$354169][PD[-h$354171][pertF[LI[1], -i]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*timevec[h$354171]*timevec[h$354172]*
  PD[-h$354172][perth[LI[1], -i, -p$354170]]*
  PD[-p$354169][PD[-h$354171][pertF[LI[1], -j]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354173]*
  metric\[Delta][p$354170, p$354174]*timevec[h$354171]*
  PD[-h$354171][perth[LI[1], -p$354173, -p$354174]]*
  PD[-p$354170][pertF[LI[1], -p$354169]] + 2*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, p$354174]*
  timevec[h$354171]*timevec[h$354172]*PD[-h$354172][
   PD[-h$354171][perth[LI[1], -p$354173, -p$354174]]]*
  PD[-p$354170][pertF[LI[1], -p$354169]] + 
 2*Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, 
   p$354174]*PD[-j][PD[-i][perth[LI[1], -p$354173, -p$354174]]]*
  PD[-p$354170][pertF[LI[1], -p$354169]] + 3*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  timevec[h$354171]*timevec[h$354172]*PD[-h$354172][
   perth[LI[1], -p$354170, -p$354174]]*
  PD[-p$354173][PD[-h$354171][pertF[LI[1], -p$354169]]] + 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-j][perth[LI[1], -p$354170, -p$354174]]*
  PD[-p$354173][PD[-i][pertF[LI[1], -p$354169]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-p$354170][perth[LI[1], -j, -p$354174]]*
  PD[-p$354173][PD[-i][pertF[LI[1], -p$354169]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, p$354174]*
  PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354173][PD[-i][perth[LI[1], -j, -p$354174]]] + 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-i][perth[LI[1], -p$354170, -p$354174]]*
  PD[-p$354173][PD[-j][pertF[LI[1], -p$354169]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-p$354170][perth[LI[1], -i, -p$354174]]*
  PD[-p$354173][PD[-j][pertF[LI[1], -p$354169]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, p$354174]*
  PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354173][PD[-j][perth[LI[1], -i, -p$354174]]] + 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-p$354173][PD[-j][pertF[LI[1], -p$354169]]]*
  PD[-p$354174][perth[LI[1], -i, -p$354170]] + 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-p$354173][PD[-i][pertF[LI[1], -p$354169]]]*
  PD[-p$354174][perth[LI[1], -j, -p$354170]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354170]*
  metric\[Delta][p$354173, p$354174]*perth[LI[1], -p$354169, -p$354173]*
  timevec[h$354171]*PD[-p$354174][PD[-h$354171][pertF[LI[1], -p$354170]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354170]*
  metric\[Delta][p$354173, p$354174]*perth[LI[1], -p$354169, -p$354173]*
  timevec[h$354171]*timevec[h$354172]*PD[-p$354174][
   PD[-h$354172][PD[-h$354171][pertF[LI[1], -p$354170]]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, p$354174]*
  PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354174][PD[-i][perth[LI[1], -j, -p$354173]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, p$354174]*
  PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354174][PD[-j][perth[LI[1], -i, -p$354173]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-i][perth[LI[1], -j, -p$354170]]*PD[-p$354174][
   PD[-p$354173][pertF[LI[1], -p$354169]]] - 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-j][perth[LI[1], -i, -p$354170]]*PD[-p$354174][
   PD[-p$354173][pertF[LI[1], -p$354169]]] + 
 Mpl^2*metric\[Delta][p$354169, p$354170]*metric\[Delta][p$354173, p$354174]*
  PD[-p$354170][perth[LI[1], -i, -j]]*PD[-p$354174][
   PD[-p$354173][pertF[LI[1], -p$354169]]] + 
 2*Mpl^2*metric\[Delta][p$354169, p$354173]*metric\[Delta][p$354170, 
   p$354174]*PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354174][PD[-p$354173][perth[LI[1], -i, -j]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354170]*
  metric\[Delta][p$354173, p$354174]*metric\[Delta][p$354175, p$354176]*
  PD[-p$354170][perth[LI[1], -p$354174, -p$354176]]*
  PD[-p$354175][PD[-p$354173][pertF[LI[1], -p$354169]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354170]*
  metric\[Delta][p$354173, p$354174]*metric\[Delta][p$354175, p$354176]*
  PD[-p$354175][PD[-p$354173][pertF[LI[1], -p$354169]]]*
  PD[-p$354176][perth[LI[1], -p$354170, -p$354174]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354169, p$354173]*
  metric\[Delta][p$354170, p$354174]*metric\[Delta][p$354175, p$354176]*
  PD[-p$354170][pertF[LI[1], -p$354169]]*
  PD[-p$354176][PD[-p$354175][perth[LI[1], -p$354173, -p$354174]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354340, p$354341]*
  PD[-p$354340][pertphi[LI[1]]]*PD[-p$354341][pertphi[LI[1]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354340, p$354341]*
  pertphi[LI[1]]*PD[-p$354341][PD[-p$354340][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*
  PD[-j][PD[-i][pertF[LI[1], -p$354935]]]*PD[-p$354934][pertphi[LI[1]]] - 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-p$354934][pertphi[LI[1]]]*
  PD[-p$354935][PD[-i][pertF[LI[1], -j]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-j][pertF[LI[1], -p$354934]]*
  PD[-p$354935][PD[-i][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-p$354934][pertF[LI[1], -j]]*
  PD[-p$354935][PD[-i][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-p$354934][pertphi[LI[1]]]*
  PD[-p$354935][PD[-j][pertF[LI[1], -i]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-i][pertF[LI[1], -p$354934]]*
  PD[-p$354935][PD[-j][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-p$354934][pertF[LI[1], -i]]*
  PD[-p$354935][PD[-j][pertphi[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-j][pertphi[LI[1]]]*
  PD[-p$354935][PD[-p$354934][pertF[LI[1], -i]]] + 
 2*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-i][pertphi[LI[1]]]*
  PD[-p$354935][PD[-p$354934][pertF[LI[1], -j]]] - 
 4*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-i][pertF[LI[1], -j]]*
  PD[-p$354935][PD[-p$354934][pertphi[LI[1]]]] - 
 4*Mpl^2*metric\[Delta][p$354934, p$354935]*PD[-j][pertF[LI[1], -i]]*
  PD[-p$354935][PD[-p$354934][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$354934, p$354935]*
  metric\[Delta][p$354936, p$354937]*PD[-p$354934][pertphi[LI[1]]]*
  PD[-p$354937][PD[-p$354936][pertF[LI[1], -p$354935]]] - 
 Mpl^2*metric\[Delta][p$355733, p$355734]*timevec[h$355735]*timevec[h$355736]*
  PD[-h$355736][PD[-j][pertF[LI[1], -p$355734]]]*
  PD[-p$355733][PD[-h$355735][pertF[LI[1], -i]]] - 
 Mpl^2*metric\[Delta][p$355733, p$355734]*timevec[h$355735]*timevec[h$355736]*
  PD[-h$355736][PD[-i][pertF[LI[1], -p$355734]]]*
  PD[-p$355733][PD[-h$355735][pertF[LI[1], -j]]] - 
 Mpl^2*metric\[Delta][p$355733, p$355734]*timevec[h$355735]*timevec[h$355736]*
  PD[-p$355733][PD[-h$355735][pertF[LI[1], -i]]]*
  PD[-p$355734][PD[-h$355736][pertF[LI[1], -j]]] + 
 (3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355734]*
   metric\[Delta][p$355737, p$355738]*timevec[h$355735]*timevec[h$355736]*
   PD[-p$355734][PD[-h$355736][pertF[LI[1], -p$355738]]]*
   PD[-p$355737][PD[-h$355735][pertF[LI[1], -p$355733]]])/2 + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355737]*
  metric\[Delta][p$355734, p$355738]*timevec[h$355735]*
  PD[-p$355734][pertF[LI[1], -p$355733]]*
  PD[-p$355737][PD[-h$355735][pertF[LI[1], -p$355738]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355737]*
  metric\[Delta][p$355734, p$355738]*timevec[h$355735]*timevec[h$355736]*
  PD[-p$355734][pertF[LI[1], -p$355733]]*
  PD[-p$355737][PD[-h$355736][PD[-h$355735][pertF[LI[1], -p$355738]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355737]*
  metric\[Delta][p$355734, p$355738]*timevec[h$355735]*
  PD[-p$355734][pertF[LI[1], -p$355733]]*
  PD[-p$355738][PD[-h$355735][pertF[LI[1], -p$355737]]] + 
 (3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355734]*
   metric\[Delta][p$355737, p$355738]*timevec[h$355735]*timevec[h$355736]*
   PD[-p$355737][PD[-h$355735][pertF[LI[1], -p$355733]]]*
   PD[-p$355738][PD[-h$355736][pertF[LI[1], -p$355734]]])/2 + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355737]*
  metric\[Delta][p$355734, p$355738]*timevec[h$355735]*timevec[h$355736]*
  PD[-p$355734][pertF[LI[1], -p$355733]]*
  PD[-p$355738][PD[-h$355736][PD[-h$355735][pertF[LI[1], -p$355737]]]] + 
 2*Mpl^2*metric\[Delta][p$355733, p$355734]*metric\[Delta][p$355737, 
   p$355738]*PD[-p$355737][PD[-i][pertF[LI[1], -p$355733]]]*
  PD[-p$355738][PD[-j][pertF[LI[1], -p$355734]]] - 
 2*Mpl^2*metric\[Delta][p$355733, p$355734]*metric\[Delta][p$355737, 
   p$355738]*PD[-j][PD[-i][pertF[LI[1], -p$355733]]]*
  PD[-p$355738][PD[-p$355737][pertF[LI[1], -p$355734]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355734]*
  metric\[Delta][p$355737, p$355738]*metric\[Delta][p$355739, p$355740]*
  PD[-p$355739][PD[-p$355737][pertF[LI[1], -p$355733]]]*
  PD[-p$355740][PD[-p$355738][pertF[LI[1], -p$355734]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$355733, p$355734]*
  metric\[Delta][p$355737, p$355738]*metric\[Delta][p$355739, p$355740]*
  PD[-p$355738][PD[-p$355737][pertF[LI[1], -p$355733]]]*
  PD[-p$355740][PD[-p$355739][pertF[LI[1], -p$355734]]] + 
 4*Mpl^2*metric\[Delta][p$356255, p$356256]*
  PD[-p$356255][PD[-j][pertE[LI[1]]]]*PD[-p$356256][PD[-i][pertphi[LI[1]]]] + 
 4*Mpl^2*metric\[Delta][p$356255, p$356256]*
  PD[-p$356255][PD[-i][pertE[LI[1]]]]*PD[-p$356256][PD[-j][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$356255, p$356256]*PD[-p$356255][pertphi[LI[1]]]*
  PD[-p$356256][PD[-j][PD[-i][pertE[LI[1]]]]] - 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$356255, p$356256]*
  timevec[h$356253]*PD[-h$356253][pertphi[LI[1]]]*
  PD[-p$356256][PD[-p$356255][pertE[LI[1]]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$356255, p$356256]*
  timevec[h$356253]*timevec[h$356254]*PD[-h$356254][
   PD[-h$356253][pertphi[LI[1]]]]*PD[-p$356256][
   PD[-p$356255][pertE[LI[1]]]] - 4*Mpl^2*metric\[Delta][p$356255, p$356256]*
  PD[-j][PD[-i][pertphi[LI[1]]]]*PD[-p$356256][PD[-p$356255][pertE[LI[1]]]] - 
 8*Mpl^2*metric\[Delta][p$356255, p$356256]*PD[-j][PD[-i][pertE[LI[1]]]]*
  PD[-p$356256][PD[-p$356255][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$356255, p$356256]*
  timevec[h$356253]*timevec[h$356254]*PD[-h$356253][pertphi[LI[1]]]*
  PD[-p$356256][PD[-p$356255][PD[-h$356254][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$356255, p$356256]*
  metric\[Delta][p$356257, p$356258]*PD[-p$356256][
   PD[-p$356255][pertE[LI[1]]]]*PD[-p$356258][
   PD[-p$356257][pertphi[LI[1]]]] + 2*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$356255, p$356256]*metric\[Delta][p$356257, p$356258]*
  PD[-p$356255][pertphi[LI[1]]]*PD[-p$356258][
   PD[-p$356257][PD[-p$356256][pertE[LI[1]]]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$356528, p$356529]*
  timevec[h$356527]*PD[-p$356528][pertB[LI[1]]]*
  PD[-p$356529][PD[-h$356527][pertphi[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$356528, p$356529]*
  timevec[h$356527]*PD[-h$356527][pertphi[LI[1]]]*
  PD[-p$356529][PD[-p$356528][pertB[LI[1]]]] + 
 (density[]*metric\[Delta][-i, -j]*metric\[Delta][p$357013, p$357014]*
   PD[-p$357014][PD[-p$357013][pertshear[LI[2]]]])/2 + 
 (metric\[Delta][-i, -j]*metric\[Delta][p$357013, p$357014]*pressure[]*
   PD[-p$357014][PD[-p$357013][pertshear[LI[2]]]])/2 - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
  timevec[h$357740]*PD[-h$357740][pertS[LI[1], -p$357739]]*
  PD[-p$357738][pertB[LI[1]]] + 2*Mpl^2*metric\[Delta][p$357738, p$357739]*
  PD[-j][PD[-i][pertS[LI[1], -p$357739]]]*PD[-p$357738][pertB[LI[1]]] - 
 4*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, 
   p$357739]*pertS[LI[1], -p$357738]*PD[-p$357739][pertB[LI[1]]] - 
 8*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
  pertS[LI[1], -p$357738]*primehubbleC[]*PD[-p$357739][pertB[LI[1]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
  pertS[LI[1], -p$357738]*timevec[h$357740]*
  PD[-p$357739][PD[-h$357740][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-j][pertS[LI[1], -p$357738]]*
  PD[-p$357739][PD[-i][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-p$357738][pertS[LI[1], -j]]*
  PD[-p$357739][PD[-i][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-p$357738][pertB[LI[1]]]*
  PD[-p$357739][PD[-i][pertS[LI[1], -j]]] + 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-i][pertS[LI[1], -p$357738]]*
  PD[-p$357739][PD[-j][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-p$357738][pertS[LI[1], -i]]*
  PD[-p$357739][PD[-j][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-p$357738][pertB[LI[1]]]*
  PD[-p$357739][PD[-j][pertS[LI[1], -i]]] - 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-i][pertS[LI[1], -j]]*
  PD[-p$357739][PD[-p$357738][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$357738, p$357739]*PD[-j][pertS[LI[1], -i]]*
  PD[-p$357739][PD[-p$357738][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357741]*
  metric\[Delta][p$357739, p$357742]*PD[-p$357739][pertS[LI[1], -p$357738]]*
  PD[-p$357742][PD[-p$357741][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
  metric\[Delta][p$357741, p$357742]*PD[-p$357738][pertB[LI[1]]]*
  PD[-p$357742][PD[-p$357741][pertS[LI[1], -p$357739]]] + 
 Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-h$358367][PD[-i][perth[LI[1], -j, -p$358366]]]*
  PD[-p$358365][pertB[LI[1]]] + Mpl^2*metric\[Delta][p$358365, p$358366]*
  timevec[h$358367]*PD[-h$358367][PD[-j][perth[LI[1], -i, -p$358366]]]*
  PD[-p$358365][pertB[LI[1]]] + 2*Mpl^2*hubbleC[]*
  metric\[Delta][p$358365, p$358366]*PD[-i][perth[LI[1], -j, -p$358366]]*
  PD[-p$358365][pertB[LI[1]]] + 2*Mpl^2*hubbleC[]*
  metric\[Delta][p$358365, p$358366]*PD[-j][perth[LI[1], -i, -p$358366]]*
  PD[-p$358365][pertB[LI[1]]] + Mpl^2*metric\[Delta][p$358365, p$358366]*
  timevec[h$358367]*PD[-i][perth[LI[1], -j, -p$358366]]*
  PD[-p$358365][PD[-h$358367][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-j][perth[LI[1], -i, -p$358366]]*PD[-p$358365][
   PD[-h$358367][pertB[LI[1]]]] + Mpl^2*metric\[Delta][p$358365, p$358366]*
  timevec[h$358367]*PD[-h$358367][perth[LI[1], -j, -p$358366]]*
  PD[-p$358365][PD[-i][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-h$358367][perth[LI[1], -i, -p$358366]]*
  PD[-p$358365][PD[-j][pertB[LI[1]]]] - 2*Mpl^2*hubbleC[]*
  metric\[Delta][p$358365, p$358366]*PD[-p$358365][pertB[LI[1]]]*
  PD[-p$358366][perth[LI[1], -i, -j]] - 
 Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-p$358365][PD[-h$358367][pertB[LI[1]]]]*
  PD[-p$358366][perth[LI[1], -i, -j]] - 
 2*Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-p$358365][pertB[LI[1]]]*PD[-p$358366][
   PD[-h$358367][perth[LI[1], -i, -j]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$358365, p$358366]*perth[LI[1], -i, -j]*
  PD[-p$358366][PD[-p$358365][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$358365, p$358366]*timevec[h$358367]*
  PD[-h$358367][perth[LI[1], -i, -j]]*PD[-p$358366][
   PD[-p$358365][pertB[LI[1]]]] + 2*Mpl^2*metric\[Delta][p$358365, p$358366]*
  perth[LI[1], -i, -j]*timevec[h$358367]*
  PD[-p$358366][PD[-p$358365][PD[-h$358367][pertB[LI[1]]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$358365, p$358366]*
  metric\[Delta][p$358368, p$358369]*timevec[h$358367]*
  PD[-h$358367][perth[LI[1], -p$358366, -p$358369]]*
  PD[-p$358368][PD[-p$358365][pertB[LI[1]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$358365, p$358366]*
  metric\[Delta][p$358368, p$358369]*perth[LI[1], -p$358365, -p$358368]*
  PD[-p$358369][PD[-p$358366][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$358365, p$358366]*
  metric\[Delta][p$358368, p$358369]*perth[LI[1], -p$358365, -p$358368]*
  timevec[h$358367]*PD[-p$358369][PD[-p$358366][
    PD[-h$358367][pertB[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-h$359261][PD[-j][PD[-i][pertF[LI[1], -p$359260]]]]*
  PD[-p$359259][pertB[LI[1]]] + 4*Mpl^2*hubbleC[]*
  metric\[Delta][p$359259, p$359260]*PD[-j][PD[-i][pertF[LI[1], -p$359260]]]*
  PD[-p$359259][pertB[LI[1]]] + 2*Mpl^2*metric\[Delta][p$359259, p$359260]*
  timevec[h$359261]*PD[-j][PD[-i][pertF[LI[1], -p$359260]]]*
  PD[-p$359259][PD[-h$359261][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-h$359261][PD[-j][pertF[LI[1], -p$359260]]]*
  PD[-p$359259][PD[-i][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-h$359261][PD[-i][pertF[LI[1], -p$359260]]]*
  PD[-p$359259][PD[-j][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-p$359259][PD[-j][pertB[LI[1]]]]*PD[-p$359260][
   PD[-h$359261][pertF[LI[1], -i]]] + 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-p$359259][PD[-i][pertB[LI[1]]]]*PD[-p$359260][
   PD[-h$359261][pertF[LI[1], -j]]] - 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-p$359259][pertB[LI[1]]]*PD[-p$359260][
   PD[-h$359261][PD[-i][pertF[LI[1], -j]]]] - 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-p$359259][pertB[LI[1]]]*PD[-p$359260][
   PD[-h$359261][PD[-j][pertF[LI[1], -i]]]] - 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-h$359261][PD[-i][pertF[LI[1], -j]]]*
  PD[-p$359260][PD[-p$359259][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-h$359261][PD[-j][pertF[LI[1], -i]]]*
  PD[-p$359260][PD[-p$359259][pertB[LI[1]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$359259, p$359260]*
  PD[-i][pertF[LI[1], -j]]*PD[-p$359260][PD[-p$359259][pertB[LI[1]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$359259, p$359260]*
  PD[-j][pertF[LI[1], -i]]*PD[-p$359260][PD[-p$359259][pertB[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-i][pertF[LI[1], -j]]*PD[-p$359260][PD[-p$359259][
    PD[-h$359261][pertB[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$359259, p$359260]*timevec[h$359261]*
  PD[-j][pertF[LI[1], -i]]*PD[-p$359260][PD[-p$359259][
    PD[-h$359261][pertB[LI[1]]]]] - 2*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$359259, p$359260]*metric\[Delta][p$359262, p$359263]*
  timevec[h$359261]*PD[-p$359262][PD[-p$359259][pertB[LI[1]]]]*
  PD[-p$359263][PD[-h$359261][pertF[LI[1], -p$359260]]] - 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$359259, p$359262]*
  metric\[Delta][p$359260, p$359263]*PD[-p$359260][pertF[LI[1], -p$359259]]*
  PD[-p$359263][PD[-p$359262][pertB[LI[1]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$359259, p$359260]*
  metric\[Delta][p$359262, p$359263]*PD[-p$359259][pertB[LI[1]]]*
  PD[-p$359263][PD[-p$359262][pertF[LI[1], -p$359260]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359259, p$359260]*
  metric\[Delta][p$359262, p$359263]*timevec[h$359261]*
  PD[-p$359259][PD[-h$359261][pertB[LI[1]]]]*
  PD[-p$359263][PD[-p$359262][pertF[LI[1], -p$359260]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359259, p$359262]*
  metric\[Delta][p$359260, p$359263]*timevec[h$359261]*
  PD[-p$359260][pertF[LI[1], -p$359259]]*
  PD[-p$359263][PD[-p$359262][PD[-h$359261][pertB[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359259, p$359260]*
  metric\[Delta][p$359262, p$359263]*timevec[h$359261]*
  PD[-p$359259][pertB[LI[1]]]*PD[-p$359263][
   PD[-p$359262][PD[-h$359261][pertF[LI[1], -p$359260]]]] + 
 2*Mpl^2*hubbleC[]^2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, 
   p$359827]*PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
  primehubbleC[]*PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
  timevec[h$359828]*PD[-p$359826][pertB[LI[1]]]*
  PD[-p$359827][PD[-h$359828][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$359826, p$359827]*
  PD[-p$359826][PD[-i][pertB[LI[1]]]]*PD[-p$359827][PD[-j][pertB[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$359826, p$359827]*PD[-j][PD[-i][pertB[LI[1]]]]*
  PD[-p$359827][PD[-p$359826][pertB[LI[1]]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
  metric\[Delta][p$359829, p$359830]*PD[-p$359829][
   PD[-p$359826][pertB[LI[1]]]]*PD[-p$359830][PD[-p$359827][pertB[LI[1]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
  metric\[Delta][p$359829, p$359830]*PD[-p$359827][
   PD[-p$359826][pertB[LI[1]]]]*PD[-p$359830][PD[-p$359829][pertB[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$360275, p$360276]*timevec[h$360277]*
  PD[-p$360275][pertS[LI[1], -j]]*PD[-p$360276][
   PD[-h$360277][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$360275, p$360276]*timevec[h$360277]*
  PD[-p$360275][pertS[LI[1], -i]]*PD[-p$360276][
   PD[-h$360277][PD[-j][pertE[LI[1]]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$360275, p$360276]*pertS[LI[1], -p$360275]*
  PD[-p$360276][PD[-j][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$360275, p$360276]*timevec[h$360277]*
  PD[-h$360277][pertS[LI[1], -p$360275]]*
  PD[-p$360276][PD[-j][PD[-i][pertE[LI[1]]]]] + 
 Mpl^2*metric\[Delta][p$360275, p$360276]*timevec[h$360277]*
  PD[-i][pertS[LI[1], -j]]*PD[-p$360276][PD[-p$360275][
    PD[-h$360277][pertE[LI[1]]]]] + Mpl^2*metric\[Delta][p$360275, p$360276]*
  timevec[h$360277]*PD[-j][pertS[LI[1], -i]]*
  PD[-p$360276][PD[-p$360275][PD[-h$360277][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360275, p$360276]*
  metric\[Delta][p$360278, p$360279]*timevec[h$360277]*
  PD[-p$360278][PD[-p$360275][pertE[LI[1]]]]*
  PD[-p$360279][PD[-h$360277][pertS[LI[1], -p$360276]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$360275, p$360278]*
  metric\[Delta][p$360276, p$360279]*PD[-p$360276][pertS[LI[1], -p$360275]]*
  PD[-p$360279][PD[-p$360278][pertE[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360275, p$360278]*
  metric\[Delta][p$360276, p$360279]*timevec[h$360277]*
  PD[-p$360276][pertS[LI[1], -p$360275]]*
  PD[-p$360279][PD[-p$360278][PD[-h$360277][pertE[LI[1]]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$360275, p$360276]*
  metric\[Delta][p$360278, p$360279]*pertS[LI[1], -p$360275]*
  PD[-p$360279][PD[-p$360278][PD[-p$360276][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360275, p$360276]*
  metric\[Delta][p$360278, p$360279]*timevec[h$360277]*
  PD[-h$360277][pertS[LI[1], -p$360275]]*
  PD[-p$360279][PD[-p$360278][PD[-p$360276][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-p$360899][PD[-j][pertB[LI[1]]]]*PD[-p$360900][
   PD[-h$360901][PD[-i][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-p$360899][PD[-i][pertB[LI[1]]]]*PD[-p$360900][
   PD[-h$360901][PD[-j][pertE[LI[1]]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$360899, p$360900]*
  PD[-p$360899][pertB[LI[1]]]*PD[-p$360900][PD[-j][PD[-i][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-p$360899][PD[-h$360901][pertB[LI[1]]]]*
  PD[-p$360900][PD[-j][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-h$360901][PD[-j][PD[-i][pertE[LI[1]]]]]*
  PD[-p$360900][PD[-p$360899][pertB[LI[1]]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][p$360899, p$360900]*
  PD[-j][PD[-i][pertE[LI[1]]]]*PD[-p$360900][PD[-p$360899][pertB[LI[1]]]] + 
 4*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-j][PD[-i][pertE[LI[1]]]]*PD[-p$360900][
   PD[-p$360899][PD[-h$360901][pertB[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$360899, p$360900]*timevec[h$360901]*
  PD[-j][PD[-i][pertB[LI[1]]]]*PD[-p$360900][
   PD[-p$360899][PD[-h$360901][pertE[LI[1]]]]] - 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*PD[-p$360902][
   PD[-p$360899][pertB[LI[1]]]]*PD[-p$360903][PD[-p$360900][pertE[LI[1]]]] - 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*timevec[h$360901]*
  PD[-p$360902][PD[-p$360899][pertE[LI[1]]]]*
  PD[-p$360903][PD[-p$360900][PD[-h$360901][pertB[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*timevec[h$360901]*
  PD[-p$360902][PD[-p$360899][pertB[LI[1]]]]*
  PD[-p$360903][PD[-p$360900][PD[-h$360901][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*timevec[h$360901]*
  PD[-p$360900][PD[-p$360899][pertB[LI[1]]]]*
  PD[-p$360903][PD[-p$360902][PD[-h$360901][pertE[LI[1]]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*PD[-p$360899][pertB[LI[1]]]*
  PD[-p$360903][PD[-p$360902][PD[-p$360900][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$360899, p$360900]*
  metric\[Delta][p$360902, p$360903]*timevec[h$360901]*
  PD[-p$360899][PD[-h$360901][pertB[LI[1]]]]*
  PD[-p$360903][PD[-p$360902][PD[-p$360900][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*timevec[h$361963]*
  timevec[h$361964]*PD[-h$361963][perth[LI[1], -j, -p$361961]]*
  PD[-p$361962][PD[-h$361964][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*timevec[h$361963]*
  timevec[h$361964]*PD[-h$361963][perth[LI[1], -i, -p$361961]]*
  PD[-p$361962][PD[-h$361964][PD[-j][pertE[LI[1]]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$361961, p$361962]*perth[LI[1], -i, -j]*
  timevec[h$361963]*PD[-p$361962][PD[-p$361961][
    PD[-h$361963][pertE[LI[1]]]]] + Mpl^2*metric\[Delta][p$361961, p$361962]*
  timevec[h$361963]*timevec[h$361964]*PD[-h$361963][perth[LI[1], -i, -j]]*
  PD[-p$361962][PD[-p$361961][PD[-h$361964][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*perth[LI[1], -i, -j]*
  timevec[h$361963]*timevec[h$361964]*PD[-p$361962][
   PD[-p$361961][PD[-h$361964][PD[-h$361963][pertE[LI[1]]]]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*timevec[h$361963]*
  PD[-h$361963][perth[LI[1], -p$361962, -p$361966]]*
  PD[-p$361965][PD[-p$361961][pertE[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*timevec[h$361963]*timevec[h$361964]*
  PD[-h$361964][PD[-h$361963][perth[LI[1], -p$361962, -p$361966]]]*
  PD[-p$361965][PD[-p$361961][pertE[LI[1]]]] + 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, 
   p$361966]*PD[-j][PD[-i][perth[LI[1], -p$361962, -p$361966]]]*
  PD[-p$361965][PD[-p$361961][pertE[LI[1]]]] - 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, 
   p$361966]*PD[-p$361965][PD[-p$361961][pertE[LI[1]]]]*
  PD[-p$361966][PD[-i][perth[LI[1], -j, -p$361962]]] - 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, 
   p$361966]*PD[-p$361965][PD[-p$361961][pertE[LI[1]]]]*
  PD[-p$361966][PD[-j][perth[LI[1], -i, -p$361962]]] + 
 2*Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, 
   p$361966]*PD[-p$361965][PD[-p$361961][pertE[LI[1]]]]*
  PD[-p$361966][PD[-p$361962][perth[LI[1], -i, -j]]] + 
 4*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*perth[LI[1], -p$361961, -p$361965]*
  timevec[h$361963]*PD[-p$361966][PD[-p$361962][
    PD[-h$361963][pertE[LI[1]]]]] + 3*Mpl^2*metric\[Delta][-i, -j]*
  metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  timevec[h$361963]*timevec[h$361964]*PD[-h$361963][
   perth[LI[1], -p$361961, -p$361965]]*
  PD[-p$361966][PD[-p$361962][PD[-h$361964][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*perth[LI[1], -p$361961, -p$361965]*
  timevec[h$361963]*timevec[h$361964]*PD[-p$361966][
   PD[-p$361962][PD[-h$361964][PD[-h$361963][pertE[LI[1]]]]]] + 
 Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  PD[-j][perth[LI[1], -p$361961, -p$361965]]*
  PD[-p$361966][PD[-p$361962][PD[-i][pertE[LI[1]]]]] + 
 Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  PD[-i][perth[LI[1], -p$361961, -p$361965]]*
  PD[-p$361966][PD[-p$361962][PD[-j][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  PD[-i][perth[LI[1], -j, -p$361961]]*PD[-p$361966][
   PD[-p$361965][PD[-p$361962][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  PD[-j][perth[LI[1], -i, -p$361961]]*PD[-p$361966][
   PD[-p$361965][PD[-p$361962][pertE[LI[1]]]]] + 
 Mpl^2*metric\[Delta][p$361961, p$361962]*metric\[Delta][p$361965, p$361966]*
  PD[-p$361961][perth[LI[1], -i, -j]]*PD[-p$361966][
   PD[-p$361965][PD[-p$361962][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*metric\[Delta][p$361967, p$361968]*
  PD[-p$361967][perth[LI[1], -p$361961, -p$361965]]*
  PD[-p$361968][PD[-p$361966][PD[-p$361962][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$361961, p$361962]*
  metric\[Delta][p$361965, p$361966]*metric\[Delta][p$361967, p$361968]*
  PD[-p$361965][PD[-p$361961][pertE[LI[1]]]]*
  PD[-p$361968][PD[-p$361967][perth[LI[1], -p$361962, -p$361966]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-h$363059][PD[-j][pertF[LI[1], -p$363057]]]*
  PD[-p$363058][PD[-h$363060][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-p$363057][PD[-h$363059][pertF[LI[1], -j]]]*
  PD[-p$363058][PD[-h$363060][PD[-i][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-h$363059][PD[-i][pertF[LI[1], -p$363057]]]*
  PD[-p$363058][PD[-h$363060][PD[-j][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-p$363057][PD[-h$363059][pertF[LI[1], -i]]]*
  PD[-p$363058][PD[-h$363060][PD[-j][pertE[LI[1]]]]] - 
 4*Mpl^2*hubbleC[]*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  PD[-i][pertF[LI[1], -j]]*PD[-p$363058][PD[-p$363057][
    PD[-h$363059][pertE[LI[1]]]]] - 4*Mpl^2*hubbleC[]*
  metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  PD[-j][pertF[LI[1], -i]]*PD[-p$363058][PD[-p$363057][
    PD[-h$363059][pertE[LI[1]]]]] + Mpl^2*metric\[Delta][p$363057, p$363058]*
  timevec[h$363059]*timevec[h$363060]*PD[-h$363059][PD[-i][pertF[LI[1], -j]]]*
  PD[-p$363058][PD[-p$363057][PD[-h$363060][pertE[LI[1]]]]] + 
 Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*timevec[h$363060]*
  PD[-h$363059][PD[-j][pertF[LI[1], -i]]]*
  PD[-p$363058][PD[-p$363057][PD[-h$363060][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-i][pertF[LI[1], -j]]*
  PD[-p$363058][PD[-p$363057][PD[-h$363060][PD[-h$363059][pertE[LI[1]]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*timevec[h$363059]*
  timevec[h$363060]*PD[-j][pertF[LI[1], -i]]*
  PD[-p$363058][PD[-p$363057][PD[-h$363060][PD[-h$363059][pertE[LI[1]]]]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363058]*
  metric\[Delta][p$363061, p$363062]*timevec[h$363059]*
  PD[-p$363061][PD[-p$363057][pertE[LI[1]]]]*
  PD[-p$363062][PD[-h$363059][pertF[LI[1], -p$363058]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363058]*
  metric\[Delta][p$363061, p$363062]*timevec[h$363059]*timevec[h$363060]*
  PD[-p$363061][PD[-p$363057][pertE[LI[1]]]]*
  PD[-p$363062][PD[-h$363060][PD[-h$363059][pertF[LI[1], -p$363058]]]] + 
 6*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363058]*
  metric\[Delta][p$363061, p$363062]*timevec[h$363059]*timevec[h$363060]*
  PD[-p$363061][PD[-h$363059][pertF[LI[1], -p$363057]]]*
  PD[-p$363062][PD[-p$363058][PD[-h$363060][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*metric\[Delta][p$363061, 
   p$363062]*PD[-p$363061][PD[-j][pertF[LI[1], -p$363057]]]*
  PD[-p$363062][PD[-p$363058][PD[-i][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*metric\[Delta][p$363061, 
   p$363062]*PD[-p$363061][PD[-i][pertF[LI[1], -p$363057]]]*
  PD[-p$363062][PD[-p$363058][PD[-j][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*metric\[Delta][p$363061, 
   p$363062]*PD[-p$363058][PD[-j][PD[-i][pertE[LI[1]]]]]*
  PD[-p$363062][PD[-p$363061][pertF[LI[1], -p$363057]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363061]*
  metric\[Delta][p$363058, p$363062]*timevec[h$363059]*
  PD[-p$363058][pertF[LI[1], -p$363057]]*
  PD[-p$363062][PD[-p$363061][PD[-h$363059][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363061]*
  metric\[Delta][p$363058, p$363062]*timevec[h$363059]*timevec[h$363060]*
  PD[-p$363058][pertF[LI[1], -p$363057]]*
  PD[-p$363062][PD[-p$363061][PD[-h$363060][PD[-h$363059][pertE[LI[1]]]]]] - 
 2*Mpl^2*metric\[Delta][p$363057, p$363058]*metric\[Delta][p$363061, 
   p$363062]*PD[-j][PD[-i][pertF[LI[1], -p$363057]]]*
  PD[-p$363062][PD[-p$363061][PD[-p$363058][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363058]*
  metric\[Delta][p$363061, p$363062]*metric\[Delta][p$363063, p$363064]*
  PD[-p$363063][PD[-p$363061][pertF[LI[1], -p$363057]]]*
  PD[-p$363064][PD[-p$363062][PD[-p$363058][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363057, p$363058]*
  metric\[Delta][p$363061, p$363062]*metric\[Delta][p$363063, p$363064]*
  PD[-p$363062][PD[-p$363061][pertF[LI[1], -p$363057]]]*
  PD[-p$363064][PD[-p$363063][PD[-p$363058][pertE[LI[1]]]]] - 
 4*Mpl^2*metric\[Delta][p$363766, p$363767]*timevec[h$363768]*
  timevec[h$363769]*PD[-p$363766][PD[-h$363768][PD[-i][pertE[LI[1]]]]]*
  PD[-p$363767][PD[-h$363769][PD[-j][pertE[LI[1]]]]] - 
 8*Mpl^2*hubbleC[]*metric\[Delta][p$363766, p$363767]*timevec[h$363768]*
  PD[-j][PD[-i][pertE[LI[1]]]]*PD[-p$363767][
   PD[-p$363766][PD[-h$363768][pertE[LI[1]]]]] + 
 2*Mpl^2*metric\[Delta][p$363766, p$363767]*timevec[h$363768]*
  timevec[h$363769]*PD[-h$363768][PD[-j][PD[-i][pertE[LI[1]]]]]*
  PD[-p$363767][PD[-p$363766][PD[-h$363769][pertE[LI[1]]]]] - 
 4*Mpl^2*metric\[Delta][p$363766, p$363767]*timevec[h$363768]*
  timevec[h$363769]*PD[-j][PD[-i][pertE[LI[1]]]]*
  PD[-p$363767][PD[-p$363766][PD[-h$363769][PD[-h$363768][pertE[LI[1]]]]]] + 
 8*Mpl^2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*timevec[h$363768]*
  PD[-p$363770][PD[-p$363766][pertE[LI[1]]]]*
  PD[-p$363771][PD[-p$363767][PD[-h$363768][pertE[LI[1]]]]] + 
 3*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*timevec[h$363768]*timevec[h$363769]*
  PD[-p$363770][PD[-p$363766][PD[-h$363768][pertE[LI[1]]]]]*
  PD[-p$363771][PD[-p$363767][PD[-h$363769][pertE[LI[1]]]]] + 
 4*Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*timevec[h$363768]*timevec[h$363769]*
  PD[-p$363770][PD[-p$363766][pertE[LI[1]]]]*
  PD[-p$363771][PD[-p$363767][PD[-h$363769][PD[-h$363768][pertE[LI[1]]]]]] + 
 2*Mpl^2*metric\[Delta][p$363766, p$363767]*metric\[Delta][p$363770, 
   p$363771]*PD[-p$363770][PD[-p$363766][PD[-i][pertE[LI[1]]]]]*
  PD[-p$363771][PD[-p$363767][PD[-j][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*timevec[h$363768]*timevec[h$363769]*
  PD[-p$363767][PD[-p$363766][PD[-h$363768][pertE[LI[1]]]]]*
  PD[-p$363771][PD[-p$363770][PD[-h$363769][pertE[LI[1]]]]] - 
 2*Mpl^2*metric\[Delta][p$363766, p$363767]*metric\[Delta][p$363770, 
   p$363771]*PD[-p$363766][PD[-j][PD[-i][pertE[LI[1]]]]]*
  PD[-p$363771][PD[-p$363770][PD[-p$363767][pertE[LI[1]]]]] - 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*metric\[Delta][p$363772, p$363773]*
  PD[-p$363772][PD[-p$363770][PD[-p$363766][pertE[LI[1]]]]]*
  PD[-p$363773][PD[-p$363771][PD[-p$363767][pertE[LI[1]]]]] + 
 Mpl^2*metric\[Delta][-i, -j]*metric\[Delta][p$363766, p$363767]*
  metric\[Delta][p$363770, p$363771]*metric\[Delta][p$363772, p$363773]*
  PD[-p$363770][PD[-p$363767][PD[-p$363766][pertE[LI[1]]]]]*
  PD[-p$363773][PD[-p$363772][PD[-p$363771][pertE[LI[1]]]]] + 
 2*perth[LI[1], -i, -j]*pertpsi[LI[1]]*primescalar[]^2*
  Derivative[0, 1][G2fun][scalar[], X[]] - 4*metric\[Delta][-i, -j]*
  pertphi[LI[1]]*pertpsi[LI[1]]*primescalar[]^2*
  Derivative[0, 1][G2fun][scalar[], X[]] - 4*metric\[Delta][-i, -j]*
  pertpsi[LI[1]]^2*primescalar[]^2*Derivative[0, 1][G2fun][scalar[], X[]] + 
 metric\[Delta][-i, -j]*pertpsi[LI[2]]*primescalar[]^2*
  Derivative[0, 1][G2fun][scalar[], X[]] + metric\[Delta][-i, -j]*
  metric\[Delta][p$347314, p$347315]*pertS[LI[1], -p$347314]*
  pertS[LI[1], -p$347315]*primescalar[]^2*Derivative[0, 1][G2fun][scalar[], 
   X[]] - 2*perth[LI[1], -i, -j]*primescalar[]*timevec[h$364063]*
  PD[-h$364063][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] + 
 4*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]*timevec[h$364378]*
  PD[-h$364378][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] + 
 4*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]*timevec[h$365084]*
  PD[-h$365084][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] - 
 metric\[Delta][-i, -j]*primescalar[]*timevec[h$365348]*
  PD[-h$365348][pertscalar[LI[2]]]*Derivative[0, 1][G2fun][scalar[], X[]] - 
 metric\[Delta][-i, -j]*timevec[h$366624]*timevec[h$366625]*
  PD[-h$366624][pertscalar[LI[1]]]*PD[-h$366625][pertscalar[LI[1]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] + 2*pertpsi[LI[1]]*primescalar[]^2*
  PD[-i][pertF[LI[1], -j]]*Derivative[0, 1][G2fun][scalar[], X[]] - 
 2*primescalar[]*timevec[h$367275]*PD[-h$367275][pertscalar[LI[1]]]*
  PD[-i][pertF[LI[1], -j]]*Derivative[0, 1][G2fun][scalar[], X[]] + 
 2*pertpsi[LI[1]]*primescalar[]^2*PD[-j][pertF[LI[1], -i]]*
  Derivative[0, 1][G2fun][scalar[], X[]] - 2*primescalar[]*timevec[h$367275]*
  PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
  Derivative[0, 1][G2fun][scalar[], X[]] - 2*PD[-i][pertscalar[LI[1]]]*
  PD[-j][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] + 
 4*pertpsi[LI[1]]*primescalar[]^2*PD[-j][PD[-i][pertE[LI[1]]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] - 4*primescalar[]*timevec[h$367638]*
  PD[-h$367638][pertscalar[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] - 2*metric\[Delta][-i, -j]*
  metric\[Delta][p$357738, p$357739]*pertS[LI[1], -p$357738]*primescalar[]^2*
  PD[-p$357739][pertB[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] + 
 metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*primescalar[]^2*
  PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] + metric\[Delta][-i, -j]*
  metric\[Delta][p$366626, p$366627]*PD[-p$366626][pertscalar[LI[1]]]*
  PD[-p$366627][pertscalar[LI[1]]]*Derivative[0, 1][G2fun][scalar[], X[]] - 
 2*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
  pertS[LI[1], -p$367889]*primescalar[]*PD[-p$367890][pertscalar[LI[1]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] + 2*metric\[Delta][-i, -j]*
  metric\[Delta][p$368180, p$368181]*primescalar[]*
  PD[-p$368180][pertB[LI[1]]]*PD[-p$368181][pertscalar[LI[1]]]*
  Derivative[0, 1][G2fun][scalar[], X[]] + 
 (perth[LI[2], -i, -j]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*pertphi[LI[2]]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (8*perth[LI[1], -i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (16*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*pprimescalar[]*
   primescalar[]^2*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (24*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*pertpsi[LI[2]]*pprimescalar[]*primescalar[]^2*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315]*pprimescalar[]*
   primescalar[]^2*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*perth[LI[2], -i, -j]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[2]]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (8*hubbleC[]*perth[LI[1], -i, -j]*pertpsi[LI[1]]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (16*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*
   primescalar[]^3*Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (24*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[2]]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315]*primescalar[]^3*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*perth[LI[1], -i, -j]*primescalar[]^3*timevec[h$345683]*
   PD[-h$345683][pertpsi[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^3*
   timevec[h$346074]*PD[-h$346074][pertpsi[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (12*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^3*timevec[h$346451]*
   PD[-h$346451][pertpsi[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (metric\[Delta][-i, -j]*primescalar[]^3*timevec[h$346694]*
   PD[-h$346694][pertpsi[LI[2]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*primescalar[]^3*timevec[h$347316]*
   PD[-h$347316][pertS[LI[1], -p$347315]]*Derivative[0, 1][G3fun][scalar[], 
    X[]])/scale[]^2 + (4*perth[LI[1], -i, -j]*pprimescalar[]*primescalar[]*
   timevec[h$364063]*PD[-h$364063][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*hubbleC[]*perth[LI[1], -i, -j]*primescalar[]^2*timevec[h$364063]*
   PD[-h$364063][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*perth[LI[1], -i, -j]*primescalar[]^2*timevec[h$364063]*
   timevec[h$364064]*PD[-h$364064][PD[-h$364063][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (8*metric\[Delta][-i, -j]*pertphi[LI[1]]*pprimescalar[]*primescalar[]*
   timevec[h$364378]*PD[-h$364378][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (12*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^2*
   timevec[h$364378]*PD[-h$364378][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^2*timevec[h$364378]*
   timevec[h$364379]*PD[-h$364379][PD[-h$364378][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (16*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (24*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^2*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$365084]*
   timevec[h$365085]*PD[-h$365084][pertpsi[LI[1]]]*
   PD[-h$365085][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (8*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^2*
   timevec[h$365084]*timevec[h$365085]*PD[-h$365085][
    PD[-h$365084][pertscalar[LI[1]]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]*
   timevec[h$365348]*PD[-h$365348][pertscalar[LI[2]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (3*hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$365348]*
   PD[-h$365348][pertscalar[LI[2]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$365348]*
   timevec[h$365349]*PD[-h$365349][PD[-h$365348][pertscalar[LI[2]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*pprimescalar[]*timevec[h$366624]*timevec[h$366625]*
   PD[-h$366624][pertscalar[LI[1]]]*PD[-h$366625][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*hubbleC[]*metric\[Delta][-i, -j]*primescalar[]*timevec[h$366624]*
   timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*metric\[Delta][-i, -j]*primescalar[]*timevec[h$366624]*
   timevec[h$366625]*timevec[h$366628]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366628][PD[-h$366625][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (8*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^2*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (8*hubbleC[]*pertpsi[LI[1]]*primescalar[]^3*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*primescalar[]^3*timevec[h$349419]*PD[-h$349419][pertpsi[LI[1]]]*
   PD[-i][pertF[LI[1], -j]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*pprimescalar[]*primescalar[]*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*hubbleC[]*primescalar[]^2*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]^2*timevec[h$367275]*timevec[h$367276]*
   PD[-h$367276][PD[-h$367275][pertscalar[LI[1]]]]*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (pprimescalar[]*primescalar[]^2*PD[-i][pertF[LI[2], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*primescalar[]^3*PD[-i][pertF[LI[2], -j]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]*timevec[h$366624]*PD[-h$366624][PD[-j][pertscalar[LI[1]]]]*
   PD[-i][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (8*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^2*
   PD[-j][pertF[LI[1], -i]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (8*hubbleC[]*pertpsi[LI[1]]*primescalar[]^3*
   PD[-j][pertF[LI[1], -i]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*primescalar[]^3*timevec[h$349419]*
   PD[-h$349419][pertpsi[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*pprimescalar[]*primescalar[]*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*hubbleC[]*primescalar[]^2*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]^2*timevec[h$367275]*timevec[h$367276]*
   PD[-h$367276][PD[-h$367275][pertscalar[LI[1]]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (pprimescalar[]*primescalar[]^2*PD[-j][pertF[LI[2], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*primescalar[]^3*PD[-j][pertF[LI[2], -i]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*primescalar[]^2*PD[-i][pertscalar[LI[1]]]*PD[-j][pertpsi[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]*timevec[h$366624]*PD[-h$366624][PD[-i][pertscalar[LI[1]]]]*
   PD[-j][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*primescalar[]^2*PD[-i][pertpsi[LI[1]]]*
   PD[-j][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*pprimescalar[]*PD[-i][pertscalar[LI[1]]]*
   PD[-j][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*hubbleC[]*primescalar[]*PD[-i][pertscalar[LI[1]]]*
   PD[-j][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (16*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^2*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (16*hubbleC[]*pertpsi[LI[1]]*primescalar[]^3*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*primescalar[]^3*timevec[h$351798]*
   PD[-h$351798][pertpsi[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (8*pprimescalar[]*primescalar[]*timevec[h$367638]*
   PD[-h$367638][pertscalar[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (12*hubbleC[]*primescalar[]^2*timevec[h$367638]*
   PD[-h$367638][pertscalar[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*primescalar[]^2*timevec[h$367638]*timevec[h$367639]*
   PD[-h$367639][PD[-h$367638][pertscalar[LI[1]]]]*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*pprimescalar[]*primescalar[]^2*PD[-j][PD[-i][pertE[LI[2]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*hubbleC[]*primescalar[]^3*PD[-j][PD[-i][pertE[LI[2]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$350661, p$350662]*
   pertS[LI[1], -p$350661]*primescalar[]^3*PD[-p$350662][pertpsi[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$350990, p$350991]*primescalar[]^3*
   PD[-p$350990][pertB[LI[1]]]*PD[-p$350991][pertpsi[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*primescalar[]^3*
   timevec[h$357740]*PD[-h$357740][pertS[LI[1], -p$357739]]*
   PD[-p$357738][pertB[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (8*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
   pertS[LI[1], -p$357738]*pprimescalar[]*primescalar[]^2*
   PD[-p$357739][pertB[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (8*hubbleC[]*metric\[Delta][-i, -j]*
   metric\[Delta][p$357738, p$357739]*pertS[LI[1], -p$357738]*primescalar[]^3*
   PD[-p$357739][pertB[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
   pertS[LI[1], -p$357738]*primescalar[]^3*timevec[h$357740]*
   PD[-p$357739][PD[-h$357740][pertB[LI[1]]]]*Derivative[0, 1][G3fun][
    scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*pprimescalar[]*
   primescalar[]^2*PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
   primescalar[]^3*PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*primescalar[]^3*
   timevec[h$359828]*PD[-p$359826][pertB[LI[1]]]*
   PD[-p$359827][PD[-h$359828][pertB[LI[1]]]]*Derivative[0, 1][G3fun][
    scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$365086, p$365087]*primescalar[]^2*
   PD[-p$365086][pertpsi[LI[1]]]*PD[-p$365087][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$366626, p$366627]*
   primescalar[]*PD[-p$366626][pertscalar[LI[1]]]*
   PD[-p$366627][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*metric\[Delta][-i, -j]*metric\[Delta][p$366626, p$366627]*
   primescalar[]*timevec[h$366624]*PD[-p$366626][pertscalar[LI[1]]]*
   PD[-p$366627][PD[-h$366624][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
   pertS[LI[1], -p$367889]*pprimescalar[]*primescalar[]*
   PD[-p$367890][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (6*hubbleC[]*metric\[Delta][-i, -j]*
   metric\[Delta][p$367889, p$367890]*pertS[LI[1], -p$367889]*primescalar[]^2*
   PD[-p$367890][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
   primescalar[]^2*timevec[h$367891]*PD[-h$367891][pertS[LI[1], -p$367889]]*
   PD[-p$367890][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
   pertS[LI[1], -p$367889]*primescalar[]^2*timevec[h$367891]*
   PD[-p$367890][PD[-h$367891][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*pprimescalar[]*
   primescalar[]*PD[-p$368180][pertB[LI[1]]]*PD[-p$368181][pertscalar[LI[1]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (6*hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*
   primescalar[]^2*PD[-p$368180][pertB[LI[1]]]*
   PD[-p$368181][pertscalar[LI[1]]]*Derivative[0, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*
   primescalar[]^2*timevec[h$368182]*PD[-p$368180][pertscalar[LI[1]]]*
   PD[-p$368181][PD[-h$368182][pertB[LI[1]]]]*Derivative[0, 1][G3fun][
    scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*primescalar[]^2*
   timevec[h$368182]*PD[-p$368180][pertB[LI[1]]]*
   PD[-p$368181][PD[-h$368182][pertscalar[LI[1]]]]*
   Derivative[0, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^4*
   Derivative[0, 2][G2fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^3*timevec[h$365084]*
   PD[-h$365084][pertscalar[LI[1]]]*Derivative[0, 2][G2fun][scalar[], X[]])/
  scale[]^2 - (metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$366624]*
   timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 2][G2fun][scalar[], X[]])/
  scale[]^2 - (2*perth[LI[1], -i, -j]*pertpsi[LI[1]]*pprimescalar[]*
   primescalar[]^4*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*pprimescalar[]*
   primescalar[]^4*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (12*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*pprimescalar[]*primescalar[]^4*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (metric\[Delta][-i, -j]*pertpsi[LI[2]]*pprimescalar[]*primescalar[]^4*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315]*pprimescalar[]*
   primescalar[]^4*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*hubbleC[]*perth[LI[1], -i, -j]*pertpsi[LI[1]]*primescalar[]^5*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (4*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*
   primescalar[]^5*Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (12*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^5*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[2]]*primescalar[]^5*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315]*primescalar[]^5*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^5*timevec[h$346451]*
   PD[-h$346451][pertpsi[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (2*perth[LI[1], -i, -j]*pprimescalar[]*primescalar[]^3*
   timevec[h$364063]*PD[-h$364063][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*hubbleC[]*perth[LI[1], -i, -j]*primescalar[]^4*timevec[h$364063]*
   PD[-h$364063][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pprimescalar[]*
   primescalar[]^3*timevec[h$364378]*PD[-h$364378][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^4*
   timevec[h$364378]*PD[-h$364378][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (16*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^3*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (18*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^4*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*metric\[Delta][-i, -j]*primescalar[]^4*timevec[h$365084]*
   timevec[h$365085]*PD[-h$365084][pertpsi[LI[1]]]*
   PD[-h$365085][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^4*
   timevec[h$365084]*timevec[h$365085]*PD[-h$365085][
    PD[-h$365084][pertscalar[LI[1]]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]^3*
   timevec[h$365348]*PD[-h$365348][pertscalar[LI[2]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^4*timevec[h$365348]*
   PD[-h$365348][pertscalar[LI[2]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (5*metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]^2*
   timevec[h$366624]*timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (7*hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^3*
   timevec[h$366624]*timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (2*metric\[Delta][-i, -j]*primescalar[]^3*timevec[h$366624]*
   timevec[h$366625]*timevec[h$366628]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366628][PD[-h$366625][pertscalar[LI[1]]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^4*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*hubbleC[]*pertpsi[LI[1]]*primescalar[]^5*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*pprimescalar[]*primescalar[]^3*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*hubbleC[]*primescalar[]^4*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-i][pertF[LI[1], -j]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^4*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*hubbleC[]*pertpsi[LI[1]]*primescalar[]^5*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*pprimescalar[]*primescalar[]^3*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*hubbleC[]*primescalar[]^4*timevec[h$367275]*
   PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (4*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^4*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (4*hubbleC[]*pertpsi[LI[1]]*primescalar[]^5*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (4*pprimescalar[]*primescalar[]^3*timevec[h$367638]*
   PD[-h$367638][pertscalar[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (4*hubbleC[]*primescalar[]^4*timevec[h$367638]*
   PD[-h$367638][pertscalar[LI[1]]]*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
   pertS[LI[1], -p$357738]*pprimescalar[]*primescalar[]^4*
   PD[-p$357739][pertB[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (2*hubbleC[]*metric\[Delta][-i, -j]*
   metric\[Delta][p$357738, p$357739]*pertS[LI[1], -p$357738]*primescalar[]^5*
   PD[-p$357739][pertB[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*
   pprimescalar[]*primescalar[]^4*PD[-p$359826][pertB[LI[1]]]*
   PD[-p$359827][pertB[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$359826, 
    p$359827]*primescalar[]^5*PD[-p$359826][pertB[LI[1]]]*
   PD[-p$359827][pertB[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (metric\[Delta][-i, -j]*metric\[Delta][p$366626, p$366627]*
   pprimescalar[]*primescalar[]^2*PD[-p$366626][pertscalar[LI[1]]]*
   PD[-p$366627][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (hubbleC[]*metric\[Delta][-i, -j]*metric\[Delta][p$366626, 
    p$366627]*primescalar[]^3*PD[-p$366626][pertscalar[LI[1]]]*
   PD[-p$366627][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (2*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
   pertS[LI[1], -p$367889]*pprimescalar[]*primescalar[]^3*
   PD[-p$367890][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (2*hubbleC[]*metric\[Delta][-i, -j]*
   metric\[Delta][p$367889, p$367890]*pertS[LI[1], -p$367889]*primescalar[]^4*
   PD[-p$367890][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 - (2*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*
   pprimescalar[]*primescalar[]^3*PD[-p$368180][pertB[LI[1]]]*
   PD[-p$368181][pertscalar[LI[1]]]*Derivative[0, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (2*hubbleC[]*metric\[Delta][-i, -j]*
   metric\[Delta][p$368180, p$368181]*primescalar[]^4*
   PD[-p$368180][pertB[LI[1]]]*PD[-p$368181][pertscalar[LI[1]]]*
   Derivative[0, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*pprimescalar[]*primescalar[]^6*
   Derivative[0, 3][G3fun][scalar[], X[]])/scale[]^6 - 
 (hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^7*
   Derivative[0, 3][G3fun][scalar[], X[]])/scale[]^6 - 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pprimescalar[]*primescalar[]^5*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 3][G3fun][scalar[], X[]])/scale[]^6 + 
 (2*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^6*
   timevec[h$365084]*PD[-h$365084][pertscalar[LI[1]]]*
   Derivative[0, 3][G3fun][scalar[], X[]])/scale[]^6 + 
 (metric\[Delta][-i, -j]*pprimescalar[]*primescalar[]^4*timevec[h$366624]*
   timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 3][G3fun][scalar[], X[]])/
  scale[]^6 - (hubbleC[]*metric\[Delta][-i, -j]*primescalar[]^5*
   timevec[h$366624]*timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[0, 3][G3fun][scalar[], X[]])/
  scale[]^6 - 2*perth[LI[1], -i, -j]*pertscalar[LI[1]]*scale[]^2*
  Derivative[1, 0][G2fun][scalar[], X[]] + 4*metric\[Delta][-i, -j]*
  pertphi[LI[1]]*pertscalar[LI[1]]*scale[]^2*Derivative[1, 0][G2fun][
   scalar[], X[]] - metric\[Delta][-i, -j]*pertscalar[LI[2]]*scale[]^2*
  Derivative[1, 0][G2fun][scalar[], X[]] - 2*pertscalar[LI[1]]*scale[]^2*
  PD[-i][pertF[LI[1], -j]]*Derivative[1, 0][G2fun][scalar[], X[]] - 
 2*pertscalar[LI[1]]*scale[]^2*PD[-j][pertF[LI[1], -i]]*
  Derivative[1, 0][G2fun][scalar[], X[]] - 4*pertscalar[LI[1]]*scale[]^2*
  PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 0][G2fun][scalar[], X[]] + 
 perth[LI[2], -i, -j]*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], 
   X[]] - 2*metric\[Delta][-i, -j]*pertphi[LI[2]]*primescalar[]^2*
  Derivative[1, 0][G3fun][scalar[], X[]] - 4*perth[LI[1], -i, -j]*
  pertpsi[LI[1]]*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]] + 
 8*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*primescalar[]^2*
  Derivative[1, 0][G3fun][scalar[], X[]] + 8*metric\[Delta][-i, -j]*
  pertpsi[LI[1]]^2*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], X[]] - 
 2*metric\[Delta][-i, -j]*pertpsi[LI[2]]*primescalar[]^2*
  Derivative[1, 0][G3fun][scalar[], X[]] - 2*metric\[Delta][-i, -j]*
  metric\[Delta][p$347314, p$347315]*pertS[LI[1], -p$347314]*
  pertS[LI[1], -p$347315]*primescalar[]^2*Derivative[1, 0][G3fun][scalar[], 
   X[]] + 4*perth[LI[1], -i, -j]*primescalar[]*timevec[h$364063]*
  PD[-h$364063][pertscalar[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] - 
 8*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]*timevec[h$364378]*
  PD[-h$364378][pertscalar[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] - 
 8*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]*timevec[h$365084]*
  PD[-h$365084][pertscalar[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 2*metric\[Delta][-i, -j]*primescalar[]*timevec[h$365348]*
  PD[-h$365348][pertscalar[LI[2]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 2*metric\[Delta][-i, -j]*timevec[h$366624]*timevec[h$366625]*
  PD[-h$366624][pertscalar[LI[1]]]*PD[-h$366625][pertscalar[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] - 4*pertpsi[LI[1]]*primescalar[]^2*
  PD[-i][pertF[LI[1], -j]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 4*primescalar[]*timevec[h$367275]*PD[-h$367275][pertscalar[LI[1]]]*
  PD[-i][pertF[LI[1], -j]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 primescalar[]^2*PD[-i][pertF[LI[2], -j]]*Derivative[1, 0][G3fun][scalar[], 
   X[]] - 4*pertpsi[LI[1]]*primescalar[]^2*PD[-j][pertF[LI[1], -i]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 4*primescalar[]*timevec[h$367275]*
  PD[-h$367275][pertscalar[LI[1]]]*PD[-j][pertF[LI[1], -i]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 
 primescalar[]^2*PD[-j][pertF[LI[2], -i]]*Derivative[1, 0][G3fun][scalar[], 
   X[]] + 4*PD[-i][pertscalar[LI[1]]]*PD[-j][pertscalar[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] - 8*pertpsi[LI[1]]*primescalar[]^2*
  PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 8*primescalar[]*timevec[h$367638]*PD[-h$367638][pertscalar[LI[1]]]*
  PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 2*primescalar[]^2*PD[-j][PD[-i][pertE[LI[2]]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 4*metric\[Delta][-i, -j]*
  metric\[Delta][p$357738, p$357739]*pertS[LI[1], -p$357738]*primescalar[]^2*
  PD[-p$357739][pertB[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] - 
 2*metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*primescalar[]^2*
  PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] - 2*metric\[Delta][-i, -j]*
  metric\[Delta][p$366626, p$366627]*PD[-p$366626][pertscalar[LI[1]]]*
  PD[-p$366627][pertscalar[LI[1]]]*Derivative[1, 0][G3fun][scalar[], X[]] + 
 4*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
  pertS[LI[1], -p$367889]*primescalar[]*PD[-p$367890][pertscalar[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] - 4*metric\[Delta][-i, -j]*
  metric\[Delta][p$368180, p$368181]*primescalar[]*
  PD[-p$368180][pertB[LI[1]]]*PD[-p$368181][pertscalar[LI[1]]]*
  Derivative[1, 0][G3fun][scalar[], X[]] + 2*metric\[Delta][-i, -j]*
  pertpsi[LI[1]]*pertscalar[LI[1]]*primescalar[]^2*
  Derivative[1, 1][G2fun][scalar[], X[]] - 2*metric\[Delta][-i, -j]*
  pertscalar[LI[1]]*primescalar[]*timevec[h$366624]*
  PD[-h$366624][pertscalar[LI[1]]]*Derivative[1, 1][G2fun][scalar[], X[]] + 
 (2*perth[LI[1], -i, -j]*pertscalar[LI[1]]*pprimescalar[]*primescalar[]^2*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertscalar[LI[1]]*pprimescalar[]*
   primescalar[]^2*Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (8*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pertscalar[LI[1]]*pprimescalar[]*
   primescalar[]^2*Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (metric\[Delta][-i, -j]*pertscalar[LI[2]]*pprimescalar[]*primescalar[]^2*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*hubbleC[]*perth[LI[1], -i, -j]*pertscalar[LI[1]]*primescalar[]^3*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*hubbleC[]*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertscalar[LI[1]]*
   primescalar[]^3*Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (8*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pertscalar[LI[1]]*
   primescalar[]^3*Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][-i, -j]*pertscalar[LI[2]]*primescalar[]^3*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*perth[LI[1], -i, -j]*pertpsi[LI[1]]*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertpsi[LI[1]]*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (8*metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*pertpsi[LI[2]]*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$347314, p$347315]*
   pertS[LI[1], -p$347314]*pertS[LI[1], -p$347315]*primescalar[]^4*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*perth[LI[1], -i, -j]*primescalar[]^3*timevec[h$364063]*
   PD[-h$364063][pertscalar[LI[1]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*metric\[Delta][-i, -j]*pertphi[LI[1]]*primescalar[]^3*
   timevec[h$364378]*PD[-h$364378][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^3*
   timevec[h$365084]*PD[-h$365084][pertpsi[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (12*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^3*timevec[h$365084]*
   PD[-h$365084][pertscalar[LI[1]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (metric\[Delta][-i, -j]*primescalar[]^3*timevec[h$365348]*
   PD[-h$365348][pertscalar[LI[2]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*metric\[Delta][-i, -j]*pertscalar[LI[1]]*pprimescalar[]*
   primescalar[]*timevec[h$366624]*PD[-h$366624][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (6*hubbleC[]*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^2*
   timevec[h$366624]*PD[-h$366624][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (5*metric\[Delta][-i, -j]*primescalar[]^2*timevec[h$366624]*
   timevec[h$366625]*PD[-h$366624][pertscalar[LI[1]]]*
   PD[-h$366625][pertscalar[LI[1]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^2*
   timevec[h$366624]*timevec[h$366625]*PD[-h$366625][
    PD[-h$366624][pertscalar[LI[1]]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*pertscalar[LI[1]]*pprimescalar[]*primescalar[]^2*
   PD[-i][pertF[LI[1], -j]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*hubbleC[]*pertscalar[LI[1]]*primescalar[]^3*
   PD[-i][pertF[LI[1], -j]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*pertpsi[LI[1]]*primescalar[]^4*PD[-i][pertF[LI[1], -j]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]^3*timevec[h$367275]*PD[-h$367275][pertscalar[LI[1]]]*
   PD[-i][pertF[LI[1], -j]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*pertscalar[LI[1]]*pprimescalar[]*primescalar[]^2*
   PD[-j][pertF[LI[1], -i]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*hubbleC[]*pertscalar[LI[1]]*primescalar[]^3*
   PD[-j][pertF[LI[1], -i]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (2*pertpsi[LI[1]]*primescalar[]^4*PD[-j][pertF[LI[1], -i]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*primescalar[]^3*timevec[h$367275]*PD[-h$367275][pertscalar[LI[1]]]*
   PD[-j][pertF[LI[1], -i]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (4*pertscalar[LI[1]]*pprimescalar[]*primescalar[]^2*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*hubbleC[]*pertscalar[LI[1]]*primescalar[]^3*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 - (4*pertpsi[LI[1]]*primescalar[]^4*PD[-j][PD[-i][pertE[LI[1]]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (4*primescalar[]^3*timevec[h$367638]*PD[-h$367638][pertscalar[LI[1]]]*
   PD[-j][PD[-i][pertE[LI[1]]]]*Derivative[1, 1][G3fun][scalar[], X[]])/
  scale[]^2 + (2*metric\[Delta][-i, -j]*metric\[Delta][p$357738, p$357739]*
   pertS[LI[1], -p$357738]*primescalar[]^4*PD[-p$357739][pertB[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$359826, p$359827]*primescalar[]^4*
   PD[-p$359826][pertB[LI[1]]]*PD[-p$359827][pertB[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (metric\[Delta][-i, -j]*metric\[Delta][p$366626, p$366627]*primescalar[]^2*
   PD[-p$366626][pertscalar[LI[1]]]*PD[-p$366627][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$367889, p$367890]*
   pertS[LI[1], -p$367889]*primescalar[]^3*PD[-p$367890][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*metric\[Delta][p$368180, p$368181]*primescalar[]^3*
   PD[-p$368180][pertB[LI[1]]]*PD[-p$368181][pertscalar[LI[1]]]*
   Derivative[1, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pertscalar[LI[1]]*pprimescalar[]*
   primescalar[]^4*Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (2*hubbleC[]*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pertscalar[LI[1]]*
   primescalar[]^5*Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (metric\[Delta][-i, -j]*pertpsi[LI[1]]^2*primescalar[]^6*
   Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*primescalar[]^5*timevec[h$365084]*
   PD[-h$365084][pertscalar[LI[1]]]*Derivative[1, 2][G3fun][scalar[], X[]])/
  scale[]^4 + (2*metric\[Delta][-i, -j]*pertscalar[LI[1]]*pprimescalar[]*
   primescalar[]^3*timevec[h$366624]*PD[-h$366624][pertscalar[LI[1]]]*
   Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 (2*hubbleC[]*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^4*
   timevec[h$366624]*PD[-h$366624][pertscalar[LI[1]]]*
   Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 + 
 (metric\[Delta][-i, -j]*primescalar[]^4*timevec[h$366624]*timevec[h$366625]*
   PD[-h$366624][pertscalar[LI[1]]]*PD[-h$366625][pertscalar[LI[1]]]*
   Derivative[1, 2][G3fun][scalar[], X[]])/scale[]^4 - 
 metric\[Delta][-i, -j]*pertscalar[LI[1]]^2*scale[]^2*
  Derivative[2, 0][G2fun][scalar[], X[]] + 2*perth[LI[1], -i, -j]*
  pertscalar[LI[1]]*primescalar[]^2*Derivative[2, 0][G3fun][scalar[], X[]] - 
 4*metric\[Delta][-i, -j]*pertphi[LI[1]]*pertscalar[LI[1]]*primescalar[]^2*
  Derivative[2, 0][G3fun][scalar[], X[]] - 4*metric\[Delta][-i, -j]*
  pertpsi[LI[1]]*pertscalar[LI[1]]*primescalar[]^2*
  Derivative[2, 0][G3fun][scalar[], X[]] + metric\[Delta][-i, -j]*
  pertscalar[LI[2]]*primescalar[]^2*Derivative[2, 0][G3fun][scalar[], X[]] + 
 4*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]*timevec[h$366624]*
  PD[-h$366624][pertscalar[LI[1]]]*Derivative[2, 0][G3fun][scalar[], X[]] + 
 2*pertscalar[LI[1]]*primescalar[]^2*PD[-i][pertF[LI[1], -j]]*
  Derivative[2, 0][G3fun][scalar[], X[]] + 
 2*pertscalar[LI[1]]*primescalar[]^2*PD[-j][pertF[LI[1], -i]]*
  Derivative[2, 0][G3fun][scalar[], X[]] + 
 4*pertscalar[LI[1]]*primescalar[]^2*PD[-j][PD[-i][pertE[LI[1]]]]*
  Derivative[2, 0][G3fun][scalar[], X[]] + 
 (metric\[Delta][-i, -j]*pertscalar[LI[1]]^2*pprimescalar[]*primescalar[]^2*
   Derivative[2, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (hubbleC[]*metric\[Delta][-i, -j]*pertscalar[LI[1]]^2*primescalar[]^3*
   Derivative[2, 1][G3fun][scalar[], X[]])/scale[]^2 - 
 (2*metric\[Delta][-i, -j]*pertpsi[LI[1]]*pertscalar[LI[1]]*primescalar[]^4*
   Derivative[2, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 (2*metric\[Delta][-i, -j]*pertscalar[LI[1]]*primescalar[]^3*
   timevec[h$366624]*PD[-h$366624][pertscalar[LI[1]]]*
   Derivative[2, 1][G3fun][scalar[], X[]])/scale[]^2 + 
 metric\[Delta][-i, -j]*pertscalar[LI[1]]^2*primescalar[]^2*
  Derivative[3, 0][G3fun][scalar[], X[]]
