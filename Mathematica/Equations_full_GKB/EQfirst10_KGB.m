(* Created with the Wolfram Language : www.wolfram.com *)
density[]*hubbleC[]*metric\[Delta][p$2790518, p$2790519]*
  PD[-p$2790519][PD[-p$2790518][pertS[LI[1], -i]]] + 
 hubbleC[]*metric\[Delta][p$2790518, p$2790519]*pressure[]*
  PD[-p$2790519][PD[-p$2790518][pertS[LI[1], -i]]] + 
 metric\[Delta][p$2790518, p$2790519]*primepressure[]*
  PD[-p$2790519][PD[-p$2790518][pertS[LI[1], -i]]] + 
 density[]*metric\[Delta][p$2790518, p$2790519]*timevec[h$2790520]*
  PD[-p$2790519][PD[-p$2790518][PD[-h$2790520][pertS[LI[1], -i]]]] + 
 metric\[Delta][p$2790518, p$2790519]*pressure[]*timevec[h$2790520]*
  PD[-p$2790519][PD[-p$2790518][PD[-h$2790520][pertS[LI[1], -i]]]] - 
 (metric\[Delta][-i, -p$2790518]*metric\[Delta][p$2790519, p$2790521]*
   primepressure[]*PD[-p$2790521][PD[-p$2790519][pertvelocityvec[LI[1], 
      p$2790518]]])/scale[] - (density[]*metric\[Delta][-i, -p$2790518]*
   metric\[Delta][p$2790519, p$2790521]*timevec[h$2790520]*
   PD[-p$2790521][PD[-p$2790519][PD[-h$2790520][pertvelocityvec[LI[1], 
       p$2790518]]]])/scale[] - (metric\[Delta][-i, -p$2790518]*
   metric\[Delta][p$2790519, p$2790521]*pressure[]*timevec[h$2790520]*
   PD[-p$2790521][PD[-p$2790519][PD[-h$2790520][pertvelocityvec[LI[1], 
       p$2790518]]]])/scale[] - 
 (density[]*metric\[Delta][p$2790518, p$2790519]*metric\[Delta][p$2790521, 
    p$2790522]*PD[-p$2790522][PD[-p$2790521][PD[-p$2790519][
      PD[-p$2790518][pertshearvec[LI[1], -i]]]]])/scale[]^2 - 
 (metric\[Delta][p$2790518, p$2790519]*metric\[Delta][p$2790521, p$2790522]*
   pressure[]*PD[-p$2790522][PD[-p$2790521][PD[-p$2790519][
      PD[-p$2790518][pertshearvec[LI[1], -i]]]]])/scale[]^2
