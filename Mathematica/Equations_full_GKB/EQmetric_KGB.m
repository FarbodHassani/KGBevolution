(* Created with the Wolfram Language : www.wolfram.com *)
-(G2fun[scalarcov[], Xcov[]]*metricg[-\[Mu], -\[Nu]]) + 
 Mpl^2*RicciCD[-\[Mu], -\[Nu]] - 
 (Mpl^2*metricg[-\[Mu], -\[Nu]]*RicciScalarCD[])/2 - 
 stressenergy[-\[Mu], -\[Nu]] - CD[-\[Mu]][scalarcov[]]*
  CD[-\[Nu]][scalarcov[]]*Derivative[0, 1][G2fun][scalarcov[], Xcov[]] + 
 metricg[-\[Mu], -\[Nu]]*CD[\[Gamma]$22053][scalarcov[]]*
  CD[-\[Gamma]$30563][CD[-\[Gamma]$22053][scalarcov[]]]*
  CD[\[Gamma]$30563][scalarcov[]]*Derivative[0, 1][G3fun][scalarcov[], 
   Xcov[]] + CD[-\[Gamma]$22053][CD[\[Gamma]$22053][scalarcov[]]]*
  CD[-\[Mu]][scalarcov[]]*CD[-\[Nu]][scalarcov[]]*
  Derivative[0, 1][G3fun][scalarcov[], Xcov[]] - 
 CD[\[Gamma]$22053][scalarcov[]]*CD[-\[Mu]][CD[-\[Gamma]$22053][scalarcov[]]]*
  CD[-\[Nu]][scalarcov[]]*Derivative[0, 1][G3fun][scalarcov[], Xcov[]] - 
 CD[\[Gamma]$22053][scalarcov[]]*CD[-\[Mu]][scalarcov[]]*
  CD[-\[Nu]][CD[-\[Gamma]$22053][scalarcov[]]]*
  Derivative[0, 1][G3fun][scalarcov[], Xcov[]] - 
 metricg[-\[Mu], -\[Nu]]*CD[-\[Gamma]$22053][scalarcov[]]*
  CD[\[Gamma]$22053][scalarcov[]]*Derivative[1, 0][G3fun][scalarcov[], 
   Xcov[]] + 2*CD[-\[Mu]][scalarcov[]]*CD[-\[Nu]][scalarcov[]]*
  Derivative[1, 0][G3fun][scalarcov[], Xcov[]]
