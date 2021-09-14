(* ::Package:: *)

L2 = G2fun[scalarcov[], Xcov[]];
L3 = -G3fun[scalarcov[], Xcov[]] CD[-\[Mu]]@CD[\[Mu]]@scalarcov[];
L4 = Mpl^2/2 RicciScalarCD[];


Lgravity = Sqrt[-Detmetricg[]] ContractMetric[ToCanonical[EinsteinToRicci[L2 + L3 + L4]]] // NoScalar // ReplaceDummies;
