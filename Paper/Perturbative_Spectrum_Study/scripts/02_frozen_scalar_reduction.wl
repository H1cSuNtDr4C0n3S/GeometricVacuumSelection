ClearAll["Global`*"];

nb = Notebook[{
  Cell["02 - Frozen Scalar Reduction", "Title"],
  Cell["Integrate out a non-dynamical scalar in a mixed quadratic action and extract effective kinetic/gradient coefficients for the remaining scalar mode.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L2 == 1/2 A zd^2 + B zd pi - 1/2 (C + F q2) pi^2 - E q2 z pi - 1/2 D q2 z^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    piRule == Solve[D[L2, pi] == 0, pi][[1]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Leff == L2 /. piRule
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Keff == A + B^2/(C + F q2) && Geff == D - (E^2 q2)/(C + F q2)
  ], "Input"]
}];

L2 = 1/2 A zd^2 + B zd pi - 1/2 (C + F q2) pi^2 - E q2 z pi - 1/2 D q2 z^2;
piRule = First[Solve[D[L2, pi] == 0, pi]];
Leff = FullSimplify[Expand[L2 /. piRule], Assumptions -> C + F q2 != 0];

Keff = FullSimplify[2 Coefficient[Leff, zd^2], Assumptions -> C + F q2 != 0];
crossEff = FullSimplify[Coefficient[Leff, z zd], Assumptions -> C + F q2 != 0];
Geff = FullSimplify[-2 Coefficient[Leff, z^2]/q2, Assumptions -> C + F q2 != 0];

KeffExpected = A + B^2/(C + F q2);
crossExpected = -(B E q2)/(C + F q2);
GeffExpected = D - (E^2 q2)/(C + F q2);

checkA = FullSimplify[Keff == KeffExpected, Assumptions -> C + F q2 != 0];
checkB = FullSimplify[crossEff == crossExpected, Assumptions -> C + F q2 != 0];
checkC = FullSimplify[Geff == GeffExpected, Assumptions -> C + F q2 != 0];

KeffHigh = FullSimplify[Limit[Keff, q2 -> Infinity], Assumptions -> F > 0];
GeffHigh = FullSimplify[Limit[Geff, q2 -> Infinity], Assumptions -> F > 0];
crossHigh = FullSimplify[Limit[crossEff, q2 -> Infinity], Assumptions -> F > 0];

checkD = FullSimplify[KeffHigh == A, Assumptions -> F > 0];
checkE = FullSimplify[GeffHigh == D - E^2/F, Assumptions -> F > 0];
checkF = FullSimplify[crossHigh == -(B E)/F, Assumptions -> F > 0];

KeffFreeze = FullSimplify[Limit[Keff, C -> Infinity], Assumptions -> F > 0 && q2 > 0];
GeffFreeze = FullSimplify[Limit[Geff, C -> Infinity], Assumptions -> F > 0 && q2 > 0];
crossFreeze = FullSimplify[Limit[crossEff, C -> Infinity], Assumptions -> F > 0 && q2 > 0];

checkG = FullSimplify[KeffFreeze == A];
checkH = FullSimplify[GeffFreeze == D];
checkI = FullSimplify[crossFreeze == 0];

stableVals = Table[
  N[
    {Keff, Geff} /. {
      A -> 2, B -> 1/5, C -> 4, D -> 13/10, E -> 1/2, F -> 11/10,
      q2 -> qq
    },
    30
  ],
  {qq, {1, 10, 100}}
];

checkNumA = Min[stableVals[[All, 1]]] > 0;
checkNumB = Min[stableVals[[All, 2]]] > 0;

unstableHigh = N[
  GeffHigh /. {D -> 1, E -> 2, F -> 1},
  30
];
checkNumC = unstableHigh < 0;

check = TrueQ[
  checkA && checkB && checkC &&
  checkD && checkE && checkF &&
  checkG && checkH && checkI &&
  checkNumA && checkNumB && checkNumC
];

nbName = "02_Frozen_Scalar_Reduction.nb";
logName = "02_frozen_scalar_reduction.log";
logLines = {
  "Notebook: " <> nbName,
  "piRule = " <> ToString[piRule, InputForm],
  "Leff = " <> ToString[Leff, InputForm],
  "Keff = " <> ToString[Keff, InputForm],
  "crossEff = " <> ToString[crossEff, InputForm],
  "Geff = " <> ToString[Geff, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "KeffHigh = " <> ToString[KeffHigh, InputForm],
  "GeffHigh = " <> ToString[GeffHigh, InputForm],
  "crossHigh = " <> ToString[crossHigh, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "checkE = " <> ToString[checkE, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "KeffFreeze = " <> ToString[KeffFreeze, InputForm],
  "GeffFreeze = " <> ToString[GeffFreeze, InputForm],
  "crossFreeze = " <> ToString[crossFreeze, InputForm],
  "checkG = " <> ToString[checkG, InputForm],
  "checkH = " <> ToString[checkH, InputForm],
  "checkI = " <> ToString[checkI, InputForm],
  "stableVals(q2={1,10,100}) = " <> ToString[stableVals, InputForm],
  "checkNumA(Keff>0) = " <> ToString[checkNumA, InputForm],
  "checkNumB(Geff>0) = " <> ToString[checkNumB, InputForm],
  "unstableHigh = " <> ToString[unstableHigh, InputForm],
  "checkNumC(unstable sample) = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "02_frozen_scalar_reduction",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];

