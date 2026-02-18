ClearAll["Global`*"];

nb = Notebook[{
  Cell["05 - Nonlocal Quadratic Mixing Kernel", "Title"],
  Cell["Two-field quadratic kernel induced by second variation of Q=I1/I0: explicit mixing coefficients and IR scaling.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dI0 == a0 x + b0 y && dI1 == a1 x + b1 y
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    d2I0 == aa0 x^2 + 2 ab0 x y + bb0 y^2 && d2I1 == aa1 x^2 + 2 ab1 x y + bb1 y^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2 == (d2I1 - q d2I0)/i0 - 2 dI0 (dI1 - q dI0)/i0^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2 == Mxx x^2 + Mxy x y + Myy y^2
  ], "Input"]
}];

dI0 = a0 x + b0 y;
dI1 = a1 x + b1 y;
d2I0 = aa0 x^2 + 2 ab0 x y + bb0 y^2;
d2I1 = aa1 x^2 + 2 ab1 x y + bb1 y^2;

q2 = FullSimplify[(d2I1 - q d2I0)/i0 - 2 dI0 (dI1 - q dI0)/i0^2];

Mxx = FullSimplify[Coefficient[q2, x^2]];
Mxy = FullSimplify[Coefficient[q2, x y]];
Myy = FullSimplify[Coefficient[q2, y^2]];

A = a1 - q a0;
B = b1 - q b0;

MxxExpected = FullSimplify[(aa1 - q aa0)/i0 - 2 a0 A/i0^2];
MxyExpected = FullSimplify[2 (ab1 - q ab0)/i0 - 2 (a0 B + b0 A)/i0^2];
MyyExpected = FullSimplify[(bb1 - q bb0)/i0 - 2 b0 B/i0^2];

checkA = FullSimplify[Mxx == MxxExpected];
checkB = FullSimplify[Mxy == MxyExpected];
checkC = FullSimplify[Myy == MyyExpected];

Kmix = {
  {Mxx, Mxy/2},
  {Mxy/2, Myy}
};

Klead = FullSimplify[Limit[i0 Kmix, i0 -> Infinity]];
checkD = FullSimplify[Klead == {{aa1 - q aa0, ab1 - q ab0}, {ab1 - q ab0, bb1 - q bb0}}];

numRule = {
  q -> 3/5,
  a0 -> 1/20, b0 -> 1/20,
  a1 -> 2/25, b1 -> 9/100,
  aa0 -> 1/10, ab0 -> 1/10, bb0 -> 1/10,
  aa1 -> 3/5, ab1 -> 1/5, bb1 -> 7/10
};

eigVals = Table[
  N[
    Eigenvalues[Kmix /. numRule /. i0 -> LL^3],
    30
  ],
  {LL, {10, 20, 40}}
];

eigPos = And @@ Flatten[Map[# > 0 &, eigVals, {2}]];
eigDecreasing = Max[Abs[eigVals[[3]]]] < Max[Abs[eigVals[[1]]]];

scaledKVals = Table[
  N[(LL^3 Kmix) /. numRule /. i0 -> LL^3, 30],
  {LL, {10, 20, 40}}
];

leadNum = N[Klead /. numRule, 30];
errLead = N[Max[Abs[Flatten[scaledKVals[[-1]] - leadNum]]], 30];
checkNum = eigPos && eigDecreasing && errLead < 1*^-3;

check = TrueQ[checkA && checkB && checkC && checkD && checkNum];

nbName = "05_Nonlocal_Quadratic_Mixing_Kernel.nb";
logName = "05_nonlocal_quadratic_mixing_kernel.log";
logLines = {
  "Notebook: " <> nbName,
  "q2 = " <> ToString[q2, InputForm],
  "Mxx = " <> ToString[Mxx, InputForm],
  "MxxExpected = " <> ToString[MxxExpected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "Mxy = " <> ToString[Mxy, InputForm],
  "MxyExpected = " <> ToString[MxyExpected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "Myy = " <> ToString[Myy, InputForm],
  "MyyExpected = " <> ToString[MyyExpected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "Klead = " <> ToString[Klead, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "eigVals(L={10,20,40}) = " <> ToString[eigVals, InputForm],
  "scaledKVals = " <> ToString[scaledKVals, InputForm],
  "leadNum = " <> ToString[leadNum, InputForm],
  "errLead = " <> ToString[errLead, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "05_nonlocal_quadratic_mixing_kernel",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
