ClearAll["Global`*"];

nb = Notebook[{
  Cell["03 - FRW Scalar Constraint Freeze Map", "Title"],
  Cell["Map from the linear FRW constraint channel (already derived in Variational_Derivation/17) to explicit scalar freezing relations.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    C1onBranch == -2 Exp[3 H0 t] (K0^2 n1 - 3 H0 up)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    H0 == K0/Sqrt[3] && n1 == (Sqrt[3]/K0) up
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    n1 == 0 -> up == 0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lkin == a1 up^2 + a2 up n1 + a3 n1^2 -> LkinRed[up]
  ], "Input"]
}];

t = Symbol["t"];
C1onBranch = -2 Exp[3 H0 t] (K0^2 n1 - 3 H0 up);
ruleH = H0 -> K0/Sqrt[3];

nRule = First[Solve[(C1onBranch /. ruleH) == 0, n1]];
checkA = FullSimplify[nRule[[1, 2]] == (Sqrt[3]/K0) up];

upGauge = FullSimplify[(up /. nRule) /. n1 -> 0];
upRuleGauge = First[Solve[(n1 /. nRule) == 0, up]];
upGauge = FullSimplify[up /. upRuleGauge];
checkB = FullSimplify[upGauge == 0];

Lkin = a1 up^2 + a2 up n1 + a3 n1^2;
Hess = {
  {D[Lkin, {up, 2}], D[D[Lkin, up], n1]},
  {D[D[Lkin, n1], up], D[Lkin, {n1, 2}]}
};
detHess = FullSimplify[Det[Hess]];

LkinRed = FullSimplify[Expand[Lkin /. nRule]];
Kred = FullSimplify[Coefficient[LkinRed, up^2]];
checkC = FreeQ[LkinRed, n1];

LkinGauge = FullSimplify[LkinRed /. upRuleGauge];
checkD = FullSimplify[LkinGauge == 0];

numKred = N[
  Kred /. {a1 -> 1, a2 -> 1/2, a3 -> 1/4, K0 -> 3/5},
  30
];
checkNum = NumericQ[numKred];

check = TrueQ[checkA && checkB && checkC && checkD && checkNum];

nbName = "03_FRW_Scalar_Constraint_Freeze_Map.nb";
logName = "03_frw_scalar_constraint_freeze_map.log";
logLines = {
  "Notebook: " <> nbName,
  "C1onBranch = " <> ToString[C1onBranch, InputForm],
  "nRule = " <> ToString[nRule, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "upRuleGauge = " <> ToString[upRuleGauge, InputForm],
  "upGauge = " <> ToString[upGauge, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "Hess = " <> ToString[Hess, InputForm],
  "detHess = " <> ToString[detHess, InputForm],
  "LkinRed = " <> ToString[LkinRed, InputForm],
  "Kred = " <> ToString[Kred, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "LkinGauge = " <> ToString[LkinGauge, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "numKred = " <> ToString[numKred, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "03_frw_scalar_constraint_freeze_map",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
