ClearAll["Global`*"];

nb = Notebook[{
  Cell["03 - Theta Chain Rule Variation", "Title"],
  Cell["Chain-rule structure for the Theta variation of the constrained leafwise integrand.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    f[th_] := lambda[th] (I1[th] - K0^2 I0[th])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dfdth = D[f[th], th] // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    expected = lambda'[th] (I1[th] - K0^2 I0[th]) + lambda[th] (I1'[th] - K0^2 I0'[th])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkA = Simplify[dfdth == expected]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaTheta = dfdth dTheta[th]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    expectedDelta = (lambda'[th] (I1[th] - K0^2 I0[th]) + lambda[th] (I1'[th] - K0^2 I0'[th])) dTheta[th]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkB = Simplify[deltaTheta == expectedDelta]
  ], "Input"]
}];

f[th_] := lambda[th] (I1[th] - K0^2 I0[th]);
dfdth = Simplify[D[f[th], th]];
expected = lambda'[th] (I1[th] - K0^2 I0[th]) + lambda[th] (I1'[th] - K0^2 I0'[th]);
checkA = Simplify[dfdth == expected];

deltaTheta = dfdth dTheta[th];
expectedDelta = (lambda'[th] (I1[th] - K0^2 I0[th]) + lambda[th] (I1'[th] - K0^2 I0'[th])) dTheta[th];
checkB = Simplify[deltaTheta == expectedDelta];

nbName = "03_Theta_ChainRule.nb";
logName = "03_theta_chainrule.log";
logLines = {
  "Notebook: " <> nbName,
  "dfdth = " <> ToString[dfdth, InputForm],
  "expected = " <> ToString[expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "deltaTheta = " <> ToString[deltaTheta, InputForm],
  "expectedDelta = " <> ToString[expectedDelta, InputForm],
  "checkB = " <> ToString[checkB, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "03_theta_chainrule",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "checkA" -> checkA,
  "checkB" -> checkB
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[checkA && checkB], Exit[0], Exit[1]];
