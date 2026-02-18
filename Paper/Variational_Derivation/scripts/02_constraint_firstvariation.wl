ClearAll["Global`*"];

nb = Notebook[{
  Cell["02 - Constraint First Variation", "Title"],
  Cell["First variation of the leafwise constrained integrand.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    integrand[th_] := lambda[th] (I1[th] - K0^2 I0[th])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    perturbed = Expand[
      (lambda[th] + eps eta[th]) ((I1[th] + eps dI1[th]) - K0^2 (I0[th] + eps dI0[th]))
    ]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaIntegrand = SeriesCoefficient[perturbed, {eps, 0, 1}] // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    expected = eta[th] (I1[th] - K0^2 I0[th]) + lambda[th] (dI1[th] - K0^2 dI0[th])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkA = Simplify[deltaIntegrand == expected]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q[th_] := I1[th]/I0[th]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    expectedQ = eta[th] I0[th] (q[th] - K0^2) + lambda[th] (dI0[th] (q[th] - K0^2) + I0[th] dQ[th])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    fromIQ = expected /. {I1[th] -> I0[th] q[th], dI1[th] -> dI0[th] q[th] + I0[th] dQ[th]} // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkB = Simplify[fromIQ == expectedQ]
  ], "Input"]
}];

integrand[th_] := lambda[th] (I1[th] - K0^2 I0[th]);
perturbed = Expand[
  (lambda[th] + eps eta[th]) ((I1[th] + eps dI1[th]) - K0^2 (I0[th] + eps dI0[th]))
];
deltaIntegrand = Simplify[SeriesCoefficient[perturbed, {eps, 0, 1}]];
expected = eta[th] (I1[th] - K0^2 I0[th]) + lambda[th] (dI1[th] - K0^2 dI0[th]);
checkA = Simplify[deltaIntegrand == expected];

q[th_] := I1[th]/I0[th];
expectedQ = eta[th] I0[th] (q[th] - K0^2) + lambda[th] (dI0[th] (q[th] - K0^2) + I0[th] dQ[th]);
fromIQ = Simplify[expected /. {I1[th] -> I0[th] q[th], dI1[th] -> dI0[th] q[th] + I0[th] dQ[th]}];
checkB = Simplify[fromIQ == expectedQ];

nbName = "02_Constraint_FirstVariation.nb";
logName = "02_constraint_firstvariation.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaIntegrand = " <> ToString[deltaIntegrand, InputForm],
  "expected = " <> ToString[expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "fromIQ = " <> ToString[fromIQ, InputForm],
  "expectedQ = " <> ToString[expectedQ, InputForm],
  "checkB = " <> ToString[checkB, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "02_constraint_firstvariation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "checkA" -> checkA,
  "checkB" -> checkB
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[checkA && checkB], Exit[0], Exit[1]];
