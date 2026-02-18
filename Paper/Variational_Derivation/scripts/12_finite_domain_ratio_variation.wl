ClearAll["Global`*"];

nb = Notebook[{
  Cell["12 - Finite-Domain Ratio Variation", "Title"],
  Cell["First variation of Q = I1/I0 on a finite domain including moving-boundary contribution.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ = ((deltaI1) I0 - I1 (deltaI0))/I0^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI0 = dI0bulk + rho f0(R),   deltaI1 = dI1bulk + rho f1(R)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    boundaryPiece = rho (f1(R) I0 - I1 f0(R))/I0^2
  ], "Input"]
}];

i0 = Symbol["I0"];
i1 = Symbol["I1"];
di0 = Symbol["dI0bulk"];
di1 = Symbol["dI1bulk"];
f0R = Symbol["f0R"];
f1R = Symbol["f1R"];

deltaI0 = di0 + rho f0R;
deltaI1 = di1 + rho f1R;

deltaQ = Simplify[SeriesCoefficient[(i1 + eps deltaI1)/(i0 + eps deltaI0), {eps, 0, 1}]];
expected = Simplify[(deltaI1 i0 - i1 deltaI0)/i0^2];
checkA = Simplify[deltaQ == expected];

boundaryPiece = Simplify[Coefficient[deltaQ, rho]];
expectedBoundary = Simplify[(f1R i0 - i1 f0R)/i0^2];
checkB = Simplify[boundaryPiece == expectedBoundary];

(* numerical consistency check against direct finite-difference on integrals *)
f0fun[s_] := 1 + (3/10) s + (1/10) s^2;
f1fun[s_] := 1/2 + s^2;
df0fun[s_] := 1/5 + (1/20) s;
df1fun[s_] := (2/5) s;

R0 = 13/10;
rho0 = 7/100;
h = 1/10^6;

qExact[ee_] := N[
  Integrate[f1fun[s] + ee df1fun[s], {s, 0, R0 + ee rho0}]/
  Integrate[f0fun[s] + ee df0fun[s], {s, 0, R0 + ee rho0}],
  40
];

deltaNum = N[(qExact[h] - qExact[-h])/(2 h), 30];

i0v = N[Integrate[f0fun[s], {s, 0, R0}], 40];
i1v = N[Integrate[f1fun[s], {s, 0, R0}], 40];
di0v = N[Integrate[df0fun[s], {s, 0, R0}], 40];
di1v = N[Integrate[df1fun[s], {s, 0, R0}], 40];

deltaPred = N[((di1v + rho0 f1fun[R0]) i0v - i1v (di0v + rho0 f0fun[R0]))/i0v^2, 30];
absErr = N[Abs[deltaNum - deltaPred], 30];
checkNum = absErr < 10^-8;

check = TrueQ[checkA && checkB && checkNum];

nbName = "12_Finite_Domain_Ratio_Variation.nb";
logName = "12_finite_domain_ratio_variation.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaQ = " <> ToString[deltaQ, InputForm],
  "expected = " <> ToString[expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "boundaryPiece = " <> ToString[boundaryPiece, InputForm],
  "expectedBoundary = " <> ToString[expectedBoundary, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "deltaNum = " <> ToString[deltaNum, InputForm],
  "deltaPred = " <> ToString[deltaPred, InputForm],
  "absErr = " <> ToString[absErr, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "12_finite_domain_ratio_variation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "absErr" -> ToString[absErr, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
