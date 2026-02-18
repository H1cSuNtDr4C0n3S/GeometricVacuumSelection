ClearAll["Global`*"];

nb = Notebook[{
  Cell["14 - Constraint Variation with Finite-Domain Boundary", "Title"],
  Cell["First variation of C = I1 - K0^2 I0 with explicit bulk+boundary split on R -> R + eps rho.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI0 = dI0bulk + rho f0(R),   deltaI1 = dI1bulk + rho f1(R)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaC = deltaI1 - K0^2 deltaI0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaC = (dI1bulk - K0^2 dI0bulk) + rho (f1(R) - K0^2 f0(R))
  ], "Input"]
}];

k0 = Symbol["K0"];
di0 = Symbol["dI0bulk"];
di1 = Symbol["dI1bulk"];
f0R = Symbol["f0R"];
f1R = Symbol["f1R"];

deltaI0 = di0 + rho f0R;
deltaI1 = di1 + rho f1R;

deltaC = Expand[deltaI1 - k0^2 deltaI0];
bulkPart = Simplify[Coefficient[deltaC, rho, 0]];
boundaryPart = Simplify[Coefficient[deltaC, rho, 1]];

expectedBulk = di1 - k0^2 di0;
expectedBoundary = f1R - k0^2 f0R;
checkA = Simplify[bulkPart == expectedBulk];
checkB = Simplify[boundaryPart == expectedBoundary];
checkC = Simplify[(boundaryPart /. f1R -> k0^2 f0R) == 0];

(* independent finite-difference check *)
f0fun[s_] := 1 + s + (2/5) s^2;
f1fun[s_] := 2 + (3/10) s + (1/2) s^2;
df0fun[s_] := 1/5 + (1/10) s;
df1fun[s_] := 1/7 + (2/5) s;

R0 = 7/5;
rho0 = 3/40;
k00 = 11/10;
h = 1/10^6;

cExact[ee_] := N[
  Integrate[f1fun[s] + ee df1fun[s], {s, 0, R0 + ee rho0}] -
  k00^2 Integrate[f0fun[s] + ee df0fun[s], {s, 0, R0 + ee rho0}],
  50
];

deltaNum = N[(cExact[h] - cExact[-h])/(2 h), 30];

di0v = N[Integrate[df0fun[s], {s, 0, R0}], 50];
di1v = N[Integrate[df1fun[s], {s, 0, R0}], 50];
deltaPred = N[
  (di1v + rho0 f1fun[R0]) - k00^2 (di0v + rho0 f0fun[R0]),
  30
];
absErr = N[Abs[deltaNum - deltaPred], 30];
checkNum = absErr < 10^-8;

check = TrueQ[checkA && checkB && checkC && checkNum];

nbName = "14_Constraint_Finite_Domain_Boundary.nb";
logName = "14_constraint_finite_domain_boundary.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaC = " <> ToString[deltaC, InputForm],
  "bulkPart = " <> ToString[bulkPart, InputForm],
  "expectedBulk = " <> ToString[expectedBulk, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "boundaryPart = " <> ToString[boundaryPart, InputForm],
  "expectedBoundary = " <> ToString[expectedBoundary, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkBoundaryNeutral = " <> ToString[checkC, InputForm],
  "deltaNum = " <> ToString[deltaNum, InputForm],
  "deltaPred = " <> ToString[deltaPred, InputForm],
  "absErr = " <> ToString[absErr, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "14_constraint_finite_domain_boundary",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "absErr" -> ToString[absErr, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
