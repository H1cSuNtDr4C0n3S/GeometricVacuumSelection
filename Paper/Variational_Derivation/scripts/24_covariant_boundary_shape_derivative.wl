ClearAll["Global`*"];

nb = Notebook[{
  Cell["24 - Covariant Boundary Shape Derivative", "Title"],
  Cell["Causal-domain indicator as Heaviside level-set: delta chi gives boundary DiracDelta term, reproducing the finite-domain boundary contribution.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    chiEps[r_] := HeavisideTheta[R + eps rho - r]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dchi == D[chiEps[r], eps] /. eps -> 0 == rho DiracDelta[R - r]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI0boundary == Integrate[mu[r] rho DiracDelta[R - r], {r, 0, Infinity}] == rho mu[R]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI1boundary == Integrate[mu[r] k2[r] rho DiracDelta[R - r], {r, 0, Infinity}] == rho mu[R] k2[R]
  ], "Input"]
}];

assum = Element[{R, r, rho}, Reals] && R > 0;

chiEps[r_] := HeavisideTheta[R + eps rho - r];
dchi = FullSimplify[(D[chiEps[r], eps] /. eps -> 0), assum];
checkA = FullSimplify[dchi == rho DiracDelta[R - r], assum];

mu[s_] := Exp[-s] (1 + s^2/5);
k2[s_] := 3/2 + s/3;

dI0Delta = FullSimplify[
  Integrate[mu[s] rho DiracDelta[R - s], {s, 0, Infinity}],
  assum
];
dI0Leib = FullSimplify[
  (D[Integrate[mu[s], {s, 0, R + eps rho}], eps] /. eps -> 0),
  assum
];
checkB = FullSimplify[dI0Delta == dI0Leib, assum];

dI1Delta = FullSimplify[
  Integrate[mu[s] k2[s] rho DiracDelta[R - s], {s, 0, Infinity}],
  assum
];
dI1Leib = FullSimplify[
  (D[Integrate[mu[s] k2[s], {s, 0, R + eps rho}], eps] /. eps -> 0),
  assum
];
checkC = FullSimplify[dI1Delta == dI1Leib, assum];

(* numeric finite-difference consistency *)
R0 = 13/10;
rho0 = 7/100;
h = 1/10^6;

i0Exact[ee_] := N[Integrate[mu[s], {s, 0, R0 + ee rho0}], 50];
i1Exact[ee_] := N[Integrate[mu[s] k2[s], {s, 0, R0 + ee rho0}], 50];

dI0Num = N[(i0Exact[h] - i0Exact[-h])/(2 h), 30];
dI1Num = N[(i1Exact[h] - i1Exact[-h])/(2 h), 30];

dI0Pred = N[rho0 mu[R0], 30];
dI1Pred = N[rho0 mu[R0] k2[R0], 30];

err0 = N[Abs[dI0Num - dI0Pred], 30];
err1 = N[Abs[dI1Num - dI1Pred], 30];

checkNum0 = err0 < 10^-8;
checkNum1 = err1 < 10^-8;

check = TrueQ[checkA && checkB && checkC && checkNum0 && checkNum1];

nbName = "24_Covariant_Boundary_Shape_Derivative.nb";
logName = "24_covariant_boundary_shape_derivative.log";
logLines = {
  "Notebook: " <> nbName,
  "dchi = " <> ToString[dchi, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "dI0Delta = " <> ToString[dI0Delta, InputForm],
  "dI0Leib = " <> ToString[dI0Leib, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "dI1Delta = " <> ToString[dI1Delta, InputForm],
  "dI1Leib = " <> ToString[dI1Leib, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "dI0Num = " <> ToString[dI0Num, InputForm],
  "dI0Pred = " <> ToString[dI0Pred, InputForm],
  "err0 = " <> ToString[err0, InputForm],
  "checkNum0 = " <> ToString[checkNum0, InputForm],
  "dI1Num = " <> ToString[dI1Num, InputForm],
  "dI1Pred = " <> ToString[dI1Pred, InputForm],
  "err1 = " <> ToString[err1, InputForm],
  "checkNum1 = " <> ToString[checkNum1, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "24_covariant_boundary_shape_derivative",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
