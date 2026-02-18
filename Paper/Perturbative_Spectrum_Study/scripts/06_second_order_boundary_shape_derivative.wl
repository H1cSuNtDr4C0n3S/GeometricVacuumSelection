ClearAll["Global`*"];

nb = Notebook[{
  Cell["06 - Second-Order Boundary Shape Derivative", "Title"],
  Cell["Second-order variation of causal-domain indicator and consistency with Leibniz derivatives for moving boundaries.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    chiEps[r_] := HeavisideTheta[R + eps rho + eps^2 sigma/2 - r]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dchi1 == D[chiEps[r], eps] /. eps -> 0 == rho DiracDelta[R - r]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dchi2 == D[chiEps[r], {eps, 2}] /. eps -> 0 == sigma DiracDelta[R - r] + rho^2 DiracDelta'[R - r]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    I2 == sigma f[R] + rho^2 f'[R]
  ], "Input"]
}];

assum = Element[{R, r, rho, sigma}, Reals] && R > 0;

chiEps[r_] := HeavisideTheta[R + eps rho + eps^2 sigma/2 - r];
dchi1 = FullSimplify[(D[chiEps[r], eps] /. eps -> 0), assum];
dchi2 = FullSimplify[(D[chiEps[r], {eps, 2}] /. eps -> 0), assum];

checkA = FullSimplify[dchi1 == rho DiracDelta[r - R], assum];

f[s_] := Exp[-s] (1 + s + s^2/5);
fp[s_] := D[f[s], s];

iExact[ee_] := Integrate[f[s], {s, 0, R + ee rho + ee^2 sigma/2}];
i1Leib = FullSimplify[D[iExact[eps], eps] /. eps -> 0, assum];
i2Leib = FullSimplify[D[iExact[eps], {eps, 2}] /. eps -> 0, assum];

i1Dist = FullSimplify[Integrate[f[s] (dchi1 /. r -> s), {s, 0, Infinity}], assum];
i2Dist = FullSimplify[Integrate[f[s] (dchi2 /. r -> s), {s, 0, Infinity}], assum];

i1Expected = FullSimplify[rho f[R], assum];
i2Expected = FullSimplify[sigma f[R] + rho^2 fp[R], assum];

checkB = FullSimplify[i1Leib == i1Dist, assum];
checkC = FullSimplify[i2Leib == i2Dist, assum];
checkD = FullSimplify[i2Leib == i2Expected, assum];

R0 = 13/10;
rho0 = 7/100;
sigma0 = 3/50;
h = 1/10^4;

inum[ee_] := N[iExact[ee] /. {R -> R0, rho -> rho0, sigma -> sigma0}, 50];

i1Num = N[(inum[h] - inum[-h])/(2 h), 30];
i1Pred = N[i1Expected /. {R -> R0, rho -> rho0}, 30];
i2Num = N[(inum[h] - 2 inum[0] + inum[-h])/h^2, 30];
i2Pred = N[i2Expected /. {R -> R0, rho -> rho0, sigma -> sigma0}, 30];

err1 = N[Abs[i1Num - i1Pred], 30];
err2 = N[Abs[i2Num - i2Pred], 30];

checkNumA = err1 < 1*^-7;
checkNumB = err2 < 1*^-5;

check = TrueQ[checkA && checkB && checkC && checkD && checkNumA && checkNumB];

nbName = "06_Second_Order_Boundary_Shape_Derivative.nb";
logName = "06_second_order_boundary_shape_derivative.log";
logLines = {
  "Notebook: " <> nbName,
  "dchi1 = " <> ToString[dchi1, InputForm],
  "dchi2 = " <> ToString[dchi2, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "i1Leib = " <> ToString[i1Leib, InputForm],
  "i1Dist = " <> ToString[i1Dist, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "i2Leib = " <> ToString[i2Leib, InputForm],
  "i2Dist = " <> ToString[i2Dist, InputForm],
  "i2Expected = " <> ToString[i2Expected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "i1Num = " <> ToString[i1Num, InputForm],
  "i1Pred = " <> ToString[i1Pred, InputForm],
  "err1 = " <> ToString[err1, InputForm],
  "i2Num = " <> ToString[i2Num, InputForm],
  "i2Pred = " <> ToString[i2Pred, InputForm],
  "err2 = " <> ToString[err2, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "06_second_order_boundary_shape_derivative",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
