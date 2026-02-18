ClearAll["Global`*"];

numK2 = (psiP^2) * ((-4) r^2 BB^2 AAp^2 + 16 AA^3 BB psiP^2 - 8 AA^4 psiP^4 +
    4 r^2 AA BB AAp (BBp + AAp psiP^2) - AA^2 (8 BB^2 + r^2 (BBp + AAp psiP^2)^2)) +
  4 r^2 AA BB psiP (AA BBp + AAp ((-2) BB + AA psiP^2)) psiPP - 4 r^2 AA^2 BB^2 psiPP^2;
denK2 = 4 r^2 AA BB (-BB + AA psiP^2)^3;
k2General = numK2/denK2;

Rexpr[x_, a_, b_, ap_, bp_, app_] :=
  (-x^2 b ap^2 - 4 a^2 ((-1 + b) b + x bp) + x a (ap (4 b - x bp) + 2 x b app)) /
  (2 x^2 a^2 b^2);
Xexpr[a_, b_, q1_] := 1/a - q1^2/b;

Lsym[x_, a_, b_, ap_, bp_, app_, q1_, q2_] := Module[
  {kappa = 1.7, alpha = 0.3, lam = 0.2, k0 = 0.4, ll = 0.05},
  4 Pi x^2 Sqrt[a b] ((Rexpr[x, a, b, ap, bp, app] - 2 ll)/(2 kappa) - (alpha/2) Xexpr[a, b, q1]) +
  lam ((k2General /. {r -> x, AA -> a, BB -> b, AAp -> ap, BBp -> bp, psiP -> q1, psiPP -> q2}) - k0^2)
];

nb = Notebook[{
  Cell["11 - Full Proxy Variational Checks", "Title"],
  Cell["Channel-by-channel first-variation consistency checks for A, B, psi on the full proxy density.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL_A = E_A deltaA + d/dr(boundary_A)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL_B = E_B deltaB + d/dr(boundary_B)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL_psi = E_psi deltaPsi + d/dr(boundary_psi)
  ], "Input"]
}];

A[x_] := 2 + (1/10) x + (1/30) x^2;
B[x_] := 3 + (1/20) x + (1/100) x^2;
psi[x_] := (1/5) x + (1/40) x^2 + (1/200) x^3;
dA[x_] := 1 + (1/12) x + (1/60) x^2;
dB[x_] := (3/2) + (1/9) x + (1/80) x^2;
dpsi[x_] := 1 + (1/13) x + (1/70) x^2;

sub = {
  a -> A[x], b -> B[x], ap -> D[A[x], x], bp -> D[B[x], x], app -> D[A[x], {x, 2}],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};

(* A channel *)
dLda = D[Lsym[x, a, b, ap, bp, app, q1, q2], a] /. sub;
dLdap = D[Lsym[x, a, b, ap, bp, app, q1, q2], ap] /. sub;
dLdapp = D[Lsym[x, a, b, ap, bp, app, q1, q2], app] /. sub;
EAexpr = dLda - D[dLdap, x] + D[dLdapp, {x, 2}];

LepsA = Lsym[
  x,
  A[x] + eps dA[x],
  B[x],
  D[A[x] + eps dA[x], x],
  D[B[x], x],
  D[A[x] + eps dA[x], {x, 2}],
  D[psi[x], x],
  D[psi[x], {x, 2}]
];
deltaAexpr = SeriesCoefficient[LepsA, {eps, 0, 1}];
boundA = dA[x] dLdap + D[dA[x], x] dLdapp - dA[x] D[dLdapp, x];
resAexpr = N[Together[deltaAexpr - (EAexpr dA[x] + D[boundA, x])], 40];

(* B channel *)
dLdb = D[Lsym[x, a, b, ap, bp, app, q1, q2], b] /. sub;
dLdbp = D[Lsym[x, a, b, ap, bp, app, q1, q2], bp] /. sub;
EBexpr = dLdb - D[dLdbp, x];

LepsB = Lsym[
  x,
  A[x],
  B[x] + eps dB[x],
  D[A[x], x],
  D[B[x] + eps dB[x], x],
  D[A[x], {x, 2}],
  D[psi[x], x],
  D[psi[x], {x, 2}]
];
deltaBexpr = SeriesCoefficient[LepsB, {eps, 0, 1}];
boundB = dB[x] dLdbp;
resBexpr = N[Together[deltaBexpr - (EBexpr dB[x] + D[boundB, x])], 40];

(* psi channel *)
dLdq1 = D[Lsym[x, a, b, ap, bp, app, q1, q2], q1] /. sub;
dLdq2 = D[Lsym[x, a, b, ap, bp, app, q1, q2], q2] /. sub;
Epsiexpr = -D[dLdq1, x] + D[dLdq2, {x, 2}];

LepsPsi = Lsym[
  x,
  A[x],
  B[x],
  D[A[x], x],
  D[B[x], x],
  D[A[x], {x, 2}],
  D[psi[x], x] + eps D[dpsi[x], x],
  D[psi[x], {x, 2}] + eps D[dpsi[x], {x, 2}]
];
deltaPsiexpr = SeriesCoefficient[LepsPsi, {eps, 0, 1}];
boundPsi = dpsi[x] dLdq1 + D[dpsi[x], x] dLdq2 - dpsi[x] D[dLdq2, x];
resPsiexpr = N[Together[deltaPsiexpr - (Epsiexpr dpsi[x] + D[boundPsi, x])], 40];

samplePts = {0.8, 1.0, 1.2};
resAvals = Table[Chop[N[resAexpr /. x -> x0, 30], 10^-16], {x0, samplePts}];
resBvals = Table[Chop[N[resBexpr /. x -> x0, 30], 10^-16], {x0, samplePts}];
resPsivals = Table[Chop[N[resPsiexpr /. x -> x0, 30], 10^-16], {x0, samplePts}];

maxResA = Max[Abs[resAvals]];
maxResB = Max[Abs[resBvals]];
maxResPsi = Max[Abs[resPsivals]];
checkA = maxResA < 10^-12;
checkB = maxResB < 10^-12;
checkPsi = maxResPsi < 10^-12;
check = TrueQ[checkA && checkB && checkPsi];

nbName = "11_Full_Proxy_Variational_Checks.nb";
logName = "11_full_proxy_variational_checks.log";
logLines = {
  "Notebook: " <> nbName,
  "samplePts = " <> ToString[samplePts, InputForm],
  "resAvals = " <> ToString[resAvals, InputForm],
  "resBvals = " <> ToString[resBvals, InputForm],
  "resPsivals = " <> ToString[resPsivals, InputForm],
  "maxResA = " <> ToString[maxResA, InputForm],
  "maxResB = " <> ToString[maxResB, InputForm],
  "maxResPsi = " <> ToString[maxResPsi, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkPsi = " <> ToString[checkPsi, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "11_full_proxy_variational_checks",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxResA" -> ToString[maxResA, InputForm],
  "maxResB" -> ToString[maxResB, InputForm],
  "maxResPsi" -> ToString[maxResPsi, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
