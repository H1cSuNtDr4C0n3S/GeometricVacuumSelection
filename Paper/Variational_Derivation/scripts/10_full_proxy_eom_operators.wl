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
  Cell["10 - Full Proxy EOM Operators", "Title"],
  Cell["Explicit reduced operators E_A, E_B, E_psi for the full local+constraint proxy density.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_A = dL/da - d/dr(dL/da') + d^2/dr^2(dL/da'')
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_B = dL/db - d/dr(dL/db')
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_psi = -d/dr(dL/dpsi') + d^2/dr^2(dL/dpsi'')
  ], "Input"]
}];

A[x_] := 2 + (1/10) x + (1/30) x^2;
B[x_] := 3 + (1/20) x + (1/100) x^2;
psi[x_] := (1/5) x + (1/40) x^2 + (1/200) x^3;

sub = {
  a -> A[x], b -> B[x], ap -> D[A[x], x], bp -> D[B[x], x], app -> D[A[x], {x, 2}],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};

dLda = D[Lsym[x, a, b, ap, bp, app, q1, q2], a] /. sub;
dLdap = D[Lsym[x, a, b, ap, bp, app, q1, q2], ap] /. sub;
dLdapp = D[Lsym[x, a, b, ap, bp, app, q1, q2], app] /. sub;
EAexpr = Together[dLda - D[dLdap, x] + D[dLdapp, {x, 2}]];

dLdb = D[Lsym[x, a, b, ap, bp, app, q1, q2], b] /. sub;
dLdbp = D[Lsym[x, a, b, ap, bp, app, q1, q2], bp] /. sub;
EBexpr = Together[dLdb - D[dLdbp, x]];

dLdq1 = D[Lsym[x, a, b, ap, bp, app, q1, q2], q1] /. sub;
dLdq2 = D[Lsym[x, a, b, ap, bp, app, q1, q2], q2] /. sub;
EpsiExpr = Together[-D[dLdq1, x] + D[dLdq2, {x, 2}]];

samplePts = {0.8, 1.0, 1.2};
EAvalues = Table[N[EAexpr /. x -> x0, 20], {x0, samplePts}];
EBvalues = Table[N[EBexpr /. x -> x0, 20], {x0, samplePts}];
Epsivalues = Table[N[EpsiExpr /. x -> x0, 20], {x0, samplePts}];

eaLeafCount = LeafCount[EAexpr];
ebLeafCount = LeafCount[EBexpr];
epsiLeafCount = LeafCount[EpsiExpr];
eaPreview = StringTake[ToString[EAexpr, InputForm], UpTo[600]];
ebPreview = StringTake[ToString[EBexpr, InputForm], UpTo[600]];
epsiPreview = StringTake[ToString[EpsiExpr, InputForm], UpTo[600]];

finiteCheck = And @@ Flatten[Map[(NumericQ[#] && # =!= Indeterminate && # =!= ComplexInfinity) &, {EAvalues, EBvalues, Epsivalues}, {2}]];
check = TrueQ[finiteCheck && eaLeafCount > 0 && ebLeafCount > 0 && epsiLeafCount > 0];

nbName = "10_Full_Proxy_EOM_Operators.nb";
logName = "10_full_proxy_eom_operators.log";
logLines = {
  "Notebook: " <> nbName,
  "samplePts = " <> ToString[samplePts, InputForm],
  "EAvalues = " <> ToString[EAvalues, InputForm],
  "EBvalues = " <> ToString[EBvalues, InputForm],
  "Epsivalues = " <> ToString[Epsivalues, InputForm],
  "EAleafCount = " <> ToString[eaLeafCount, InputForm],
  "EBleafCount = " <> ToString[ebLeafCount, InputForm],
  "EpsileafCount = " <> ToString[epsiLeafCount, InputForm],
  "EApreview = " <> eaPreview,
  "EBpreview = " <> ebPreview,
  "Epsipreview = " <> epsiPreview,
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "10_full_proxy_eom_operators",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "EAleafCount" -> eaLeafCount,
  "EBleafCount" -> ebLeafCount,
  "EpsileafCount" -> epsiLeafCount
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
