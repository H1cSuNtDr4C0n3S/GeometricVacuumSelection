ClearAll["Global`*"];

numK2 = (psiP^2) * ((-4) r^2 BB^2 AAp^2 + 16 AA^3 BB psiP^2 - 8 AA^4 psiP^4 +
    4 r^2 AA BB AAp (BBp + AAp psiP^2) - AA^2 (8 BB^2 + r^2 (BBp + AAp psiP^2)^2)) +
  4 r^2 AA BB psiP (AA BBp + AAp ((-2) BB + AA psiP^2)) psiPP - 4 r^2 AA^2 BB^2 psiPP^2;
denK2 = 4 r^2 AA BB (-BB + AA psiP^2)^3;
k2General = numK2/denK2;

nb = Notebook[{
  Cell["09 - Metric Kernels A,B (Numeric)", "Title"],
  Cell["Numerical verification of metric variation kernels E_A and E_B for the K2_general block.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL = E_A deltaA + E_B deltaB + d/dr(boundary)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_A = ∂L/∂A - d/dr(∂L/∂A'),   E_B = ∂L/∂B - d/dr(∂L/∂B')
  ], "Input"]
}];

Aexpr[x_] := 1 + (1/10) x + (1/30) x^2;
Bexpr[x_] := 2 + (1/7) x + (1/25) x^2;
psi[x_] := (1/3) x + (1/20) x^2 + (1/100) x^3;
dA[x_] := 1 + (1/12) x + (1/60) x^2;
dB[x_] := (3/2) + (1/9) x + (1/80) x^2;

Lsym[x_, a_, b_, ap_, bp_, q1_, q2_] := k2General /. {
  r -> x, AA -> a, BB -> b, AAp -> ap, BBp -> bp, psiP -> q1, psiPP -> q2
};

L0 = Lsym[x, Aexpr[x], Bexpr[x], D[Aexpr[x], x], D[Bexpr[x], x], D[psi[x], x], D[psi[x], {x, 2}]];
Leps = Lsym[
  x,
  Aexpr[x] + eps dA[x],
  Bexpr[x] + eps dB[x],
  D[Aexpr[x] + eps dA[x], x],
  D[Bexpr[x] + eps dB[x], x],
  D[psi[x], x],
  D[psi[x], {x, 2}]
];
deltaL = SeriesCoefficient[Leps, {eps, 0, 1}];

dLdA = D[Lsym[x, a, b, ap, bp, q1, q2], a] /. {
  a -> Aexpr[x], b -> Bexpr[x], ap -> D[Aexpr[x], x], bp -> D[Bexpr[x], x],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};
dLdAp = D[Lsym[x, a, b, ap, bp, q1, q2], ap] /. {
  a -> Aexpr[x], b -> Bexpr[x], ap -> D[Aexpr[x], x], bp -> D[Bexpr[x], x],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};
dLdB = D[Lsym[x, a, b, ap, bp, q1, q2], b] /. {
  a -> Aexpr[x], b -> Bexpr[x], ap -> D[Aexpr[x], x], bp -> D[Bexpr[x], x],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};
dLdBp = D[Lsym[x, a, b, ap, bp, q1, q2], bp] /. {
  a -> Aexpr[x], b -> Bexpr[x], ap -> D[Aexpr[x], x], bp -> D[Bexpr[x], x],
  q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]
};

EAexpr = Together[dLdA - D[dLdAp, x]];
EBexpr = Together[dLdB - D[dLdBp, x]];
boundary = dA[x] dLdAp + dB[x] dLdBp;

residueExpr = N[Together[deltaL - (EAexpr dA[x] + EBexpr dB[x] + D[boundary, x])], 40];
samplePts = {0.7, 0.9, 1.1, 1.3};
residualValues = Table[Chop[N[residueExpr /. x -> x0, 30], 10^-20], {x0, samplePts}];
EAvalues = Table[N[EAexpr /. x -> x0, 30], {x0, samplePts}];
EBvalues = Table[N[EBexpr /. x -> x0, 30], {x0, samplePts}];
maxResidual = Max[Abs[residualValues]];
check = maxResidual < 10^-14;

eaPreview = StringTake[ToString[EAexpr, InputForm], UpTo[700]];
ebPreview = StringTake[ToString[EBexpr, InputForm], UpTo[700]];
eaLeafCount = LeafCount[EAexpr];
ebLeafCount = LeafCount[EBexpr];

nbName = "09_Metric_Kernels_AB_Numeric.nb";
logName = "09_metric_kernels_AB_numeric.log";
logLines = {
  "Notebook: " <> nbName,
  "samplePts = " <> ToString[samplePts, InputForm],
  "EAvalues = " <> ToString[EAvalues, InputForm],
  "EBvalues = " <> ToString[EBvalues, InputForm],
  "EAleafCount = " <> ToString[eaLeafCount, InputForm],
  "EBleafCount = " <> ToString[ebLeafCount, InputForm],
  "EApreview = " <> eaPreview,
  "EBpreview = " <> ebPreview,
  "residualValues = " <> ToString[residualValues, InputForm],
  "maxResidual = " <> ToString[maxResidual, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "09_metric_kernels_AB_numeric",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxResidual" -> ToString[maxResidual, InputForm],
  "EAleafCount" -> eaLeafCount,
  "EBleafCount" -> ebLeafCount
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
