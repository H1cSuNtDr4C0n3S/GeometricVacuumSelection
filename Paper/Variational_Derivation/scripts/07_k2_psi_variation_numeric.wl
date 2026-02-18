ClearAll["Global`*"];

numK2 = (psiP^2) * ((-4) r^2 BB^2 AAp^2 + 16 AA^3 BB psiP^2 - 8 AA^4 psiP^4 +
    4 r^2 AA BB AAp (BBp + AAp psiP^2) - AA^2 (8 BB^2 + r^2 (BBp + AAp psiP^2)^2)) +
  4 r^2 AA BB psiP (AA BBp + AAp ((-2) BB + AA psiP^2)) psiPP - 4 r^2 AA^2 BB^2 psiPP^2;
denK2 = 4 r^2 AA BB (-BB + AA psiP^2)^3;
k2General = numK2/denK2;

nb = Notebook[{
  Cell["07 - K2 Psi Variation (Numeric Check)", "Title"],
  Cell["Numerical verification of second-order EL identity for the K2_general block with nontrivial A(r), B(r).", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lsym[x_, q1_, q2_] := K2general /. {
      r -> x, AA -> Aexpr[x], BB -> Bexpr[x], AAp -> Aexpr'[x], BBp -> Bexpr'[x], psiP -> q1, psiPP -> q2
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL = coeff[eps, 1] of Lsym[x, psi'[x] + eps eta'[x], psi''[x] + eps eta''[x]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    ELpsi = -d/dx(∂L/∂psi') + d^2/dx^2(∂L/∂psi'')
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    residue = deltaL - (ELpsi eta + d/dx boundary)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    check = Max[Abs[residue(x_i)]] < 10^-14 on sample points
  ], "Input"]
}];

Aexpr[x_] := 1 + (1/10) x + (1/30) x^2;
Bexpr[x_] := 2 + (1/7) x + (1/25) x^2;
psi[x_] := (1/3) x + (1/20) x^2 + (1/100) x^3;
eta[x_] := 1 + (1/11) x + (1/50) x^2;

Lsym[x_, q1_, q2_] := k2General /. {
  r -> x,
  AA -> Aexpr[x],
  BB -> Bexpr[x],
  AAp -> D[Aexpr[x], x],
  BBp -> D[Bexpr[x], x],
  psiP -> q1,
  psiPP -> q2
};

L0 = Lsym[x, D[psi[x], x], D[psi[x], {x, 2}]];
Leps = Lsym[x, D[psi[x], x] + eps D[eta[x], x], D[psi[x], {x, 2}] + eps D[eta[x], {x, 2}]];
deltaL = SeriesCoefficient[Leps, {eps, 0, 1}];

dLdq1 = D[Lsym[x, q1, q2], q1] /. {q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]};
dLdq2 = D[Lsym[x, q1, q2], q2] /. {q1 -> D[psi[x], x], q2 -> D[psi[x], {x, 2}]};
ELpsi = -D[dLdq1, x] + D[dLdq2, {x, 2}];

boundary = eta[x] dLdq1 + D[eta[x], x] dLdq2 - eta[x] D[dLdq2, x];
residueExpr = N[Together[deltaL - (ELpsi eta[x] + D[boundary, x])], 40];

samplePts = {0.7, 0.9, 1.1, 1.3};
residualValues = Table[Chop[N[residueExpr /. x -> x0, 30], 10^-20], {x0, samplePts}];
maxResidual = Max[Abs[residualValues]];
check = maxResidual < 10^-14;

nbName = "07_K2_Psi_Variation_Numeric.nb";
logName = "07_k2_psi_variation_numeric.log";
logLines = {
  "Notebook: " <> nbName,
  "samplePts = " <> ToString[samplePts, InputForm],
  "residualValues = " <> ToString[residualValues, InputForm],
  "maxResidual = " <> ToString[maxResidual, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "07_k2_psi_variation_numeric",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxResidual" -> ToString[maxResidual, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
