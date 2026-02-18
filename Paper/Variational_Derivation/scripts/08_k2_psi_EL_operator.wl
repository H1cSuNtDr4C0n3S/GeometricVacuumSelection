ClearAll["Global`*"];

numK2 = (psiP^2) * ((-4) r^2 BB^2 AAp^2 + 16 AA^3 BB psiP^2 - 8 AA^4 psiP^4 +
    4 r^2 AA BB AAp (BBp + AAp psiP^2) - AA^2 (8 BB^2 + r^2 (BBp + AAp psiP^2)^2)) +
  4 r^2 AA BB psiP (AA BBp + AAp ((-2) BB + AA psiP^2)) psiPP - 4 r^2 AA^2 BB^2 psiPP^2;
denK2 = 4 r^2 AA BB (-BB + AA psiP^2)^3;
k2General = numK2/denK2;

nb = Notebook[{
  Cell["08 - K2 Psi EL Operator", "Title"],
  Cell["Explicit EL_psi operator for the K2_general block on a nontrivial spherical background.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lsym[x_, q1_, q2_] := K2general /. {
      r -> x, AA -> Aexpr[x], BB -> Bexpr[x], AAp -> Aexpr'[x], BBp -> Bexpr'[x], psiP -> q1, psiPP -> q2
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    ELpsi[x] = -d/dx(∂L/∂psi') + d^2/dx^2(∂L/∂psi'')
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    check = deltaL - (ELpsi eta + d/dx boundary) == 0 on sample points
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
ELpsiExpr = Together[-D[dLdq1, x] + D[dLdq2, {x, 2}]];

boundary = eta[x] dLdq1 + D[eta[x], x] dLdq2 - eta[x] D[dLdq2, x];
residueExpr = N[Together[deltaL - (ELpsiExpr eta[x] + D[boundary, x])], 40];

samplePts = {0.7, 0.9, 1.1, 1.3};
residualValues = Table[Chop[N[residueExpr /. x -> x0, 30], 10^-20], {x0, samplePts}];
maxResidual = Max[Abs[residualValues]];
ELvalues = Table[N[ELpsiExpr /. x -> x0, 30], {x0, samplePts}];
check = maxResidual < 10^-14;

elPreview = StringTake[ToString[ELpsiExpr, InputForm], UpTo[900]];
elLeafCount = LeafCount[ELpsiExpr];

nbName = "08_K2_Psi_EL_Operator.nb";
logName = "08_k2_psi_EL_operator.log";
logLines = {
  "Notebook: " <> nbName,
  "samplePts = " <> ToString[samplePts, InputForm],
  "ELvalues = " <> ToString[ELvalues, InputForm],
  "ELleafCount = " <> ToString[elLeafCount, InputForm],
  "ELpreview = " <> elPreview,
  "residualValues = " <> ToString[residualValues, InputForm],
  "maxResidual = " <> ToString[maxResidual, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "08_k2_psi_EL_operator",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxResidual" -> ToString[maxResidual, InputForm],
  "ELleafCount" -> elLeafCount
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
