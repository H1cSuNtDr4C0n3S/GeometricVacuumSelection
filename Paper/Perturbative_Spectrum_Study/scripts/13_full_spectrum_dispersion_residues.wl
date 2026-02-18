ClearAll["Global`*"];

nb = Notebook[{
  Cell["13 - Full Spectrum Dispersion and Residue Map", "Title"],
  Cell["Unified tensor/vector/scalar dispersion analysis with scalar mixing and mass terms. Includes UV-speed recovery, positivity checks, and IR-flow scaling.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    det[x Ks - (q2 Gs + Ms)] == 0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    xTensor == (q2 g2 + mt2)/c2 && xVector == (q2 g1 + mv2)/c1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    cScalarUV == Eigenvalues[Inverse[Ks].Gs]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    coeffs -> coeffs0 + ncoeffs/I0
  ], "Input"]
}];

Ks = {{ks, km}, {km, kw}};
Gs = {{gs, gm}, {gm, gw}};
Ms = {{ms, mm}, {mm, mw}};

x = xw2;
A = q2 Gs + Ms;

poly = FullSimplify[Expand[Det[x Ks - A]]];
a2 = FullSimplify[Coefficient[poly, x, 2]];
a1 = FullSimplify[Coefficient[poly, x, 1]];
a0 = FullSimplify[Coefficient[poly, x, 0]];

adjKs = {{kw, -km}, {-km, ks}};
a2Expected = FullSimplify[Det[Ks]];
a1Expected = FullSimplify[-Tr[adjKs.A]];
a0Expected = FullSimplify[Det[A]];

checkA = FullSimplify[a2 == a2Expected];
checkB = FullSimplify[a1 == a1Expected];
checkC = FullSimplify[a0 == a0Expected];

disc = FullSimplify[a1^2 - 4 a2 a0];
xPlus = FullSimplify[(-a1 + Sqrt[disc])/(2 a2)];
xMinus = FullSimplify[(-a1 - Sqrt[disc])/(2 a2)];

checkD = FullSimplify[xPlus + xMinus == -a1/a2];
checkE = FullSimplify[xPlus*xMinus == a0/a2];

decRule = {km -> 0, gm -> 0, mm -> 0};
checkF1 = FullSimplify[((xPlus + xMinus) /. decRule) == (q2 gs + ms)/ks + (q2 gw + mw)/kw];
checkF2 = FullSimplify[((xPlus*xMinus) /. decRule) == ((q2 gs + ms)/ks) ((q2 gw + mw)/kw)];
checkF = TrueQ[checkF1 && checkF2];

xTensor = FullSimplify[(q2 g2 + mt2)/c2];
xVector = FullSimplify[(q2 g1 + mv2)/c1];
cScalarUV = FullSimplify[Eigenvalues[Inverse[Ks].Gs]];

checkG = FreeQ[xTensor, {ks, km, kw, gs, gm, gw, ms, mm, mw}] &&
  FreeQ[xVector, {ks, km, kw, gs, gm, gw, ms, mm, mw}];
checkH1 = FullSimplify[Total[cScalarUV /. decRule] == gs/ks + gw/kw];
checkH2 = FullSimplify[Times @@ (cScalarUV /. decRule) == (gs*gw)/(ks*kw)];
checkH = TrueQ[checkH1 && checkH2];

(* IR-flow parametrization *)
ksIR = ks0 + nks/i0;
kmIR = km0 + nkm/i0;
kwIR = kw0 + nkw/i0;

gsIR = gs0 + ngs/i0;
gmIR = gm0 + ngm/i0;
gwIR = gw0 + ngw/i0;

msIR = ms0 + nms/i0;
mmIR = mm0 + nmm/i0;
mwIR = mw0 + nmw/i0;

c2IR = c20 + nc2/i0;
g2IR = g20 + ng2/i0;
mt2IR = mt20 + nmt2/i0;

c1IR = c10 + nc1/i0;
g1IR = g10 + ng1/i0;
mv2IR = mv20 + nmv2/i0;

KsIR = {{ksIR, kmIR}, {kmIR, kwIR}};
GsIR = {{gsIR, gmIR}, {gmIR, gwIR}};
MsIR = {{msIR, mmIR}, {mmIR, mwIR}};

numRule = {
  ks0 -> 6/5, km0 -> 1/12, kw0 -> 7/5,
  gs0 -> 4/5, gm0 -> 1/20, gw0 -> 9/10,
  ms0 -> 1/200, mm0 -> 1/500, mw0 -> 1/150,
  nks -> 1/3, nkm -> 1/40, nkw -> -1/4,
  ngs -> 1/5, ngm -> 1/50, ngw -> -1/6,
  nms -> 1/300, nmm -> 1/700, nmw -> 1/250,
  c20 -> 1, g20 -> 1, mt20 -> 1/300,
  c10 -> 19/20, g10 -> 9/10, mv20 -> 1/250,
  nc2 -> 1/4, ng2 -> -1/6, nmt2 -> 1/400,
  nc1 -> -1/5, ng1 -> 1/8, nmv2 -> 1/350
};

Lvals = {10, 20, 40};
q2Vals = {1/10, 1, 4};
q2UV = 10^6;

scalarRootsNum[q2val_, i0val_] := Module[{ksN, gsN, msN},
  ksN = N[KsIR /. numRule /. i0 -> i0val, 60];
  gsN = N[GsIR /. numRule /. i0 -> i0val, 60];
  msN = N[MsIR /. numRule /. i0 -> i0val, 60];
  Sort[N[Eigenvalues[LinearSolve[ksN, q2val gsN + msN]], 40]]
];

scalarUVNum[i0val_] := Module[{ksN, gsN},
  ksN = N[KsIR /. numRule /. i0 -> i0val, 60];
  gsN = N[GsIR /. numRule /. i0 -> i0val, 60];
  Sort[N[Eigenvalues[LinearSolve[ksN, gsN]], 40]]
];

tensorRootNum[q2val_, i0val_] := N[((q2val g2IR + mt2IR)/c2IR) /. numRule /. i0 -> i0val, 40];
vectorRootNum[q2val_, i0val_] := N[((q2val g1IR + mv2IR)/c1IR) /. numRule /. i0 -> i0val, 40];

ksEigVals = Table[
  Sort[N[Eigenvalues[KsIR /. numRule /. i0 -> LL^3], 40]],
  {LL, Lvals}
];

gsEigVals = Table[
  Sort[N[Eigenvalues[GsIR /. numRule /. i0 -> LL^3], 40]],
  {LL, Lvals}
];

scalarRootVals = Table[
  scalarRootsNum[qv, LL^3],
  {LL, Lvals}, {qv, q2Vals}
];

tensorRootVals = Table[
  Table[tensorRootNum[qv, LL^3], {qv, q2Vals}],
  {LL, Lvals}
];

vectorRootVals = Table[
  Table[vectorRootNum[qv, LL^3], {qv, q2Vals}],
  {LL, Lvals}
];

scalarUVFromRoots = Table[
  Sort[N[scalarRootsNum[q2UV, LL^3]/q2UV, 40]],
  {LL, Lvals}
];
scalarUVEigs = Table[
  scalarUVNum[LL^3],
  {LL, Lvals}
];
uvErr = Table[
  N[Max[Abs[scalarUVFromRoots[[ii]] - scalarUVEigs[[ii]]]], 40],
  {ii, 1, Length[Lvals]}
];

KsLoc = N[{{ks0, km0}, {km0, kw0}} /. numRule, 60];
GsLoc = N[{{gs0, gm0}, {gm0, gw0}} /. numRule, 60];
MsLoc = N[{{ms0, mm0}, {mm0, mw0}} /. numRule, 60];

scalarLoc[qv_] := Sort[N[Eigenvalues[LinearSolve[KsLoc, qv GsLoc + MsLoc]], 40]];
tensorLoc[qv_] := N[((qv g20 + mt20)/c20) /. numRule, 40];
vectorLoc[qv_] := N[((qv g10 + mv20)/c10) /. numRule, 40];

q2Ref = 1;
idxRef = 2;

scalarErr = Table[
  N[Max[Abs[scalarRootVals[[ii, idxRef]] - scalarLoc[q2Ref]]], 40],
  {ii, 1, Length[Lvals]}
];

tensorErr = Table[
  N[Abs[tensorRootVals[[ii, idxRef]] - tensorLoc[q2Ref]], 40],
  {ii, 1, Length[Lvals]}
];

vectorErr = Table[
  N[Abs[vectorRootVals[[ii, idxRef]] - vectorLoc[q2Ref]], 40],
  {ii, 1, Length[Lvals]}
];

ratioScalar = N[{scalarErr[[1]]/scalarErr[[2]], scalarErr[[2]]/scalarErr[[3]]}, 30];
ratioTensor = N[{tensorErr[[1]]/tensorErr[[2]], tensorErr[[2]]/tensorErr[[3]]}, 30];
ratioVector = N[{vectorErr[[1]]/vectorErr[[2]], vectorErr[[2]]/vectorErr[[3]]}, 30];

checkNumA = And @@ Join[
  Flatten[Map[# > 0 &, ksEigVals, {2}]],
  Flatten[Map[# > 0 &, gsEigVals, {2}]],
  Flatten[Map[# > 0 &, scalarRootVals, {3}]],
  Flatten[Map[# > 0 &, tensorRootVals, {2}]],
  Flatten[Map[# > 0 &, vectorRootVals, {2}]]
];

checkNumB = Max[uvErr] < 1*^-6;
checkNumC = Max[Abs[ratioScalar - {8, 8}]] < 6*^-1;
checkNumD = Max[Abs[ratioTensor - {8, 8}]] < 3*^-1;
checkNumE = Max[Abs[ratioVector - {8, 8}]] < 3*^-1;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF &&
   checkG && checkH && checkNumA && checkNumB && checkNumC && checkNumD && checkNumE
];

nbName = "13_Full_Spectrum_Dispersion_Residues.nb";
logName = "13_full_spectrum_dispersion_residues.log";
logLines = {
  "Notebook: " <> nbName,
  "poly = " <> ToString[poly, InputForm],
  "a2 = " <> ToString[a2, InputForm],
  "a2Expected = " <> ToString[a2Expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "a1 = " <> ToString[a1, InputForm],
  "a1Expected = " <> ToString[a1Expected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "a0 = " <> ToString[a0, InputForm],
  "a0Expected = " <> ToString[a0Expected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "xPlus = " <> ToString[xPlus, InputForm],
  "xMinus = " <> ToString[xMinus, InputForm],
  "checkD(sum rule) = " <> ToString[checkD, InputForm],
  "checkE(product rule) = " <> ToString[checkE, InputForm],
  "checkF1(dec sum) = " <> ToString[checkF1, InputForm],
  "checkF2(dec prod) = " <> ToString[checkF2, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "xTensor = " <> ToString[xTensor, InputForm],
  "xVector = " <> ToString[xVector, InputForm],
  "cScalarUV = " <> ToString[cScalarUV, InputForm],
  "checkG(channel separation) = " <> ToString[checkG, InputForm],
  "checkH1(dec UV sum) = " <> ToString[checkH1, InputForm],
  "checkH2(dec UV product) = " <> ToString[checkH2, InputForm],
  "checkH(dec UV speeds) = " <> ToString[checkH, InputForm],
  "ksEigVals = " <> ToString[ksEigVals, InputForm],
  "gsEigVals = " <> ToString[gsEigVals, InputForm],
  "scalarRootVals = " <> ToString[scalarRootVals, InputForm],
  "tensorRootVals = " <> ToString[tensorRootVals, InputForm],
  "vectorRootVals = " <> ToString[vectorRootVals, InputForm],
  "scalarUVFromRoots = " <> ToString[scalarUVFromRoots, InputForm],
  "scalarUVEigs = " <> ToString[scalarUVEigs, InputForm],
  "uvErr = " <> ToString[uvErr, InputForm],
  "scalarErr = " <> ToString[scalarErr, InputForm],
  "tensorErr = " <> ToString[tensorErr, InputForm],
  "vectorErr = " <> ToString[vectorErr, InputForm],
  "ratioScalar = " <> ToString[ratioScalar, InputForm],
  "ratioTensor = " <> ToString[ratioTensor, InputForm],
  "ratioVector = " <> ToString[ratioVector, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "checkNumE = " <> ToString[checkNumE, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "13_full_spectrum_dispersion_residues",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
