ClearAll["Global`*"];

nb = Notebook[{
  Cell["09 - Vector Sector and Full Stability Map", "Title"],
  Cell["Complete stability map with explicit vector channel and scalar diagonalization from K^-1 G, including IR-suppressed non-local corrections.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    cT2 == g2/c2 && cV2 == g1/c1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Ks == {{cs, cm}, {cm, cw}} && Gs == {{gs, gm}, {gm, gw}}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    eigScalar == Eigenvalues[Inverse[Ks].Gs]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    coeffs -> coeffs0 + ncoeffs/I0
  ], "Input"]
}];

Ks = {{cs, cm}, {cm, cw}};
Gs = {{gs, gm}, {gm, gw}};

detKs = FullSimplify[Det[Ks]];
traceSG = FullSimplify[Tr[Inverse[Ks].Gs]];
detSG = FullSimplify[Det[Inverse[Ks].Gs]];

traceExpected = FullSimplify[(gw cs - 2 gm cm + gs cw)/(cs cw - cm^2)];
detExpected = FullSimplify[(gs gw - gm^2)/(cs cw - cm^2)];

checkA = FullSimplify[detKs == cs cw - cm^2];
checkB = FullSimplify[traceSG == traceExpected];
checkC = FullSimplify[detSG == detExpected];

disc = FullSimplify[traceSG^2 - 4 detSG];
sPlus = FullSimplify[(traceSG + Sqrt[disc])/2];
sMinus = FullSimplify[(traceSG - Sqrt[disc])/2];

checkD = FullSimplify[sPlus + sMinus == traceSG];
checkE = FullSimplify[sPlus sMinus == detSG];

decRule = {cm -> 0, gm -> 0};
checkF = FullSimplify[(traceSG /. decRule) == gs/cs + gw/cw];
checkG = FullSimplify[(detSG /. decRule) == gs gw/(cs cw)];

cT2 = FullSimplify[g2/c2];
cV2 = FullSimplify[g1/c1];
checkH = FreeQ[cT2, {cs, cw, cm, gs, gw, gm}] && FreeQ[cV2, {cs, cw, cm, gs, gw, gm}];

(* IR-suppressed parametrization *)
c2IR = c20 + nc2/i0;
c1IR = c10 + nc1/i0;
csIR = cs0 + ncs/i0;
cwIR = cw0 + ncw/i0;
cmIR = cm0 + ncm/i0;

g2IR = g20 + ng2/i0;
g1IR = g10 + ng1/i0;
gsIR = gs0 + ngs/i0;
gwIR = gw0 + ngw/i0;
gmIR = gm0 + ngm/i0;

KsIR = {{csIR, cmIR}, {cmIR, cwIR}};
GsIR = {{gsIR, gmIR}, {gmIR, gwIR}};

cT2IR = FullSimplify[g2IR/c2IR];
cV2IR = FullSimplify[g1IR/c1IR];
scalarIR = FullSimplify[Eigenvalues[Inverse[KsIR].GsIR]];

cT2Inf = FullSimplify[Limit[cT2IR, i0 -> Infinity]];
cV2Inf = FullSimplify[Limit[cV2IR, i0 -> Infinity]];
scalarInf = FullSimplify[Limit[scalarIR, i0 -> Infinity]];
scalarLoc = FullSimplify[Eigenvalues[Inverse[{{cs0, cm0}, {cm0, cw0}}].{{gs0, gm0}, {gm0, gw0}}]];

checkI = FullSimplify[cT2Inf == g20/c20];
checkJ = FullSimplify[cV2Inf == g10/c10];
checkK = FullSimplify[Total[scalarInf] == Total[scalarLoc] && Times @@ scalarInf == Times @@ scalarLoc];

numRule = {
  c20 -> 1, c10 -> 9/10, cs0 -> 4/5, cw0 -> 6/5, cm0 -> 1/10,
  nc2 -> 1/3, nc1 -> 1/4, ncs -> 1/2, ncw -> -1/5, ncm -> 1/25,
  g20 -> 1, g10 -> 19/20, gs0 -> 7/10, gw0 -> 9/10, gm0 -> 1/20,
  ng2 -> -1/6, ng1 -> 1/10, ngs -> 1/5, ngw -> 1/8, ngm -> 1/30
};

cT2Vals = Table[N[cT2IR /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
cV2Vals = Table[N[cV2IR /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];

scalarVals = Table[
  Sort[N[scalarIR /. numRule /. i0 -> LL^3, 30]],
  {LL, {10, 20, 40}}
];
scalarLocNum = Sort[N[scalarLoc /. numRule, 30]];

detKsVals = Table[N[Det[KsIR /. numRule /. i0 -> LL^3], 30], {LL, {10, 20, 40}}];
detGsVals = Table[N[Det[GsIR /. numRule /. i0 -> LL^3], 30], {LL, {10, 20, 40}}];

cT2Err = Table[N[Abs[cT2Vals[[ii]] - (g20/c20 /. numRule)], 30], {ii, 1, 3}];
cV2Err = Table[N[Abs[cV2Vals[[ii]] - (g10/c10 /. numRule)], 30], {ii, 1, 3}];
scalarErr = Table[
  N[Max[Abs[scalarVals[[ii]] - scalarLocNum]], 30],
  {ii, 1, 3}
];

ratioT = N[{cT2Err[[1]]/cT2Err[[2]], cT2Err[[2]]/cT2Err[[3]]}, 30];
ratioV = N[{cV2Err[[1]]/cV2Err[[2]], cV2Err[[2]]/cV2Err[[3]]}, 30];
ratioS = N[{scalarErr[[1]]/scalarErr[[2]], scalarErr[[2]]/scalarErr[[3]]}, 30];

checkNumA = And @@ Join[
   Map[# > 0 &, cT2Vals],
   Map[# > 0 &, cV2Vals],
   Flatten[Map[# > 0 &, scalarVals, {2}]],
   Map[# > 0 &, detKsVals],
   Map[# > 0 &, detGsVals]
];
checkNumB = Max[Abs[ratioT - {8, 8}]] < 2*^-1;
checkNumC = Max[Abs[ratioV - {8, 8}]] < 2*^-1;
checkNumD = Max[Abs[ratioS - {8, 8}]] < 5*^-1;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF && checkG &&
   checkH && checkI && checkJ && checkK &&
   checkNumA && checkNumB && checkNumC && checkNumD
];

nbName = "09_Vector_Sector_Full_Stability_Map.nb";
logName = "09_vector_sector_full_stability_map.log";
logLines = {
  "Notebook: " <> nbName,
  "detKs = " <> ToString[detKs, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "traceSG = " <> ToString[traceSG, InputForm],
  "traceExpected = " <> ToString[traceExpected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "detSG = " <> ToString[detSG, InputForm],
  "detExpected = " <> ToString[detExpected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "sPlus = " <> ToString[sPlus, InputForm],
  "sMinus = " <> ToString[sMinus, InputForm],
  "checkD(sum rule) = " <> ToString[checkD, InputForm],
  "checkE(product rule) = " <> ToString[checkE, InputForm],
  "checkF(dec trace) = " <> ToString[checkF, InputForm],
  "checkG(dec det) = " <> ToString[checkG, InputForm],
  "cT2 = " <> ToString[cT2, InputForm],
  "cV2 = " <> ToString[cV2, InputForm],
  "checkH(channel separation) = " <> ToString[checkH, InputForm],
  "cT2Inf = " <> ToString[cT2Inf, InputForm],
  "cV2Inf = " <> ToString[cV2Inf, InputForm],
  "scalarInf = " <> ToString[scalarInf, InputForm],
  "scalarLoc = " <> ToString[scalarLoc, InputForm],
  "checkI = " <> ToString[checkI, InputForm],
  "checkJ = " <> ToString[checkJ, InputForm],
  "checkK = " <> ToString[checkK, InputForm],
  "cT2Vals(L={10,20,40}) = " <> ToString[cT2Vals, InputForm],
  "cV2Vals(L={10,20,40}) = " <> ToString[cV2Vals, InputForm],
  "scalarVals = " <> ToString[scalarVals, InputForm],
  "scalarLocNum = " <> ToString[scalarLocNum, InputForm],
  "detKsVals = " <> ToString[detKsVals, InputForm],
  "detGsVals = " <> ToString[detGsVals, InputForm],
  "cT2Err = " <> ToString[cT2Err, InputForm],
  "cV2Err = " <> ToString[cV2Err, InputForm],
  "scalarErr = " <> ToString[scalarErr, InputForm],
  "ratioT = " <> ToString[ratioT, InputForm],
  "ratioV = " <> ToString[ratioV, InputForm],
  "ratioS = " <> ToString[ratioS, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "09_vector_sector_full_stability_map",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
