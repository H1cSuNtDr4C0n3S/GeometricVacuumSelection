ClearAll["Global`*"];

nb = Notebook[{
  Cell["12 - Full 4D Spin-Projector Algebra and Inversion", "Title"],
  Cell["Barnes-Rivers decomposition in full 4D (not only spatial proxy), with explicit inverse operator and IR-suppressed coefficient flow.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Op4D == c2 P2 + c1 P1 + cs P0s + cw P0w + cm (P0sw + P0ws)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Op4DInv == P2/c2 + P1/c1 + (cw P0s + cs P0w - cm (P0sw + P0ws))/(cs*cw - cm^2)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    P2 + P1 + P0s + P0w == Isym
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {c2, c1, cs, cw, cm} -> {c20 + nc2/I0, c10 + nc1/I0, cs0 + ncs/I0, cw0 + ncw/I0, cm0 + ncm/I0}
  ], "Input"]
}];

idx = Range[4];
d = 4;
dm1 = d - 1;
zero4 = ConstantArray[0, {4, 4, 4, 4}];

delta = Table[KroneckerDelta[mu, nu], {mu, idx}, {nu, idx}];
omega = Table[KroneckerDelta[mu, 4] KroneckerDelta[nu, 4], {mu, idx}, {nu, idx}];
theta = delta - omega;

P2 = Table[
  (theta[[mu, al]] theta[[nu, be]] + theta[[mu, be]] theta[[nu, al]])/2 -
    theta[[mu, nu]] theta[[al, be]]/dm1,
  {mu, idx}, {nu, idx}, {al, idx}, {be, idx}
];

P1 = Table[
  (theta[[mu, al]] omega[[nu, be]] + theta[[mu, be]] omega[[nu, al]] +
     theta[[nu, al]] omega[[mu, be]] + theta[[nu, be]] omega[[mu, al]])/2,
  {mu, idx}, {nu, idx}, {al, idx}, {be, idx}
];

P0s = Table[theta[[mu, nu]] theta[[al, be]]/dm1, {mu, idx}, {nu, idx}, {al, idx}, {be, idx}];
P0w = Table[omega[[mu, nu]] omega[[al, be]], {mu, idx}, {nu, idx}, {al, idx}, {be, idx}];
P0sw = Table[theta[[mu, nu]] omega[[al, be]]/Sqrt[dm1], {mu, idx}, {nu, idx}, {al, idx}, {be, idx}];
P0ws = Table[omega[[mu, nu]] theta[[al, be]]/Sqrt[dm1], {mu, idx}, {nu, idx}, {al, idx}, {be, idx}];

Isym = Table[
  (KroneckerDelta[mu, al] KroneckerDelta[nu, be] + KroneckerDelta[mu, be] KroneckerDelta[nu, al])/2,
  {mu, idx}, {nu, idx}, {al, idx}, {be, idx}
];

ComposeOp[A_, B_] := Table[
  Sum[A[[mu, nu, ro, si]] B[[ro, si, al, be]], {ro, idx}, {si, idx}],
  {mu, idx}, {nu, idx}, {al, idx}, {be, idx}
];

TensorEqQ[A_, B_, ass_: True] := TrueQ[Simplify[A == B, Assumptions -> ass]];

checkA = TensorEqQ[ComposeOp[P2, P2], P2];
checkB = TensorEqQ[ComposeOp[P1, P1], P1];
checkC = TensorEqQ[ComposeOp[P0s, P0s], P0s];
checkD = TensorEqQ[ComposeOp[P0w, P0w], P0w];
checkE = TensorEqQ[ComposeOp[P2, P1], zero4];
checkF = TensorEqQ[ComposeOp[P2, P0s], zero4];
checkG = TensorEqQ[ComposeOp[P1, P0s], zero4];
checkH = TensorEqQ[ComposeOp[P0sw, P0ws], P0s];
checkI = TensorEqQ[ComposeOp[P0ws, P0sw], P0w];
checkJ = TensorEqQ[P2 + P1 + P0s + P0w, Isym];

Op4D = c2 P2 + c1 P1 + cs P0s + cw P0w + cm (P0sw + P0ws);
detS = FullSimplify[cs*cw - cm^2];
Op4DInv = P2/c2 + P1/c1 + (cw P0s + cs P0w - cm (P0sw + P0ws))/detS;

invAss = c1 != 0 && c2 != 0 && detS != 0;
checkK = TensorEqQ[ComposeOp[Op4D, Op4DInv], Isym, invAss];
checkL = TensorEqQ[ComposeOp[Op4DInv, Op4D], Isym, invAss];

checkM = TensorEqQ[ComposeOp[P2, ComposeOp[Op4D, P0s]], zero4] &&
  TensorEqQ[ComposeOp[P2, ComposeOp[Op4D, P0w]], zero4] &&
  TensorEqQ[ComposeOp[P2, ComposeOp[Op4D, P0sw]], zero4];

c2IR = c20 + nc2/i0;
c1IR = c10 + nc1/i0;
csIR = cs0 + ncs/i0;
cwIR = cw0 + ncw/i0;
cmIR = cm0 + ncm/i0;

detSIR = FullSimplify[csIR*cwIR - cmIR^2];
resTIR = FullSimplify[1/c2IR];
resVIR = FullSimplify[1/c1IR];
KsIR = {{csIR, cmIR}, {cmIR, cwIR}};
resSIR = FullSimplify[Inverse[KsIR]];

checkN = FullSimplify[Det[resSIR] == 1/detSIR];

c2Inf = FullSimplify[Limit[c2IR, i0 -> Infinity]];
c1Inf = FullSimplify[Limit[c1IR, i0 -> Infinity]];
ksInf = FullSimplify[Limit[KsIR, i0 -> Infinity]];
resTInf = FullSimplify[Limit[resTIR, i0 -> Infinity]];
resVInf = FullSimplify[Limit[resVIR, i0 -> Infinity]];

checkO = FullSimplify[c2Inf == c20 && c1Inf == c10];
checkP = FullSimplify[ksInf == {{cs0, cm0}, {cm0, cw0}}];
checkQ = FullSimplify[resTInf == 1/c20 && resVInf == 1/c10];

numRule = {
  c20 -> 1,
  c10 -> 9/10,
  cs0 -> 6/5,
  cw0 -> 7/5,
  cm0 -> 1/10,
  nc2 -> 1/3,
  nc1 -> -1/5,
  ncs -> 1/2,
  ncw -> -1/4,
  ncm -> 1/30
};

Lvals = {10, 20, 40};

c2Vals = Table[N[c2IR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
c1Vals = Table[N[c1IR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
detVals = Table[N[detSIR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];

resTVals = Table[N[resTIR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
resVVals = Table[N[resVIR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];

ksVals = Table[N[KsIR /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
ksEigVals = Table[Sort[Eigenvalues[ksVals[[ii]]]], {ii, 1, Length[Lvals]}];
resSEigVals = Table[Sort[Eigenvalues[N[resSIR /. numRule /. i0 -> LL^3, 30]]], {LL, Lvals}];

c2Err = Table[N[Abs[c2Vals[[ii]] - (c20 /. numRule)], 30], {ii, 1, Length[Lvals]}];
c1Err = Table[N[Abs[c1Vals[[ii]] - (c10 /. numRule)], 30], {ii, 1, Length[Lvals]}];
detErr = Table[N[Abs[detVals[[ii]] - (cs0*cw0 - cm0^2 /. numRule)], 30], {ii, 1, Length[Lvals]}];

resTErr = Table[N[Abs[resTVals[[ii]] - (1/c20 /. numRule)], 30], {ii, 1, Length[Lvals]}];
resVErr = Table[N[Abs[resVVals[[ii]] - (1/c10 /. numRule)], 30], {ii, 1, Length[Lvals]}];

ratioC2 = N[{c2Err[[1]]/c2Err[[2]], c2Err[[2]]/c2Err[[3]]}, 30];
ratioC1 = N[{c1Err[[1]]/c1Err[[2]], c1Err[[2]]/c1Err[[3]]}, 30];
ratioDet = N[{detErr[[1]]/detErr[[2]], detErr[[2]]/detErr[[3]]}, 30];
ratioResT = N[{resTErr[[1]]/resTErr[[2]], resTErr[[2]]/resTErr[[3]]}, 30];
ratioResV = N[{resVErr[[1]]/resVErr[[2]], resVErr[[2]]/resVErr[[3]]}, 30];

checkNumA = And @@ Join[
  Map[# > 0 &, c2Vals],
  Map[# > 0 &, c1Vals],
  Map[# > 0 &, detVals],
  Flatten[Map[# > 0 &, ksEigVals, {2}]],
  Flatten[Map[# > 0 &, resSEigVals, {2}]],
  Map[# > 0 &, resTVals],
  Map[# > 0 &, resVVals]
];

checkNumB = Max[Abs[ratioC2 - {8, 8}]] < 2*^-1 && Max[Abs[ratioC1 - {8, 8}]] < 2*^-1;
checkNumC = Max[Abs[ratioDet - {8, 8}]] < 3*^-1;
checkNumD = Max[Abs[ratioResT - {8, 8}]] < 3*^-1 && Max[Abs[ratioResV - {8, 8}]] < 3*^-1;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF && checkG &&
   checkH && checkI && checkJ && checkK && checkL && checkM && checkN &&
   checkO && checkP && checkQ && checkNumA && checkNumB && checkNumC && checkNumD
];

nbName = "12_Full4D_Spin_Projector_Inversion.nb";
logName = "12_full4d_spin_projector_inversion.log";
logLines = {
  "Notebook: " <> nbName,
  "checkA(P2*P2=P2) = " <> ToString[checkA, InputForm],
  "checkB(P1*P1=P1) = " <> ToString[checkB, InputForm],
  "checkC(P0s*P0s=P0s) = " <> ToString[checkC, InputForm],
  "checkD(P0w*P0w=P0w) = " <> ToString[checkD, InputForm],
  "checkE(P2*P1=0) = " <> ToString[checkE, InputForm],
  "checkF(P2*P0s=0) = " <> ToString[checkF, InputForm],
  "checkG(P1*P0s=0) = " <> ToString[checkG, InputForm],
  "checkH(P0sw*P0ws=P0s) = " <> ToString[checkH, InputForm],
  "checkI(P0ws*P0sw=P0w) = " <> ToString[checkI, InputForm],
  "checkJ(completeness) = " <> ToString[checkJ, InputForm],
  "detS = " <> ToString[detS, InputForm],
  "checkK(Op*OpInv=I) = " <> ToString[checkK, InputForm],
  "checkL(OpInv*Op=I) = " <> ToString[checkL, InputForm],
  "checkM(no tensor-scalar leakage) = " <> ToString[checkM, InputForm],
  "detSIR = " <> ToString[detSIR, InputForm],
  "resSIR = " <> ToString[resSIR, InputForm],
  "checkN(det inverse) = " <> ToString[checkN, InputForm],
  "c2Inf = " <> ToString[c2Inf, InputForm],
  "c1Inf = " <> ToString[c1Inf, InputForm],
  "ksInf = " <> ToString[ksInf, InputForm],
  "resTInf = " <> ToString[resTInf, InputForm],
  "resVInf = " <> ToString[resVInf, InputForm],
  "checkO = " <> ToString[checkO, InputForm],
  "checkP = " <> ToString[checkP, InputForm],
  "checkQ = " <> ToString[checkQ, InputForm],
  "c2Vals(L={10,20,40}) = " <> ToString[c2Vals, InputForm],
  "c1Vals(L={10,20,40}) = " <> ToString[c1Vals, InputForm],
  "detVals = " <> ToString[detVals, InputForm],
  "ksEigVals = " <> ToString[ksEigVals, InputForm],
  "resSEigVals = " <> ToString[resSEigVals, InputForm],
  "resTVals = " <> ToString[resTVals, InputForm],
  "resVVals = " <> ToString[resVVals, InputForm],
  "ratioC2 = " <> ToString[ratioC2, InputForm],
  "ratioC1 = " <> ToString[ratioC1, InputForm],
  "ratioDet = " <> ToString[ratioDet, InputForm],
  "ratioResT = " <> ToString[ratioResT, InputForm],
  "ratioResV = " <> ToString[ratioResV, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "12_full4d_spin_projector_inversion",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
