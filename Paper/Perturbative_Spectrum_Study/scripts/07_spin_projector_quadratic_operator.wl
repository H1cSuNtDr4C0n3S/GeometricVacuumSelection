ClearAll["Global`*"];

nb = Notebook[{
  Cell["07 - Spin-Projector Quadratic Operator (4D Proxy)", "Title"],
  Cell["General tensor/vector/scalar block decomposition for a symmetric spatial tensor kernel using Barnes-Rivers projectors; includes IR-suppressed non-local scalar mixing.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Op == c2 P2 + c1 P1 + cs P0s + cw P0w + cm (P0sw + P0ws)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    P2 + P1 + P0s + P0w == Isym
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Kscalar == {{cs, cm}, {cm, cw}} && detKscalar == cs*cw - cm^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {cs, cw, cm} -> {cs0 + ns/I0, cw0 + nw/I0, cm0 + nm/I0}
  ], "Input"]
}];

idx = Range[3];
zero4 = ConstantArray[0, {3, 3, 3, 3}];

delta = Table[KroneckerDelta[i, j], {i, idx}, {j, idx}];
omega = Table[KroneckerDelta[i, 3] KroneckerDelta[j, 3], {i, idx}, {j, idx}];
theta = delta - omega;

P2 = Table[
  (theta[[i, k]] theta[[j, l]] + theta[[i, l]] theta[[j, k]])/2 -
    theta[[i, j]] theta[[k, l]]/2,
  {i, idx}, {j, idx}, {k, idx}, {l, idx}
];

P1 = Table[
  (theta[[i, k]] omega[[j, l]] + theta[[i, l]] omega[[j, k]] +
     theta[[j, k]] omega[[i, l]] + theta[[j, l]] omega[[i, k]])/2,
  {i, idx}, {j, idx}, {k, idx}, {l, idx}
];

P0s = Table[theta[[i, j]] theta[[k, l]]/2, {i, idx}, {j, idx}, {k, idx}, {l, idx}];
P0w = Table[omega[[i, j]] omega[[k, l]], {i, idx}, {j, idx}, {k, idx}, {l, idx}];
P0sw = Table[theta[[i, j]] omega[[k, l]]/Sqrt[2], {i, idx}, {j, idx}, {k, idx}, {l, idx}];
P0ws = Table[omega[[i, j]] theta[[k, l]]/Sqrt[2], {i, idx}, {j, idx}, {k, idx}, {l, idx}];

Isym = Table[
  (KroneckerDelta[i, k] KroneckerDelta[j, l] + KroneckerDelta[i, l] KroneckerDelta[j, k])/2,
  {i, idx}, {j, idx}, {k, idx}, {l, idx}
];

ComposeOp[A_, B_] := Table[
  Sum[A[[i, j, m, n]] B[[m, n, k, l]], {m, idx}, {n, idx}],
  {i, idx}, {j, idx}, {k, idx}, {l, idx}
];

TensorEqQ[A_, B_] := TrueQ[FullSimplify[A == B]];

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

Op = c2 P2 + c1 P1 + cs P0s + cw P0w + cm (P0sw + P0ws);

checkK = TensorEqQ[ComposeOp[P2, ComposeOp[Op, P0s]], zero4] &&
  TensorEqQ[ComposeOp[P2, ComposeOp[Op, P0w]], zero4] &&
  TensorEqQ[ComposeOp[P2, ComposeOp[Op, P0sw]], zero4];

checkL = TensorEqQ[ComposeOp[P2, ComposeOp[Op, P2]], c2 P2];

Kscalar = {{cs, cm}, {cm, cw}};
traceKscalar = FullSimplify[Tr[Kscalar]];
detKscalar = FullSimplify[Det[Kscalar]];
checkM = TrueQ[FullSimplify[detKscalar == cs*cw - cm^2]];

csIR = cs0 + ns/i0;
cwIR = cw0 + nw/i0;
cmIR = cm0 + nm/i0;
c2IR = c20 + n2/i0;

traceIR = FullSimplify[(csIR + cwIR)];
detIR = FullSimplify[csIR*cwIR - cmIR^2];

traceIRInf = FullSimplify[Limit[traceIR, i0 -> Infinity]];
detIRInf = FullSimplify[Limit[detIR, i0 -> Infinity]];
c2IRInf = FullSimplify[Limit[c2IR, i0 -> Infinity]];

checkN = TrueQ[FullSimplify[traceIRInf == cs0 + cw0]];
checkO = TrueQ[FullSimplify[detIRInf == cs0*cw0 - cm0^2]];
checkP = TrueQ[FullSimplify[c2IRInf == c20]];

numRule = {
  cs0 -> 4/5, cw0 -> 6/5, cm0 -> 1/10,
  ns -> 1/2, nw -> -1/5, nm -> 1/20,
  c20 -> 1, n2 -> 1/3
};

scalarVals = Table[
  N[{{csIR, cmIR}, {cmIR, cwIR}} /. numRule /. i0 -> LL^3, 30],
  {LL, {10, 20, 40}}
];

scalarEigVals = Table[N[Eigenvalues[scalarVals[[ii]]], 30], {ii, 1, 3}];
tensorVals = Table[N[c2IR /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];

scalarPositive = And @@ Flatten[Map[# > 0 &, scalarEigVals, {2}]];
tensorPositive = And @@ Map[# > 0 &, tensorVals];

scalarLim = N[{{cs0, cm0}, {cm0, cw0}} /. numRule, 30];
scalarErr = Table[
  N[Max[Abs[Flatten[scalarVals[[ii]] - scalarLim]]], 30],
  {ii, 1, 3}
];
scalarErrRatios = N[{scalarErr[[1]]/scalarErr[[2]], scalarErr[[2]]/scalarErr[[3]]}, 30];

tensorErr = Table[N[Abs[tensorVals[[ii]] - (c20 /. numRule)], 30], {ii, 1, 3}];
tensorErrRatios = N[{tensorErr[[1]]/tensorErr[[2]], tensorErr[[2]]/tensorErr[[3]]}, 30];

checkNumA = scalarPositive && tensorPositive;
checkNumB = Max[Abs[scalarErrRatios - {8, 8}]] < 1*^-2;
checkNumC = Max[Abs[tensorErrRatios - {8, 8}]] < 1*^-2;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF && checkG &&
   checkH && checkI && checkJ && checkK && checkL && checkM && checkN &&
   checkO && checkP && checkNumA && checkNumB && checkNumC
];

nbName = "07_Spin_Projector_Quadratic_Operator.nb";
logName = "07_spin_projector_quadratic_operator.log";
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
  "checkK(no tensor-scalar mixing in isotropic Op) = " <> ToString[checkK, InputForm],
  "checkL(P2 O P2 = c2 P2) = " <> ToString[checkL, InputForm],
  "traceKscalar = " <> ToString[traceKscalar, InputForm],
  "detKscalar = " <> ToString[detKscalar, InputForm],
  "checkM(det formula) = " <> ToString[checkM, InputForm],
  "traceIR = " <> ToString[traceIR, InputForm],
  "detIR = " <> ToString[detIR, InputForm],
  "c2IR = " <> ToString[c2IR, InputForm],
  "traceIRInf = " <> ToString[traceIRInf, InputForm],
  "detIRInf = " <> ToString[detIRInf, InputForm],
  "c2IRInf = " <> ToString[c2IRInf, InputForm],
  "checkN = " <> ToString[checkN, InputForm],
  "checkO = " <> ToString[checkO, InputForm],
  "checkP = " <> ToString[checkP, InputForm],
  "scalarVals(L={10,20,40}) = " <> ToString[scalarVals, InputForm],
  "scalarEigVals = " <> ToString[scalarEigVals, InputForm],
  "tensorVals(L={10,20,40}) = " <> ToString[tensorVals, InputForm],
  "scalarErr = " <> ToString[scalarErr, InputForm],
  "scalarErrRatios = " <> ToString[scalarErrRatios, InputForm],
  "tensorErr = " <> ToString[tensorErr, InputForm],
  "tensorErrRatios = " <> ToString[tensorErrRatios, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "07_spin_projector_quadratic_operator",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
