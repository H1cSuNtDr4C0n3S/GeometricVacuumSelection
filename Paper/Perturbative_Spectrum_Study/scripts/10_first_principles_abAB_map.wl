ClearAll["Global`*"];

nb = Notebook[{
  Cell["10 - First-Principles abAB Map (Covariant Expansion)", "Title"],
  Cell["Derivation of first/second-order coefficient blocks {a,b,A,B} directly from covariant epsilon-expansion of j0=chi mu and j1=chi mu K2.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    muEps == mu + eps mu1 + eps^2 mu2/2 && chiEps == chi + eps chi1 + eps^2 chi2/2 && kEps == k0 + eps k1 + eps^2 k2/2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    j0eps == chiEps muEps && j1eps == chiEps muEps kEps
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dI0 == a.phi && dI1 == b.phi && d2I0 == phi.A.phi && d2I1 == phi.B.phi
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    qKernel == (B - q A)/I0 - (a (b-q a)^T + (b-q a) a^T)/I0^2
  ], "Input"]
}];

n = 4;
phi = {ht, hv, xs, xw};

symMat[prefix_] := Table[
  Symbol[prefix <> ToString[Min[i, j]] <> ToString[Max[i, j]]],
  {i, 1, n}, {j, 1, n}
];

mVec = Table[Symbol["m" <> ToString[i]], {i, 1, n}];
cVec = Table[Symbol["c" <> ToString[i]], {i, 1, n}];
pVec = Table[Symbol["p" <> ToString[i]], {i, 1, n}];

M2 = symMat["m2"];
C2 = symMat["c2"];
K2 = symMat["k2"];

mu1 = mVec.phi;
chi1 = cVec.phi;
k1 = pVec.phi;

mu2 = Sum[M2[[i, j]] phi[[i]] phi[[j]], {i, 1, n}, {j, 1, n}];
chi2 = Sum[C2[[i, j]] phi[[i]] phi[[j]], {i, 1, n}, {j, 1, n}];
k2 = Sum[K2[[i, j]] phi[[i]] phi[[j]], {i, 1, n}, {j, 1, n}];

muEps = mu + eps mu1 + eps^2 mu2/2;
chiEps = chi + eps chi1 + eps^2 chi2/2;
kEps = k0 + eps k1 + eps^2 k2/2;

j0eps = Expand[chiEps muEps];
j1eps = Expand[chiEps muEps kEps];

dj0 = FullSimplify[D[j0eps, eps] /. eps -> 0];
d2j0 = FullSimplify[D[j0eps, {eps, 2}] /. eps -> 0];
dj1 = FullSimplify[D[j1eps, eps] /. eps -> 0];
d2j1 = FullSimplify[D[j1eps, {eps, 2}] /. eps -> 0];

aVec = Table[FullSimplify[Coefficient[dj0, phi[[i]]]], {i, 1, n}];
bVec = Table[FullSimplify[Coefficient[dj1, phi[[i]]]], {i, 1, n}];

halfHessian[expr_] := Table[
  FullSimplify[D[expr, phi[[i]], phi[[j]]]/2],
  {i, 1, n}, {j, 1, n}
];

AKernel = halfHessian[d2j0];
BKernel = halfHessian[d2j1];

aExpected = FullSimplify[chi mVec + mu cVec];
bExpected = FullSimplify[chi k0 mVec + chi mu pVec + mu k0 cVec];

AExpected = FullSimplify[chi M2 + mu C2 + Outer[Times, mVec, cVec] + Outer[Times, cVec, mVec]];
BExpected = FullSimplify[
  chi k0 M2 + mu k0 C2 + chi mu K2 +
   k0 (Outer[Times, mVec, cVec] + Outer[Times, cVec, mVec]) +
   mu (Outer[Times, cVec, pVec] + Outer[Times, pVec, cVec]) +
   chi (Outer[Times, mVec, pVec] + Outer[Times, pVec, mVec])
];

vecEqQ[v1_, v2_] := TrueQ[And @@ Table[FullSimplify[v1[[i]] == v2[[i]]], {i, 1, Length[v1]}]];
matEqQ[m1_, m2_] := TrueQ[And @@ Flatten@Table[FullSimplify[m1[[i, j]] == m2[[i, j]]], {i, 1, Length[m1]}, {j, 1, Length[m1[[1]]]}]];

checkA = vecEqQ[aVec, aExpected];
checkB = vecEqQ[bVec, bExpected];
checkC = matEqQ[AKernel, AExpected];
checkD = matEqQ[BKernel, BExpected];

qKernelFromAB = FullSimplify[
  (BKernel - q AKernel)/i0 -
   (Outer[Times, aVec, bVec - q aVec] + Outer[Times, bVec - q aVec, aVec])/i0^2
];

q2Poly = FullSimplify[(d2j1 - q d2j0)/i0 - 2 dj0 (dj1 - q dj0)/i0^2];
qKernelFromPoly = halfHessian[q2Poly];

checkE = matEqQ[qKernelFromAB, qKernelFromPoly];

qKernelLead = FullSimplify[Limit[i0 qKernelFromAB, i0 -> Infinity]];
qKernelLeadExpected = FullSimplify[BKernel - q AKernel];
checkF = matEqQ[qKernelLead, qKernelLeadExpected];

mkRule[prefix_, fun_] := Flatten@Table[
  Symbol[prefix <> ToString[i] <> ToString[j]] -> fun[i, j],
  {i, 1, n}, {j, i, n}
];

numRule = Join[
  {
    mu -> 1,
    chi -> 3/2,
    k0 -> 2/3,
    q -> 3/5
  },
  Table[Symbol["m" <> ToString[i]] -> 1/(18 + 2 i), {i, 1, n}],
  Table[Symbol["c" <> ToString[i]] -> 1/(20 + 2 i), {i, 1, n}],
  Table[Symbol["p" <> ToString[i]] -> 1/(16 + 2 i), {i, 1, n}],
  mkRule["m2", Function[{i, j}, 1/(8 + i + j)]],
  mkRule["c2", Function[{i, j}, 1/(10 + i + j)]],
  mkRule["k2", Function[{i, j}, 1/(12 + i + j)]]
];

Lvals = {10, 20, 40};
kernelVals = Table[
  N[qKernelFromAB /. numRule /. i0 -> LL^3, 30],
  {LL, Lvals}
];
normVals = Table[N[Max[Abs[Flatten[kernelVals[[ii]]]]], 30], {ii, 1, 3}];
ratios = N[{normVals[[1]]/normVals[[2]], normVals[[2]]/normVals[[3]]}, 30];

leadNum = N[qKernelLead /. numRule, 30];
scaledVals = Table[N[(Lvals[[ii]]^3) kernelVals[[ii]], 30], {ii, 1, 3}];
leadErr = Table[N[Max[Abs[Flatten[scaledVals[[ii]] - leadNum]]], 30], {ii, 1, 3}];
leadRatios = N[{leadErr[[1]]/leadErr[[2]], leadErr[[2]]/leadErr[[3]]}, 30];

checkNumA = normVals[[1]] > normVals[[2]] > normVals[[3]] > 0;
checkNumB = Max[Abs[ratios - {8, 8}]] < 2*^-1;
checkNumC = Max[Abs[leadRatios - {8, 8}]] < 3*^-1;

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkNumA && checkNumB && checkNumC];

nbName = "10_First_Principles_abAB_Map.nb";
logName = "10_first_principles_abAB_map.log";
logLines = {
  "Notebook: " <> nbName,
  "aVec = " <> ToString[aVec, InputForm],
  "aExpected = " <> ToString[aExpected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "bVec = " <> ToString[bVec, InputForm],
  "bExpected = " <> ToString[bExpected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "AKernel = " <> ToString[AKernel, InputForm],
  "AExpected = " <> ToString[AExpected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "BKernel = " <> ToString[BKernel, InputForm],
  "BExpected = " <> ToString[BExpected, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "checkE(qKernel AB vs poly) = " <> ToString[checkE, InputForm],
  "qKernelLead = " <> ToString[qKernelLead, InputForm],
  "qKernelLeadExpected = " <> ToString[qKernelLeadExpected, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "kernelVals(L={10,20,40}) = " <> ToString[kernelVals, InputForm],
  "normVals = " <> ToString[normVals, InputForm],
  "ratios = " <> ToString[ratios, InputForm],
  "leadNum = " <> ToString[leadNum, InputForm],
  "scaledVals = " <> ToString[scaledVals, InputForm],
  "leadErr = " <> ToString[leadErr, InputForm],
  "leadRatios = " <> ToString[leadRatios, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "10_first_principles_abAB_map",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
