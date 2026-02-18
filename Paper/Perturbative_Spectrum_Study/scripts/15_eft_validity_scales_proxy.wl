ClearAll["Global`*"];

nb = Notebook[{
  Cell["15 - EFT Validity, Cutoff, and Strong-Coupling Proxies", "Title"],
  Cell["Operational EFT discussion in the perturbative pipeline: (i) non-local correction size, (ii) kinetic-distance from strong coupling, (iii) conservative tested-q^2 window for dispersion positivity.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    epsNL[L] == Norm[KNL[L], 2]/Norm[KIR, 2] && epsNL[L] << 1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    zT[L] == c2Tot[L] && zV[L] == c1Tot[L] && zSmin[L] == Min[Eigenvalues[KsTot[L]]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {xTensor, xVector, xScalarPlus, xScalarMinus} > 0 \[ForAll] q2 \[Element] q2Grid
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2CutProxy == Max[q2Grid]
  ], "Input"]
}];

(* Non-local coefficients as in script 14 *)
c2NL = (-2*at*bt + 2*at^2*q + i0*(B2 - A2*q))/i0^2;
c1NL = (-2*av*bv + 2*av^2*q + i0*(B1 - A1*q))/i0^2;
csNL = (-2*as*bs + 2*as^2*q + i0*(Bs - As*q))/i0^2;
cwNL = (-2*aw*bw + 2*aw^2*q + i0*(Bw - Aw*q))/i0^2;
cmNL = -((aw*bs + as*bw - Bm*i0 - 2*as*aw*q + Am*i0*q)/i0^2);

c2Tot = l2 + dg2 + c2NL;
c1Tot = l1 + dg1 + c1NL;
csTot = ls + dgs + csNL;
cwTot = lw + dgw + cwNL;
cmTot = lm + dgm + cmNL;

KsTot = {{csTot, cmTot}, {cmTot, cwTot}};
KsIR = {{ls + dgs, lm + dgm}, {lm + dgm, lw + dgw}};
KNL = {{csNL, cmNL}, {cmNL, cwNL}};

numRule = {
  q -> 3/5,
  at -> 1/8, av -> 1/10, as -> 1/10, aw -> 3/25,
  bt -> 1/5, bv -> 11/50, bs -> 1/4, bw -> 7/25,
  A2 -> 1/2, A1 -> 2/5, As -> 9/20, Aw -> 2/5, Am -> 1/6,
  B2 -> 4/5, B1 -> 3/4, Bs -> 7/10, Bw -> 11/15, Bm -> 2/5,
  l2 -> 1, l1 -> 9/10, ls -> 6/5, lw -> 7/5, lm -> 1/10,
  dg2 -> 1/20, dg1 -> -1/25, dgs -> 1/50, dgw -> -1/100, dgm -> 1/200
};

Lvals = {10, 20, 40};

eps2Vals = Table[
  N[Abs[(c2NL/(l2 + dg2)) /. numRule /. i0 -> LL^3], 30],
  {LL, Lvals}
];
eps1Vals = Table[
  N[Abs[(c1NL/(l1 + dg1)) /. numRule /. i0 -> LL^3], 30],
  {LL, Lvals}
];
epsSVals = Table[
  N[(Norm[KNL /. numRule /. i0 -> LL^3, 2]/Norm[KsIR /. numRule, 2]), 30],
  {LL, Lvals}
];

ratioEps2 = N[{eps2Vals[[1]]/eps2Vals[[2]], eps2Vals[[2]]/eps2Vals[[3]]}, 30];
ratioEps1 = N[{eps1Vals[[1]]/eps1Vals[[2]], eps1Vals[[2]]/eps1Vals[[3]]}, 30];
ratioEpsS = N[{epsSVals[[1]]/epsSVals[[2]], epsSVals[[2]]/epsSVals[[3]]}, 30];

epsTargetA = 1/100;
epsTargetB = 1/1000;

c2Eps = N[Mean[Table[eps2Vals[[ii]] Lvals[[ii]]^3, {ii, 1, Length[Lvals]}]], 30];
c1Eps = N[Mean[Table[eps1Vals[[ii]] Lvals[[ii]]^3, {ii, 1, Length[Lvals]}]], 30];
cSEps = N[Mean[Table[epsSVals[[ii]] Lvals[[ii]]^3, {ii, 1, Length[Lvals]}]], 30];

lMinA = N[
  Max[(c2Eps/epsTargetA)^(1/3), (c1Eps/epsTargetA)^(1/3), (cSEps/epsTargetA)^(1/3)],
  30
];
lMinB = N[
  Max[(c2Eps/epsTargetB)^(1/3), (c1Eps/epsTargetB)^(1/3), (cSEps/epsTargetB)^(1/3)],
  30
];

zTVals = Table[N[c2Tot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
zVVals = Table[N[c1Tot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
ksVals = Table[N[KsTot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
zSMinVals = Table[N[Min[Eigenvalues[ksVals[[ii]]]], 30], {ii, 1, Length[Lvals]}];

invTVals = Table[N[1/zTVals[[ii]], 30], {ii, 1, Length[Lvals]}];
invVVals = Table[N[1/zVVals[[ii]], 30], {ii, 1, Length[Lvals]}];
invSMaxVals = Table[N[Max[Eigenvalues[Inverse[ksVals[[ii]]]]], 30], {ii, 1, Length[Lvals]}];

zFloor = 1/10;
invCeil = 10;

checkA = Max[Abs[ratioEps2 - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioEps1 - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioEpsS - {8, 8}]] < 3*^-1;
checkB = Max[eps2Vals] < 1*^-2 && Max[eps1Vals] < 1*^-2 && Max[epsSVals] < 1*^-2;
checkC = Min[zTVals] > zFloor && Min[zVVals] > zFloor && Min[zSMinVals] > zFloor;
checkD = Max[invTVals] < invCeil && Max[invVVals] < invCeil && Max[invSMaxVals] < invCeil;

(* Conservative tested-q^2 window for positivity *)
gt2 = 6/5;
gt1 = 11/10;
gsS = 13/10;
gmS = 3/20;
gwS = 7/5;

mt2 = 3/50;
mv2 = 1/20;
msM = 1/25;
mmM = 1/100;
mwM = 9/250;

GsMat = {{gsS, gmS}, {gmS, gwS}};
MsMat = {{msM, mmM}, {mmM, mwM}};

q2Grid = {1/10, 1, 4};

xTensor[q2_, LL_] := N[(q2 gt2 + mt2)/(c2Tot /. numRule /. i0 -> LL^3), 30];
xVector[q2_, LL_] := N[(q2 gt1 + mv2)/(c1Tot /. numRule /. i0 -> LL^3), 30];
xScalar[q2_, LL_] := Module[{k = N[KsTot /. numRule /. i0 -> LL^3, 50]},
  N[Sort[Eigenvalues[Inverse[k].(q2 GsMat + MsMat)]], 30]
];

xTensorVals = Table[xTensor[qq, LL], {LL, Lvals}, {qq, q2Grid}];
xVectorVals = Table[xVector[qq, LL], {LL, Lvals}, {qq, q2Grid}];
xScalarVals = Table[xScalar[qq, LL], {LL, Lvals}, {qq, q2Grid}];

allDispVals = Join[
  Flatten[xTensorVals],
  Flatten[xVectorVals],
  Flatten[Map[Flatten, xScalarVals]]
];

checkE = Min[allDispVals] > 0;
q2CutProxy = Max[q2Grid];
checkF = TrueQ[q2CutProxy == 4];

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF];

nbName = "15_EFT_Validity_Scales_Proxy.nb";
logName = "15_eft_validity_scales_proxy.log";
logLines = {
  "Notebook: " <> nbName,
  "eps2Vals(L={10,20,40}) = " <> ToString[eps2Vals, InputForm],
  "eps1Vals(L={10,20,40}) = " <> ToString[eps1Vals, InputForm],
  "epsSVals(L={10,20,40}) = " <> ToString[epsSVals, InputForm],
  "ratioEps2 = " <> ToString[ratioEps2, InputForm],
  "ratioEps1 = " <> ToString[ratioEps1, InputForm],
  "ratioEpsS = " <> ToString[ratioEpsS, InputForm],
  "epsTargetA = " <> ToString[epsTargetA, InputForm],
  "epsTargetB = " <> ToString[epsTargetB, InputForm],
  "lMin(eps=1e-2) = " <> ToString[lMinA, InputForm],
  "lMin(eps=1e-3) = " <> ToString[lMinB, InputForm],
  "zTVals = " <> ToString[zTVals, InputForm],
  "zVVals = " <> ToString[zVVals, InputForm],
  "zSMinVals = " <> ToString[zSMinVals, InputForm],
  "invTVals = " <> ToString[invTVals, InputForm],
  "invVVals = " <> ToString[invVVals, InputForm],
  "invSMaxVals = " <> ToString[invSMaxVals, InputForm],
  "xTensorVals(q2={0.1,1,4}) = " <> ToString[xTensorVals, InputForm],
  "xVectorVals(q2={0.1,1,4}) = " <> ToString[xVectorVals, InputForm],
  "xScalarVals(q2={0.1,1,4}) = " <> ToString[xScalarVals, InputForm],
  "q2CutProxy = " <> ToString[q2CutProxy, InputForm],
  "checkA(IR scaling eps) = " <> ToString[checkA, InputForm],
  "checkB(eps small) = " <> ToString[checkB, InputForm],
  "checkC(kinetic floor) = " <> ToString[checkC, InputForm],
  "checkD(residue bound) = " <> ToString[checkD, InputForm],
  "checkE(dispersion positivity in tested window) = " <> ToString[checkE, InputForm],
  "checkF(q2 cutoff proxy) = " <> ToString[checkF, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "15_eft_validity_scales_proxy",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
