ClearAll["Global`*"];

nb = Notebook[{
  Cell["14 - Full 4D Covariant Gauge-Fixed Coefficients", "Title"],
  Cell["Direct 4D derivation map from the non-local action variation to quadratic metric+Theta coefficients, with fully covariant gauge-fixing shifts and IR checks.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2 == (d2I1 - q d2I0)/I0 - 2 dI0 (dI1 - q dI0)/I0^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {c2NL, c1NL, csNL, cwNL, cmNL} == Coefficients[q2, {ht^2, hv^2, xs^2, xw^2, 2 xs xw}]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    cA == cA_loc + cA_NL + DeltaCgfA
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Op4D == c2 P2 + c1 P1 + cs P0s + cw P0w + cm (P0sw + P0ws)
  ], "Input"]
}];

dI0 = at ht + av hv + as xs + aw xw;
dI1 = bt ht + bv hv + bs xs + bw xw;
d2I0 = A2 ht^2 + A1 hv^2 + As xs^2 + Aw xw^2 + 2 Am xs xw;
d2I1 = B2 ht^2 + B1 hv^2 + Bs xs^2 + Bw xw^2 + 2 Bm xs xw;

q2Expr = FullSimplify[(d2I1 - q d2I0)/i0 - 2 dI0 (dI1 - q dI0)/i0^2];

c2NL = FullSimplify[Coefficient[q2Expr, ht^2]];
c1NL = FullSimplify[Coefficient[q2Expr, hv^2]];
csNL = FullSimplify[Coefficient[q2Expr, xs^2]];
cwNL = FullSimplify[Coefficient[q2Expr, xw^2]];
cmNL = FullSimplify[Coefficient[q2Expr, xs xw]/2];

c2Expected = FullSimplify[(B2 - q A2)/i0 - 2 at (bt - q at)/i0^2];
c1Expected = FullSimplify[(B1 - q A1)/i0 - 2 av (bv - q av)/i0^2];
csExpected = FullSimplify[(Bs - q As)/i0 - 2 as (bs - q as)/i0^2];
cwExpected = FullSimplify[(Bw - q Aw)/i0 - 2 aw (bw - q aw)/i0^2];
cmExpected = FullSimplify[(Bm - q Am)/i0 - (as (bw - q aw) + aw (bs - q as))/i0^2];

checkA = FullSimplify[c2NL == c2Expected];
checkB = FullSimplify[c1NL == c1Expected];
checkC = FullSimplify[csNL == csExpected];
checkD = FullSimplify[cwNL == cwExpected];
checkE = FullSimplify[cmNL == cmExpected];

qBranch = q -> k0^2;
checkF = FullSimplify[(c2NL /. qBranch) == ((B2 - k0^2 A2)/i0 - 2 at (bt - k0^2 at)/i0^2)];
checkG = FullSimplify[(c1NL /. qBranch) == ((B1 - k0^2 A1)/i0 - 2 av (bv - k0^2 av)/i0^2)];
checkH = FullSimplify[(csNL /. qBranch) == ((Bs - k0^2 As)/i0 - 2 as (bs - k0^2 as)/i0^2)];
checkI = FullSimplify[(cwNL /. qBranch) == ((Bw - k0^2 Aw)/i0 - 2 aw (bw - k0^2 aw)/i0^2)];
checkJ = FullSimplify[(cmNL /. qBranch) == ((Bm - k0^2 Am)/i0 - (as (bw - k0^2 aw) + aw (bs - k0^2 as))/i0^2)];

(* local + gauge-fixing shifts *)
c2Tot = l2 + c2NL + g2;
c1Tot = l1 + c1NL + g1;
csTot = ls + csNL + gs;
cwTot = lw + cwNL + gw;
cmTot = lm + cmNL + gm;

checkK = FullSimplify[
  c2Tot == l2 + g2 + c2NL &&
   c1Tot == l1 + g1 + c1NL &&
   csTot == ls + gs + csNL &&
   cwTot == lw + gw + cwNL &&
   cmTot == lm + gm + cmNL
];

KsTot = {{csTot, cmTot}, {cmTot, cwTot}};
detSTot = FullSimplify[Det[KsTot]];
checkL = FullSimplify[detSTot == csTot*cwTot - cmTot^2];

resT = FullSimplify[1/c2Tot];
resV = FullSimplify[1/c1Tot];
resS = FullSimplify[Inverse[KsTot]];
checkM = FullSimplify[Det[resS] == 1/detSTot];

c2Lead = FullSimplify[Limit[i0 c2NL, i0 -> Infinity]];
c1Lead = FullSimplify[Limit[i0 c1NL, i0 -> Infinity]];
csLead = FullSimplify[Limit[i0 csNL, i0 -> Infinity]];
cwLead = FullSimplify[Limit[i0 cwNL, i0 -> Infinity]];
cmLead = FullSimplify[Limit[i0 cmNL, i0 -> Infinity]];

checkN = FullSimplify[c2Lead == B2 - q A2];
checkO = FullSimplify[c1Lead == B1 - q A1];
checkP = FullSimplify[csLead == Bs - q As];
checkQ = FullSimplify[cwLead == Bw - q Aw];
checkR = FullSimplify[cmLead == Bm - q Am];

c2TotInf = FullSimplify[Limit[c2Tot, i0 -> Infinity]];
c1TotInf = FullSimplify[Limit[c1Tot, i0 -> Infinity]];
ksTotInf = FullSimplify[Limit[KsTot, i0 -> Infinity]];
checkS = FullSimplify[c2TotInf == l2 + g2];
checkT = FullSimplify[c1TotInf == l1 + g1];
checkU = FullSimplify[ksTotInf == {{ls + gs, lm + gm}, {lm + gm, lw + gw}}];

numRule = {
  q -> 3/5,
  at -> 1/8, av -> 1/10, as -> 1/10, aw -> 3/25,
  bt -> 1/5, bv -> 11/50, bs -> 1/4, bw -> 7/25,
  A2 -> 1/2, A1 -> 2/5, As -> 9/20, Aw -> 2/5, Am -> 1/6,
  B2 -> 4/5, B1 -> 3/4, Bs -> 7/10, Bw -> 11/15, Bm -> 2/5,
  l2 -> 1, l1 -> 9/10, ls -> 6/5, lw -> 7/5, lm -> 1/10,
  g2 -> 1/20, g1 -> -1/25, gs -> 1/50, gw -> -1/100, gm -> 1/200,
  k0 -> 2/3
};

Lvals = {10, 20, 40};

c2NLVals = Table[N[c2NL /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
c1NLVals = Table[N[c1NL /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
csNLVals = Table[N[csNL /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
cwNLVals = Table[N[cwNL /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
cmNLVals = Table[N[cmNL /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];

ratioC2NL = N[{Abs[c2NLVals[[1]]/c2NLVals[[2]]], Abs[c2NLVals[[2]]/c2NLVals[[3]]]}, 30];
ratioC1NL = N[{Abs[c1NLVals[[1]]/c1NLVals[[2]]], Abs[c1NLVals[[2]]/c1NLVals[[3]]]}, 30];
ratioCsNL = N[{Abs[csNLVals[[1]]/csNLVals[[2]]], Abs[csNLVals[[2]]/csNLVals[[3]]]}, 30];
ratioCwNL = N[{Abs[cwNLVals[[1]]/cwNLVals[[2]]], Abs[cwNLVals[[2]]/cwNLVals[[3]]]}, 30];
ratioCmNL = N[{Abs[cmNLVals[[1]]/cmNLVals[[2]]], Abs[cmNLVals[[2]]/cmNLVals[[3]]]}, 30];

c2TotVals = Table[N[c2Tot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
c1TotVals = Table[N[c1Tot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
ksTotVals = Table[N[KsTot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
detTotVals = Table[N[detSTot /. numRule /. i0 -> LL^3, 30], {LL, Lvals}];
ksEigVals = Table[Sort[N[Eigenvalues[ksTotVals[[ii]]], 30]], {ii, 1, Length[Lvals]}];

c2InfNum = N[c2TotInf /. numRule, 30];
c1InfNum = N[c1TotInf /. numRule, 30];
detInfNum = N[(Det[ksTotInf] /. numRule), 30];

c2Err = Table[N[Abs[c2TotVals[[ii]] - c2InfNum], 30], {ii, 1, Length[Lvals]}];
c1Err = Table[N[Abs[c1TotVals[[ii]] - c1InfNum], 30], {ii, 1, Length[Lvals]}];
detErr = Table[N[Abs[detTotVals[[ii]] - detInfNum], 30], {ii, 1, Length[Lvals]}];

ratioC2Err = N[{c2Err[[1]]/c2Err[[2]], c2Err[[2]]/c2Err[[3]]}, 30];
ratioC1Err = N[{c1Err[[1]]/c1Err[[2]], c1Err[[2]]/c1Err[[3]]}, 30];
ratioDetErr = N[{detErr[[1]]/detErr[[2]], detErr[[2]]/detErr[[3]]}, 30];

decompErr = Table[
  N[Max[Abs[{
      (c2Tot - (l2 + g2 + c2NL)),
      (c1Tot - (l1 + g1 + c1NL)),
      (csTot - (ls + gs + csNL)),
      (cwTot - (lw + gw + cwNL)),
      (cmTot - (lm + gm + cmNL))
    } /. numRule /. i0 -> LL^3]], 30],
  {LL, Lvals}
];

checkNumA = Max[Abs[ratioC2NL - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioC1NL - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioCsNL - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioCwNL - {8, 8}]] < 2*^-1 &&
  Max[Abs[ratioCmNL - {8, 8}]] < 3*^-1;

checkNumB = And @@ Join[
   Map[# > 0 &, c2TotVals],
   Map[# > 0 &, c1TotVals],
   Map[# > 0 &, detTotVals],
   Flatten[Map[# > 0 &, ksEigVals, {2}]]
];

checkNumC = Max[Abs[ratioC2Err - {8, 8}]] < 3*^-1 &&
  Max[Abs[ratioC1Err - {8, 8}]] < 3*^-1 &&
  Max[Abs[ratioDetErr - {8, 8}]] < 4*^-1;

checkNumD = Max[decompErr] < 10^-25;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE &&
   checkF && checkG && checkH && checkI && checkJ &&
   checkK && checkL && checkM &&
   checkN && checkO && checkP && checkQ && checkR &&
   checkS && checkT && checkU &&
   checkNumA && checkNumB && checkNumC && checkNumD
];

nbName = "14_Full4D_Covariant_GaugeFixed_Coefficients.nb";
logName = "14_full4d_covariant_gauge_fixed_coefficients.log";
logLines = {
  "Notebook: " <> nbName,
  "q2Expr = " <> ToString[q2Expr, InputForm],
  "c2NL = " <> ToString[c2NL, InputForm],
  "c1NL = " <> ToString[c1NL, InputForm],
  "csNL = " <> ToString[csNL, InputForm],
  "cwNL = " <> ToString[cwNL, InputForm],
  "cmNL = " <> ToString[cmNL, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "checkE = " <> ToString[checkE, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "checkG = " <> ToString[checkG, InputForm],
  "checkH = " <> ToString[checkH, InputForm],
  "checkI = " <> ToString[checkI, InputForm],
  "checkJ = " <> ToString[checkJ, InputForm],
  "c2Tot = " <> ToString[c2Tot, InputForm],
  "c1Tot = " <> ToString[c1Tot, InputForm],
  "csTot = " <> ToString[csTot, InputForm],
  "cwTot = " <> ToString[cwTot, InputForm],
  "cmTot = " <> ToString[cmTot, InputForm],
  "detSTot = " <> ToString[detSTot, InputForm],
  "checkK = " <> ToString[checkK, InputForm],
  "checkL = " <> ToString[checkL, InputForm],
  "checkM = " <> ToString[checkM, InputForm],
  "c2Lead = " <> ToString[c2Lead, InputForm],
  "c1Lead = " <> ToString[c1Lead, InputForm],
  "csLead = " <> ToString[csLead, InputForm],
  "cwLead = " <> ToString[cwLead, InputForm],
  "cmLead = " <> ToString[cmLead, InputForm],
  "checkN = " <> ToString[checkN, InputForm],
  "checkO = " <> ToString[checkO, InputForm],
  "checkP = " <> ToString[checkP, InputForm],
  "checkQ = " <> ToString[checkQ, InputForm],
  "checkR = " <> ToString[checkR, InputForm],
  "c2TotInf = " <> ToString[c2TotInf, InputForm],
  "c1TotInf = " <> ToString[c1TotInf, InputForm],
  "ksTotInf = " <> ToString[ksTotInf, InputForm],
  "checkS = " <> ToString[checkS, InputForm],
  "checkT = " <> ToString[checkT, InputForm],
  "checkU = " <> ToString[checkU, InputForm],
  "c2NLVals(L={10,20,40}) = " <> ToString[c2NLVals, InputForm],
  "c1NLVals(L={10,20,40}) = " <> ToString[c1NLVals, InputForm],
  "csNLVals(L={10,20,40}) = " <> ToString[csNLVals, InputForm],
  "cwNLVals(L={10,20,40}) = " <> ToString[cwNLVals, InputForm],
  "cmNLVals(L={10,20,40}) = " <> ToString[cmNLVals, InputForm],
  "ratioC2NL = " <> ToString[ratioC2NL, InputForm],
  "ratioC1NL = " <> ToString[ratioC1NL, InputForm],
  "ratioCsNL = " <> ToString[ratioCsNL, InputForm],
  "ratioCwNL = " <> ToString[ratioCwNL, InputForm],
  "ratioCmNL = " <> ToString[ratioCmNL, InputForm],
  "c2TotVals = " <> ToString[c2TotVals, InputForm],
  "c1TotVals = " <> ToString[c1TotVals, InputForm],
  "detTotVals = " <> ToString[detTotVals, InputForm],
  "ksEigVals = " <> ToString[ksEigVals, InputForm],
  "ratioC2Err = " <> ToString[ratioC2Err, InputForm],
  "ratioC1Err = " <> ToString[ratioC1Err, InputForm],
  "ratioDetErr = " <> ToString[ratioDetErr, InputForm],
  "decompErr = " <> ToString[decompErr, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "14_full4d_covariant_gauge_fixed_coefficients",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
