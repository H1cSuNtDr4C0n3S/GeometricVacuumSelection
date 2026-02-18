ClearAll["Global`*"];

nb = Notebook[{
  Cell["08 - Covariant to Spin-Coefficient Map", "Title"],
  Cell["Coefficient identification from the covariant second variation of Q=I1/I0 to spin blocks {c2,c1,cs,cw,cm}.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2 == (d2I1 - q d2I0)/I0 - 2 dI0 (dI1 - q dI0)/I0^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dI0 == at ht + av hv + as xs + aw xw && dI1 == bt ht + bv hv + bs xs + bw xw
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    d2I0 == A2 ht^2 + A1 hv^2 + As xs^2 + Aw xw^2 + 2 Am xs xw
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    d2I1 == B2 ht^2 + B1 hv^2 + Bs xs^2 + Bw xw^2 + 2 Bm xs xw
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q2 == c2 ht^2 + c1 hv^2 + cs xs^2 + cw xw^2 + 2 cm xs xw
  ], "Input"]
}];

dI0 = at ht + av hv + as xs + aw xw;
dI1 = bt ht + bv hv + bs xs + bw xw;
d2I0 = A2 ht^2 + A1 hv^2 + As xs^2 + Aw xw^2 + 2 Am xs xw;
d2I1 = B2 ht^2 + B1 hv^2 + Bs xs^2 + Bw xw^2 + 2 Bm xs xw;

q2Expr = FullSimplify[(d2I1 - q d2I0)/i0 - 2 dI0 (dI1 - q dI0)/i0^2];

c2 = FullSimplify[Coefficient[q2Expr, ht^2]];
c1 = FullSimplify[Coefficient[q2Expr, hv^2]];
cs = FullSimplify[Coefficient[q2Expr, xs^2]];
cw = FullSimplify[Coefficient[q2Expr, xw^2]];
cm = FullSimplify[Coefficient[q2Expr, xs xw]/2];

c2Expected = FullSimplify[(B2 - q A2)/i0 - 2 at (bt - q at)/i0^2];
c1Expected = FullSimplify[(B1 - q A1)/i0 - 2 av (bv - q av)/i0^2];
csExpected = FullSimplify[(Bs - q As)/i0 - 2 as (bs - q as)/i0^2];
cwExpected = FullSimplify[(Bw - q Aw)/i0 - 2 aw (bw - q aw)/i0^2];
cmExpected = FullSimplify[(Bm - q Am)/i0 - (as (bw - q aw) + aw (bs - q as))/i0^2];

checkA = FullSimplify[c2 == c2Expected];
checkB = FullSimplify[c1 == c1Expected];
checkC = FullSimplify[cs == csExpected];
checkD = FullSimplify[cw == cwExpected];
checkE = FullSimplify[cm == cmExpected];

Kscalar = {{cs, cm}, {cm, cw}};
detKscalar = FullSimplify[Det[Kscalar]];
detExpected = FullSimplify[cs cw - cm^2];
checkF = FullSimplify[detKscalar == detExpected];

c2Lead = FullSimplify[Limit[i0 c2, i0 -> Infinity]];
c1Lead = FullSimplify[Limit[i0 c1, i0 -> Infinity]];
csLead = FullSimplify[Limit[i0 cs, i0 -> Infinity]];
cwLead = FullSimplify[Limit[i0 cw, i0 -> Infinity]];
cmLead = FullSimplify[Limit[i0 cm, i0 -> Infinity]];

checkG = FullSimplify[c2Lead == B2 - q A2];
checkH = FullSimplify[c1Lead == B1 - q A1];
checkI = FullSimplify[csLead == Bs - q As];
checkJ = FullSimplify[cwLead == Bw - q Aw];
checkK = FullSimplify[cmLead == Bm - q Am];

numRule = {
  q -> 3/5,
  at -> 1/8, av -> 1/10, as -> 1/10, aw -> 3/25,
  bt -> 1/5, bv -> 11/50, bs -> 1/4, bw -> 7/25,
  A2 -> 1/2, A1 -> 2/5, As -> 9/20, Aw -> 2/5, Am -> 1/6,
  B2 -> 4/5, B1 -> 3/4, Bs -> 7/10, Bw -> 11/15, Bm -> 2/5
};

c2Vals = Table[N[c2 /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
c1Vals = Table[N[c1 /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
csVals = Table[N[cs /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
cwVals = Table[N[cw /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
cmVals = Table[N[cm /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];

Kvals = Table[N[Kscalar /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
eigs = Table[N[Eigenvalues[Kvals[[ii]]], 30], {ii, 1, 3}];

ratioC2 = N[{c2Vals[[1]]/c2Vals[[2]], c2Vals[[2]]/c2Vals[[3]]}, 30];
ratioC1 = N[{c1Vals[[1]]/c1Vals[[2]], c1Vals[[2]]/c1Vals[[3]]}, 30];
ratioCs = N[{csVals[[1]]/csVals[[2]], csVals[[2]]/csVals[[3]]}, 30];
ratioCm = N[{cmVals[[1]]/cmVals[[2]], cmVals[[2]]/cmVals[[3]]}, 30];

checkNumA = And @@ Join[
   Map[# > 0 &, c2Vals],
   Map[# > 0 &, c1Vals],
   Flatten[Map[# > 0 &, eigs, {2}]]
];
checkNumB = Max[Abs[ratioC2 - {8, 8}]] < 5*^-2 && Max[Abs[ratioC1 - {8, 8}]] < 5*^-2;
checkNumC = Max[Abs[ratioCs - {8, 8}]] < 8*^-2 && Max[Abs[ratioCm - {8, 8}]] < 8*^-2;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF &&
   checkG && checkH && checkI && checkJ && checkK &&
   checkNumA && checkNumB && checkNumC
];

nbName = "08_Covariant_To_Spin_Coefficients_Map.nb";
logName = "08_covariant_to_spin_coefficients_map.log";
logLines = {
  "Notebook: " <> nbName,
  "q2Expr = " <> ToString[q2Expr, InputForm],
  "c2 = " <> ToString[c2, InputForm],
  "c2Expected = " <> ToString[c2Expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "c1 = " <> ToString[c1, InputForm],
  "c1Expected = " <> ToString[c1Expected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "cs = " <> ToString[cs, InputForm],
  "csExpected = " <> ToString[csExpected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "cw = " <> ToString[cw, InputForm],
  "cwExpected = " <> ToString[cwExpected, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "cm = " <> ToString[cm, InputForm],
  "cmExpected = " <> ToString[cmExpected, InputForm],
  "checkE = " <> ToString[checkE, InputForm],
  "detKscalar = " <> ToString[detKscalar, InputForm],
  "detExpected = " <> ToString[detExpected, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "c2Lead = " <> ToString[c2Lead, InputForm],
  "c1Lead = " <> ToString[c1Lead, InputForm],
  "csLead = " <> ToString[csLead, InputForm],
  "cwLead = " <> ToString[cwLead, InputForm],
  "cmLead = " <> ToString[cmLead, InputForm],
  "checkG = " <> ToString[checkG, InputForm],
  "checkH = " <> ToString[checkH, InputForm],
  "checkI = " <> ToString[checkI, InputForm],
  "checkJ = " <> ToString[checkJ, InputForm],
  "checkK = " <> ToString[checkK, InputForm],
  "c2Vals(L={10,20,40}) = " <> ToString[c2Vals, InputForm],
  "c1Vals(L={10,20,40}) = " <> ToString[c1Vals, InputForm],
  "csVals(L={10,20,40}) = " <> ToString[csVals, InputForm],
  "cwVals(L={10,20,40}) = " <> ToString[cwVals, InputForm],
  "cmVals(L={10,20,40}) = " <> ToString[cmVals, InputForm],
  "Kvals = " <> ToString[Kvals, InputForm],
  "eigs = " <> ToString[eigs, InputForm],
  "ratioC2 = " <> ToString[ratioC2, InputForm],
  "ratioC1 = " <> ToString[ratioC1, InputForm],
  "ratioCs = " <> ToString[ratioCs, InputForm],
  "ratioCm = " <> ToString[ratioCm, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "08_covariant_to_spin_coefficients_map",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
