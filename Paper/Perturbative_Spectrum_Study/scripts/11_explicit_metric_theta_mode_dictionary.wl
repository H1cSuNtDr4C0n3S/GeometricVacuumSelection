ClearAll["Global`*"];

nb = Notebook[{
  Cell["11 - Explicit Metric+Theta Mode Dictionary", "Title"],
  Cell["Extraction of {m,c,p,m2,c2,k2} from an explicit ADM-inspired perturbative ansatz for gamma_ij, K_ij, and Theta.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    gamma == a^2 Exp[2 eps xs] DiagonalMatrix[{Exp[2 eps ht], Exp[-eps ht], Exp[-eps ht]}]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Kcov == H0 gamma + eps hv S && ThetaDot == 1 + eps xw
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    mu == Sqrt[Det[gamma]] Sqrt[X] && X == ThetaDot^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    chiProxy == Exp[3 eps (-xs + xw)]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    k2inv == gammaInv.gammaInv.Kcov.Kcov
  ], "Input"]
}];

phi = {ht, hv, xs, xw};
n = Length[phi];

halfHessian[expr_] := Table[
  FullSimplify[D[expr, phi[[i]], phi[[j]]]/2],
  {i, 1, n}, {j, 1, n}
];

vecFromLinear[expr_] := Table[FullSimplify[Coefficient[expr, phi[[i]]]], {i, 1, n}];
matEqQ[m1_, m2_] := TrueQ[And @@ Flatten@Table[FullSimplify[m1[[i, j]] == m2[[i, j]]], {i, 1, Length[m1]}, {j, 1, Length[m1[[1]]]}]];
vecEqQ[v1_, v2_] := TrueQ[And @@ Table[FullSimplify[v1[[i]] == v2[[i]]], {i, 1, Length[v1]}]];

(* Explicit ADM-inspired ansatz *)
gamma = a^2 Exp[2 eps xs] DiagonalMatrix[{Exp[2 eps ht], Exp[-eps ht], Exp[-eps ht]}];
gammaInv = FullSimplify[Inverse[gamma]];

S = {{0, sv, 0}, {sv, 0, 0}, {0, 0, 0}};
Kcov = FullSimplify[H0 gamma + eps hv S];

thetaDot = 1 + eps xw;
sqrtGamma = FullSimplify[Sqrt[Det[gamma]]];
sqrtX = thetaDot; (* positive-branch around eps=0 *)
muExpr = FullSimplify[sqrtGamma sqrtX];

chiExpr = FullSimplify[Exp[3 eps (-xs + xw)]];

K2Expr = FullSimplify[
  Sum[gammaInv[[i, k]] gammaInv[[j, l]] Kcov[[i, j]] Kcov[[k, l]], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}, {l, 1, 3}]
];

(* Build dictionaries *)
mu0 = FullSimplify[muExpr /. eps -> 0];
mu1lin = FullSimplify[D[muExpr, eps] /. eps -> 0];
mu2sec = FullSimplify[D[muExpr, {eps, 2}] /. eps -> 0];
mVec = vecFromLinear[mu1lin];
M2 = halfHessian[mu2sec];

chi0 = FullSimplify[chiExpr /. eps -> 0];
chi1lin = FullSimplify[D[chiExpr, eps] /. eps -> 0];
chi2sec = FullSimplify[D[chiExpr, {eps, 2}] /. eps -> 0];
cVec = vecFromLinear[chi1lin];
C2 = halfHessian[chi2sec];

k0 = FullSimplify[K2Expr /. eps -> 0];
k1lin = FullSimplify[D[K2Expr, eps] /. eps -> 0];
k2sec = FullSimplify[D[K2Expr, {eps, 2}] /. eps -> 0];
pVec = vecFromLinear[k1lin];
K2 = halfHessian[k2sec];

(* Reconstruction checks up to O(eps^2) *)
muRec = FullSimplify[mu0 + eps (mVec.phi) + eps^2/2 (phi.M2.phi)];
chiRec = FullSimplify[chi0 + eps (cVec.phi) + eps^2/2 (phi.C2.phi)];
kRec = FullSimplify[k0 + eps (pVec.phi) + eps^2/2 (phi.K2.phi)];

muSeries2 = Normal[Series[muExpr, {eps, 0, 2}]];
chiSeries2 = Normal[Series[chiExpr, {eps, 0, 2}]];
kSeries2 = Normal[Series[K2Expr, {eps, 0, 2}]];

checkA = FullSimplify[Expand[muRec - muSeries2] == 0];
checkB = FullSimplify[Expand[chiRec - chiSeries2] == 0];
checkC = FullSimplify[Expand[kRec - kSeries2] == 0];

(* Structural checks *)
checkD = vecEqQ[mVec, {0, 0, 3 mu0, mu0}];
checkE = vecEqQ[cVec, {0, 0, -3, 3}];
checkF = FullSimplify[k0 == 3 H0^2];

(* K1/K2 channel consistency *)
checkG = FullSimplify[Coefficient[k1lin, hv] == 0];
checkH = FullSimplify[Coefficient[k2sec, hv^2] > 0 /. {a -> 2, sv -> 1, H0 -> 1/3}];

(* Numerical truncation residual ~ O(eps^3) *)
numRule = {
  a -> 2,
  H0 -> 1/3,
  sv -> 1/5,
  ht -> 1/20,
  hv -> 1/25,
  xs -> 1/30,
  xw -> 1/35
};

muErr[e_] := N[Abs[(muExpr - muRec) /. numRule /. eps -> e], 40];
chiErr[e_] := N[Abs[(chiExpr - chiRec) /. numRule /. eps -> e], 40];
kErr[e_] := N[Abs[(K2Expr - kRec) /. numRule /. eps -> e], 40];

eVals = {1/20, 1/40, 1/80};
muErrVals = Table[muErr[e], {e, eVals}];
chiErrVals = Table[chiErr[e], {e, eVals}];
kErrVals = Table[kErr[e], {e, eVals}];

muRat = N[{muErrVals[[1]]/muErrVals[[2]], muErrVals[[2]]/muErrVals[[3]]}, 30];
chiRat = N[{chiErrVals[[1]]/chiErrVals[[2]], chiErrVals[[2]]/chiErrVals[[3]]}, 30];
kRat = N[{kErrVals[[1]]/kErrVals[[2]], kErrVals[[2]]/kErrVals[[3]]}, 30];

checkNumA = And @@ Join[Map[# > 0 &, muErrVals], Map[# > 0 &, chiErrVals], Map[# > 0 &, kErrVals]];
checkNumB = Max[Abs[muRat - {8, 8}]] < 1*^-1;
checkNumC = Max[Abs[chiRat - {8, 8}]] < 1*^-1;
checkNumD = Max[Abs[kRat - {8, 8}]] < 2*^-1;

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF &&
   checkG && checkH && checkNumA && checkNumB && checkNumC && checkNumD
];

nbName = "11_Explicit_Metric_Theta_Mode_Dictionary.nb";
logName = "11_explicit_metric_theta_mode_dictionary.log";
logLines = {
  "Notebook: " <> nbName,
  "mu0 = " <> ToString[mu0, InputForm],
  "mVec = " <> ToString[mVec, InputForm],
  "M2 = " <> ToString[M2, InputForm],
  "chi0 = " <> ToString[chi0, InputForm],
  "cVec = " <> ToString[cVec, InputForm],
  "C2 = " <> ToString[C2, InputForm],
  "k0 = " <> ToString[k0, InputForm],
  "pVec = " <> ToString[pVec, InputForm],
  "K2 = " <> ToString[K2, InputForm],
  "checkA(mu reconstruction) = " <> ToString[checkA, InputForm],
  "checkB(chi reconstruction) = " <> ToString[checkB, InputForm],
  "checkC(k reconstruction) = " <> ToString[checkC, InputForm],
  "checkD(mVec explicit) = " <> ToString[checkD, InputForm],
  "checkE(cVec explicit) = " <> ToString[checkE, InputForm],
  "checkF(k0=3H0^2) = " <> ToString[checkF, InputForm],
  "checkG(no linear hv in K2) = " <> ToString[checkG, InputForm],
  "checkH(positive hv^2 channel sample) = " <> ToString[checkH, InputForm],
  "muErrVals = " <> ToString[muErrVals, InputForm],
  "chiErrVals = " <> ToString[chiErrVals, InputForm],
  "kErrVals = " <> ToString[kErrVals, InputForm],
  "muRat = " <> ToString[muRat, InputForm],
  "chiRat = " <> ToString[chiRat, InputForm],
  "kRat = " <> ToString[kRat, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "11_explicit_metric_theta_mode_dictionary",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
