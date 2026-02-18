ClearAll["Global`*"];

nb = Notebook[{
  Cell["19 - FRW Linearized Full EOM Reduction", "Title"],
  Cell["Linearized system (Ea, EN, Elam) around de Sitter: eliminate constraint variables and obtain the reduced scalar-mode ODE.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    a = Exp[H0 t] (1 + eps u[t]), N = 1 + eps n[t], lambda = lam0 + eps l1[t]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Elam^(1) == 0 -> n[t] == u'[t]/H0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    EN^(1) == 0 -> l1[t] = l1[u,u']
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Ea^(1) -> reduced ODE for u[t]
  ], "Input"]
}];

t = Symbol["t"];
aT = aa[t];
nT = NN[t];
lamT = ll[t];

Lfrw = -3 aT D[aT, t]^2/(kappa nT) - (Lambda aT^3 nT)/kappa +
  lamT (3 aT D[aT, t]^2/nT^2 - K0^2 aT^3);

Ea = Simplify[D[Lfrw, aT] - D[D[Lfrw, D[aT, t]], t]];
EN = Simplify[D[Lfrw, nT]];
Elam = Simplify[D[Lfrw, lamT]];

aPert = Exp[H0 t] (1 + eps u[t]);
nPert = 1 + eps n[t];
lPert = lam0 + eps l1[t];

subPert = {
  aa[t] -> aPert,
  Derivative[1][aa][t] -> D[aPert, t],
  Derivative[2][aa][t] -> D[aPert, {t, 2}],
  NN[t] -> nPert,
  Derivative[1][NN][t] -> D[nPert, t],
  ll[t] -> lPert,
  Derivative[1][ll][t] -> D[lPert, t]
};

eaSeries = Normal[Series[Simplify[Ea /. subPert], {eps, 0, 1}]];
enSeries = Normal[Series[Simplify[EN /. subPert], {eps, 0, 1}]];
elSeries = Normal[Series[Simplify[Elam /. subPert], {eps, 0, 1}]];

ea1 = Simplify[Coefficient[eaSeries, eps, 1]];
en1 = Simplify[Coefficient[enSeries, eps, 1]];
el1 = Simplify[Coefficient[elSeries, eps, 1]];

branchRule = {
  H0 -> K0/Sqrt[3],
  lam0 -> (K0^2 - Lambda)/(2 K0^2 kappa)
};

ea1b = Simplify[ea1 /. branchRule];
en1b = Simplify[en1 /. branchRule];
el1b = Simplify[el1 /. branchRule];

nRule = First[Solve[el1b == 0, n[t]]];
expectedN = Simplify[(Sqrt[3] Derivative[1][u][t])/K0];
dnRule = Simplify[Derivative[1][n][t] -> D[n[t] /. nRule, t]];
checkA = Simplify[(n[t] /. nRule) == expectedN, Assumptions -> K0 != 0];

enReduced = Simplify[en1b /. nRule /. dnRule];
lRule = First[Solve[enReduced == 0, l1[t]]];
dlRule = Simplify[Derivative[1][l1][t] -> D[l1[t] /. lRule, t]];
checkB = Simplify[(enReduced /. lRule) == 0];

eaReduced = Simplify[ea1b /. nRule /. dnRule /. lRule /. dlRule];
u2Rules = Solve[eaReduced == 0, Derivative[2][u][t]];
checkC = Length[u2Rules] >= 1;
u2Rule = First[u2Rules];
checkD = Simplify[(eaReduced /. u2Rule) == 0];

u2Expr = Simplify[Derivative[2][u][t] /. u2Rule];

(* numeric algebraic consistency on random local jet values *)
eaJet = eaReduced /. {
  u[t] -> u0,
  Derivative[1][u][t] -> u1,
  Derivative[2][u][t] -> u2
};
u2Jet = u2Expr /. {
  u[t] -> u0,
  Derivative[1][u][t] -> u1
};

paramSub = {K0 -> 6/5, kappa -> 3/2, Lambda -> 1/10, H0 -> Sqrt[(6/5)^2/3]};
samples = {
  {u0 -> 1/10, u1 -> 1/20},
  {u0 -> -2/15, u1 -> 1/8},
  {u0 -> 3/25, u1 -> -1/12}
};
resVals = Table[
  N[eaJet /. u2 -> u2Jet /. paramSub /. ss, 30],
  {ss, samples}
];
maxAbs = N[Max[Abs[resVals]], 30];
checkNum = maxAbs < 10^-10;

check = TrueQ[checkA && checkB && checkC && checkD && checkNum];

nbName = "19_FRW_Linearized_Full_EOM_Reduction.nb";
logName = "19_frw_linearized_full_eom_reduction.log";
logLines = {
  "Notebook: " <> nbName,
  "ea1_branch = " <> ToString[ea1b, InputForm],
  "en1_branch = " <> ToString[en1b, InputForm],
  "el1_branch = " <> ToString[el1b, InputForm],
  "nRule = " <> ToString[nRule, InputForm],
  "dnRule = " <> ToString[dnRule, InputForm],
  "expectedN = " <> ToString[expectedN, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "lRule = " <> ToString[lRule, InputForm],
  "dlRule = " <> ToString[dlRule, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "eaReduced = " <> ToString[eaReduced, InputForm],
  "u2Rule = " <> ToString[u2Rule, InputForm],
  "u2Expr = " <> ToString[u2Expr, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "resVals = " <> ToString[resVals, InputForm],
  "maxAbs = " <> ToString[maxAbs, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "19_frw_linearized_full_eom_reduction",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxAbs" -> ToString[maxAbs, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
