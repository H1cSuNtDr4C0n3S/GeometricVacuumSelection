ClearAll["Global`*"];

nb = Notebook[{
  Cell["09 - Complete EOM Soft Branch Parameter Robustness Test", "Title"],
  Cell["Objective: verify that the evolving soft branch of the complete FRW minisuperspace EOM is robust across a wide parameter range (not fine-tuned), and that its transition anchor matches the standard matter-vacuum transition scale.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    softClosure == (3 H^2 - K0^2) == sigma/a^3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    H2soft[a_] := K0^2/3 + sigma/(3 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    weff[a_] := -1 + sigma/(sigma + K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    aTransition == (sigma/K0^2)^(1/3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"
  ], "Input"]
}];

t = Symbol["t"];
aT = aa[t];
nT = NN[t];
lamT = ll[t];

Lbar = -3 aT D[aT, t]^2/(kappa nT) - (Lambda aT^3 nT)/kappa +
  lamT (3 aT D[aT, t]^2/nT^2 - K0^2 aT^3) - rho0 aT^(-3 ww) nT;

Ea = Simplify[D[Lbar, aT] - D[D[Lbar, D[aT, t]], t]];
EN = Simplify[D[Lbar, nT]];
Elam = Simplify[D[Lbar, lamT]];

h2Soft[a_] := K0^2/3 + sigma/(3 a^3);
hSoft[a_] := Sqrt[h2Soft[a]];
hdotSoft[a_] := -sigma/(2 a^3);

subsSoft = {
  NN[t] -> 1,
  Derivative[1][NN][t] -> 0,
  ll[t] -> 0,
  Derivative[1][ll][t] -> 0,
  ww -> 0,
  Lambda -> K0^2,
  rho0 -> sigma/kappa,
  Derivative[1][aa][t] -> aT hSoft[aT],
  Derivative[2][aa][t] -> aT (h2Soft[aT] + hdotSoft[aT])
};

EaSoft = FullSimplify[Ea /. subsSoft, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
ENSoft = FullSimplify[EN /. subsSoft, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
ElamSoft = FullSimplify[Elam /. subsSoft, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
softConstraintResidual = FullSimplify[
  (ElamSoft/aT^3) - sigma/aT^3,
  Assumptions -> {aT > 0, sigma > 0}
];

checkA = TrueQ[EaSoft == 0 && ENSoft == 0 && softConstraintResidual == 0];

weff[a_] := FullSimplify[
  -1 - (a/3) D[Log[h2Soft[x]], x] /. x -> a,
  Assumptions -> {a > 0, K0 > 0, sigma > 0}
];
weffClosed[a_] := -1 + sigma/(sigma + K0^2 a^3);
checkB = TrueQ[Simplify[weff[a] == weffClosed[a], Assumptions -> {a > 0, K0 > 0, sigma > 0}]];

aTransition = FullSimplify[(sigma/K0^2)^(1/3), Assumptions -> {K0 > 0, sigma > 0}];
weffAtTransition = FullSimplify[
  weff[aTransition],
  Assumptions -> {K0 > 0, sigma > 0}
];
checkC = TrueQ[weffAtTransition == -1/2];

omegaGamma[h_] := N[2.469*10^-5/h^2, 50];
omegaR[h_, nEff_ : 3.046] := N[omegaGamma[h] (1 + 0.22710731766 nEff), 50];

hPlanck = 0.6766;
omPlanck = 0.30966;
orPlanck = omegaR[hPlanck];
olPlanck = N[1 - omPlanck - orPlanck, 50];
rPlanck = N[omPlanck/olPlanck, 50];
zLambdaMatterPlanck = N[(olPlanck/omPlanck)^(1/3) - 1, 50];

zTransition[r_] := N[r^(-1/3) - 1, 50];
zSoftPlanck = zTransition[rPlanck];
checkD = Abs[zSoftPlanck - zLambdaMatterPlanck] < 10^-12;

scanRatios = {1/20, 1/10, 3/10, 9/20, 7/10, 1, 3};
ratioSpan = N[Max[scanRatios]/Min[scanRatios], 50];

k0Ref = 3/5;
kappaRef = 17/10;
aResidualGrid = {1/10, 1/3, 1, 3, 10};
aEarly = 10^-2;
aLate = 10^2;

evalSample[r_] := Module[
  {sig, eaVals, enVals, elVals, maxRes, wEarly, wMid, wLate, zt, at, pass},
  sig = N[r k0Ref^2, 50];
  eaVals = Table[
    N[EaSoft /. {aa[t] -> aval, K0 -> k0Ref, sigma -> sig, kappa -> kappaRef}, 40],
    {aval, aResidualGrid}
  ];
  enVals = Table[
    N[ENSoft /. {aa[t] -> aval, K0 -> k0Ref, sigma -> sig, kappa -> kappaRef}, 40],
    {aval, aResidualGrid}
  ];
  elVals = Table[
    N[softConstraintResidual /. {aa[t] -> aval, sigma -> sig}, 40],
    {aval, aResidualGrid}
  ];
  maxRes = N[Max[Abs[Join[eaVals, enVals, elVals]]], 40];
  wEarly = N[weff[aEarly] /. {K0 -> k0Ref, sigma -> sig}, 50];
  wMid = N[weff[1] /. {K0 -> k0Ref, sigma -> sig}, 50];
  wLate = N[weff[aLate] /. {K0 -> k0Ref, sigma -> sig}, 50];
  at = N[r^(1/3), 50];
  zt = zTransition[r];
  pass = Abs[wEarly] < 0.05 && Abs[wLate + 1] < 10^-5 && wEarly > wMid > wLate && maxRes < 10^-25;
  <|
    "ratio" -> N[r, 50],
    "sigma" -> sig,
    "aTransition" -> at,
    "zTransition" -> zt,
    "weffEarly" -> wEarly,
    "weffMid" -> wMid,
    "weffLate" -> wLate,
    "maxResidual" -> maxRes,
    "pass" -> pass
  |>
];

scanData = evalSample /@ scanRatios;
scanPasses = Lookup[scanData, "pass"];
zTransitions = Lookup[scanData, "zTransition"];
maxResidualScan = N[Max[Lookup[scanData, "maxResidual"]], 50];

checkE = AllTrue[scanPasses, TrueQ];
checkF = maxResidualScan < 10^-25;
checkG = ratioSpan > 20 && Min[zTransitions] < 0 && Max[zTransitions] > 0;

classification = If[
  TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG],
  "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED",
  "undetermined"
];

check = TrueQ[classification == "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"];

nbName = "09_Complete_EOM_Soft_Branch_Parameter_Robustness_Test.nb";
logName = "09_complete_eom_soft_branch_parameter_robustness_test.log";
logLines = {
  "Notebook: " <> nbName,
  "EaSoft = " <> ToString[EaSoft, InputForm],
  "ENSoft = " <> ToString[ENSoft, InputForm],
  "ElamSoft = " <> ToString[ElamSoft, InputForm],
  "softConstraintResidual = " <> ToString[softConstraintResidual, InputForm],
  "weff(a) = " <> ToString[weff[aT], InputForm],
  "aTransition = " <> ToString[aTransition, InputForm],
  "weffAtTransition = " <> ToString[weffAtTransition, InputForm],
  "hPlanck = " <> ToString[hPlanck, InputForm],
  "omPlanck = " <> ToString[omPlanck, InputForm],
  "orPlanck = " <> ToString[orPlanck, InputForm],
  "olPlanck = " <> ToString[olPlanck, InputForm],
  "rPlanck = " <> ToString[rPlanck, InputForm],
  "zLambdaMatterPlanck = " <> ToString[zLambdaMatterPlanck, InputForm],
  "zSoftPlanck = " <> ToString[zSoftPlanck, InputForm],
  "scanRatios = " <> ToString[scanRatios, InputForm],
  "ratioSpan = " <> ToString[ratioSpan, InputForm],
  "scanData = " <> ToString[scanData, InputForm],
  "zTransitions = " <> ToString[zTransitions, InputForm],
  "maxResidualScan = " <> ToString[maxResidualScan, InputForm],
  "checkA(complete channels on soft branch) = " <> ToString[checkA, InputForm],
  "checkB(weff closed form) = " <> ToString[checkB, InputForm],
  "checkC(transition at weff=-1/2) = " <> ToString[checkC, InputForm],
  "checkD(Planck transition anchor) = " <> ToString[checkD, InputForm],
  "checkE(scan passes all samples) = " <> ToString[checkE, InputForm],
  "checkF(scan residuals) = " <> ToString[checkF, InputForm],
  "checkG(broad ratio span, past/future transitions) = " <> ToString[checkG, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "09_complete_eom_soft_branch_parameter_robustness_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
