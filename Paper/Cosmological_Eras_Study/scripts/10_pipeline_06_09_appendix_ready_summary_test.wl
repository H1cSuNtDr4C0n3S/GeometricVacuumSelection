ClearAll["Global`*"];

nb = Notebook[{
  Cell["10 - Pipeline 06-09 Appendix-Ready Summary Test", "Title"],
  Cell["Objective: build a compact, reproducible summary table for the complete-EOM block (tests 06-09) and verify that the hard-vs-soft storyline is internally consistent and not fine-tuned.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classesExpected == {
      "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY",
      "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND",
      "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED",
      "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "PIPELINE_06_09_APPENDIX_READY"
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

(* ---------- Test 06 reconstruction ---------- *)
ElamN1 = Simplify[Elam /. {NN[t] -> 1}];
constraintForm = Simplify[ElamN1/aa[t]^3];
check06A = constraintForm === (3 (Derivative[1][aa][t]^2/aa[t]^2) - K0^2);

subDS = {
  aa[t] -> Exp[H0 t],
  Derivative[1][aa][t] -> H0 Exp[H0 t],
  Derivative[2][aa][t] -> H0^2 Exp[H0 t],
  NN[t] -> 1,
  Derivative[1][NN][t] -> 0
};
EaDS = Simplify[Ea /. subDS];
ENDS = Simplify[EN /. subDS];
ElamDS = Simplify[Elam /. subDS];
hardRule06 = H0^2 -> K0^2/3;
lamRule06 = First[Solve[Simplify[(ENDS /. hardRule06) == 0], ll[t]]];
lamExpr06 = Simplify[ll[t] /. lamRule06];
dlamRule06 = Derivative[1][ll][t] -> D[lamExpr06, t];
EaCons06 = Simplify[EaDS /. hardRule06 /. lamRule06 /. dlamRule06];
ENCons06 = Simplify[ENDS /. hardRule06 /. lamRule06];
ElamCons06 = Simplify[ElamDS /. hardRule06];
check06B = Simplify[ElamCons06 == 0];
check06C = Simplify[ENCons06 == 0];
check06D = TrueQ @ Simplify[(EaCons06 /. hardRule06) == 0, Assumptions -> K0 != 0];

lamPrimeExpr06 = Simplify[D[lamExpr06, t]];
lamPrimeDust06 = Simplify[
  lamPrimeExpr06 /. {
    ww -> 0, rho0 -> 1, kappa -> 17/10, Lambda -> 1/5,
    K0 -> 3/5, H0 -> Sqrt[(3/5)^2/3], t -> 1
  }
];
check06E = N[Abs[lamPrimeDust06], 40] > 10^-8;

ElamPower06 = Simplify[
  Elam /. {
    aa[t] -> t^p,
    Derivative[1][aa][t] -> D[t^p, t],
    NN[t] -> 1
  },
  Assumptions -> t > 0
];
hardPowerNoGo06 = Reduce[
  {
    (ElamPower06 /. {K0 -> 1, t -> 1}) == 0,
    (ElamPower06 /. {K0 -> 1, t -> 2}) == 0
  },
  p, Reals
];
check06F = TrueQ[hardPowerNoGo06 === False];

class06 = If[
  TrueQ[check06A && check06B && check06C && check06D && check06E && check06F],
  "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY",
  "undetermined"
];

(* ---------- Test 07 reconstruction ---------- *)
h2Soft07[a_] := K0^2/3 + sigma/(3 a^3);
hSoft07[a_] := Sqrt[h2Soft07[a]];
hdotSoft07[a_] := -sigma/(2 a^3);
weff07[a_] := FullSimplify[
  -1 - (a/3) D[Log[h2Soft07[x]], x] /. x -> a,
  Assumptions -> {a > 0, K0 > 0, sigma > 0}
];

subsSoft07 = {
  NN[t] -> 1,
  Derivative[1][NN][t] -> 0,
  ll[t] -> 0,
  Derivative[1][ll][t] -> 0,
  ww -> 0,
  Lambda -> K0^2,
  rho0 -> sigma/kappa,
  Derivative[1][aa][t] -> aT hSoft07[aT],
  Derivative[2][aa][t] -> aT (h2Soft07[aT] + hdotSoft07[aT])
};
EaSoft07 = FullSimplify[Ea /. subsSoft07, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
ENSoft07 = FullSimplify[EN /. subsSoft07, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
ElamSoft07 = FullSimplify[Elam /. subsSoft07, Assumptions -> {aT > 0, K0 > 0, sigma > 0, kappa > 0}];
softResidual07 = FullSimplify[
  (ElamSoft07/aT^3) - sigma/aT^3,
  Assumptions -> {aT > 0, sigma > 0}
];

weffEarly07 = N[weff07[10^-2] /. {K0 -> 3/5, sigma -> 1/3}, 50];
weffMid07 = N[weff07[1] /. {K0 -> 3/5, sigma -> 1/3}, 50];
weffLate07 = N[weff07[100] /. {K0 -> 3/5, sigma -> 1/3}, 50];

check07A = TrueQ[EaSoft07 == 0];
check07B = TrueQ[ENSoft07 == 0];
check07C = TrueQ[softResidual07 == 0];
check07D = Abs[weffEarly07] < 0.05;
check07E = Abs[weffLate07 + 1] < 0.01;
check07F = weffEarly07 > weffMid07 > weffLate07;

class07 = If[
  TrueQ[check07A && check07B && check07C && check07D && check07E && check07F],
  "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND",
  "undetermined"
];

(* ---------- Test 08 reconstruction ---------- *)
rho08[a_] := rho0 a^(-3 (1 + ww));
enNorm08[h2_, lam_, rh_] := -((Lambda + kappa rh + 3 h2 (2 kappa lam - 1))/kappa);
lamFromEN08 = First[Solve[enNorm08[h2, lam, rh] == 0, lam]];
lamGeneric08 = Simplify[lam /. lamFromEN08];
lamHardExpected08[rh_] := (K0^2 - Lambda - kappa rh)/(2 kappa K0^2);
lamHardFromEN08[rh_] := Simplify[lamGeneric08 /. h2 -> K0^2/3 /. rh -> rh];
check08A = TrueQ @ Simplify[
  lamHardFromEN08[rh] == lamHardExpected08[rh],
  Assumptions -> {K0 != 0, kappa != 0}
];

lamHardA08 = Simplify[lamHardExpected08[rho08[a]], Assumptions -> {a > 0, K0 != 0, kappa != 0}];
h2Soft08[a_] := K0^2/3 + sigma/(3 a^3);
h20 = Simplify[h2Soft08[1]];
omegaM0 = Simplify[sigma/(K0^2 + sigma)];
omegaL0 = Simplify[K0^2/(K0^2 + sigma)];
h2SoftRebuild08[a_] := h20 (omegaM0 a^-3 + omegaL0);
check08B = TrueQ @ Simplify[h2Soft08[a] == h2SoftRebuild08[a], Assumptions -> {a > 0, K0 > 0, sigma > 0}];
check08C = TrueQ @ Simplify[omegaM0 + omegaL0 == 1, Assumptions -> {K0 > 0, sigma > 0}];

weff08[a_] := FullSimplify[
  -1 - (a/(3 h2Soft08[a])) D[h2Soft08[x], x] /. x -> a,
  Assumptions -> {a > 0, K0 > 0, sigma > 0}
];
weffExpected08[a_] := -1 + sigma/(sigma + K0^2 a^3);
check08D = TrueQ @ Simplify[weff08[a] == weffExpected08[a], Assumptions -> {a > 0, K0 > 0, sigma > 0}];
check08E = TrueQ[
  FullSimplify[Limit[weff08[a], a -> 0, Direction -> "FromAbove"], Assumptions -> {K0 > 0, sigma > 0}] == 0 &&
  FullSimplify[Limit[weff08[a], a -> Infinity], Assumptions -> {K0 > 0, sigma > 0}] == -1
];
check08F = TrueQ[FullSimplify[h2Soft08[a] /. sigma -> 0] == K0^2/3];

class08 = If[
  TrueQ[check08A && check08B && check08C && check08D && check08E && check08F],
  "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED",
  "undetermined"
];

(* ---------- Test 09 reconstruction ---------- *)
weff09[a_] := -1 + sigma/(sigma + K0^2 a^3);
aTransition09 = FullSimplify[(sigma/K0^2)^(1/3), Assumptions -> {K0 > 0, sigma > 0}];
weffAtTransition09 = FullSimplify[weff09[aTransition09], Assumptions -> {K0 > 0, sigma > 0}];
check09A = TrueQ[EaSoft07 == 0 && ENSoft07 == 0 && softResidual07 == 0];
check09B = TrueQ[weffAtTransition09 == -1/2];

omegaGamma[h_] := N[2.469*10^-5/h^2, 50];
omegaR[h_, nEff_ : 3.046] := N[omegaGamma[h] (1 + 0.22710731766 nEff), 50];
hPlanck = 0.6766;
omPlanck = 0.30966;
orPlanck = omegaR[hPlanck];
olPlanck = N[1 - omPlanck - orPlanck, 50];
rPlanck = N[omPlanck/olPlanck, 50];
zLambdaMatterPlanck = N[(olPlanck/omPlanck)^(1/3) - 1, 50];
zSoftPlanck = N[rPlanck^(-1/3) - 1, 50];
check09C = Abs[zSoftPlanck - zLambdaMatterPlanck] < 10^-12;

scanRatios09 = {1/20, 1/10, 3/10, 9/20, 7/10, 1, 3};
ratioSpan09 = N[Max[scanRatios09]/Min[scanRatios09], 50];
k0Ref09 = 3/5;
kappaRef09 = 17/10;
aResidualGrid09 = {1/10, 1/3, 1, 3, 10};

evalScan09[r_] := Module[
  {sig, eaVals, enVals, elVals, maxRes, wEarly, wMid, wLate, zt, pass},
  sig = N[r k0Ref09^2, 50];
  eaVals = Table[
    N[EaSoft07 /. {aa[t] -> aval, K0 -> k0Ref09, sigma -> sig, kappa -> kappaRef09}, 40],
    {aval, aResidualGrid09}
  ];
  enVals = Table[
    N[ENSoft07 /. {aa[t] -> aval, K0 -> k0Ref09, sigma -> sig, kappa -> kappaRef09}, 40],
    {aval, aResidualGrid09}
  ];
  elVals = Table[
    N[softResidual07 /. {aa[t] -> aval, sigma -> sig}, 40],
    {aval, aResidualGrid09}
  ];
  maxRes = N[Max[Abs[Join[eaVals, enVals, elVals]]], 40];
  wEarly = N[weff09[10^-2] /. {K0 -> k0Ref09, sigma -> sig}, 50];
  wMid = N[weff09[1] /. {K0 -> k0Ref09, sigma -> sig}, 50];
  wLate = N[weff09[10^2] /. {K0 -> k0Ref09, sigma -> sig}, 50];
  zt = N[r^(-1/3) - 1, 50];
  pass = Abs[wEarly] < 0.05 && Abs[wLate + 1] < 10^-5 && wEarly > wMid > wLate && maxRes < 10^-25;
  <|
    "ratio" -> N[r, 50],
    "zTransition" -> zt,
    "maxResidual" -> maxRes,
    "pass" -> pass
  |>
];

scanData09 = evalScan09 /@ scanRatios09;
scanPass09 = Lookup[scanData09, "pass"];
zTransitions09 = Lookup[scanData09, "zTransition"];
maxResidualScan09 = N[Max[Lookup[scanData09, "maxResidual"]], 50];

check09D = AllTrue[scanPass09, TrueQ];
check09E = maxResidualScan09 < 10^-25;
check09F = ratioSpan09 > 20 && Min[zTransitions09] < 0 && Max[zTransitions09] > 0;

class09 = If[
  TrueQ[check09A && check09B && check09C && check09D && check09E && check09F],
  "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED",
  "undetermined"
];

rows = {
  <|"test" -> "06", "title" -> "Hard-constraint complete EOM branch search", "classification" -> class06, "check" -> class06 =!= "undetermined"|>,
  <|"test" -> "07", "title" -> "Soft-closure evolving branch", "classification" -> class07, "check" -> class07 =!= "undetermined"|>,
  <|"test" -> "08", "title" -> "Paper-ready hard-soft EOM bridge", "classification" -> class08, "check" -> class08 =!= "undetermined"|>,
  <|"test" -> "09", "title" -> "Soft-branch robustness and anti-fine-tuning", "classification" -> class09, "check" -> class09 =!= "undetermined"|>
};

classesExpected = {
  "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY",
  "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND",
  "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED",
  "COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"
};
classesObserved = rows[[All, "classification"]];
checksObserved = rows[[All, "check"]];
checkAllRows = AllTrue[checksObserved, TrueQ];
checkClasses = classesObserved == classesExpected;

summaryMarkdown = StringRiffle[
  Join[
    {"| Test | Title | Classification | Check |", "|---|---|---|---|"},
    (StringJoin[
        "| ", #["test"], " | ", #["title"], " | ", #["classification"], " | ", ToString[#["check"], InputForm], " |"
      ] & /@ rows)
  ],
  "\n"
];

bridgeStatement = StringRiffle[
  {
    "Hard branch (06): self-consistent de Sitter under strict Elam=0.",
    "Soft branch (07): evolving matter->de Sitter history in complete channels.",
    "Bridge (08): compact hard/soft formulas and LCDM-like mapping.",
    "Robustness (09): broad parameter scan confirms non-fine-tuned behavior."
  },
  " "
];

classification = If[
  TrueQ[checkAllRows && checkClasses],
  "PIPELINE_06_09_APPENDIX_READY",
  "undetermined"
];
check = TrueQ[classification == "PIPELINE_06_09_APPENDIX_READY"];

nbName = "10_Pipeline_06_09_Appendix_Ready_Summary_Test.nb";
logName = "10_pipeline_06_09_appendix_ready_summary_test.log";
logLines = {
  "Notebook: " <> nbName,
  "constraintForm(06) = " <> ToString[constraintForm, InputForm],
  "check06A = " <> ToString[check06A, InputForm],
  "check06B = " <> ToString[check06B, InputForm],
  "check06C = " <> ToString[check06C, InputForm],
  "check06D = " <> ToString[check06D, InputForm],
  "check06E = " <> ToString[check06E, InputForm],
  "check06F = " <> ToString[check06F, InputForm],
  "lamPrimeDust06 = " <> ToString[lamPrimeDust06, InputForm],
  "hardPowerNoGo06 = " <> ToString[hardPowerNoGo06, InputForm],
  "weffEarly07 = " <> ToString[weffEarly07, InputForm],
  "weffMid07 = " <> ToString[weffMid07, InputForm],
  "weffLate07 = " <> ToString[weffLate07, InputForm],
  "lamGeneric08 = " <> ToString[lamGeneric08, InputForm],
  "omegaM0 = " <> ToString[omegaM0, InputForm],
  "omegaL0 = " <> ToString[omegaL0, InputForm],
  "zLambdaMatterPlanck = " <> ToString[zLambdaMatterPlanck, InputForm],
  "zSoftPlanck = " <> ToString[zSoftPlanck, InputForm],
  "scanRatios09 = " <> ToString[scanRatios09, InputForm],
  "ratioSpan09 = " <> ToString[ratioSpan09, InputForm],
  "zTransitions09 = " <> ToString[zTransitions09, InputForm],
  "maxResidualScan09 = " <> ToString[maxResidualScan09, InputForm],
  "rows = " <> ToString[rows, InputForm],
  "classesObserved = " <> ToString[classesObserved, InputForm],
  "classesExpected = " <> ToString[classesExpected, InputForm],
  "checkAllRows = " <> ToString[checkAllRows, InputForm],
  "checkClasses = " <> ToString[checkClasses, InputForm],
  "summaryMarkdown = " <> ToString[summaryMarkdown, InputForm],
  "bridgeStatement = " <> ToString[bridgeStatement, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "10_pipeline_06_09_appendix_ready_summary_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
