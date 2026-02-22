ClearAll["Global`*"];

nb = Notebook[{
  Cell["08 - Paper-Ready Hard vs Soft EOM Bridge Test", "Title"],
  Cell["Objective: produce compact hard-vs-soft EOM formulas ready for main text/appendix and verify internal algebraic consistency against the complete-channel structure.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    hardBranch == {
      3 H^2 - K0^2 == 0,
      lambdaHard[a_] == (K0^2 - Lambda - kappa rho[a])/(2 kappa K0^2)
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    softBranch == {
      3 H^2 - K0^2 == sigma/a^3,
      H2soft[a_] == K0^2/3 + sigma/(3 a^3),
      weff[a_] == -1 + sigma/(sigma + K0^2 a^3)
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED"
  ], "Input"]
}];

rho[a_] := rho0 a^(-3 (1 + ww));

(* Normalized EN channel in FRW gauge N=1 *)
enNorm[h2_, lam_, rh_] := -((Lambda + kappa rh + 3 h2 (2 kappa lam - 1))/kappa);

lamFromEN = First[Solve[enNorm[h2, lam, rh] == 0, lam]];
lamGeneric = Simplify[lam /. lamFromEN];
lamHardExpected[rh_] := (K0^2 - Lambda - kappa rh)/(2 kappa K0^2);
lamHardFromEN[rh_] := Simplify[lamGeneric /. h2 -> K0^2/3 /. rh -> rh];

checkA = TrueQ @ Simplify[lamHardFromEN[rh] == lamHardExpected[rh], Assumptions -> {K0 != 0, kappa != 0}];

lamHardA = Simplify[lamHardExpected[rho[a]], Assumptions -> {a > 0, K0 != 0, kappa != 0}];
lamHardVac = Simplify[lamHardA /. {rho0 -> 0, Lambda -> K0^2}, Assumptions -> {K0 != 0, kappa != 0}];
checkB = TrueQ[lamHardVac == 0];

(* Soft closure branch from finite-domain motivated channel *)
h2Soft[a_] := K0^2/3 + sigma/(3 a^3);
h20 = Simplify[h2Soft[1]];
omegaM0 = Simplify[sigma/(K0^2 + sigma)];
omegaL0 = Simplify[K0^2/(K0^2 + sigma)];
h2SoftRebuild[a_] := h20 (omegaM0 a^-3 + omegaL0);

checkC = TrueQ @ Simplify[h2Soft[a] == h2SoftRebuild[a], Assumptions -> {a > 0, K0 > 0, sigma > 0}];
checkD = TrueQ @ Simplify[omegaM0 + omegaL0 == 1, Assumptions -> {K0 > 0, sigma > 0}];

weffSoft[a_] := FullSimplify[-1 - (a/(3 h2Soft[a])) D[h2Soft[x], x] /. x -> a, Assumptions -> {a > 0, K0 > 0, sigma > 0}];
weffExpected[a_] := -1 + sigma/(sigma + K0^2 a^3);
checkE = TrueQ @ Simplify[weffSoft[a] == weffExpected[a], Assumptions -> {a > 0, K0 > 0, sigma > 0}];

weffEarlyLim = FullSimplify[Limit[weffSoft[a], a -> 0, Direction -> "FromAbove"], Assumptions -> {K0 > 0, sigma > 0}];
weffLateLim = FullSimplify[Limit[weffSoft[a], a -> Infinity], Assumptions -> {K0 > 0, sigma > 0}];
checkF = TrueQ[weffEarlyLim == 0 && weffLateLim == -1];

(* Hard-soft bridge checks *)
softToHard = FullSimplify[h2Soft[a] /. sigma -> 0];
checkG = TrueQ[softToHard == K0^2/3];

(* Numeric sample for ready-to-cite values *)
numRules = {K0 -> 3/5, sigma -> 1/3, kappa -> 17/10, Lambda -> 1/5, rho0 -> 1/3, ww -> 0};
aGrid = {10^-2, 1, 10^2};
h2SoftVals = N[(h2Soft[#] /. numRules) & /@ aGrid, 50];
weffVals = N[(weffSoft[#] /. numRules) & /@ aGrid, 50];
lamHardVals = N[(lamHardA /. numRules /. a -> #) & /@ aGrid, 50];

checkH = TrueQ[weffVals[[1]] > -0.05 && Abs[weffVals[[3]] + 1] < 10^-5];
checkI = TrueQ[lamHardVals[[1]] < lamHardVals[[2]] < lamHardVals[[3]]];

paperSnippet = StringRiffle[
  {
    "Hard branch (complete channels):",
    "  3 H^2 - K0^2 = 0,   lambda_hard(a) = (K0^2 - Lambda - kappa rho(a)) / (2 kappa K0^2).",
    "Soft finite-domain closure:",
    "  3 H^2 - K0^2 = sigma a^-3,   H^2(a) = K0^2/3 + sigma/(3 a^3).",
    "  This is equivalent to H^2(a)=H0^2[Omega_m0 a^-3 + Omega_L0] with",
    "  Omega_m0 = sigma/(K0^2+sigma), Omega_L0 = K0^2/(K0^2+sigma).",
    "  Effective equation of state: w_eff(a) = -1 + sigma/(sigma + K0^2 a^3),",
    "  hence w_eff -> 0 (early) and w_eff -> -1 (late)."
  },
  "\n"
];

classification = If[
  TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI],
  "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED",
  "undetermined"
];

check = TrueQ[classification == "PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED"];

nbName = "08_Paper_Ready_Hard_Soft_EOM_Bridge_Test.nb";
logName = "08_paper_ready_hard_soft_eom_bridge_test.log";
logLines = {
  "Notebook: " <> nbName,
  "lamGeneric = " <> ToString[lamGeneric, InputForm],
  "lamHardA = " <> ToString[lamHardA, InputForm],
  "h2Soft(a) = " <> ToString[h2Soft[a], InputForm],
  "h2SoftRebuild(a) = " <> ToString[h2SoftRebuild[a], InputForm],
  "omegaM0 = " <> ToString[omegaM0, InputForm],
  "omegaL0 = " <> ToString[omegaL0, InputForm],
  "weffSoft(a) = " <> ToString[weffSoft[a], InputForm],
  "weffEarlyLim = " <> ToString[weffEarlyLim, InputForm],
  "weffLateLim = " <> ToString[weffLateLim, InputForm],
  "softToHard = " <> ToString[softToHard, InputForm],
  "aGrid = " <> ToString[aGrid, InputForm],
  "h2SoftVals = " <> ToString[h2SoftVals, InputForm],
  "weffVals = " <> ToString[weffVals, InputForm],
  "lamHardVals = " <> ToString[lamHardVals, InputForm],
  "checkA(lam hard from EN) = " <> ToString[checkA, InputForm],
  "checkB(lam hard vacuum limit) = " <> ToString[checkB, InputForm],
  "checkC(soft H2 LCDM form) = " <> ToString[checkC, InputForm],
  "checkD(omegas sum to one) = " <> ToString[checkD, InputForm],
  "checkE(weff closed form) = " <> ToString[checkE, InputForm],
  "checkF(weff limits early/late) = " <> ToString[checkF, InputForm],
  "checkG(soft->hard at sigma=0) = " <> ToString[checkG, InputForm],
  "checkH(weff numeric anchors) = " <> ToString[checkH, InputForm],
  "checkI(lamHard monotonic sample) = " <> ToString[checkI, InputForm],
  "paperSnippet = " <> ToString[paperSnippet, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "08_paper_ready_hard_soft_eom_bridge_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
