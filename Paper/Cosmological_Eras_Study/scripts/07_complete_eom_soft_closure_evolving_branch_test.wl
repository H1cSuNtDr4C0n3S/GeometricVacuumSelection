ClearAll["Global`*"];

nb = Notebook[{
  Cell["07 - Complete EOM Soft-Closure Evolving Branch Test", "Title"],
  Cell["Objective: test a finite-domain inspired soft closure for the constraint channel and verify existence of an evolving matter->de Sitter branch that is self-consistent in the complete EOM channels.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    softConstraintNormalized == (3 H^2 - K0^2) == sigma/a^3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    H2soft[a_] := K0^2/3 + sigma/(3 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {Ea, EN} == 0 && (Elam/a^3) == sigma/a^3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    weff[a_] := -1 - (1/3) D[Log[H2soft[a]], Log[a]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND"
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

(* Soft closure: normalized constraint residual decays as a^-3 *)
h2Soft[a_] := K0^2/3 + sigma/(3 a^3);
hSoft[a_] := Sqrt[h2Soft[a]];
hdotSoft[a_] := -sigma/(2 a^3); (* from d(H^2)/dt = -sigma H/a^3 *)
weff[a_] := FullSimplify[-1 - (a/3) D[Log[h2Soft[x]], x] /. x -> a, Assumptions -> a > 0];

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

checkA = TrueQ[EaSoft == 0];
checkB = TrueQ[ENSoft == 0];
checkC = TrueQ[softConstraintResidual == 0];

(* Asymptotic behavior of effective equation of state *)
weffEarly = N[weff[10^-2] /. {K0 -> 3/5, sigma -> 1/3}, 50];
weffMid = N[weff[1] /. {K0 -> 3/5, sigma -> 1/3}, 50];
weffLate = N[weff[100] /. {K0 -> 3/5, sigma -> 1/3}, 50];
checkD = Abs[weffEarly] < 0.05;
checkE = Abs[weffLate + 1] < 0.01;
checkF = weffEarly > weffMid > weffLate;

(* Direct numerical residual spot checks *)
eaNumVals = Table[
  N[
    EaSoft /. {
      aa[t] -> aval,
      kappa -> 17/10,
      K0 -> 3/5,
      sigma -> 1/3
    },
    40
  ],
  {aval, {1/10, 1/3, 1, 3, 10}}
];
enNumVals = Table[
  N[
    ENSoft /. {
      aa[t] -> aval,
      kappa -> 17/10,
      K0 -> 3/5,
      sigma -> 1/3
    },
    40
  ],
  {aval, {1/10, 1/3, 1, 3, 10}}
];
elamNumVals = Table[
  N[
    softConstraintResidual /. {
      aa[t] -> aval,
      sigma -> 1/3
    },
    40
  ],
  {aval, {1/10, 1/3, 1, 3, 10}}
];
maxResidual = N[
  Max[Abs[Join[eaNumVals, enNumVals, elamNumVals]]],
  40
];
checkG = maxResidual < 10^-25;

(* Hard branch comparison: under hard channel, power-law fails *)
ElamPowerHard = Simplify[
  Elam /. {
    aa[t] -> t^p,
    Derivative[1][aa][t] -> D[t^p, t],
    NN[t] -> 1
  },
  Assumptions -> t > 0
];
hardPowerNoGo = Reduce[
  {
    (ElamPowerHard /. {K0 -> 1, t -> 1}) == 0,
    (ElamPowerHard /. {K0 -> 1, t -> 2}) == 0
  },
  p, Reals
];
checkH = TrueQ[hardPowerNoGo === False];

classification = If[
  TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH],
  "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND",
  "undetermined"
];

check = TrueQ[classification == "COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND"];

nbName = "07_Complete_EOM_Soft_Closure_Evolving_Branch_Test.nb";
logName = "07_complete_eom_soft_closure_evolving_branch_test.log";
logLines = {
  "Notebook: " <> nbName,
  "Ea = " <> ToString[Ea, InputForm],
  "EN = " <> ToString[EN, InputForm],
  "Elam = " <> ToString[Elam, InputForm],
  "h2Soft(a) = " <> ToString[h2Soft[aT], InputForm],
  "EaSoft = " <> ToString[EaSoft, InputForm],
  "ENSoft = " <> ToString[ENSoft, InputForm],
  "ElamSoft = " <> ToString[ElamSoft, InputForm],
  "softConstraintResidual = " <> ToString[softConstraintResidual, InputForm],
  "weff(a) = " <> ToString[weff[aT], InputForm],
  "weffEarly(a=1e-2) = " <> ToString[weffEarly, InputForm],
  "weffMid(a=1) = " <> ToString[weffMid, InputForm],
  "weffLate(a=1e2) = " <> ToString[weffLate, InputForm],
  "eaNumVals = " <> ToString[eaNumVals, InputForm],
  "enNumVals = " <> ToString[enNumVals, InputForm],
  "elamNumVals = " <> ToString[elamNumVals, InputForm],
  "maxResidual = " <> ToString[maxResidual, InputForm],
  "ElamPowerHard = " <> ToString[ElamPowerHard, InputForm],
  "hardPowerNoGo = " <> ToString[hardPowerNoGo, InputForm],
  "checkA(Ea channel) = " <> ToString[checkA, InputForm],
  "checkB(EN channel) = " <> ToString[checkB, InputForm],
  "checkC(soft constraint channel) = " <> ToString[checkC, InputForm],
  "checkD(early matter-like weff) = " <> ToString[checkD, InputForm],
  "checkE(late deSitter-like weff) = " <> ToString[checkE, InputForm],
  "checkF(monotonic transition weff) = " <> ToString[checkF, InputForm],
  "checkG(numerical residuals) = " <> ToString[checkG, InputForm],
  "checkH(hard-branch power-law no-go) = " <> ToString[checkH, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "07_complete_eom_soft_closure_evolving_branch_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
