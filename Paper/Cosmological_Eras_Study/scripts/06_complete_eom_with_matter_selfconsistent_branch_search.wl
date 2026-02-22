ClearAll["Global`*"];

nb = Notebook[{
  Cell["06 - Complete EOM With Matter: Self-Consistent Branch Search", "Title"],
  Cell["Objective: solve the complete minisuperspace EOM channels (Ea, EN, Elam) with barotropic matter, identify an exact self-consistent branch, and test whether a standard power-law expanding branch can satisfy the hard constraint channel.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L == -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3) - rho0 a^(-3 w) N
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {Ea == 0, EN == 0, Elam == 0}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Elam == 0 \[Implies] H^2 == K0^2/3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY"
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

ElamN1 = Simplify[Elam /. {NN[t] -> 1}];
constraintForm = Simplify[ElamN1/aa[t]^3];
checkA = constraintForm === (3 (Derivative[1][aa][t]^2/aa[t]^2) - K0^2);

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

hardRule = H0^2 -> K0^2/3;
lamRule = First[Solve[Simplify[(ENDS /. hardRule) == 0], ll[t]]];
lamExpr = Simplify[ll[t] /. lamRule];
dlamRule = Derivative[1][ll][t] -> D[lamExpr, t];

EaCons = Simplify[EaDS /. hardRule /. lamRule /. dlamRule];
ENCons = Simplify[ENDS /. hardRule /. lamRule];
ElamCons = Simplify[ElamDS /. hardRule];

checkB = Simplify[ElamCons == 0];
checkC = Simplify[ENCons == 0];
checkD = TrueQ @ Simplify[(EaCons /. hardRule) == 0, Assumptions -> K0 != 0];

rhoSol[t_] := rho0 Exp[-3 H0 (1 + ww) t];
rhoConservation = Simplify[D[rhoSol[t], t] + 3 H0 (1 + ww) rhoSol[t]];
checkE = Simplify[rhoConservation == 0];

lamPrimeExpr = Simplify[D[lamExpr, t]];
lamPrimeDust = Simplify[
  lamPrimeExpr /. {
    ww -> 0, rho0 -> 1, kappa -> 17/10, Lambda -> 1/5,
    K0 -> 3/5, H0 -> Sqrt[(3/5)^2/3], t -> 1
  }
];
checkF = N[Abs[lamPrimeDust], 40] > 10^-8;

lamPrimeDE = Simplify[lamPrimeExpr /. ww -> -1];
checkG = Simplify[lamPrimeDE == 0];

aPower = t^p;
ElamPower = Simplify[
  Elam /. {
    aa[t] -> aPower,
    Derivative[1][aa][t] -> D[aPower, t],
    NN[t] -> 1
  },
  Assumptions -> t > 0
];

powerNoGoReduce = Reduce[
  {
    (ElamPower /. {K0 -> 1, t -> 1}) == 0,
    (ElamPower /. {K0 -> 1, t -> 2}) == 0
  },
  p, Reals
];
checkH = TrueQ[powerNoGoReduce === False];

sampleResiduals = Table[
  N[
    {
      EaCons,
      ENCons,
      ElamCons
    } /. {
      ww -> 0,
      rho0 -> 1,
      kappa -> 17/10,
      Lambda -> 1/5,
      K0 -> 3/5,
      H0 -> Sqrt[(3/5)^2/3],
      t -> tt
    },
    40
  ],
  {tt, {0, 1/2, 1, 2}}
];
maxResidual = N[Max[Abs[Flatten[sampleResiduals]]], 40];
checkI = maxResidual < 10^-25;

classification = If[
  TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI],
  "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY",
  "undetermined"
];

check = TrueQ[classification == "COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY"];

nbName = "06_Complete_EOM_With_Matter_Selfconsistent_Branch_Search.nb";
logName = "06_complete_eom_with_matter_selfconsistent_branch_search.log";
logLines = {
  "Notebook: " <> nbName,
  "Ea = " <> ToString[Ea, InputForm],
  "EN = " <> ToString[EN, InputForm],
  "Elam = " <> ToString[Elam, InputForm],
  "ElamN1 = " <> ToString[ElamN1, InputForm],
  "constraintForm = " <> ToString[constraintForm, InputForm],
  "EaDS = " <> ToString[EaDS, InputForm],
  "ENDS = " <> ToString[ENDS, InputForm],
  "ElamDS = " <> ToString[ElamDS, InputForm],
  "hardRule = " <> ToString[hardRule, InputForm],
  "lamRule = " <> ToString[lamRule, InputForm],
  "lamExpr = " <> ToString[lamExpr, InputForm],
  "lamPrimeExpr = " <> ToString[lamPrimeExpr, InputForm],
  "EaCons = " <> ToString[EaCons, InputForm],
  "ENCons = " <> ToString[ENCons, InputForm],
  "ElamCons = " <> ToString[ElamCons, InputForm],
  "rhoConservation = " <> ToString[rhoConservation, InputForm],
  "lamPrimeDust = " <> ToString[lamPrimeDust, InputForm],
  "lamPrimeDE = " <> ToString[lamPrimeDE, InputForm],
  "ElamPower = " <> ToString[ElamPower, InputForm],
  "powerNoGoReduce = " <> ToString[powerNoGoReduce, InputForm],
  "sampleResiduals = " <> ToString[sampleResiduals, InputForm],
  "maxResidual = " <> ToString[maxResidual, InputForm],
  "checkA(constraint channel form) = " <> ToString[checkA, InputForm],
  "checkB(Elam on hard branch) = " <> ToString[checkB, InputForm],
  "checkC(EN solved by lambda) = " <> ToString[checkC, InputForm],
  "checkD(Ea consistency on branch) = " <> ToString[checkD, InputForm],
  "checkE(matter continuity eq) = " <> ToString[checkE, InputForm],
  "checkF(lambda time dependence for dust) = " <> ToString[checkF, InputForm],
  "checkG(lambda static for w=-1) = " <> ToString[checkG, InputForm],
  "checkH(power-law no-go under hard constraint) = " <> ToString[checkH, InputForm],
  "checkI(numerical residuals) = " <> ToString[checkI, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "06_complete_eom_with_matter_selfconsistent_branch_search",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
