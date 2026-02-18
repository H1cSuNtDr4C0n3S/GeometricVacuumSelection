ClearAll["Global`*"];

nb = Notebook[{
  Cell["00 - Barotropic Hard-Constraint Consistency", "Title"],
  Cell["FRW minisuperspace with barotropic matter: test consistency of a hard geometric constraint H^2=K0^2/3.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    rho[a_] == rho0 a^(-3 (1 + w))
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L == -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3) - rho0 a^(-3 w) N
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    hardConstraint == Elam == 0 \[Implies] H^2 == K0^2/3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    with EN==0, lambda[t] is fixed; Ea checks whether a genuine dust/radiation era survives under hard constraint
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
dlamRule = Derivative[1][ll][t] -> D[(ll[t] /. lamRule), t];

lamExpr = Simplify[ll[t] /. lamRule];
lamPrimeExpr = Simplify[D[lamExpr, t]];
EaCons = Together @ Simplify[EaDS /. hardRule /. lamRule /. dlamRule];

checkA = Simplify[(ElamDS /. hardRule) == 0];
checkB = Simplify[(EaCons /. rho0 -> 0) == 0];
checkC = Simplify[(EaCons /. ww -> -1) == 0];

matterValsDust = Table[
  N[
    EaCons /. {
      ww -> 0,
      rho0 -> 1,
      kappa -> 17/10,
      Lambda -> 1/5,
      K0 -> 3/5,
      H0 -> Sqrt[(3/5)^2/3],
      t -> tt
    },
    30
  ],
  {tt, {0, 1/2, 1, 3/2}}
];
maxDust = N[Max[Abs[matterValsDust]], 30];
checkD = maxDust < 10^-8;

lamPrimeDustVals = Table[
  N[
    lamPrimeExpr /. {
      ww -> 0,
      rho0 -> 1,
      kappa -> 17/10,
      Lambda -> 1/5,
      K0 -> 3/5,
      H0 -> Sqrt[(3/5)^2/3],
      t -> tt
    },
    30
  ],
  {tt, {0, 1/2, 1}}
];
maxLamPrimeDust = N[Max[Abs[lamPrimeDustVals]], 30];
checkE = maxLamPrimeDust > 10^-8;
lamPrimeDE = Simplify[lamPrimeExpr /. ww -> -1];
checkF = Simplify[lamPrimeDE == 0];

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF];

nbName = "00_Barotropic_Hard_Constraint_Consistency.nb";
logName = "00_barotropic_hard_constraint_consistency.log";
logLines = {
  "Notebook: " <> nbName,
  "EaDS = " <> ToString[EaDS, InputForm],
  "ENDS = " <> ToString[ENDS, InputForm],
  "ElamDS = " <> ToString[ElamDS, InputForm],
  "hardRule = " <> ToString[hardRule, InputForm],
  "lamRule = " <> ToString[lamRule, InputForm],
  "lamExpr = " <> ToString[lamExpr, InputForm],
  "lamPrimeExpr = " <> ToString[lamPrimeExpr, InputForm],
  "EaCons = " <> ToString[EaCons, InputForm],
  "checkA(Elam) = " <> ToString[checkA, InputForm],
  "checkB(vacuum) = " <> ToString[checkB, InputForm],
  "checkC(w=-1) = " <> ToString[checkC, InputForm],
  "matterValsDust = " <> ToString[matterValsDust, InputForm],
  "maxDust = " <> ToString[maxDust, InputForm],
  "checkD(Ea redundancy on hard branch) = " <> ToString[checkD, InputForm],
  "lamPrimeDustVals = " <> ToString[lamPrimeDustVals, InputForm],
  "maxLamPrimeDust = " <> ToString[maxLamPrimeDust, InputForm],
  "checkE(lambda time dependence for dust) = " <> ToString[checkE, InputForm],
  "lamPrimeDE(w=-1) = " <> ToString[lamPrimeDE, InputForm],
  "checkF(lambda static for w=-1) = " <> ToString[checkF, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "00_barotropic_hard_constraint_consistency",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
