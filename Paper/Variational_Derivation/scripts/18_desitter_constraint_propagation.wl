ClearAll["Global`*"];

nb = Notebook[{
  Cell["18 - de Sitter Constraint Propagation", "Title"],
  Cell["Minisuperspace consistency check: once (Ea, EN, Elam)=0 on the de Sitter branch, their time derivatives vanish as well.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_a = dL/da - d/dt(dL/da'),   E_N = dL/dN,   E_lambda = dL/dlambda
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deSitter: a = Exp[H0 t], N = 1, H0^2 = K0^2/3, lambda = (K0^2 - Lambda)/(2 K0^2 kappa)
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

subDS = {
  aa[t] -> Exp[H0 t],
  Derivative[1][aa][t] -> H0 Exp[H0 t],
  Derivative[2][aa][t] -> H0^2 Exp[H0 t],
  NN[t] -> 1,
  Derivative[1][NN][t] -> 0,
  ll[t] -> lam0,
  Derivative[1][ll][t] -> 0
};

EaDS = Simplify[Ea /. subDS];
ENDS = Simplify[EN /. subDS];
ElamDS = Simplify[Elam /. subDS];

branchRule = H0^2 -> K0^2/3;
lamRule = First[Solve[Simplify[(ENDS /. branchRule) == 0], lam0]];

EaOn = Simplify[EaDS /. branchRule /. lamRule];
ENOn = Simplify[ENDS /. branchRule /. lamRule];
ElamOn = Simplify[ElamDS /. branchRule];

dEaOn = Simplify[D[EaOn, t]];
dENOn = Simplify[D[ENOn, t]];
dElamOn = Simplify[D[ElamOn, t]];

checkA = Simplify[EaOn == 0];
checkB = Simplify[ENOn == 0];
checkC = Simplify[ElamOn == 0];
checkD = Simplify[dEaOn == 0];
checkE = Simplify[dENOn == 0];
checkF = Simplify[dElamOn == 0];

(* numerical sanity *)
kappa0 = 13/10;
lambda0 = 7/100;
k00 = 4/5;
h0 = N[Sqrt[k00^2/3], 30];
lam0num = N[lam0 /. lamRule /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00}, 30];

eaVals = Table[
  N[EaDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> tt, 30],
  {tt, {0, 1/2, 1}}
];
enVals = Table[
  N[ENDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> tt, 30],
  {tt, {0, 1/2, 1}}
];
elVals = Table[
  N[ElamDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> tt, 30],
  {tt, {0, 1/2, 1}}
];

maxAbsNum = N[Max[Join[Abs[eaVals], Abs[enVals], Abs[elVals]]], 30];
checkNum = maxAbsNum < 10^-10;

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkNum];

nbName = "18_deSitter_Constraint_Propagation.nb";
logName = "18_desitter_constraint_propagation.log";
logLines = {
  "Notebook: " <> nbName,
  "EaDS = " <> ToString[EaDS, InputForm],
  "ENDS = " <> ToString[ENDS, InputForm],
  "ElamDS = " <> ToString[ElamDS, InputForm],
  "branchRule = " <> ToString[branchRule, InputForm],
  "lamRule = " <> ToString[lamRule, InputForm],
  "EaOn = " <> ToString[EaOn, InputForm],
  "ENOn = " <> ToString[ENOn, InputForm],
  "ElamOn = " <> ToString[ElamOn, InputForm],
  "dEaOn = " <> ToString[dEaOn, InputForm],
  "dENOn = " <> ToString[dENOn, InputForm],
  "dElamOn = " <> ToString[dElamOn, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "checkE = " <> ToString[checkE, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "eaVals = " <> ToString[eaVals, InputForm],
  "enVals = " <> ToString[enVals, InputForm],
  "elVals = " <> ToString[elVals, InputForm],
  "maxAbsNum = " <> ToString[maxAbsNum, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "18_desitter_constraint_propagation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxAbsNum" -> ToString[maxAbsNum, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
