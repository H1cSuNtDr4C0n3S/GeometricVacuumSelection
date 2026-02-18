ClearAll["Global`*"];

nb = Notebook[{
  Cell["16 - FRW Minisuperspace EOM", "Title"],
  Cell["Euler-Lagrange equations for (a,N,lambda) in a flat FRW minisuperspace proxy with leafwise constraint term.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E_a = dL/da - d/dt(dL/da'),   E_N = dL/dN,   E_lambda = dL/dlambda
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deSitter branch: a=e^(H t), N=1, H^2=K0^2/3
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

constraintRule = H0^2 -> K0^2/3;
lamRule = First[Solve[Simplify[(ENDS /. constraintRule) == 0], lam0]];

EaOn = Simplify[EaDS /. constraintRule /. lamRule];
ENOn = Simplify[ENDS /. constraintRule /. lamRule];
ElamOn = Simplify[ElamDS /. constraintRule];

checkA = Simplify[ElamOn == 0];
checkB = Simplify[ENOn == 0];
checkC = Simplify[EaOn == 0];

(* numerical sanity *)
kappa0 = 2;
lambda0 = 1/5;
k00 = 3/5;
h0 = N[Sqrt[k00^2/3], 30];
lam0num = N[lam0 /. lamRule /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00}, 30];

eaNum = N[EaDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> 0.4, 30];
enNum = N[ENDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> 0.4, 30];
elamNum = N[ElamDS /. {kappa -> kappa0, Lambda -> lambda0, K0 -> k00, H0 -> h0, lam0 -> lam0num} /. t -> 0.4, 30];
maxAbsNum = N[Max[Abs[eaNum], Abs[enNum], Abs[elamNum]], 30];
checkNum = maxAbsNum < 10^-10;

check = TrueQ[checkA && checkB && checkC && checkNum];

nbName = "16_FRW_Minisuperspace_EOM.nb";
logName = "16_frw_minisuperspace_eom.log";
logLines = {
  "Notebook: " <> nbName,
  "EaDS = " <> ToString[EaDS, InputForm],
  "ENDS = " <> ToString[ENDS, InputForm],
  "ElamDS = " <> ToString[ElamDS, InputForm],
  "constraintRule = " <> ToString[constraintRule, InputForm],
  "lamRule = " <> ToString[lamRule, InputForm],
  "checkA(Elam) = " <> ToString[checkA, InputForm],
  "checkB(EN) = " <> ToString[checkB, InputForm],
  "checkC(Ea) = " <> ToString[checkC, InputForm],
  "eaNum = " <> ToString[eaNum, InputForm],
  "enNum = " <> ToString[enNum, InputForm],
  "elamNum = " <> ToString[elamNum, InputForm],
  "maxAbsNum = " <> ToString[maxAbsNum, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "16_frw_minisuperspace_eom",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxAbsNum" -> ToString[maxAbsNum, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
