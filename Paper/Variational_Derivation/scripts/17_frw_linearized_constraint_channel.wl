ClearAll["Global`*"];

nb = Notebook[{
  Cell["17 - FRW Linearized Constraint Channel", "Title"],
  Cell["First-order perturbation of C = 3 a a'^2/N^2 - K0^2 a^3 around de Sitter.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {a[t] == Exp[H t] (1 + eps u[t]), N[t] == 1 + eps n[t]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    C1 == 6 H Exp[3 H t] (Derivative[1][u][t] - H n[t])
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    n[t] == 0 \[Implies] Derivative[1][u][t] == 0
  ], "Input"]
}];

t = Symbol["t"];
a0[t_] := Exp[H0 t];

aPert = a0[t] (1 + eps u[t]);
nPert = 1 + eps n[t];

constraintExpr = Simplify[3 aPert D[aPert, t]^2/nPert^2 - K0^2 aPert^3];
c0 = Simplify[SeriesCoefficient[constraintExpr, {eps, 0, 0}]];
c1 = Simplify[SeriesCoefficient[constraintExpr, {eps, 0, 1}]];

branchRule = H0^2 -> K0^2/3;
c0On = Simplify[c0 /. branchRule];
c1On = Simplify[c1 /. branchRule];

expectedC1 = Simplify[6 H0 Exp[3 H0 t] (Derivative[1][u][t] - H0 n[t])];
checkA = Simplify[c0On == 0];
checkBExpr = Simplify[(c1On - expectedC1) /. K0^2 -> 3 H0^2];
checkB = Simplify[checkBExpr == 0];

c1Gauge = Simplify[c1On /. n[t] -> 0];
expectedGauge = Simplify[6 H0 Exp[3 H0 t] Derivative[1][u][t]];
checkC = Simplify[c1Gauge == expectedGauge];

c1Const = Simplify[c1Gauge /. {u[t] -> u0, Derivative[1][u][t] -> 0}];
checkD = Simplify[c1Const == 0];

(* numerical sanity on constrained linear mode: n = u'/H0 *)
hval = 7/10;
k0val = Sqrt[3] hval;
samplePts = {0, 1/2, 1, 3/2};
c1Vals = Table[
  N[
    c1On /. {
      H0 -> hval, K0 -> k0val,
      u[t] -> Sin[t], Derivative[1][u][t] -> Cos[t],
      n[t] -> Cos[t]/hval,
      t -> tt
    },
    30
  ],
  {tt, samplePts}
];
maxAbs = N[Max[Abs[c1Vals]], 30];
checkNum = maxAbs < 10^-12;

check = TrueQ[checkA && checkB && checkC && checkD && checkNum];

nbName = "17_FRW_Linearized_Constraint_Channel.nb";
logName = "17_frw_linearized_constraint_channel.log";
logLines = {
  "Notebook: " <> nbName,
  "C0 = " <> ToString[c0, InputForm],
  "C1 = " <> ToString[c1, InputForm],
  "C0_onBranch = " <> ToString[c0On, InputForm],
  "C1_onBranch = " <> ToString[c1On, InputForm],
  "expectedC1 = " <> ToString[expectedC1, InputForm],
  "checkBExpr = " <> ToString[checkBExpr, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "C1_gauge_n0 = " <> ToString[c1Gauge, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "samplePts = " <> ToString[samplePts, InputForm],
  "c1Vals(constrained) = " <> ToString[c1Vals, InputForm],
  "maxAbs = " <> ToString[maxAbs, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "17_frw_linearized_constraint_channel",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxAbs" -> ToString[maxAbs, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
