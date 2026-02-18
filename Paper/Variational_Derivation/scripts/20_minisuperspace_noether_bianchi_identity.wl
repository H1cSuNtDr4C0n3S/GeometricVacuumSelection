ClearAll["Global`*"];

nb = Notebook[{
  Cell["20 - Minisuperspace Noether-Bianchi Identity", "Title"],
  Cell["Off-shell gauge identity for time-reparametrization in FRW minisuperspace: structural relation among Ea, EN, Elam.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    L = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Ea == dL/da - d/dt(dL/da'), EN == dL/dN, Elam == dL/dlambda
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    noetherId == Ea a' + EN N' + Elam lambda' - d/dt(N EN) == 0
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

noetherRaw = Together[Simplify[Ea D[aT, t] + EN D[nT, t] + Elam D[lamT, t] - D[nT EN, t]]];
auxTerm = Together[Simplify[D[lamT Elam, t]]];
noetherId = Together[Simplify[noetherRaw - auxTerm]];
checkA = Simplify[noetherId == 0];

(* equivalent solved form for EN propagation *)
propEq = Simplify[D[nT EN + lamT Elam, t] - (Ea D[aT, t] + EN D[nT, t] + Elam D[lamT, t])];
checkB = Simplify[propEq == 0];

(* de Sitter branch implication: if Ea=Elam=EN=0 then d/dt(N EN)=0 trivially *)
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

idDS = Simplify[noetherId /. subDS /. branchRule /. lamRule];
checkC = Simplify[idDS == 0];

(* numeric sanity with nontrivial profiles *)
aNum[t_] := 2 + (1/10) t + (1/50) t^2 + (1/200) t^3;
nNum[t_] := 1 + (1/20) t + (1/100) t^2;
lNum[t_] := 1/5 + (1/30) t + (1/150) t^2;

numSub = {
  aa[t] -> aNum[t],
  Derivative[1][aa][t] -> D[aNum[t], t],
  Derivative[2][aa][t] -> D[aNum[t], {t, 2}],
  NN[t] -> nNum[t],
  Derivative[1][NN][t] -> D[nNum[t], t],
  ll[t] -> lNum[t],
  Derivative[1][ll][t] -> D[lNum[t], t],
  kappa -> 7/5,
  Lambda -> 3/20,
  K0 -> 2/3
};

samplePts = {0, 1/2, 1, 3/2};
idVals = Table[N[noetherId /. numSub /. t -> tt, 30], {tt, samplePts}];
maxAbs = N[Max[Abs[idVals]], 30];
checkNum = maxAbs < 10^-12;

check = TrueQ[checkA && checkB && checkC && checkNum];

nbName = "20_Minisuperspace_Noether_Bianchi_Identity.nb";
logName = "20_minisuperspace_noether_bianchi_identity.log";
logLines = {
  "Notebook: " <> nbName,
  "Ea = " <> ToString[Ea, InputForm],
  "EN = " <> ToString[EN, InputForm],
  "Elam = " <> ToString[Elam, InputForm],
  "noetherRaw = " <> ToString[noetherRaw, InputForm],
  "auxTerm = " <> ToString[auxTerm, InputForm],
  "noetherId = " <> ToString[noetherId, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "propEq = " <> ToString[propEq, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "idDS = " <> ToString[idDS, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "samplePts = " <> ToString[samplePts, InputForm],
  "idVals = " <> ToString[idVals, InputForm],
  "maxAbs = " <> ToString[maxAbs, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "20_minisuperspace_noether_bianchi_identity",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "maxAbs" -> ToString[maxAbs, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
