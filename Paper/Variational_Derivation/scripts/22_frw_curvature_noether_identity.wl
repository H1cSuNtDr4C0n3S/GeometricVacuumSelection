ClearAll["Global`*"];

nb = Notebook[{
  Cell["22 - FRW Curvature Noether Identity", "Title"],
  Cell["Noether/Bianchi-type identity in FRW minisuperspace with spatial curvature k.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lk = -3 a a'^2/(kappa N) + 3 k a N/kappa - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    noetherRaw = Ea a' + EN N' + Elam lambda' - d/dt(N EN),  noetherId = noetherRaw - d/dt(lambda Elam)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    noetherId == 0  (off-shell)
  ], "Input"]
}];

t = Symbol["t"];
aT = aa[t];
nT = NN[t];
lamT = ll[t];

Lk = -3 aT D[aT, t]^2/(kappa nT) + 3 kk aT nT/kappa - (Lambda aT^3 nT)/kappa +
  lamT (3 aT D[aT, t]^2/nT^2 - K0^2 aT^3);

Ea = Simplify[D[Lk, aT] - D[D[Lk, D[aT, t]], t]];
EN = Simplify[D[Lk, nT]];
Elam = Simplify[D[Lk, lamT]];

noetherRaw = Together[Simplify[Ea D[aT, t] + EN D[nT, t] + Elam D[lamT, t] - D[nT EN, t]]];
auxTerm = Together[Simplify[D[lamT Elam, t]]];
noetherId = Together[Simplify[noetherRaw - auxTerm]];

checkA = Simplify[noetherId == 0];
checkB = Simplify[(noetherId /. kk -> 0) == 0];
propEq = Simplify[D[nT EN + lamT Elam, t] - (Ea D[aT, t] + EN D[nT, t] + Elam D[lamT, t])];
checkC = Simplify[propEq == 0];

(* numeric sanity *)
aNum[t_] := 2 + (1/8) t + (1/40) t^2;
nNum[t_] := 1 + (1/25) t + (1/250) t^2;
lNum[t_] := 3/10 + (1/30) t + (1/400) t^2;

numSub = {
  aa[t] -> aNum[t],
  Derivative[1][aa][t] -> D[aNum[t], t],
  Derivative[2][aa][t] -> D[aNum[t], {t, 2}],
  NN[t] -> nNum[t],
  Derivative[1][NN][t] -> D[nNum[t], t],
  ll[t] -> lNum[t],
  Derivative[1][ll][t] -> D[lNum[t], t],
  kappa -> 7/5,
  Lambda -> 1/10,
  K0 -> 3/5,
  kk -> 1/4
};

samplePts = {0, 1/2, 1, 3/2};
idVals = Table[N[noetherId /. numSub /. t -> tt, 30], {tt, samplePts}];
maxAbs = N[Max[Abs[idVals]], 30];
checkNum = maxAbs < 10^-12;

check = TrueQ[checkA && checkB && checkC && checkNum];

nbName = "22_FRW_Curvature_Noether_Identity.nb";
logName = "22_frw_curvature_noether_identity.log";
logLines = {
  "Notebook: " <> nbName,
  "Ea = " <> ToString[Ea, InputForm],
  "EN = " <> ToString[EN, InputForm],
  "Elam = " <> ToString[Elam, InputForm],
  "noetherRaw = " <> ToString[noetherRaw, InputForm],
  "auxTerm = " <> ToString[auxTerm, InputForm],
  "noetherId = " <> ToString[noetherId, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB(k->0) = " <> ToString[checkB, InputForm],
  "propEq = " <> ToString[propEq, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "samplePts = " <> ToString[samplePts, InputForm],
  "idVals = " <> ToString[idVals, InputForm],
  "maxAbs = " <> ToString[maxAbs, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "22_frw_curvature_noether_identity",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
