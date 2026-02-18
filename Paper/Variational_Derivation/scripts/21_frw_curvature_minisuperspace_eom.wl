ClearAll["Global`*"];

nb = Notebook[{
  Cell["21 - FRW Curvature Minisuperspace EOM", "Title"],
  Cell["FRW minisuperspace extension with spatial-curvature term +3 k a N/kappa and consistency checks vs k->0 limit.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lk = -3 a a'^2/(kappa N) + 3 k a N/kappa - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Ea == dLk/da - d/dt(dLk/da'), EN == dLk/dN, Elam == dLk/dlambda
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    k -> 0  \[Implies]  recovery of the flat-FRW EOM used in notebook 16
  ], "Input"]
}];

t = Symbol["t"];
aT = aa[t];
nT = NN[t];
lamT = ll[t];

Lk = -3 aT D[aT, t]^2/(kappa nT) + 3 kk aT nT/kappa - (Lambda aT^3 nT)/kappa +
  lamT (3 aT D[aT, t]^2/nT^2 - K0^2 aT^3);

L0 = -3 aT D[aT, t]^2/(kappa nT) - (Lambda aT^3 nT)/kappa +
  lamT (3 aT D[aT, t]^2/nT^2 - K0^2 aT^3);

EaK = Simplify[D[Lk, aT] - D[D[Lk, D[aT, t]], t]];
ENK = Simplify[D[Lk, nT]];
ElamK = Simplify[D[Lk, lamT]];

Ea0 = Simplify[D[L0, aT] - D[D[L0, D[aT, t]], t]];
EN0 = Simplify[D[L0, nT]];
Elam0 = Simplify[D[L0, lamT]];

checkA = Simplify[(EaK /. kk -> 0) == Ea0];
checkB = Simplify[(ENK /. kk -> 0) == EN0];
checkC = Simplify[(ElamK /. kk -> 0) == Elam0];
checkD = Simplify[D[EaK, kk] == 3 nT/kappa];
checkE = Simplify[D[ENK, kk] == 3 aT/kappa];
checkF = FreeQ[ElamK, kk];

(* numerical sanity: finite values on a nontrivial profile *)
aNum[t_] := 2 + (1/10) t + (1/30) t^2 + (1/300) t^3;
nNum[t_] := 1 + (1/20) t + (1/200) t^2;
lNum[t_] := 1/4 + (1/25) t + (1/200) t^2;

numSub = {
  aa[t] -> aNum[t],
  Derivative[1][aa][t] -> D[aNum[t], t],
  Derivative[2][aa][t] -> D[aNum[t], {t, 2}],
  NN[t] -> nNum[t],
  Derivative[1][NN][t] -> D[nNum[t], t],
  ll[t] -> lNum[t],
  Derivative[1][ll][t] -> D[lNum[t], t],
  kappa -> 11/10,
  Lambda -> 1/8,
  K0 -> 2/3,
  kk -> 1/5
};

samplePts = {0, 1/2, 1, 3/2};
eaVals = Table[N[EaK /. numSub /. t -> tt, 30], {tt, samplePts}];
enVals = Table[N[ENK /. numSub /. t -> tt, 30], {tt, samplePts}];
elVals = Table[N[ElamK /. numSub /. t -> tt, 30], {tt, samplePts}];

finiteCheck = And @@ Flatten[Map[(NumericQ[#] && # =!= Indeterminate && # =!= ComplexInfinity) &, {eaVals, enVals, elVals}, {2}]];
checkNum = TrueQ[finiteCheck];

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkNum];

nbName = "21_FRW_Curvature_Minisuperspace_EOM.nb";
logName = "21_frw_curvature_minisuperspace_eom.log";
logLines = {
  "Notebook: " <> nbName,
  "EaK = " <> ToString[EaK, InputForm],
  "ENK = " <> ToString[ENK, InputForm],
  "ElamK = " <> ToString[ElamK, InputForm],
  "checkA(k->0 Ea) = " <> ToString[checkA, InputForm],
  "checkB(k->0 EN) = " <> ToString[checkB, InputForm],
  "checkC(k->0 Elam) = " <> ToString[checkC, InputForm],
  "checkD(dEa/dk) = " <> ToString[checkD, InputForm],
  "checkE(dEN/dk) = " <> ToString[checkE, InputForm],
  "checkF(Elam independent of k) = " <> ToString[checkF, InputForm],
  "samplePts = " <> ToString[samplePts, InputForm],
  "eaVals = " <> ToString[eaVals, InputForm],
  "enVals = " <> ToString[enVals, InputForm],
  "elVals = " <> ToString[elVals, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "21_frw_curvature_minisuperspace_eom",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
