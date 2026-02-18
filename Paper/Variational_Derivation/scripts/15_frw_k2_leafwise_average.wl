ClearAll["Global`*"];

nb = Notebook[{
  Cell["15 - FRW K2 and Leafwise Average", "Title"],
  Cell["Homogeneous FRW reduction: K^2, leafwise average Q=I1/I0, and de Sitter selection rule.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Kij == (a ad/N) deltaij && Kmixed == (ad/(a N)) deltaij
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    K2 == 3 (ad/(a N))^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Q[R, t] == I1/I0 == K2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    {N == 1, a == Exp[H t], Q == 3 H^2, H^2 == K0^2/3}
  ], "Input"]
}];

a = Symbol["a"];
ad = Symbol["ad"];
n = Symbol["N"];
r = Symbol["R"];

kMat = (ad/(a n)) IdentityMatrix[3];
k2 = Simplify[Tr[kMat . kMat]];
expectedK2 = 3 ad^2/(a^2 n^2);
checkA = Simplify[k2 == expectedK2];

i0 = Simplify[(4 Pi/3) r^3 a^3];
i1 = Simplify[k2 i0];
q = Simplify[i1/i0];
checkB = Simplify[q == k2];

qdS = Simplify[q /. {a -> Exp[H t], ad -> H Exp[H t], n -> 1}];
dqdt = Simplify[D[qdS, t]];
checkC = Simplify[qdS == 3 H^2 && dqdt == 0];

selectionRule = First[Solve[3 h2 == K0^2, h2]];
checkD = Simplify[(3 h2 /. selectionRule) == K0^2];

k0val = 6/5;
hval = N[Sqrt[k0val^2/3], 30];
qNum = N[qdS /. H -> hval /. t -> 2, 30];
targetNum = N[k0val^2, 30];
numErr = N[Abs[qNum - targetNum], 30];
checkNum = numErr < 10^-12;

check = TrueQ[checkA && checkB && checkC && checkD && checkNum];

nbName = "15_FRW_K2_Leafwise_Average.nb";
logName = "15_frw_k2_leafwise_average.log";
logLines = {
  "Notebook: " <> nbName,
  "K2 = " <> ToString[k2, InputForm],
  "expectedK2 = " <> ToString[expectedK2, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "I0 = " <> ToString[i0, InputForm],
  "I1 = " <> ToString[i1, InputForm],
  "Q = " <> ToString[q, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "QdS = " <> ToString[qdS, InputForm],
  "dQdt_dS = " <> ToString[dqdt, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "selectionRule(H^2) = " <> ToString[selectionRule, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "qNum = " <> ToString[qNum, InputForm],
  "targetNum = " <> ToString[targetNum, InputForm],
  "numErr = " <> ToString[numErr, InputForm],
  "checkNum = " <> ToString[checkNum, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "15_frw_k2_leafwise_average",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check,
  "numErr" -> ToString[numErr, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
