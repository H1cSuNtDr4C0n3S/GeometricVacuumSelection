ClearAll["Global`*"];

nb = Notebook[{
  Cell["04 - Covariant Ratio Second Variation", "Title"],
  Cell["Second variation of Q = I1/I0 in covariant finite-domain form, including IR scaling of first and second variation channels.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    qeps == (i1 + eps di1 + eps^2 d2i1/2)/(i0 + eps di0 + eps^2 d2i0/2)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q1 == D[qeps, eps] /. eps -> 0 && q2 == D[qeps, {eps, 2}] /. eps -> 0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q1Branch == (di1 - q di0)/i0 && q2Branch == (d2i1 - q d2i0)/i0 - 2 di0 (di1 - q di0)/i0^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    i0 -> L^3  \[Implies] q1, q2 are IR-suppressed
  ], "Input"]
}];

qeps = (i1 + eps di1 + eps^2 d2i1/2)/(i0 + eps di0 + eps^2 d2i0/2);
q1 = FullSimplify[D[qeps, eps] /. eps -> 0];
q2 = FullSimplify[D[qeps, {eps, 2}] /. eps -> 0];

q1Expected = FullSimplify[(di1 i0 - di0 i1)/i0^2];
q2Expected = FullSimplify[(d2i1 i0 - d2i0 i1)/i0^2 - 2 di0 (di1 i0 - di0 i1)/i0^3];

checkA = FullSimplify[q1 == q1Expected];
checkB = FullSimplify[q2 == q2Expected];

q1Branch = FullSimplify[q1 /. i1 -> q i0];
q2Branch = FullSimplify[q2 /. i1 -> q i0];

q1BranchExpected = FullSimplify[(di1 - q di0)/i0];
q2BranchExpected = FullSimplify[(d2i1 - q d2i0)/i0 - 2 di0 (di1 - q di0)/i0^2];

checkC = FullSimplify[q1Branch == q1BranchExpected];
checkD = FullSimplify[q2Branch == q2BranchExpected];

numRule = {
  q -> 3/5,
  di0 -> 2/7,
  di1 -> 5/11,
  d2i0 -> 3/10,
  d2i1 -> 4/9
};

q1Vals = Table[N[q1Branch /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];
q2Vals = Table[N[q2Branch /. numRule /. i0 -> LL^3, 30], {LL, {10, 20, 40}}];

ratioQ1 = N[{q1Vals[[1]]/q1Vals[[2]], q1Vals[[2]]/q1Vals[[3]]}, 30];
ratioQ2 = N[{q2Vals[[1]]/q2Vals[[2]], q2Vals[[2]]/q2Vals[[3]]}, 30];

checkNumA = Min[Abs[q1Vals]] > 0 && Abs[q1Vals[[-1]]] < Abs[q1Vals[[1]]];
checkNumB = Min[Abs[q2Vals]] > 0 && Abs[q2Vals[[-1]]] < Abs[q2Vals[[1]]];
checkNumC = Max[Abs[ratioQ1 - {8, 8}]] < 1*^-2;
checkNumD = Max[Abs[ratioQ2 - {8, 8}]] < 1*^-2;

check = TrueQ[checkA && checkB && checkC && checkD && checkNumA && checkNumB && checkNumC && checkNumD];

nbName = "04_Covariant_Ratio_Second_Variation.nb";
logName = "04_covariant_ratio_second_variation.log";
logLines = {
  "Notebook: " <> nbName,
  "q1 = " <> ToString[q1, InputForm],
  "q1Expected = " <> ToString[q1Expected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "q2 = " <> ToString[q2, InputForm],
  "q2Expected = " <> ToString[q2Expected, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "q1Branch = " <> ToString[q1Branch, InputForm],
  "q1BranchExpected = " <> ToString[q1BranchExpected, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "q2Branch = " <> ToString[q2Branch, InputForm],
  "q2BranchExpected = " <> ToString[q2BranchExpected, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "q1Vals(L={10,20,40}) = " <> ToString[q1Vals, InputForm],
  "q2Vals(L={10,20,40}) = " <> ToString[q2Vals, InputForm],
  "ratioQ1 = " <> ToString[ratioQ1, InputForm],
  "ratioQ2 = " <> ToString[ratioQ2, InputForm],
  "checkNumA = " <> ToString[checkNumA, InputForm],
  "checkNumB = " <> ToString[checkNumB, InputForm],
  "checkNumC = " <> ToString[checkNumC, InputForm],
  "checkNumD = " <> ToString[checkNumD, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "04_covariant_ratio_second_variation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];

