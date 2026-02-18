ClearAll["Global`*"];

nb = Notebook[{
  Cell["00 - Variational Seed", "Title"],
  Cell["Primo controllo simbolico per il termine con moltiplicatore di Lagrange leafwise.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[lagDensity[th_] := lambda[th] (I1[th] - K0^2 I0[th])], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[eqLambda = D[lagDensity[th], lambda[th]] // Simplify], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[check1 = Simplify[eqLambda == I1[th] - K0^2 I0[th]]], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[perturbed = Expand[(lambda[th] + eps eta[th]) (I1[th] - K0^2 I0[th])]], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[lin = SeriesCoefficient[perturbed, {eps, 0, 1}]], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[check2 = Simplify[lin == eta[th] (I1[th] - K0^2 I0[th])]], "Input"]
}];

lagDensity[th_] := lambda[th] (I1[th] - K0^2 I0[th]);
eqLambda = Simplify[D[lagDensity[th], lambda[th]]];
check1 = Simplify[eqLambda == I1[th] - K0^2 I0[th]];

Clear[eps];
perturbed = Expand[(lambda[th] + eps eta[th]) (I1[th] - K0^2 I0[th])];
lin = SeriesCoefficient[perturbed, {eps, 0, 1}];
check2 = Simplify[lin == eta[th] (I1[th] - K0^2 I0[th])];

nbName = "00_Variational_Seed.nb";
logName = "00_variational_seed.log";
logLines = {
  "Notebook: " <> nbName,
  "eqLambda = " <> ToString[eqLambda, InputForm],
  "check1 = " <> ToString[check1, InputForm],
  "lin = " <> ToString[lin, InputForm],
  "check2 = " <> ToString[check2, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "00_variational_seed",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check1" -> check1,
  "check2" -> check2
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check1 && check2], Exit[0], Exit[1]];
