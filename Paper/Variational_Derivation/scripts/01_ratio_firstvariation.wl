ClearAll["Global`*"];

nb = Notebook[{
  Cell["01 - Ratio First Variation", "Title"],
  Cell["Linearizzazione simbolica della media normalizzata <K^2> = I1/I0.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[ratioPerturbed = (i1 + eps di1)/(i0 + eps di0)], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[series1 = Normal @ Series[ratioPerturbed, {eps, 0, 1}] // Simplify], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[deltaRatio = Coefficient[series1 - i1/i0, eps] // Simplify], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[expected = di1/i0 - (i1 di0)/i0^2], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[check = Simplify[deltaRatio == expected, i0 != 0]], "Input"]
}];

$Assumptions = i0 != 0;
ratioPerturbed = (i1 + eps di1)/(i0 + eps di0);
series1 = Simplify[Normal @ Series[ratioPerturbed, {eps, 0, 1}]];
deltaRatio = Simplify[Coefficient[series1 - i1/i0, eps]];
expected = Simplify[di1/i0 - (i1 di0)/i0^2];
check = Simplify[deltaRatio == expected];

nbName = "01_Ratio_FirstVariation.nb";
logName = "01_ratio_firstvariation.log";
logLines = {
  "Notebook: " <> nbName,
  "series1 = " <> ToString[series1, InputForm],
  "deltaRatio = " <> ToString[deltaRatio, InputForm],
  "expected = " <> ToString[expected, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "01_ratio_firstvariation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
