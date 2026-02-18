ClearAll["Global`*"];

nb = Notebook[{
  Cell["00 - Tensor Channel IR Stability", "Title"],
  Cell["Quadratic tensor channel around FRW/de Sitter: ghost/gradient conditions and IR-suppressed nonlocal corrections consistent with the paper.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    ST2 == (1/8) Integrate[a^3 (GT hdot^2 - FT grad[h]^2/a^2), d4x]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    GT == 1/kappa + dGT && FT == 1/kappa + dFT
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    cT2 == FT/GT && ghostCondition == (GT > 0) && gradientCondition == (FT > 0)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dGT == aGT/I0 && dFT == aFT/I0 && cT2 -> 1 as I0 -> Infinity
  ], "Input"]
}];

GT = 1/kappa + dGT;
FT = 1/kappa + dFT;
cT2 = FullSimplify[FT/GT];

checkA = FullSimplify[(cT2 /. {dGT -> 0, dFT -> 0}) == 1];

irRule = {dGT -> aGT/I0, dFT -> aFT/I0};
cT2IR = FullSimplify[cT2 /. irRule];
deltaCT = FullSimplify[cT2IR - 1];
deltaAsym = FullSimplify[Limit[I0 deltaCT, I0 -> Infinity], Assumptions -> kappa != 0];
checkB = FullSimplify[deltaAsym == kappa (aFT - aGT), Assumptions -> kappa != 0];

ghostCondition = GT > 0;
gradientCondition = FT > 0;

numVals = Table[
  N[
    ({GT, FT} /. irRule)~Join~{cT2IR} /. {
      kappa -> 2,
      aGT -> 3/50,
      aFT -> 2/25,
      I0 -> LL
    },
    30
  ],
  {LL, {10, 20, 40}}
];

checkNumA = Min[numVals[[All, 1]]] > 0;
checkNumB = Min[numVals[[All, 2]]] > 0;
checkNumC = Abs[numVals[[-1, 3]] - 1] < Abs[numVals[[1, 3]] - 1];

check = TrueQ[checkA && checkB && checkNumA && checkNumB && checkNumC];

nbName = "00_Tensor_Channel_IR_Stability.nb";
logName = "00_tensor_channel_ir_stability.log";
logLines = {
  "Notebook: " <> nbName,
  "GT = " <> ToString[GT, InputForm],
  "FT = " <> ToString[FT, InputForm],
  "cT2 = " <> ToString[cT2, InputForm],
  "checkA(GR limit) = " <> ToString[checkA, InputForm],
  "cT2IR = " <> ToString[cT2IR, InputForm],
  "deltaCT = " <> ToString[deltaCT, InputForm],
  "deltaAsym = " <> ToString[deltaAsym, InputForm],
  "checkB(IR asymptotic) = " <> ToString[checkB, InputForm],
  "ghostCondition = " <> ToString[ghostCondition, InputForm],
  "gradientCondition = " <> ToString[gradientCondition, InputForm],
  "numVals(I0={10,20,40}) = " <> ToString[numVals, InputForm],
  "checkNumA(GT>0) = " <> ToString[checkNumA, InputForm],
  "checkNumB(FT>0) = " <> ToString[checkNumB, InputForm],
  "checkNumC(cT2->1) = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "00_tensor_channel_ir_stability",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
