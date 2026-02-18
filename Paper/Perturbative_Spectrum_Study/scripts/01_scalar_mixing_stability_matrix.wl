ClearAll["Global`*"];

nb = Notebook[{
  Cell["01 - Scalar Mixing Stability Matrix", "Title"],
  Cell["General 2-field scalar quadratic system: kinetic/gradient matrices, ghost criteria, and high-k sound-speed eigenvalues.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Kmat == {{k11, k12}, {k12, k22}} && Gmat == {{g11, g12}, {g12, g22}}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    ghostFree == (k11 > 0 && Det[Kmat] > 0)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    c2eig == Eigenvalues[Inverse[Kmat].Gmat]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    gradientStable == (Tr[Inverse[Kmat].Gmat] > 0 && Det[Inverse[Kmat].Gmat] > 0)
  ], "Input"]
}];

Kmat = {{k11, k12}, {k12, k22}};
Gmat = {{g11, g12}, {g12, g22}};

detK = FullSimplify[Det[Kmat]];
traceKG = FullSimplify[Tr[LinearSolve[Kmat, Gmat]], Assumptions -> detK != 0];
detKG = FullSimplify[Det[Gmat]/detK, Assumptions -> detK != 0];

traceExpected = FullSimplify[(g11 k22 + g22 k11 - 2 g12 k12)/detK, Assumptions -> detK != 0];
detExpected = FullSimplify[(g11 g22 - g12^2)/detK, Assumptions -> detK != 0];

checkA = FullSimplify[traceKG == traceExpected, Assumptions -> detK != 0];
checkB = FullSimplify[detKG == detExpected, Assumptions -> detK != 0];

disc = FullSimplify[traceKG^2 - 4 detKG, Assumptions -> detK != 0];
c2plus = FullSimplify[(traceKG + Sqrt[disc])/2, Assumptions -> detK != 0];
c2minus = FullSimplify[(traceKG - Sqrt[disc])/2, Assumptions -> detK != 0];

stableRule = {
  k11 -> 2, k12 -> 3/10, k22 -> 11/10,
  g11 -> 3/2, g12 -> 1/5, g22 -> 9/10
};

Kstable = Kmat /. stableRule;
Gstable = Gmat /. stableRule;
eigStable = N[Eigenvalues[LinearSolve[Kstable, Gstable]], 30];
checkNumA = Min[Eigenvalues[Kstable]] > 0;
checkNumB = Min[eigStable] > 0;

unstableRule = {
  k11 -> 2, k12 -> 3/10, k22 -> 11/10,
  g11 -> 3/2, g12 -> 1/5, g22 -> -2/5
};
Kunstable = Kmat /. unstableRule;
Gunstable = Gmat /. unstableRule;
eigUnstable = N[Eigenvalues[LinearSolve[Kunstable, Gunstable]], 30];
checkNumC = Min[eigUnstable] < 0;

check = TrueQ[checkA && checkB && checkNumA && checkNumB && checkNumC];

nbName = "01_Scalar_Mixing_Stability_Matrix.nb";
logName = "01_scalar_mixing_stability_matrix.log";
logLines = {
  "Notebook: " <> nbName,
  "detK = " <> ToString[detK, InputForm],
  "traceKG = " <> ToString[traceKG, InputForm],
  "detKG = " <> ToString[detKG, InputForm],
  "traceExpected = " <> ToString[traceExpected, InputForm],
  "detExpected = " <> ToString[detExpected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "disc = " <> ToString[disc, InputForm],
  "c2plus = " <> ToString[c2plus, InputForm],
  "c2minus = " <> ToString[c2minus, InputForm],
  "eigStable = " <> ToString[eigStable, InputForm],
  "checkNumA(K positive definite) = " <> ToString[checkNumA, InputForm],
  "checkNumB(stable eigenvalues positive) = " <> ToString[checkNumB, InputForm],
  "eigUnstable = " <> ToString[eigUnstable, InputForm],
  "checkNumC(unstable sample has negative c^2) = " <> ToString[checkNumC, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "01_scalar_mixing_stability_matrix",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];

