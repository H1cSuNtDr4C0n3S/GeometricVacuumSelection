ClearAll["Global`*"];

nb = Notebook[{
  Cell["04 - Ratio Kernel Form", "Title"],
  Cell["Kernel form of the first variation of q = I1/I0 under metric perturbations.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q = I1/I0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ = dI1/I0 - (I1 dI0)/I0^2 // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    rhs = (1/I0) dI1 - (q/I0) dI0 // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkA = Simplify[deltaQ == rhs]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dI0 = Inactive[Integrate][k0[x] dg[x], {x, xmin, xmax}]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dI1 = Inactive[Integrate][k1[x] dg[x], {x, xmin, xmax}]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    rhsKernel = rhs /. q -> I1/I0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    expectedKernelExpanded = (1/I0) Inactive[Integrate][k1[x] dg[x], {x, xmin, xmax}] -
      (I1/I0^2) Inactive[Integrate][k0[x] dg[x], {x, xmin, xmax}]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkB = Simplify[rhsKernel == expectedKernelExpanded]
  ], "Input"]
}];

q = I1/I0;
deltaQ = Simplify[dI1/I0 - (I1 dI0)/I0^2];
rhs = Simplify[(1/I0) dI1 - (q/I0) dI0];
checkA = Simplify[deltaQ == rhs];

dI0 = Inactive[Integrate][k0[x] dg[x], {x, xmin, xmax}];
dI1 = Inactive[Integrate][k1[x] dg[x], {x, xmin, xmax}];
rhsKernel = rhs /. q -> I1/I0;
expectedKernelExpanded =
  (1/I0) Inactive[Integrate][k1[x] dg[x], {x, xmin, xmax}] -
  (I1/I0^2) Inactive[Integrate][k0[x] dg[x], {x, xmin, xmax}];
checkB = Simplify[rhsKernel == expectedKernelExpanded];

nbName = "04_Ratio_Kernel_Form.nb";
logName = "04_ratio_kernel_form.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaQ = " <> ToString[deltaQ, InputForm],
  "rhs = " <> ToString[rhs, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "rhsKernel = " <> ToString[rhsKernel, InputForm],
  "expectedKernelExpanded = " <> ToString[expectedKernelExpanded, InputForm],
  "checkB = " <> ToString[checkB, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "04_ratio_kernel_form",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "checkA" -> checkA,
  "checkB" -> checkB
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[checkA && checkB], Exit[0], Exit[1]];
