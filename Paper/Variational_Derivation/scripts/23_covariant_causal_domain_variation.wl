ClearAll["Global`*"];

nb = Notebook[{
  Cell["23 - Covariant Causal-Domain Variation (4D)", "Title"],
  Cell["Full 4D covariant representation of leafwise functionals I0/I1 using delta(Theta-theta) and causal-domain indicator chi_Theta, then first-variation split.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    I0[theta] == Integrate[Sqrt[-g] W Sqrt[X] DiracDelta[Theta - theta] chiTheta, d4x]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    I1[theta] == Integrate[Sqrt[-g] W Sqrt[X] DiracDelta[Theta - theta] chiTheta K2, d4x]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI0 == Int4[chiTheta deltaMu + mu deltaChi] && deltaI1 == Int4[chiTheta K2 deltaMu + chiTheta mu deltaK2 + mu K2 deltaChi]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ == ((deltaI1) I0 - I1 (deltaI0))/I0^2
  ], "Input"]
}];

(* integrand-level perturbations in 4D covariant form *)
j0eps = Expand[(chi + eps dchi) (mu + eps dmu)];
j1eps = Expand[(chi + eps dchi) (mu + eps dmu) (k2 + eps dk2)];

dj0 = Expand[SeriesCoefficient[j0eps, {eps, 0, 1}]];
dj1 = Expand[SeriesCoefficient[j1eps, {eps, 0, 1}]];

expectedDj0 = chi dmu + mu dchi;
expectedDj1 = chi k2 dmu + chi mu dk2 + mu k2 dchi;

checkA = Simplify[dj0 == expectedDj0];
checkB = Simplify[dj1 == expectedDj1];

(* split into bulk vs moving-boundary channel: dchi -> db *)
dj0Split = Collect[dj0 /. dchi -> db, db];
dj1Split = Collect[dj1 /. dchi -> db, db];

bulk0 = Simplify[Coefficient[dj0Split, db, 0]];
boundary0 = Simplify[Coefficient[dj0Split, db, 1]];
bulk1 = Simplify[Coefficient[dj1Split, db, 0]];
boundary1 = Simplify[Coefficient[dj1Split, db, 1]];

checkC = Simplify[bulk0 == chi dmu];
checkD = Simplify[boundary0 == mu];
checkE = Simplify[bulk1 == chi (k2 dmu + mu dk2)];
checkF = Simplify[boundary1 == mu k2];

(* integrated proxy form for Q-variation; xi tracks boundary displacement amplitude *)
i0 = Symbol["I0"];
i1 = Symbol["I1"];
di0bulk = Symbol["dI0bulk"];
di1bulk = Symbol["dI1bulk"];
b0 = Symbol["b0"];
b1 = Symbol["b1"];

dI0 = di0bulk + xi b0;
dI1 = di1bulk + xi b1;
deltaQ = Simplify[(dI1 i0 - i1 dI0)/i0^2];
boundaryQ = Simplify[Coefficient[deltaQ, xi, 1]];
expectedBoundaryQ = Simplify[(b1 i0 - i1 b0)/i0^2];
checkG = Simplify[boundaryQ == expectedBoundaryQ];

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG];

nbName = "23_Covariant_Causal_Domain_Variation.nb";
logName = "23_covariant_causal_domain_variation.log";
logLines = {
  "Notebook: " <> nbName,
  "dj0 = " <> ToString[dj0, InputForm],
  "expectedDj0 = " <> ToString[expectedDj0, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "dj1 = " <> ToString[dj1, InputForm],
  "expectedDj1 = " <> ToString[expectedDj1, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "bulk0 = " <> ToString[bulk0, InputForm],
  "boundary0 = " <> ToString[boundary0, InputForm],
  "bulk1 = " <> ToString[bulk1, InputForm],
  "boundary1 = " <> ToString[boundary1, InputForm],
  "checkC = " <> ToString[checkC, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "checkE = " <> ToString[checkE, InputForm],
  "checkF = " <> ToString[checkF, InputForm],
  "deltaQ = " <> ToString[deltaQ, InputForm],
  "boundaryQ = " <> ToString[boundaryQ, InputForm],
  "expectedBoundaryQ = " <> ToString[expectedBoundaryQ, InputForm],
  "checkG = " <> ToString[checkG, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "23_covariant_causal_domain_variation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];

