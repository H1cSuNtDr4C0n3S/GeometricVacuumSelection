ClearAll["Global`*"];

nb = Notebook[{
  Cell["25 - Covariant Nonlocal Action Variation", "Title"],
  Cell["Leafwise nonlocal action variation in the paper form S_NL = Integral dTheta lambda(Theta) [I1 - K0^2 I0], with full bulk+boundary split and IR suppression in the normalized-average channel.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaSleaf == eta (I1 - K0^2 I0) + lambda (deltaI1 - K0^2 deltaI0)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaI0 == A0 + B0 && deltaI1 == A1 + C1 + B1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaC == deltaI1 - K0^2 deltaI0 == (A1 - K0^2 A0 + C1) + (B1 - K0^2 B0)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ == (deltaI1 I0 - I1 deltaI0)/I0^2
  ], "Input"]
}];

(* single-leaf variation of nonlocal action term *)
sEps = Expand[(lam + eps eta) ((i1 + eps di1) - k0^2 (i0 + eps di0))];
deltaSleaf = Expand[SeriesCoefficient[sEps, {eps, 0, 1}]];
expectedDeltaSleaf = Expand[eta (i1 - k0^2 i0) + lam (di1 - k0^2 di0)];
checkA = Simplify[deltaSleaf == expectedDeltaSleaf];

(* covariant split of dI0,dI1:
   A0 = Int4[chi dmu], B0 = Int4[mu dchi]
   A1 = Int4[chi K2 dmu], C1 = Int4[chi mu dK2], B1 = Int4[mu K2 dchi] *)
A0 = Symbol["A0"];
B0 = Symbol["B0"];
A1 = Symbol["A1"];
C1 = Symbol["C1"];
B1 = Symbol["B1"];

di0Cov = A0 + B0;
di1Cov = A1 + C1 + B1;

deltaC = Expand[di1Cov - k0^2 di0Cov];
bulkPart = Expand[A1 - k0^2 A0 + C1];
boundaryPart = Expand[B1 - k0^2 B0];

checkB = Simplify[deltaC == bulkPart + boundaryPart];
checkC = Simplify[(boundaryPart /. B1 -> k0^2 B0) == 0];

deltaSgeom = Expand[lam deltaC];
deltaSfull = Expand[eta (i1 - k0^2 i0) + deltaSgeom];

(* normalized-average channel (used in variational IR argument in the paper) *)
deltaQ = Simplify[(di1Cov i0 - i1 di0Cov)/i0^2];
deltaQIR = Simplify[deltaQ /. i1 -> q i0];
checkD = Simplify[deltaQIR == (di1Cov - q di0Cov)/i0];

(* explicit scaling test: i0 ~ L^3 *)
numSub = {
  A0 -> 3/4,
  B0 -> 2/5,
  A1 -> 7/5,
  C1 -> 1/3,
  B1 -> 1/6,
  q -> 3/10
};

vals = Table[N[Abs[deltaQIR /. numSub /. i0 -> LL^3], 30], {LL, {10, 20, 40}}];
ratios = N[{vals[[1]]/vals[[2]], vals[[2]]/vals[[3]]}, 30];
checkE = Max[Abs[ratios - {8, 8}]] < 10^-10;

check = TrueQ[checkA && checkB && checkC && checkD && checkE];

nbName = "25_Covariant_Nonlocal_Action_Variation.nb";
logName = "25_covariant_nonlocal_action_variation.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaSleaf = " <> ToString[deltaSleaf, InputForm],
  "expectedDeltaSleaf = " <> ToString[expectedDeltaSleaf, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "di0Cov = " <> ToString[di0Cov, InputForm],
  "di1Cov = " <> ToString[di1Cov, InputForm],
  "deltaC = " <> ToString[deltaC, InputForm],
  "bulkPart = " <> ToString[bulkPart, InputForm],
  "boundaryPart = " <> ToString[boundaryPart, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "checkC(boundary neutral when B1=K0^2 B0) = " <> ToString[checkC, InputForm],
  "deltaSgeom = " <> ToString[deltaSgeom, InputForm],
  "deltaSfull = " <> ToString[deltaSfull, InputForm],
  "deltaQ = " <> ToString[deltaQ, InputForm],
  "deltaQIR = " <> ToString[deltaQIR, InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "vals(L={10,20,40}) = " <> ToString[vals, InputForm],
  "ratios = " <> ToString[ratios, InputForm],
  "checkE(IR decay) = " <> ToString[checkE, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "25_covariant_nonlocal_action_variation",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];

