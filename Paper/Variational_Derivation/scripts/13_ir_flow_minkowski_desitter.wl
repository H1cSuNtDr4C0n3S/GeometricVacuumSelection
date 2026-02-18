ClearAll["Global`*"];

nb = Notebook[{
  Cell["13 - IR Flow: Minkowski vs de Sitter", "Title"],
  Cell["Finite-domain flow dQ/dR for Q(R)=I1(R)/I0(R) and its contrasting behavior in Minkowski/de Sitter toy sectors.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dQ/dR = (f1(R) I0(R) - I1(R) f0(R))/I0(R)^2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Q_M(R) = 3 c/R^2,   dQ_M/dR = -6 c/R^3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    f1 = KdS f0  ->  Q_dS(R)=KdS, dQ_dS/dR=0
  ], "Input"]
}];

(* general flow identity *)
i0 = Symbol["I0"];
i1 = Symbol["I1"];
f0R = Symbol["f0R"];
f1R = Symbol["f1R"];
qR = i1/i0;
qPrime = Simplify[D[(i1 + eps f1R)/(i0 + eps f0R), eps] /. eps -> 0];
flowExpected = Simplify[(f1R i0 - i1 f0R)/i0^2];
checkA = Simplify[qPrime == flowExpected];

(* Minkowski toy: K^2 ~ c/r^2, volume kernel ~ r^2 -> numerator kernel ~ c *)
f0M[r_] := r^2;
f1M[r_] := c;
i0M[R_] := Integrate[f0M[r], {r, 0, R}, Assumptions -> R > 0];
i1M[R_] := Integrate[f1M[r], {r, 0, R}, Assumptions -> R > 0];
qM[R_] := Simplify[i1M[R]/i0M[R], R > 0];
dqM[R_] := Simplify[D[qM[R], R], R > 0];
flowM[R_] := Simplify[(f1M[R] i0M[R] - i1M[R] f0M[R])/i0M[R]^2, R > 0];
checkM = Simplify[qM[R] == 3 c/R^2 && dqM[R] == -6 c/R^3 && dqM[R] == flowM[R], R > 0];

(* de Sitter toy: constant ratio sector f1 = KdS f0 *)
f0D[r_] := 1 + r + r^2;
f1D[r_] := KdS f0D[r];
i0D[R_] := Integrate[f0D[r], {r, 0, R}, Assumptions -> R > 0];
i1D[R_] := Integrate[f1D[r], {r, 0, R}, Assumptions -> R > 0];
qD[R_] := Simplify[i1D[R]/i0D[R], R > 0];
dqD[R_] := Simplify[D[qD[R], R], R > 0];
flowD[R_] := Simplify[(f1D[R] i0D[R] - i1D[R] f0D[R])/i0D[R]^2, R > 0];
checkD = Simplify[qD[R] == KdS && dqD[R] == 0 && dqD[R] == flowD[R], R > 0];

(* numerical sanity *)
c0 = 0.8;
kd0 = 1.2;
sampleRs = {2, 4, 8};
qMvals = N[Table[qM[rr] /. c -> c0, {rr, sampleRs}], 20];
qDvals = N[Table[qD[rr] /. KdS -> kd0, {rr, sampleRs}], 20];

check = TrueQ[checkA && checkM && checkD];

nbName = "13_IR_Flow_Minkowski_DeSitter.nb";
logName = "13_ir_flow_minkowski_desitter.log";
logLines = {
  "Notebook: " <> nbName,
  "qPrime = " <> ToString[qPrime, InputForm],
  "flowExpected = " <> ToString[flowExpected, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "qM(R) = " <> ToString[qM[R], InputForm],
  "dqM(R) = " <> ToString[dqM[R], InputForm],
  "checkM = " <> ToString[checkM, InputForm],
  "qD(R) = " <> ToString[qD[R], InputForm],
  "dqD(R) = " <> ToString[dqD[R], InputForm],
  "checkD = " <> ToString[checkD, InputForm],
  "sampleRs = " <> ToString[sampleRs, InputForm],
  "qMvals(c=0.8) = " <> ToString[qMvals, InputForm],
  "qDvals(KdS=1.2) = " <> ToString[qDvals, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "13_ir_flow_minkowski_desitter",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
