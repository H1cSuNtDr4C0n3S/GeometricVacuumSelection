ClearAll["Global`*"];

nb = Notebook[{
  Cell["06 - Second-Order Euler-Lagrange Identity", "Title"],
  Cell["Exact variational identity for L(x, q, q', q'').", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Lsym[x_, q0_, q1_, q2_] := c0 q0^2 + c1 q1^2 + c2 q2^2 + c3 q0 q1 + c4 q1 q2 + c5 x q0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    q[x_] := q0c + q1c x + q2c x^2 + q3c x^3 + q4c x^4
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    eta[x_] := e0 + e1 x + e2 x^2 + e3 x^3
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaL = SeriesCoefficient[Lsym[x, q[x] + eps eta[x], D[q[x] + eps eta[x], x], D[q[x] + eps eta[x], {x, 2}]], {eps, 0, 1}] // Expand
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dLdq0 = D[Lsym[x, q0, q1, q2], q0] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dLdq1 = D[Lsym[x, q0, q1, q2], q1] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    dLdq2 = D[Lsym[x, q0, q1, q2], q2] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    EL = dLdq0 - D[dLdq1, x] + D[dLdq2, {x, 2}] // Expand
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    boundary = eta[x] dLdq1 + D[eta[x], x] dLdq2 - eta[x] D[dLdq2, x]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    residue = Simplify[Expand[deltaL - (EL eta[x] + D[boundary, x])]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    check = Simplify[residue == 0]
  ], "Input"]
}];

Lsym[x_, q0_, q1_, q2_] := c0 q0^2 + c1 q1^2 + c2 q2^2 + c3 q0 q1 + c4 q1 q2 + c5 x q0;
q[x_] := q0c + q1c x + q2c x^2 + q3c x^3 + q4c x^4;
eta[x_] := e0 + e1 x + e2 x^2 + e3 x^3;

deltaL =
  Expand[
    SeriesCoefficient[
      Lsym[x, q[x] + eps eta[x], D[q[x] + eps eta[x], x], D[q[x] + eps eta[x], {x, 2}]],
      {eps, 0, 1}
    ]
  ];

dLdq0 = D[Lsym[x, q0, q1, q2], q0] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]};
dLdq1 = D[Lsym[x, q0, q1, q2], q1] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]};
dLdq2 = D[Lsym[x, q0, q1, q2], q2] /. {q0 -> q[x], q1 -> D[q[x], x], q2 -> D[q[x], {x, 2}]};

EL = Expand[dLdq0 - D[dLdq1, x] + D[dLdq2, {x, 2}]];
boundary = eta[x] dLdq1 + D[eta[x], x] dLdq2 - eta[x] D[dLdq2, x];
residue = Simplify[Expand[deltaL - (EL eta[x] + D[boundary, x])]];
check = Simplify[residue == 0];

nbName = "06_SecondOrder_EL_Identity.nb";
logName = "06_second_order_EL_identity.log";
logLines = {
  "Notebook: " <> nbName,
  "EL = " <> ToString[EL, InputForm],
  "residue = " <> ToString[residue, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "06_second_order_EL_identity",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
