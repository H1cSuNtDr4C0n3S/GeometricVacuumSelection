ClearAll["Global`*"];

numK2 =
  (psiP^2) * ((-4) r^2 BB^2 AAp^2 + 16 AA^3 BB psiP^2 - 8 AA^4 psiP^4 +
      4 r^2 AA BB AAp (BBp + AAp psiP^2) - AA^2 (8 BB^2 + r^2 (BBp + AAp psiP^2)^2)) +
    4 r^2 AA BB psiP (AA BBp + AAp ((-2) BB + AA psiP^2)) psiPP - 4 r^2 AA^2 BB^2 psiPP^2;

denK2 = 4 r^2 AA BB (-BB + AA psiP^2)^3;
k2General = numK2/denK2;

nb = Notebook[{
  Cell["05 - Metric Variation (Spherical)", "Title"],
  Cell["First variation structure from the general closed-form K^2 expression.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    K2general = NumK2/DenK2
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    pertRules = {
      AA -> AA + eps dAA, BB -> BB + eps dBB,
      AAp -> AAp + eps dAAp, BBp -> BBp + eps dBBp,
      psiP -> psiP + eps dpsiP, psiPP -> psiPP + eps dpsiPP
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaSeries = SeriesCoefficient[(k2General /. pertRules) // Expand, {eps, 0, 1}] // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaGateaux =
      (D[k2General, AA] dAA + D[k2General, BB] dBB +
       D[k2General, AAp] dAAp + D[k2General, BBp] dBBp +
       D[k2General, psiP] dpsiP + D[k2General, psiPP] dpsiPP) // Simplify
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkA = Simplify[
      deltaSeries == deltaGateaux,
      r != 0 && AA != 0 && BB != 0 && (-BB + AA psiP^2) != 0
    ]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    k2Minkowski =
      Simplify[
        k2General /. {AA -> 1, BB -> 1, AAp -> 0, BBp -> 0, psiP -> v, psiPP -> 0},
        r != 0 && v^2 != 1
      ]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkB = Simplify[k2Minkowski == 2 v^2/(r^2 (1 - v^2)), r != 0 && v^2 != 1]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    k2dS =
      FullSimplify[
        k2General /. {
          AA -> 1 - H^2 r^2,
          BB -> 1/(1 - H^2 r^2),
          AAp -> -2 H^2 r,
          BBp -> (2 H^2 r)/(1 - H^2 r^2)^2,
          psiP -> (H r)/(1 - H^2 r^2),
          psiPP -> H/(1 - H^2 r^2) + (2 H^3 r^2)/(1 - H^2 r^2)^2
        },
        H > 0 && r > 0 && r < 1/H
      ]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    checkC = FullSimplify[k2dS == 3 H^2, H > 0 && r > 0 && r < 1/H]
  ], "Input"]
}];

pertRules = {
  AA -> AA + eps dAA, BB -> BB + eps dBB,
  AAp -> AAp + eps dAAp, BBp -> BBp + eps dBBp,
  psiP -> psiP + eps dpsiP, psiPP -> psiPP + eps dpsiPP
};

deltaSeries = Simplify[SeriesCoefficient[Expand[k2General /. pertRules], {eps, 0, 1}]];
deltaGateaux = Simplify[
  D[k2General, AA] dAA + D[k2General, BB] dBB +
  D[k2General, AAp] dAAp + D[k2General, BBp] dBBp +
  D[k2General, psiP] dpsiP + D[k2General, psiPP] dpsiPP
];
checkAExpr = Simplify[
  deltaSeries == deltaGateaux,
  r != 0 && AA != 0 && BB != 0 && (-BB + AA psiP^2) != 0
];
checkA = TrueQ[checkAExpr];

k2Minkowski =
  Simplify[
    k2General /. {AA -> 1, BB -> 1, AAp -> 0, BBp -> 0, psiP -> v, psiPP -> 0},
    r != 0 && v^2 != 1
  ];
checkBExpr = Simplify[k2Minkowski == 2 v^2/(r^2 (1 - v^2)), r != 0 && v^2 != 1];
checkB = TrueQ[checkBExpr];

k2dS =
  FullSimplify[
    k2General /. {
      AA -> 1 - H^2 r^2,
      BB -> 1/(1 - H^2 r^2),
      AAp -> -2 H^2 r,
      BBp -> (2 H^2 r)/(1 - H^2 r^2)^2,
      psiP -> (H r)/(1 - H^2 r^2),
      psiPP -> H/(1 - H^2 r^2) + (2 H^3 r^2)/(1 - H^2 r^2)^2
    },
    H > 0 && r > 0 && r < 1/H
  ];
checkCExpr = FullSimplify[k2dS == 3 H^2, H > 0 && r > 0 && r < 1/H];
checkC = TrueQ[checkCExpr];

nbName = "05_Metric_Variation_Spherical.nb";
logName = "05_metric_variation_spherical.log";
logLines = {
  "Notebook: " <> nbName,
  "checkAExpr = " <> ToString[checkAExpr, InputForm],
  "checkA = " <> ToString[checkA, InputForm],
  "k2Minkowski = " <> ToString[k2Minkowski, InputForm],
  "checkBExpr = " <> ToString[checkBExpr, InputForm],
  "checkB = " <> ToString[checkB, InputForm],
  "k2dS = " <> ToString[k2dS, InputForm],
  "checkCExpr = " <> ToString[checkCExpr, InputForm],
  "checkC = " <> ToString[checkC, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "05_metric_variation_spherical",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "checkA" -> checkA,
  "checkB" -> checkB,
  "checkC" -> checkC,
  "checkAExpr" -> ToString[checkAExpr, InputForm],
  "checkBExpr" -> ToString[checkBExpr, InputForm],
  "checkCExpr" -> ToString[checkCExpr, InputForm]
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[checkA && checkB && checkC], Exit[0], Exit[1]];
