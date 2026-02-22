ClearAll["Global`*"];

nb = Notebook[{
  Cell["07 - Domain Proxy IR Class Equivalence Test", "Title"],
  Cell["Objective: compare two reasonable causal-domain proxies (redshift-threshold vs apparent-horizon) and verify they lead to the same IR behavior class, showing no fine-tuning.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    DThetaRed[M_, H_, z2_] == {r > 0 | f[r, M, H] >= z2}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    DThetaApp[M_, H_] == {r > 0 | f[r, M, H] >= 0}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Q[M_, H_, proxy_] == I1[M, H, proxy]/I0[M, H, proxy]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ[H_, proxy_] == Q[M, H, proxy] - Q[0, H, proxy]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    proxyClass[proxy] == "IR_SUPPRESSED_SUBLEADING" && proxyClass["redshift"] == proxyClass["apparent_horizon"]
  ], "Input"]
}];

f[r_, m_, h_] := 1 - 2 m/r - h^2 r^2;

positiveRoots[m_?NumericQ, h_?NumericQ, z2_?NumericQ] := Module[{sol},
  sol = r /. NSolve[h^2 r^3 - (1 - z2) r + 2 m == 0, r, Reals, WorkingPrecision -> 80];
  sol = Chop[N[sol, 40], 10^-25];
  sol = Select[sol, NumericQ[#] && # > 0 &];
  Sort[DeleteDuplicates[sol, Abs[#1 - #2] < 10^-20 &]]
];

horizonRoots[m_?NumericQ, h_?NumericQ] := Module[{sol},
  sol = r /. NSolve[f[r, m, h] == 0, r, Reals, WorkingPrecision -> 80];
  sol = Chop[N[sol, 40], 10^-25];
  sol = Select[sol, NumericQ[#] && # > 0 &];
  Sort[DeleteDuplicates[sol, Abs[#1 - #2] < 10^-20 &]]
];

domainBounds[proxy_String, m_?NumericQ, h_?NumericQ, z2_?NumericQ] := Module[{roots},
  Which[
    proxy == "redshift_threshold",
    If[m <= 10^-14,
      {0, N[Sqrt[1 - z2]/h, 50]},
      roots = positiveRoots[m, h, z2];
      If[Length[roots] < 2, $Failed, {First[roots], Last[roots]}]
    ],
    proxy == "apparent_horizon",
    If[m <= 10^-14,
      {0, N[1/h, 50]},
      roots = horizonRoots[m, h];
      If[Length[roots] < 2, $Failed, {First[roots], Last[roots]}]
    ],
    True,
    $Failed
  ]
];

i0Domain[bounds_List] := N[(4 Pi/3) (bounds[[2]]^3 - bounds[[1]]^3), 50];
k2bg[h_?NumericQ] := N[3 h^2, 50];

localI1[m_?NumericQ, bounds_List, alpha_?NumericQ, beta_?NumericQ, rFloor_?NumericQ] := Module[
  {rs, a, ff},
  If[m <= 10^-14, Return[0]];
  rs = 2 m;
  a = N[Sqrt[(beta rs)^2 + rFloor^2], 50];
  ff[x_] := N[(1/(2 a)) ArcTan[x/a] - x/(2 (x^2 + a^2)), 50];
  N[4 Pi alpha rs^2 (ff[bounds[[2]]] - ff[bounds[[1]]]), 50]
];

qValue[m_?NumericQ, h_?NumericQ, bounds_List, alpha_?NumericQ, beta_?NumericQ, rFloor_?NumericQ] := Module[
  {i0, loc},
  i0 = i0Domain[bounds];
  loc = localI1[m, bounds, alpha, beta, rFloor];
  N[k2bg[h] + loc/i0, 50]
];

evaluateProxy[proxy_String, mBH_?NumericQ, hRef_?NumericQ, hHalf_?NumericQ, z2_?NumericQ, alpha_?NumericQ, beta_?NumericQ, rFloor_?NumericQ] := Module[
  {
    bRef, bHalf, b0Ref, b0Half, domainsOK,
    i0Ref, i0Half, locRef, locHalf,
    qMassRef, qMassHalf, q0Ref, q0Half,
    deltaQRef, deltaQHalf, ratioDelta, ratioI0, ratioLoc, epsRef, epsHalf,
    checkA, checkB, checkC, checkD, checkE, checkF, checkG, checkProxy, classification
  },

  bRef = domainBounds[proxy, mBH, hRef, z2];
  bHalf = domainBounds[proxy, mBH, hHalf, z2];
  b0Ref = domainBounds[proxy, 0, hRef, z2];
  b0Half = domainBounds[proxy, 0, hHalf, z2];
  domainsOK = bRef =!= $Failed && bHalf =!= $Failed && b0Ref =!= $Failed && b0Half =!= $Failed;

  i0Ref = Indeterminate;
  i0Half = Indeterminate;
  locRef = Indeterminate;
  locHalf = Indeterminate;
  qMassRef = Indeterminate;
  qMassHalf = Indeterminate;
  q0Ref = Indeterminate;
  q0Half = Indeterminate;
  deltaQRef = Indeterminate;
  deltaQHalf = Indeterminate;
  ratioDelta = Indeterminate;
  ratioI0 = Indeterminate;
  ratioLoc = Indeterminate;
  epsRef = Indeterminate;
  epsHalf = Indeterminate;

  checkA = False;
  checkB = False;
  checkC = False;
  checkD = False;
  checkE = False;
  checkF = False;
  checkG = False;
  checkProxy = False;
  classification = "undetermined";

  If[TrueQ[domainsOK],
    i0Ref = i0Domain[bRef];
    i0Half = i0Domain[bHalf];
    locRef = localI1[mBH, bRef, alpha, beta, rFloor];
    locHalf = localI1[mBH, bHalf, alpha, beta, rFloor];

    qMassRef = qValue[mBH, hRef, bRef, alpha, beta, rFloor];
    qMassHalf = qValue[mBH, hHalf, bHalf, alpha, beta, rFloor];
    q0Ref = qValue[0, hRef, b0Ref, alpha, beta, rFloor];
    q0Half = qValue[0, hHalf, b0Half, alpha, beta, rFloor];

    deltaQRef = N[qMassRef - q0Ref, 50];
    deltaQHalf = N[qMassHalf - q0Half, 50];
    ratioDelta = N[deltaQHalf/deltaQRef, 40];
    ratioI0 = N[i0Half/i0Ref, 40];
    ratioLoc = N[locHalf/locRef, 40];
    epsRef = N[deltaQRef/q0Ref, 40];
    epsHalf = N[deltaQHalf/q0Half, 40];

    checkA = Abs[q0Ref - k2bg[hRef]] < 10^-16 && Abs[q0Half - k2bg[hHalf]] < 10^-16;
    checkB = deltaQRef > 0 && deltaQHalf > 0;
    checkC = Abs[ratioDelta/(1/8) - 1] < 0.2;
    checkD = Abs[ratioI0/8 - 1] < 0.2;
    checkE = Abs[ratioLoc - 1] < 0.1;
    checkF = epsRef < 0.05 && epsHalf < 0.05;
    checkG = deltaQHalf < deltaQRef/4;

    checkProxy = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG];
    classification = If[TrueQ[checkProxy], "IR_SUPPRESSED_SUBLEADING", "undetermined"];
  ];

  <|
    "proxy" -> proxy,
    "boundsRef" -> bRef,
    "boundsHalf" -> bHalf,
    "bounds0Ref" -> b0Ref,
    "bounds0Half" -> b0Half,
    "i0Ref" -> i0Ref,
    "i0Half" -> i0Half,
    "locRef" -> locRef,
    "locHalf" -> locHalf,
    "qMassRef" -> qMassRef,
    "qMassHalf" -> qMassHalf,
    "q0Ref" -> q0Ref,
    "q0Half" -> q0Half,
    "deltaQRef" -> deltaQRef,
    "deltaQHalf" -> deltaQHalf,
    "ratioDelta" -> ratioDelta,
    "ratioI0" -> ratioI0,
    "ratioLoc" -> ratioLoc,
    "epsRef" -> epsRef,
    "epsHalf" -> epsHalf,
    "checkA" -> checkA,
    "checkB" -> checkB,
    "checkC" -> checkC,
    "checkD" -> checkD,
    "checkE" -> checkE,
    "checkF" -> checkF,
    "checkG" -> checkG,
    "checkProxy" -> checkProxy,
    "classification" -> classification
  |>
];

mBH = 1;
hRef = 1/100;
hHalf = hRef/2;
z2Ref = 1/20;
alpha = 1;
beta = 1;
rFloor = 10^-12;

redRes = evaluateProxy["redshift_threshold", mBH, hRef, hHalf, z2Ref, alpha, beta, rFloor];
appRes = evaluateProxy["apparent_horizon", mBH, hRef, hHalf, z2Ref, alpha, beta, rFloor];

z2Scan = N[{1/100, 1/50, 1/20, 1/10}, 50];
scanData = Table[
  Module[{res},
    res = evaluateProxy["redshift_threshold", mBH, hRef, hHalf, zz, alpha, beta, rFloor];
    <|
      "z2" -> zz,
      "classification" -> res["classification"],
      "ratioDelta" -> res["ratioDelta"],
      "epsRef" -> res["epsRef"],
      "epsHalf" -> res["epsHalf"],
      "checkProxy" -> res["checkProxy"]
    |>
  ],
  {zz, z2Scan}
];

scanClasses = Lookup[scanData, "classification"];
scanRatios = Lookup[scanData, "ratioDelta"];
scanEpsRef = Lookup[scanData, "epsRef"];
scanEpsHalf = Lookup[scanData, "epsHalf"];
scanChecks = Lookup[scanData, "checkProxy"];

checkA = TrueQ[redRes["checkProxy"]];
checkB = TrueQ[appRes["checkProxy"]];
checkC = redRes["classification"] == appRes["classification"] && redRes["classification"] == "IR_SUPPRESSED_SUBLEADING";
checkD = Abs[redRes["ratioDelta"]/appRes["ratioDelta"] - 1] < 0.02;
checkE = Abs[redRes["epsRef"] - appRes["epsRef"]] < 0.01 && Abs[redRes["epsHalf"] - appRes["epsHalf"]] < 0.01;
checkF = redRes["boundsRef"] =!= appRes["boundsRef"] && redRes["boundsHalf"] =!= appRes["boundsHalf"];
checkG = AllTrue[scanChecks, TrueQ];
checkH = AllTrue[scanClasses, # == "IR_SUPPRESSED_SUBLEADING" &];
checkI = Max[Abs[scanRatios - 1/8]] < 5*10^-3;
checkJ = Max[scanEpsRef] < 0.05 && Max[scanEpsHalf] < 0.05;

classification = If[TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI && checkJ],
  "PROXY_CLASS_EQUIVALENT_IR",
  "undetermined"
];

check = TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI && checkJ &&
   classification == "PROXY_CLASS_EQUIVALENT_IR"];

nbName = "07_Domain_Proxy_IR_Class_Equivalence_Test.nb";
logName = "07_domain_proxy_ir_class_equivalence_test.log";
logLines = {
  "Notebook: " <> nbName,
  "mBH = " <> ToString[mBH, InputForm],
  "hRef = " <> ToString[hRef, InputForm],
  "hHalf = " <> ToString[hHalf, InputForm],
  "z2Ref = " <> ToString[z2Ref, InputForm],
  "alpha = " <> ToString[alpha, InputForm],
  "beta = " <> ToString[beta, InputForm],
  "rFloor = " <> ToString[rFloor, InputForm],
  "redRes = " <> ToString[redRes, InputForm],
  "appRes = " <> ToString[appRes, InputForm],
  "z2Scan = " <> ToString[z2Scan, InputForm],
  "scanData = " <> ToString[scanData, InputForm],
  "scanClasses = " <> ToString[scanClasses, InputForm],
  "scanRatios = " <> ToString[scanRatios, InputForm],
  "scanEpsRef = " <> ToString[scanEpsRef, InputForm],
  "scanEpsHalf = " <> ToString[scanEpsHalf, InputForm],
  "checkA(red proxy IR class) = " <> ToString[checkA, InputForm],
  "checkB(apparent-horizon proxy IR class) = " <> ToString[checkB, InputForm],
  "checkC(class equality across proxies) = " <> ToString[checkC, InputForm],
  "checkD(ratioDelta agreement across proxies) = " <> ToString[checkD, InputForm],
  "checkE(subleading eps agreement) = " <> ToString[checkE, InputForm],
  "checkF(domains are genuinely different) = " <> ToString[checkF, InputForm],
  "checkG(redshift scan passes all proxy checks) = " <> ToString[checkG, InputForm],
  "checkH(redshift scan keeps same class) = " <> ToString[checkH, InputForm],
  "checkI(redshift scan ratioDelta remains ~1/8) = " <> ToString[checkI, InputForm],
  "checkJ(redshift scan remains subleading) = " <> ToString[checkJ, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "07_domain_proxy_ir_class_equivalence_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
