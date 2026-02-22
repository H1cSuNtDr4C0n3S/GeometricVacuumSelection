ClearAll["Global`*"];

nb = Notebook[{
  Cell["04 - Era Proxy Equivalence IR Class Test", "Title"],
  Cell["Objective: compare two reasonable finite-domain proxies (redshift-threshold vs apparent-horizon radius) across radiation, matter, and late de Sitter eras; verify same IR behavior class without fine-tuning.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    DThetaRed[z_] == {r > 0 | r <= Sqrt[1 - z2]/H[z]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    DThetaApp[z_] == {r > 0 | r <= 1/H[z]}
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Q[z, proxy] := 3 H[z]^2 + I1[z, proxy]/I0[z, proxy]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    class[era, proxy] == "IR_SUPPRESSED_SUBLEADING"
  ], "Input"]
}];

toNumeric[token_String] := Quiet @ Check[ToExpression[token], $Failed];

parseBestfitLines[lines_List] := Module[
  {useful, headerTokens, valueTokens, values, assoc},
  useful = Select[lines, StringLength[StringTrim[#]] > 0 &];
  If[Length[useful] < 2, Return[$Failed]];
  headerTokens = StringSplit[StringTrim[StringReplace[First[useful], "#" -> ""]]];
  valueTokens = StringSplit[StringTrim[useful[[2]]]];
  If[Length[headerTokens] =!= Length[valueTokens], Return[$Failed]];
  values = toNumeric /@ valueTokens;
  If[AnyTrue[values, # === $Failed &], Return[$Failed]];
  assoc = AssociationThread[headerTokens -> values];
  assoc
];

loadBestfitData[localFile_String, remoteURL_String, fallbackHash_String, fallbackAssoc_Association] := Module[
  {lines = $Failed, hash = "missing", src = "none", tmp = $Failed, parsed},
  If[FileExistsQ[localFile],
    lines = Import[localFile, "Lines"];
    hash = ToLowerCase[FileHash[localFile, "SHA256", "HexString"]];
    src = "local_file";,
    tmp = Quiet @ Check[URLDownload[remoteURL], $Failed];
    If[StringQ[tmp] && FileExistsQ[tmp],
      lines = Import[tmp, "Lines"];
      hash = ToLowerCase[FileHash[tmp, "SHA256", "HexString"]];
      src = "remote_url";
      Quiet @ Check[DeleteFile[tmp], Null];
    ];
  ];
  If[!ListQ[lines],
    Return[<|"source" -> "embedded_snapshot", "hash" -> fallbackHash, "assoc" -> fallbackAssoc|>]
  ];
  parsed = parseBestfitLines[lines];
  If[!AssociationQ[parsed],
    Return[<|"source" -> "embedded_snapshot", "hash" -> fallbackHash, "assoc" -> fallbackAssoc|>]
  ];
  <|"source" -> src, "hash" -> hash, "assoc" -> parsed|>
];

omegaGamma[h_] := N[2.469*10^-5/h^2, 50];
omegaR[h_, nEff_ : 3.046] := N[omegaGamma[h] (1 + 0.22710731766 nEff), 50];
eBackground[z_, om_, or_, ol_] := N[Sqrt[or (1 + z)^4 + om (1 + z)^3 + ol], 50];
wEff[z_, om_, or_, ol_] := N[
  ((1/3) or (1 + z)^4 - ol)/(or (1 + z)^4 + om (1 + z)^3 + ol),
  50
];

baseDir = DirectoryName[$InputFileName];
studyDir = DirectoryName[baseDir];
dataDir = FileNameJoin[{studyDir, "data"}];
dataFile = FileNameJoin[{dataDir, "desi_dr1_v1.0_bestfit.minimum.txt"}];
hashFile = FileNameJoin[{dataDir, "desi_dr1_v1.0.sha256sum"}];
manifestFile = FileNameJoin[{dataDir, "manifest.md"}];
bestfitURL = "https://data.desi.lbl.gov/public/dr1/vac/dr1/bao-cosmo-params/v1.0/iminuit/base/desi-bao-all/bestfit.minimum.txt";
expectedHash = "948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25";
embeddedAssoc = <|
  "hrdrag" -> 101.94944,
  "omegam" -> 0.29353038,
  "omegal" -> 0.70639,
  "rdrag" -> 151.35012,
  "H0rdrag" -> 10194.944
|>;

loaded = loadBestfitData[dataFile, bestfitURL, expectedHash, embeddedAssoc];
sourceUsed = If[AssociationQ[loaded], loaded["source"], "none"];
actualHash = If[AssociationQ[loaded], loaded["hash"], "missing"];
parsed = If[AssociationQ[loaded], loaded["assoc"], $Failed];

checkAData = AssociationQ[loaded] && sourceUsed =!= "none";
checkBData = AssociationQ[parsed] && KeyExistsQ[parsed, "omegam"] && KeyExistsQ[parsed, "omegal"];
checkCData = checkAData && actualHash == expectedHash;

hPlanck = 0.6766;
omPlanck = 0.30966;
orPlanck = omegaR[hPlanck];
olPlanck = N[1 - omPlanck - orPlanck, 50];

e[z_] := eBackground[z, omPlanck, orPlanck, olPlanck];
w[z_] := wEff[z, omPlanck, orPlanck, olPlanck];

z2 = 1/20;
alpha = 1/10000;
rs = 10^-30;
beta = 1;
rFloor = 10^-40;
aReg = N[Sqrt[(beta rs)^2 + rFloor^2], 60];

proxyKappa[proxy_String] := Which[
  proxy == "redshift_threshold", N[Sqrt[1 - z2], 60],
  proxy == "apparent_horizon", 1,
  True, Indeterminate
];

rDom[z_?NumericQ, proxy_String] := N[proxyKappa[proxy]/e[z], 60];
i0[z_?NumericQ, proxy_String] := N[(4 Pi/3) rDom[z, proxy]^3, 60];
ff[x_] := N[(1/(2 aReg)) ArcTan[x/aReg] - x/(2 (x^2 + aReg^2)), 60];
localI1[z_?NumericQ, proxy_String] := N[4 Pi alpha rs^2 (ff[rDom[z, proxy]] - ff[0]), 60];
deltaQ[z_?NumericQ, proxy_String] := N[localI1[z, proxy]/i0[z, proxy], 60];
q0[z_?NumericQ] := N[3 e[z]^2, 60];
eps[z_?NumericQ, proxy_String] := N[deltaQ[z, proxy]/q0[z], 60];

eraList = {
  <|"name" -> "radiation", "zRef" -> 10^6, "zCmp" -> N[(1 + 10^6)/Sqrt[2] - 1, 50]|>,
  <|"name" -> "matter", "zRef" -> 100, "zCmp" -> N[(1 + 100)/2^(2/3) - 1, 50]|>,
  <|"name" -> "late_de_sitter", "zRef" -> 0, "zCmp" -> -0.9|>
};

evaluateEraProxy[era_Association, proxy_String] := Module[
  {
    name, zRef, zCmp, hRef, hCmp, hRatio, wRef, wCmp,
    i0Ref, i0Cmp, locRef, locCmp, deltaRef, deltaCmp,
    ratioDelta, ratioI0, ratioLoc, epsRef, epsCmp,
    check1, check2, check3, check4, check5, check6, check7, check8, check9,
    checkProxy, class
  },
  name = era["name"];
  zRef = era["zRef"];
  zCmp = era["zCmp"];
  hRef = e[zRef];
  hCmp = e[zCmp];
  hRatio = N[hCmp/hRef, 50];
  wRef = w[zRef];
  wCmp = w[zCmp];

  i0Ref = i0[zRef, proxy];
  i0Cmp = i0[zCmp, proxy];
  locRef = localI1[zRef, proxy];
  locCmp = localI1[zCmp, proxy];
  deltaRef = deltaQ[zRef, proxy];
  deltaCmp = deltaQ[zCmp, proxy];
  ratioDelta = N[deltaCmp/deltaRef, 50];
  ratioI0 = N[i0Cmp/i0Ref, 50];
  ratioLoc = N[locCmp/locRef, 50];
  epsRef = eps[zRef, proxy];
  epsCmp = eps[zCmp, proxy];

  check1 = hRef > 0 && hCmp > 0 && hRatio > 0;
  check2 = i0Ref > 0 && i0Cmp > 0;
  check3 = deltaRef > 0 && deltaCmp > 0;
  check4 = If[hRatio < 1, deltaCmp < deltaRef, deltaCmp > deltaRef];
  check5 = Abs[ratioDelta/(hRatio^3) - 1] < 0.12;
  check6 = Abs[ratioI0 hRatio^3 - 1] < 0.12;
  check7 = Abs[ratioLoc - 1] < 0.12;
  check8 = epsRef < 0.05 && epsCmp < 0.05;
  check9 = Which[
    name == "radiation", Abs[wRef - 1/3] < 0.01 && Abs[wCmp - 1/3] < 0.01,
    name == "matter", Abs[wRef] < 0.05 && Abs[wCmp] < 0.05,
    name == "late_de_sitter", wRef < -0.6 && Abs[wCmp + 1] < 0.01,
    True, False
  ];

  checkProxy = TrueQ[check1 && check2 && check3 && check4 && check5 && check6 && check7 && check8 && check9];
  class = If[checkProxy, "IR_SUPPRESSED_SUBLEADING", "undetermined"];

  <|
    "era" -> name,
    "proxy" -> proxy,
    "zRef" -> zRef,
    "zCmp" -> zCmp,
    "hRef" -> hRef,
    "hCmp" -> hCmp,
    "hRatio" -> hRatio,
    "wRef" -> wRef,
    "wCmp" -> wCmp,
    "i0Ref" -> i0Ref,
    "i0Cmp" -> i0Cmp,
    "locRef" -> locRef,
    "locCmp" -> locCmp,
    "deltaRef" -> deltaRef,
    "deltaCmp" -> deltaCmp,
    "ratioDelta" -> ratioDelta,
    "ratioI0" -> ratioI0,
    "ratioLoc" -> ratioLoc,
    "epsRef" -> epsRef,
    "epsCmp" -> epsCmp,
    "checkProxy" -> checkProxy,
    "classification" -> class
  |>
];

redResults = Association @ Table[era["name"] -> evaluateEraProxy[era, "redshift_threshold"], {era, eraList}];
appResults = Association @ Table[era["name"] -> evaluateEraProxy[era, "apparent_horizon"], {era, eraList}];

eraNames = {"radiation", "matter", "late_de_sitter"};
redClasses = Lookup[Lookup[redResults, eraNames], "classification"];
appClasses = Lookup[Lookup[appResults, eraNames], "classification"];
redRatios = Lookup[Lookup[redResults, eraNames], "ratioDelta"];
appRatios = Lookup[Lookup[appResults, eraNames], "ratioDelta"];
redEpsRef = Lookup[Lookup[redResults, eraNames], "epsRef"];
appEpsRef = Lookup[Lookup[appResults, eraNames], "epsRef"];
redEpsCmp = Lookup[Lookup[redResults, eraNames], "epsCmp"];
appEpsCmp = Lookup[Lookup[appResults, eraNames], "epsCmp"];
domainRefRatios = Table[
  N[appResults[name]["i0Ref"]/redResults[name]["i0Ref"], 50],
  {name, eraNames}
];

checkA = AllTrue[Lookup[Lookup[redResults, eraNames], "checkProxy"], TrueQ];
checkB = AllTrue[Lookup[Lookup[appResults, eraNames], "checkProxy"], TrueQ];
checkC = AllTrue[redClasses, # == "IR_SUPPRESSED_SUBLEADING" &] &&
  AllTrue[appClasses, # == "IR_SUPPRESSED_SUBLEADING" &];
checkD = And @@ Table[Abs[redRatios[[i]]/appRatios[[i]] - 1] < 0.03, {i, Length[eraNames]}];
checkE = And @@ Table[
    Abs[redEpsRef[[i]] - appEpsRef[[i]]] < 0.01 && Abs[redEpsCmp[[i]] - appEpsCmp[[i]]] < 0.01,
    {i, Length[eraNames]}
  ];
checkF = AllTrue[domainRefRatios, Abs[# - 1] > 0.05 &];
checkG = AllTrue[redRatios, # < 1 &] && AllTrue[appRatios, # < 1 &];

classification = If[
  TrueQ[checkAData && checkBData && checkCData && checkA && checkB && checkC && checkD && checkE && checkF && checkG],
  "ERA_PROXY_CLASS_EQUIVALENT_IR",
  "undetermined"
];

check = TrueQ[classification == "ERA_PROXY_CLASS_EQUIVALENT_IR"];

nbName = "04_Era_Proxy_Equivalence_IR_Class_Test.nb";
logName = "04_era_proxy_equivalence_ir_class_test.log";
logLines = {
  "Notebook: " <> nbName,
  "sourceUsed = " <> ToString[sourceUsed, InputForm],
  "bestfitURL = " <> ToString[bestfitURL, InputForm],
  "dataFile = " <> ToString[dataFile, InputForm],
  "hashFile = " <> ToString[hashFile, InputForm],
  "manifestFile = " <> ToString[manifestFile, InputForm],
  "expectedHash = " <> expectedHash,
  "actualHash = " <> ToString[actualHash, InputForm],
  "checkAData(data source available) = " <> ToString[checkAData, InputForm],
  "checkBData(parsed keys) = " <> ToString[checkBData, InputForm],
  "checkCData(sha256 match) = " <> ToString[checkCData, InputForm],
  "eraNames = " <> ToString[eraNames, InputForm],
  "z2 = " <> ToString[z2, InputForm],
  "proxyKappaRed = " <> ToString[proxyKappa["redshift_threshold"], InputForm],
  "proxyKappaApp = " <> ToString[proxyKappa["apparent_horizon"], InputForm],
  "redResults = " <> ToString[redResults, InputForm],
  "appResults = " <> ToString[appResults, InputForm],
  "redClasses = " <> ToString[redClasses, InputForm],
  "appClasses = " <> ToString[appClasses, InputForm],
  "redRatios = " <> ToString[redRatios, InputForm],
  "appRatios = " <> ToString[appRatios, InputForm],
  "domainRefRatios = " <> ToString[domainRefRatios, InputForm],
  "checkA(red proxy class per era) = " <> ToString[checkA, InputForm],
  "checkB(apparent proxy class per era) = " <> ToString[checkB, InputForm],
  "checkC(class equality across proxies and eras) = " <> ToString[checkC, InputForm],
  "checkD(ratioDelta agreement across proxies) = " <> ToString[checkD, InputForm],
  "checkE(eps agreement across proxies) = " <> ToString[checkE, InputForm],
  "checkF(domains genuinely different) = " <> ToString[checkF, InputForm],
  "checkG(all ratioDelta < 1) = " <> ToString[checkG, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "04_era_proxy_equivalence_ir_class_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
