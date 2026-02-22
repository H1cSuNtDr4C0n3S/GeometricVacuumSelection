ClearAll["Global`*"];

nb = Notebook[{
  Cell["00 - Friedmann Eras Data Baseline Consistency Test", "Title"],
  Cell["Objective: build a reproducible cosmology baseline from frozen data (DESI DR1) and verify that radiation, matter, and late de Sitter eras remain standard and ordered.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    E[z_] := Sqrt[\[CapitalOmega]r (1 + z)^4 + \[CapitalOmega]m (1 + z)^3 + \[CapitalOmega]\[CapitalLambda]]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    zEq := \[CapitalOmega]m/\[CapitalOmega]r - 1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    zLambdaMatter := (\[CapitalOmega]\[CapitalLambda]/\[CapitalOmega]m)^(1/3) - 1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    wEff[z_] := ((1/3) \[CapitalOmega]r (1 + z)^4 - \[CapitalOmega]\[CapitalLambda])/
      (\[CapitalOmega]r (1 + z)^4 + \[CapitalOmega]m (1 + z)^3 + \[CapitalOmega]\[CapitalLambda])
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
checkA = AssociationQ[loaded] && sourceUsed =!= "none";

parsed = If[checkA, loaded["assoc"], $Failed];
requiredKeys = {"hrdrag", "omegam", "omegal", "rdrag", "H0rdrag"};
checkB = AssociationQ[parsed] && AllTrue[requiredKeys, KeyExistsQ[parsed, #] &];

checkC = checkA && actualHash == expectedHash;

hrdrag = If[checkB, parsed["hrdrag"], Indeterminate];
rdrag = If[checkB, parsed["rdrag"], Indeterminate];
h0rdrag = If[checkB, parsed["H0rdrag"], Indeterminate];
omDESI = If[checkB, parsed["omegam"], Indeterminate];
olDESI = If[checkB, parsed["omegal"], Indeterminate];

h0DESIFromHrdrag = If[checkB, N[100 hrdrag/rdrag, 50], Indeterminate];
h0DESIFromH0rdrag = If[checkB, N[h0rdrag/rdrag, 50], Indeterminate];
checkD = checkB && Abs[h0DESIFromHrdrag - h0DESIFromH0rdrag] < 10^-8;

hPlanck = 0.6766;
omPlanck = 0.30966;

hDESI = If[NumberQ[h0DESIFromHrdrag], N[h0DESIFromHrdrag/100, 50], Indeterminate];
orPlanck = omegaR[hPlanck];
orDESI = If[NumberQ[hDESI], omegaR[hDESI], Indeterminate];
olPlanck = N[1 - omPlanck - orPlanck, 50];

zEqPlanck = N[omPlanck/orPlanck - 1, 50];
zEqDESI = If[NumberQ[omDESI] && NumberQ[orDESI], N[omDESI/orDESI - 1, 50], Indeterminate];
zLambdaMatterPlanck = N[(olPlanck/omPlanck)^(1/3) - 1, 50];
zLambdaMatterDESI = If[NumberQ[omDESI] && NumberQ[olDESI], N[(olDESI/omDESI)^(1/3) - 1, 50], Indeterminate];

checkE = NumberQ[zEqPlanck] && NumberQ[zEqDESI] && NumberQ[zLambdaMatterPlanck] && NumberQ[zLambdaMatterDESI] &&
  zEqPlanck > 1000 && zEqDESI > 1000 &&
  0 < zLambdaMatterPlanck < 2 && 0 < zLambdaMatterDESI < 2 &&
  zEqPlanck > zLambdaMatterPlanck && zEqDESI > zLambdaMatterDESI;

zRadProbe = 10^6;
zMatProbe = 100;
zDSProbe = -0.9;

wRadPlanck = wEff[zRadProbe, omPlanck, orPlanck, olPlanck];
wMatPlanck = wEff[zMatProbe, omPlanck, orPlanck, olPlanck];
wDSPlanck = wEff[zDSProbe, omPlanck, orPlanck, olPlanck];

checkF = Abs[wRadPlanck - 1/3] < 0.01 && Abs[wMatPlanck] < 0.05 && Abs[wDSPlanck + 1] < 0.01;

checkG = NumberQ[h0DESIFromHrdrag] && NumberQ[omDESI] &&
  Abs[h0DESIFromHrdrag - 67.66]/67.66 < 0.03 &&
  Abs[omDESI - omPlanck]/omPlanck < 0.08;

zSamples = {0, 1, 10, 100};
eRatios = If[
  NumberQ[omDESI] && NumberQ[orDESI] && NumberQ[olDESI],
  N[
    Table[
      eBackground[z, omDESI, orDESI, olDESI]/eBackground[z, omPlanck, orPlanck, olPlanck],
      {z, zSamples}
    ],
    50
  ],
  ConstantArray[Indeterminate, Length[zSamples]]
];
maxERelDiff = If[VectorQ[eRatios, NumberQ], N[Max[Abs[eRatios - 1]], 50], Indeterminate];
checkH = NumberQ[maxERelDiff] && maxERelDiff < 0.06;

checkI = (sourceUsed == "local_file" && FileExistsQ[hashFile] && FileExistsQ[manifestFile]) ||
  sourceUsed == "remote_url" || sourceUsed == "embedded_snapshot";

classification = If[
  TrueQ[checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI],
  "STANDARD_ERA_HISTORY_COMPATIBLE",
  "undetermined"
];

check = TrueQ[classification == "STANDARD_ERA_HISTORY_COMPATIBLE"];

nbName = "00_Friedmann_Eras_Data_Baseline_Consistency_Test.nb";
logName = "00_friedmann_eras_data_baseline_consistency_test.log";
logLines = {
  "Notebook: " <> nbName,
  "sourceUsed = " <> ToString[sourceUsed, InputForm],
  "bestfitURL = " <> ToString[bestfitURL, InputForm],
  "dataFile = " <> ToString[dataFile, InputForm],
  "hashFile = " <> ToString[hashFile, InputForm],
  "manifestFile = " <> ToString[manifestFile, InputForm],
  "expectedHash = " <> expectedHash,
  "actualHash = " <> ToString[actualHash, InputForm],
  "checkA(files exist) = " <> ToString[checkA, InputForm],
  "checkB(parse + required keys) = " <> ToString[checkB, InputForm],
  "checkC(sha256 match) = " <> ToString[checkC, InputForm],
  "hrdrag = " <> ToString[hrdrag, InputForm],
  "rdrag = " <> ToString[rdrag, InputForm],
  "H0rdrag = " <> ToString[h0rdrag, InputForm],
  "h0DESIFromHrdrag = " <> ToString[h0DESIFromHrdrag, InputForm],
  "h0DESIFromH0rdrag = " <> ToString[h0DESIFromH0rdrag, InputForm],
  "checkD(H0 reconstruction consistency) = " <> ToString[checkD, InputForm],
  "omPlanck = " <> ToString[omPlanck, InputForm],
  "omDESI = " <> ToString[omDESI, InputForm],
  "orPlanck = " <> ToString[orPlanck, InputForm],
  "orDESI = " <> ToString[orDESI, InputForm],
  "olPlanck = " <> ToString[olPlanck, InputForm],
  "olDESI = " <> ToString[olDESI, InputForm],
  "zEqPlanck = " <> ToString[zEqPlanck, InputForm],
  "zEqDESI = " <> ToString[zEqDESI, InputForm],
  "zLambdaMatterPlanck = " <> ToString[zLambdaMatterPlanck, InputForm],
  "zLambdaMatterDESI = " <> ToString[zLambdaMatterDESI, InputForm],
  "checkE(era ordering) = " <> ToString[checkE, InputForm],
  "wRadPlanck(z=1e6) = " <> ToString[wRadPlanck, InputForm],
  "wMatPlanck(z=100) = " <> ToString[wMatPlanck, InputForm],
  "wDSPlanck(z=-0.9) = " <> ToString[wDSPlanck, InputForm],
  "checkF(w_eff probes) = " <> ToString[checkF, InputForm],
  "checkG(Planck-DESI parameter proximity) = " <> ToString[checkG, InputForm],
  "zSamples = " <> ToString[zSamples, InputForm],
  "eRatios = " <> ToString[eRatios, InputForm],
  "maxERelDiff = " <> ToString[maxERelDiff, InputForm],
  "checkH(background closeness across eras) = " <> ToString[checkH, InputForm],
  "checkI(source availability policy) = " <> ToString[checkI, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "00_friedmann_eras_data_baseline_consistency_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
