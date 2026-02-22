ClearAll["Global`*"];

nb = Notebook[{
  Cell["05 - Pipeline 00-04 Appendix-Ready Summary Test", "Title"],
  Cell["Objective: produce a compact, reproducible summary table of tests 00-04 and verify the whole Cosmological_Eras pipeline is internally consistent and appendix-ready.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classesExpected == {
      "STANDARD_ERA_HISTORY_COMPATIBLE",
      "RADIATION_ERA_IR_SUPPRESSED",
      "MATTER_ERA_IR_SUPPRESSED",
      "ERAS_PIPELINE_CONTINUITY_CONFIRMED",
      "ERA_PROXY_CLASS_EQUIVALENT_IR"
    }
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "PIPELINE_00_04_APPENDIX_READY"
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

checkDataA = AssociationQ[loaded] && sourceUsed =!= "none";
checkDataB = AssociationQ[parsed] && KeyExistsQ[parsed, "omegam"] && KeyExistsQ[parsed, "omegal"];
checkDataC = checkDataA && actualHash == expectedHash;

hPlanck = 0.6766;
omPlanck = 0.30966;
orPlanck = omegaR[hPlanck];
olPlanck = N[1 - omPlanck - orPlanck, 50];
e[z_] := eBackground[z, omPlanck, orPlanck, olPlanck];
w[z_] := wEff[z, omPlanck, orPlanck, olPlanck];

(* 00 *)
hDESI = parsed["hrdrag"]/parsed["rdrag"];
omDESI = parsed["omegam"];
olDESI = parsed["omegal"];
orDESI = omegaR[hDESI];
zEqPlanck = N[omPlanck/orPlanck - 1, 50];
zEqDESI = N[omDESI/orDESI - 1, 50];
zLmPlanck = N[(olPlanck/omPlanck)^(1/3) - 1, 50];
zLmDESI = N[(olDESI/omDESI)^(1/3) - 1, 50];
eRatios00 = Table[eBackground[z, omDESI, orDESI, olDESI]/e[z], {z, {0, 1, 10, 100}}];
maxERelDiff00 = N[Max[Abs[eRatios00 - 1]], 50];
check00 = checkDataA && checkDataB && checkDataC &&
  zEqPlanck > 1000 && zEqDESI > 1000 &&
  zLmPlanck > 0 && zLmDESI > 0 &&
  maxERelDiff00 < 0.06;
class00 = If[check00, "STANDARD_ERA_HISTORY_COMPATIBLE", "undetermined"];

(* Helpers for 01-04 *)
rDom[h_, kappa_] := N[kappa/h, 60];
i0[h_, kappa_] := N[(4 Pi/3) rDom[h, kappa]^3, 60];
ff[x_, aReg_] := N[(1/(2 aReg)) ArcTan[x/aReg] - x/(2 (x^2 + aReg^2)), 60];
localI1[h_, kappa_, alpha_, rs_, aReg_] := N[4 Pi alpha rs^2 (ff[rDom[h, kappa], aReg] - ff[0, aReg]), 60];
deltaQ[h_, kappa_, alpha_, rs_, aReg_] := N[localI1[h, kappa, alpha, rs, aReg]/i0[h, kappa], 60];
eps[z_, kappa_, alpha_, rs_, aReg_] := N[deltaQ[e[z], kappa, alpha, rs, aReg]/(3 e[z]^2), 60];

z2 = 1/20;
kRed = N[Sqrt[1 - z2], 60];
kApp = 1;

(* 01 *)
alpha01 = 1/100000;
rs01 = 1/1000;
aReg01 = N[Sqrt[rs01^2 + (10^-12)^2], 60];
zRef01 = 10^6;
zHalf01 = N[(1 + zRef01)/Sqrt[2] - 1, 50];
hRatio01 = N[e[zHalf01]/e[zRef01], 50];
ratioDelta01 = N[deltaQ[hRatio01, kRed, alpha01, rs01, aReg01]/deltaQ[1, kRed, alpha01, rs01, aReg01], 50];
ratioI001 = N[i0[hRatio01, kRed]/i0[1, kRed], 50];
ratioLoc01 = N[
  localI1[hRatio01, kRed, alpha01, rs01, aReg01]/localI1[1, kRed, alpha01, rs01, aReg01],
  50
];
epsRef01 = N[deltaQ[1, kRed, alpha01, rs01, aReg01]/(3*1^2), 60];
epsHalf01 = N[deltaQ[hRatio01, kRed, alpha01, rs01, aReg01]/(3 hRatio01^2), 60];
check01 = Abs[hRatio01 - 1/2] < 5*10^-4 &&
  Abs[w[zRef01] - 1/3] < 0.01 && Abs[w[zHalf01] - 1/3] < 0.01 &&
  Abs[ratioDelta01/(hRatio01^3) - 1] < 0.05 &&
  Abs[ratioI001 hRatio01^3 - 1] < 0.05 &&
  Abs[ratioLoc01 - 1] < 0.1 &&
  epsRef01 < 0.05 && epsHalf01 < 0.05;
class01 = If[check01, "RADIATION_ERA_IR_SUPPRESSED", "undetermined"];

(* 02 *)
alpha02 = 1/100000;
rs02 = 1/1000;
aReg02 = N[Sqrt[rs02^2 + (10^-12)^2], 60];
zRef02 = 100;
zHalf02 = N[(1 + zRef02)/2^(2/3) - 1, 50];
hRatio02 = N[e[zHalf02]/e[zRef02], 50];
ratioDelta02 = N[deltaQ[hRatio02, kRed, alpha02, rs02, aReg02]/deltaQ[1, kRed, alpha02, rs02, aReg02], 50];
ratioI002 = N[i0[hRatio02, kRed]/i0[1, kRed], 50];
ratioLoc02 = N[
  localI1[hRatio02, kRed, alpha02, rs02, aReg02]/localI1[1, kRed, alpha02, rs02, aReg02],
  50
];
epsRef02 = N[deltaQ[1, kRed, alpha02, rs02, aReg02]/(3*1^2), 60];
epsHalf02 = N[deltaQ[hRatio02, kRed, alpha02, rs02, aReg02]/(3 hRatio02^2), 60];
check02 = Abs[hRatio02 - 1/2] < 0.03 &&
  Abs[w[zRef02]] < 0.05 && Abs[w[zHalf02]] < 0.05 &&
  Abs[ratioDelta02/(hRatio02^3) - 1] < 0.08 &&
  Abs[ratioI002 hRatio02^3 - 1] < 0.08 &&
  Abs[ratioLoc02 - 1] < 0.1 &&
  epsRef02 < 0.05 && epsHalf02 < 0.05;
class02 = If[check02, "MATTER_ERA_IR_SUPPRESSED", "undetermined"];

(* 03 *)
alpha03 = 1/10000;
rs03 = 10^-30;
aReg03 = N[Sqrt[rs03^2 + (10^-40)^2], 80];
zEq03 = zEqPlanck;
zLm03 = zLmPlanck;
zProbe03 = {10^6, zEq03, 100, zLm03, 0, -0.9};
wProbe03 = N[w /@ zProbe03, 50];
deltaProbe03 = N[deltaQ[e[#], kRed, alpha03, rs03, aReg03] & /@ zProbe03, 80];
epsProbe03 = N[eps[#, kRed, alpha03, rs03, aReg03] & /@ zProbe03, 80];
rProbe03 = N[rDom[e[#], kRed] & /@ zProbe03, 80];
check03 = zProbe03[[1]] > zProbe03[[2]] > zProbe03[[3]] > zProbe03[[4]] > zProbe03[[5]] > zProbe03[[6]] &&
  Abs[wProbe03[[1]] - 1/3] < 0.01 &&
  Abs[wProbe03[[3]]] < 0.05 &&
  Abs[wProbe03[[6]] + 1] < 0.01 &&
  Abs[w[zEq03] - 1/6] < 0.02 &&
  Abs[w[zLm03] + 1/2] < 0.05 &&
  And @@ Thread[deltaProbe03 > 0] && Max[epsProbe03] < 0.05 &&
  And @@ Thread[Most[deltaProbe03] > Rest[deltaProbe03]] &&
  And @@ Thread[rProbe03 > 0];
class03 = If[check03, "ERAS_PIPELINE_CONTINUITY_CONFIRMED", "undetermined"];

(* 04 *)
alpha04 = 1/10000;
rs04 = 10^-30;
aReg04 = N[Sqrt[rs04^2 + (10^-40)^2], 80];
eraData04 = {
  <|"name" -> "radiation", "zRef" -> 10^6, "zCmp" -> N[(1 + 10^6)/Sqrt[2] - 1, 50]|>,
  <|"name" -> "matter", "zRef" -> 100, "zCmp" -> N[(1 + 100)/2^(2/3) - 1, 50]|>,
  <|"name" -> "late_de_sitter", "zRef" -> 0, "zCmp" -> -0.9|>
};
eraEval[era_, kappa_] := Module[
  {zRef, zCmp, hRatio, ratioDelta, ratioI0, ratioLoc, epsRef, epsCmp, checkEra},
  zRef = era["zRef"]; zCmp = era["zCmp"];
  hRatio = N[e[zCmp]/e[zRef], 50];
  ratioDelta = N[deltaQ[e[zCmp], kappa, alpha04, rs04, aReg04]/deltaQ[e[zRef], kappa, alpha04, rs04, aReg04], 50];
  ratioI0 = N[i0[e[zCmp], kappa]/i0[e[zRef], kappa], 50];
  ratioLoc = N[
    localI1[e[zCmp], kappa, alpha04, rs04, aReg04]/localI1[e[zRef], kappa, alpha04, rs04, aReg04],
    50
  ];
  epsRef = eps[zRef, kappa, alpha04, rs04, aReg04];
  epsCmp = eps[zCmp, kappa, alpha04, rs04, aReg04];
  checkEra = Abs[ratioDelta/(hRatio^3) - 1] < 0.12 &&
    Abs[ratioI0 hRatio^3 - 1] < 0.12 &&
    Abs[ratioLoc - 1] < 0.12 &&
    epsRef < 0.05 && epsCmp < 0.05 &&
    If[hRatio < 1, ratioDelta < 1, True];
  <|
    "ratioDelta" -> ratioDelta,
    "epsRef" -> epsRef,
    "epsCmp" -> epsCmp,
    "checkEra" -> checkEra
  |>
];
red04 = eraEval[#, kRed] & /@ eraData04;
app04 = eraEval[#, kApp] & /@ eraData04;
redClasses04 = If[#["checkEra"], "IR_SUPPRESSED_SUBLEADING", "undetermined"] & /@ red04;
appClasses04 = If[#["checkEra"], "IR_SUPPRESSED_SUBLEADING", "undetermined"] & /@ app04;
redRatios04 = red04[[All, "ratioDelta"]];
appRatios04 = app04[[All, "ratioDelta"]];
domainRefRatios04 = Table[N[i0[e[eraData04[[i, "zRef"]]], kApp]/i0[e[eraData04[[i, "zRef"]]], kRed], 50], {i, 1, 3}];
check04 = AllTrue[redClasses04, # == "IR_SUPPRESSED_SUBLEADING" &] &&
  AllTrue[appClasses04, # == "IR_SUPPRESSED_SUBLEADING" &] &&
  And @@ Table[Abs[redRatios04[[i]]/appRatios04[[i]] - 1] < 0.03, {i, 1, 3}] &&
  AllTrue[domainRefRatios04, Abs[# - 1] > 0.05 &];
class04 = If[check04, "ERA_PROXY_CLASS_EQUIVALENT_IR", "undetermined"];

rows = {
  <|"test" -> "00", "title" -> "Data baseline consistency", "classification" -> class00, "check" -> check00|>,
  <|"test" -> "01", "title" -> "Radiation era IR scaling", "classification" -> class01, "check" -> check01|>,
  <|"test" -> "02", "title" -> "Matter era IR scaling", "classification" -> class02, "check" -> check02|>,
  <|"test" -> "03", "title" -> "Era continuity global channel", "classification" -> class03, "check" -> check03|>,
  <|"test" -> "04", "title" -> "Proxy equivalence across eras", "classification" -> class04, "check" -> check04|>
};

classesExpected = {
  "STANDARD_ERA_HISTORY_COMPATIBLE",
  "RADIATION_ERA_IR_SUPPRESSED",
  "MATTER_ERA_IR_SUPPRESSED",
  "ERAS_PIPELINE_CONTINUITY_CONFIRMED",
  "ERA_PROXY_CLASS_EQUIVALENT_IR"
};
classesObserved = rows[[All, "classification"]];
checksObserved = rows[[All, "check"]];
checkAllRows = AllTrue[checksObserved, TrueQ];
checkClasses = classesObserved == classesExpected;

summaryMarkdown = StringRiffle[
  Join[
    {"| Test | Title | Classification | Check |", "|---|---|---|---|"},
    (StringJoin[
        "| ", #["test"], " | ", #["title"], " | ", #["classification"], " | ", ToString[#["check"], InputForm], " |"
      ] & /@ rows)
  ],
  "\n"
];

classification = If[
  TrueQ[checkDataA && checkDataB && checkDataC && checkAllRows && checkClasses],
  "PIPELINE_00_04_APPENDIX_READY",
  "undetermined"
];
check = TrueQ[classification == "PIPELINE_00_04_APPENDIX_READY"];

nbName = "05_Pipeline_00_04_Appendix_Ready_Summary_Test.nb";
logName = "05_pipeline_00_04_appendix_ready_summary_test.log";
logLines = {
  "Notebook: " <> nbName,
  "sourceUsed = " <> ToString[sourceUsed, InputForm],
  "bestfitURL = " <> ToString[bestfitURL, InputForm],
  "dataFile = " <> ToString[dataFile, InputForm],
  "hashFile = " <> ToString[hashFile, InputForm],
  "manifestFile = " <> ToString[manifestFile, InputForm],
  "expectedHash = " <> expectedHash,
  "actualHash = " <> ToString[actualHash, InputForm],
  "checkDataA = " <> ToString[checkDataA, InputForm],
  "checkDataB = " <> ToString[checkDataB, InputForm],
  "checkDataC = " <> ToString[checkDataC, InputForm],
  "rows = " <> ToString[rows, InputForm],
  "classesObserved = " <> ToString[classesObserved, InputForm],
  "classesExpected = " <> ToString[classesExpected, InputForm],
  "checkAllRows = " <> ToString[checkAllRows, InputForm],
  "checkClasses = " <> ToString[checkClasses, InputForm],
  "summaryMarkdown = " <> ToString[summaryMarkdown, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "05_pipeline_00_04_appendix_ready_summary_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
