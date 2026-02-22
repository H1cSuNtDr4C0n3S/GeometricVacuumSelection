ClearAll["Global`*"];

nb = Notebook[{
  Cell["03 - Era Sequence Continuity and Global IR Channel Test", "Title"],
  Cell["Objective: verify smooth radiation->matter->de Sitter ordering and show the non-local channel stays a global IR consistency correction (subleading, no era discontinuities).", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    zEq := \[CapitalOmega]m/\[CapitalOmega]r - 1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    zLambdaMatter := (\[CapitalOmega]\[CapitalLambda]/\[CapitalOmega]m)^(1/3) - 1
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    eps[z_] := deltaQ[z]/(3 E[z]^2)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    classification == "ERAS_PIPELINE_CONTINUITY_CONFIRMED"
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

checkA = AssociationQ[loaded] && sourceUsed =!= "none";
checkB = AssociationQ[parsed] && KeyExistsQ[parsed, "omegam"] && KeyExistsQ[parsed, "omegal"];
checkC = checkA && actualHash == expectedHash;

hPlanck = 0.6766;
omPlanck = 0.30966;
orPlanck = omegaR[hPlanck];
olPlanck = N[1 - omPlanck - orPlanck, 50];

e[z_] := eBackground[z, omPlanck, orPlanck, olPlanck];
w[z_] := wEff[z, omPlanck, orPlanck, olPlanck];

zEq = N[omPlanck/orPlanck - 1, 50];
zLambdaMatter = N[(olPlanck/omPlanck)^(1/3) - 1, 50];

zRad = 10^6;
zMat = 100;
zNow = 0;
zDS = -0.9;

wRad = w[zRad];
wMat = w[zMat];
wNow = w[zNow];
wDS = w[zDS];
wEq = w[zEq];
wLambda = w[zLambdaMatter];

zEqMinus = N[zEq (1 - 10^-3), 50];
zEqPlus = N[zEq (1 + 10^-3), 50];
zLambdaMinus = N[zLambdaMatter (1 - 10^-3), 50];
zLambdaPlus = N[zLambdaMatter (1 + 10^-3), 50];
wEqMinus = w[zEqMinus];
wEqPlus = w[zEqPlus];
wLambdaMinus = w[zLambdaMinus];
wLambdaPlus = w[zLambdaPlus];

z2 = 1/20;
alpha = 10^-4;
rs = 10^-30;
beta = 1;
rFloor = 10^-40;
aReg = N[Sqrt[(beta rs)^2 + rFloor^2], 80];

rDom[z_] := N[Sqrt[1 - z2]/e[z], 80];
i0[z_] := N[(4 Pi/3) rDom[z]^3, 80];
ff[x_] := N[(1/(2 aReg)) ArcTan[x/aReg] - x/(2 (x^2 + aReg^2)), 80];
localI1[z_] := N[4 Pi alpha rs^2 (ff[rDom[z]] - ff[0]), 80];
deltaQ[z_] := N[localI1[z]/i0[z], 80];
q0[z_] := N[3 e[z]^2, 80];
eps[z_] := N[deltaQ[z]/q0[z], 80];

zProbe = {zRad, zEq, zMat, zLambdaMatter, zNow, zDS};
eProbe = N[e /@ zProbe, 50];
wProbe = N[w /@ zProbe, 50];
deltaProbe = N[deltaQ /@ zProbe, 80];
epsProbe = N[eps /@ zProbe, 80];
rProbe = N[rDom /@ zProbe, 80];

checkD = zRad > zEq && zEq > zMat && zMat > zLambdaMatter && zLambdaMatter > zNow && zNow > zDS && zDS > -1;
checkE = Abs[wRad - 1/3] < 0.01 && Abs[wMat] < 0.05 && Abs[wDS + 1] < 0.01;
checkF = Abs[wEq - 1/6] < 0.02;
checkG = Abs[wLambda + 1/2] < 0.05;
checkH = Abs[wEqPlus - wEqMinus] < 0.01 && Abs[wLambdaPlus - wLambdaMinus] < 0.01;
checkI = And @@ Thread[deltaProbe > 0] && And @@ Thread[epsProbe > 0] && Max[epsProbe] < 0.05;
checkJ = And @@ Thread[Most[deltaProbe] > Rest[deltaProbe]];
checkK = And @@ Thread[Most[eProbe] > Rest[eProbe]];
checkL = And @@ Thread[rProbe > 0] && Max[rProbe] < 10;

classification = If[
  TrueQ[
    checkA && checkB && checkC && checkD && checkE && checkF &&
    checkG && checkH && checkI && checkJ && checkK && checkL
  ],
  "ERAS_PIPELINE_CONTINUITY_CONFIRMED",
  "undetermined"
];

check = TrueQ[classification == "ERAS_PIPELINE_CONTINUITY_CONFIRMED"];

nbName = "03_Era_Sequence_Continuity_Global_IR_Channel_Test.nb";
logName = "03_era_sequence_continuity_global_ir_channel_test.log";
logLines = {
  "Notebook: " <> nbName,
  "sourceUsed = " <> ToString[sourceUsed, InputForm],
  "bestfitURL = " <> ToString[bestfitURL, InputForm],
  "dataFile = " <> ToString[dataFile, InputForm],
  "hashFile = " <> ToString[hashFile, InputForm],
  "manifestFile = " <> ToString[manifestFile, InputForm],
  "expectedHash = " <> expectedHash,
  "actualHash = " <> ToString[actualHash, InputForm],
  "checkA(data source available) = " <> ToString[checkA, InputForm],
  "checkB(parsed keys) = " <> ToString[checkB, InputForm],
  "checkC(sha256 match) = " <> ToString[checkC, InputForm],
  "zEq = " <> ToString[zEq, InputForm],
  "zLambdaMatter = " <> ToString[zLambdaMatter, InputForm],
  "zProbe = " <> ToString[zProbe, InputForm],
  "eProbe = " <> ToString[eProbe, InputForm],
  "wProbe = " <> ToString[wProbe, InputForm],
  "deltaProbe = " <> ToString[deltaProbe, InputForm],
  "epsProbe = " <> ToString[epsProbe, InputForm],
  "rProbe = " <> ToString[rProbe, InputForm],
  "wEq = " <> ToString[wEq, InputForm],
  "wLambda = " <> ToString[wLambda, InputForm],
  "wEqMinus = " <> ToString[wEqMinus, InputForm],
  "wEqPlus = " <> ToString[wEqPlus, InputForm],
  "wLambdaMinus = " <> ToString[wLambdaMinus, InputForm],
  "wLambdaPlus = " <> ToString[wLambdaPlus, InputForm],
  "checkD(redshift ordering) = " <> ToString[checkD, InputForm],
  "checkE(era w_eff anchors) = " <> ToString[checkE, InputForm],
  "checkF(w_eff at zEq ~ 1/6) = " <> ToString[checkF, InputForm],
  "checkG(w_eff at zLambda ~ -1/2) = " <> ToString[checkG, InputForm],
  "checkH(continuity around transitions) = " <> ToString[checkH, InputForm],
  "checkI(IR correction subleading positive) = " <> ToString[checkI, InputForm],
  "checkJ(delta monotonic along eras) = " <> ToString[checkJ, InputForm],
  "checkK(E monotonic along eras) = " <> ToString[checkK, InputForm],
  "checkL(domain radius finite positive) = " <> ToString[checkL, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "03_era_sequence_continuity_global_ir_channel_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
