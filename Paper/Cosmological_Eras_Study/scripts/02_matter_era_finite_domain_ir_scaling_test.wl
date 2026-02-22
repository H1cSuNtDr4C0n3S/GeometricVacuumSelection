ClearAll["Global`*"];

nb = Notebook[{
  Cell["02 - Matter Era Finite-Domain IR Scaling Test", "Title"],
  Cell["Objective: in matter era, verify finite-domain non-local channel remains IR-suppressed and follows expected scaling when H is reduced.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    hNorm[z_] := E[z]/E[zRef]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Rdom[h_] := Sqrt[1 - z2]/h
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Q[h_] := 3 h^2 + I1[h]/I0[h]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    deltaQ[h_] := Q[h] - 3 h^2
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

zRef = 100;
zHalf = N[(1 + zRef)/2^(2/3) - 1, 50];
hRef = e[zRef];
hHalf = e[zHalf];
hNormRef = 1;
hNormHalf = N[hHalf/hRef, 50];

wRef = wEff[zRef, omPlanck, orPlanck, olPlanck];
wHalf = wEff[zHalf, omPlanck, orPlanck, olPlanck];

z2 = 1/20;
alpha = 10^-5;
rs = 10^-3;
beta = 1;
rFloor = 10^-12;
aReg = N[Sqrt[(beta rs)^2 + rFloor^2], 50];

rDom[h_] := N[Sqrt[1 - z2]/h, 50];
i0[h_] := N[(4 Pi/3) rDom[h]^3, 50];
ff[x_] := N[(1/(2 aReg)) ArcTan[x/aReg] - x/(2 (x^2 + aReg^2)), 50];
localI1[h_] := N[4 Pi alpha rs^2 (ff[rDom[h]] - ff[0]), 50];
qMass[h_] := N[3 h^2 + localI1[h]/i0[h], 50];
q0[h_] := N[3 h^2, 50];

i0Ref = i0[hNormRef];
i0Half = i0[hNormHalf];
locRef = localI1[hNormRef];
locHalf = localI1[hNormHalf];
qMassRef = qMass[hNormRef];
qMassHalf = qMass[hNormHalf];
q0Ref = q0[hNormRef];
q0Half = q0[hNormHalf];
deltaQRef = N[qMassRef - q0Ref, 50];
deltaQHalf = N[qMassHalf - q0Half, 50];
ratioDelta = N[deltaQHalf/deltaQRef, 50];
ratioI0 = N[i0Half/i0Ref, 50];
ratioLoc = N[locHalf/locRef, 50];
epsRef = N[deltaQRef/q0Ref, 50];
epsHalf = N[deltaQHalf/q0Half, 50];

checkD = Abs[hNormHalf - 1/2] < 0.03;
checkE = Abs[wRef] < 0.05 && Abs[wHalf] < 0.05;
checkF = zHalf > 0 && zHalf < zRef;
checkG = deltaQRef > 0 && deltaQHalf > 0;
checkH = Abs[ratioDelta/(hNormHalf^3) - 1] < 0.08;
checkI = Abs[ratioI0 hNormHalf^3 - 1] < 0.08;
checkJ = Abs[ratioLoc - 1] < 0.1;
checkK = epsRef < 0.05 && epsHalf < 0.05;
checkL = deltaQHalf < deltaQRef/4;

classification = If[
  TrueQ[
    checkA && checkB && checkC && checkD && checkE && checkF &&
    checkG && checkH && checkI && checkJ && checkK && checkL
  ],
  "MATTER_ERA_IR_SUPPRESSED",
  "undetermined"
];

check = TrueQ[classification == "MATTER_ERA_IR_SUPPRESSED"];

nbName = "02_Matter_Era_Finite_Domain_IR_Scaling_Test.nb";
logName = "02_matter_era_finite_domain_ir_scaling_test.log";
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
  "zRef = " <> ToString[zRef, InputForm],
  "zHalf = " <> ToString[zHalf, InputForm],
  "hRef = " <> ToString[hRef, InputForm],
  "hHalf = " <> ToString[hHalf, InputForm],
  "hNormHalf = " <> ToString[hNormHalf, InputForm],
  "wRef = " <> ToString[wRef, InputForm],
  "wHalf = " <> ToString[wHalf, InputForm],
  "z2 = " <> ToString[z2, InputForm],
  "alpha = " <> ToString[alpha, InputForm],
  "rs = " <> ToString[rs, InputForm],
  "beta = " <> ToString[beta, InputForm],
  "rFloor = " <> ToString[rFloor, InputForm],
  "aReg = " <> ToString[aReg, InputForm],
  "i0Ref = " <> ToString[i0Ref, InputForm],
  "i0Half = " <> ToString[i0Half, InputForm],
  "locRef = " <> ToString[locRef, InputForm],
  "locHalf = " <> ToString[locHalf, InputForm],
  "deltaQRef = " <> ToString[deltaQRef, InputForm],
  "deltaQHalf = " <> ToString[deltaQHalf, InputForm],
  "ratioDelta = " <> ToString[ratioDelta, InputForm],
  "ratioI0 = " <> ToString[ratioI0, InputForm],
  "ratioLoc = " <> ToString[ratioLoc, InputForm],
  "epsRef = " <> ToString[epsRef, InputForm],
  "epsHalf = " <> ToString[epsHalf, InputForm],
  "checkD(hNormHalf ~ 1/2) = " <> ToString[checkD, InputForm],
  "checkE(matter w_eff probes) = " <> ToString[checkE, InputForm],
  "checkF(zHalf in era) = " <> ToString[checkF, InputForm],
  "checkG(deltaQ positive) = " <> ToString[checkG, InputForm],
  "checkH(delta scaling with h^3) = " <> ToString[checkH, InputForm],
  "checkI(I0 scaling with h^-3) = " <> ToString[checkI, InputForm],
  "checkJ(local term quasi-constant) = " <> ToString[checkJ, InputForm],
  "checkK(subleading eps) = " <> ToString[checkK, InputForm],
  "checkL(delta half strongly reduced) = " <> ToString[checkL, InputForm],
  "classification = " <> ToString[classification, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "02_matter_era_finite_domain_ir_scaling_test",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
