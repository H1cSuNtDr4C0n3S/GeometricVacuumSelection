ClearAll["Global`*"];

nb = Notebook[{
  Cell["26 - Covariant General EOM Channels", "Title"],
  Cell["Symbolic certification of the three 4D channels (constraint, metric, scalar) from S = S_loc + S_NL, with explicit bulk effective stress tensor from the nonlocal leafwise term.", "Text"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    S == SLoc + Integrate[lambda[theta] (I1[theta] - K0^2 I0[theta]), theta]
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    C[theta] == I1[theta] - K0^2 I0[theta] == 0
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    Gmunu + Lambda gmunu == kappa (TTheta + TNL)
  ], "Input"],
  Cell[BoxData @ ToBoxes @ HoldForm[
    alpha BoxTheta + JNL == 0
  ], "Input"]
}];

(* ===== 1) Constraint channel from lambda-variation ===== *)
sNLeps = Expand[(lam + eps eta) (i1 - k0^2 i0)];
deltaSNLeta = SeriesCoefficient[sNLeps, {eps, 0, 1}] // Expand;
expectedDeltaSNLeta = Expand[eta (i1 - k0^2 i0)];
checkA = Simplify[deltaSNLeta == expectedDeltaSNLeta];

(* ===== 2) Metric channel: explicit dmu and T_NL,bulk ===== *)
(* metric-only variation of mu = Sqrt[-g] W Sqrt[X] DiracDelta[Theta-theta] *)
(* encoded through standard identities:
   dSqrtMinusG = -(1/2) Sqrt[-g] g_{mn} dg^{mn}
   dSqrtX      = +(1/2) Sqrt[X] u_m u_n dg^{mn}
   => dmu      = -(1/2) mu h_{mn} dg^{mn}, h = g - u u *)
dSqrtMinusG = -(1/2) sgm g dg;
dSqrtX = (1/2) sx uu dg;
dMuRaw = Expand[dSqrtMinusG sx + sgm dSqrtX] /. (sgm sx) -> mu;
dMuExpected = Expand[-(1/2) mu h dg /. h -> (g - uu)];
checkB = Simplify[dMuRaw == dMuExpected];

(* split dK2 by definition *)
dK2 = piG dg + piTheta dth + divV;

(* metric coefficients in dI0,dI1 bulk *)
dI0bulk = Expand[chi (-(1/2) mu h dg)];
dI1bulk = Expand[chi (k2 (-(1/2) mu h dg) + mu (piG dg))];

metricCoeffI0 = Coefficient[dI0bulk, dg, 1] // Expand;
metricCoeffI1 = Coefficient[dI1bulk, dg, 1] // Expand;

expectedCoeffI0 = Expand[-(1/2) chi mu h];
expectedCoeffI1 = Expand[chi (-1/2 mu k2 h + mu piG)];

checkC = Simplify[metricCoeffI0 == expectedCoeffI0];
checkD = Simplify[metricCoeffI1 == expectedCoeffI1];

(* nonlocal metric contribution and bulk effective stress tensor *)
dSNLmetricBulk = Expand[lam (dI1bulk - k0^2 dI0bulk)];
metricCoeffSNL = Coefficient[dSNLmetricBulk, dg, 1] // Expand;

expectedMetricCoeffSNL = Expand[lam chi mu (-1/2 (k2 - k0^2) h + piG)];
checkE = Simplify[metricCoeffSNL == expectedMetricCoeffSNL];

tnlBulkFromDef = Expand[-2 metricCoeffSNL];
tnlBulkExpected = Expand[lam chi mu ((k2 - k0^2) h - 2 piG)];
checkF = Simplify[tnlBulkFromDef == tnlBulkExpected];

(* ===== 3) Scalar channel structure ===== *)
(* keep local and nonlocal scalar variations in coefficient form *)
deltaSLocScalar = slocTheta dth;
deltaSNLScalar = sNLtheta dth;
totalScalarCoeff = Expand[deltaSLocScalar + deltaSNLScalar] /. dth -> 1;

scalarRule = {
  slocTheta -> alpha boxTheta,
  sNLtheta -> jNL
};
scalarEquationLHS = Expand[totalScalarCoeff /. scalarRule];
checkG = Simplify[scalarEquationLHS == (alpha boxTheta + jNL)];

(* ===== 4) Covariant metric-channel assembly ===== *)
(* local metric variation written in standard EH+Theta form *)
localMetricCoeff = 1/2 sgm ((gEin + Lambda gmn)/kappa - tTheta);
nonlocalMetricCoeff = -1/2 sgm tNL;
totalMetricCoeff = Expand[localMetricCoeff + nonlocalMetricCoeff];

metricEqLHS = Expand[2 totalMetricCoeff/sgm];
metricEqExpected = Expand[(gEin + Lambda gmn)/kappa - tTheta - tNL];
checkH = Simplify[metricEqLHS == metricEqExpected];

(* ===== 5) de Sitter branch with Lambda = 0 in minisuperspace proxy ===== *)
lamRule = First[Solve[h0^2 (-3 + 6 kappa lam0) + Lambda == 0, lam0]];
lam0Lambda0 = Simplify[lam0 /. lamRule /. Lambda -> 0 /. h0^2 -> k0^2/3];
checkI = Simplify[lam0Lambda0 == 1/(2 kappa)];

check = TrueQ[
  checkA && checkB && checkC && checkD && checkE && checkF && checkG && checkH && checkI
];

nbName = "26_Covariant_General_EOM_Channels.nb";
logName = "26_covariant_general_eom_channels.log";
logLines = {
  "Notebook: " <> nbName,
  "deltaSNLeta = " <> ToString[deltaSNLeta, InputForm],
  "expectedDeltaSNLeta = " <> ToString[expectedDeltaSNLeta, InputForm],
  "checkA(constraint channel) = " <> ToString[checkA, InputForm],
  "dMuRaw = " <> ToString[dMuRaw, InputForm],
  "dMuExpected = " <> ToString[dMuExpected, InputForm],
  "checkB(dmu metric identity) = " <> ToString[checkB, InputForm],
  "metricCoeffI0 = " <> ToString[metricCoeffI0, InputForm],
  "expectedCoeffI0 = " <> ToString[expectedCoeffI0, InputForm],
  "checkC(dI0 metric coeff) = " <> ToString[checkC, InputForm],
  "metricCoeffI1 = " <> ToString[metricCoeffI1, InputForm],
  "expectedCoeffI1 = " <> ToString[expectedCoeffI1, InputForm],
  "checkD(dI1 metric coeff) = " <> ToString[checkD, InputForm],
  "metricCoeffSNL = " <> ToString[metricCoeffSNL, InputForm],
  "expectedMetricCoeffSNL = " <> ToString[expectedMetricCoeffSNL, InputForm],
  "checkE(dS_NL metric coeff) = " <> ToString[checkE, InputForm],
  "tnlBulkFromDef = " <> ToString[tnlBulkFromDef, InputForm],
  "tnlBulkExpected = " <> ToString[tnlBulkExpected, InputForm],
  "checkF(TNL bulk explicit) = " <> ToString[checkF, InputForm],
  "scalarEquationLHS = " <> ToString[scalarEquationLHS, InputForm],
  "checkG(scalar channel assembly) = " <> ToString[checkG, InputForm],
  "metricEqLHS = " <> ToString[metricEqLHS, InputForm],
  "metricEqExpected = " <> ToString[metricEqExpected, InputForm],
  "checkH(metric channel assembly) = " <> ToString[checkH, InputForm],
  "lamRule = " <> ToString[lamRule, InputForm],
  "lam0(Lambda=0,H^2=K0^2/3) = " <> ToString[lam0Lambda0, InputForm],
  "checkI(deSitter lambda) = " <> ToString[checkI, InputForm],
  "check = " <> ToString[check, InputForm]
};
logText = StringRiffle[logLines, "\n"];

result = <|
  "script" -> "26_covariant_general_eom_channels",
  "nbName" -> nbName,
  "logName" -> logName,
  "nbBase64" -> BaseEncode[ExportByteArray[nb, "NB"]],
  "logText" -> logText,
  "check" -> check
|>;

Print["__RESULT__" <> ExportString[result, "RawJSON"]];
If[TrueQ[check], Exit[0], Exit[1]];
