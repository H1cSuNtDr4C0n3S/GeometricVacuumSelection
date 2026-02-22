# Variational Derivation - Detailed Calculations and Results

Questo documento entra nel dettaglio dei calcoli eseguiti in `Paper/Variational_Derivation/scripts/` e dei risultati registrati in `Paper/Variational_Derivation/logs/`.

Convenzioni usate nei notebook:
1. `I0` = integrale di misura/volume sulla foglia.
2. `I1` = integrale del contenuto geometrico (nel flusso interpretato come media di `K^2` pesata).
3. `Q = I1/I0`.
4. `lambda` = moltiplicatore di Lagrange del vincolo.
5. `theta` = etichetta di foliazione (tempo geometrico).

## 00 - Variational Seed

Riferimento: `scripts/00_variational_seed.wl`, log `logs/00_variational_seed.log`.

Calcolo fatto:
1. Si definisce la densita variazionale minima
   `lagDensity[th] = lambda[th] ( I1[th] - K0^2 I0[th] )`.
2. Si deriva rispetto a `lambda[th]`.
3. Si controlla il termine lineare in una perturbazione `lambda -> lambda + eps eta`.

Output:
1. `eqLambda = -(K0^2*I0[th]) + I1[th]`.
2. `lin = -(K0^2*eta[th]*I0[th]) + eta[th]*I1[th]`.
3. `check1=True`, `check2=True`.

Conclusione:
1. La variazione in `lambda` impone direttamente il vincolo globale leafwise `I1 - K0^2 I0 = 0`.
2. Il settore `lambda` non introduce dinamica propria locale: enforce puro del vincolo.

## 01 - Ratio First Variation

Riferimento: `scripts/01_ratio_firstvariation.wl`, log `logs/01_ratio_firstvariation.log`.

Calcolo fatto:
1. Si perturba `(i1 + eps di1)/(i0 + eps di0)`.
2. Si estrae il coefficiente lineare in `eps`.
3. Si confronta con la formula analitica del quoziente.

Output:
1. `series1 = (di1*eps*i0 - di0*eps*i1 + i0*i1)/i0^2`.
2. `deltaRatio = (di1*i0 - di0*i1)/i0^2`.
3. `expected = (di1*i0 - di0*i1)/i0^2`.
4. `check=True`.

Conclusione:
1. La variazione di `Q=I1/I0` contiene sempre differenza tra contributo numeratore e contributo denominatore.
2. E' la formula base che ricompare in tutti i passaggi successivi.

## 02 - Constraint First Variation

Riferimento: `scripts/02_constraint_firstvariation.wl`, log `logs/02_constraint_firstvariation.log`.

Calcolo fatto:
1. Si varia l integrando vincolato tenendo separati i canali `eta`, `dI0`, `dI1`.
2. Si riscrive anche in termini di `dQ`.
3. Si confronta con forma attesa.

Output:
1. `deltaIntegrand = eta[th]*(-(K0^2*I0[th]) + I1[th]) + (-(K0^2*dI0[th]) + dI1[th])*lambda[th]`.
2. `fromIQ = eta[th]*(-(K0^2*I0[th]) + I1[th]) + (dQ[th]*I0[th] + dI0[th]*(-K0^2 + I1[th]/I0[th]))*lambda[th]`.
3. `checkA=True`, `checkB=True`.

Conclusione:
1. Separazione pulita tra enforcement (`eta` contro il vincolo) e risposta geometrica (`dI0,dI1`).
2. Conferma che usare `Q` o usare direttamente `I0,I1` e algebraicamente equivalente.

## 03 - Theta Chain Rule

Riferimento: `scripts/03_theta_chainrule.wl`, log `logs/03_theta_chainrule.log`.

Calcolo fatto:
1. Si deriva `f(th) = lambda(th)(I1(th)-K0^2 I0(th))` rispetto a `th`.
2. Si moltiplica per `dTheta(th)` per ottenere `deltaTheta`.
3. Si verifica identita completa della chain rule.

Output:
1. `dfdth = lambda[th]*(-(K0^2*Derivative[1][I0][th]) + Derivative[1][I1][th]) + (-(K0^2*I0[th]) + I1[th])*Derivative[1][lambda][th]`.
2. `deltaTheta = dTheta[th]*(...)` con gli stessi termini.
3. `checkA=True`, `checkB=True`.

Conclusione:
1. Il contributo di `theta` non e trascurabile: entra sia tramite `I0,I1` sia tramite `lambda`.
2. Struttura coerente con una foliazione dinamica.

## 04 - Ratio Kernel Form

Riferimento: `scripts/04_ratio_kernel_form.wl`, log `logs/04_ratio_kernel_form.log`.

Calcolo fatto:
1. Si assume forma a kernel per le variazioni:
   `dI0 = Int[dg(x) k0(x)]`, `dI1 = Int[dg(x) k1(x)]`.
2. Si inserisce nella formula del quoziente.
3. Si riscrive in forma espansa per isolare i pesi.

Output:
1. `deltaQ = (-(I1*Int[dg k0]) + I0*Int[dg k1])/I0^2`.
2. `expectedKernelExpanded = -((I1*Int[dg k0])/I0^2) + Int[dg k1]/I0`.
3. `checkA=True`, `checkB=True`.

Conclusione:
1. La risposta variazionale ha la forma differenza tra kernel locale e media globale.
2. Questa e la base per leggere il regime IR come filtro globale.

## 05 - Metric Variation Spherical

Riferimento: `scripts/05_metric_variation_spherical.wl`, log `logs/05_metric_variation_spherical.log`.

Calcolo fatto:
1. Si usa la forma sferica di `K2_general`.
2. Si valutano i limiti/background noti (Minkowski e de Sitter).
3. Si confronta con le espressioni attese.

Output:
1. `k2Minkowski = (-2*v^2)/(r^2*(-1 + v^2))`.
2. `k2dS = 3*H^2`.
3. `checkA=True`, `checkB=True`, `checkC=True`.

Conclusione:
1. In settore Minkowski il contributo effettivo scala come `1/r^2` e non resta costante in IR.
2. In settore de Sitter il valore e costante positivo (`3H^2`), quindi compatibile con un vincolo a media costante.

## 06 - Second-Order EL Identity

Riferimento: `scripts/06_second_order_EL_identity.wl`, log `logs/06_second_order_EL_identity.log`.

Calcolo fatto:
1. Si costruisce un test polinomiale con dipendenza da derivate seconde.
2. Si applica identita di Eulero-Lagrange con termini di bordo.
3. Si controlla che il residuo sia nullo.

Output:
1. `EL = 2*c0*q0c - 4*c1*q2c + 48*c2*q4c + ... + 2*c0*q4c*x^4`.
2. `residue = 0`.
3. `check=True`.

Conclusione:
1. Il formalismo con derivate seconde e implementato correttamente.
2. Nessuna perdita di termini di bordo nel passaggio a EL.

## 07 - K2 Psi Variation Numeric

Riferimento: `scripts/07_k2_psi_variation_numeric.wl`, log `logs/07_k2_psi_variation_numeric.log`.

Calcolo fatto:
1. Si applica la struttura variazionale al blocco reale `K2_general` nel canale `psi`.
2. Si testa su punti campione radiali.
3. Si confronta differenza numerica vs formula analitica.

Output:
1. `samplePts = {0.7, 0.9, 1.1, 1.3}`.
2. `residualValues = {0, 0, 0, 0}`.
3. `maxResidual = 0`, `check=True`.

Conclusione:
1. Identita variazionale soddisfatta esattamente nei campioni testati.
2. Validazione operativa del canale `psi`.

## 08 - K2 Psi EL Operator

Riferimento: `scripts/08_k2_psi_EL_operator.wl`, log `logs/08_k2_psi_EL_operator.log`.

Calcolo fatto:
1. Si costruisce esplicitamente l operatore `EL_psi` per `K2_general`.
2. Si campiona l operatore su quattro raggi.
3. Si verifica anche la consistenza differenziale (`residualValues`).

Output:
1. `ELvalues = {2.4306327483237995, 1.1798451310566078, 0.6699496788795586, 0.4263114587001963}`.
2. `ELleafCount = 209`.
3. `residualValues = {0, 0, 0, 0}`, `maxResidual = 0`, `check=True`.
4. `ELpreview` e un polinomio razionale lungo (troncato nel log).

Conclusione:
1. L operatore dinamico di `psi` e non banale e decresce radialmente nel background testato.
2. La forma lunga (`ELpreview`) indica struttura nonlineare ricca, non riducibile a pochi termini.

## 09 - Metric Kernels A/B Numeric

Riferimento: `scripts/09_metric_kernels_AB_numeric.wl`, log `logs/09_metric_kernels_AB_numeric.log`.

Calcolo fatto:
1. Si estraggono i kernel metrici efficaci nei canali `A` e `B`.
2. Si valutano numericamente su punti campione.
3. Si verifica residuo variazionale nullo.

Output:
1. `EAvalues = {0.164714498102404, 0.10136972652898028, 0.06597158587489217, 0.04249281815249714}`.
2. `EBvalues = {-0.18750122813491985, -0.1296360616888789, -0.09978467003444905, -0.08258512679849422}`.
3. `EAleafCount = 199`, `EBleafCount = 191`.
4. `residualValues = {0, 0, 0, 0}`, `maxResidual = 0`, `check=True`.
5. `EApreview` e `EBpreview` mostrano operatori polinomiali lunghi (troncati nel log).

Conclusione:
1. I canali metrici hanno risposta differenziata (segni e scale diverse).
2. Il vincolo non agisce come riscalatura uniforme della metrica.

## 10 - Full Proxy EOM Operators

Riferimento: `scripts/10_full_proxy_eom_operators.wl`, log `logs/10_full_proxy_eom_operators.log`.

Calcolo fatto:
1. Si assembla una densita proxy completa (locale + termine di vincolo).
2. Si costruiscono gli operatori EOM dei tre canali `A`, `B`, `psi`.
3. Si campionano i valori su `r = {0.8,1.0,1.2}`.

Output:
1. `EAvalues = {-2.7927302932294946, -2.676808775860709, -2.539069081983486}`.
2. `EBvalues = {-2.356199205071267, -2.5444503443257043, -2.7757012519239197}`.
3. `Epsivalues = {-1.247200893332045, -1.8565167768369895, -2.517655494263335}`.
4. `EAleafCount = 1309`, `EBleafCount = 1117`, `EpsileafCount = 1127`.
5. `check=True`.

Conclusione:
1. Il sistema proxy attiva simultaneamente tutti i canali dinamici.
2. Le EOM risultanti sono estese (conteggio termini elevato), coerente con una teoria efficace nonlineare.

## 11 - Full Proxy Variational Checks

Riferimento: `scripts/11_full_proxy_variational_checks.wl`, log `logs/11_full_proxy_variational_checks.log`.

Calcolo fatto:
1. Si verifica canale per canale la decomposizione
   `deltaL = E_field * delta(field) + total derivative`.
2. Si campionano i residui numerici su `r = {0.8,1.0,1.2}`.
3. Si controlla soglia di consistenza.

Output:
1. `resAvals = {-1.3361305999372705*^-15, -1.2302681317945482*^-15, -1.1229065778314945*^-15}`.
2. `resBvals = {3.0888261575352174*^-16, 4.106390994612664*^-16, 5.248909667811335*^-16}`.
3. `resPsivals = {1.8789636663694201*^-16, 2.556792130909114*^-16, 3.3869698479312834*^-16}`.
4. `maxResA = 1.3361305999372705*^-15`, `maxResB = 5.248909667811335*^-16`, `maxResPsi = 3.3869698479312834*^-16`.
5. `checkA=True`, `checkB=True`, `checkPsi=True`, `check=True`.

Conclusione:
1. Residui a livello macchina (`10^-15`-`10^-16`): coerenza variazionale piena del proxy.
2. Nessun mismatch operativo tra canali metrici e canale `psi`.

## 12 - Finite-Domain Ratio Variation

Riferimento: `scripts/12_finite_domain_ratio_variation.wl`, log `logs/12_finite_domain_ratio_variation.log`.

Calcolo fatto:
1. Si introduce dominio finito con bordo mobile: `R -> R + eps rho`.
2. Si definiscono
   `deltaI0 = dI0bulk + rho f0(R)`,
   `deltaI1 = dI1bulk + rho f1(R)`.
3. Si espande `Q = I1/I0` al primo ordine.
4. Si isola il pezzo di bordo.
5. Si fa check numerico indipendente con derivata centrale da integrali espliciti.

Output:
1. `deltaQ = (dI1bulk*I0 - dI0bulk*I1 + f1R*I0*rho - f0R*I1*rho)/I0^2`.
2. `boundaryPiece = (f1R*I0 - f0R*I1)/I0^2`.
3. `deltaNum = 0.0871228388174273043170618346112253574359857178106621173916`.
4. `deltaPred = 0.0871228388174262140273603379268826091806459542071834440351`.
5. `absErr = 1.0902897014966843427482553397636034786733565*10^-15`.
6. `checkA=True`, `checkB=True`, `checkNum=True`, `check=True`.

Conclusione:
1. La variazione a dominio finito separa in modo controllato contributi bulk e contributi bordo.
2. Il bordo pesa secondo contrasto tra valore locale e media globale.
3. L accordo numerico a `10^-15` conferma formula analitica robusta.

## 13 - IR Flow Minkowski vs de Sitter

Riferimento: `scripts/13_ir_flow_minkowski_desitter.wl`, log `logs/13_ir_flow_minkowski_desitter.log`.

Calcolo fatto:
1. Si deriva la legge di flusso radiale
   `dQ/dR = (f1(R) I0(R) - f0(R) I1(R))/I0(R)^2`.
2. Si applica a due toy sectors:
   `Minkowski-like`: `f0=r^2`, `f1=c`.
   `de Sitter-like`: `f1=KdS f0`.
3. Si verifica formula analitica e campionamento numerico.

Output:
1. `qPrime = (f1R*I0 - f0R*I1)/I0^2`.
2. `qM(R) = (3*c)/R^2`, `dqM(R) = (-6*c)/R^3`, `checkM=True`.
3. `qD(R) = KdS`, `dqD(R) = 0`, `checkD=True`.
4. `sampleRs = {2, 4, 8}`.
5. `qMvals(c=0.8) = {0.6000000000000001, 0.15000000000000002, 0.037500000000000006}`.
6. `qDvals(KdS=1.2) = {1.2, 1.2, 1.2}`.
7. `checkA=True`, `check=True`.

Conclusione:
1. In settore Minkowski-like il rapporto medio decade come `1/R^2` e tende a zero in IR.
2. In settore de Sitter-like il rapporto resta costante.
3. Questa differenza realizza in forma operativa il criterio di selezione IR tra ramo statico e ramo de Sitter.

## 14 - Constraint Finite-Domain Boundary

Riferimento: `scripts/14_constraint_finite_domain_boundary.wl`, log `logs/14_constraint_finite_domain_boundary.log`.

Calcolo fatto:
1. Si considera il funzionale lineare del vincolo:
   `C = I1 - K0^2 I0`.
2. Su dominio finito con bordo mobile:
   `deltaI0 = dI0bulk + rho f0(R)`,
   `deltaI1 = dI1bulk + rho f1(R)`.
3. Si separa `deltaC` in parte bulk e parte bordo.
4. Si verifica numericamente con derivata centrale indipendente.

Output:
1. `deltaC = dI1bulk - dI0bulk*K0^2 + f1R*rho - f0R*K0^2*rho`.
2. `bulkPart = dI1bulk - dI0bulk*K0^2`.
3. `boundaryPart = f1R - f0R*K0^2`.
4. `checkA=True`, `checkB=True`, `checkBoundaryNeutral=True`.
5. `deltaNum = 0.1006720000000007869375...`.
6. `deltaPred = 0.100672`.
7. `absErr = 7.869375...*10^-16`.
8. `checkNum=True`, `check=True`.

Conclusione:
1. Anche nel canale lineare del vincolo il contributo di bordo e separabile in modo esatto.
2. Il bordo e inattivo quando soddisfa localmente il matching `f1(R)=K0^2 f0(R)`.
3. Consistenza numerica a livello macchina.

## 15 - FRW K2 Leafwise Average

Riferimento: `scripts/15_frw_k2_leafwise_average.wl`, log `logs/15_frw_k2_leafwise_average.log`.

Calcolo fatto:
1. Riduzione FRW omogenea:
   `K2 = 3 (ad/(a N))^2`.
2. Calcolo dei funzionali su palla comovente:
   `I0 = (4 Pi/3) a^3 R^3`,
   `I1 = K2 I0`.
3. Verifica che `Q=I1/I0=K2`.
4. Specializzazione de Sitter `a=e^(H t), N=1`.

Output:
1. `K2 = (3*ad^2)/(a^2*N^2)`, `checkA=True`.
2. `I0 = (4*a^3*Pi*R^3)/3`.
3. `I1 = (4*a*ad^2*Pi*R^3)/N^2`.
4. `Q = (3*ad^2)/(a^2*N^2)`, `checkB=True`.
5. `QdS = 3*H^2`, `dQdt_dS = 0`, `checkC=True`.
6. `selectionRule(H^2) = {h2 -> K0^2/3}`, `checkD=True`.
7. `numErr = 0`, `checkNum=True`, `check=True`.

Conclusione:
1. In FRW omogeneo il vincolo leafwise coincide con il valore locale di `K2`.
2. Il ramo de Sitter soddisfa in modo diretto la selezione `H^2=K0^2/3`.

## 16 - FRW Minisuperspace EOM

Riferimento: `scripts/16_frw_minisuperspace_eom.wl`, log `logs/16_frw_minisuperspace_eom.log`.

Calcolo fatto:
1. Si definisce il lagrangiano ridotto:
   `L = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3)`.
2. Si derivano le tre EOM:
   `E_a`, `E_N`, `E_lambda`.
3. Si testa il ramo de Sitter:
   `a=e^(H0 t), N=1, H0^2=K0^2/3`.
4. Si risolve il valore compatibile di `lambda0` dal canale `E_N`.

Output:
1. `EaDS = (-3*E^(2*H0*t)*(K0^2*kappa*lam0 + 3*H0^2*(-1 + kappa*lam0) + Lambda))/kappa`.
2. `ENDS = -((E^(3*H0*t)*(H0^2*(-3 + 6*kappa*lam0) + Lambda))/kappa)`.
3. `ElamDS = E^(3*H0*t)*(3*H0^2 - K0^2)`.
4. `lamRule = {lam0 -> (K0^2 - Lambda)/(2*K0^2*kappa)}`.
5. `checkA(Elam)=True`, `checkB(EN)=True`, `checkC(Ea)=True`.
6. `eaNum=0`, `enNum=0`, `elamNum=0`, `maxAbsNum=0`, `checkNum=True`, `check=True`.

Conclusione:
1. Le EOM minisuperspace sono chiuse e coerenti.
2. Il ramo de Sitter e una soluzione simultanea dei tre canali, con `lambda0` fissato da `K0`, `Lambda`, `kappa`.

## 17 - FRW Linearized Constraint Channel

Riferimento: `scripts/17_frw_linearized_constraint_channel.wl`, log `logs/17_frw_linearized_constraint_channel.log`.

Calcolo fatto:
1. Si linearizza attorno a de Sitter:
   `a = e^(H0 t) (1 + eps u(t))`, `N = 1 + eps n(t)`.
2. Si espande il vincolo:
   `C = 3 a a'^2/N^2 - K0^2 a^3 = C0 + eps C1 + ...`.
3. Si impone il ramo `H0^2=K0^2/3`.
4. Si verifica il canale linearizzato anche numericamente su campioni.

Output:
1. `C0 = E^(3*H0*t)*(3*H0^2 - K0^2)`.
2. `C1 = -3*E^(3*H0*t)*(2*H0^2*n[t] + (-3*H0^2 + K0^2)*u[t] - 2*H0*Derivative[1][u][t])`.
3. `C1_onBranch = -2*E^(3*H0*t)*(K0^2*n[t] - 3*H0*Derivative[1][u][t])`.
4. `C1_gauge_n0 = 6*E^(3*H0*t)*H0*Derivative[1][u][t]`.
5. `checkA=True`, `checkB=True`, `checkC=True`, `checkD=True`.
6. `c1Vals(constrained) = {0, 0, 0, 0}`, `maxAbs=0`, `checkNum=True`, `check=True`.

Conclusione:
1. Al primo ordine il vincolo lega direttamente perturbazioni metriche/scalari:
   `u'(t) = H0 n(t)`.
2. In gauge cosmica (`n=0`) il canale vincolato e congelato (`u'=0`) a livello lineare.

## 18 - de Sitter Constraint Propagation

Riferimento: `scripts/18_desitter_constraint_propagation.wl`, log `logs/18_desitter_constraint_propagation.log`.

Calcolo fatto:
1. Si riparte dal lagrangiano minisuperspace FRW del notebook 16.
2. Si valutano `Ea`, `EN`, `Elam` sul ramo de Sitter con
   `H0^2=K0^2/3` e `lam0=(K0^2-Lambda)/(2*K0^2*kappa)`.
3. Si controllano anche le derivate temporali dei tre canali on-shell.
4. Si fa verifica numerica su punti campione.

Output:
1. `EaOn = 0`, `ENOn = 0`, `ElamOn = 0`.
2. `dEaOn = 0`, `dENOn = 0`, `dElamOn = 0`.
3. `checkA=True`, `checkB=True`, `checkC=True`, `checkD=True`, `checkE=True`, `checkF=True`.
4. `eaVals`, `enVals`, `elVals` numericamente nulli.
5. `maxAbsNum = 0`, `checkNum=True`, `check=True`.

Conclusione:
1. I canali di vincolo e dinamica restano consistenti lungo il flusso temporale sul ramo de Sitter.
2. Non compaiono drift temporali dei vincoli nel settore omogeneo.

## 19 - FRW Linearized Full EOM Reduction

Riferimento: `scripts/19_frw_linearized_full_eom_reduction.wl`, log `logs/19_frw_linearized_full_eom_reduction.log`.

Calcolo fatto:
1. Si linearizzano simultaneamente `Ea`, `EN`, `Elam` con
   `a=Exp[H0 t](1+eps u)`, `N=1+eps n`, `lambda=lam0+eps l1`.
2. Si impone il ramo de Sitter (`H0=K0/Sqrt[3]`) e il valore background di `lam0`.
3. Da `Elam^(1)=0` si ricava `nRule`.
4. Da `EN^(1)=0` si ricava `lRule`.
5. Si sostituiscono anche `dnRule` e `dlRule` in `Ea^(1)` per la riduzione completa.

Output:
1. `nRule = {n[t] -> (Sqrt[3]*u'[t])/K0}`, `checkA=True`.
2. `lRule = {l1[t] -> (Sqrt[3]*(K0^2-Lambda)*u'[t])/(2*K0^3*kappa)}`, `checkB=True`.
3. `eaReduced = 0`.
4. `u2Rule = {}`, `u2Expr = u''[t]` (equazione ridotta identicamente soddisfatta).
5. `resVals = {0,0,0}`, `maxAbs=0`, `checkNum=True`, `check=True`.

Conclusione:
1. Nel proxy lineare FRW, i canali vincolati eliminano completamente le variabili ausiliarie al primo ordine.
2. Sul ramo de Sitter omogeneo non rimane una dinamica lineare residua nel canale scalare di questo settore ridotto.

## 20 - Minisuperspace Noether-Bianchi Identity

Riferimento: `scripts/20_minisuperspace_noether_bianchi_identity.wl`, log `logs/20_minisuperspace_noether_bianchi_identity.log`.

Calcolo fatto:
1. Si considerano le EOM minisuperspace complete `Ea`, `EN`, `Elam`.
2. Si costruisce la combinazione:
   `noetherRaw = Ea a' + EN N' + Elam lambda' - d/dt(N EN)`.
3. Si identifica il termine di bilancio del moltiplicatore:
   `auxTerm = d/dt(lambda Elam)`.
4. Si verifica l'identita' off-shell:
   `noetherId = noetherRaw - auxTerm = 0`.
5. Si verifica la forma equivalente:
   `d/dt(N EN + lambda Elam) - (Ea a' + EN N' + Elam lambda') = 0`.

Output:
1. `noetherRaw = ...` (espressione completa nel log).
2. `auxTerm = ...` (uguale a `noetherRaw`).
3. `noetherId = 0`, `checkA=True`.
4. `propEq = 0`, `checkB=True`.
5. `idDS = 0`, `checkC=True`.
6. `idVals = {0,0,0,0}`, `maxAbs=0`, `checkNum=True`, `check=True`.

Conclusione:
1. Esiste una relazione off-shell esatta che vincola i tre canali EOM minisuperspace.
2. Il settore del moltiplicatore non rompe la consistenza: compare come derivata totale `d(lambda Elam)/dt`.
3. Questo giustifica in modo strutturale la propagazione dei vincoli osservata nei notebook 18-19.

## 21 - FRW Curvature Minisuperspace EOM

Riferimento: `scripts/21_frw_curvature_minisuperspace_eom.wl`, log `logs/21_frw_curvature_minisuperspace_eom.log`.

Calcolo fatto:
1. Si estende il lagrangiano minisuperspace con curvatura spaziale:
   `Lk = Lflat + 3 k a N/kappa`.
2. Si ricalcolano `EaK`, `ENK`, `ElamK`.
3. Si verifica il limite `k->0` rispetto alle EOM del notebook 16.
4. Si controlla la sensibilita' differenziale ai termini di curvatura.

Output:
1. `checkA(k->0 Ea)=True`, `checkB(k->0 EN)=True`, `checkC(k->0 Elam)=True`.
2. `checkD(dEa/dk)=True`, `checkE(dEN/dk)=True`.
3. `checkF(Elam independent of k)=True`.
4. Campioni numerici finiti (`eaVals`, `enVals`, `elVals`), `checkNum=True`, `check=True`.

Conclusione:
1. L'estensione con curvatura e' coerente e recupera esattamente il caso piatto come limite.
2. La curvatura entra soltanto nei canali metrici, non nel vincolo `Elam`.

## 22 - FRW Curvature Noether Identity

Riferimento: `scripts/22_frw_curvature_noether_identity.wl`, log `logs/22_frw_curvature_noether_identity.log`.

Calcolo fatto:
1. Si ripete il test di identita' strutturale del notebook 20 nel caso `k != 0`.
2. Si definiscono:
   `noetherRaw = Ea a' + EN N' + Elam lambda' - d/dt(N EN)`,
   `noetherId = noetherRaw - d/dt(lambda Elam)`.
3. Si verifica off-shell e nel limite `k->0`.

Output:
1. `noetherId = 0`, `checkA=True`.
2. `checkB(k->0)=True`.
3. `propEq = 0`, `checkC=True`.
4. `idVals={0,0,0,0}`, `maxAbs=0`, `checkNum=True`, `check=True`.

Conclusione:
1. L'identita' Noether/Bianchi minisuperspace resta valida anche con curvatura spaziale.
2. La propagazione dei vincoli mantiene la stessa struttura del caso piatto.

## 23 - Covariant Causal-Domain Variation (4D)

Riferimento: `scripts/23_covariant_causal_domain_variation.wl`, log `logs/23_covariant_causal_domain_variation.log`.

Calcolo fatto:
1. Si implementa la forma 4D covariante dei funzionali leafwise del paper:
   `I0[theta] = Int4[mu chi]`, `I1[theta] = Int4[mu chi K2]`,
   con `mu = sqrt(-g) W sqrt(X) delta(Theta-theta)`.
2. Si espande al primo ordine:
   `mu -> mu + eps dmu`, `chi -> chi + eps dchi`, `K2 -> K2 + eps dK2`.
3. Si ottengono i contributi integrand-level:
   `dj0 = chi dmu + mu dchi`,
   `dj1 = chi K2 dmu + chi mu dK2 + mu K2 dchi`.
4. Si separa formalmente bulk/boundary sostituendo `dchi -> db`.
5. Si verifica la forma del termine di bordo nel canale `deltaQ`.

Output:
1. `dj0 = chi*dmu + dchi*mu`, `checkA=True`.
2. `dj1 = chi*dmu*k2 + chi*dk2*mu + dchi*k2*mu`, `checkB=True`.
3. `bulk0 = chi*dmu`, `boundary0 = mu`, `checkC=True`, `checkD=True`.
4. `bulk1 = chi*(dmu*k2 + dk2*mu)`, `boundary1 = k2*mu`, `checkE=True`, `checkF=True`.
5. `deltaQ = (-(I1*(dI0bulk + b0*xi)) + I0*(dI1bulk + b1*xi))/I0^2`.
6. `boundaryQ = (b1*I0 - b0*I1)/I0^2`, `checkG=True`.
7. `check=True`.

Conclusione:
1. La decomposizione covariante 4D dei funzionali con indicatore causale `chi_Theta` e' algebricamente consistente.
2. Il termine di bordo nel rapporto `Q` coincide con la struttura gia' usata nei notebook finite-domain, ma ora con significato covariante esplicito.

## 24 - Covariant Boundary Shape Derivative

Riferimento: `scripts/24_covariant_boundary_shape_derivative.wl`, log `logs/24_covariant_boundary_shape_derivative.log`.

Calcolo fatto:
1. Si rappresenta il dominio mobile tramite indicatore:
   `chiEps(r) = HeavisideTheta(R + eps rho - r)`.
2. Si deriva rispetto a `eps` al primo ordine:
   `dchi = rho DiracDelta(R-r)`.
3. Si verifica equivalenza tra:
   derivata distribuzionale (integrale con DiracDelta) e derivata di Leibniz dell'integrale a limite mobile.
4. Si fa controllo numerico indipendente con differenze finite.

Output:
1. `dchi = rho*DiracDelta[-r + R]`, `checkA=True`.
2. `dI0Delta = ((5 + R^2)*rho)/(5*E^R)`, `dI0Leib = ((5 + R^2)*rho)/(5*E^R)`, `checkB=True`.
3. `dI1Delta = ((9 + 2*R)*(5 + R^2)*rho)/(30*E^R)`, `dI1Leib = ((9 + 2*R)*(5 + R^2)*rho)/(30*E^R)`, `checkC=True`.
4. `err0 = 1.0874654...*10^-17`, `err1 = 1.2528183...*10^-17`.
5. `checkNum0=True`, `checkNum1=True`, `check=True`.

Conclusione:
1. Il contributo di bordo non e' postulato: emerge direttamente dalla variazione dell'indicatore causale.
2. L'accordo analitico e numerico a precisione macchina conferma la forma distribuzionale del termine di shape derivative.

## 25 - Covariant Nonlocal Action Variation

Riferimento: `scripts/25_covariant_nonlocal_action_variation.wl`, log `logs/25_covariant_nonlocal_action_variation.log`.

Calcolo fatto:
1. Si varia il termine di azione non-locale del paper su una foglia:
   `Sleaf = lambda (I1 - K0^2 I0)`,
   con `lambda -> lambda + eps eta`, `I0 -> I0 + eps dI0`, `I1 -> I1 + eps dI1`.
2. Si ottiene:
   `deltaSleaf = eta (I1-K0^2 I0) + lambda (dI1-K0^2 dI0)`.
3. Si inserisce la decomposizione covariante:
   `dI0 = A0 + B0`,
   `dI1 = A1 + C1 + B1`,
   con `A`/`C` bulk e `B` boundary.
4. Si separa:
   `deltaC = (A1 + C1 - K0^2 A0) + (B1 - K0^2 B0)`.
5. Si verifica il canale normalizzato `deltaQ` e la soppressione IR al crescere del volume `I0 ~ L^3`.

Output:
1. `deltaSleaf = eta*i1 - eta*i0*k0^2 + di1*lam - di0*k0^2*lam`, `checkA=True`.
2. `di0Cov = A0 + B0`, `di1Cov = A1 + B1 + C1`.
3. `deltaC = A1 + B1 + C1 - A0*k0^2 - B0*k0^2`.
4. `bulkPart = A1 + C1 - A0*k0^2`, `boundaryPart = B1 - B0*k0^2`, `checkB=True`.
5. `checkC(boundary neutral when B1=K0^2 B0)=True`.
6. `deltaQIR = (A1 + B1 + C1 - A0*q - B0*q)/i0`, `checkD=True`.
7. `vals(L={10,20,40}) = {1.555e-3, 1.94375e-4, 2.4296875e-5}`.
8. `ratios = {~8, ~8}`, `checkE=True`.
9. `check=True`.

Conclusione:
1. La variazione dell'azione non-locale e' chiusa in forma bulk+boundary coerente con `S = Sloc + integral dTheta lambda(Theta)[I1-K0^2 I0]`.
2. Nel canale normalizzato, il contributo variazionale scala come `1/I0` (qui `~L^-3`), cioe' resta IR-dominato e localmente soppresso.

## 26 - Covariant General EOM Channels

Riferimento: `scripts/26_covariant_general_eom_channels.wl`, log `logs/26_covariant_general_eom_channels.log`.

Calcolo fatto:
1. Si parte dall'azione completa `S = SLoc + integral dtheta lambda(theta) (I1-K0^2 I0)`.
2. Si certifica il canale vincolo da variazione in `lambda`:
   `deltaS_eta = eta (I1-K0^2 I0)`.
3. Si deriva in forma simbolica l'identita' di misura covariante per variazione metrica:
   `delta mu = -(1/2) mu h_{mn} delta g^{mn}`.
4. Si estraggono i coefficienti metrici bulk di `delta I0` e `delta I1`, introducendo il kernel variazionale di `K2`:
   `delta K2 = PiG delta g + PiTheta deltaTheta + divV`.
5. Si ottiene il tensore efficace non-locale bulk dalla definizione variazionale:
   `TNLbulk = lambda chi mu [ (K2-K0^2) h - 2 PiG ]`.
6. Si assembla il canale metrico covariante:
   `(G+Lambda g)/kappa - TTheta - TNL = 0`.
7. Si assembla il canale scalare:
   `alpha BoxTheta + JNL = 0`.
8. Si aggiunge check minisuperspace del ramo de Sitter:
   `lambda0(Lambda=0,H^2=K0^2/3)=1/(2 kappa)`.

Output:
1. `checkA(constraint channel)=True`.
2. `checkB(dmu metric identity)=True`.
3. `checkC(dI0 metric coeff)=True`, `checkD(dI1 metric coeff)=True`.
4. `checkE(dS_NL metric coeff)=True`, `checkF(TNL bulk explicit)=True`.
5. `checkG(scalar channel assembly)=True`.
6. `checkH(metric channel assembly)=True`.
7. `checkI(deSitter lambda)=True`.
8. `check=True`.

Conclusione:
1. Le EOM generali (vincolo/metrica/scalare) sono ora certificate in pipeline Mathematica, non solo riportate in forma testuale.
2. Il contributo non-locale metrico e' esplicitato come tensore efficace bulk con origine variazionale tracciabile.
3. Il ramo de Sitter con `Lambda=0` resta coerente nella chiusura EOM del proxy omogeneo.

## Quadro complessivo (00-26)

Risultato tecnico consolidato:
1. Tutti gli script `00..26` hanno check positivi nei log; nei notebook iniziali compaiono come `checkA/checkB/...`, mentre nei blocchi successivi anche come `check=True` unificato.
2. Le identita variazionali simboliche sono consistenti.
3. Le verifiche numeriche indipendenti danno residui a livello macchina o nulli.
4. La riduzione FRW produce EOM chiuse, propagazione coerente dei vincoli, riduzione lineare completa e identita' Noether/Bianchi off-shell, anche con estensione a curvatura spaziale.
5. E' ora esplicitata anche la catena 4D covariante del dominio causale (`chi_Theta`, `delta(Theta-theta)`) con variazione bulk+boundary, derivata di forma e canali EOM generali certificati.
6. Coerentemente con il testo del paper, la parte qui validata riguarda la variazione covariante dato un indicatore causale `chi_Theta`; la costruzione geometrica universale di `chi_Theta` da mutua raggiungibilita' causale resta un blocco separato.

Risultato fisico consolidato:
1. Il vincolo globale su `Q=I1/I0` e matematicamente ben posto in formulazione variazionale.
2. A dominio finito compaiono contributi di bordo controllabili e fisicamente interpretabili.
3. Il flusso IR distingue naturalmente comportamento tipo Minkowski (decadimento) e tipo de Sitter (costanza).
4. In cosmologia omogenea il vincolo seleziona `H^2=K0^2/3` e, al primo ordine, irrigidisce il canale perturbativo vincolato.
5. Nel formalismo covariante generale, i termini da variazione del dominio causale entrano come contributi distribuzionali di bordo e mantengono la soppressione IR nel canale normalizzato.

## Riproduzione rapida

Comando unico:

```powershell
cd c:\lean\StaticSanityCheck\Paper\Variational_Derivation
.\run_all.ps1
```

Artefatti:
1. Notebook: `Paper/Variational_Derivation/*.nb`.
2. Log: `Paper/Variational_Derivation/logs/*.log`.

## Nota EFT / scale di validita
1. La catena `00..26` valida la parte variazionale e mostra soppressione IR controllata, ma non fissa una scala UV fondamentale della teoria.
2. La discussione quantitativa di cutoff operativo, strong-coupling proxy e regime di validita e' implementata nel blocco perturbativo (`Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`).
3. Quindi questa cartella va letta come fondazione strutturale EFT (operatori e variazioni), non come determinazione completa del dominio UV.
