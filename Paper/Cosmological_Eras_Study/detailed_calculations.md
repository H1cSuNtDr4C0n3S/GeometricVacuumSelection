# Cosmological Eras Study - Detailed Calculations and Results

Questo documento raccoglie i calcoli eseguiti in `Paper/Cosmological_Eras_Study/scripts/`
e i risultati registrati in `Paper/Cosmological_Eras_Study/logs/`.

Convenzioni:
1. `a(t)`: fattore di scala FRW.
2. `H(t)`: parametro di Hubble.
3. `w`: parametro barotropico dell'era considerata.
4. `D_Theta`: dominio causale/finito sulla foglia `Theta=const`.
5. `Q = I1/I0`: canale normalizzato del vincolo non-locale.

## Notebook 00 - Friedmann Eras Data Baseline Consistency Test

Script: `Paper/Cosmological_Eras_Study/scripts/00_friedmann_eras_data_baseline_consistency_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/00_friedmann_eras_data_baseline_consistency_test.log`

### Obiettivo

Costruire un baseline FRW riproducibile con parametri Planck18 + DESI DR1 e verificare:
1. consistenza del parsing/best-fit;
2. consistenza interna della ricostruzione di `H0`;
3. ordine corretto delle transizioni tra ere;
4. comportamento atteso di `w_eff` nelle tre regioni (radiazione, materia, de Sitter tardo);
5. vicinanza Planck-DESI sul background `E(z)`.

### Setup verificato dal log

1. `sourceUsed="embedded_snapshot"` (kernel cloud), con hash verificato:
   - `expectedHash=948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25`
   - `actualHash=948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25`
2. Parametri DESI letti:
   - `hrdrag=101.94944`
   - `rdrag=151.35012`
   - `H0rdrag=10194.944`
   - `omegam=0.29353038`
   - `omegal=0.70639`
3. Ricostruzione:
   - `h0DESIFromHrdrag=67.35999945028124`
   - `h0DESIFromH0rdrag=67.35999945028124`

### Derivati cosmologici

1. Densita di radiazione:
   - `orPlanck=0.0000912427019293302`
   - `orDESI=0.0000920572451620294`
2. Densita di vuoto:
   - `olPlanck=0.6902487572980707`
   - `olDESI=0.70639`
3. Redshift di uguaglianza:
   - `zEqPlanck=3392.805679273281`
   - `zEqDESI=3187.563588703517`
4. Transizione materia-vuoto:
   - `zLambdaMatterPlanck=0.306291721795789`
   - `zLambdaMatterDESI=0.34007928913383556`

### Sonde fisiche di era

1. `wRadPlanck(z=1e6)=0.33220589221708013` (vicino a `1/3`)
2. `wMatPlanck(z=100)=0.009631223177569712` (vicino a `0`)
3. `wDSPlanck(z=-0.9)=-0.9995515626863059` (vicino a `-1`)

### Confronto background Planck vs DESI

Campioni `zSamples={0,1,10,100}`:
1. `eRatios={1.0000062186032455,0.9820279811418215,0.9737732857167074,0.9745128732102355}`
2. `maxERelDiff=0.026226714283292596`

### Check e classificazione

1. `checkA` ... `checkI` = `True`
2. `classification="STANDARD_ERA_HISTORY_COMPATIBLE"`
3. `check=True`

Interpretazione sintetica: la baseline congelata e consistente, e la sequenza standard delle ere cosmologiche resta preservata; questo abilita i prossimi test di dominio finito/non-localita era-per-era.

## Notebook 01 - Radiation Era Finite-Domain IR Scaling Test

Script: `Paper/Cosmological_Eras_Study/scripts/01_radiation_era_finite_domain_ir_scaling_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/01_radiation_era_finite_domain_ir_scaling_test.log`

### Obiettivo

Validare nel regime di radiazione che:
1. il dominio finito causale resti ben definito;
2. il contributo non-locale resti subleading;
3. lo scaling al variare di `H` segua la legge attesa (`deltaQ ~ h^3`, `I0 ~ h^-3`).

### Parametri principali dal log

1. Pivot era:
   - `zRef=1000000`
   - `zHalf=707106.4882933287...`
   - `hNormHalf=0.5003501284609077`
2. Sonde `w_eff`:
   - `wRef=0.33220589221708013`
   - `wHalf=0.3317411215114061`
3. Proxy e regularization:
   - `z2=1/20`
   - `alpha=1/100000`
   - `rs=1/1000`
   - `aReg=0.0010000000000000000005...`

### Quantita IR

1. `i0Ref=3.8785912845623018...`
2. `i0Half=30.96363721205441`
3. `locRef=9.856711585922334...*10^-8`
4. `locHalf=9.86315347597068*10^-8`
5. `deltaQRef=2.54131226075672...*10^-8`
6. `deltaQHalf=3.185398833593922*10^-9`
7. `ratioDelta=0.12534464507897247`
8. `ratioI0=7.983217343703347`
9. `ratioLoc=1.0006535536717485`
10. `epsRef=8.47104086918...*10^-9`
11. `epsHalf=4.24125642671...*10^-9`

### Esito

1. Tutti i check `D...L` risultano `True`
2. `classification="RADIATION_ERA_IR_SUPPRESSED"`
3. `check=True`

Lettura fisica: nel regime dominato da radiazione il canale non-locale rimane una correzione piccola e scala correttamente con il dominio IR.

## Notebook 02 - Matter Era Finite-Domain IR Scaling Test

Script: `Paper/Cosmological_Eras_Study/scripts/02_matter_era_finite_domain_ir_scaling_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/02_matter_era_finite_domain_ir_scaling_test.log`

### Obiettivo

Ripetere il test IR nel regime materia e verificare la stessa classe di comportamento del caso radiazione.

### Parametri principali dal log

1. Pivot era:
   - `zRef=100`
   - `zHalf=62.62601301969109...`
   - `hNormHalf=0.49732086530589376`
2. Sonde `w_eff`:
   - `wRef=0.009631223177569712`
   - `wHalf=0.006125680279393943`
3. Proxy e regularization:
   - `z2=1/20`
   - `alpha=1/100000`
   - `rs=1/1000`
   - `aReg=0.0010000000000000000005...`

### Quantita IR

1. `i0Ref=3.8785912845623018...`
2. `i0Half=31.53290448964249`
3. `locRef=9.856711585922334...*10^-8`
4. `locHalf=9.863192531707576*10^-8`
5. `deltaQRef=2.54131226075672...*10^-8`
6. `deltaQHalf=3.1279048240406837*10^-9`
7. `ratioDelta=0.12308227022480485`
8. `ratioI0=8.129989002747159`
9. `ratioLoc=1.0006575160213167`
10. `epsRef=8.47104086918...*10^-9`
11. `epsHalf=4.21559532194...*10^-9`

### Esito

1. Tutti i check `D...L` risultano `True`
2. `classification="MATTER_ERA_IR_SUPPRESSED"`
3. `check=True`

Lettura fisica: anche in era materia il termine non-locale resta subleading e con scaling coerente; stessa classe IR del caso radiazione.

## Notebook 03 - Era Sequence Continuity and Global IR Channel Test

Script: `Paper/Cosmological_Eras_Study/scripts/03_era_sequence_continuity_global_ir_channel_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/03_era_sequence_continuity_global_ir_channel_test.log`

### Obiettivo

Verificare con un singolo test globale:
1. ordinamento corretto delle ere;
2. continuita di `w_eff` attorno a `zEq` e `zLambdaMatter`;
3. carattere subleading e monotono del canale IR lungo la sequenza cosmologica.

### Valori chiave

1. `zEq=3392.805679273281`
2. `zLambdaMatter=0.306291721795789`
3. `zProbe={1000000,3392.805679273281,100,0.306291721795789,0,-0.9}`
4. `eProbe={9.568318965715258*10^9,155592.26776155384,573.1827370567321,1.1750588047327108,1.,0.8309984515162112}`
5. `wProbe={0.33220589221708013,0.1666666666334025,0.009631223177569712,-0.49983965390495094,-0.6902183430640942,-0.9995515626863059}`
6. `deltaProbe={2.2291159593946854*10^-4,9.584947981371325*10^-19,4.7918700731453395*10^-26,4.1286163794239125*10^-34,2.5446363581470136*10^-34,1.4602471648707273*10^-34}`
7. `epsProbe={8.1159641088355545*10^-25,1.319752471975318*10^-29,4.8618054419226003*10^-32,9.966991191612092*10^-35,8.482121193823378*10^-35,7.048629577640062*10^-35}`

### Esito

1. Check `D...L` tutti `True`
2. `classification="ERAS_PIPELINE_CONTINUITY_CONFIRMED"`
3. `check=True`

Lettura fisica: la sequenza radiazione->materia->de Sitter e continua e il canale non-locale resta una correzione IR globale, positiva e fortemente subleading.

## Notebook 04 - Era Proxy Equivalence IR Class Test

Script: `Paper/Cosmological_Eras_Study/scripts/04_era_proxy_equivalence_ir_class_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/04_era_proxy_equivalence_ir_class_test.log`

### Obiettivo

Confrontare due proxy quasi-locali per `D_Theta`:
1. `redshift_threshold`: `r <= sqrt(1-z2)/H(z)`
2. `apparent_horizon`: `r <= 1/H(z)`

e verificare che producano la stessa classe IR nelle ere:
1. radiazione
2. materia
3. de Sitter tardo

### Setup dal log

1. `sourceUsed="embedded_snapshot"`
2. `z2=1/20`
3. `proxyKappaRed=0.974679434480896...`
4. `proxyKappaApp=1`

### Risultati per era

Per ciascuna era, entrambi i proxy danno `classification="IR_SUPPRESSED_SUBLEADING"`.

1. Radiazione:
   - `ratioDelta(red)=0.12526278027351173`
   - `ratioDelta(app)=0.12526278027351173`
2. Materia:
   - `ratioDelta(red)=0.12300139639329095`
   - `ratioDelta(app)=0.12300139639329098`
3. De Sitter tardo:
   - `ratioDelta(red)=0.5738529830384365`
   - `ratioDelta(app)=0.5738529830384366`

### Robustezza proxy

1. Classi coincidenti:
   - `redClasses={"IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING"}`
   - `appClasses={"IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING"}`
2. Domini realmente diversi:
   - `domainRefRatios={1.0799772127212148,1.0799772127212148,1.0799772127212148}`
3. Check di equivalenza:
   - `checkA=True`
   - `checkB=True`
   - `checkC=True`
   - `checkD=True`
   - `checkE=True`
   - `checkF=True`
   - `checkG=True`

### Esito finale

1. `classification="ERA_PROXY_CLASS_EQUIVALENT_IR"`
2. `check=True`

Lettura fisica: la classe IR resta invariata cambiando proxy quasi-locali ragionevoli; il comportamento non e un artefatto di tuning della frontiera di dominio.

## Notebook 05 - Pipeline 00-04 Appendix-Ready Summary Test

Script: `Paper/Cosmological_Eras_Study/scripts/05_pipeline_00_04_appendix_ready_summary_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/05_pipeline_00_04_appendix_ready_summary_test.log`

### Obiettivo

Produrre una sintesi unica, riproducibile e citabile dei risultati `00-04`:
1. classi osservate;
2. check aggregato;
3. tabella markdown pronta per appendice.

### Output principali dal log

1. Righe aggregate:
   - `00 -> STANDARD_ERA_HISTORY_COMPATIBLE`
   - `01 -> RADIATION_ERA_IR_SUPPRESSED`
   - `02 -> MATTER_ERA_IR_SUPPRESSED`
   - `03 -> ERAS_PIPELINE_CONTINUITY_CONFIRMED`
   - `04 -> ERA_PROXY_CLASS_EQUIVALENT_IR`
2. `classesObserved=classesExpected`
3. `checkAllRows=True`
4. `checkClasses=True`
5. `summaryMarkdown` generato automaticamente (tabella completa)

### Esito

1. `classification="PIPELINE_00_04_APPENDIX_READY"`
2. `check=True`

Lettura fisica/editoriale: la pipeline ere cosmologiche e consolidata in forma compatta e direttamente traslabile nel materiale supplementare.

## Notebook 06 - Complete EOM With Matter Selfconsistent Branch Search

Script: `Paper/Cosmological_Eras_Study/scripts/06_complete_eom_with_matter_selfconsistent_branch_search.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/06_complete_eom_with_matter_selfconsistent_branch_search.log`

### Obiettivo

Cercare una soluzione cosmologica auto-consistente delle EOM complete con materia nel settore minisuperspace FRW, verificando simultaneamente i canali:
1. `Ea=0`
2. `EN=0`
3. `Elam=0`

### Equazioni e riduzione

1. Lagrangiano usato:
   - `L = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3) - rho0 a^(-3 w) N`
2. Canale vincolo in gauge `N=1`:
   - `constraintForm = -K0^2 + 3 (a'/a)^2`
3. Sul ramo hard (`H0^2=K0^2/3`) si trova:
   - `lamRule` da `EN=0`
   - `EaCons=0`, `ENCons=0`, `ElamCons=0`

### Test materia e ricerca ramo alternativo

1. Continuita materia:
   - `rhoConservation=0`
2. Dinamica di `lambda`:
   - `lamPrimeDust = 5/(2*Sqrt[3]*E^((3*Sqrt[3])/5))` (dust, non nullo)
   - `lamPrimeDE = 0` per `w=-1`
3. Test no-go power-law sotto vincolo hard:
   - `ElamPower = t^(-2+3 p) (3 p^2 - K0^2 t^2)`
   - `powerNoGoReduce=False`
4. Residui numerici canali:
   - `sampleResiduals={{0,0,0},{0,0,0},{0,0,0},{0,0,0}}`
   - `maxResidual=0`

### Esito

1. `checkA..checkI=True`
2. `classification="COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY"`
3. `check=True`

Lettura fisica: nel setup hard il sistema completo ammette un ramo auto-consistente de Sitter con `lambda(t)` che assorbe la materia; non emerge un ramo power-law standard soddisfacente il canale di vincolo.

## Notebook 07 - Complete EOM Soft-Closure Evolving Branch Test

Script: `Paper/Cosmological_Eras_Study/scripts/07_complete_eom_soft_closure_evolving_branch_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/07_complete_eom_soft_closure_evolving_branch_test.log`

### Obiettivo

Verificare se una chiusura soft (ispirata alla formulazione finite-domain) del canale di vincolo consente una soluzione cosmologica evolutiva auto-consistente nelle EOM complete.

### Chiusura soft testata

1. Vincolo normalizzato:
   - `3 H^2 - K0^2 = sigma/a^3`
2. Equivalente:
   - `H2soft(a) = K0^2/3 + sigma/(3 a^3)`

### Coerenza canali EOM completi

Con sostituzioni:
1. `N=1`, `w=0`, `Lambda=K0^2`, `rho0=sigma/kappa`
2. `lambda=0`
3. `a_dot = a H`, `a_ddot = a (H^2 + H_dot)`, `H_dot = -sigma/(2 a^3)`

si ottiene:
1. `EaSoft=0`
2. `ENSoft=0`
3. `ElamSoft=sigma`
4. `softConstraintResidual=0`

### Transizione materia -> de Sitter

1. `weff(a) = -1 + sigma/(sigma + K0^2 a^3)`
2. Valori campione (`K0=3/5`, `sigma=1/3`):
   - `weffEarly(a=1e-2) = -1.07999883360...*10^-6` (vicino a `0`, materia-like)
   - `weffMid(a=1) = -0.51923076923...`
   - `weffLate(a=1e2) = -0.999999074074...` (vicino a `-1`, de Sitter-like)

### Check numerici

1. `eaNumVals={0,0,0,0,0}`
2. `enNumVals={0,0,0,0,0}`
3. `elamNumVals={0,0,0,0,0}`
4. `maxResidual=0`

### Confronto con ramo hard

1. `ElamPowerHard = t^(-2+3 p) (3 p^2 - K0^2 t^2)`
2. `hardPowerNoGo=False`
3. Il no-go power-law del canale hard resta confermato.

### Esito

1. `checkA..checkH=True`
2. `classification="COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND"`
3. `check=True`

Lettura fisica: la chiusura soft del canale di vincolo apre un ramo cosmologico evolutivo auto-consistente nelle EOM complete, compatibile con una dinamica che interpola tra regime materia-like e de Sitter tardivo.

## Notebook 08 - Paper-Ready Hard/Soft EOM Bridge Test

Script: `Paper/Cosmological_Eras_Study/scripts/08_paper_ready_hard_soft_eom_bridge_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/08_paper_ready_hard_soft_eom_bridge_test.log`

### Obiettivo

Costruire una sintesi paper-ready, algebraicamente chiusa, del ponte tra:
1. ramo hard (vincolo rigido da `Elam=0`);
2. ramo soft (chiusura finite-domain efficace);
3. forma LCDM-like emergente nel ramo soft.

### Risultati simbolici principali

1. Formula generale da `EN`:
   - `lamGeneric = -1/6*(-3*h2 + Lambda + kappa*rh)/(h2*kappa)`
2. Ramo hard:
   - `3 H^2 - K0^2 = 0`
   - `lamHardA = -1/2*(-K0^2 + Lambda + (kappa*rho0)/a^(3*(1 + ww)))/(K0^2*kappa)`
   - equivalente: `lambda_hard(a) = (K0^2 - Lambda - kappa rho(a)) / (2 kappa K0^2)`
3. Ramo soft:
   - `h2Soft(a) = K0^2/3 + sigma/(3*a^3)`
   - `3 H^2 - K0^2 = sigma/a^3`
4. Ricostruzione LCDM-like:
   - `h2SoftRebuild(a) = ((K0^2 + sigma)*(K0^2/(K0^2 + sigma) + sigma/(a^3*(K0^2 + sigma))))/3`
   - `omegaM0 = sigma/(K0^2 + sigma)`
   - `omegaL0 = K0^2/(K0^2 + sigma)`
5. Equazione di stato efficace:
   - `weffSoft(a) = -((a^3*K0^2)/(a^3*K0^2 + sigma))`
   - forma equivalente: `weff(a) = -1 + sigma/(sigma + K0^2 a^3)`
   - limiti: `weffEarlyLim=0`, `weffLateLim=-1`
6. Limite di raccordo:
   - `softToHard = K0^2/3` (per `sigma -> 0`).

### Anchor numeriche dal log

Con griglia `aGrid={1/100, 1, 100}`:
1. `h2SoftVals={111111.231111..., 0.231111..., 0.120000111111...}`
2. `weffVals={-1.07999883360...*10^-6, -0.51923076923..., -0.999999074074...}`
3. `lamHardVals={-462962.832244..., -0.332244008714..., 0.130718491285...}`

### Check formali

1. `checkA(lam hard from EN) = True`
2. `checkB(lam hard vacuum limit) = True`
3. `checkC(soft H2 LCDM form) = True`
4. `checkD(omegas sum to one) = True`
5. `checkE(weff closed form) = True`
6. `checkF(weff limits early/late) = True`
7. `checkG(soft->hard at sigma=0) = True`
8. `checkH(weff numeric anchors) = True`
9. `checkI(lamHard monotonic sample) = True`

### Esito

1. `classification="PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED"`
2. `check=True`

Lettura fisica/editoriale: il ramo hard e il ramo soft vengono ricondotti a un quadro unico e presentabile nel main text; il limite hard resta de Sitter, mentre la chiusura soft conserva la consistenza EOM e produce un comportamento materia-like a piccoli `a` con asintotica de Sitter a grandi `a`.

## Notebook 09 - Complete EOM Soft Branch Parameter Robustness Test

Script: `Paper/Cosmological_Eras_Study/scripts/09_complete_eom_soft_branch_parameter_robustness_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/09_complete_eom_soft_branch_parameter_robustness_test.log`

### Obiettivo

Verificare che il ramo soft delle EOM complete:
1. resti auto-consistente nei canali dinamici completi;
2. non dipenda da tuning stretto dei parametri;
3. riproduca la scala di transizione materia-vuoto in accordo con il riferimento Planck.

### Coerenza canali completi (ramo soft)

Con la chiusura:
1. `3 H^2 - K0^2 = sigma/a^3`
2. `H2soft(a)=K0^2/3 + sigma/(3 a^3)`

si ottiene dal log:
1. `EaSoft = 0`
2. `ENSoft = 0`
3. `ElamSoft = sigma`
4. `softConstraintResidual = 0`

### Dinamica efficace e scala di transizione

1. `weff(a) = -1 + sigma/(sigma + K0^2 a^3)`
2. `aTransition = (sigma/K0^2)^(1/3)`
3. `weffAtTransition = -1/2`

Questo identifica in modo univoco il passaggio materia-like -> de Sitter-like del ramo.

### Anchor Planck

Dal log:
1. `hPlanck=0.6766`
2. `omPlanck=0.30966`
3. `orPlanck=0.0000912427019293302`
4. `olPlanck=0.6902487572980707`
5. `rPlanck=omPlanck/olPlanck=0.4486208728750804`
6. `zLambdaMatterPlanck=0.306291721795789`
7. `zSoftPlanck=0.306291721795789`

La coincidenza `zSoftPlanck = zLambdaMatterPlanck` verifica l'ancoraggio del ramo soft alla transizione standard.

### Scan di robustezza parametrica

Range testato:
1. `scanRatios = {1/20, 1/10, 3/10, 9/20, 7/10, 1, 3}`
2. `ratioSpan = 60`

Per tutti i campioni:
1. `pass=True`
2. `maxResidual=0` per campione
3. `maxResidualScan=0`
4. `zTransitions` copre sia passato sia futuro:
   - massimo `1.71441761659...`
   - minimo `-0.306638725649...`

Esempi di `weff` dal log:
1. `ratio=0.05`: `weffEarly=-1.9999600008...*10^-5`, `weffMid=-0.95238...`, `weffLate=-0.99999995...`
2. `ratio=0.45`: `weffEarly=-2.22221728396...*10^-6`, `weffMid=-0.689655...`, `weffLate=-0.99999955...`
3. `ratio=3`: `weffEarly=-3.33333222222...*10^-7`, `weffMid=-0.25`, `weffLate=-0.999997...`

### Check finali

1. `checkA(complete channels on soft branch)=True`
2. `checkB(weff closed form)=True`
3. `checkC(transition at weff=-1/2)=True`
4. `checkD(Planck transition anchor)=True`
5. `checkE(scan passes all samples)=True`
6. `checkF(scan residuals)=True`
7. `checkG(broad ratio span, past/future transitions)=True`

### Esito

1. `classification="COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"`
2. `check=True`

Lettura fisica: il ramo soft non e una soluzione fragile; resta coerente nelle EOM complete su un intervallo parametrico ampio e mantiene una scala di transizione compatibile con il benchmark cosmologico standard.

## Notebook 10 - Pipeline 06-09 Appendix-Ready Summary Test

Script: `Paper/Cosmological_Eras_Study/scripts/10_pipeline_06_09_appendix_ready_summary_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/10_pipeline_06_09_appendix_ready_summary_test.log`

### Obiettivo

Consolidare in un singolo test il blocco EOM `06-09`:
1. ramo hard completo;
2. ramo soft evolutivo;
3. ponte hard/soft paper-ready;
4. robustezza parametrica anti-fine-tuning.

### Check del blocco hard (06)

Dal log:
1. `constraintForm(06) = -K0^2 + (3*Derivative[1][aa][t]^2)/aa[t]^2`
2. `check06A=True`
3. `check06B=True`
4. `check06C=True`
5. `check06D=True`
6. `check06E=True`
7. `check06F=True`
8. `lamPrimeDust06 = 5/(2*Sqrt[3]*E^((3*Sqrt[3])/5))`
9. `hardPowerNoGo06 = False`

### Anchor del blocco soft e del bridge (07-08)

1. `weffEarly07=-1.07999883360...*10^-6`
2. `weffMid07=-0.51923076923...`
3. `weffLate07=-0.999999074074...`
4. `lamGeneric08 = -1/6*(-3*h2 + Lambda + kappa*rh)/(h2*kappa)`
5. `omegaM0 = sigma/(K0^2 + sigma)`
6. `omegaL0 = K0^2/(K0^2 + sigma)`

### Anchor della robustezza (09)

1. `zLambdaMatterPlanck = 0.306291721795789`
2. `zSoftPlanck = 0.306291721795789`
3. `scanRatios09 = {1/20, 1/10, 3/10, 9/20, 7/10, 1, 3}`
4. `ratioSpan09 = 60`
5. `zTransitions09` copre:
   - massimo `1.71441761659...`
   - minimo `-0.306638725649...`
6. `maxResidualScan09 = 0`

### Tabella consolidata 06-09

Righe aggregate:
1. `06 -> COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY`
2. `07 -> COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND`
3. `08 -> PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED`
4. `09 -> COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED`

Controlli finali:
1. `checkAllRows = True`
2. `checkClasses = True`
3. `summaryMarkdown` generato nel log con tabella pronta da appendice.

### Esito

1. `classification="PIPELINE_06_09_APPENDIX_READY"`
2. `check=True`

Lettura fisica/editoriale: la catena hard -> soft -> robustezza e chiusa e internamente coerente; il blocco EOM e pronto per essere citato nel testo principale/appendice come evidenza unica non frammentata.

## Notebook 11 - Section Closure Referee-Ready Test

Script: `Paper/Cosmological_Eras_Study/scripts/11_section_closure_referee_ready_test.wl`  
Log: `Paper/Cosmological_Eras_Study/logs/11_section_closure_referee_ready_test.log`

### Obiettivo

Chiudere la sezione con un test unico e sintetico che validi, in forma referee-facing:
1. claim hard;
2. claim soft;
3. claim di robustezza (non fine-tuning) con ancoraggio cosmologico standard.

### Claim hard (vincolo rigido)

Dal log:
1. `constraintForm = -K0^2 + (3*Derivative[1][aa][t]^2)/aa[t]^2`
2. `hardPowerNoGo = False`
3. `checkA=True`
4. `checkB=True`
5. `checkC=True`
6. `checkD=True`
7. `checkE=True`

Interpretazione: sotto `Elam=0` il ramo auto-consistente resta de Sitter e il ramo power-law non passa.

### Claim soft (chiusura finite-domain)

Dal log:
1. `weffEarly = -1.07999883360...*10^-6`
2. `weffMid = -0.51923076923...`
3. `weffLate = -0.999999074074...`
4. `checkF=True`
5. `checkG=True`
6. `checkH=True`
7. `checkI=True`

Interpretazione: nei canali completi il ramo soft e coerente e realizza la transizione materia-like -> de Sitter-like.

### Ancoraggio standard + anti-fine-tuning

Dal log:
1. `zLambdaMatterPlanck = 0.306291721795789`
2. `zSoftPlanck = 0.306291721795789`
3. `checkJ=True`
4. `scanRatios = {1/20, 1/5, 1, 3}`
5. `ratioSpan = 60`
6. `scanData = { ... "pass" -> True ... }` per tutti i punti
7. `maxResidualScan = 0`
8. `checkK=True`

Interpretazione: il ramo soft non e fragile, resta coerente su uno spazio parametrico ampio e mantiene la scala di transizione corretta.

### Output editoriali pronti

Il log include direttamente:
1. `refereeClosureSnippet` (testo corto per main text/response to referee)
2. `closureChecklistMarkdown` (tabella compatta claim -> status)

### Esito

1. `classification="COSMOLOGICAL_ERAS_SECTION_CLOSURE_READY"`
2. `check=True`

Lettura finale: la sezione cosmologica e formalmente chiudibile; i risultati principali sono consolidati in un pacchetto sintetico, verificato e citabile.
