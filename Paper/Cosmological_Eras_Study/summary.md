# Cosmological Eras Study Summary

Mappa notebook `Paper/Cosmological_Eras_Study` -> passaggi matematici, output verificati e lettura fisica.

## Notebook

1. `00_Friedmann_Eras_Data_Baseline_Consistency_Test.nb`
`Riferimento`: fissare baseline cosmologica riproducibile (Planck18 + DESI DR1) e verificare che la sequenza radiazione -> materia -> de Sitter tardo resti standard.
`Output`:
1. sorgente dati usata nel run cloud: `sourceUsed="embedded_snapshot"` con hash coerente:
   `expectedHash=948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25`,
   `actualHash=948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25`;
2. ricostruzione `H0` consistente dal best-fit DESI:
   `h0DESIFromHrdrag=67.35999945028124`,
   `h0DESIFromH0rdrag=67.35999945028124`;
3. transizioni tra ere:
   `zEqPlanck=3392.805679273281`, `zEqDESI=3187.563588703517`,
   `zLambdaMatterPlanck=0.306291721795789`, `zLambdaMatterDESI=0.34007928913383556`;
4. sonde `w_eff`:
   `wRadPlanck(z=1e6)=0.33220589221708013`,
   `wMatPlanck(z=100)=0.009631223177569712`,
   `wDSPlanck(z=-0.9)=-0.9995515626863059`;
5. vicinanza background Planck-DESI:
   `eRatios={1.0000062186032455,0.9820279811418215,0.9737732857167074,0.9745128732102355}`,
   `maxERelDiff=0.026226714283292596`;
6. classificazione finale:
   `classification="STANDARD_ERA_HISTORY_COMPATIBLE"`.
`Perche questo passaggio`: stabilisce una baseline FRW congelata e verificata su cui costruire i test delle diverse ere senza dipendere da setup esterni del kernel.
`Interpretazione`: la storia cosmologica standard resta intatta (ordine ere corretto, `w_eff` coerente) e Planck/DESI restano nella stessa classe di comportamento di background.
`Esito log`: `check=True`.

2. `01_Radiation_Era_Finite_Domain_IR_Scaling_Test.nb`
`Riferimento`: testare in era di radiazione che il canale non-locale su dominio finito resti IR-soppresso e segua lo scaling atteso.
`Output`:
1. pivot era:
   `zRef=1000000`, `zHalf=707106.4882933287...`,
   `hNormHalf=0.5003501284609077`;
2. sonde radiazione:
   `wRef=0.33220589221708013`,
   `wHalf=0.3317411215114061`;
3. scaling canale:
   `ratioDelta=0.12534464507897247`,
   `ratioI0=7.983217343703347`,
   `ratioLoc=1.0006535536717485`;
4. subleading:
   `epsRef=8.47104086918...*10^-9`,
   `epsHalf=4.24125642671...*10^-9`;
5. classificazione:
   `classification="RADIATION_ERA_IR_SUPPRESSED"`.
`Perche questo passaggio`: verifica che il ramo hard-constraint non distrugga la fisica di era standard quando il dominio e reso operativo/finito.
`Interpretazione`: anche in regime radiazione il termine non-locale resta correzione piccola e con scaling coerente (`deltaQ ~ h^3`).
`Esito log`: `check=True`.

3. `02_Matter_Era_Finite_Domain_IR_Scaling_Test.nb`
`Riferimento`: ripetere il test precedente in era di materia per controllare robustezza della classe IR.
`Output`:
1. pivot era:
   `zRef=100`,
   `zHalf=62.62601301969109...`,
   `hNormHalf=0.49732086530589376`;
2. sonde materia:
   `wRef=0.009631223177569712`,
   `wHalf=0.006125680279393943`;
3. scaling canale:
   `ratioDelta=0.12308227022480485`,
   `ratioI0=8.129989002747159`,
   `ratioLoc=1.0006575160213167`;
4. subleading:
   `epsRef=8.47104086918...*10^-9`,
   `epsHalf=4.21559532194...*10^-9`;
5. classificazione:
   `classification="MATTER_ERA_IR_SUPPRESSED"`.
`Perche questo passaggio`: mostra che il comportamento IR non dipende da un singolo regime cosmologico.
`Interpretazione`: in era materia il meccanismo resta nella stessa classe qualitativa del caso radiazione (nessun fine-tuning di era).
`Esito log`: `check=True`.

4. `03_Era_Sequence_Continuity_Global_IR_Channel_Test.nb`
`Riferimento`: test globale di continuita tra ere con controlli su `zEq`, `zLambdaMatter` e sul canale IR lungo la sequenza cosmologica.
`Output`:
1. transizioni:
   `zEq=3392.805679273281`,
   `zLambdaMatter=0.306291721795789`;
2. ancore fisiche:
   `wProbe={0.33220589221708013, 0.1666666666334025, 0.009631223177569712, -0.49983965390495094, -0.6902183430640942, -0.9995515626863059}`;
3. canale IR:
   `deltaProbe={2.2291159593946854*10^-4, 9.584947981371325*10^-19, 4.7918700731453395*10^-26, 4.1286163794239125*10^-34, 2.5446363581470136*10^-34, 1.4602471648707273*10^-34}`;
   `epsProbe={8.1159641088355545*10^-25, 1.319752471975318*10^-29, 4.8618054419226003*10^-32, 9.966991191612092*10^-35, 8.482121193823378*10^-35, 7.048629577640062*10^-35}`;
4. classificazione:
   `classification="ERAS_PIPELINE_CONTINUITY_CONFIRMED"`.
`Perche questo passaggio`: chiude il requisito di pipeline: non solo test locali per era, ma continuita globale tra regimi.
`Interpretazione`: la storia standard e liscia tra le ere e il canale non-locale rimane un vincolo IR subleading su tutta la sequenza.
`Esito log`: `check=True`.

5. `04_Era_Proxy_Equivalence_IR_Class_Test.nb`
`Riferimento`: confronto diretto tra due proxy quasi-locali ragionevoli per il dominio finito (`redshift-threshold` vs `apparent-horizon`) nelle tre ere.
`Output`:
1. classi per era/proxy:
   `redClasses={"IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING"}`,
   `appClasses={"IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING","IR_SUPPRESSED_SUBLEADING"}`;
2. rapporti di scaling per era:
   `redRatios={0.12526278027351173,0.12300139639329095,0.5738529830384365}`,
   `appRatios={0.12526278027351173,0.12300139639329098,0.5738529830384366}`;
3. proxy effettivamente distinti:
   `domainRefRatios={1.0799772127212148,1.0799772127212148,1.0799772127212148}`;
4. check globali:
   `checkA..checkG=True`;
5. classificazione finale:
   `classification="ERA_PROXY_CLASS_EQUIVALENT_IR"`.
`Perche questo passaggio`: risponde esplicitamente al punto "non fine-tuning": il risultato IR non dipende da una scelta unica di proxy.
`Interpretazione`: con domini diversi ma causalmente ragionevoli, la pipeline conserva la stessa classe IR in radiazione, materia e de Sitter tardo.
`Esito log`: `check=True`.

6. `05_Pipeline_00_04_Appendix_Ready_Summary_Test.nb`
`Riferimento`: creare una tabella finale unica (test `00-04`) pronta da inserire in appendice con classi e check consolidati.
`Output`:
1. righe aggregate:
   `rows={("00","STANDARD_ERA_HISTORY_COMPATIBLE"),("01","RADIATION_ERA_IR_SUPPRESSED"),("02","MATTER_ERA_IR_SUPPRESSED"),("03","ERAS_PIPELINE_CONTINUITY_CONFIRMED"),("04","ERA_PROXY_CLASS_EQUIVALENT_IR")}`;
2. confronto classi:
   `classesObserved=classesExpected`;
3. tabella markdown pronta:
   `summaryMarkdown="| Test | Title | Classification | Check | ..."` con tutte le righe `True`;
4. classificazione finale:
   `classification="PIPELINE_00_04_APPENDIX_READY"`.
`Perche questo passaggio`: trasforma i risultati sparsi in una sintesi verificabile e direttamente citabile nel testo/appendice.
`Interpretazione`: la catena `00-04` e internamente coerente in blocco unico, non solo notebook-per-notebook.
`Esito log`: `check=True`.

7. `06_Complete_EOM_With_Matter_Selfconsistent_Branch_Search.nb`
`Riferimento`: avviare la ricerca di una soluzione cosmologica auto-consistente delle EOM complete con materia nel minisuperspace FRW.
`Output`:
1. canali completi derivati:
   `Ea`, `EN`, `Elam` (con materia barotropica);
2. vincolo hard dal canale `Elam`:
   `constraintForm = -K0^2 + 3 (a'/a)^2`;
3. ramo auto-consistente trovato:
   de Sitter (`H0^2=K0^2/3`) + `lambda(t)` da `EN=0`, con `EaCons=0`, `ENCons=0`, `ElamCons=0`;
4. test materia:
   `rhoConservation=0`, `lamPrimeDust=5/(2*Sqrt[3]*E^((3*Sqrt[3])/5))` (dinamica), `lamPrimeDE=0` per `w=-1`;
5. test no-go ramo power-law sotto hard-constraint:
   `ElamPower=t^(-2+3 p) (3 p^2-K0^2 t^2)`, `powerNoGoReduce=False`;
6. residui numerici:
   `sampleResiduals={{0,0,0},...}`, `maxResidual=0`;
7. classificazione finale:
   `classification="COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY"`.
`Perche questo passaggio`: e il primo step concreto della ricerca EOM complete con materia, con verifica simultanea dei tre canali dinamici.
`Interpretazione`: nel setup hard, esiste un ramo auto-consistente completo (de Sitter con `lambda(t)` che assorbe la materia), mentre un ramo cosmologico power-law standard non soddisfa il canale di vincolo.
`Esito log`: `check=True`.

8. `07_Complete_EOM_Soft_Closure_Evolving_Branch_Test.nb`
`Riferimento`: testare una chiusura soft ispirata al dominio finito per il canale di vincolo e verificare l'esistenza di un ramo evolutivo materia -> de Sitter nei canali completi.
`Output`:
1. chiusura soft normalizzata:
   `3 H^2 - K0^2 = sigma/a^3`,
   `H2soft = K0^2/3 + sigma/(3 a^3)`;
2. coerenza canali completi sul ramo trovato:
   `EaSoft=0`, `ENSoft=0`, `ElamSoft=sigma`, `softConstraintResidual=0`;
3. transizione cosmologica effettiva:
   `weff(a)=-1+sigma/(sigma+K0^2 a^3)`,
   `weffEarly(a=1e-2)=-1.07999883360...*10^-6` (materia-like),
   `weffMid(a=1)=-0.51923076923...`,
   `weffLate(a=1e2)=-0.999999074074...` (de Sitter-like);
4. residui numerici:
   `eaNumVals={0,0,0,0,0}`, `enNumVals={0,0,0,0,0}`, `elamNumVals={0,0,0,0,0}`, `maxResidual=0`;
5. confronto hard:
   `ElamPowerHard=t^(-2+3 p) (3 p^2-K0^2 t^2)`, `hardPowerNoGo=False`;
6. classificazione finale:
   `classification="COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND"`.
`Perche questo passaggio`: e la prima evidenza costruttiva di un ramo evolutivo auto-consistente quando il canale di vincolo e reso soft/finitario invece che hard.
`Interpretazione`: la difficolta trovata nel ramo hard non chiude il framework: con chiusura soft coerente i canali completi ammettono una storia tipo materia->de Sitter.
`Esito log`: `check=True`.

9. `08_Paper_Ready_Hard_Soft_EOM_Bridge_Test.nb`
`Riferimento`: fissare in forma compatta da main-text il ponte tra ramo hard (vincolo rigido) e ramo soft (chiusura finite-domain) nelle EOM complete con materia.
`Output`:
1. forma generale da `EN`:
   `lamGeneric = -( -3 h2 + Lambda + kappa rh)/(6 h2 kappa)`;
2. ramo hard:
   `3 H^2 - K0^2 = 0`,
   `lambda_hard(a)= (K0^2 - Lambda - kappa rho(a))/(2 kappa K0^2)`;
3. ramo soft:
   `3 H^2 - K0^2 = sigma/a^3`,
   `H^2(a)=K0^2/3 + sigma/(3 a^3)`;
4. equivalenza LCDM-like:
   `H^2(a)=H0^2[Omega_m0 a^-3 + Omega_L0]`,
   `Omega_m0=sigma/(K0^2+sigma)`, `Omega_L0=K0^2/(K0^2+sigma)`;
5. EoS efficace:
   `weff(a)=-1+sigma/(sigma+K0^2 a^3)`,
   limiti `weffEarly=0`, `weffLate=-1`;
6. ancore numeriche (sample):
   `aGrid={1/100,1,100}`,
   `weffVals={-1.07999883360...*10^-6,-0.51923076923...,-0.999999074074...}`;
7. check globali:
   `checkA..checkI=True`;
8. classificazione finale:
   `classification="PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED"`.
`Perche questo passaggio`: chiude il rischio di contraddizione narrativa nel paper mostrando che hard e soft non sono due storie scollegate, ma due limiti della stessa struttura EOM.
`Interpretazione`: hard seleziona il ramo de Sitter auto-consistente; soft mantiene la consistenza EOM e riapre una storia evolutiva materia->de Sitter, in forma direttamente presentabile nel testo principale.
`Esito log`: `check=True`.

10. `09_Complete_EOM_Soft_Branch_Parameter_Robustness_Test.nb`
`Riferimento`: verificare che il ramo soft auto-consistente delle EOM complete non sia fine-tuned, ma stabile su un ampio intervallo di parametri.
`Output`:
1. canali completi sul ramo soft:
   `EaSoft=0`, `ENSoft=0`, `ElamSoft=sigma`, `softConstraintResidual=0`;
2. transizione intrinseca del ramo:
   `weff(a)=-1+sigma/(sigma+K0^2 a^3)`,
   `aTransition=(sigma/K0^2)^(1/3)`,
   `weffAtTransition=-1/2`;
3. ancoraggio Planck:
   `rPlanck=omPlanck/olPlanck=0.4486208728750804`,
   `zLambdaMatterPlanck=0.306291721795789`,
   `zSoftPlanck=0.306291721795789`;
4. scan parametrico ampio:
   `scanRatios={1/20,1/10,3/10,9/20,7/10,1,3}`,
   `ratioSpan=60`;
5. risultati scan:
   tutti i sample con `pass=True`,
   `maxResidualScan=0`,
   `zTransitions` da `1.7144...` a `-0.3066...`;
6. check globali:
   `checkA..checkG=True`;
7. classificazione finale:
   `classification="COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED"`.
`Perche questo passaggio`: chiude il punto fisico/metodologico "non e fine-tuning" sul ramo evolutivo trovato in `07-08`.
`Interpretazione`: il ramo soft resta auto-consistente nei canali completi su un intervallo parametrico largo e riproduce in modo naturale la scala di transizione materia-vuoto.
`Esito log`: `check=True`.

11. `10_Pipeline_06_09_Appendix_Ready_Summary_Test.nb`
`Riferimento`: consolidare in un unico test il blocco EOM `06-09` per ottenere una tabella appendix-ready della narrativa hard -> soft -> robustezza.
`Output`:
1. validazione blocco hard (`06`):
   `check06A..check06F=True`,
   `constraintForm(06) = -K0^2 + 3 (a'/a)^2`,
   `hardPowerNoGo06=False`;
2. ancore blocco soft/bridge (`07-08`):
   `weffEarly07=-1.07999883360...*10^-6`,
   `weffMid07=-0.51923076923...`,
   `weffLate07=-0.999999074074...`,
   `lamGeneric08=-( -3 h2 + Lambda + kappa rh)/(6 h2 kappa)`,
   `omegaM0=sigma/(K0^2+sigma)`, `omegaL0=K0^2/(K0^2+sigma)`;
3. ancore robustezza (`09`):
   `zLambdaMatterPlanck=0.306291721795789`,
   `zSoftPlanck=0.306291721795789`,
   `scanRatios09={1/20,1/10,3/10,9/20,7/10,1,3}`,
   `ratioSpan09=60`,
   `maxResidualScan09=0`;
4. tabella aggregata:
   `rows={("06","COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY"),("07","COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND"),("08","PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED"),("09","COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED")}`;
5. check globali:
   `checkAllRows=True`, `checkClasses=True`;
6. classificazione finale:
   `classification="PIPELINE_06_09_APPENDIX_READY"`.
`Perche questo passaggio`: chiude in modo editoriale il blocco EOM con una singola evidenza riproducibile da citare in appendice.
`Interpretazione`: il percorso logico regge end-to-end: hard seleziona il ramo de Sitter, soft riapre l'evoluzione materia->de Sitter, e la soluzione soft resta robusta senza fine-tuning.
`Esito log`: `check=True`.

12. `11_Section_Closure_Referee_Ready_Test.nb`
`Riferimento`: test di chiusura sezione in formato referee-facing, con checklist sintetica dei claim fisici principali.
`Output`:
1. claim hard verificato:
   `constraintForm = -K0^2 + 3 (a'/a)^2`,
   `hardPowerNoGo=False`,
   `checkA..checkE=True`;
2. claim soft verificato:
   `weffEarly=-1.07999883360...*10^-6`,
   `weffMid=-0.51923076923...`,
   `weffLate=-0.999999074074...`,
   `checkF..checkI=True`;
3. ancoraggio standard:
   `zLambdaMatterPlanck=0.306291721795789`,
   `zSoftPlanck=0.306291721795789`,
   `checkJ=True`;
4. anti-fine-tuning compatto:
   `scanRatios={1/20,1/5,1,3}`,
   `ratioSpan=60`,
   `scanData pass=True` per tutti i punti,
   `maxResidualScan=0`,
   `checkK=True`;
5. snippet e checklist pronti:
   `refereeClosureSnippet="Complete-EOM closure test: ..."` ,
   `closureChecklistMarkdown="| Claim | Status | ..."` ;
6. classificazione finale:
   `classification="COSMOLOGICAL_ERAS_SECTION_CLOSURE_READY"`.
`Perche questo passaggio`: fornisce una chiusura breve e inattaccabile da integrare direttamente nel testo/appendice.
`Interpretazione`: i tre messaggi chiave sono consolidati in un solo test: hard coerente, soft evolutivo, robustezza non fine-tuned.
`Esito log`: `check=True`.

## Nota operativa

Generazione e verifica:

```powershell
cd c:\lean\StaticSanityCheck\Paper\Cosmological_Eras_Study
.\run_all.ps1
```

Risultati:
1. notebook in `Paper/Cosmological_Eras_Study/*.nb`
2. log in `Paper/Cosmological_Eras_Study/logs/*.log`
