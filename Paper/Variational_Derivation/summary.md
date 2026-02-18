# Variational Derivation Summary

Mappa notebook `Paper/Variational_Derivation` -> passaggi matematici, output verificati e lettura fisica.

## Notebook 00-25

1. `00_Variational_Seed.nb`
`Riferimento`: azione con moltiplicatore di Lagrange, `eq:lag_mult` (`Paper/main.tex:301`).
`Output`: `eqLambda = I1[th] - K0^2 I0[th]`, `check1=True`, `check2=True`.
`Perche questo passaggio`: fissa subito la condizione di vincolo corretta da `delta S / delta lambda`.
`Interpretazione fisica dell'output`: il vincolo e' leafwise e globale; non impone una legge locale puntuale su ogni punto dello spazio.

2. `01_Ratio_FirstVariation.nb`
`Riferimento`: media leafwise `I1/I0` (sezione 4, `Paper/main.tex:269`).
`Output`: `deltaRatio = (di1*i0 - di0*i1)/i0^2`, `check=True`.
`Perche questo passaggio`: serve la formula esatta di variazione del rapporto per tutte le EOM successive.
`Interpretazione fisica dell'output`: la risposta dinamica dipende sia da variazione del numeratore (curvatura) sia del denominatore (volume leafwise), con pesatura `1/I0^2`.

3. `02_Constraint_FirstVariation.nb`
`Riferimento`: termine vincolato dell'azione, `Paper/main.tex:292`-`Paper/main.tex:300`.
`Output`: decomposizione corretta della variazione in parte `eta` e parte `dI0,dI1`, `checkA=True`, `checkB=True`.
`Perche questo passaggio`: separa algebraicamente "enforcement del vincolo" da "risposta geometrica".
`Interpretazione fisica dell'output`: il meccanismo ha due canali: uno che impone la regola globale e uno che modifica la geometria media della foglia.

4. `03_Theta_ChainRule.nb`
`Riferimento`: dipendenza da `Theta` di `lambda`, `I0`, `I1` (sezione 4).
`Output`: chain rule completa, `checkA=True`, `checkB=True`.
`Perche questo passaggio`: evita errori nel contributo `delta Theta` quando i funzionali dipendono dalla foliazione.
`Interpretazione fisica dell'output`: il campo `Theta` agisce come coordinata d'ordine temporale della foliazione; le variazioni scorrono lungo le foglie.

5. `04_Ratio_Kernel_Form.nb`
`Riferimento`: struttura IR del vincolo (sezione 4.4).
`Output`: forma kernel
`deltaQ = (I0*Int[k1 dg] - I1*Int[k0 dg])/I0^2`, `checkA=True`, `checkB=True`.
`Perche questo passaggio`: mette la variazione in forma integrale pronta per stime IR/UV.
`Interpretazione fisica dell'output`: il termine efficace e' il contrasto tra kernel locale e media globale (`k1 - q k0`), quindi le correzioni locali sono naturalmente filtrate dalla media cosmica.

6. `05_Metric_Variation_Spherical.nb`
`Riferimento`: `K2_general = NumK2/DenK2` (`Paper/main.tex:710`) + test sezione 6.
`Output`: `checkA=True`, `k2Minkowski = 2 v^2/(r^2 (1-v^2))`, `k2dS = 3 H^2`, `checkB=True`, `checkC=True`.
`Perche questo passaggio`: aggancia la pipeline variazionale alla stessa forma chiusa usata in Lean.
`Interpretazione fisica dell'output`: conferma esplicita della selezione geometrica: Minkowski non mantiene `K^2` costante positivo nel limite IR, de Sitter si'.

7. `06_SecondOrder_EL_Identity.nb`
`Riferimento`: identita' generale per lagrangiane con derivate seconde.
`Output`: `residue = 0`, `check=True`.
`Perche questo passaggio`: fondazione formale per non perdere termini di bordo quando compare `psi''`.
`Interpretazione fisica dell'output`: la dinamica bulk (`EL`) e i flussi di bordo sono separabili in modo coerente; senza questa separazione la fisica estratta sarebbe ambigua.

8. `07_K2_Psi_Variation_Numeric.nb`
`Riferimento`: applicazione dell'identita' precedente al blocco reale `K2_general`.
`Output`: `samplePts={0.7,0.9,1.1,1.3}`, `residualValues={0,0,0,0}`, `maxResidual=0`, `check=True`.
`Perche questo passaggio`: test operativo robusto entro il limite temporale del kernel cloud.
`Interpretazione fisica dell'output`: la formula variazionale usata per `psi` e' consistente nel regime nonlineare del modello, non solo in casi toy lineari.

9. `08_K2_Psi_EL_Operator.nb`
`Riferimento`: operatore `EL_psi` esplicito del blocco `K2_general`.
`Output`: `ELvalues={2.4306327483, 1.1798451311, 0.6699496789, 0.4263114587}`, `ELleafCount=209`, `maxResidual=0`, `check=True`.
`Perche questo passaggio`: oltre al solo check, espone l'operatore dinamico che guida `psi`.
`Interpretazione fisica dell'output`: nel background testato la forza efficace su `psi` e' positiva e decrescente con `r`; indica che il settore di foliazione non e' inerte ma ha una dinamica IR ben definita.

10. `09_Metric_Kernels_AB_Numeric.nb`
`Riferimento`: primi kernel metrici nel settore sferico.
`Output`: `EAvalues={0.1647144981,0.1013697265,0.0659715859,0.0424928182}`, `EBvalues={-0.1875012281,-0.1296360617,-0.09978467,-0.0825851268}`, `EAleafCount=199`, `EBleafCount=191`, `maxResidual=0`, `check=True`.
`Perche questo passaggio`: e' l'ingresso diretto verso le EOM metriche (canali `A` e `B`) con decomposizione variazionale verificata.
`Interpretazione fisica dell'output`: i due canali metrici rispondono in modo diverso (segno opposto nel test), quindi il vincolo non agisce come semplice riscalatura uniforme della metrica ma come sorgente strutturata sul settore radiale.

11. `10_Full_Proxy_EOM_Operators.nb`
`Riferimento`: assemblaggio operativo della densita' ridotta locale+vincolo proxy nel settore sferico.
`Output`: `EAvalues={-2.7927302932,-2.6768087759,-2.5390690820}`, `EBvalues={-2.3561992051,-2.5444503443,-2.7757012519}`, `Epsivalues={-1.2472008933,-1.8565167768,-2.5176554943}`, `EAleafCount=1309`, `EBleafCount=1117`, `EpsileafCount=1127`, `check=True`.
`Perche questo passaggio`: produce esplicitamente gli operatori EOM dei tre canali (`A`,`B`,`psi`) in un unico blocco coerente.
`Interpretazione fisica dell'output`: nel background testato tutti e tre i canali sono dinamicamente attivi e comparabili in scala; la dinamica non e' confinata al solo settore `psi`.

12. `11_Full_Proxy_Variational_Checks.nb`
`Riferimento`: verifica di consistenza variazionale canale-per-canale sulla densita' completa proxy.
`Output`: `resAvals~O(10^-15)`, `resBvals~O(10^-16)`, `resPsivals~O(10^-16)`, `maxResA=1.3361305999*10^-15`, `maxResB=5.2489096678*10^-16`, `maxResPsi=3.3869698479*10^-16`, `checkA=True`, `checkB=True`, `checkPsi=True`, `check=True`.
`Perche questo passaggio`: certifica che la decomposizione `deltaL = EOM*delta(field) + total-derivative` regge simultaneamente per tutti i campi del proxy.
`Interpretazione fisica dell'output`: il framework variazionale completo e' internamente consistente a precisione numerica molto alta; i termini di bordo e bulk sono bilanciati come richiesto da una teoria ben posta.

13. `12_Finite_Domain_Ratio_Variation.nb`
`Riferimento`: formulazione a dominio finito del rapporto `Q=I1/I0` con bordo mobile `R -> R + eps rho`.
`Output`: `deltaQ = (dI1bulk*I0 - dI0bulk*I1 + f1R*I0*rho - f0R*I1*rho)/I0^2`, `boundaryPiece = (f1R*I0 - f0R*I1)/I0^2`, `deltaNum=0.08712283881742730...`, `deltaPred=0.08712283881742621...`, `absErr=1.0902897015*10^-15`, `checkA=True`, `checkB=True`, `checkNum=True`, `check=True`.
`Perche questo passaggio`: chiude il punto critico emerso nel dibattito ("cosa succede a epoche/domini finiti?") verificando che la variazione includa correttamente il contributo di bordo oltre al bulk.
`Interpretazione fisica dell'output`: il flusso del vincolo su domini finiti e' controllato da due pezzi distinti: variazione interna della foglia e "iniezione" dal bordo. L'errore numerico ~`10^-15` mostra che la formula analitica di bordo non e' un'approssimazione qualitativa ma coincide con la derivata effettiva del rapporto.

14. `13_IR_Flow_Minkowski_DeSitter.nb`
`Riferimento`: legge di flusso in raggio `dQ/dR` e confronto IR tra settore tipo Minkowski e settore tipo de Sitter.
`Output`: `qPrime = (f1R*I0 - f0R*I1)/I0^2`, `qM(R)=3*c/R^2`, `dqM(R)=-6*c/R^3`, `qD(R)=KdS`, `dqD(R)=0`, `qMvals(c=0.8)={0.6,0.15,0.0375}`, `qDvals(KdS=1.2)={1.2,1.2,1.2}`, `checkA=True`, `checkM=True`, `checkD=True`, `check=True`.
`Perche questo passaggio`: serve a rendere esplicita la dinamica IR del vincolo con una variabile fisicamente trasparente (`R`) e a separare due comportamenti-limite del rapporto medio.
`Interpretazione fisica dell'output`: nel toy Minkowski il rapporto decade come `1/R^2` (quindi tende a zero in IR), mentre nel toy de Sitter resta costante. Questo formalizza il meccanismo di selezione: l'IR spinge naturalmente fuori dal settore statico e privilegia un settore a rapporto costante positivo.

15. `14_Constraint_Finite_Domain_Boundary.nb`
`Riferimento`: variazione del funzionale vincolato lineare `C=I1-K0^2 I0` su dominio finito con bordo mobile.
`Output`: `deltaC = dI1bulk - dI0bulk*K0^2 + f1R*rho - f0R*K0^2*rho`, `bulkPart = dI1bulk - dI0bulk*K0^2`, `boundaryPart = f1R - f0R*K0^2`, `deltaNum=0.1006720000...`, `deltaPred=0.100672`, `absErr=7.869375...*10^-16`, `checkA=True`, `checkB=True`, `checkBoundaryNeutral=True`, `checkNum=True`, `check=True`.
`Perche questo passaggio`: rende espliciti i termini di bordo nel canale del moltiplicatore di Lagrange (non solo nel rapporto `I1/I0`).
`Interpretazione fisica dell'output`: il bordo contribuisce tramite mismatch locale `f1R-K0^2 f0R`; se il bordo e' gia' allineato al vincolo, il flusso di bordo si spegne.

16. `15_FRW_K2_Leafwise_Average.nb`
`Riferimento`: formulazione cosmologica omogenea (riduzione FRW del vincolo leafwise).
`Output`: `K2 = 3*ad^2/(a^2*N^2)`, `I0=(4*a^3*Pi*R^3)/3`, `Q=I1/I0=3*ad^2/(a^2*N^2)`, `QdS=3*H^2`, `dQdt_dS=0`, `selectionRule(H^2)={h2->K0^2/3}`, `numErr=0`, `check=True`.
`Perche questo passaggio`: costruisce la base matematica per scrivere EOM cosmologiche coerenti con il vincolo.
`Interpretazione fisica dell'output`: in un background omogeneo il vincolo globale coincide con il valore locale, e il ramo de Sitter realizza naturalmente `H^2=K0^2/3`.

17. `16_FRW_Minisuperspace_EOM.nb`
`Riferimento`: EOM complete nel settore minisuperspace FRW per i campi `(a,N,lambda)`.
`Output`: `EaDS = (-3*E^(2*H0*t)*(K0^2*kappa*lam0 + 3*H0^2*(-1 + kappa*lam0) + Lambda))/kappa`, `ENDS = -(E^(3*H0*t)*(H0^2*(-3 + 6*kappa*lam0) + Lambda))/kappa`, `ElamDS=E^(3*H0*t)*(3*H0^2-K0^2)`, `lamRule={lam0->(K0^2-Lambda)/(2*K0^2*kappa)}`, `eaNum=0`, `enNum=0`, `elamNum=0`, `check=True`.
`Perche questo passaggio`: fornisce una derivazione variazionale cosmologica completa (non solo qualitativa) nel proxy FRW.
`Interpretazione fisica dell'output`: il ramo de Sitter soddisfa simultaneamente i tre canali EOM; `lambda` fissa il bilanciamento tra settore locale e vincolo globale.

18. `17_FRW_Linearized_Constraint_Channel.nb`
`Riferimento`: linearizzazione del vincolo cosmologico attorno a de Sitter.
`Output`: `C0=E^(3*H0*t)*(3*H0^2-K0^2)`, `C1=-3*E^(3*H0*t)*(2*H0^2*n[t]+(-3*H0^2+K0^2)*u[t]-2*H0*u'[t])`, `C1_onBranch=-2*E^(3*H0*t)*(K0^2*n[t]-3*H0*u'[t])`, `C1_gauge_n0=6*E^(3*H0*t)*H0*u'[t]`, `c1Vals(constrained)={0,0,0,0}`, `maxAbs=0`, `check=True`.
`Perche questo passaggio`: testa la propagazione del vincolo nel canale lineare cosmologico.
`Interpretazione fisica dell'output`: al primo ordine il vincolo impone `u' = H n`; in gauge cosmica (`n=0`) il modo scalare del canale vincolato e' congelato (`u'=0`), quindi non emerge una dinamica lineare libera in quel settore.

19. `18_deSitter_Constraint_Propagation.nb`
`Riferimento`: coerenza dinamica dei vincoli nel minisuperspace FRW sul ramo de Sitter.
`Output`: `EaOn=0`, `ENOn=0`, `ElamOn=0`, `dEaOn=0`, `dENOn=0`, `dElamOn=0`, `lamRule={lam0->(K0^2-Lambda)/(2*K0^2*kappa)}`, `checkA..F=True`, `maxAbsNum=0`, `check=True`.
`Perche questo passaggio`: verifica che una volta imposto il ramo de Sitter con il valore consistente di `lambda`, i vincoli non solo sono soddisfatti ma restano propagati nel tempo.
`Interpretazione fisica dell'output`: il ramo selezionato e' dinamicamente autosufficiente nel settore omogeneo; non emergono derive temporali spurie dei canali vincolati.

20. `19_FRW_Linearized_Full_EOM_Reduction.nb`
`Riferimento`: linearizzazione completa del sistema `(Ea,EN,Elam)` e riduzione per eliminazione dei campi vincolati.
`Output`: `nRule={n[t]->(Sqrt[3] u'[t])/K0}`, `lRule={l1[t]->(Sqrt[3](K0^2-Lambda)u'[t])/(2*K0^3*kappa)}`, `dnRule`, `dlRule`, `eaReduced=0`, `resVals={0,0,0}`, `maxAbs=0`, `check=True`.
`Perche questo passaggio`: completa la chiusura EOM lineare mostrando esplicitamente come i canali di vincolo eliminano le variabili ausiliarie nel settore FRW.
`Interpretazione fisica dell'output`: al primo ordine, sul ramo de Sitter omogeneo, l'equazione dinamica residua nel canale scalare si riduce a identita' (nessuna forza lineare residua in questo proxy), coerente con il carattere vincolato del settore.

21. `20_Minisuperspace_Noether_Bianchi_Identity.nb`
`Riferimento`: identita' strutturale da riparametrizzazione temporale nel minisuperspace FRW.
`Output`: `noetherRaw = d/dt(lambda*Elam)`, `auxTerm = d/dt(lambda*Elam)`, `noetherId = 0`, `propEq = 0`, `idDS = 0`, `idVals={0,0,0,0}`, `maxAbs=0`, `checkA=True`, `checkB=True`, `checkC=True`, `checkNum=True`, `check=True`.
`Perche questo passaggio`: lega formalmente i canali EOM in una relazione di consistenza (analogo minisuperspace di una identita' di tipo Bianchi/Noether).
`Interpretazione fisica dell'output`: la coerenza dei vincoli non e' solo numerica su un branch, ma discende da una identita' off-shell; il settore del moltiplicatore entra come termine di bilancio totale `d/dt(lambda*Elam)`.

22. `21_FRW_Curvature_Minisuperspace_EOM.nb`
`Riferimento`: estensione FRW minisuperspace con curvatura spaziale `k/a^2` tramite termine locale `+3 k a N/kappa`.
`Output`: `checkA..C=True` (recupero limite `k->0` delle EOM del notebook 16), `checkD=True` (`dEa/dk=3N/kappa`), `checkE=True` (`dEN/dk=3a/kappa`), `checkF=True` (`Elam` indipendente da `k`), `checkNum=True`, `check=True`.
`Perche questo passaggio`: testa che la catena EOM non si rompe aggiungendo la curvatura spaziale standard GR.
`Interpretazione fisica dell'output`: il termine di curvatura modifica solo i canali metrici (`Ea`,`EN`) come atteso, mentre il canale di vincolo (`Elam`) resta invariato.

23. `22_FRW_Curvature_Noether_Identity.nb`
`Riferimento`: identita' Noether/Bianchi minisuperspace con curvatura spaziale non nulla.
`Output`: `noetherId=0`, `checkA=True`, `checkB(k->0)=True`, `propEq=0`, `checkC=True`, `idVals={0,0,0,0}`, `maxAbs=0`, `checkNum=True`, `check=True`.
`Perche questo passaggio`: verifica che la relazione strutturale tra EOM e vincoli (notebook 20) sopravvive anche con `k != 0`.
`Interpretazione fisica dell'output`: la coerenza di gauge temporale del minisuperspace e' robusta rispetto alla curvatura spaziale; non compaiono anomalie nella propagazione dei vincoli.

24. `23_Covariant_Causal_Domain_Variation.nb`
`Riferimento`: definizione del dominio causale leafwise e funzionali `I0`,`I1` (`Paper/main.tex:242-263`) + struttura variazionale del rapporto (`Paper/main.tex:317-324`).
`Output`: `dj0 = chi*dmu + dchi*mu`, `dj1 = chi*dmu*k2 + chi*dk2*mu + dchi*k2*mu`, `bulk0=chi*dmu`, `boundary0=mu`, `bulk1=chi*(dmu*k2 + dk2*mu)`, `boundary1=k2*mu`, `boundaryQ=(b1*I0-b0*I1)/I0^2`, `checkA..G=True`, `check=True`.
`Perche questo passaggio`: porta i blocchi a una formulazione 4D covariante esplicita con indicatore di dominio `chi_Theta` (non solo proxy minisuperspace).
`Interpretazione fisica dell'output`: la variazione si separa in modo pulito tra contributi bulk (variazioni metriche/di `K^2`) e contributi di bordo (variazione del dominio causale), mantenendo la stessa struttura del flusso `deltaQ` usata nei notebook 12-14.

25. `24_Covariant_Boundary_Shape_Derivative.nb`
`Riferimento`: domini finiti/causali e saturazione a scala di orizzonte (`Paper/main.tex:278-281`, `Paper/main.tex:495-497`).
`Output`: `dchi = rho*DiracDelta[-r+R]`, `dI0Delta=dI0Leib=((5+R^2)rho)/(5 E^R)`, `dI1Delta=dI1Leib=((9+2R)(5+R^2)rho)/(30 E^R)`, `err0~1.09e-17`, `err1~1.25e-17`, `checkA=True`, `checkB=True`, `checkC=True`, `checkNum0=True`, `checkNum1=True`, `check=True`.
`Perche questo passaggio`: dimostra in modo esplicito che la variazione dell'indicatore di dominio genera il termine di bordo distribuzionale corretto.
`Interpretazione fisica dell'output`: il bordo non e' un'aggiunta ad hoc, ma emerge geometricamente come termine di forma (`DiracDelta`) e riproduce numericamente la stessa correzione finite-domain gia' osservata nella catena proxy.

26. `25_Covariant_Nonlocal_Action_Variation.nb`
`Riferimento`: azione non-locale con moltiplicatore (`Paper/main.tex:292-301`) + dominanza IR del canale normalizzato (`Paper/main.tex:321-329`).
`Output`: `deltaSleaf = eta*(i1-k0^2 i0) + lam*(di1-k0^2 di0)`, `di0Cov=A0+B0`, `di1Cov=A1+B1+C1`, `deltaC=(A1+C1-k0^2 A0) + (B1-k0^2 B0)`, `deltaQIR=(A1+B1+C1-A0 q-B0 q)/i0`, `vals(L={10,20,40})={1.555e-3,1.94375e-4,2.4296875e-5}`, `ratios~{8,8}`, `checkA..E=True`, `check=True`.
`Perche questo passaggio`: chiude la derivazione variazionale in forma coerente con l'azione del paper, distinguendo enforcement (`delta lambda`) da risposta geometrica bulk+boundary.
`Interpretazione fisica dell'output`: nel canale della media normalizzata le correzioni decadono come `1/I0` (qui `~L^-3`), quindi la modifica non-locale resta IR-dominata e localmente soppressa, in linea con l'argomento di recupero GR locale nel paper.

Nota di scope:
1. La catena `23-25` copre la variazione covariante 4D una volta introdotto l'indicatore `chi_Theta` del dominio causale.
2. Resta separato (come gia' dichiarato nel paper) il problema della dimostrazione generale di equivalenza tra `chi_Theta` e una costruzione esplicita via mutua raggiungibilita' causale in background arbitrari.

## Nota operativa

Generazione e verifica (Wolfram Cloud):

```powershell
cd c:\lean\StaticSanityCheck\Paper\Variational_Derivation
.\run_all.ps1
```

I risultati sono in `Paper/Variational_Derivation/logs/*.log` e i notebook in `Paper/Variational_Derivation/*.nb`.

## Nota EFT / scale di validita
1. Questa pipeline fissa la struttura variazionale (bulk/boundary) e la soppressione IR (`~1/I0`), ma non determina da sola un cutoff UV fondamentale.
2. La discussione operativa su cutoff, strong coupling e regime di validita e' stata implementata in `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`.
3. In termini EFT, questa cartella fornisce i building blocks; la stima quantitativa del dominio affidabile avviene nel blocco perturbativo completo.
