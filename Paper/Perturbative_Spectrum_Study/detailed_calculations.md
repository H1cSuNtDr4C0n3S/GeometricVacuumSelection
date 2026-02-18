# Perturbative Spectrum - Detailed Calculations and Results

Questo documento dettaglia i calcoli in `Paper/Perturbative_Spectrum_Study/scripts/` e i risultati in `Paper/Perturbative_Spectrum_Study/logs/`.

Convenzioni:
1. `Kmat`: matrice cinetica scalare.
2. `Gmat`: matrice dei termini gradiente.
3. `c^2`: autovalori di `Kmat^-1 Gmat` ad alta frequenza.
4. `I0`: scala IR del dominio causale leafwise (come nel paper).

## 00 - Tensor Channel IR Stability

Riferimento:
1. Script: `scripts/00_tensor_channel_ir_stability.wl`.
2. Log: `logs/00_tensor_channel_ir_stability.log`.

Calcolo:
1. Si usa la forma quadratica tensoriale:
   `S_T^(2) ~ Integral a^3 [GT hdot^2 - FT (grad h)^2/a^2]`.
2. Si definisce:
   `GT = 1/kappa + dGT`, `FT = 1/kappa + dFT`.
3. Si ottiene:
   `cT2 = FT/GT`.
4. Si testa il limite GR (`dGT=dFT=0`).
5. Si impone una parametrizzazione IR-soppressa:
   `dGT=aGT/I0`, `dFT=aFT/I0`.

Output:
1. `cT2 = (1 + dFT*kappa)/(1 + dGT*kappa)`.
2. `cT2IR = (I0 + aFT*kappa)/(I0 + aGT*kappa)`.
3. `deltaCT = ((aFT-aGT)*kappa)/(I0 + aGT*kappa)`.
4. `deltaAsym = (aFT-aGT)*kappa`.
5. Test numerici su `I0={10,20,40}`: `cT2` si avvicina a `1`.
6. `check=True`.

Lettura fisica:
1. Ghost-free tensor: `GT>0`.
2. Gradient-stable tensor: `FT>0`.
3. La deviazione di velocita dei tensori e' naturalmente soppressa in IR se le correzioni entrano come `1/I0`.

## 01 - Scalar Mixing Stability Matrix

Riferimento:
1. Script: `scripts/01_scalar_mixing_stability_matrix.wl`.
2. Log: `logs/01_scalar_mixing_stability_matrix.log`.

Calcolo:
1. Sistema scalare quadratico a due campi con:
   `Kmat={{k11,k12},{k12,k22}}`,
   `Gmat={{g11,g12},{g12,g22}}`.
2. Si calcolano:
   `detK`, `traceKG = Tr[K^-1 G]`, `detKG = Det[K^-1 G]`.
3. Si scrivono gli autovalori:
   `c2plus`, `c2minus`.
4. Si testa un set stabile e un set instabile.

Output:
1. `detK = -k12^2 + k11*k22`.
2. `traceKG = (g22*k11 - 2*g12*k12 + g11*k22)/detK`.
3. `detKG = (g11*g22 - g12^2)/detK`.
4. Esempio stabile: autovalori positivi `{0.8318..., 0.7463...}`.
5. Esempio instabile: un autovalore negativo `{0.7502..., -0.4042...}`.
6. `check=True`.

Lettura fisica:
1. Ghost-free: `Kmat` positiva definita (`k11>0`, `detK>0`).
2. Gradient stability (alta frequenza): autovalori di `K^-1 G` positivi.
3. Il mixing entra esplicitamente nelle soglie tramite `k12`, `g12`.

## 02 - Frozen Scalar Reduction

Riferimento:
1. Script: `scripts/02_frozen_scalar_reduction.wl`.
2. Log: `logs/02_frozen_scalar_reduction.log`.

Calcolo:
1. Si considera un Lagrangiano quadratico con modo scalare non dinamico `pi`:
   `L2 = 1/2 A zd^2 + B zd pi - 1/2 (C+F q2) pi^2 - E q2 z pi - 1/2 D q2 z^2`.
2. Si integra `pi` tramite equazione algebrica:
   `piRule = (-(E q2 z) + B zd)/(C+F q2)`.
3. Si ottiene l'azione efficace:
   `Keff`, `crossEff`, `Geff`.
4. Si studiano i limiti:
   `q2 -> Infinity` (UV) e `C -> Infinity` (freeze hard).

Output:
1. `Keff = A + B^2/(C+F q2)`.
2. `crossEff = -(B E q2)/(C+F q2)`.
3. `Geff = D - (E^2 q2)/(C+F q2)`.
4. UV: `Keff->A`, `Geff->D-E^2/F`.
5. Freeze hard: `Keff->A`, `Geff->D`, `cross->0`.
6. `check=True`.

Lettura fisica:
1. Integrare il modo non dinamico e' ben posto.
2. Il criterio gradient UV del modo residuo e' `D-E^2/F>0`.
3. Nel limite hard (`C` grande) si recupera una dinamica scalare ridotta senza mixing residuo.

## 03 - FRW Scalar Constraint Freeze Map

Riferimento:
1. Script: `scripts/03_frw_scalar_constraint_freeze_map.wl`.
2. Log: `logs/03_frw_scalar_constraint_freeze_map.log`.

Calcolo:
1. Si riparte dalla struttura del canale lineare FRW:
   `C1onBranch = -2 Exp[3 H0 t](K0^2 n1 - 3 H0 up)`.
2. Sul branch `H0=K0/Sqrt[3]` si ottiene:
   `n1 = (Sqrt[3]/K0) up`.
3. Impostando gauge `n1=0`, si ricava `up=0`.
4. Su un kinetic toy `Lkin=a1 up^2 + a2 up n1 + a3 n1^2`,
   la riduzione vincolata elimina `n1` e in gauge cosmico annulla il contributo.

Output:
1. `nRule = {n1 -> (Sqrt[3]*up)/K0}`.
2. `upRuleGauge = {up -> 0}`, `upGauge = 0`.
3. `LkinRed = ((3*a3 + K0*(Sqrt[3]*a2 + a1*K0))*up^2)/K0^2`.
4. `LkinGauge = 0`.
5. `check=True`.

Lettura fisica:
1. Il freeze scalare osservato nella catena EOM FRW e' coerente anche a livello di riduzione perturbativa.
2. In gauge cosmico il modo scalare del canale vincolato non propaga.

## 04 - Covariant Ratio Second Variation

Riferimento:
1. Script: `scripts/04_covariant_ratio_second_variation.wl`.
2. Log: `logs/04_covariant_ratio_second_variation.log`.

Calcolo:
1. Si definisce l'espansione a due ordini del rapporto non-locale:
   `Q(eps) = (i1 + eps di1 + eps^2 d2i1/2)/(i0 + eps di0 + eps^2 d2i0/2)`.
2. Si derivano:
   `q1 = dQ/deps|0`, `q2 = d2Q/deps2|0`.
3. Si passa al branch `i1 = q i0` (coerente con vincolo leafwise su soluzione).
4. Si testa la scalatura IR con `i0=L^3` su `L={10,20,40}`.

Output:
1. `q1 = (di1*i0 - di0*i1)/i0^2`, `checkA=True`.
2. `q2 = (-2*di0*di1*i0 + 2*di0^2*i1 + i0*(d2i1*i0 - d2i0*i1))/i0^3`, `checkB=True`.
3. `q1Branch = (di1 - di0*q)/i0`, `checkC=True`.
4. `q2Branch = (-2*di0*di1 + 2*di0^2*q + i0*(d2i1 - d2i0*q))/i0^2`, `checkD=True`.
5. `q1Vals = {2.8311688e-4, 3.5389610e-5, 4.4237013e-6}`.
6. `q2Vals = {2.6428266e-4, 3.3053028e-5, 4.1319049e-6}`.
7. `ratioQ1 ~ {8,8}`, `ratioQ2 ~ {7.9957,7.9995}`.
8. `check=True`.

Lettura fisica:
1. I termini al primo e secondo ordine del rapporto non-locale restano entrambi soppressi dal volume IR (`~1/i0` con `i0~L^3` nel setup testato).
2. Questo e' il blocco algebrico che serve per passare dalle EOM lineari a un kernel quadratico non-locale controllato.

## 05 - Nonlocal Quadratic Mixing Kernel

Riferimento:
1. Script: `scripts/05_nonlocal_quadratic_mixing_kernel.wl`.
2. Log: `logs/05_nonlocal_quadratic_mixing_kernel.log`.

Calcolo:
1. Si parametrizzano due perturbazioni scalari `x,y`:
   `dI0 = a0 x + b0 y`, `dI1 = a1 x + b1 y`,
   `d2I0 = aa0 x^2 + 2 ab0 x y + bb0 y^2`,
   `d2I1 = aa1 x^2 + 2 ab1 x y + bb1 y^2`.
2. Si usa la formula covariante di `q2`:
   `q2 = (d2I1 - q d2I0)/i0 - 2 dI0 (dI1 - q dI0)/i0^2`.
3. Si estraggono i coefficienti del kernel:
   `Mxx`, `Mxy`, `Myy`.
4. Si studia il limite IR della matrice `Kmix`.

Output:
1. Forme chiuse dei coefficienti verificate:
   `checkA=True`, `checkB=True`, `checkC=True`.
2. Limite leading:
   `Klead = {{aa1-aa0*q, ab1-ab0*q},{ab1-ab0*q, bb1-bb0*q}}`, `checkD=True`.
3. Test numerico (`i0=L^3`, `L={10,20,40}`):
   autovalori positivi e decrescenti in modulo.
4. Confronto `L^3 Kmix` vs `Klead`:
   `errLead = 9.375e-8`, `checkNum=True`.
5. `check=True`.

Lettura fisica:
1. Il mixing non-locale quadratico ha una struttura matriciale esplicita e regolare.
2. In IR il kernel decresce come `1/i0`, mentre la matrice riscalata `i0 Kmix` tende a un limite finito: questo e' compatibile con un effetto non-locale globalmente rilevante ma localmente soppresso.

## 06 - Second-Order Boundary Shape Derivative

Riferimento:
1. Script: `scripts/06_second_order_boundary_shape_derivative.wl`.
2. Log: `logs/06_second_order_boundary_shape_derivative.log`.

Calcolo:
1. Dominio causale mobile descritto da:
   `chiEps(r)=HeavisideTheta[R + eps rho + eps^2 sigma/2 - r]`.
2. Derivata prima e seconda dell'indicatore:
   `dchi1`, `dchi2`.
3. Confronto tra:
   derivazione Leibniz dell'integrale a bordo mobile e formulazione distribuzionale con DiracDelta/derivata di DiracDelta.
4. Verifica numerica con differenze finite.

Output:
1. `dchi1 = rho*DiracDelta[r-R]`, `checkA=True`.
2. `dchi2 = sigma*DiracDelta[r-R] + 4 rho^2 Derivative[1][DiracDelta][-2 r + 2 R]`.
3. Primo ordine: `i1Leib=i1Dist`, `checkB=True`.
4. Secondo ordine: `i2Leib=i2Dist=i2Expected`, `checkC=True`, `checkD=True`.
5. Test numerico: `err1 ~ 6.4e-12`, `err2 ~ 2.7e-12`, `checkNumA=True`, `checkNumB=True`.
6. `check=True`.

Lettura fisica:
1. I contributi di bordo al secondo ordine emergono direttamente dalla variazione del dominio causale (non sono ipotesi aggiunte a mano).
2. Questo chiude il pezzo tecnico piu delicato per trattare in modo consistente il cutoff causale variabile nel blocco perturbativo quadratico.

## 07 - Spin-Projector Quadratic Operator (4D Proxy)

Riferimento:
1. Script: `scripts/07_spin_projector_quadratic_operator.wl`.
2. Log: `logs/07_spin_projector_quadratic_operator.log`.

Calcolo:
1. Si costruiscono i proiettori di Barnes-Rivers (spaziali) per un tensore simmetrico:
   `P2`, `P1`, `P0s`, `P0w`, `P0sw`, `P0ws`.
2. Si verifica l'algebra:
   idempotenza, ortogonalita e completezza `P2+P1+P0s+P0w=Isym`.
3. Si definisce un operatore isotropo generale:
   `Op = c2 P2 + c1 P1 + cs P0s + cw P0w + cm(P0sw+P0ws)`.
4. Si verifica che il settore tensoriale non miscela con lo scalare:
   `P2 Op P0* = 0`, e che `P2 Op P2 = c2 P2`.
5. Si estrae il blocco scalare equivalente:
   `Kscalar={{cs,cm},{cm,cw}}`, `detKscalar=cs*cw-cm^2`.
6. Si introducono correzioni IR-soppresse:
   `cs=cs0+ns/I0`, `cw=cw0+nw/I0`, `cm=cm0+nm/I0`, `c2=c20+n2/I0`.
7. Si controllano limiti `I0->Infinity` e scaling numerico con `I0=L^3`.

Output:
1. `checkA..checkJ = True` (algebra proiettori e completezza).
2. `checkK=True` (assenza di mixing tensor-scalare nell'operatore isotropo).
3. `checkL=True` (`P2 Op P2 = c2 P2`).
4. `traceKscalar = cs + cw`, `detKscalar = -cm^2 + cs*cw`, `checkM=True`.
5. Limiti IR: `traceIRInf=cs0+cw0`, `detIRInf=-cm0^2+cs0*cw0`, `c2IRInf=c20`,
   con `checkN=True`, `checkO=True`, `checkP=True`.
6. Esempio numerico stabile:
   autovalori scalari positivi per `L={10,20,40}` e `tensorVals>0`.
7. Errori rispetto al limite IR:
   `scalarErrRatios={8,8}`, `tensorErrRatios={8,8}` (entro tolleranza), quindi scaling `~L^-3`.
8. `check=True`.

Lettura fisica:
1. Questo blocco fornisce la forma "general case" del kernel quadratico a livello operatoriale: tensor, vector e scalar sono separati dalla simmetria isotropa, con mixing confinato al blocco scalare.
2. Il termine non-locale entra come deformazione IR-soppressa dei coefficienti del blocco; nel limite di dominio causale grande si recupera il kernel locale.
3. Non e' ancora la derivazione completa dei coefficienti dall'azione 4D totale, ma e' il ponte tecnico necessario per chiudere il passo finale senza tornare al minisuperspace.

## 08 - Covariant to Spin-Coefficients Map

Riferimento:
1. Script: `scripts/08_covariant_to_spin_coefficients_map.wl`.
2. Log: `logs/08_covariant_to_spin_coefficients_map.log`.

Calcolo:
1. Si parte dalla forma covariante su branch:
   `q2 = (d2I1 - q d2I0)/I0 - 2 dI0 (dI1 - q dI0)/I0^2`.
2. Si parametrizzano i blocchi lineari e quadratici per canali:
   tensor (`ht`), vector (`hv`), scalari (`xs`,`xw`).
3. Si espande `q2` e si estraggono i coefficienti:
   `c2`, `c1`, `cs`, `cw`, `cm`.
4. Si verifica la forma analitica attesa coefficiente per coefficiente.
5. Si verifica il limite leading `I0*c -> B-qA`.
6. Si controlla numericamente stabilita del blocco scalare e scaling IR.

Output:
1. `checkA..checkE=True` (formule chiuse di `c2,c1,cs,cw,cm`).
2. `checkF=True` (`detKscalar = cs*cw-cm^2`).
3. `checkG..checkK=True` (limiti leading `I0*c`).
4. Esempio numerico:
   `c2Vals,c1Vals,csVals,cwVals,cmVals` tutti decrescenti con `L`.
5. `ratioC2, ratioC1, ratioCs, ratioCm ~ {8,8}` per `L={10,20,40}`.
6. `eigs(Kscalar)>0` nel campione scelto.
7. `check=True`.

Lettura fisica:
1. Questa e' la mappa operativa mancante tra variazione covariante e coefficienti di spin del kernel quadratico.
2. Rende esplicito come i contributi non-locali (`a,b,A,B`) modificano canali tensor/vector/scalar.
3. Conferma che l'effetto resta IR-dominato (`~1/I0`) con subleading `~1/I0^2`.

## 09 - Vector Sector Full Stability Map

Riferimento:
1. Script: `scripts/09_vector_sector_full_stability_map.wl`.
2. Log: `logs/09_vector_sector_full_stability_map.log`.

Calcolo:
1. Si definiscono:
   `cT2=g2/c2`, `cV2=g1/c1`,
   `Ks={{cs,cm},{cm,cw}}`, `Gs={{gs,gm},{gm,gw}}`.
2. Si ricavano forme chiuse:
   `trace = Tr[Ks^-1 Gs]`, `det = Det[Ks^-1 Gs]`,
   autovalori scalari `sPlus,sMinus`.
3. Si verificano:
   regole somma/prodotto degli autovalori;
   limite di disaccoppiamento `cm=gm=0`;
   separazione dei canali tensor/vector dal blocco scalare.
4. Si inseriscono correzioni IR:
   `coeff = coeff0 + ncoeff/I0`.
5. Si testa convergenza a limite locale e stabilita numerica.

Output:
1. `checkA..checkH=True` (determinanti/tracce/formule autovalori/disaccoppiamento/separazione canali).
2. `checkI=True`, `checkJ=True`, `checkK=True` (limiti `I0->Infinity` corretti).
3. Valori numerici su `L={10,20,40}`:
   `cT2Vals`, `cV2Vals`, `scalarVals` tutti positivi nel campione.
4. `detKsVals>0`, `detGsVals>0`.
5. Errori verso limite locale con rapporti:
   `ratioT ~ 8`, `ratioV ~ 8`, `ratioS ~ 8`.
6. `check=True`.

Lettura fisica:
1. Il settore vettoriale entra esplicitamente nella stessa infrastruttura di stabilita dei canali tensoriale/scalare.
2. La gerarchia IR e' uniforme: tutti i canali convergono al limite locale con correzioni `~1/I0`.
3. Questo chiude il "full stability map" al livello di kernel parametrico completo.

## 10 - First-Principles abAB Map (Covariant Expansion)

Riferimento:
1. Script: `scripts/10_first_principles_abAB_map.wl`.
2. Log: `logs/10_first_principles_abAB_map.log`.

Calcolo:
1. Si espandono in `eps`:
   `muEps = mu + eps mu1 + eps^2 mu2/2`,
   `chiEps = chi + eps chi1 + eps^2 chi2/2`,
   `kEps = k0 + eps k1 + eps^2 k2/2`.
2. Si definiscono:
   `j0 = chi mu`, `j1 = chi mu K2`, e si estraggono `dj0,d2j0,dj1,d2j1`.
3. Si proietta su base di campi `phi={ht,hv,xs,xw}`:
   `dI0 = a.phi`, `dI1 = b.phi`,
   `d2I0 = phi.A.phi`, `d2I1 = phi.B.phi`.
4. Si ricavano `a,b,A,B` via coefficienti/Hessiano e si confrontano con le forme attese da regola di prodotto.
5. Si costruisce il kernel quadratico:
   `qKernel = (B-qA)/I0 - (a(b-q a)^T + (b-q a)a^T)/I0^2`.
6. Si verifica che coincide con il kernel ottenuto direttamente da
   `q2 = (d2I1-q d2I0)/I0 - 2 dI0(dI1-q dI0)/I0^2`.

Output:
1. `aVec = {chi*m1 + c1*mu, ...}`, `checkA=True`.
2. `bVec = {chi*k0*m1 + c1*k0*mu + chi*mu*p1, ...}`, `checkB=True`.
3. `AKernel`/`BKernel` uguali alle forme attese, `checkC=True`, `checkD=True`.
4. `checkE=True` (`qKernel` da `A,B,a,b` coincide con quello da espansione diretta).
5. `qKernelLead = Limit[I0*qKernel, I0->Infinity] = B-qA`, `checkF=True`.
6. Test numerico:
   `normVals = {1.3636e-4, 1.7048e-5, 2.1310e-6}`,
   `ratios ~ {7.9989, 7.9999}`,
   `leadErr` decresce con `leadRatios = {8,8}`.
7. `check=True`.

Lettura fisica:
1. Il blocco intermedio `{a,b,A,B}` e' ora derivato direttamente dalla struttura covariante non-locale, non postulado.
2. La gerarchia IR e' consistente: `qKernel ~ 1/I0`, con termine leading controllato da `B-qA`.
3. Questo riduce in modo sostanziale il gap tecnico segnalato dal referee su "parametrizzazione non derivata".

## 11 - Explicit Metric+Theta Mode Dictionary

Riferimento:
1. Script: `scripts/11_explicit_metric_theta_mode_dictionary.wl`.
2. Log: `logs/11_explicit_metric_theta_mode_dictionary.log`.

Calcolo:
1. Si imposta un ansatz ADM esplicito:
   `gamma_ij = a^2 exp(2 eps xs) diag(exp(2 eps ht), exp(-eps ht), exp(-eps ht))`,
   `K_ij = H0 gamma_ij + eps hv S_ij`,
   `ThetaDot = 1 + eps xw`.
2. Si definiscono:
   `mu = sqrt(det gamma) sqrt(X)` con `X = ThetaDot^2`,
   `chiProxy = exp(3 eps(-xs + xw))`,
   `K2 = gamma^-1 gamma^-1 K K`.
3. Si estraggono i dizionari:
   `m,m2` da `mu`, `c,c2` da `chiProxy`, `p,k2` da `K2`.
4. Si verifica ricostruzione fino a ordine `eps^2` e residuo `O(eps^3)`.

Output:
1. `mVec = {0,0,3*sqrt(a^6),sqrt(a^6)}`, `checkD=True`.
2. `cVec = {0,0,-3,3}`, `checkE=True`.
3. `k0 = 3 H0^2`, `checkF=True`.
4. `pVec = {0,0,0,0}`.
5. `K2` ha unico canale quadratico non nullo:
   `K2_22 = 4 sv^2 / a^4`.
6. Ricostruzioni: `checkA=True`, `checkB=True`, `checkC=True`.
7. Residui numerici:
   `muRat ~ {8.0058, 8.0029}`,
   `chiRat ~ {7.9993, 7.9996}`,
   `kRat ~ {7.9817, 7.9908}`,
   con tutti i check numerici `True`.
8. `check=True`.

Lettura fisica:
1. I tensori elementari (`m,c,p,m2,c2,k2`) sono ora estratti da un modello geometrico esplicito, non solo introdotti come simboli.
2. Nel setup scelto, il canale vettoriale entra in `K2` al secondo ordine (`hv^2`) senza termine lineare, in coerenza con la struttura di stabilita osservata nei notebook precedenti.
3. La coerenza di ordine (`O(eps^3)` nei residui) conferma la consistenza del dizionario perturbativo usato in `10`.

## 12 - Full 4D Spin-Projector Inversion

Riferimento:
1. Script: `scripts/12_full4d_spin_projector_inversion.wl`.
2. Log: `logs/12_full4d_spin_projector_inversion.log`.

Calcolo:
1. Si costruiscono i proiettori Barnes-Rivers full 4D:
   `P2, P1, P0s, P0w, P0sw, P0ws` con `d-1=3`.
2. Si verifica l'algebra:
   idempotenza, ortogonalita e completezza su `Isym`.
3. Si definisce:
   `Op4D = c2 P2 + c1 P1 + cs P0s + cw P0w + cm(P0sw+P0ws)`.
4. Si costruisce l'inversa chiusa:
   `Op4DInv = P2/c2 + P1/c1 + (cw P0s + cs P0w - cm(P0sw+P0ws))/(cs*cw-cm^2)`.
5. Si verifica esplicitamente:
   `Op4D.Op4DInv = Isym` e `Op4DInv.Op4D = Isym`.
6. Si testa il flusso IR con:
   `c* = c*0 + nc*/I0`, e con i residui canale per canale.

Output:
1. `checkA..checkJ=True` (algebra proiettori + completezza).
2. `detS = cs*cw-cm^2`.
3. `checkK=True`, `checkL=True` (inversione bilatera verificata).
4. `checkM=True` (assenza leakage tensor-scalare spurio).
5. `checkN=True` (`Det[resSIR]=1/detSIR`).
6. Limiti IR:
   `c2Inf=c20`, `c1Inf=c10`, `ksInf={{cs0,cm0},{cm0,cw0}}`, `resTInf=1/c20`, `resVInf=1/c10`,
   con `checkO=True`, `checkP=True`, `checkQ=True`.
7. Test numerici (`L={10,20,40}`):
   coefficienti e residui positivi; rapporti `ratioC2, ratioC1, ratioDet, ratioResT, ratioResV ~ {8,8}`.
8. `check=True`.

Lettura fisica:
1. Il kernel perturbativo 4D non e' solo decomposto: e' anche invertibile in forma chiusa, con condizioni trasparenti (`c2!=0`, `c1!=0`, `detS!=0`).
2. I residui effettivi dei canali mantengono la gerarchia IR (`~1/I0`) senza introdurre mixing patologico tensor-scalare.

## 13 - Full Spectrum Dispersion and Residues

Riferimento:
1. Script: `scripts/13_full_spectrum_dispersion_residues.wl`.
2. Log: `logs/13_full_spectrum_dispersion_residues.log`.

Calcolo:
1. Settore scalare misto:
   `det[x Ks - (q2 Gs + Ms)] = 0`.
2. Si estraggono coefficienti del polinomio:
   `a2, a1, a0` e si verificano contro le forme attese.
3. Si ottengono le radici `xPlus,xMinus` e si verificano regole somma/prodotto.
4. Si testa il limite di disaccoppiamento (`km=gm=mm=0`).
5. Si includono anche canali:
   `xTensor=(q2 g2 + mt2)/c2`, `xVector=(q2 g1 + mv2)/c1`.
6. Si verifica il ponte UV scalare:
   radici ad alto `q2` / `q2` contro `Eigenvalues[K^-1 G]`.
7. Si applica flusso IR:
   `coeff = coeff0 + ncoeff/I0` e si controlla scaling numerico.

Output:
1. `checkA=True`, `checkB=True`, `checkC=True` (`a2,a1,a0` corretti).
2. `checkD=True`, `checkE=True` (somma/prodotto radici).
3. `checkF=True` (limite di disaccoppiamento).
4. `checkG=True` (separazione tensor/vector dal blocco scalare).
5. `checkH=True` (limite UV scalare da `K^-1 G`).
6. `uvErr ~ 5e-9`.
7. Positivita numerica dei canali su `L={10,20,40}`, `q2={0.1,1,4}`.
8. Scaling IR:
   `ratioScalar ~ {7.998, 7.9998}`,
   `ratioTensor ~ {7.9983, 7.9998}`,
   `ratioVector ~ {8.0015, 8.0002}`.
9. `check=True`.

Lettura fisica:
1. Il quadro ghost/gradient non resta limitato a condizioni statiche su coefficienti: si traduce in dispersioni esplicite dei modi propaganti.
2. Il limite UV recupera correttamente le velocita da `K^-1 G`, mentre il settore IR resta controllato e coerente con la soppressione non-locale del paper.

## 14 - Passaggio 4D generale (metrica+Theta) con chiusura gauge-fixing parametrica

Riferimento:
1. Script: `scripts/14_full4d_covariant_gauge_fixed_coefficients.wl`.
2. Log: `logs/14_full4d_covariant_gauge_fixed_coefficients.log`.
3. Blocco teorico sorgente nel paper: `S_loc + integral dTheta lambda(Theta)[I1-K0^2 I0]` (`Paper/main.tex`, eq. `eq:lag_mult`) e definizioni leafwise (`I0`,`I1`).

Calcolo:
1. Si parte dalla forma 4D generale:
   `S = integral d4x sqrt(-g)[(R-2Lambda)/(2kappa) - (alpha/2) X] + integral dtheta lambda(theta)[I1[theta]-K0^2 I0[theta]]`,
   con `X = g^(mu nu) partial_mu Theta partial_nu Theta`.
2. Si usa la rappresentazione covariante dei funzionali leafwise:
   `I0[theta] = integral d4x mu chi_theta`,
   `I1[theta] = integral d4x mu chi_theta K^2`,
   `mu = sqrt(-g) W sqrt(X) DiracDelta[Theta-theta]`.
3. Si espande attorno a background generico (`Xbar>0`, nessun ansatz simmetrico):
   `g_mu nu = gbar_mu nu + eps h_mu nu`,
   `Theta = Thetabar + eps vartheta`,
   `lambda = lambdabar + eps ell`.
4. Si introduce gauge-fixing pienamente covariante al secondo ordine:
   `S_gf = -(1/(2 xi_g)) integral d4x sqrt(-gbar) F_mu F^mu - (1/(2 xi_Theta)) integral d4x sqrt(-gbar) G^2`,
   con scelta minimale
   `F_mu = nablaBar^nu h_mu nu - (1/2) nablaBar_mu h`,
   `G = ubar^mu nablaBar_mu vartheta + nu h`,
   `ubar_mu = partial_mu Thetabar/sqrt(Xbar)`.
5. Dalla variazione covariante 4D (script `23`/`25`) si ottiene la decomposizione bulk+boundary:
   `deltaI0 = Int4[chi deltaMu] + Int4[mu deltaChi]`,
   `deltaI1 = Int4[chi K^2 deltaMu] + Int4[chi mu deltaK2] + Int4[mu K^2 deltaChi]`.
6. Sul branch `q = I1/I0` (on-shell del vincolo: `q=K0^2`) la variazione seconda del rapporto e':
   `delta2Q = (delta2I1 - q delta2I0)/I0 - 2 deltaI0 (deltaI1 - q deltaI0)/I0^2`.
7. Proiettando su base di modi `phi={ht,hv,xs,xw}`:
   `deltaI0 = a.phi`, `deltaI1 = b.phi`,
   `delta2I0 = phi.A.phi`, `delta2I1 = phi.B.phi`,
   la parte non-locale dei coefficienti e':
   `c2NL = (B2-q A2)/I0 - 2 at(bt-q at)/I0^2`,
   `c1NL = (B1-q A1)/I0 - 2 av(bv-q av)/I0^2`,
   `csNL = (Bs-q As)/I0 - 2 as(bs-q as)/I0^2`,
   `cwNL = (Bw-q Aw)/I0 - 2 aw(bw-q aw)/I0^2`,
   `cmNL = (Bm-q Am)/I0 - (as(bw-q aw)+aw(bs-q as))/I0^2`.
8. I coefficienti totali metrici/`Theta` sono implementati in forma operativa come:
   `cA = cA_loc + cA_NL + Delta cA_gf`, con `A in {2,1,s,w,m}`,
   dove nel notebook `Delta cA_gf` e' rappresentato tramite shift parametrici (`g2,g1,gs,gw,gm`), lasciando aperta la derivazione esplicita completa da `F_mu,G`.
9. L'operatore quadratico 4D mantiene la forma:
   `Op4D = c2 P2 + c1 P1 + cs P0s + cw P0w + cm(P0sw+P0ws)`,
   con inversa chiusa (se `c2!=0`, `c1!=0`, `cs*cw-cm^2!=0`):
   `Op4DInv = P2/c2 + P1/c1 + (cw P0s + cs P0w - cm(P0sw+P0ws))/(cs*cw-cm^2)`.
10. Nel limite IR (`I0->Infinity`) la parte non-locale decade (`~1/I0`) e si recupera il kernel locale gauge-fissato del background scelto.

Output:
1. Coefficienti non-locali estratti in forma chiusa:
   `c2NL,c1NL,csNL,cwNL,cmNL` come in `q2Expr`,
   con check simbolici `checkA..checkE=True`.
2. Verifica del branch vincolato `q -> K0^2`:
   `checkF..checkJ=True`.
3. Coefficienti totali gauge-fissati:
   `c2Tot = l2 + g2 + c2NL`, `c1Tot = l1 + g1 + c1NL`,
   `csTot = ls + gs + csNL`, `cwTot = lw + gw + cwNL`, `cmTot = lm + gm + cmNL`,
   con `checkK=True`.
4. Blocco scalare e residui:
   `detSTot = csTot*cwTot-cmTot^2`, `Det[resS]=1/detSTot`,
   con `checkL=True`, `checkM=True`.
5. Limiti IR:
   `c2Lead=B2-A2 q`, `c1Lead=B1-A1 q`, `csLead=Bs-As q`, `cwLead=Bw-Aw q`, `cmLead=Bm-Am q`,
   con `checkN..checkR=True`;
   e per il totale:
   `c2TotInf=l2+g2`, `c1TotInf=l1+g1`, `ksTotInf={{ls+gs,lm+gm},{lm+gm,lw+gw}}`,
   con `checkS=True`, `checkT=True`, `checkU=True`.
6. Test numerici (`L={10,20,40}`):
   `ratioC2NL,ratioC1NL,ratioCsNL,ratioCwNL,ratioCmNL ~ {8,8}`,
   `ratioC2Err,ratioC1Err,ratioDetErr ~ {8,8}`,
   `ksEigVals>0`, `decompErr={0,0,0}`.
7. `checkNumA=True`, `checkNumB=True`, `checkNumC=True`, `checkNumD=True`, `check=True`.

Lettura fisica:
1. Il "gap referee" sul passaggio diretto da azione 4D generale a coefficienti perturbativi e' ridotto in modo sostanziale a livello analitico-operativo, con chiusura attuale parametrica/proxy del blocco gauge-fixing.
2. Il gauge-fixing covariante rende invertibile il kernel senza alterare la gerarchia IR indotta dal termine leafwise non-locale.
3. Il contenuto del paper resta rispettato: la derivazione esplicita completa di `Delta c_gf` da `F_mu,G` resta un'estensione aperta, coerente con lo scope del testo principale.

## 15 - EFT Interpretation: Cutoff, Strong-Coupling, Validity Regime

Riferimento:
1. Script: `scripts/15_eft_validity_scales_proxy.wl`.
2. Log: `logs/15_eft_validity_scales_proxy.log`.
3. Blocchi usati: coefficienti gauge-fissati da `14`, dispersioni da struttura `13`.

Calcolo:
1. Si definiscono indicatori EFT di piccolezza delle correzioni non-locali:
   `eps2 = |c2NL/(l2+dg2)|`,
   `eps1 = |c1NL/(l1+dg1)|`,
   `epsS = ||KNL||2/||KIR||2`.
2. Si testa lo scaling IR su `L={10,20,40}` (con `I0=L^3`) e si controlla la legge `~1/L^3`.
3. Si stima una soglia di validita in lunghezza IR:
   `Lmin(eps<=1e-2)` e `Lmin(eps<=1e-3)` da fit `eps ~ C/L^3`.
4. Si introducono proxy di strong coupling dai canali cinetici:
   `zT=c2Tot`, `zV=c1Tot`, `zSmin=min eig(KsTot)`,
   e proxy di residuo:
   `invT=1/zT`, `invV=1/zV`, `invSmax=max eig(KsTot^-1)`.
5. Si definisce una finestra conservativa di cutoff operativo nel settore testato:
   `q2Grid={0.1,1,4}`, `q2CutProxy=max(q2Grid)=4`,
   con controllo di positivita di `xTensor`, `xVector`, `xScalarPlus`, `xScalarMinus` su tutta la griglia.

Output:
1. Piccolezza non-locale:
   `eps2Vals={4.7616e-4,5.9523e-5,7.4405e-6}`,
   `eps1Vals={5.9299e-4,7.4127e-5,9.2660e-6}`,
   `epsSVals={5.3000e-4,6.6256e-5,8.2821e-6}`.
2. Scaling:
   `ratioEps2~{7.9996,7.9999}`,
   `ratioEps1~{7.9996,7.9999}`,
   `ratioEpsS~{7.9992,7.9999}`,
   `checkA=True`, `checkB=True`.
3. Stima regime:
   `Lmin(eps=1e-2)=3.8996...`,
   `Lmin(eps=1e-3)=8.4014...`.
4. Strong-coupling proxy:
   `zTVals~{1.0505,1.05006,1.05001}`,
   `zVVals~{0.86051,0.86006,0.86001}`,
   `zSMinVals~{1.17012,1.16993,1.16991}`,
   `invTVals,invVVals,invSMaxVals` tutti `O(1)`,
   con `checkC=True`, `checkD=True`.
5. Cutoff operativo su finestra testata:
   tutti i modi di dispersione sono positivi su `q2={0.1,1,4}`,
   `checkE=True`, `q2CutProxy=4`, `checkF=True`.
6. `check=True`.

Lettura fisica:
1. Regime di validita perturbativo: nel campione validato (`L>=10`) le correzioni non-locali sono piccole (`~10^-4`-`10^-6`) rispetto ai blocchi IR locali.
2. Nel dominio testato non emerge segnale di strong coupling nei proxy adottati (`zFloor=0.1`, `invCeil=10`): i coefficienti cinetici restano lontani da zero e i residui restano finiti.
3. `q2CutProxy=4` va letto come massimo `q^2` testato con positivita delle dispersioni nella griglia corrente, non come cutoff UV fondamentale della teoria.
4. Il completamento UV resta aperto (coerente con il paper): qui si fornisce una diagnostica EFT interna e riproducibile del regime affidabile della pipeline.

## Quadro complessivo

Risultato tecnico:
1. Tutti gli script `00..15` hanno `check=True`.
2. I criteri tensoriali, di mixing scalare e di riduzione del modo non dinamico sono espliciti e verificati.
3. La variazione seconda covariante di `Q=I1/I0` e il kernel quadratico a due campi sono derivati in forma chiusa.
4. I termini di bordo fino al secondo ordine sono verificati analiticamente e numericamente.
5. La decomposizione operatoriale generale in blocchi di spin e' verificata con algebra dei proiettori e limiti IR controllati.
6. La mappa covariante -> coefficienti di spin e la stabilita completa con settore vettoriale sono ora esplicite.
7. La mappa first-principles `eps-expansion -> {a,b,A,B} -> qKernel` e' verificata simbolicamente.
8. Il dizionario `metric+Theta -> {m,c,p,m2,c2,k2}` e' verificato in ansatz ADM esplicito.
9. Il kernel in proiettori full 4D e' invertito esplicitamente con controlli di residuo e limiti IR.
10. Le dispersioni complete tensor/vector/scalar (con mixing e massa) sono verificate con ponte UV e scaling IR.
11. Il passaggio 4D generale ai coefficienti metrici/`Theta` e' ora chiuso in forma parametrica/proxy compatibile con gauge-fixing covariante (sezione `14`), mentre la derivazione esplicita completa di `Delta c_gf` resta aperta.
12. L'interpretazione EFT (cutoff operativo, proxy strong-coupling, e regime di validita) e' ora esplicitata e verificata (sezione `15`).

Risultato fisico:
1. Tensor mode: condizioni ghost/gradient chiare (`GT>0`, `FT>0`) e controllo IR.
2. Vector mode: velocita `cV^2=g1/c1` trattata esplicitamente con condizioni di positivita.
3. Mixing scalare: stabilita letta dagli autovalori di `K^-1 G`.
4. Canale scalare FRW: conferma del congelamento nel branch/gauge rilevante.
5. Settore non-locale quadratico: contributi IR-soppressi e dipendenza dal dominio causale resa esplicita.
6. A livello operatoriale generale, il mixing non-locale non forza accoppiamenti tensor-scalari spurii in isotropia.
7. Tutti i canali considerati convergono al limite locale con legge IR coerente.
8. Le frequenze fisiche dei modi (tensor/vector/scalar) restano positive nel campione stabile, con transizione UV corretta verso le velocita `K^-1 G`.
9. Il gauge-fixing covariante modifica i canali puramente di gauge ma lascia invariata la struttura fisica IR/UV dei modi propaganti.
10. Nel dominio testato la teoria funziona come EFT controllata: correzioni non-locali piccole, nessun segnale di strong coupling nei proxy adottati, e finestra conservativa `q2<=4` intesa come massimo `q^2` testato con dispersioni positive.

Limite attuale:
1. La discussione EFT qui introdotta usa cutoff/proxy operativi interni alla pipeline; la determinazione di un cutoff UV fondamentale richiede un completamento UV esplicito.
2. Resta da estendere la chiusura BRST/Faddeev-Popov esplicita di indipendenza dal gauge.
3. Resta aperto il passaggio finale a predizioni osservative complete (spettri/transfer) in setup cosmologico realistico.

## Riproduzione rapida

```powershell
cd c:\lean\StaticSanityCheck\Paper\Perturbative_Spectrum_Study
.\run_all.ps1
```
