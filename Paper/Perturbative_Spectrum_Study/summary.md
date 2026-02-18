# Perturbative Spectrum Summary

Mappa notebook `Paper/Perturbative_Spectrum_Study` -> output verificati (`check=True`) e lettura fisica.

## Notebook 00-15

1. `00_Tensor_Channel_IR_Stability.nb`
`Riferimento`: canale tensoriale su background FRW/de Sitter in presenza di correzioni IR-soppresse, coerente con `S_loc` (`Paper/main.tex:209-218`) e dominanza IR (`Paper/main.tex:321-329`).
`Output`: `cT2=(1+dFT*kappa)/(1+dGT*kappa)`, condizioni `GT>0`, `FT>0`, forma IR `cT2IR=(I0+aFT*kappa)/(I0+aGT*kappa)`, `deltaAsym=(aFT-aGT)kappa`.
`Interpretazione`: il tensore resta ghost/gradient-stable per `GT,FT>0`; le deviazioni da GR vanno a zero con crescita del dominio causale.

2. `01_Scalar_Mixing_Stability_Matrix.nb`
`Riferimento`: matrice scalare generale a due campi.
`Output`: `detK=-k12^2+k11*k22`, `trace(K^-1G)`, `det(K^-1G)`, autovalori `c2plus,c2minus`; benchmark stabile (entrambi positivi) e instabile (uno negativo).
`Interpretazione`: la stabilita UV si riduce a positivita di `K` e autovalori di `K^-1G`.

3. `02_Frozen_Scalar_Reduction.nb`
`Riferimento`: integrazione di un modo scalare non dinamico.
`Output`: `piRule={pi->(-(E q2 z)+B zd)/(C+F q2)}`, `Keff=A+B^2/(C+F q2)`, `Geff=D-(E^2 q2)/(C+F q2)`, limiti UV e freeze hard coerenti.
`Interpretazione`: il modo ausiliario si elimina in modo controllato; il criterio gradient UV diventa `D-E^2/F>0`.

4. `03_FRW_Scalar_Constraint_Freeze_Map.nb`
`Riferimento`: mappa del vincolo lineare FRW gia ottenuto nella parte EOM.
`Output`: `n1=(Sqrt[3]/K0)up`, gauge `n1=0 => up=0`, `LkinGauge=0`.
`Interpretazione`: il canale scalare FRW del branch vincolato risulta congelato nel gauge cosmico.

5. `04_Covariant_Ratio_Second_Variation.nb`
`Riferimento`: variazione seconda covariante di `Q=I1/I0`, ponte diretto verso l'azione quadratica non-locale.
`Output`:  
`q1=(di1*i0-di0*i1)/i0^2`,  
`q2=(-2*di0*di1*i0+2*di0^2*i1+i0*(d2i1*i0-d2i0*i1))/i0^3`,  
forma su branch `i1=q i0`:  
`q1Branch=(di1-di0*q)/i0`,  
`q2Branch=(-2*di0*di1+2*di0^2*q+i0*(d2i1-d2i0*q))/i0^2`.  
Test numerico con `i0=L^3`: rapporti `~8` tra `L` e `2L` per `L={10,20,40}`.
`Interpretazione`: primo e secondo canale variazionale sono entrambi IR-soppressi (`~1/L^3` nel setup testato).

6. `05_Nonlocal_Quadratic_Mixing_Kernel.nb`
`Riferimento`: kernel quadratico a due campi generato da `q2`.
`Output`: coefficienti espliciti `Mxx,Mxy,Myy`; limite `Klead=Limit[i0 Kmix, i0->Infinity]` uguale alla matrice dei coefficienti di secondo ordine `(aa1-q aa0, ab1-q ab0, bb1-q bb0)`; autovalori numerici positivi e decrescenti con `L`.
`Interpretazione`: il mixing non-locale quadratico resta controllato in IR; la parte leading e' ben definita e stabile nel campione scelto.

7. `06_Second_Order_Boundary_Shape_Derivative.nb`
`Riferimento`: variazione del dominio causale mobile al secondo ordine.
`Output`:  
`dchi1=rho DiracDelta[r-R]`,  
`dchi2=sigma DiracDelta[r-R]+4 rho^2 Derivative[1][DiracDelta][-2 r+2 R]`;  
equivalenza Leibniz/distribuzionale per derivata prima e seconda; errore numerico `~10^-12`.
`Interpretazione`: i termini di bordo al secondo ordine non sono assunti ma derivati; la forma distribuzionale e' verificata analiticamente e numericamente.

8. `07_Spin_Projector_Quadratic_Operator.nb`
`Riferimento`: decomposizione generale del kernel quadratico isotropo in proiettori di Barnes-Rivers (settori spin-2, spin-1, spin-0), come proxy tecnico del "caso generale" oltre minisuperspace.
`Output`: algebra dei proiettori verificata (`P2^2=P2`, `P1^2=P1`, completezza `P2+P1+P0s+P0w=Isym`), assenza di mixing tensor-scalare in operatore isotropo (`P2 Op P0*=0`), blocco scalare `Kscalar={{cs,cm},{cm,cw}}` con `det=-cm^2+cs*cw`; con coefficenti IR `cs,cw,cm ~ c0 + n/I0`, limiti `I0->Infinity` corretti e scaling numerico `~1/L^3`.
`Interpretazione`: il canale tensoriale resta separato e il mixing non-locale entra naturalmente nel blocco scalare 2x2; il limite IR recupera il kernel locale di riferimento.

9. `08_Covariant_To_Spin_Coefficients_Map.nb`
`Riferimento`: identificazione esplicita dei coefficienti di spin a partire dalla formula covariante di `delta^2Q` su branch.
`Output`:  
`c2=(-2 at bt + 2 at^2 q + I0(B2-A2 q))/I0^2`,  
`c1=(-2 av bv + 2 av^2 q + I0(B1-A1 q))/I0^2`,  
`cs=(-2 as bs + 2 as^2 q + I0(Bs-As q))/I0^2`,  
`cw=(-2 aw bw + 2 aw^2 q + I0(Bw-Aw q))/I0^2`,  
`cm=(Bm-Am q)/I0 - (as(bw-q aw)+aw(bs-q as))/I0^2`;  
limiti `I0*c* -> {B*-q A*}` verificati; test numerici positivi e scaling `~1/L^3`.
`Interpretazione`: la dipendenza dei coefficienti dinamici dai blocchi variazionali covarianti e' ora esplicita, non piu solo parametrica.

10. `09_Vector_Sector_Full_Stability_Map.nb`
`Riferimento`: mappa di stabilita completa con canali tensoriale, vettoriale e scalare.
`Output`: `cT2=g2/c2`, `cV2=g1/c1`, settore scalare da autovalori di `Ks^-1 Gs`; formule chiuse di `trace`/`det` verificate, limiti di disaccoppiamento (`cm=gm=0`) verificati, canali separati, convergenza IR ai valori locali con rapporti errori `~8`.
`Interpretazione`: il settore vettoriale e' trattato allo stesso livello degli altri canali; la stabilita completa puo' essere letta da condizioni algebriche uniformi.

11. `10_First_Principles_abAB_Map.nb`
`Riferimento`: derivazione first-principles dei blocchi intermedi `{a,b,A,B}` dalla espansione covariante in `eps` di `j0=chi mu` e `j1=chi mu K2`.
`Output`:  
`a = chi m + mu c`,  
`b = chi k0 m + mu k0 c + chi mu p`,  
`A = chi M2 + mu C2 + m c^T + c m^T`,  
`B = chi k0 M2 + mu k0 C2 + chi mu K2 + k0(mc^T+cm^T) + mu(cp^T+pc^T) + chi(mp^T+pm^T)`;  
identita `qKernel(AB)=qKernel(direttamente da q2)` verificata; limiti IR `I0 qKernel -> B-qA`; scaling numerico `~L^-3` con rapporti `~8`.
`Interpretazione`: il collegamento tra variazione covariante e coefficienti intermedi non e' piu assunto: e' derivato simbolicamente dalla struttura dell'azione non-locale.

12. `11_Explicit_Metric_Theta_Mode_Dictionary.nb`
`Riferimento`: dizionario esplicito `metric+Theta -> {m,c,p,m2,c2,k2}` da ansatz ADM.
`Output`:  
`mVec={0,0,3*sqrt(a^6),sqrt(a^6)}`,  
`cVec={0,0,-3,3}`,  
`k0=3 H0^2`,  
`pVec={0,0,0,0}`,  
`K2` con unico canale quadratico non nullo `hv^2` (`K2_22=4 sv^2/a^4`);  
ricostruzioni ordine-2 (`mu`, `chi`, `K2`) verificate; residui numerici di troncamento con ratio `~8`.
`Interpretazione`: i tensori elementari non sono piu solo oggetti astratti; in un ansatz geometrico esplicito si estraggono direttamente e mostrano una gerarchia IR coerente.

13. `12_Full4D_Spin_Projector_Inversion.nb`
`Riferimento`: estensione del blocco a proiettori Barnes-Rivers in full 4D con inversione esplicita dell'operatore quadratico.
`Output`: algebra dei proiettori `P2,P1,P0s,P0w,P0sw,P0ws` verificata (`checkA..checkJ=True`), completezza `P2+P1+P0s+P0w=Isym`, inversa chiusa `OpInv` verificata su entrambi i lati (`checkK=True`, `checkL=True`), assenza leakage tensor-scalare (`checkM=True`), limiti IR corretti, scaling numerico `~1/L^3` per coefficienti e residui (`ratio* ~ 8`).
`Interpretazione`: il controllo non e' piu solo su decomposizione parametrica; ora l'inversione del kernel e i residui dei canali sono consistenti nel formalismo 4D.

14. `13_Full_Spectrum_Dispersion_Residues.nb`
`Riferimento`: mappa completa delle dispersioni in canali tensoriale/vettoriale/scalare con mixing e termini di massa.
`Output`: polinomio scalare `det[x Ks-(q2 Gs+Ms)]` con coefficienti chiusi (`checkA..checkC=True`), regole somma/prodotto radici (`checkD`, `checkE`), limiti di disaccoppiamento (`checkF=True`), separazione tensor/vector (`checkG=True`), limite UV scalare coerente con `Eigenvalues[K^-1 G]` (`checkH=True`, `uvErr ~ 5e-9`), positivita dei canali numerici e scaling IR `~1/L^3` (`ratioScalar`, `ratioTensor`, `ratioVector ~ 8`).
`Interpretazione`: il quadro ghost/gradient e' ora chiuso al livello di spettro di dispersione completo, con verifica esplicita del ponte UV e della soppressione IR.

15. `14_Full4D_Covariant_GaugeFixed_Coefficients.nb`
`Riferimento`: passaggio 4D generale ai coefficienti quadratici metrici/`Theta` da `q2`, con split locale + non-locale e chiusura gauge-fixing implementata in forma parametrica (`g2,g1,gs,gw,gm`).
`Output`: formule chiuse `c2NL,c1NL,csNL,cwNL,cmNL` (`checkA..checkE=True`), verifica branch `q->K0^2` (`checkF..checkJ=True`), costruzione coefficienti totali `cTot = l + g + cNL` (`checkK=True`), consistenza inversa del blocco scalare (`checkL=True`, `checkM=True`), limiti IR (`checkN..checkU=True`), scaling numerico `~1/L^3` con rapporti `~8` e `checkNumA..D=True`.
`Interpretazione`: il passaggio da variazione covariante 4D a coefficienti effettivi gauge-fissati e' implementato end-to-end in forma operativa; resta aperta la derivazione esplicita completa di `Delta c_gf` da `F_mu,G`.

16. `15_EFT_Validity_Scales_Proxy.nb`
`Riferimento`: interpretazione EFT operativa del blocco perturbativo con tre diagnostiche: piccolezza delle correzioni non-locali, distanza da strong coupling, finestra di cutoff conservativa nel dominio testato.
`Output`: `eps2Vals,eps1Vals,epsSVals` piccoli (`~10^-4`-`10^-6`) e con scaling `~L^-3` (`ratioEps* ~ 8`), stime `Lmin(eps=1e-2)=3.8996...`, `Lmin(eps=1e-3)=8.4014...`, proxy cinetici `zT,zV,zSmin` tutti `O(1)>0`, dispersioni positive su `q2={0.1,1,4}`, `q2CutProxy=4`, `checkA..F=True`, `check=True`.
`Interpretazione`: il framework funziona come EFT controllata nel dominio testato; `q2CutProxy=4` e' il massimo `q^2` testato con dispersioni positive (non un cutoff UV fondamentale) e l'assenza di strong coupling va letta nei proxy/soglie adottati.

## Quadro fisico attuale
Coperto:
1. criteri ghost/gradient tensor e scalar mixing.
2. collegamento con freeze scalare FRW della catena EOM.
3. estensione quadratica non-locale via `delta^2(I1/I0)` e boundary shape derivative di ordine 2.
4. decomposizione operatoriale generale in blocchi di spin con separazione tensor/scalar verificata.
5. mappa esplicita covariante -> coefficienti di spin (`08`).
6. stabilita completa con settore vettoriale esplicito (`09`).
7. derivazione first-principles dei blocchi intermedi `a,b,A,B` (`10`).
8. estrazione esplicita dei tensori elementari in ansatz ADM (`11`).
9. algebra+inversione del kernel in formalismo full 4D (`12`).
10. dispersioni complete tensor/vector/scalar con controllo UV e residui (`13`).
11. passaggio 4D generale ai coefficienti metrici/`Theta` con chiusura gauge-fixing parametrica/proxy (`14`), con derivazione esplicita completa di `Delta c_gf` ancora aperta.
12. discussione EFT esplicita di cutoff operativo, strong-coupling proxy e regime di validita (`15`).

Ancora aperto:
1. prova BRST/Faddeev-Popov di indipendenza dal gauge sui residui fisici nella pipeline numerica.
2. cutoff UV fondamentale e completamento UV della teoria oltre il cutoff operativo della pipeline.
3. estendere il dizionario del punto `11` da ansatz controllato a decomposizione 4D completamente generale senza restrizioni di simmetria.
4. passaggio finale a predizioni osservabili del settore perturbativo in setup cosmologico completo.
