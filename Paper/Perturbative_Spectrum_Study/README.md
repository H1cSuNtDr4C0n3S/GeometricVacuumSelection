# Perturbative_Spectrum_Study

Pipeline Mathematica/Wolfram dedicata al blocco "spettro perturbativo (ghost/gradient)" della teoria, in continuita con la parte EOM gia chiusa in `Paper/Variational_Derivation`.

Focus operativo:
1. canale tensoriale e stabilita IR;
2. mixing scalare generale (2 campi);
3. riduzione del canale scalare congelato in FRW;
4. estensione al livello quadratico del termine non-locale via variazioni seconde di `Q=I1/I0`;
5. termini di bordo di shape derivative fino al secondo ordine;
6. decomposizione generale a proiettori (spin-2/spin-1/spin-0) del kernel quadratico isotropo con mixing non-locale IR;
7. mappa esplicita dai blocchi variazionali covarianti ai coefficienti di spin `{c2,c1,cs,cw,cm}`;
8. chiusura del settore vettoriale e mappa di stabilita completa (`tensor + vector + scalar`);
9. derivazione first-principles dei blocchi intermedi `{a,b,A,B}` dalla espansione covariante di `j0=chi mu`, `j1=chi mu K2`;
10. estrazione esplicita del dizionario modi `metric+Theta -> {m,c,p,m2,c2,k2}` in un ansatz ADM verificato;
11. algebra e inversione del kernel in proiettori Barnes-Rivers full 4D;
12. mappa completa delle dispersioni tensor/vector/scalar (con mixing e masse) e controllo dei residui/limite UV;
13. passaggio 4D generale ai coefficienti metrici/`Theta` con chiusura gauge-fixing parametrica/proxy in script/notebook dedicato (`14`);
14. interpretazione EFT operativa: cutoff conservativo nel dominio testato, proxy strong-coupling, e regime di validita quantitativo (`15`).

## Attinenza al paper
Il paper marca la stabilita perturbativa completa e la formulazione fully non-local come lavoro futuro (`Paper/main.tex:102-107`, `Paper/main.tex:230-231`, `Paper/main.tex:315-317`).
Questa cartella resta aderente a quel perimetro: costruisce e verifica i blocchi necessari senza introdurre claim oltre il testo.

Allineamento esplicito con il paper:
1. settore locale `S_loc` (`Paper/main.tex:209-218`);
2. termine leafwise con moltiplicatore `lambda(Theta)` (`Paper/main.tex:292-301`);
3. dominanza IR e soppressione locale (`Paper/main.tex:321-329`);
4. definizione causale leafwise di dominio (`Paper/main.tex:242-258`).

## Struttura
- `scripts/00_tensor_channel_ir_stability.wl`
- `scripts/01_scalar_mixing_stability_matrix.wl`
- `scripts/02_frozen_scalar_reduction.wl`
- `scripts/03_frw_scalar_constraint_freeze_map.wl`
- `scripts/04_covariant_ratio_second_variation.wl`
- `scripts/05_nonlocal_quadratic_mixing_kernel.wl`
- `scripts/06_second_order_boundary_shape_derivative.wl`
- `scripts/07_spin_projector_quadratic_operator.wl`
- `scripts/08_covariant_to_spin_coefficients_map.wl`
- `scripts/09_vector_sector_full_stability_map.wl`
- `scripts/10_first_principles_abAB_map.wl`
- `scripts/11_explicit_metric_theta_mode_dictionary.wl`
- `scripts/12_full4d_spin_projector_inversion.wl`
- `scripts/13_full_spectrum_dispersion_residues.wl`
- `scripts/14_full4d_covariant_gauge_fixed_coefficients.wl`
- `scripts/15_eft_validity_scales_proxy.wl`
- `logs/`
- `summary.md`
- `detailed_calculations.md`

## Stato di verifica
Tutti gli script `00..15` risultano con `check=True` nei rispettivi log:
1. `logs/00_tensor_channel_ir_stability.log`
2. `logs/01_scalar_mixing_stability_matrix.log`
3. `logs/02_frozen_scalar_reduction.log`
4. `logs/03_frw_scalar_constraint_freeze_map.log`
5. `logs/04_covariant_ratio_second_variation.log`
6. `logs/05_nonlocal_quadratic_mixing_kernel.log`
7. `logs/06_second_order_boundary_shape_derivative.log`
8. `logs/07_spin_projector_quadratic_operator.log`
9. `logs/08_covariant_to_spin_coefficients_map.log`
10. `logs/09_vector_sector_full_stability_map.log`
11. `logs/10_first_principles_abAB_map.log`
12. `logs/11_explicit_metric_theta_mode_dictionary.log`
13. `logs/12_full4d_spin_projector_inversion.log`
14. `logs/13_full_spectrum_dispersion_residues.log`
15. `logs/14_full4d_covariant_gauge_fixed_coefficients.log`
16. `logs/15_eft_validity_scales_proxy.log`

## Cosa copre gia la richiesta del "referee cattivo"
Coperto:
1. criteri ghost/gradient tensor/scalar in forma chiusa;
2. struttura quadratica non-locale derivata da variazioni di rapporto covarianti;
3. termini di bordo causali trattati con derivata distribuzionale (1o e 2o ordine);
4. decomposizione operatoriale generale del kernel quadratico in blocchi tensor/vector/scalar con algebre dei proiettori verificate;
5. identificazione esplicita dei coefficienti di spin dalla variazione seconda covariante (`08`);
6. mappa di stabilita completa con settore vettoriale esplicito e diagonalizzazione scalare (`09`);
7. mappa first-principles `eps-expansion -> {a,b,A,B} -> qKernel` verificata (`10`);
8. dizionario esplicito `metric+Theta -> {m,c,p,m2,c2,k2}` su ansatz ADM (`11`);
9. algebra+inversione del kernel perturbativo in proiettori full 4D (`12`);
10. dispersioni complete tensor/vector/scalar con controllo UV, mixing scalare e scaling IR (`13`).
11. passaggio 4D generale (senza ansatz simmetrico) ai coefficienti perturbativi metrici/`Theta`, con gauge-fixing implementato in forma parametrica/proxy (`scripts/14_full4d_covariant_gauge_fixed_coefficients.wl`); la derivazione esplicita completa di `Delta c_gf` da `F_mu,G` resta aperta.
12. discussione EFT esplicita:
   cutoff operativo come massimo `q^2` testato (`q2CutProxy=4` nel dominio testato),
   proxy strong-coupling (assenza di segnale nelle soglie adottate),
   regime di validita (`epsNL << 1` con stime `Lmin`) da `scripts/15_eft_validity_scales_proxy.wl`.

Non ancora coperto:
1. prova BRST/Faddeev-Popov di indipendenza dal gauge sui residui fisici nella pipeline numerica;
2. cutoff UV fondamentale e completamento UV oltre il cutoff operativo della pipeline;
3. estendere il dizionario del punto `11` da ansatz ADM controllato a decomposizione perturbativa 4D completamente generale (senza scelte di simmetria/fissaggi parziali);
4. confronto diretto con osservabili perturbative cosmologiche (spettri/transfer) in un setup fenomenologico completo.

## Esecuzione (Wolfram Cloud trial)
Autenticazione una tantum:

```powershell
& "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" -authenticate
```

Esecuzione completa:

```powershell
cd c:\lean\StaticSanityCheck\Paper\Perturbative_Spectrum_Study
.\run_all.ps1
```
