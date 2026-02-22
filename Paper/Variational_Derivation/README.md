# Variational_Derivation

Questa cartella contiene la pipeline per la derivazione variazionale completa,
con verifica obbligatoria via Mathematica (`wolframscript`) prima della consegna.

## Struttura
- `scripts/00_variational_seed.wl`
- `scripts/01_ratio_firstvariation.wl`
- `scripts/02_constraint_firstvariation.wl`
- `scripts/03_theta_chainrule.wl`
- `scripts/04_ratio_kernel_form.wl`
- `scripts/05_metric_variation_spherical.wl`
- `scripts/06_second_order_EL_identity.wl`
- `scripts/07_k2_psi_variation_numeric.wl`
- `scripts/08_k2_psi_EL_operator.wl`
- `scripts/09_metric_kernels_AB_numeric.wl`
- `scripts/10_full_proxy_eom_operators.wl`
- `scripts/11_full_proxy_variational_checks.wl`
- `scripts/12_finite_domain_ratio_variation.wl`
- `scripts/13_ir_flow_minkowski_desitter.wl`
- `scripts/14_constraint_finite_domain_boundary.wl`
- `scripts/15_frw_k2_leafwise_average.wl`
- `scripts/16_frw_minisuperspace_eom.wl`
- `scripts/17_frw_linearized_constraint_channel.wl`
- `scripts/18_desitter_constraint_propagation.wl`
- `scripts/19_frw_linearized_full_eom_reduction.wl`
- `scripts/20_minisuperspace_noether_bianchi_identity.wl`
- `scripts/21_frw_curvature_minisuperspace_eom.wl`
- `scripts/22_frw_curvature_noether_identity.wl`
- `scripts/23_covariant_causal_domain_variation.wl`
- `scripts/24_covariant_boundary_shape_derivative.wl`
- `scripts/25_covariant_nonlocal_action_variation.wl`
- `scripts/26_covariant_general_eom_channels.wl`
- `logs/` (riempita dopo esecuzione)
- `summary.md` (mappa notebook <-> formule del paper)
- `detailed_calculations.md` (dettaglio passaggi matematici + output log)

## Esecuzione (Wolfram Cloud trial)
Autenticati una volta:

```powershell
& "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" -authenticate
```

Poi, da PowerShell nella cartella `Paper/Variational_Derivation`:

```powershell
.\run_all.ps1
```

`run_all.ps1` esegue gli script in cloud (`-o`), raccoglie il payload JSON
restituito da Mathematica e salva localmente:
- notebook `.nb` nella cartella corrente;
- log testuali nella cartella `logs/`.

Nota: `run_all.ps1` esegue automaticamente tutti i file `*.wl` presenti
in `scripts/`, ordinati per nome.

## Copertura EOM rispetto al paper
La catena `23-26` estende i blocchi minisuperspace con una formulazione 4D covariante dei funzionali leafwise:
1. rappresentazione covariante di `I0`, `I1` con `delta(Theta-theta)` e indicatore causale `chi_Theta` (aderente a `Paper/main.tex:242-263`);
2. variazione del bordo tramite derivata di forma (`Heaviside -> DiracDelta`) coerente con la discussione su domini finiti/causali (`Paper/main.tex:278-281`, `Paper/main.tex:524-528`);
3. variazione dell'azione non-locale `S_NL = integral dTheta lambda(Theta) [I1-K0^2 I0]` in split bulk+boundary, in linea con `Paper/main.tex:292-301` e con l'argomento di dominanza IR (`Paper/main.tex:314-329`).
4. derivazione simbolica certificata dei canali EOM generali (vincolo, metrica, scalare), con forma bulk esplicita del tensore efficace non-locale e verifica del ramo de Sitter con `Lambda=0`.

## Nota EFT / validita
Questa cartella chiude la parte variazionale (bulk+boundary) ma non fissa da sola cutoff UV o scala di strong coupling.
L'interpretazione EFT operativa (proxy di cutoff, criteri di strong coupling e finestra di validita) e' implementata in
`Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`,
usando i coefficienti quadratici estratti dalla pipeline perturbativa.
