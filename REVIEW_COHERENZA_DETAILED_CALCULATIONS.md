# Revisione Coerenza Interna - detailed_calculations.md

Data revisione: 2026-02-18

## Perimetro analizzato
1. `Paper/Perturbative_Spectrum_Study/detailed_calculations.md`
2. `Paper/Variational_Derivation/detailed_calculations.md`
3. `Paper/Barotropic_Cosmology_Study/detailed_calculations.md`

## Check eseguiti
1. Coerenza struttura documento: ogni sezione ha riferimento a script/log e conclusione fisica.
2. Coerenza dichiarazioni vs log:
   - `Paper/Perturbative_Spectrum_Study/logs/*.log`: nessun `False`.
   - `Paper/Variational_Derivation/logs/*.log`: nessun `False`.
   - `Paper/Barotropic_Cosmology_Study/logs/*.log`: nessun `False`.
3. Coerenza inter-cartella:
   - blocco variazionale usato come base del blocco perturbativo;
   - blocco barotropico coerente col limite hard dichiarato.
4. Audit assunzioni/overclaim:
   - verificata distinzione tra risultati dimostrati e proxy numerici;
   - cercate affermazioni piu forti del supporto effettivo.

## Esito sintetico
1. Non emergono contraddizioni logiche gravi tra i tre `detailed_calculations.md`.
2. Sono presenti alcuni punti con assunzioni forti o formulazione troppo assertiva rispetto a quanto implementato nei notebook/script.

## Finding (ordinati per severita)

### 1. ALTA - Overclaim su "derivazione diretta 4D con gauge-fixing pienamente covariante"
1. In `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:491`-`Paper/Perturbative_Spectrum_Study/detailed_calculations.md:513` si presenta il gauge-fixing come derivato a livello d'azione 4D.
2. In `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:546` si dichiara il "gap referee" chiuso.
3. Pero lo script operativo `Paper/Perturbative_Spectrum_Study/scripts/14_full4d_covariant_gauge_fixed_coefficients.wl:53`-`Paper/Perturbative_Spectrum_Study/scripts/14_full4d_covariant_gauge_fixed_coefficients.wl:57` usa shift parametrico (`g2,g1,gs,gw,gm`), senza derivazione esplicita di `Delta c_gf` da `F_mu,G`.
4. Quindi il risultato e una chiusura parametrica/proxy, non ancora una derivazione completa da gauge-fixing covariante esplicito.

### 2. MEDIA - "Cutoff" definito in modo tautologico
1. In `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:570` e `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:600` si introduce `q2CutProxy`.
2. In `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl:117`, `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl:136`, `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl:137` il valore e fissato da griglia (`q2Grid={0.1,1,4}`, `q2CutProxy=Max[q2Grid]`, `checkF==4`).
3. Questo e coerente come "massimo testato", ma non e una stima dinamica del cutoff EFT.

### 3. MEDIA - Strong coupling trattata con proxy e soglie arbitrarie
1. In `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:599` si legge "Strong coupling non osservata".
2. In `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl:91`-`Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl:99` i criteri dipendono da soglie convenzionali (`zFloor=0.1`, `invCeil=10`).
3. La conclusione va qualificata esplicitamente come "assenza di segnale nel proxy scelto".

### 4. BASSA - Artefatto di encoding nel testo EFT
1. In `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:571` compare `xScalarÂ±`.
2. E un problema editoriale, non fisico, ma riduce chiarezza e puo creare ambiguita.

### 5. BASSA - Criterio "tutti check=True" non uniforme nella cartella Variational
1. `Paper/Variational_Derivation/detailed_calculations.md:599` dichiara "`00..25` check=True".
2. I log iniziali (es. `00`, `02`, `03`, `04`, `05`) usano `checkA/checkB` o `check1/check2` senza riga finale unica `check=True`.
3. Non e un errore matematico, ma indebolisce la verificabilita automatica uniforme.

## Punti da correggere (proposti)
1. Riformulare in `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:546`:
   - da "gap chiuso" a "chiusura parametrica/proxy del passaggio 4D; derivazione esplicita di `Delta c_gf` ancora aperta".
2. Riformulare in `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:600`:
   - da "cutoff operativo" a "massimo `q^2` testato con stabilita nella griglia corrente".
3. Riformulare in `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:599`:
   - esplicitare "no strong-coupling signal nel proxy con soglie adottate".
4. Correggere encoding in `Paper/Perturbative_Spectrum_Study/detailed_calculations.md:571`:
   - sostituire `xScalarÂ±` con `xScalarPlus,xScalarMinus` oppure `xScalar+-` in ASCII.
5. Standardizzare i log variational:
   - aggiungere riga finale `check=True` anche agli script `00`-`05` per omogeneita con il resto della pipeline.

## Nota finale
1. La pipeline e complessivamente coerente e solida sul piano algebrico/numerico.
2. I punti da correggere sono principalmente di livello interpretativo (forza delle affermazioni) e qualita di reporting/verificabilita.
