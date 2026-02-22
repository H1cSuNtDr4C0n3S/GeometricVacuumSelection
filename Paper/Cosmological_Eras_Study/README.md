# Cosmological_Eras_Study

Questa cartella e dedicata alla nuova pipeline sull'applicazione della teoria nelle varie ere cosmologiche.

## Scope iniziale
1. Definire una prescrizione operativa di dominio finito/causale in FRW.
2. Costruire test separati per era di radiazione, materia e de Sitter tardivo.
3. Verificare continuita tra ere e consistenza del canale `Q=I1/I0`.
4. Misurare quando i contributi non-locali restano IR-soppressi.
5. Preparare output riproducibili (notebook + log testuali).

## Stato attuale (00-11)
1. `00`: baseline dati/cosmologia consistente (`STANDARD_ERA_HISTORY_COMPATIBLE`).
2. `01`: classe IR in radiazione (`RADIATION_ERA_IR_SUPPRESSED`).
3. `02`: classe IR in materia (`MATTER_ERA_IR_SUPPRESSED`).
4. `03`: continuita globale tra ere (`ERAS_PIPELINE_CONTINUITY_CONFIRMED`).
5. `04`: equivalenza di classe tra proxy (`ERA_PROXY_CLASS_EQUIVALENT_IR`).
6. `05`: tabella finale `00-04` pronta per appendice (`PIPELINE_00_04_APPENDIX_READY`).
7. `06`: ricerca EOM complete con materia: ramo auto-consistente trovato sotto hard constraint (`COMPLETE_EOM_MATTER_SELFCONSISTENT_DS_BRANCH_ONLY`).
8. `07`: ricerca EOM complete con chiusura soft/finite-domain: ramo evolutivo trovato (`COMPLETE_EOM_SOFT_CLOSURE_EVOLVING_BRANCH_FOUND`).
9. `08`: ponte paper-ready hard/soft nelle EOM complete con mappa LCDM-like (`PAPER_READY_HARD_SOFT_EOM_BRIDGE_CONFIRMED`).
10. `09`: robustezza parametrica del ramo soft EOM + test anti-fine-tuning (`COMPLETE_EOM_SOFT_BRANCH_ROBUST_NOT_FINE_TUNED`).
11. `10`: sintesi appendix-ready del blocco EOM `06-09` (`PIPELINE_06_09_APPENDIX_READY`).
12. `11`: chiusura sezione in formato referee-ready (`COSMOLOGICAL_ERAS_SECTION_CLOSURE_READY`).

## Struttura
1. `scripts/`
2. `logs/`
3. `data/`
4. `summary.md`
5. `detailed_calculations.md`
6. `run_all.ps1`

## Data baseline (frozen)
1. Baseline cosmologica primaria: `Planck18` (Astropy realization, riferimento stabile).
2. Cross-check di robustezza: `DESI DR1 BAO VAC v1.0`.
3. File scaricati localmente:
   - `data/desi_dr1_v1.0_bestfit.minimum.txt`
   - `data/desi_dr1_v1.0.sha256sum`
4. Verifica integrita:
   - SHA256 `data/desi_dr1_v1.0_bestfit.minimum.txt`:
     `948784a698e81e31ef496e212bde038604b71ff4fed13ff7bda0a4ee2db12b25`
   - Match con entry ufficiale `iminuit/base/desi-bao-all/bestfit.minimum.txt` nel file checksum DR1 v1.0.
5. Nota esecuzione cloud:
   - se il kernel Wolfram Cloud non ha accesso a filesystem/rete locale, gli script usano una snapshot embedded con gli stessi valori e lo stesso hash, mantenendo output riproducibile.

## Fonti ufficiali
1. Planck18 (Astropy): https://docs.astropy.org/en/stable/api/astropy.cosmology.realizations.Planck18.html
2. Planck Legacy cosmological parameters: https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Cosmological_Parameters
3. DESI DR1 BAO VAC v1.0: https://data.desi.lbl.gov/public/dr1/vac/dr1/bao-cosmo-params/v1.0/

## Avvio rapido
```powershell
cd c:\lean\StaticSanityCheck\Paper\Cosmological_Eras_Study
.\run_all.ps1
```
