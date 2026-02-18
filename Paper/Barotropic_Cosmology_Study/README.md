# Barotropic_Cosmology_Study

Questa cartella e' separata da `Paper/Variational_Derivation` e contiene il lavoro dedicato al settore cosmologico con materia barotropica.

## Scope
1. Non modifica la pipeline variazionale `00..22` gia' validata in `Paper/Variational_Derivation`.
2. Isola il problema cosmologico FRW con materia barotropica in una linea di sviluppo indipendente.
3. Prepara il passaggio successivo a una prescrizione finita/causale per ottenere una storia cosmologica realistica.

## Struttura
- `scripts/00_barotropic_hard_constraint_consistency.wl`
- `logs/00_barotropic_hard_constraint_consistency.log`
- `00_Barotropic_Hard_Constraint_Consistency.nb`
- `summary.md`
- `detailed_calculations.md`

## Cosa verifica il notebook 00
1. Deriva le EOM minisuperspace FRW del modello con materia barotropica e moltiplicatore `lambda`.
2. Impone il vincolo hard `H^2 = K0^2/3` (canale `Elam=0`) e risolve `EN=0` per `lambda(t)`.
3. Controlla coerenza in vuoto e in caso `w=-1`.
4. Mostra che, per dust (`w=0`), `lambda` resta dinamico (`lambda'(t) != 0`) pur mantenendo ridondanza del canale `Ea` sul ramo hard.

## Esito attuale
Dal log `logs/00_barotropic_hard_constraint_consistency.log`:
- `checkA=True`
- `checkB=True`
- `checkC=True`
- `checkD=True`
- `checkE=True`
- `checkF=True`
- `check=True`

Quindi il test implementato nel notebook 00 e' validato.

## Esecuzione (Wolfram Cloud trial)
Autenticati una volta:

```powershell
& "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" -authenticate
```

Poi, da PowerShell nella cartella `Paper/Barotropic_Cosmology_Study`:

```powershell
.\run_all.ps1
```

`run_all.ps1` esegue tutti gli script `.wl` in `scripts/`, salva i notebook `.nb` nella cartella corrente e i log testuali in `logs/`.

## Nota metodologica
Nel settore FRW omogeneo con vincolo hard su `Q=I1/I0=K^2`, il branch selezionato e' direttamente de Sitter (`H^2=K0^2/3`).
Per recuperare ere materia/radiazione realistiche serve il passo successivo: prescrizione a dominio finito/causale con correzioni finite-volume.

## Nota EFT / validita
Il notebook 00 verifica coerenza minisuperspace sul ramo hard ma non determina cutoff UV o scala di strong coupling.
Per una discussione EFT operativa (cutoff testato, strong-coupling proxy, regime di validita) si usa
`Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`.
