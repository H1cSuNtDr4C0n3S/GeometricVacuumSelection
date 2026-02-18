# Barotropic Cosmology Study Summary

Mappa notebook `Paper/Barotropic_Cosmology_Study` -> passaggi matematici, output verificati e lettura fisica.

## Notebook 00

1. `00_Barotropic_Hard_Constraint_Consistency.nb`
`Riferimento`: minisuperspace FRW con materia barotropica e vincolo geometrico hard.
`Output`:
- `ElamDS = Exp[3 H0 t] (3 H0^2 - K0^2)`.
- `hardRule = H0^2 -> K0^2/3`.
- `lamExpr = -(1/2) (-K0^2 + Lambda + (kappa rho0)/(Exp[3 H0 t] (Exp[H0 t])^(3 ww))) / (K0^2 kappa)`.
- `lamPrimeExpr = (3 H0 rho0 (1 + ww)) / (2 Exp[3 H0 t] (Exp[H0 t])^(3 ww) K0^2)`.
- `EaCons = (3 (Exp[H0 t])^(-1 - 3 ww) (-3 H0^2 + K0^2) rho0 (1 + ww)) / K0^2`.
- `checkA=True`, `checkB=True`, `checkC=True`, `checkD=True`, `checkE=True`, `checkF=True`, `check=True`.
`Perche questo passaggio`: chiarisce cosa succede quando il vincolo `H^2=K0^2/3` e' imposto in modo hard in presenza di materia barotropica.
`Interpretazione fisica dell'output`: sul ramo hard, il canale del vincolo forza de Sitter; la materia non riapre una dinamica standard di `H(t)`, ma viene assorbita nel profilo temporale di `lambda(t)`. Per `w=-1`, `lambda` torna statico (`lamPrimeDE=0`), coerente con un settore effettivamente tipo costante cosmologica.

## Nota operativa

Generazione e verifica:

```powershell
cd c:\lean\StaticSanityCheck\Paper\Barotropic_Cosmology_Study
.\run_all.ps1
```

Risultati:
- notebook in `Paper/Barotropic_Cosmology_Study/*.nb`
- log in `Paper/Barotropic_Cosmology_Study/logs/*.log`

## Nota EFT / scale di validita
1. Il blocco barotropico corrente (`00`) e' un test di coerenza sul vincolo hard, non una determinazione completa EFT di cutoff e strong coupling.
2. In particolare, il ramo hard resta de Sitter e non fornisce da solo una finestra cosmologica realistica materia/radiazione.
3. La discussione EFT operativa del settore perturbativo e' in `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`.
