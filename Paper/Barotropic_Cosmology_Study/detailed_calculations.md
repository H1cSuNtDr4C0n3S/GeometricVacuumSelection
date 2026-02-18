# Barotropic Cosmology Study - Detailed Calculations and Results

Questo documento dettaglia i calcoli eseguiti in `Paper/Barotropic_Cosmology_Study/scripts/` e i risultati registrati in `Paper/Barotropic_Cosmology_Study/logs/`.

Convenzioni:
1. `a(t)` fattore di scala.
2. `N(t)` lapse.
3. `lambda(t)` moltiplicatore di Lagrange del vincolo geometrico.
4. `rho0`, `w` parametri della materia barotropica (`rho(a) ~ a^{-3(1+w)}`).
5. `kappa`, `Lambda`, `K0` parametri del modello.

## 00 - Barotropic Hard-Constraint Consistency

Riferimenti:
1. Script: `scripts/00_barotropic_hard_constraint_consistency.wl`.
2. Log: `logs/00_barotropic_hard_constraint_consistency.log`.
3. Notebook generato: `00_Barotropic_Hard_Constraint_Consistency.nb`.

### 1) Lagrangiano minisuperspace usato

Lo script imposta:

`Lbar = -3 a a'^2/(kappa N) - (Lambda a^3 N)/kappa + lambda (3 a a'^2/N^2 - K0^2 a^3) - rho0 a^(-3 w) N`.

Commento tecnico:
1. Primo e secondo termine: parte GR ridotta FRW (kinetic + `Lambda`).
2. Terzo termine: canale vincolato con `lambda`.
3. Quarto termine: settore materia barotropica efficace.

### 2) Equazioni del moto simboliche

Si derivano le tre EOM:
1. `Ea = dL/da - d/dt(dL/da')`.
2. `EN = dL/dN`.
3. `Elam = dL/dlambda`.

Poi si specializza su ansatz de Sitter:
1. `a(t)=Exp[H0 t]`.
2. `a'(t)=H0 Exp[H0 t]`.
3. `a''(t)=H0^2 Exp[H0 t]`.
4. `N(t)=1`.

Output dal log:
1. `EaDS = (-3*(-3*E^(2*H0*t)*H0^2 + E^(2*H0*t)*Lambda - (E^(H0*t))^(-1 - 3*ww)*kappa*rho0*ww + E^(2*H0*t)*(3*H0^2 + K0^2)*kappa*ll[t] + 2*E^(2*H0*t)*H0*kappa*Derivative[1][ll][t]))/kappa`.
2. `ENDS = -((E^(3*H0*t)*Lambda + (kappa*rho0)/(E^(H0*t))^(3*ww) + 3*E^(3*H0*t)*H0^2*(-1 + 2*kappa*ll[t]))/kappa)`.
3. `ElamDS = E^(3*H0*t)*(3*H0^2 - K0^2)`.

Interpretazione:
1. `ElamDS=0` impone direttamente `H0^2 = K0^2/3` (branch hard).
2. Sul branch hard, `EN` non determina `H(t)` ma determina `lambda(t)`.

### 3) Imposizione vincolo hard e soluzione per lambda(t)

Regola usata:
1. `hardRule = H0^2 -> K0^2/3`.

Da `ENDS=0` si ricava:
1. `lamRule = {ll[t] -> (E^(3*H0*t)*(E^(H0*t))^(3*ww)*K0^2 - E^(3*H0*t)*(E^(H0*t))^(3*ww)*Lambda - kappa*rho0)/(2*E^(3*H0*t)*(E^(H0*t))^(3*ww)*K0^2*kappa)}`.
2. `lamExpr = -1/2*(-K0^2 + Lambda + (kappa*rho0)/(E^(3*H0*t)*(E^(H0*t))^(3*ww)))/(K0^2*kappa)`.
3. `lamPrimeExpr = (3*H0*rho0*(1 + ww))/(2*E^(3*H0*t)*(E^(H0*t))^(3*ww)*K0^2)`.

Interpretazione:
1. Con `w != -1`, `lambda'(t)` e' in generale non nullo.
2. La dipendenza temporale di `lambda` codifica la presenza della materia sul branch geometricamente rigido.

### 4) Riduzione del canale Ea sul branch hard

Dopo sostituzione di `hardRule`, `lamRule`, `dlamRule`:

1. `EaCons = (3*(E^(H0*t))^(-1 - 3*ww)*(-3*H0^2 + K0^2)*rho0*(1 + ww))/K0^2`.

Poiche sul branch hard vale `-3 H0^2 + K0^2 = 0`, segue:
1. `EaCons = 0` (ridondanza del canale `Ea`).

Interpretazione:
1. Nel setup hard, i canali dinamici non danno evoluzione indipendente di `H`.
2. La struttura e' consistente con un settore fortemente vincolato, non con una dinamica FRW standard materia/radiazione.

### 5) Verifiche automatiche implementate

Check dal log:
1. `checkA(Elam) = True`.
2. `checkB(vacuum) = True`.
3. `checkC(w=-1) = True`.
4. `matterValsDust = {0, 0, 0, 0}`.
5. `maxDust = 0`.
6. `checkD(Ea redundancy on hard branch) = True`.
7. `lamPrimeDustVals = {1.4433756729740644..., 0.8584467265031961..., 0.5105606226032707...}`.
8. `maxLamPrimeDust = 1.4433756729740644...`.
9. `checkE(lambda time dependence for dust) = True`.
10. `lamPrimeDE(w=-1) = 0`.
11. `checkF(lambda static for w=-1) = True`.
12. `check = True`.

Interpretazione fisica dei check:
1. `checkA`: il ramo hard e' effettivamente imposto dal canale vincolato.
2. `checkD`: una volta imposto il ramo, il canale `Ea` diventa identita.
3. `checkE`: per dust/radiazione la materia non muove `H`, ma muove `lambda`.
4. `checkF`: per `w=-1` il settore materia si comporta come energia di vuoto, e `lambda` puo' essere costante.

## Conclusione tecnica del blocco 00

1. Il notebook e' internamente consistente (`check=True`).
2. Nel regime FRW omogeneo con vincolo hard, il ramo de Sitter e' fissato direttamente da `Elam`.
3. La materia barotropica entra in `lambda(t)` ma non produce, in questo setup, una fase standard con `H(t)` evolutivo.
4. Questo conferma che per una cosmologia completa (materia/radiazione + late-time de Sitter) serve il passo successivo a dominio finito/causale, non il vincolo hard applicato uniformemente a ogni tempo.

## Riproduzione rapida

```powershell
cd c:\lean\StaticSanityCheck\Paper\Barotropic_Cosmology_Study
.\run_all.ps1
```

Artefatti:
1. Notebook: `Paper/Barotropic_Cosmology_Study/*.nb`.
2. Log: `Paper/Barotropic_Cosmology_Study/logs/*.log`.

## Nota EFT / scale di validita
1. Questo blocco controlla coerenza analitica del caso FRW barotropico con vincolo hard, ma non fissa cutoff UV e non testa direttamente strong coupling.
2. In termini EFT, il risultato indica che il branch hard assorbe la materia in `lambda(t)`; per un regime cosmologico realistico serve il passo finite-domain/causal gia segnalato.
3. La discussione quantitativa di cutoff operativo e strong-coupling proxy e' stata implementata nel blocco perturbativo: `Paper/Perturbative_Spectrum_Study/scripts/15_eft_validity_scales_proxy.wl`.
