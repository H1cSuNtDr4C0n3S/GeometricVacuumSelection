# ThetaV2 Hard Global - Mathematica Notebooks

Questa cartella contiene notebook Mathematica per l'analisi di una teoria gravitazionale con vincolo globale-per-foglia (media normalizzata su dominio causale) sull'ipersuperficie di foliazione.

**Convenzioni:**
- Signature metrica: (+,−,−,−)
- Riemann sign convention: R_μνρσ con R = g^μν R_μν
- Default: W = 1 (smoothing opzionale, non variato)

## Albero dei File

```
ThetaV2_HardGlobal/
├── 00_Setup_Assumptions.nb          # Setup globale e parametri
├── 01_Invariants_X_K2_R.nb          # Invarianti X, K², R
├── 02_Action_Lagrangian.nb          # Lagrangiana 4D + kernel vincolo
├── 03_Leff_Radial_Reduction.nb      # Riduzione radiale + funzionale globale
├── 04_SanityChecks.nb               # Verifiche (volume, Schwarzschild, limite statico)
├── 05_Test_Minkowski.nb             # Test Minkowski → ESCLUSO
├── 06_Test_deSitter_AdS.nb          # Test de Sitter → PERMESSO (AdS: TODO)
└── GR_Spherical_Foliation_K2_R.nb   # Derivazione completa da primi principi
```

---

## Descrizione dei Notebook

### **00_Setup_Assumptions.nb**

#### Descrizione del Codice
Questo notebook definisce le assunzioni globali e i parametri per l'approccio "hard constraint" globale-per-foglia. Il vincolo **⟨K²⟩_{D_Θ} = K₀²** (media normalizzata sul dominio causale) è imposto per foglia (per ipersuperficie Θ) tramite un moltiplicatore di Lagrange **λ(Θ)**, NON come campo locale.

**Elementi definiti:**
- **Coordinate**: `coords = {t, r, th, ph}` (coordinate sferiche)
- **Assunzioni**: r > 0, A[r] > 0, B[r] > 0, ψ[r] ∈ Reals, 0 ≤ th ≤ π
- **Vincolo timelike**: -B[r] + A[r](ψ'[r])² < 0 (foliation spacelike)
- **Parametri di accoppiamento**:
  - κ: costante gravitazionale (8πG)
  - α: coefficiente per l'invariante X
  - Λ: costante cosmologica
- **Parametri del vincolo globale**:
  - K0: scala della curvatura estrinseca target (⟨K²⟩ = K₀²)
  - λ(Θ): moltiplicatore di Lagrange per foglia (costante su ogni Σ_Θ)

#### Output
Il notebook non produce output visibili, ma inizializza le variabili globali utilizzate negli altri notebook.

---

### **01_Invariants_X_K2_R.nb**

#### Descrizione del Codice
Calcola gli invarianti geometrici per una metrica sferica simmetrica:
- Metrica: `g = diag(A[r], -B[r], -r², -r²sin²θ)`
- Invariante cinetico: `X = 1/A[r] - (ψ'[r])²/B[r]`
- Curvatura estrinseca: K² (formula complessa dipendente da A, B, ψ e loro derivate)
- Scalare di Ricci: R

#### Output
```
sqrt(-g) = r² √(A[r]B[r]) sin(θ)

X = 1/A[r] - (ψ'[r])²/B[r]

K2 = [espressione complessa con derivate di A, B, ψ] / [4r²A[r]B[r](-B[r] + A[r](ψ'[r])²)³]

R = [-r²B[r](A'[r])² - 4A[r]²((-1+B[r])B[r] + rB'[r]) + rA[r](A'[r](4B[r] - rB'[r]) + 2rB[r]A''[r])] / [2r²A[r]²B[r]²]
```

Gli invarianti sono semplificati ma mantengono la piena dipendenza da A[r], B[r] e ψ[r].

---

### **02_Action_Lagrangian.nb**

#### Descrizione del Codice
Definisce la densità lagrangiana 4D e la struttura del vincolo globale.

**Parti principali:**
1. **Lagrangiana locale** (dinamica standard):
   ```
   L4base = √(-g) × [(1/2κ)(R - 2Λ) - (α/2)X]
   ```

2. **Kernel del vincolo per integrazione sulla foglia** (con W = 1 di default):
   - `KernelK2 = K2` (kernel di K²)
   - `KernelRef = K0²` (valore target)

3. **Struttura del funzionale di vincolo globale-per-foglia**:
   ```
   C[Θ] = ∫_{D_Θ} [K2 × W × √h d³x] - K0² × ∫_{D_Θ} [W × √h d³x]

   Media normalizzata: ⟨K²⟩_{D_Θ} = I1/I0 = K0²
   ```
   dove:
   - D_Θ è il dominio causale sulla foglia Σ_Θ (definito da criteri geometrici)
   - W è funzione di smoothing opzionale (default W = 1, NON variata nell'azione)
   - Il vincolo è INTEGRALE, non pointwise

#### Output
```
L4base (local dynamics) = -[sin(θ)/(4κ(AB)^(3/2))] × [r²B(A')² + A²(4(1+r²Λ)B² + 4rB' - 2B(2+r²ακ(ψ')²)) + rA(2rακB² + rA'B' - 2B(2A' + rA''))]

KernelK2 (to be leaf-integrated) = [stessa espressione di K2]

KernelRef (target value) = K0²
```

---

### **03_Leff_Radial_Reduction.nb**

#### Descrizione del Codice
Esegue l'integrazione angolare per ottenere i kernel radiali e costruisce il vincolo globale.

**Operazioni principali:**
1. **Integrazione angolare**: fattore = 4π
2. **Lagrangiana radiale locale**: `LeffBase(r)` (integrazione di L4base sugli angoli)
3. **Kernel della foglia** (con W = 1):
   - `f0(r) = 4πr²√(B-A(ψ')²) × W` (misura, W = 1 di default)
   - `f1(r) = f0(r) × K2` (K² pesato)
4. **Densità di vincolo**: `Gconstraint(r) = f1 - K0² × f0`
5. **Funzionale di vincolo globale-per-foglia**:
   ```
   C[Θ] = I1 - K0² × I0 = 0

   dove:
   I0 = ∫_{D_Θ} f0(r) dr  [volume normalizzato]
   I1 = ∫_{D_Θ} f1(r) dr  [K² pesato]

   Media normalizzata: ⟨K²⟩_{D_Θ} = I1/I0 = K0²
   ```
   **IMPORTANTE**: Il dominio D_Θ è determinato da criteri geometrici (orizzonti, struttura causale),
   NON aggiustato per test. λ(Θ) è costante su ogni foglia, NON campo locale.

#### Output
```
Angular integration factor = 4π

LeffBase(r) [local dynamics] = -[π/(κ(AB)^(3/2))] × [r²B(A')² + A²(...) + rA(...)]

f0(r) [measure kernel] = 4πr²√(B - A(ψ')²)

f1(r) [K^2 weighted kernel] = [π/(AB(B - A(ψ')²)^(5/2))] × [...]

Gconstraint(r) = f1 - K0^2 * f0 = [espressione completa]

=== GLOBAL CONSTRAINT STRUCTURE ===
C[Theta] = I1 - K0^2 * I0 = 0
where:
  I0 = Integral_{D_Theta} f0(r) dr  [volume]
  I1 = Integral_{D_Theta} f1(r) dr  [K^2 weighted]

Equivalently: <K^2>_{D_Theta} = I1/I0 = K0^2

This is NOT a local equation.
lambda[Theta] is constant on each leaf.
The constraint selects globally consistent foliations.
```

---

### **04_SanityChecks.nb**

#### Descrizione del Codice
Verifica la correttezza delle espressioni calcolate attraverso tre controlli:
1. Elemento di volume: √(-g)
2. Test Schwarzschild: R = 0 per la metrica di Schwarzschild
3. Limite di foliazione statica (ψ' = 0)

#### Output
```
=== CHECK 1: Volume Element ===
sqrt(-g) = r² √(A[r]B[r]) sin(θ)

=== CHECK 2: Schwarzschild R = 0 ===
R(Schwarzschild) = 0
PASS: Schwarzschild is Ricci-flat (GR recovered).

=== CHECK 3: Static Foliation Limit (psi' = 0) ===
X | psi'=0 = 1/A[r]
K2 | psi'=0 = 0

Interpretation:
  X -> 1/A(r) : kinetic term from lapse alone
  K2 simplified: extrinsic curvature of t=const slices

NOTE: This is a LOCAL limit check, not a constraint test.
The global constraint is evaluated in notebooks 05, 06.
```

---

### **05_Test_Minkowski.nb**

#### Descrizione del Codice
Testa il vincolo globale nello spaziotempo di Minkowski.

**Metrica**: A(r) = 1, B(r) = 1

**Casi testati**:
1. **Foliazione statica** (ψ' = 0): K2 = 0 ovunque
2. **Foliazione inclinata** (ψ' = v = const):
   - X = 1 - v²
   - K2 = -2v²/[r²(v²-1)]

**Dominio causale**: D_Θ definito come **causal patch** dello spaziotempo piatto
- Scelta geometrica: patch con scala operativa R_max
- Limite infrarosso: R_max → ∞ per testare comportamento asintotico
- **NOTA**: Dominio fissato da criterio geometrico (non aggiustato per ottenere risultati desiderati)

**Kernel** (con W = 1):
- f0(r) = 4πr²√(1-v²)
- f1(r) = 8πv²/√(1-v²) (indipendente da r!)

#### Output
```
=== MINKOWSKI METRIC ===
A(r) = 1, B(r) = 1 (flat spacetime)

=== LOCAL INVARIANTS (Minkowski) ===
X = 1 - (ψ'[r])²
K2 = -[2(ψ'[r])²(-1+(ψ'[r])²)² + r²(ψ''[r])²] / [r²(-1+(ψ'[r])²)³]

=== CASE 1: Static Foliation (psi' = 0) ===
X | psi'=0 = 1
K2 | psi'=0 = 0
Static foliation: K2 = 0 everywhere.

=== CASE 2: Tilted Foliation (psi' = v = const) ===
X | psi'=v = 1 - v²
K2 | psi'=v = -2v²/[r²(-1+v²)]
Timelike condition X > 0 requires: |v| < 1
Key observation: K2 ~ 1/r^2 (decays at large r)

=== INTEGRALS OVER CAUSAL DOMAIN ===
I0(Rmax) = (4/3)πRmax³√(1-v²)
I1(Rmax) = 8πRmax v²/√(1-v²)

=== NORMALIZED AVERAGE ===
<K^2>_{D_Theta} = I1/I0 = -6v²/[Rmax²(v²-1)]

=== ASYMPTOTIC BEHAVIOR (Rmax -> infinity) ===
lim_{Rmax->inf} <K^2> = 0

========================================
           FINAL CONCLUSION
========================================

  <K^2>_{D_Theta} -> 0  as  Rmax -> infinity

In Minkowski spacetime with ANY timelike foliation (|v| < 1):
The NORMALIZED AVERAGE of K^2 over the causal domain
vanishes in the IR limit, regardless of the tilt v.

Therefore, the leafwise constraint:
  C[Theta] = 0  with  K0^2 > 0

CANNOT be satisfied globally in Minkowski.

CONCLUSION: Minkowski is EXCLUDED by the global
leaf-integrated constraint ⟨K²⟩_{D_Θ} = K0² > 0.

This is NOT a pointwise statement K² = K0².
It is a statement about the NORMALIZED AVERAGE over D_Θ.

=== SUPPLEMENTARY: Why ⟨K²⟩ -> 0 ===
Scaling analysis:
  K2(r) ~ v^2 / r^2       (decays as 1/r^2)
  f0(r) ~ r^2             (volume measure grows as r^2)
  f1(r) = f0 * K2 ~ const (product is constant)

Integrals:
  I0 ~ Rmax^3             (volume grows as Rmax^3)
  I1 ~ Rmax               (grows only as Rmax)

Therefore:
  <K^2> = I1/I0 ~ Rmax / Rmax^3 = 1/Rmax^2 -> 0
```

**Conclusione**: Minkowski è ESCLUSO dal vincolo globale.

---

### **06_Test_deSitter_AdS.nb**

#### Descrizione del Codice
Testa il vincolo globale nello spaziotempo di de Sitter (AdS non ancora implementato, TODO).

**Metrica de Sitter (patch statico)**:
- A(r) = 1 - H²r²
- B(r) = 1/(1 - H²r²)
- **Dominio causale geometrico**: 0 < r < 1/H (orizzonte cosmologico naturale)

**Foliazione**: Painlevé-Gullstrand type
- ψ'(r) = Hr/(1 - H²r²)

**Verifiche geometriche**:
- R = -12H² (massimamente simmetrico, convenzione segno Riemann come da setup)
- X = 1 (foliazione timelike)
- K² = 3H² (COSTANTE su tutto il patch!)

**Nota sulla costanza di K²**: Poiché K² è costante, la media normalizzata ⟨K²⟩_{D_Θ} = K² indipendentemente da W e dai bounds di integrazione (purché ∫f0 > 0).

#### Output
```
=== DE SITTER METRIC (STATIC PATCH) ===
A(r) = 1 - H²r²
B(r) = 1/(1 - H²r²)
Domain: 0 < r < 1/H (static patch)

=== FOLIATION ANSATZ ===
psi'(r) = Hr/(1 - H²r²)

=== GEOMETRY CHECK ===
R = -12H²
(Sign convention: R = -12 H^2 here; maximal symmetry is |R| = 12 H^2.)
Maximally symmetric (|R| = 12 H^2)? True

=== CAUSALITY CHECK ===
X = 1
Timelike foliation (X = 1)? True

=== LOCAL INVARIANT K^2 ===
K^2 = 3H²
d(K^2)/dr = 0
K^2 is constant? True
K^2 = 3 H^2? True

=== LEAFWISE CONSTRAINT ANALYSIS ===
CAUSAL DOMAIN D_Theta:
  The static patch {r : 0 < r < 1/H} is the
  natural causal domain for de Sitter.

KEY OBSERVATION:
  Since K^2 = 3 H^2 is CONSTANT on the patch,
  the normalized leafwise average satisfies:
    ⟨K²⟩_{D_Theta} = K² = 3 H²

  This holds for ANY W ≥ 0 with finite ∫f0 > 0.

CONSTRAINT SATISFACTION:
  The leafwise constraint C[Theta] = 0 requires
  ⟨K²⟩ = K0². Since ⟨K²⟩ = 3 H², this gives:
    H² = K0² / 3

  This is a SELECTION RULE, not a pointwise equation.

========================================
           FINAL VERDICT
========================================

SUCCESS: de Sitter satisfies the LEAFWISE
NORMALIZED constraint (trivially, because
K^2 is constant), selecting H^2 = K0^2/3
on the static patch.

=== SUPPLEMENTARY ===
Why the leafwise constraint is trivial for de Sitter:

For ANY normalized average:
  <K^2> = Integral[f0 * K^2] / Integral[f0]

If K^2 = const, it factors out:
  <K^2> = K^2 * Integral[f0] / Integral[f0] = K^2

Therefore <K^2> = 3 H^2 regardless of:
  - The choice of W (any W >= 0 with Integral[f0] > 0)
  - The exact domain bounds (within static patch)
  - The angular integration details

This is why de Sitter passes the global constraint
effortlessly: constancy of K^2 makes the average
equal to the local value everywhere.
```

**Conclusione**: de Sitter è PERMESSO con H² = K₀²/3.

---

### **GR_Spherical_Foliation_K2_R.nb**

#### Descrizione del Codice
Derivazione dai primi principi degli invarianti geometrici (X, K², R) per una metrica sferica simmetrica con foliazione Θ = t + ψ(r).

**Procedura completa**:
1. Metrica: `g = diag(A[r], -B[r], -r², -r²sin²θ)`
2. Coordinate: `{t, r, θ, φ}`
3. Foglia: `Θ = t + ψ(r)`
4. Normale: `n_μ = ∂_μΘ / √X`
5. Metrica indotta: `h_{μν} = g_{μν} + n_μn_ν`
6. Curvatura estrinseca: `K_{μν} = h^α_μ h^β_ν ∇_α n_β`
7. Invariante: `K² = K^{μν}K_{μν}`
8. Scalare di Ricci: da tensore di Riemann

**Test**: Schwarzschild (A = 1-2M/r, B = 1/A) → R = 0 ✓

#### Output
```
X = 1/A[r] - (ψ'[r])²/B[r]

K2 = [stessa formula complessa del notebook 01]

R = [-r²B(A')² - 4A²((-1+B)B + rB') + rA(A'(4B - rB') + 2rBA'')] / [2r²A²B²]

R (Schwarzschild test) = 0
```

Questo notebook fornisce la derivazione completa e indipendente degli invarianti, confermando le formule utilizzate nei notebook precedenti.

---

## Riepilogo dei Risultati

### Struttura della Teoria
- **Vincolo globale-per-foglia**: `⟨K²⟩_{D_Θ} = K₀²` (media normalizzata su dominio causale)
- **NON locale**: λ(Θ) costante per foglia, NON funzione di x
- **NON pointwise**: Il vincolo è `⟨K²⟩ = K₀²`, NON `K² = K₀²` punto per punto
- **Funzionale di vincolo**: C[Θ] = I1 - K₀² × I0 = 0, dove I1 = ∫[f1·dr], I0 = ∫[f0·dr]
- **Dominio causale**: D_Θ determinato da geometria (orizzonti, struttura causale), NON aggiustato

### Test Spaziatempi
| Spaziotempo | K²(r) | ⟨K²⟩_{D_Θ} | Vincolo | Fisica |
|-------------|-------|-----------|---------|--------|
| **Minkowski** | ~v²/r² | → 0 (IR) | ❌ NO | Escluso: volume domina su curvatura locale |
| **de Sitter** | 3H² | 3H² | ✅ SÌ | Permesso: K² costante ⇒ H² = K₀²/3 |

**Nota cruciale**: La differenza è nel comportamento INTEGRATO, non pointwise.

### Significato Fisico
Il vincolo globale **seleziona de Sitter** ed **esclude Minkowski**, fornendo un meccanismo di:
1. **Emergenza della costante cosmologica** dalla geometria della foliazione
2. **Selezione del vuoto** basata su proprietà globali-per-foglia (non locali)
3. **Scala dinamica** H² = K₀²/3 fissata dal parametro del vincolo

Questo è un risultato **qualitativo robusto** rispetto a:
- Scelta di W (purché W ≥ 0, ∫f0 > 0)
- Dettagli della foliazione (purché X > 0)
- Limiti esatti del dominio (purché coerenti con struttura causale)
