import Mathlib.Data.Real.Basic
import Mathlib.Tactic

import StaticSanityCheck.K2General

namespace StaticSanityCheck.FromGeneral

set_option autoImplicit false
noncomputable section

/-- de Sitter static patch: A(r)=1-H^2 r^2. -/
def dS_A (H r : Real) : Real := 1 - (H^2) * (r^2)

/-- de Sitter static patch: B(r)=1/A(r). -/
def dS_B (H r : Real) : Real := (dS_A H r)⁻¹

/-- PG slicing: ψ′(r) = Hr / A(r). -/
def dS_PG_psi' (H r : Real) : Real := (H * r) / (dS_A H r)

/-- X invariant (paper eq. (9)): X = 1/A - (ψ′)^2 / B. -/
def X_of (A B psip : Real) : Real := (1 / A) - (psip^2) / B

theorem deSitterPG_X_eq_one {H r : Real} (hA : dS_A H r ≠ 0) :
    X_of (dS_A H r) (dS_B H r) (dS_PG_psi' H r) = 1 := by
  unfold X_of dS_B dS_PG_psi'
  have hA0 : (1 - H ^ 2 * r ^ 2) ≠ (0 : Real) := by
    simpa [dS_A] using hA
  unfold dS_A
  field_simp [hA0]

/-- Derivatives needed for the substitution into K2_general. -/
def dS_A' (H r : Real) : Real := -2 * (H^2) * r
def dS_B' (H r : Real) : Real := (2 * (H^2) * r) / ((dS_A H r) ^ 2)

/-- ψ″ for PG slicing: d/dr (Hr/A) = H/A + 2H^3 r^2/A^2. -/
def dS_PG_psi'' (H r : Real) : Real :=
  (H / (dS_A H r)) + (2 * (H^3) * (r^2)) / ((dS_A H r) ^ 2)

/-- Core check: de Sitter static patch + PG slicing gives K^2 = 3H^2.
    (Proved by substituting into the *general* closed-form expression.) -/
theorem deSitterPG_K2_eq_threeH2_from_general
    {H r : Real} (hA : dS_A H r ≠ 0) (hr : r ≠ 0) :
    K2_general r
      (dS_A H r) (dS_B H r)
      (dS_A' H r) (dS_B' H r)
      (dS_PG_psi' H r) (dS_PG_psi'' H r)
    = 3 * H^2 := by
  unfold K2_general NumK2 DenK2

  -- r ≠ 0 side conditions
  have hr2 : r ^ 2 ≠ (0 : Real) := by
    exact pow_ne_zero 2 hr
  have hDen : (-4 : Real) * r ^ 2 ≠ 0 := by
    exact mul_ne_zero (by norm_num) hr2

  -- A ≠ 0 in expanded form
  have hA0 : (1 - H ^ 2 * r ^ 2) ≠ (0 : Real) := by
    simpa [dS_A] using hA

  -- core := -BB + AA*(psiP^2) = -1
  have hCoreEq :
      (- (dS_B H r) + (dS_A H r) * (dS_PG_psi' H r) ^ 2) = (-1 : Real) := by
    unfold dS_B dS_PG_psi'
    unfold dS_A
    field_simp [hA0]
    ring_nf

  -- simplify only the goal
  simp [hCoreEq] at ⊢

  -- clear denominators from A and r^2
  field_simp [hA0, hDen] at ⊢

  /-
    After the field_simp, Lean produces an inverse of a large polynomial P.
    That P is actually - (dS_A H r)^6.
  -/
  set P : Real :=
    (H ^ 2 * dS_A H r ^ 2 * r ^ 2 * 3 +
        (-(H ^ 4 * dS_A H r * r ^ 4 * 3) - H ^ 4 * dS_A H r ^ 2 * r ^ 4 * 3) +
      H ^ 6 * dS_A H r * r ^ 6 * 6 +
    H ^ 6 * r ^ 6 +
      (-(H ^ 8 * dS_A H r * r ^ 8 * 3) - H ^ 8 * r ^ 8 * 3) +
    H ^ 10 * r ^ 10 * 3 +
  (-(H ^ 12 * r ^ 12) - dS_A H r ^ 3)) with hPdef

  -- ✅ Correct exponent: P = - (dS_A)^6
  have hPsimp : P = - (dS_A H r) ^ 6 := by
    -- normalize P into a binomial expansion
    simp [P, dS_A]
    ring_nf

  have hPne : P ≠ 0 := by
    have hA6 : (dS_A H r) ^ 6 ≠ (0 : Real) := by
      exact pow_ne_zero 6 hA
    simpa [hPsimp] using (neg_ne_zero.mpr hA6)

  -- expose P in the goal and clear its inverse
  unfold dS_B dS_A' dS_B' dS_PG_psi' dS_PG_psi'' dS_A
  field_simp [hPne, hA0]
  ring_nf

/-- Leafwise constraint K^2 = K0^2 implies H^2 = K0^2/3 (selection rule). -/
theorem deSitter_selection_rule_from_general
    {H r K0 : Real} (hA : dS_A H r ≠ 0) (hr : r ≠ 0)
    (hConstraint :
      K2_general r
        (dS_A H r) (dS_B H r)
        (dS_A' H r) (dS_B' H r)
        (dS_PG_psi' H r) (dS_PG_psi'' H r)
      = K0^2) :
    H^2 = K0^2 / 3 := by
  have hK :
      K2_general r
        (dS_A H r) (dS_B H r)
        (dS_A' H r) (dS_B' H r)
        (dS_PG_psi' H r) (dS_PG_psi'' H r)
      = 3 * H^2 :=
    deSitterPG_K2_eq_threeH2_from_general (H := H) (r := r) hA hr
  have hk : 3 * H^2 = K0^2 := by
    nlinarith [hK, hConstraint]
  nlinarith [hk]

end
