import Mathlib.Topology.Algebra.Ring.Real
import Mathlib.Algebra.Order.Ring.Abs

namespace StaticSanityCheck

open Filter
set_option autoImplicit false
noncomputable section

/-- de Sitter static patch metric functions: A(r)=1-H^2 r^2, B(r)=1/A(r). -/
def dS_A (H r : Real) : Real := 1 - (H^2) * (r^2)
def dS_B (H r : Real) : Real := (dS_A H r)⁻¹

/-- Painlevé–Gullstrand slicing in de Sitter static patch (paper eq. (27)):
    ψ′(r) = Hr / (1 - H^2 r^2) = Hr / A(r). -/
def dS_PG_psi' (H r : Real) : Real := (H * r) / (dS_A H r)

/-- Kinetic invariant for Θ = t + ψ(r) (paper eq. (9)):
    X = 1/A - (ψ′)^2 / B. -/
def X_of (A B psip : Real) : Real := (1 / A) - (psip^2) / B

/-- In de Sitter static patch with PG slicing, X = 1 on the admissible region A ≠ 0. -/
theorem deSitterPG_X_eq_one {H r : Real} (hA : dS_A H r ≠ 0) :
    X_of (dS_A H r) (dS_B H r) (dS_PG_psi' H r) = 1 := by
  unfold X_of dS_A dS_B dS_PG_psi'
  have hA' : (1 - H ^ 2 * r ^ 2) ≠ 0 := by
    simpa using hA
  field_simp [hA']
  -- crucial: expand A inside the remaining polynomial goal
  simp [dS_A] at *
  ring_nf

/-- Closed-form result stated in the paper (eq. (28)): K^2 = 3 H^2. -/
def deSitterPG_K2 (H : Real) : Real := 3 * H^2

/-- Since K^2 is constant on the causal domain of the static patch, the normalized
    leafwise average equals the local value (paper eq. (29)). -/
def deSitterPG_leafwiseAverage (H : Real) : Real := deSitterPG_K2 H

theorem deSitterPG_leafwiseAverage_eq (H : Real) :
    deSitterPG_leafwiseAverage H = 3 * H^2 := by
  rfl

/-- Leafwise constraint selection rule (paper eq. (30)):
    if <K^2> = K0^2 then H^2 = K0^2 / 3. -/
theorem deSitter_selection_rule {H K0 : Real}
    (h : deSitterPG_leafwiseAverage H = K0^2) :
    H^2 = K0^2 / 3 := by
  have : 3 * H^2 = K0^2 := by
    simpa [deSitterPG_leafwiseAverage, deSitterPG_K2] using h
  nlinarith

end
