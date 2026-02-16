import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Measure.OpenPos
import Mathlib.MeasureTheory.Measure.Typeclasses.Finite
import Mathlib.Tactic

namespace StaticSanityCheck.DeSitterWellPosed

open Filter MeasureTheory
open scoped Topology

set_option autoImplicit false
noncomputable section

abbrev E3 := EuclideanSpace ℝ (Fin 3)

/-- Horizon scale of the de Sitter static patch: `r_H = H⁻¹`. -/
def dSHorizonRadius (H : Real) : Real := H⁻¹

/-- Nested causal subdomains used in the IR definition:
`D_Θ(R) = ball(0, min(R, H⁻¹))`. -/
def dSDomain (H R : Real) : Set E3 :=
  Metric.ball (0 : E3) (min R (dSHorizonRadius H))

/-- Denominator functional of the leafwise average (weight `W = 1`). -/
def I0 (H R : Real) : Real := volume.real (dSDomain H R)

/-- Numerator functional of the leafwise average (weight `W = 1`). -/
def I1 (H R : Real) (f : E3 → Real) : Real :=
  ∫ x in dSDomain H R, f x

/-- Leafwise normalized average on `D_Θ(R)`. -/
def leafAverage (H R : Real) (f : E3 → Real) : Real :=
  I1 H R f / I0 H R

theorem dSHorizonRadius_pos {H : Real} (hH : 0 < H) :
    0 < dSHorizonRadius H := by
  exact inv_pos.mpr hH

/-- Well-posedness of the denominator: `I0 > 0` for positive `H` and positive cutoff `R`. -/
theorem I0_pos_of_pos {H R : Real} (hH : 0 < H) (hR : 0 < R) :
    0 < I0 H R := by
  have hMinPos : 0 < min R (dSHorizonRadius H) := by
    exact lt_min hR (dSHorizonRadius_pos hH)
  unfold I0 dSDomain
  have hMeasurePos :
      0 < volume (Metric.ball (0 : E3) (min R (dSHorizonRadius H))) :=
    Metric.measure_ball_pos (μ := volume) (x := (0 : E3)) hMinPos
  have hMeasureTop :
      volume (Metric.ball (0 : E3) (min R (dSHorizonRadius H))) ≠ ⊤ := by
    exact (MeasureTheory.measure_ball_lt_top
      (μ := volume) (x := (0 : E3)) (r := min R (dSHorizonRadius H))).ne
  exact ENNReal.toReal_pos hMeasurePos.ne' hMeasureTop

/-- In de Sitter static patch, each `D_Θ(R)` has finite volume (horizon-bounded domain). -/
theorem dSDomain_measure_lt_top (H R : Real) :
    volume (dSDomain H R) < ⊤ := by
  unfold dSDomain
  exact MeasureTheory.measure_ball_lt_top
    (μ := volume) (x := (0 : E3)) (r := min R (dSHorizonRadius H))

/-- Saturation of the causal domain at the horizon scale `H⁻¹`. -/
theorem dSDomain_eq_horizon_of_ge {H R : Real}
    (hR : dSHorizonRadius H ≤ R) :
    dSDomain H R = Metric.ball (0 : E3) (dSHorizonRadius H) := by
  unfold dSDomain
  simp [min_eq_right hR]

theorem I0_eq_horizon_of_ge {H R : Real}
    (hR : dSHorizonRadius H ≤ R) :
    I0 H R = I0 H (dSHorizonRadius H) := by
  unfold I0
  rw [dSDomain_eq_horizon_of_ge (H := H) (R := R) hR]
  unfold dSDomain
  simp

theorem I1_eq_horizon_of_ge {H R : Real} {f : E3 → Real}
    (hR : dSHorizonRadius H ≤ R) :
    I1 H R f = I1 H (dSHorizonRadius H) f := by
  unfold I1
  rw [dSDomain_eq_horizon_of_ge (H := H) (R := R) hR]
  unfold dSDomain
  simp

/-- Saturation of the normalized average at finite radius `R = H⁻¹`. -/
theorem leafAverage_eq_horizon_of_ge {H R : Real} {f : E3 → Real}
    (hR : dSHorizonRadius H ≤ R) :
    leafAverage H R f = leafAverage H (dSHorizonRadius H) f := by
  unfold leafAverage
  rw [I1_eq_horizon_of_ge (H := H) (R := R) (f := f) hR]
  rw [I0_eq_horizon_of_ge (H := H) (R := R) hR]

/-- Filter formulation of saturation: beyond `H⁻¹`, the average is eventually constant. -/
theorem leafAverage_eventuallyEq_horizon (H : Real) (f : E3 → Real) :
    (fun R : Real => leafAverage H R f) =ᶠ[atTop]
      fun _ : Real => leafAverage H (dSHorizonRadius H) f := by
  exact (eventually_ge_atTop (dSHorizonRadius H)).mono
    (fun R hR => leafAverage_eq_horizon_of_ge (H := H) (R := R) (f := f) hR)

/-- IR limit in de Sitter is saturated at finite radius: the `R → ∞` limit equals the horizon value. -/
theorem leafAverage_tendsto_horizon (H : Real) (f : E3 → Real) :
    Tendsto (fun R : Real => leafAverage H R f) atTop
      (nhds (leafAverage H (dSHorizonRadius H) f)) := by
  refine Tendsto.congr' ?_ tendsto_const_nhds
  exact (leafAverage_eventuallyEq_horizon H f).symm

end
