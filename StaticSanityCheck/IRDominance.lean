import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Measure.Lebesgue.VolumeOfBalls
import Mathlib.Order.Filter.AtTopBot.Field
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

namespace StaticSanityCheck

open Filter MeasureTheory
open scoped Topology

set_option autoImplicit false
noncomputable section

abbrev E3 := EuclideanSpace ℝ (Fin 3)

/--
Abstract IR lemma: if a numerator becomes eventually constant while the denominator
tends to `+∞`, the normalized ratio goes to `0`.
-/
theorem tendsto_div_atTop_zero_of_eventuallyEq_const
    {N D : ℝ → ℝ} {C : ℝ}
    (hN : ∀ᶠ R in atTop, N R = C)
    (hD : Tendsto D atTop atTop) :
    Tendsto (fun R => N R / D R) atTop (nhds 0) := by
  have hInv : Tendsto (fun R => (D R)⁻¹) atTop (nhds 0) :=
    tendsto_inv_atTop_zero.comp hD
  have hConstDiv : Tendsto (fun R => C / D R) atTop (nhds 0) := by
    have hMul : Tendsto (fun R => C * (D R)⁻¹) atTop (nhds (C * 0)) :=
      tendsto_const_nhds.mul hInv
    simpa [div_eq_mul_inv] using hMul
  refine Tendsto.congr' ?_ hConstDiv
  exact hN.mono (fun R hR => by simp [hR])

/-- In 3D Euclidean space, the volume of `ball 0 R` diverges as `R → +∞`. -/
theorem tendsto_volume_real_ball_fin_three_atTop :
    Tendsto (fun R : ℝ => volume.real (Metric.ball (0 : E3) R)) atTop atTop := by
  have hpow : Tendsto (fun R : ℝ => R ^ 3) atTop atTop :=
    tendsto_pow_atTop (by norm_num : (3 : ℕ) ≠ 0)
  have hconst_pos : 0 < (Real.pi * 4 / 3) := by
    nlinarith [Real.pi_pos]
  have hconst_nonneg : 0 ≤ (Real.pi * 4 / 3) := le_of_lt hconst_pos
  have hmul : Tendsto (fun R : ℝ => (R ^ 3) * (Real.pi * 4 / 3)) atTop atTop :=
    hpow.atTop_mul_const hconst_pos
  refine Tendsto.congr' ?_ hmul
  filter_upwards [eventually_gt_atTop (1 : ℝ)] with R hR
  have hRpos : 0 < R := lt_trans zero_lt_one hR
  have hRle : 0 ≤ R := le_of_lt hRpos
  simp [Measure.real, EuclideanSpace.volume_ball_fin_three, hRle, hconst_nonneg]

/--
IR dominance (compact/local source, 3D version): if `f` is supported inside a
fixed finite ball, then its normalized integral over `ball 0 R` tends to `0`.
-/
theorem tendsto_ballIntegral_div_volume_real_zero_of_tsupport_subset_closedBall
    {f : E3 → ℝ} {R0 : ℝ}
    (hSupp : tsupport f ⊆ Metric.closedBall (0 : E3) R0) :
    Tendsto
      (fun R : ℝ =>
        (∫ x in Metric.ball (0 : E3) R, f x) /
          volume.real (Metric.ball (0 : E3) R))
      atTop (nhds 0) := by
  let C : ℝ := ∫ x, f x
  have hNumEventuallyConst :
      ∀ᶠ R in atTop, (∫ x in Metric.ball (0 : E3) R, f x) = C := by
    filter_upwards [eventually_gt_atTop R0] with R hR
    have hZeroOutside : ∀ x, x ∉ Metric.ball (0 : E3) R → f x = 0 := by
      intro x hxBall
      have hge : R ≤ dist x 0 := by
        have hnot : ¬ dist x 0 < R := by
          simpa [Metric.mem_ball] using hxBall
        exact le_of_not_gt hnot
      have hR0lt : R0 < dist x 0 := lt_of_lt_of_le hR hge
      have hxNotClosed : x ∉ Metric.closedBall (0 : E3) R0 := by
        simpa [Metric.mem_closedBall, not_le] using hR0lt
      have hxNotTsupport : x ∉ tsupport f := fun hxT => hxNotClosed (hSupp hxT)
      exact image_eq_zero_of_notMem_tsupport hxNotTsupport
    calc
      (∫ x in Metric.ball (0 : E3) R, f x)
          = ∫ x, f x := by
              exact MeasureTheory.setIntegral_eq_integral_of_forall_compl_eq_zero hZeroOutside
      _ = C := by rfl
  exact tendsto_div_atTop_zero_of_eventuallyEq_const
    (N := fun R => ∫ x in Metric.ball (0 : E3) R, f x)
    (D := fun R => volume.real (Metric.ball (0 : E3) R))
    (C := C)
    hNumEventuallyConst
    tendsto_volume_real_ball_fin_three_atTop

end
