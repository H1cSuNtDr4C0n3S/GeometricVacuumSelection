import Mathlib.Topology.Algebra.Ring.Real
import Mathlib.Algebra.Order.Ring.Abs

namespace StaticSanityCheck

open Filter

set_option autoImplicit false

noncomputable section

/--
`minkowskiInfraredAverage v R` stores the closed-form expression obtained by
evaluating the normalized radial average of `K^2` for the linear tilt
`psi(r) = v * r` inside Minkowski spacetime, truncated at radius `R`.
The numerator scales like the radial volume measure while the denominator
captures the `R^3` growth of the domain, leaving an overall `R^-2` decay.
-/
def minkowskiInfraredAverage (v R : Real) : Real :=
  6 * v ^ 2 / ((1 - v ^ 2) * R ^ 2)

lemma minkowskiInfraredAverage_eq_const_div (v R : Real) :
    minkowskiInfraredAverage v R =
      (6 * v ^ 2 / (1 - v ^ 2)) / R ^ 2 := by
  unfold minkowskiInfraredAverage
  calc
    6 * v ^ 2 / ((1 - v ^ 2) * R ^ 2)
        = 6 * v ^ 2 * ((1 - v ^ 2) * R ^ 2)⁻¹ := by
          simp [div_eq_mul_inv]
    _ = 6 * v ^ 2 * ((1 - v ^ 2)⁻¹ * (R ^ 2)⁻¹) := by
          simp [mul_comm, mul_left_comm, mul_assoc]
    _ = (6 * v ^ 2 / (1 - v ^ 2)) / R ^ 2 := by
          simp [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc]

/-- The infrared average in Minkowski spacetime decays quadratically with the
cutoff radius, forcing the limit to vanish for any subluminal tilt. -/
theorem minkowskiInfraredAverage_tendsto_zero {v : Real}
    (hv : |v| < 1) :
    Tendsto (fun R : Real => minkowskiInfraredAverage v R) atTop (nhds 0) := by
  have hv_sq : v ^ 2 < 1 := by
    have := (sq_lt_sq (a := v) (b := (1 : Real))).2 (by simpa [abs_one] using hv)
    simpa using this
  set C := 6 * v ^ 2 / (1 - v ^ 2)
  have h_inv : Tendsto (fun R : Real => (1 : Real) / R) atTop (nhds 0) := by
    simpa [one_div] using tendsto_inv_atTop_zero
  have h_inv_sq :
      Tendsto (fun R : Real => (1 : Real) / R ^ 2) atTop (nhds 0) := by
    have := h_inv.mul h_inv
    simpa [pow_two, one_div, mul_comm, mul_left_comm, mul_assoc] using this
  have hconst : Tendsto (fun _ : Real => C) atTop (nhds C) := tendsto_const_nhds
  have hdiv : Tendsto (fun R : Real => C / R ^ 2) atTop (nhds 0) := by
    have hmul := hconst.mul h_inv_sq
    simpa [div_eq_mul_inv] using hmul
  have : Tendsto (fun R : Real => (6 * v ^ 2 / (1 - v ^ 2)) / R ^ 2) atTop (nhds 0) := by
    simpa [C] using hdiv
  simpa [minkowskiInfraredAverage_eq_const_div] using this

end
