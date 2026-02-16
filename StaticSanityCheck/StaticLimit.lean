import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

namespace StaticSanityCheck

set_option autoImplicit false
noncomputable section

def XInvariant (A B psi' : Real) : Real :=
  1 / A - psi' ^ 2 / B

@[simp] lemma XInvariant_static (A B : Real) :
    XInvariant A B 0 = 1 / A := by
  simp [XInvariant]

def K2Denominator (A B r : Real) : Real :=
  r ^ 2 * (-B) ^ 3 * (A * B)

lemma K2Denominator_ne_zero
    (A B r : Real) (hA : A ≠ 0) (hB : B ≠ 0) (hr : r ≠ 0) :
    K2Denominator A B r ≠ 0 := by
  have hr2 : r ^ 2 ≠ 0 := by simpa using (pow_ne_zero 2 hr)
  have hBpow : (-B) ^ 3 ≠ 0 := by
    simpa using pow_ne_zero 3 (neg_ne_zero.mpr hB)
  have hAB : A * B ≠ 0 := mul_ne_zero hA hB
  have hprod : r ^ 2 * (-B) ^ 3 * (A * B) ≠ 0 :=
    mul_ne_zero (mul_ne_zero hr2 hBpow) hAB
  simpa [K2Denominator] using hprod

def K2Invariant (A B Aprime Bprime psi' psi'' r : Real) : Real :=
  (psi' ^ 2 * (Aprime ^ 2 + Bprime ^ 2) + psi'' ^ 2 * r ^ 2) /
    K2Denominator A B r

lemma K2Invariant_static
    (A B Aprime Bprime r : Real)
    (hA : A ≠ 0) (hB : B ≠ 0) (hr : r ≠ 0) :
    K2Invariant A B Aprime Bprime 0 0 r = 0 := by
  have _ : K2Denominator A B r ≠ 0 := K2Denominator_ne_zero A B r hA hB hr
  simp [K2Invariant, K2Denominator]

lemma static_sanity_check
    (A B Aprime Bprime r : Real)
    (hA : A ≠ 0) (hB : B ≠ 0) (hr : r ≠ 0) :
    XInvariant A B 0 = 1 / A ∧
      K2Invariant A B Aprime Bprime 0 0 r = 0 :=
  ⟨XInvariant_static A B, K2Invariant_static A B Aprime Bprime r hA hB hr⟩

end
