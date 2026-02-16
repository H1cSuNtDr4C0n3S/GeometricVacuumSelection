import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace StaticSanityCheck

set_option autoImplicit false
noncomputable section

/-- Numerator of the general closed-form K^2 expression (from `06_Test_deSitter.nb`). -/
def NumK2 (r AA BB AAp BBp psiP psiPP : Real) : Real :=
  (psiP ^ 2) *
      ( (-4) * (r ^ 2) * (BB ^ 2) * (AAp ^ 2)
        + 16 * (AA ^ 3) * BB * (psiP ^ 2)
        - 8 * (AA ^ 4) * (psiP ^ 4)
        + 4 * (r ^ 2) * AA * BB * AAp * (BBp + AAp * (psiP ^ 2))
        - (AA ^ 2) * ( 8 * (BB ^ 2) + (r ^ 2) * ((BBp + AAp * (psiP ^ 2)) ^ 2) )
      )
    + 4 * (r ^ 2) * AA * BB * psiP
        * ( AA * BBp + AAp * ( (-2) * BB + AA * (psiP ^ 2) ) )
        * psiPP
    - 4 * (r ^ 2) * (AA ^ 2) * (BB ^ 2) * (psiPP ^ 2)

/-- Denominator of the general closed-form K^2 expression (from `06_Test_deSitter.nb`). -/
def DenK2 (r AA BB psiP : Real) : Real :=
  4 * (r ^ 2) * AA * BB * ((-BB + AA * (psiP ^ 2)) ^ 3)

/-- General closed-form K^2 = NumK2 / DenK2. -/
def K2_general (r AA BB AAp BBp psiP psiPP : Real) : Real :=
  (NumK2 r AA BB AAp BBp psiP psiPP) / (DenK2 r AA BB psiP)

end
