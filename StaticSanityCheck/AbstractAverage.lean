import Mathlib.MeasureTheory.Integral.Average
import Mathlib.MeasureTheory.Integral.IntegrableOn
import Mathlib.Tactic

namespace StaticSanityCheck

open Filter MeasureTheory
open scoped Topology

set_option autoImplicit false
noncomputable section

variable {alpha iota : Type*} [MeasurableSpace alpha] {mu : Measure alpha}

/-- Normalized average on a domain `s`, written without notation. -/
def leafAverage (mu : Measure alpha) (s : Set alpha) (f : alpha -> Real) : Real :=
  MeasureTheory.average (Measure.restrict mu s) f

/-- If a real-valued quantity is constant on a measurable leaf/domain, then its
normalized average equals that constant. -/
theorem leafAverage_eq_const_of_forall_eq
    {s : Set alpha} {f : alpha -> Real} {c : Real}
    (hs : MeasurableSet s)
    (hs0 : Not (mu s = 0))
    (hFin : IsFiniteMeasure (Measure.restrict mu s))
    (hconst : forall x, s x -> f x = c) :
    leafAverage mu s f = c := by
  let nu : Measure alpha := Measure.restrict mu s
  haveI : IsFiniteMeasure nu := by simpa [nu] using hFin
  have hnu0 : Not (nu = 0) := by
    intro h0
    apply hs0
    have h0u : nu Set.univ = 0 := by simp [h0]
    simpa [nu] using h0u
  letI : NeZero nu := by exact { out := hnu0 }
  have hconstAe : Filter.EventuallyEq (MeasureTheory.ae nu) f (fun _ : alpha => c) := by
    have hmem : Filter.Eventually (fun x => s x) (MeasureTheory.ae nu) := by
      simpa [nu] using (MeasureTheory.ae_restrict_mem hs)
    have hforall : Filter.Eventually (fun x => s x -> f x = c) (MeasureTheory.ae nu) :=
      Filter.Eventually.of_forall (fun x hx => hconst x hx)
    exact (hmem.and hforall).mono (fun x hx => hx.2 hx.1)
  have hcongr :
      MeasureTheory.average nu f = MeasureTheory.average nu (fun _ : alpha => c) :=
    MeasureTheory.average_congr hconstAe
  have hconstAvg : MeasureTheory.average nu (fun _ : alpha => c) = c := by
    exact MeasureTheory.average_const nu c
  simpa [leafAverage, nu] using hcongr.trans hconstAvg

/-- Quantitative stability of a normalized average:
if `f` stays within `eps` from `L` on `s`, then the normalized average on `s`
also stays within `eps` from `L`. -/
theorem norm_leafAverage_sub_le_of_forall_norm_sub_le
    {s : Set alpha} {f : alpha -> Real} {L eps : Real}
    (hs : MeasurableSet s)
    (hs0 : Not (mu s = 0))
    (hFin : IsFiniteMeasure (Measure.restrict mu s))
    (hInt : IntegrableOn f s mu)
    (hBound : forall x, s x -> norm (f x - L) <= eps) :
    norm (leafAverage mu s f - L) <= eps := by
  let nu : Measure alpha := Measure.restrict mu s
  haveI : IsFiniteMeasure nu := by simpa [nu] using hFin
  have hIntNu : Integrable f nu := by simpa [nu] using hInt
  have hnuUniv0 : Not (nu Set.univ = 0) := by
    intro h0
    apply hs0
    simpa [nu] using h0
  have hnuUnivTop : Not (nu Set.univ = ⊤) := by
    exact MeasureTheory.measure_ne_top nu Set.univ
  have hnuRealPos : 0 < nu.real Set.univ := by
    exact ENNReal.toReal_pos hnuUniv0 hnuUnivTop
  have hAvgMul :
      nu.real Set.univ * MeasureTheory.average nu f = MeasureTheory.integral nu f := by
    simpa [smul_eq_mul] using (MeasureTheory.measure_smul_average nu f)
  have hConstMul :
      nu.real Set.univ * MeasureTheory.average nu (fun _ : alpha => L) =
        MeasureTheory.integral nu (fun _ : alpha => L) := by
    simpa [smul_eq_mul] using (MeasureTheory.measure_smul_average nu (fun _ : alpha => L))
  have hAvgConst : MeasureTheory.average nu (fun _ : alpha => L) = L := by
    have hnu0 : Not (nu = 0) := by
      intro h0
      apply hs0
      have h0u : nu Set.univ = 0 := by simp [h0]
      simpa [nu] using h0u
    letI : NeZero nu := by exact { out := hnu0 }
    exact MeasureTheory.average_const nu L
  have hSmul :
      nu.real Set.univ * (leafAverage mu s f - L) =
        MeasureTheory.integral nu (fun x => f x - L) := by
    calc
      nu.real Set.univ * (leafAverage mu s f - L)
          = nu.real Set.univ * (MeasureTheory.average nu f - L) := by
              rfl
      _ = nu.real Set.univ * MeasureTheory.average nu f - nu.real Set.univ * L := by ring
      _ = nu.real Set.univ * MeasureTheory.average nu f -
            nu.real Set.univ * MeasureTheory.average nu (fun _ : alpha => L) := by
            simp [hAvgConst]
      _ = MeasureTheory.integral nu f - MeasureTheory.integral nu (fun _ : alpha => L) := by
            rw [hAvgMul, hConstMul]
      _ = MeasureTheory.integral nu (fun x => f x - L) := by
            simpa using (MeasureTheory.integral_sub hIntNu (integrable_const L)).symm
  have hBoundAe : Filter.Eventually (fun x => norm (f x - L) <= eps) (MeasureTheory.ae nu) := by
    have hmem : Filter.Eventually (fun x => s x) (MeasureTheory.ae nu) := by
      simpa [nu] using (MeasureTheory.ae_restrict_mem hs)
    have hforall : Filter.Eventually (fun x => s x -> norm (f x - L) <= eps) (MeasureTheory.ae nu) :=
      Filter.Eventually.of_forall (fun x hx => hBound x hx)
    exact (hmem.and hforall).mono (fun x hx => hx.2 hx.1)
  have hNormInt :
      norm (MeasureTheory.integral nu (fun x => f x - L)) <= eps * nu.real Set.univ := by
    simpa using (MeasureTheory.norm_integral_le_of_norm_le_const
      (f := fun x => f x - L) (C := eps) hBoundAe)
  have hNormSmul :
      norm (nu.real Set.univ * (leafAverage mu s f - L)) <= eps * nu.real Set.univ := by
    simpa [hSmul] using hNormInt
  have hMul :
      norm (leafAverage mu s f - L) * nu.real Set.univ <= eps * nu.real Set.univ := by
    have h1 : abs (nu.real Set.univ) * norm (leafAverage mu s f - L) <= eps * nu.real Set.univ := by
      simpa [norm_mul] using hNormSmul
    have habs : abs (nu.real Set.univ) = nu.real Set.univ := abs_of_nonneg ENNReal.toReal_nonneg
    have h2 : nu.real Set.univ * norm (leafAverage mu s f - L) <= eps * nu.real Set.univ := by
      simpa [habs] using h1
    simpa [mul_comm, mul_left_comm, mul_assoc] using h2
  exact le_of_mul_le_mul_right hMul hnuRealPos

/-- Abstract average lemma:
if `K2 i` converges uniformly to `Klim` on each domain `D i` (eventually),
then the normalized averages converge to `Klim`. -/
theorem tendsto_leafAverage_of_uniform_on_domains
    {l : Filter iota} {D : iota -> Set alpha} {K2 : iota -> alpha -> Real} {Klim : Real}
    (hMeas : Filter.Eventually (fun i => MeasurableSet (D i)) l)
    (hNonzero : Filter.Eventually (fun i => Not (mu (D i) = 0)) l)
    (hFinite : Filter.Eventually (fun i => IsFiniteMeasure (Measure.restrict mu (D i))) l)
    (hInt : Filter.Eventually (fun i => IntegrableOn (K2 i) (D i) mu) l)
    (hUniform :
      forall eps, eps > 0 ->
        Filter.Eventually (fun i => forall x, D i x -> norm (K2 i x - Klim) <= eps) l) :
    Tendsto (fun i => leafAverage mu (D i) (K2 i)) l (nhds Klim) := by
  refine Metric.tendsto_nhds.2 ?_
  intro eps heps
  have heps2 : 0 < eps / 2 := by linarith
  filter_upwards [hMeas, hNonzero, hFinite, hInt, hUniform (eps / 2) heps2]
    with i hMeasI h0 hFinI hIntI hUnifI
  have hle :
      norm (leafAverage mu (D i) (K2 i) - Klim) <= eps / 2 :=
    norm_leafAverage_sub_le_of_forall_norm_sub_le
      (mu := mu) (s := D i) (f := K2 i) (L := Klim) (eps := eps / 2)
      hMeasI h0 hFinI hIntI (fun x hx => hUnifI x hx)
  have hlt : norm (leafAverage mu (D i) (K2 i) - Klim) < eps := by
    linarith
  simpa [Real.dist_eq] using hlt

/-- Constant-on-domain corollary: if `K2` is eventually constant on each domain,
then the normalized averages converge to that constant. -/
theorem tendsto_leafAverage_of_eventually_const
    {l : Filter iota} {D : iota -> Set alpha} {K2 : iota -> alpha -> Real} {Klim : Real}
    (hMeas : Filter.Eventually (fun i => MeasurableSet (D i)) l)
    (hNonzero : Filter.Eventually (fun i => Not (mu (D i) = 0)) l)
    (hFinite : Filter.Eventually (fun i => IsFiniteMeasure (Measure.restrict mu (D i))) l)
    (hConst : Filter.Eventually (fun i => forall x, D i x -> K2 i x = Klim) l) :
    Tendsto (fun i => leafAverage mu (D i) (K2 i)) l (nhds Klim) := by
  refine Tendsto.congr' ?_ tendsto_const_nhds
  filter_upwards [hMeas, hNonzero, hFinite, hConst] with i hMeasI h0 hFinI hConstI
  exact (leafAverage_eq_const_of_forall_eq
    (mu := mu) (s := D i) (f := K2 i) (c := Klim) hMeasI h0 hFinI hConstI).symm

end
end StaticSanityCheck
