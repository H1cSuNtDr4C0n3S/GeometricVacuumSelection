#!/usr/bin/env python3
"""
Quantitative analysis of causal domain sensitivity for vacuum selection.
Based on "Geometric Vacuum Selection from Constrained Spacetime Foliations"
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Set up plotting style
plt.rcParams['font.size'] = 12
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# ==============================================================================
# CASE 1: MINKOWSKI WITH LINEAR TILT (Section 4.1)
# ==============================================================================

def K_squared_minkowski(r, v):
    """
    K² for Minkowski with linear tilt (Eq. 23)
    K²(r) = 2v²/[r²(1-v²)]
    """
    if r == 0:
        return np.inf
    return 2 * v**2 / (r**2 * (1 - v**2))

def volume_element_minkowski(r):
    """Volume element for spherical shell: 4πr²"""
    return 4 * np.pi * r**2

def average_K_squared_minkowski(R_max, v):
    """
    Compute ⟨K²⟩ for domain [0, R_max]
    ⟨K²⟩ = ∫ K² dV / ∫ dV
    """
    # I₁ = ∫₀^R_max K²(r) · 4πr² dr
    def integrand_I1(r):
        if r < 1e-10:  # Avoid singularity at r=0
            return 0
        return K_squared_minkowski(r, v) * volume_element_minkowski(r)
    
    I1, _ = quad(integrand_I1, 1e-10, R_max, limit=100)
    
    # I₀ = ∫₀^R_max 4πr² dr = (4π/3)R_max³
    I0 = (4 * np.pi / 3) * R_max**3
    
    return I1 / I0

# Calculate for different domain sizes
v = 0.5  # velocity parameter
R_values = np.logspace(-1, 3, 50)  # 0.1 to 1000
avg_K2_values = [average_K_squared_minkowski(R, v) for R in R_values]

# ==============================================================================
# CASE 2: DE SITTER (Section 4.3)
# ==============================================================================

def K_squared_de_sitter(H):
    """
    K² for de Sitter in Painlevé-Gullstrand coordinates (Eq. 28)
    K² = 3H² (constant everywhere)
    """
    return 3 * H**2

def average_K_squared_de_sitter(R_max, H):
    """
    For de Sitter, K² is constant, so ⟨K²⟩ = K² regardless of domain
    """
    return K_squared_de_sitter(H)

# ==============================================================================
# CASE 3: HYBRID MODEL (Minkowski core + de Sitter exterior)
# ==============================================================================

def K_squared_hybrid(r, v, H, r_transition, smoothness=0.1):
    """
    Hybrid model: Minkowski-like for r < r_transition, de Sitter for r > r_transition
    Use smooth interpolation
    """
    K2_minkowski = K_squared_minkowski(r, v) if r > 1e-10 else 0
    K2_dS = K_squared_de_sitter(H)
    
    # Smooth transition using tanh
    transition_factor = 0.5 * (1 + np.tanh((r - r_transition) / (smoothness * r_transition)))
    
    return (1 - transition_factor) * K2_minkowski + transition_factor * K2_dS

def average_K_squared_hybrid(R_max, v, H, r_transition):
    """Compute ⟨K²⟩ for hybrid model"""
    def integrand_I1(r):
        return K_squared_hybrid(r, v, H, r_transition) * volume_element_minkowski(r)
    
    I1, _ = quad(integrand_I1, 0, R_max, limit=100)
    I0 = (4 * np.pi / 3) * R_max**3
    
    return I1 / I0

# ==============================================================================
# PLOTTING (MODIFICATO PER SALVARE 4 FILE DISTINTI)
# ==============================================================================

print("Generating plots...")

# ----- Plot 1: Minkowski - Domain Size Dependence -----
plt.figure(figsize=(10, 7))
plt.loglog(R_values, avg_K2_values, 'b-', linewidth=2, label=f'v = {v}')
plt.loglog(R_values, 1/R_values**2, 'r--', linewidth=1.5, label='$1/R^2$ (scaling prediction)')
plt.xlabel('Domain Radius $R_{max}$ (arbitrary units)')
plt.ylabel('$\\langle K^2 \\rangle$ (arbitrary units)')
plt.title('Minkowski with Linear Tilt: Domain Sensitivity')
plt.legend()
plt.text(0.05, 0.95, 'Infrared limit: $\\langle K^2 \\rangle \\sim 1/R^2$', 
         transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('plot_1_minkowski_sensitivity.png', dpi=300)
plt.close() # Chiude la figura per liberare memoria
print("✓ Saved: plot_1_minkowski_sensitivity.png")

# ----- Plot 2: Different velocities -----
plt.figure(figsize=(10, 7))
for v_test in [0.1, 0.3, 0.5, 0.7, 0.9]:
    avg_K2_v = [average_K_squared_minkowski(R, v_test) for R in R_values]
    plt.loglog(R_values, avg_K2_v, linewidth=2, label=f'v = {v_test}')

plt.xlabel('Domain Radius $R_{max}$ (arbitrary units)')
plt.ylabel('$\\langle K^2 \\rangle$ (arbitrary units)')
plt.title('Velocity Dependence in Minkowski')
plt.legend()

plt.tight_layout()
plt.savefig('plot_2_velocity_dependence.png', dpi=300)
plt.close()
print("✓ Saved: plot_2_velocity_dependence.png")

# ----- Plot 3: Hybrid Model - Transition Effects -----
plt.figure(figsize=(10, 7))
H = 1.0  # Hubble parameter
v_hybrid = 0.5
r_trans_values = [1.0, 5.0, 10.0, 50.0]
R_hybrid = np.logspace(0, 2.5, 40)

for r_trans in r_trans_values:
    avg_K2_hybrid = [average_K_squared_hybrid(R, v_hybrid, H, r_trans) for R in R_hybrid]
    plt.semilogx(R_hybrid, avg_K2_hybrid, linewidth=2, label=f'$r_{{trans}}$ = {r_trans}')

plt.axhline(y=3*H**2, color='k', linestyle='--', linewidth=1.5, label='Pure de Sitter limit')
plt.xlabel('Domain Radius $R_{max}$ (arbitrary units)')
plt.ylabel('$\\langle K^2 \\rangle$ (arbitrary units)')
plt.title('Hybrid Model: Transition Scale Effects')
plt.legend()

plt.tight_layout()
plt.savefig('plot_3_hybrid_transition.png', dpi=300)
plt.close()
print("✓ Saved: plot_3_hybrid_transition.png")

# ----- Plot 4: Convergence Analysis -----
plt.figure(figsize=(10, 7))

# Calculate relative change in ⟨K²⟩ as domain increases
R_convergence = np.logspace(0, 3, 100)
avg_K2_conv = np.array([average_K_squared_minkowski(R, 0.5) for R in R_convergence])

# Relative change: |⟨K²⟩(R) - ⟨K²⟩(R_max)| / ⟨K²⟩(R_max)
reference_value = avg_K2_conv[-1]
relative_change = np.abs(avg_K2_conv - reference_value) / reference_value

plt.loglog(R_convergence, relative_change * 100, 'g-', linewidth=2)
plt.axhline(y=10, color='r', linestyle='--', linewidth=1.5, label='10% threshold')
plt.axhline(y=1, color='orange', linestyle='--', linewidth=1.5, label='1% threshold')
plt.xlabel('Domain Radius $R_{max}$ (arbitrary units)')
plt.ylabel('Relative Change (%)')
plt.title('Convergence to Infrared Limit')
plt.legend()
plt.text(0.05, 0.95, 'Minkowski case (v=0.5)', 
         transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

plt.tight_layout()
plt.savefig('plot_4_convergence_analysis.png', dpi=300)
plt.close()
print("✓ Saved: plot_4_convergence_analysis.png")

# ==============================================================================
# QUANTITATIVE RESULTS TABLE (INVARIATO)
# ==============================================================================

print("\n" + "="*80)
print("QUANTITATIVE SENSITIVITY ANALYSIS")
print("="*80)

print("\n1. MINKOWSKI WITH LINEAR TILT (v = 0.5)")
print("-" * 80)
print(f"{'Domain Radius':<20} {'⟨K²⟩':<20} {'Relative to R=1000':<25}")
print("-" * 80)

R_test = [1, 10, 100, 1000]
K2_values = [average_K_squared_minkowski(R, 0.5) for R in R_test]
K2_ref = K2_values[-1]

for R, K2 in zip(R_test, K2_values):
    rel_diff = abs(K2 - K2_ref) / K2_ref * 100
    print(f"{R:<20} {K2:<20.6f} {rel_diff:<25.2f}%")

print(f"\nScaling behavior: ⟨K²⟩ ∝ R^α where α ≈ {np.log(K2_values[-1]/K2_values[0]) / np.log(R_test[-1]/R_test[0]):.3f}")
print("Expected: α = -2 for infrared limit")

print("\n2. DE SITTER CASE (H = 1.0)")
print("-" * 80)
print(f"K² = 3H² = {3 * 1.0**2:.6f} (constant, independent of domain)")
print("Domain sensitivity: NONE (homogeneous spacetime)")

print("\n3. HYBRID MODEL SENSITIVITY")
print("-" * 80)
print(f"{'Transition r':<15} {'R_max=10':<15} {'R_max=50':<15} {'R_max=100':<15}")
print("-" * 80)

for r_trans in [1.0, 5.0, 10.0]:
    K2_10 = average_K_squared_hybrid(10, 0.5, 1.0, r_trans)
    K2_50 = average_K_squared_hybrid(50, 0.5, 1.0, r_trans)
    K2_100 = average_K_squared_hybrid(100, 0.5, 1.0, r_trans)
    print(f"{r_trans:<15.1f} {K2_10:<15.6f} {K2_50:<15.6f} {K2_100:<15.6f}")

print("\n4. CONVERGENCE CRITERIA")
print("-" * 80)
print("For Minkowski case (v=0.5), domain radius required for:")

for threshold in [10, 5, 1, 0.1]:
    # Find R where relative change < threshold%
    for R in R_convergence:
        K2 = average_K_squared_minkowski(R, 0.5)
        rel_change = abs(K2 - reference_value) / reference_value * 100
        if rel_change < threshold:
            print(f"  {threshold}% accuracy: R > {R:.1f}")
            break

print("\n5. VELOCITY PARAMETER SENSITIVITY")
print("-" * 80)
print(f"{'Velocity v':<15} {'⟨K²⟩ at R=10':<20} {'⟨K²⟩ at R=100':<20}")
print("-" * 80)

for v_test in [0.1, 0.3, 0.5, 0.7, 0.9]:
    K2_10 = average_K_squared_minkowski(10, v_test)
    K2_100 = average_K_squared_minkowski(100, v_test)
    print(f"{v_test:<15.1f} {K2_10:<20.6f} {K2_100:<20.6f}")

print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)
print("""
1. MINKOWSKI SENSITIVITY: HIGH
   - ⟨K²⟩ ∝ 1/R^2 in infrared limit
   - Requires R >> transition scale for convergence
   - Domain choice CRITICAL for vacuum selection

2. DE SITTER SENSITIVITY: NONE
   - K² constant everywhere
   - Result independent of domain size
   - Vacuum selection ROBUST

3. HYBRID MODELS: MODERATE
   - Sensitive to R_max/r_transition ratio
   - Converges to dominant component when R >> r_transition
   - Physical interpretation: dominant geometry determines vacuum

4. PRACTICAL IMPLICATION:
   - For inhomogeneous spacetimes: domain choice affects vacuum selection
   - For homogeneous spacetimes: domain choice irrelevant
   - Physical cutoff (e.g., cosmological horizon) may naturally define D_Θ
""")

print("\n✓ Analysis complete! 4 Plots saved in the current directory.")