"""
Benchmark & Comparison: Real-Gas CO2 Shock Tube (1000 Grid Points)
Compares exact analytical Helmholtz EoS vs Look-Up Table (LuT) 2D Bicubic Splines.
"""

import time
import pickle
import numpy as np
from pyshockflow import Driver, Config

def run_simulation(input_file, case_label):
    print(f"\n{'='*80}")
    print(f"  RUNNING: {case_label} ({input_file})")
    print(f"{'='*80}")
    
    config = Config(input_file)
    driver = Driver(config)
    
    start_time = time.perf_counter()
    driver.solve()
    elapsed_time = time.perf_counter() - start_time
    
    return elapsed_time

if __name__ == '__main__':
    print("=" * 80)
    print("  PYSHOCKFLOW PERFORMANCE COMPARISON: EXACT HELMHOLTZ vs 2D LuT SPLINES")
    print("  Test Case: Supercritical CO2 Shock Tube (1000 mesh points, Arabi Roe solver)")
    print("=" * 80)
    
    # 1. Run Exact Helmholtz
    t_exact = run_simulation('input_exact.ini', 'Case 1: Exact Analytical Helmholtz EoS (USE_LUT=no)')
    
    # 2. Run Look-Up Table
    t_lut = run_simulation('input_lut.ini', 'Case 2: 2D Bicubic Spline Look-Up Table (USE_LUT=yes)')
    
    # 3. Accuracy check from saved results
    with open('Results/real_exact_NX_500/Results.pik', 'rb') as f:
        res_exact = pickle.load(f)
        
    with open('Results/real_lut_NX_500/Results.pik', 'rb') as f:
        res_lut = pickle.load(f)
        
    p_exact = res_exact['Primitive']['Pressure'][:, -1]
    p_lut = res_lut['Primitive']['Pressure'][:, -1]
    rho_exact = res_exact['Primitive']['Density'][:, -1]
    rho_lut = res_lut['Primitive']['Density'][:, -1]
    u_exact = res_exact['Primitive']['Velocity'][:, -1]
    u_lut = res_lut['Primitive']['Velocity'][:, -1]
    
    err_p = np.max(np.abs(p_lut - p_exact)) / np.max(np.abs(p_exact))
    err_rho = np.max(np.abs(rho_lut - rho_exact)) / np.max(np.abs(rho_exact))
    err_u = np.max(np.abs(u_lut - u_exact)) / (np.max(np.abs(u_exact)) + 1e-12)
    
    speedup = t_exact / t_lut
    
    # 4. Summary Table
    print("\n" + "=" * 80)
    print("                           PERFORMANCE & ACCURACY SUMMARY")
    print("=" * 80)
    print(f"  Exact Helmholtz EoS Runtime (1000 pts) : {t_exact:8.4f} s")
    print(f"  LuT 2D Splines Runtime      (1000 pts) : {t_lut:8.4f} s")
    print(f"  Total Speedup                          : {speedup:8.2f} x")
    print("-" * 80)
    print("  Solution Accuracy (Max Relative Error vs Exact at t_final):")
    print(f"    • Pressure Error                     : {err_p:8.4e} ({err_p*100:.3f}%)")
    print(f"    • Density Error                      : {err_rho:8.4e} ({err_rho*100:.3f}%)")
    print(f"    • Velocity Error                     : {err_u:8.4e} ({err_u*100:.3f}%)")
    print("=" * 80)
    print("  To visualize comparison plots, run: python plot_comparison.py\n")
