"""
Plot comparison between Exact Helmholtz EoS and Look-Up Table (LuT) 2D Bicubic Splines.
"""

import pickle
import os
import matplotlib.pyplot as plt
import numpy as np

def plot_comparison():
    file_exact = 'Results/real_exact_NX_500/Results.pik'
    file_lut = 'Results/real_lut_NX_500/Results.pik'
    
    if not os.path.exists(file_exact) or not os.path.exists(file_lut):
        print("Results files not found! Please run 'python main.py' first.")
        return
        
    with open(file_exact, 'rb') as f:
        res_exact = pickle.load(f)
    with open(file_lut, 'rb') as f:
        res_lut = pickle.load(f)
        
    x_nodes = res_exact['X Coords']
    times = res_exact['Time']
    last_idx = -1
    t_final = times[last_idx]
    
    # In Results.pik, arrays have shape (nNodes, nTimes)
    p_exact = res_exact['Primitive']['Pressure'][:, last_idx]
    p_lut = res_lut['Primitive']['Pressure'][:, last_idx]
    
    rho_exact = res_exact['Primitive']['Density'][:, last_idx]
    rho_lut = res_lut['Primitive']['Density'][:, last_idx]
    
    u_exact = res_exact['Primitive']['Velocity'][:, last_idx]
    u_lut = res_lut['Primitive']['Velocity'][:, last_idx]
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    
    # 1. Pressure
    axs[0].plot(x_nodes, p_exact * 1e-5, 'k-', linewidth=2, label='Exact Helmholtz EoS')
    axs[0].plot(x_nodes[::25], p_lut[::25] * 1e-5, 'r--', markersize=4, label='LuT Splines (250x250)')
    axs[0].set_title(f'Real CO2 Shock Tube (1000 pts) at t = {t_final*1e3:.2f} ms')
    axs[0].set_ylabel('Pressure [bar]')
    axs[0].grid(True, linestyle='--', alpha=0.6)
    axs[0].legend(loc='best')
    
    # 2. Density
    axs[1].plot(x_nodes, rho_exact, 'k-', linewidth=2, label='Exact Helmholtz EoS')
    axs[1].plot(x_nodes[::25], rho_lut[::25], 'r--', markersize=4, label='LuT Splines (250x250)')
    axs[1].set_ylabel('Density [kg/m³]')
    axs[1].grid(True, linestyle='--', alpha=0.6)
    axs[1].legend(loc='best')
    
    # 3. Velocity
    axs[2].plot(x_nodes, u_exact, 'k-', linewidth=2, label='Exact Helmholtz EoS')
    axs[2].plot(x_nodes[::25], u_lut[::25], 'r--', markersize=4, label='LuT Splines (250x250)')
    axs[2].set_xlabel('x [m]')
    axs[2].set_ylabel('Velocity [m/s]')
    axs[2].grid(True, linestyle='--', alpha=0.6)
    axs[2].legend(loc='best')
    
    plt.tight_layout()
    out_img = 'comparison_plot.png'
    plt.savefig(out_img, dpi=200)
    print(f"Comparison plot successfully saved to: {out_img}")

if __name__ == '__main__':
    plot_comparison()
    plt.show()