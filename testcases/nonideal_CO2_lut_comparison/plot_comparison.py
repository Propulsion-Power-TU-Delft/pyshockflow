"""
Plot comparison between Exact Helmholtz EoS and Look-Up Table (LuT) 2D Bicubic Splines.
"""

import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from pyshockflow.thesis_plots import *
set_thesis_style()

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
    
    fig, axs = create_figure(fraction=1.0, aspect_ratio=0.7, subplots=(1, 3), is_print=False)
    
    # 1. Pressure
    axs[0].plot(x_nodes, p_exact * 1e-5, 'k-', label='Exact Helmholtz EoS')
    axs[0].plot(x_nodes[::25], p_lut[::25] * 1e-5, 'r--', label='LuT (250x250)')
    axs[0].set_ylabel(r'$P$ [bar]')
    axs[0].set_xlabel(r'$x$ [m]')
    axs[0].grid(True, linestyle='--', alpha=0.6)
    
    # 2. Density
    axs[1].plot(x_nodes, rho_exact, 'k-')
    axs[1].plot(x_nodes[::25], rho_lut[::25], 'r--')
    axs[1].set_ylabel(r'$\rho$ [kg/m$^3$]')
    axs[1].set_xlabel(r'$x$ [m]')
    axs[1].grid(True, linestyle='--', alpha=0.6)
    
    # 3. Velocity
    axs[2].plot(x_nodes, u_exact, 'k-')
    axs[2].plot(x_nodes[::25], u_lut[::25], 'r--')
    axs[2].set_xlabel(r'$x$ [m]')
    axs[2].set_ylabel(r'$u$ [m/s]')
    axs[2].grid(True, linestyle='--', alpha=0.6)
    
    fig.legend(ncol=2, loc='outside upper center')
    out_img = 'comparison_plot.pdf'
    plt.savefig(out_img)
    print(f"Comparison plot successfully saved to: {out_img}")

if __name__ == '__main__':
    plot_comparison()
    plt.show()