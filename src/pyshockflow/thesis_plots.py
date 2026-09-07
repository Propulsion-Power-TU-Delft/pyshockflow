import matplotlib.pyplot as plt
import matplotlib as mpl

"""For a standalone single plot, 
setting your width fraction to 0.6 with an aspect ratio of 1.3 provides the most objective balance between data legibility and page economy.
"""


def set_thesis_style():
    latex_preamble = r"""
    \usepackage[T1]{fontenc}
    \usepackage[utf8]{inputenc}
    \usepackage{lmodern}
    \usepackage{fourier}
    \usepackage{amsmath}
    """

    mpl.rcParams.update({
        # Existing typography settings
        "text.usetex": True,
        "text.latex.preamble": latex_preamble,
        "font.family": "serif",
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "savefig.pad_inches": 0.05,
        "figure.dpi": 300,
        
        # New geometry and thickness settings
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "axes.linewidth": 0.6,
        "grid.linewidth": 0.4,
        "patch.linewidth": 0.6,
        "lines.markeredgewidth": 0.6,
    })
    

def get_fig_dim(fraction=1.0, aspect_ratio=1.618, subplots=(1, 1), is_print=False):
    """
    Computes canvas dimensions locked to the TU Delft text width.
    """
    # 125.25 mm for print, 127.5 mm for digital
    text_width_in = 4.931 if is_print else 5.0196 
    
    fig_width_in = text_width_in * fraction
    
    # Calculate width per column, then apply aspect ratio for height
    ax_width = fig_width_in / subplots[1]
    ax_height = ax_width / aspect_ratio
    fig_height_in = ax_height * subplots[0]
    
    return (fig_width_in, fig_height_in)

def create_figure(fraction=1.0, aspect_ratio=1.618, subplots=(1, 1), is_print=False):
    """
    Shortcut to create a properly sized figure and axes.
    """
    figsize = get_fig_dim(fraction=fraction, aspect_ratio=aspect_ratio, subplots=subplots, is_print=is_print)
    fig, axes = plt.subplots(subplots[0], subplots[1], figsize=figsize, layout='constrained')
    return fig, axes

# example of use
if __name__ == "__main__":
    
    import numpy as np
    # import thesis_plots as tp

    # 1. Apply the TU Delft pdfLaTeX typography rules globally
    set_thesis_style()

    # 2. Create a full-width figure with a standard 1.618 aspect ratio
    fig, ax = create_figure(fraction=1.00, aspect_ratio=1.618, subplots=(1, 1))

    # 3. Generate some dummy data (mock density profile across a shock)
    x = np.linspace(0, 10, 200)
    y = 1.2 + 0.8 * np.tanh(x - 5)

    # 4. Plot the data
    ax.plot(x, y, color='black', linewidth=1.2, label="Numerical solution")
    ax.axvline(x=5, color='red', linestyle='--', linewidth=1.0, label="Exact shock location")

    # 5. Add LaTeX-formatted labels to test the Fourier math font
    ax.set_xlabel(r"Axial coordinate $x$ (m)")
    ax.set_ylabel(r"Density $\rho~\mathrm{[kg/m^3]}$")
    ax.set_title(r"Validation of typography")

    # 6. Add legend and grid
    ax.legend(loc="upper left")
    ax.grid(True, linestyle=':', alpha=0.6)

    # 7. Export the figure to PDF
    fig.savefig("typography_test_fullpage.pdf")
    print("Plot successfully saved as typography_test_fullpage.pdf")




    # 2. Create a half-width figure with a standard 1.618 aspect ratio
    fig, ax = create_figure(fraction=0.49, aspect_ratio=1.3, subplots=(1, 1), is_print=False)

    # 3. Generate some dummy data (mock density profile across a shock)
    x = np.linspace(0, 10, 200)
    y = 1.2 + 0.8 * np.tanh(x - 5)

    # 4. Plot the data
    ax.plot(x, y, color='black', linewidth=1.2, label="Numerical solution")
    ax.axvline(x=5, color='red', linestyle='--', linewidth=1.0, label="Exact shock location")

    # 5. Add LaTeX-formatted labels to test the Fourier math font
    ax.set_xlabel(r"Axial coordinate $x$ (m)")
    ax.set_ylabel(r"Density $\rho~\mathrm{[kg/m^3]}$")
    ax.set_title(r"Validation of typography")

    # 6. Add legend and grid
    ax.legend(loc="upper left")
    ax.grid(True, linestyle=':', alpha=0.6)

    # 7. Export the figure to PDF
    fig.savefig("typography_test_halfpage.pdf")
    print("Plot successfully saved as typography_test_halfpage.pdf")