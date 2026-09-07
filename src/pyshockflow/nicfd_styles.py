import matplotlib as mpl
import matplotlib.pyplot as plt

def set_nicfd_style():
    latex_preamble = r"""
    \usepackage[T1]{fontenc}
    \usepackage[utf8]{inputenc}
    \usepackage{times}
    \usepackage{mathptmx}
    \usepackage{amsmath}
    \usepackage{amssymb}
    """

    mpl.rcParams.update({
        # Typography matching NICFD 2026 template
        "text.usetex": True,
        "text.latex.preamble": latex_preamble,
        "font.family": "serif",
        "font.serif": ["Times", "Times New Roman"],
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        
        # Geometry & rendering
        "savefig.pad_inches": 0.03,
        "figure.dpi": 300,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "axes.linewidth": 0.6,
        "grid.linewidth": 0.4,
        "patch.linewidth": 0.6,
        "lines.markeredgewidth": 0.6,
    })


def get_fig_dim(fraction=1.0, aspect_ratio=1.618, subplots=(1, 1), two_column=False):
    """
    Computes dimensions locked to the NICFD 2026 text width.
    - Full text width: 170 mm (~6.693 in)
    - Column width (if in multicols with ~0.25 in gutter): ~3.22 in
    """
    # 170 mm total text width
    full_text_width_in = 170.0 / 25.4  # 6.6929 in
    col_sep_in = 0.25                  # Standard LaTeX multicol separation
    column_width_in = (full_text_width_in - col_sep_in) / 2.0  # ~3.221 in

    base_width = column_width_in if two_column else full_text_width_in
    fig_width_in = base_width * fraction

    # Calculate width per subplot column, apply aspect ratio for individual subplot height
    ax_width = fig_width_in / subplots[1]
    ax_height = ax_width / aspect_ratio
    fig_height_in = ax_height * subplots[0]

    return (fig_width_in, fig_height_in)


def create_figure(fraction=1.0, aspect_ratio=1.618, subplots=(1, 1), two_column=False):
    """
    Creates properly sized figures matching NICFD 2026 page layout.
    Set `two_column=True` when fitting inside a \begin{multicols}{2} column.
    """
    figsize = get_fig_dim(fraction=fraction, aspect_ratio=aspect_ratio, subplots=subplots, two_column=two_column)
    fig, axes = plt.subplots(subplots[0], subplots[1], figsize=figsize, layout='constrained')
    return fig, axes