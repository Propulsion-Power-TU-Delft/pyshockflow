"""Read solver residuals from a log file and plot them vs iteration number."""

import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from pyshockflow.thesis_plots import *


def read_residuals(filename):
    """Parse iterations and residuals from the log file.

    Returns the list of iteration numbers and a list of residual tuples,
    one value per equation (e.g. continuity, momentum, energy).
    """
    iterations = []
    residuals = []

    line_re = re.compile(
        r"Iteration\s+(\d+).*Residuals:\s*(.+)$"
    )

    with open(filename, "r") as f:
        for line in f:
            match = line_re.search(line)
            if match:
                it = int(match.group(1))
                values = [float(v) for v in match.group(2).split(",")]
                iterations.append(it)
                residuals.append(values)

    return iterations, residuals


def plot_residuals(iterations, residuals, labels=None):
    """Plot each residual component as a function of the iteration number."""
    n_components = len(residuals[0])
    if labels is None:
        labels = [f"Residual {i + 1}" for i in range(n_components)]

    set_thesis_style()
    fig, ax = create_figure(fraction=0.5, aspect_ratio=1.3, subplots=(1, 1), is_print=False)
    for i in range(n_components):
        component = [r[i] for r in residuals]
        ax.plot(iterations, component-component[0], label=labels[i])

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residuals drop")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig("residuals_nozzle.pdf")
    plt.show()
    


def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "residuals.log"
    iterations, residuals = read_residuals(filename)

    if not iterations:
        print(f"No residuals found in '{filename}'.")
        return

    print(f"Read {len(iterations)} iterations from '{filename}'.")
    labels = ["Continuity", "Momentum", "Energy"]
    plot_residuals(iterations, np.array(residuals), labels)


if __name__ == "__main__":
    main()
