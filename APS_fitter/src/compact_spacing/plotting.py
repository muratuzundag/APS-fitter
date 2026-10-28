"""Publication-style diagnostic plots."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_result(k, result, output, color="tab:blue", show=False):
    k = np.asarray(k, dtype=float)
    order = np.argsort(k)
    k = k[order]
    p = result.corrected_period_s[order]
    e = result.corrected_error_s[order]
    residual = result.residual_s[order]

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.0),
                             gridspec_kw={"height_ratios": [2.0, 1.0, 1.25]})
    ax0, ax1, ax2 = axes
    ax0.errorbar(k, p, yerr=e, fmt="o", color=color, capsize=3,
                 label=result.sequence)
    grid = np.linspace(k.min(), k.max(), 300)
    ax0.plot(grid, result.intercept_s + result.spacing_gls_s * grid,
             "--", color="0.15", label="GLS fit")
    ax0.set_ylabel("Corrected period (s)")
    ax0.legend(frameon=False)

    ax1.axhline(0, color="0.15", lw=1)
    ax1.errorbar(k, residual, yerr=e, fmt="o", color=color, capsize=3)
    ax1.set_xlabel("Relative radial order, k")
    ax1.set_ylabel("Residual (s)")

    ax2.hist(result.samples_s, bins="fd", density=True, color=color, alpha=0.60)
    ax2.axvline(result.spacing_s, color=color, lw=2,
                label=(rf"$\Delta\Pi={result.spacing_s:.3f}^{{+{result.error_plus_s:.3f}}}"
                       rf"_{{-{result.error_minus_s:.3f}}}\,$s"))
    ax2.axvspan(result.p16_s, result.p84_s, color="0.2", alpha=0.12,
                label="16th–84th percentiles")
    ax2.set_xlabel(r"Mean period spacing, $\Delta\Pi$ (s)")
    ax2.set_ylabel("Density")
    ax2.legend(frameon=False)

    for ax in axes:
        ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
