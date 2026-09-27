"""Shared Matplotlib style for publication figures.

Usage::

    from rydberg_chain import plotstyle
    plotstyle.use()
    fig, axes = plotstyle.subplots(1, 3)
    ...
    plotstyle.panel_labels(axes)
    plotstyle.save(fig, "paper/figures/fig01_example")

Colours are the Okabe-Ito blue, vermillion and bluish green plus a near-black
ink.  Checked on a white background: every pair stays distinguishable under
simulated protanopia and deuteranopia (OKLab distance >= 11), and every colour
has >= 3:1 contrast.  Colour is always paired with a second cue (line style or
marker) so that figures survive greyscale printing.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BLUE = "#0072b2"
VERMILLION = "#d55e00"
GREEN = "#009e73"
INK = "#222222"
MUTED = "#6e6e6e"
GRID = "#e4e4e4"

#: Categorical order.  Assign in this order and keep an entity's colour fixed.
CATEGORICAL = (BLUE, VERMILLION, GREEN)

#: Interaction-range models, keyed as in the study CSVs.
CUTOFF_COLORS = {"nearest1": VERMILLION, "nearest2": BLUE, "all": INK}
CUTOFF_LABELS = {"nearest1": "nearest neighbours",
                 "nearest2": "up to 2nd neighbours",
                 "all": "all pairs"}

#: Chain-length parity: filled circles / solid lines for even N,
#: open squares / dashed lines for odd N.
PARITY_STYLE = {"even": {"marker": "o", "linestyle": "-"},
                "odd": {"marker": "s", "linestyle": "--", "markerfacecolor": "white"}}

#: Line styles used alongside colour for up to three series.
LINESTYLES = ("-", "--", ":")

#: Figure widths in inches (typical two-column journal layout).
SINGLE_COLUMN = 3.4
DOUBLE_COLUMN = 7.0

#: Single-hue sequential ramp (light to dark blue) for maps.
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "rydberg_blue",
    ["#f5f9fe", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

RC = {
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": 8.5,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "legend.handlelength": 2.4,
    "legend.borderaxespad": 0.3,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "axes.prop_cycle": mpl.cycler(color=list(CATEGORICAL)),
    "grid.color": GRID,
    "grid.linewidth": 0.5,
    "grid.linestyle": "-",
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.8,
    "ytick.minor.size": 1.8,
    "lines.linewidth": 1.4,
    "lines.markersize": 4.5,
    "lines.markeredgewidth": 0.9,
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round",
    "figure.dpi": 100,
    "figure.facecolor": "white",
    "figure.constrained_layout.use": True,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "image.cmap": "rydberg_blue",
}

# Metadata without timestamps keeps regenerated files byte-stable.
_METADATA = {"png": {"Software": None},
             "pdf": {"Creator": None, "Producer": None, "CreationDate": None}}


def use():
    """Apply the publication style to the current Matplotlib session."""
    if "rydberg_blue" not in mpl.colormaps:
        mpl.colormaps.register(SEQUENTIAL)
    mpl.rcParams.update(RC)


def subplots(nrows=1, ncols=1, width=DOUBLE_COLUMN, height=None, **kwargs):
    """plt.subplots with a journal width; default height is 0.8 of a panel width."""
    if height is None:
        height = 0.8 * width / ncols * nrows + 0.25
    kwargs.setdefault("squeeze", False)
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), **kwargs)
    return fig, (axes.ravel() if axes.size > 1 else axes[0, 0])


def panel_labels(axes, start="a"):
    """Label panels (a), (b), ... at the top left, above the axes."""
    axes = [axes] if hasattr(axes, "set_title") else list(axes)
    for k, ax in enumerate(axes):
        ax.set_title(f"({chr(ord(start) + k)})", loc="left", fontweight="bold",
                     pad=4)


def save(fig, stem, formats=("png", "pdf"), dpi=300):
    """Write stem.png (300 dpi) and stem.pdf; return the written paths."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = stem.with_suffix("." + fmt)
        fig.savefig(path, dpi=dpi, metadata=_METADATA.get(fmt))
        paths.append(path)
    plt.close(fig)
    return paths
