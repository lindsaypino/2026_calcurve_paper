"""Shared figure style for the calibration-curves paper — the ONE place the
Talus look is defined. Recycles Lindsay's original vizta 'talk' theme + talusbio
palette so every panel matches.

    from src.style import set_style, LIGHT, DARK, PAL, PRIMARY
    set_style()
"""
import seaborn as sns
import vizta

# vizta.mpl.set_theme returns (PRIMARY navy #0C015B, ACCENT light-blue #64C0CA)
PRIMARY = "#0C015B"
ACCENT = "#64C0CA"          # light Talus blue
PAL = sns.color_palette("talusbio")  # [0]#0086bb blue [1]#ee8156 orange [2]#66c2a5 teal [3]#eb98b9 pink [4]#fed766 yellow

# two-series comparison colors (light = first condition, dark = second)
LIGHT = ACCENT              # #64C0CA
DARK = PAL[0]               # #0086bb

# figures-of-merit accent colors (LOD / LOQ / ULOQ markers on example curves)
FOM_COLORS = {"LOD": PAL[1], "LOQ": PAL[2], "ULOQ": PAL[3]}  # orange / teal / pink


def set_style(context="talk"):
    """Apply the vizta theme; returns (PRIMARY, ACCENT)."""
    primary, accent = vizta.mpl.set_theme(context=context)
    return primary, accent


def kde_pair(ax, va, vb, label_a, label_b, fill=True):
    """Draw a two-series filled KDE (light vs dark) — the paper's standard comparison."""
    if len(va) > 2:
        sns.kdeplot(va, ax=ax, color=LIGHT, fill=fill, alpha=0.3, linewidth=2, label=label_a)
    if len(vb) > 2:
        sns.kdeplot(vb, ax=ax, color=DARK, fill=fill, alpha=0.3, linewidth=2, label=label_b)
    ax.grid(True, alpha=0.3)
