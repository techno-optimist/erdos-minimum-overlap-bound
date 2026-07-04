#!/usr/bin/env python3
"""Regenerate Figure 1 of the note: the Hyra step construction F and its
overlap function M, plotted from the lattice values (float display only; every
certified quantity in the note is exact). Reads certs/erdos_hyra_current.json
and writes fig_hyra_construction.pdf next to this script."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "certs", "erdos_hyra_current.json")
OUT = os.path.join(HERE, "fig_hyra_construction.pdf")

v = np.array(json.load(open(DATA))["values"], dtype=float)
n = len(v)
assert n == 2400
h = 2.0 / n
f = v * ((n / 2) / v.sum())          # board rescale (float mirror of the exact one)
g = 1.0 - f

# lattice correlations C[m] = sum_i f_i g_{i+m} for m = -(n-1)..(n-1)
C = np.correlate(g, f, mode="full")  # C[k] = sum_i f_i g_{i + (k-(n-1))}
lags = np.arange(-(n - 1), n)
M_lattice = (2.0 / n) * C
k_max = int(np.argmax(M_lattice))
m_star, M_star = lags[k_max], M_lattice[k_max]
assert m_star == -92, f"lag convention broken: argmax lag {m_star}"
assert abs(M_star - 0.3808669097979876) < 1e-15, M_star

# extend to m = +-n where M vanishes (supports disjoint)
x_M = np.concatenate(([-2.0], lags * h, [2.0]))
y_M = np.concatenate(([0.0], M_lattice, [0.0]))

plt.rcParams.update({
    "font.size": 9,
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.1, 4.3))

x_F = -1.0 + np.arange(n + 1) * h
ax1.stairs(f, x_F, color="#1f4e79", linewidth=0.5, fill=True, alpha=0.85)
ax1.set_xlim(-1, 1)
ax1.set_ylim(0, 1.02)
ax1.set_xlabel(r"$x$", labelpad=1)
ax1.set_ylabel(r"$F(x)$")

ax2.plot(x_M, y_M, color="black", linewidth=0.9)
ax2.axhline(0.25, color="gray", linewidth=0.7, linestyle=":", zorder=0)
ax2.text(-1.97, 0.262, "mean value $1/4$", color="gray", fontsize=8)
ax2.plot([m_star * h], [M_star], "o", color="#b22222", markersize=4, zorder=5)
ax2.annotate(
    r"$\sup_x M(x) = M(-184/2400) = Q$",
    xy=(m_star * h, M_star), xytext=(0.28, 0.345),
    fontsize=8.5, color="#b22222",
    arrowprops=dict(arrowstyle="-", color="#b22222", linewidth=0.6),
)
ax2.set_xlim(-2, 2)
ax2.set_ylim(0, 0.42)
ax2.set_xlabel(r"$x$", labelpad=1)
ax2.set_ylabel(r"$M(x)$")

fig.tight_layout(h_pad=1.4)
fig.savefig(OUT)
print(f"wrote {OUT}; max at m={m_star}, x*={m_star*h:.6f}, M={M_star:.16f}")
