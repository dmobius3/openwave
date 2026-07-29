"""M8.1.1 figures: the blind coexact-gap measurement, drawn from the solver's own JSON.

Designer-side plotting only. Every number drawn here was produced by the blind solver
(``m8_1_1_coexact_solver.py``); nothing is recomputed and nothing is added.

Figures
    m8_1_1_gap_by_connection.png : the adjoint coexact bottom for every irreducible flat
        SU(2) connection found across the ADE family, showing the single outlier.
    m8_1_1_first_occurrence.png  : first symmetric-power level and first coexact level
        against graph distance, over every irreducible of every group.
    m8_1_1_2i_graph.png          : the measured adjacency graph of the order-120 group,
        laid out by distance, with the two adjoint nodes marked.

Run from this directory:  python3 m8_1_1_plots.py
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "m8_1_1_coexact.json"
PLOTS = HERE.parent / "plots"

GROUPS = json.loads(DATA.read_text())["groups"]


def fig_gap_by_connection():
    """One bar per (group, 2-dimensional determinant-one representation)."""
    labels, gaps, exceptional = [], [], []
    for name, block in GROUPS.items():
        for entry in block["T4"]:
            labels.append(f"{name}\nrho{entry['rho_sigma']}")
            gaps.append(entry["q_squared"])
            exceptional.append(entry["q_squared"] != 4)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    colors = ["#c1121f" if e else "#4a6fa5" for e in exceptional]
    ax.bar(range(len(gaps)), gaps, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel(r"adjoint coexact bottom  $\lambda_{min} R^2$")
    ax.set_title(
        "Blind measurement: adjoint coexact bottom for every irreducible flat SU(2) "
        f"connection ({len(gaps)} connections over 18 groups)",
        fontsize=10,
    )
    for i, g in enumerate(gaps):
        ax.text(i, g + 0.6, str(g), ha="center", fontsize=8)
    ax.set_ylim(0, max(gaps) * 1.15)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = PLOTS / "m8_1_1_gap_by_connection.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_first_occurrence():
    """least symmetric-power level and first coexact level against graph distance."""
    d_all, a_all, e_all = [], [], []
    for block in GROUPS.values():
        dist = block["distance_vector"]
        for row in block["T6"]:
            idx = row["sigma"]
            d_all.append(dist[idx])
            a_all.append(row["least_sym_level"])
            e_all.append(row["e"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    lim = max(d_all + a_all + e_all) + 1

    axes[0].scatter(d_all, a_all, s=45, alpha=0.45, color="#4a6fa5")
    axes[0].plot([0, lim], [0, lim], color="#c1121f", lw=1, ls="--", label="a = d")
    axes[0].set_xlabel("graph distance d(sigma)")
    axes[0].set_ylabel("least a with sigma in V_a restricted")
    axes[0].set_title(f"first symmetric-power level ({len(d_all)} irreducibles)", fontsize=10)
    axes[0].legend(fontsize=8)

    axes[1].scatter(d_all, e_all, s=45, alpha=0.45, color="#2a9d8f")
    axes[1].plot([0, lim], [0, lim], color="#c1121f", lw=1, ls="--", label="e = d")
    axes[1].set_xlabel("graph distance d(sigma)")
    axes[1].set_ylabel("first coexact level e(sigma)")
    axes[1].set_title("first coexact level (note the d = 0 and d = 1 branches)", fontsize=10)
    axes[1].legend(fontsize=8)

    for ax in axes:
        ax.set_xlim(-0.5, lim)
        ax.set_ylim(-0.5, lim)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = PLOTS / "m8_1_1_first_occurrence.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_2i_graph():
    """The measured adjacency graph of the order-120 group, laid out by distance."""
    block = GROUPS["2I"]
    A = block["adjacency_A"]
    dist = block["distance_vector"]
    dims = [row["dim"] for row in block["T6"]]
    adjoints = {entry["constituents"][0]["sigma"] for entry in block["T4"]}

    # x = distance, y spreads nodes that share a distance
    by_d = {}
    pos = {}
    for i, d in enumerate(dist):
        by_d.setdefault(d, []).append(i)
    for d, nodes in by_d.items():
        for j, n in enumerate(nodes):
            pos[n] = (d, 0 if len(nodes) == 1 else (j - (len(nodes) - 1) / 2) * 1.1)

    fig, ax = plt.subplots(figsize=(9, 3.4))
    for i in range(len(dist)):
        for j in range(i + 1, len(dist)):
            if A[i][j]:
                ax.plot(*zip(pos[i], pos[j]), color="#999999", lw=1, zorder=1)
    for i, (x, y) in pos.items():
        marked = i in adjoints
        ax.scatter([x], [y], s=560, zorder=2,
                   color="#c1121f" if marked else "#4a6fa5",
                   edgecolor="black", linewidth=0.6)
        ax.text(x, y, str(dims[i]), ha="center", va="center", color="white",
                fontsize=9, zorder=3, fontweight="bold")
        ax.text(x, y - 0.62, f"d={dist[i]}", ha="center", fontsize=7, color="#444444")
    ax.set_title(
        "Measured adjacency graph of the order-120 group: node label = dimension, "
        "red = the two adjoint nodes",
        fontsize=10,
    )
    ax.set_xlabel("graph distance from the trivial node")
    ax.set_yticks([])
    ax.set_xlim(-0.5, max(dist) + 0.5)
    ax.set_ylim(-1.3, 1.3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    out = PLOTS / "m8_1_1_2i_graph.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


DEFECT = json.loads((HERE.parent / "data" / "m8_1_1_defect.json").read_text())


def _frac(s):
    """'−67/720' style exact string from the solver JSON to a float, for plotting only."""
    num, den = s.split("/")
    return int(num) / int(den)


def fig_affine_relation():
    """Measured defect sum against the affine combination whose coefficients were fitted."""
    rows = DEFECT["U3"]["residuals"]
    a, b, c = (_frac(DEFECT["U3"]["universal_triple"][k]) for k in ("a", "b", "c"))
    d1 = {g: _frac(blk["D_trivial"]) for g, blk in DEFECT["U3"]["per_group"].items()}

    meas, pred = [], []
    for r in rows:
        meas.append(_frac(r["S"]))
        pred.append(a * r["dim"] + b * _frac(r["D"]) + c * r["dim"] * d1[r["group"]])

    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    lo, hi = min(meas + pred) - 1, max(meas + pred) + 1
    ax.plot([lo, hi], [lo, hi], color="#c1121f", lw=1, ls="--", label="exact agreement")
    ax.scatter(pred, meas, s=42, alpha=0.5, color="#4a6fa5")
    ax.set_xlabel("fitted affine combination  a.dim + b.D + c.dim.D(1)")
    ax.set_ylabel("measured defect sum  S(alpha)")
    ax.set_title(
        f"Blind fit: {len(rows)} twists over 17 groups, "
        f"triple (a, b, c) = ({int(a)}, {int(b)}, {int(c)})\n"
        f"nonzero residuals: {DEFECT['U3']['pooled_nonzero_residuals']}",
        fontsize=10,
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = PLOTS / "m8_1_1_affine_relation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_golden_support():
    """Per-class contributions to the character sum for the two adjoint twists."""
    rows = DEFECT["U8"]["rows"]

    def val(s):
        """Exact string possibly carrying sqrt(5), evaluated numerically for the bar height."""
        expr = s.replace("sqrt(5)", str(5 ** 0.5))
        return eval(expr, {"__builtins__": {}})  # noqa: S307, solver-generated literals only

    labels = [f"order {r['element_order']}\nsize {r['size']}" for r in rows]
    p = [val(r["contrib_D_S2P"]) for r in rows]
    pp = [val(r["contrib_D_S2Pprime"]) for r in rows]
    differ = [not r["agree"] for r in rows]

    x = range(len(rows))
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    ax.bar([i - 0.2 for i in x], p, width=0.4, label="adjoint of the defining connection",
           color="#4a6fa5")
    ax.bar([i + 0.2 for i in x], pp, width=0.4, label="adjoint of the conjugate connection",
           color="#e8a33d")
    for i, dif in enumerate(differ):
        if dif:
            ax.axvspan(i - 0.5, i + 0.5, color="#c1121f", alpha=0.10)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("contribution to the character sum D")
    ax.set_title(
        "Class-by-class support of the asymmetry: shaded classes are the only ones where "
        "the two adjoints differ",
        fontsize=10,
    )
    ax.axhline(0, color="black", lw=0.8)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = PLOTS / "m8_1_1_golden_support.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


if __name__ == "__main__":
    for maker in (fig_gap_by_connection, fig_first_occurrence, fig_2i_graph,
                  fig_affine_relation, fig_golden_support):
        print("wrote", maker())
