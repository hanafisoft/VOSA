"""Plots (matplotlib): saliences, components U, diagonal terms D_k, loadings V_k,
level/contrast plane and relative bounds. Labels follow :func:`vosa_langue`.

Every function draws on ``ax`` (a matplotlib Axes) or, if ``ax`` is None, on a
new figure, and returns the plotted values.
"""

import numpy as np

from ._outils import _fusion_etiquettes
from .langue import _tr
from .resultats import VosaFactorisation, VosaFit


def _axe(ax):
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots()
    return ax


def _couleurs(m):
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab10")
    return [cmap(i % 10) for i in range(m)]


def plot_saliences(fit, dim=1, ax=None, main=None):
    """Bar chart of the saliences of one dimension, with the ceilings and the level."""
    ax = _axe(ax)
    a = fit["saliences"].iloc[:, dim - 1]
    pos = np.arange(len(a))
    ax.bar(pos, a.values, color="0.75", edgecolor="black")
    ax.scatter(pos, fit["plafonds"].values, marker="v", color="0.25", zorder=3, label=_tr("plafond"))
    ax.axhline(a.mean(), linestyle="--", color="black", label=_tr("niveau"))
    ax.set_xticks(pos)
    ax.set_xticklabels(a.index, rotation=90)
    ax.set_ylim(0, max(a.max(), fit["plafonds"].max()) * 1.3)
    ax.set_ylabel(_tr("salience"))
    ax.set_title(main if main is not None else f"{fit['methode']} - {_tr('dimension')} {dim}")
    ax.legend(loc="upper center", ncol=2, frameon=False)
    return a


def plot_U(x, dims=(1, 2), groupes=None, etiquettes=None, ax=None, main=None):
    """Individuals on the components U (index plot for one dimension, scatter for two)."""
    ax = _axe(ax)
    U = x["U"]
    dims = [dims] if np.isscalar(dims) else list(dims)
    if max(dims) > U.shape[1]:
        raise ValueError("dims exceeds the number of components.")
    n = U.shape[0]
    noms = [str(i + 1) for i in range(n)]
    if etiquettes is None:
        etiquettes = n <= 50
    if groupes is None:
        col = ["0.15"] * n
        niveaux = []
    else:
        groupes = np.asarray(groupes)
        if len(groupes) != n:
            raise ValueError("groupes must have one value per individual.")
        niveaux = list(dict.fromkeys(groupes.tolist()))
        pal = _couleurs(len(niveaux))
        col = [pal[niveaux.index(g)] for g in groupes]
    titre = main if main is not None else x.get("methode", "")
    if len(dims) == 1:
        s = U[:, dims[0] - 1]
        xs, ys = np.arange(1, n + 1), s
        ax.scatter(xs, ys, c=col, s=18)
        ax.axhline(0, linestyle=":", color="black")
        ax.set_xlabel(_tr("individu"))
        ax.set_ylabel(f"{_tr('composante')} {dims[0]}")
        out = s
    else:
        xs, ys = U[:, dims[0] - 1], U[:, dims[1] - 1]
        ax.scatter(xs, ys, c=col, s=18)
        ax.axhline(0, linestyle=":", color="black")
        ax.axvline(0, linestyle=":", color="black")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel(f"{_tr('composante')} {dims[0]}")
        ax.set_ylabel(f"{_tr('composante')} {dims[1]}")
        out = U[:, [dims[0] - 1, dims[1] - 1]]
    if etiquettes:
        for xi, yi, t, c in zip(xs, ys, noms, col):
            ax.annotate(t, (xi, yi), textcoords="offset points", xytext=(0, 4),
                        ha="center", fontsize=6, color=c)
    for g, c in zip(niveaux, _couleurs(len(niveaux))):
        ax.scatter([], [], color=c, label=str(g))
    if niveaux:
        ax.legend(loc="upper right", frameon=False)
    ax.set_title(titre)
    return out


def plot_D(x, carre=True, ax=None, main=None):
    """Bar chart of the d_kr (or their squares, the saliences): blocks by dimensions."""
    ax = _axe(ax)
    if isinstance(x, VosaFit):
        D = np.sqrt(np.clip(x["saliences"], 0, None))
    else:
        D = x["D"]
    M = D ** 2 if carre else D
    K, R = M.shape
    largeur = 0.8 / R
    pos = np.arange(K)
    for r in range(R):
        ax.bar(pos - 0.4 + largeur * (r + 0.5), M.iloc[:, r].values, width=largeur,
               color=_couleurs(R)[r], edgecolor="black", label=M.columns[r])
    ax.set_xticks(pos)
    ax.set_xticklabels(M.index, rotation=90)
    ax.set_ylim(0, M.values.max() * 1.15)
    ax.set_ylabel(_tr("d2") if carre else _tr("d"))
    ax.set_title(main if main is not None else x.get("methode", ""))
    ax.legend(loc="upper right", frameon=False)
    return M


def plot_V(x, bloc=0, dims=1, axe=None, etiquettes=True, n_etiquettes=10, ax=None, main=None):
    """Loadings V_k of one block: profile for one dimension, scatter for two.

    ``bloc`` is the block name or its index (from 0). When the block has more
    than 40 variables, only the ``n_etiquettes`` variables with the largest
    loadings are labelled; ``axe`` gives an optional abscissa (wavelengths...).
    """
    if not isinstance(x, VosaFactorisation):
        raise TypeError("x must be a VosaFactorisation (see factorisation_canonique).")
    ax = _axe(ax)
    noms_blocs = list(x["V"].keys())
    nom_bloc = bloc if isinstance(bloc, str) else noms_blocs[bloc]
    Vk = x["V"][nom_bloc]
    p = Vk.shape[0]
    noms = [str(j + 1) for j in range(p)]
    dims = [dims] if np.isscalar(dims) else list(dims)
    ax.set_title(main if main is not None else f"{nom_bloc} - {x['methode']}")
    if len(dims) == 1:
        v = Vk[:, dims[0] - 1]
        ab = np.arange(1, p + 1) if axe is None else np.asarray(axe)
        if len(ab) != p:
            raise ValueError("axe must have one value per variable of the block.")
        if p > 40:
            ax.plot(ab, v, color="black", linewidth=1)
        else:
            ax.vlines(ab, 0, v, color="black", linewidth=3)
        ax.axhline(0, linestyle=":", color="black")
        ax.set_xlabel(_tr("variable") if axe is None else "")
        ax.set_ylabel(f"{_tr('chargement')} {dims[0]}")
        if etiquettes:
            k = range(p) if p <= 40 else np.argsort(-np.abs(v))[:min(n_etiquettes, p)]
            for j in k:
                ax.annotate(noms[j], (ab[j], v[j]), textcoords="offset points",
                            xytext=(0, 4 if v[j] >= 0 else -9), ha="center", fontsize=7)
        return v
    v = Vk[:, [dims[0] - 1, dims[1] - 1]]
    lim = np.abs(v).max() * 1.15
    ax.scatter(v[:, 0], v[:, 1], color="0.2", s=14)
    ax.axhline(0, linestyle=":", color="black")
    ax.axvline(0, linestyle=":", color="black")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{_tr('chargement')} {dims[0]}")
    ax.set_ylabel(f"{_tr('chargement')} {dims[1]}")
    if etiquettes:
        k = range(p) if p <= 40 else np.argsort(-np.sum(v ** 2, axis=1))[:min(n_etiquettes, p)]
        for j in k:
            ax.annotate(noms[j], (v[j, 0], v[j, 1]), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=7)
    return v


def plot_niveau_contraste(fits, ax=None, main=None):
    """Mean squared level against mean contrast for several fits (dict name -> fit)."""
    import pandas as pd
    ax = _axe(ax)
    if not isinstance(fits, dict):
        fits = {f["methode"]: f for f in fits}
    tab = pd.DataFrame({"methode": list(fits.keys()),
                        "niveau2": [f["resume"]["niveau2"].mean() for f in fits.values()],
                        "contraste": [f["resume"]["contraste"].mean() for f in fits.values()]})
    ax.scatter(tab["niveau2"], tab["contraste"], color="black")
    for xe, ye, t in _fusion_etiquettes(tab["niveau2"], tab["contraste"], list(tab["methode"])):
        ax.annotate(t, (xe, ye), textcoords="offset points", xytext=(5, 0), ha="left",
                    va="center", fontsize=8)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_xlim(0, tab["niveau2"].max() * 1.45)
    ax.set_ylim(0, tab["contraste"].max() * 1.25)
    ax.set_xlabel(_tr("niveau2"))
    ax.set_ylabel(_tr("contraste"))
    if main is not None:
        ax.set_title(main)
    return tab


def plot_borne_relative(b, ax=None, main=None):
    """Relative contrast gain along the path, inside the band of corollary 2."""
    ax = _axe(ax)
    ax.fill_between(b["alpha"], b["borne_inf"], b["borne_sup"], color="steelblue", alpha=0.25)
    ax.plot(b["alpha"], b["delta_contraste"], color="black", linewidth=2)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\delta_V$")
    if main is not None:
        ax.set_title(main)
    return b
