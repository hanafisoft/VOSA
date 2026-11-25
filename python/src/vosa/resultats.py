"""Result objects: printing in French or English, summaries and plot shortcuts."""

import numpy as np
import pandas as pd

from ._outils import Resultat
from .langue import _tr, traduire


def _fmt(x, digits):
    if isinstance(x, (pd.DataFrame, pd.Series)):
        return x.round(digits).to_string()
    return str(np.round(x, digits))


class VosaFit(Resultat):
    """Fitted model (CPCA, INDORT, SPA, VOSA, STATIS).

    Main items: ``methode``, ``alpha``, ``U`` (n x R components), ``saliences``
    (K x R DataFrame), ``resume`` (level and contrast by dimension), ``plafonds``
    (first eigenvalue of each block), ``n``, ``K``, ``blocs`` (normalized data
    blocks, if available).
    """

    def texte(self, digits=3):
        R = self["U"].shape[1]
        lignes = [f"{self['methode']} : n = {self['n']}, K = {self['K']}, R = {R}", "",
                  _tr("saliences"), _fmt(self["saliences"], digits), "",
                  _tr("resume"), traduire(self["resume"].round(digits)).to_string(index=False)]
        return "\n".join(lignes)

    def __repr__(self):
        return self.texte()

    __str__ = __repr__

    def summary(self, afficher=True):
        """Saliences with the ceilings of the blocks, and the floor 1/(n-1)."""
        tab = self["saliences"].copy()
        tab.insert(0, _tr("plafond"), self["plafonds"].values)
        plancher = 1 / (self["n"] - 1)
        if afficher:
            print(tab.round(3).to_string())
            print(_tr("plancher"), round(plancher, 3))
        return Resultat(methode=self["methode"], saliences=tab, resume=self["resume"],
                        plancher=plancher)

    def plot(self, dim=1, ax=None, main=None):
        """Bar chart of the saliences of one dimension, with the ceilings of the blocks."""
        from .graphiques import plot_saliences
        return plot_saliences(self, dim=dim, ax=ax, main=main)


class VosaChemin(Resultat):
    """Path of the solutions of g_alpha (items ``chemin``, ``solutions``, ``R``)."""

    def __repr__(self):
        return f"{_tr('chemin')} {self['R']}\n" + \
            traduire(self["chemin"].round(4)).to_string(index=False)

    __str__ = __repr__


class VosaFactorisation(Resultat):
    """Canonical factorization X_k = U D_k V_k' + E_k.

    Items: ``U``, ``D`` (K x R DataFrame of the d_kr), ``V`` (dict of p_k x R
    loadings), ``E`` (dict of residuals), ``part_restituee``,
    ``part_residuelle``, ``erreur_reconstruction``, ``methode``.
    """

    def __repr__(self, digits=3):
        lignes = [f"{_tr('fc_titre')} ({self['methode']})", "",
                  _tr("fc_D"), _fmt(self["D"], digits), "",
                  _tr("fc_part"), _fmt(self["part_restituee"], digits), "",
                  _tr("fc_resid"), _fmt(self["part_residuelle"], digits), "",
                  f"{_tr('fc_err')} {self['erreur_reconstruction']:.3g}"]
        return "\n".join(lignes)

    __str__ = __repr__

    def plot(self, type="D", **kwargs):
        """``type="D"``, ``"U"`` or ``"V"``: calls plot_D, plot_U or plot_V."""
        from .graphiques import plot_D, plot_U, plot_V
        fonctions = {"D": plot_D, "U": plot_U, "V": plot_V}
        if type not in fonctions:
            raise ValueError("type must be 'D', 'U' or 'V'.")
        return fonctions[type](self, **kwargs)
