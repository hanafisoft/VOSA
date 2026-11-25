"""vosa: Variance-Optimized Saliences Analysis for multiblock data (Python version).

For a list of scalar product matrices W_k = X_k X_k' (one per block, scaled to
unit trace), components are described by their saliences a_k(u) = u'W_k u.
The family g_alpha = (1 - alpha) level^2 + alpha contrast contains CPCA
(:func:`cpca`, alpha = 0), INDORT (:func:`indort`, alpha = 1/2) and VOSA
(:func:`vosa`, alpha = 1); :func:`spa` solves any alpha with the Shifted Polar
Ascent algorithm and :func:`chemin_spa` the whole path. :func:`statis` is
provided for comparison, :func:`factorisation_canonique` writes any fit in the
canonical factorization X_k = U D_k V_k' + E_k, and :func:`plafonds`,
:func:`diagnostic_blocs`, :func:`heterogeneite` and :func:`borne_relative`
give the diagnostics of the article. Results are printed in French or English
(:func:`vosa_langue`).

Same functions and outputs as the R package vosa (same repository).
"""

from .diagnostics import borne_relative, diagnostic_blocs, heterogeneite, plafonds
from .factorisation import factorisation_canonique
from .graphiques import (plot_borne_relative, plot_D, plot_niveau_contraste,
                         plot_saliences, plot_U, plot_V)
from .langue import traduire, vosa_langue
from .methodes import chemin_spa, cpca, indort, spa, statis, vosa
from .preparation import bloc_W, blocs_vers_S
from .resultats import VosaChemin, VosaFactorisation, VosaFit
from .simulation import simuler_blocs

__version__ = "1.0.1"

__all__ = [
    "blocs_vers_S", "bloc_W", "cpca", "indort", "spa", "vosa", "chemin_spa", "statis",
    "factorisation_canonique", "plafonds", "diagnostic_blocs", "heterogeneite",
    "borne_relative", "simuler_blocs", "plot_saliences", "plot_U", "plot_D", "plot_V",
    "plot_niveau_contraste", "plot_borne_relative", "vosa_langue", "traduire",
    "VosaFit", "VosaChemin", "VosaFactorisation",
]
