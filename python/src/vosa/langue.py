"""Language of the results (French or English).

The language governs printed output, plot labels and the column names of the
tables shown. Attribute and key names of the result objects stay in French
(``fit.saliences``, ``fit.resume``, ``fc.part_restituee``...), so that a script
works whatever the language; :func:`traduire` returns a copy with translated
names, for export.
"""

import os

import pandas as pd

_LANGUES = ("fr", "en")
_etat = {"langue": os.environ.get("VOSA_LANGUE", "fr") if os.environ.get("VOSA_LANGUE", "fr") in _LANGUES else "fr"}


def vosa_langue(langue=None):
    """Read or set the language of the results.

    Parameters
    ----------
    langue : {"fr", "en"} or None
        ``None`` returns the current language. The default language is French;
        it can also be set with the environment variable ``VOSA_LANGUE``.

    Returns
    -------
    str
        The current (or newly set) language.
    """
    if langue is None:
        return _etat["langue"]
    if langue not in _LANGUES:
        raise ValueError("langue must be 'fr' or 'en'.")
    _etat["langue"] = langue
    return langue


DICO_NOMS = {
    "methode": "method", "alpha": "alpha", "U": "U", "saliences": "saliences",
    "resume": "summary", "dimension": "dimension", "niveau": "level",
    "niveau2": "level2", "contraste": "contrast", "L2": "L2", "rho": "rho",
    "g_alpha": "g_alpha", "G_alpha": "G_alpha", "valeurs_propres": "eigenvalues",
    "critere": "criterion", "Wbar": "Wbar", "iterations": "iterations",
    "converge": "converged", "decalage": "shift", "criteres_departs": "start_criteria",
    "n": "n", "K": "K", "plafonds": "ceilings", "blocs": "blocks",
    "poids": "weights", "poids_somme1": "weights_sum1",
    "interstructure": "interstructure", "RV": "RV", "part_RV": "RV_share",
    "compromis": "compromise", "valeurs_propres_compromis": "compromise_eigenvalues",
    "chemin": "path", "solutions": "solutions", "R": "R",
    "D": "D", "V": "V", "E": "E", "part_restituee": "reproduced_share",
    "part_residuelle": "residual_share", "erreur_reconstruction": "reconstruction_error",
    "bloc": "block", "plafond": "ceiling", "plancher": "floor", "salience": "salience",
    "fraction_plafond": "ceiling_fraction", "cos2_cp1": "cos2_pc1",
    "au_dessus_niveau": "above_level", "H": "H", "H_relative": "H_relative",
    "indice_RV": "RV_index", "part_cumulee": "cumulative_share",
    "delta_contraste": "contrast_gain", "delta_niveau": "level_loss",
    "borne_inf": "lower_bound", "borne_sup": "upper_bound",
    "X": "X", "t_d": "t_d", "T": "T", "jeu": "dataset", "variable": "variable",
    "chargement": "loading",
}

DICO_MESSAGES = {
    "saliences": {"fr": "Saliences :", "en": "Saliences:"},
    "resume": {"fr": "Niveau et contraste par dimension :", "en": "Level and contrast by dimension:"},
    "plancher": {"fr": "Plancher 1/(n-1) =", "en": "Floor 1/(n-1) ="},
    "chemin": {"fr": "Chemin des solutions de g_alpha, R =", "en": "Path of solutions of g_alpha, R ="},
    "fc_titre": {"fr": "Factorisation canonique X_k = U D_k V_k' + E_k",
                 "en": "Canonical factorization X_k = U D_k V_k' + E_k"},
    "fc_D": {"fr": "Termes diagonaux d_kr des D_k :", "en": "Diagonal terms d_kr of D_k:"},
    "fc_part": {"fr": "Part d'inertie restituée (d_kr^2 / ||X_k||^2) :",
                "en": "Share of inertia reproduced (d_kr^2 / ||X_k||^2):"},
    "fc_resid": {"fr": "Part résiduelle ||E_k||^2 / ||X_k||^2 :",
                 "en": "Residual share ||E_k||^2 / ||X_k||^2:"},
    "fc_err": {"fr": "Erreur de reconstruction :", "en": "Reconstruction error:"},
    "salience": {"fr": "salience", "en": "salience"},
    "plafond": {"fr": "plafond", "en": "ceiling"},
    "niveau": {"fr": "niveau", "en": "level"},
    "niveau2": {"fr": "niveau au carré", "en": "squared level"},
    "contraste": {"fr": "contraste", "en": "contrast"},
    "dimension": {"fr": "dimension", "en": "dimension"},
    "individu": {"fr": "individu", "en": "individual"},
    "composante": {"fr": "composante", "en": "component"},
    "variable": {"fr": "variable", "en": "variable"},
    "chargement": {"fr": "chargement, dim.", "en": "loading, dim."},
    "d2": {"fr": "d_kr² (salience)", "en": "d_kr² (salience)"},
    "d": {"fr": "d_kr", "en": "d_kr"},
    "score": {"fr": "score x racine(n)", "en": "score x sqrt(n)"},
}


def _tr(cle, langue=None):
    langue = langue or vosa_langue()
    m = DICO_MESSAGES.get(cle)
    return cle if m is None else m[langue]


def _tr_noms(noms, langue=None):
    langue = langue or vosa_langue()
    if langue == "fr":
        return list(noms)
    return [DICO_NOMS.get(str(x), x) for x in noms]


def traduire(x, langue=None):
    """Copy of a result object, dict or DataFrame with names in the chosen language.

    Values are unchanged; only dictionary keys and DataFrame column names are
    translated. With ``langue="fr"`` the names are returned unchanged.
    """
    langue = langue or vosa_langue()
    if langue not in _LANGUES:
        raise ValueError("langue must be 'fr' or 'en'.")
    if isinstance(x, pd.DataFrame):
        out = x.copy()
        out.columns = _tr_noms(out.columns, langue)
        return out
    if isinstance(x, dict):
        cles = _tr_noms(list(x.keys()), langue)
        return {c: (traduire(v, langue) if isinstance(v, (dict, pd.DataFrame)) else v)
                for c, v in zip(cles, x.values())}
    return x
