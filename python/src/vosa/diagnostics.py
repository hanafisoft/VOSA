"""Diagnostics: ceilings, VOSA versus block principal components, heterogeneity, bounds."""

import numpy as np
import pandas as pd

from ._outils import Resultat, _eigh_decroissant, _plafonds, _preparer_S, _saliences
from .resultats import VosaChemin, VosaFit


def plafonds(S):
    """Ceiling (first eigenvalue) of each block and floor 1/(n-1) of the saliences."""
    S = _preparer_S(S, True)
    n = S[0].shape[0]
    return pd.DataFrame({"bloc": S.noms, "plafond": _plafonds(S), "plancher": 1 / (n - 1)})


def diagnostic_blocs(fit, S):
    """Compare each dimension of a fit with the first principal component of each block.

    Returns a DataFrame with, for each block and dimension: the salience, the
    ceiling, the fraction of the ceiling reached, the squared cosine with the
    first principal component of the block, and whether the salience is above
    the level of the dimension.
    """
    S = _preparer_S(S, True)
    U = fit["U"]
    A = np.asarray(fit["saliences"])
    lignes = []
    for k, Sk in enumerate(S):
        val, vec = _eigh_decroissant(Sk)
        v1 = vec[:, 0]
        for r in range(U.shape[1]):
            a = A[k, r]
            lignes.append({"bloc": S.noms[k], "dimension": r + 1, "salience": a,
                           "plafond": val[0], "fraction_plafond": a / val[0],
                           "cos2_cp1": float(np.dot(U[:, r], v1) ** 2),
                           "au_dessus_niveau": bool(a > A[:, r].mean())})
    return pd.DataFrame(lignes)


def heterogeneite(S, U=None):
    """Heterogeneity of the block structures.

    H = (1/K) sum_k ||W_k - W_bar||_F^2, its relative value, the RV index
    1 - lambda_1(RV)/K and, if components U (or a fit) are given, the share of
    H reproduced by the contrast of each dimension (V(u_r) / H) and its
    cumulative sum.
    """
    S = _preparer_S(S, True)
    K = len(S)
    Wbar = sum(S) / K
    H = np.mean([np.sum((Sk - Wbar) ** 2) for Sk in S])
    Hrel = H / np.mean([np.sum(Sk ** 2) for Sk in S])
    Sf = [Sk / np.sqrt(np.sum(Sk ** 2)) for Sk in S]
    RV = np.array([[np.sum(Sf[k] * Sf[l]) for l in range(K)] for k in range(K)])
    out = Resultat(H=H, H_relative=Hrel, indice_RV=1 - np.linalg.eigvalsh(RV)[-1] / K)
    if U is not None:
        if isinstance(U, VosaFit):
            U = U["U"]
        A = _saliences(S, np.asarray(U, float))
        V = ((A - A.mean(axis=0)) ** 2).mean(axis=0)
        out["part_restituee"] = V / H
        out["part_cumulee"] = np.cumsum(V) / H
    return out


def borne_relative(chemin):
    """Relative contrast gain along the path and the bounds of corollary 2.

    For alpha >= 1/2: delta_V (relative gain of contrast with respect to
    INDORT), delta_level (relative loss of squared level), and the bounds
    (1-alpha)/alpha * rho * delta_level <= delta_V <= rho * delta_level.
    """
    tab = chemin["chemin"] if isinstance(chemin, VosaChemin) else chemin
    tab = tab.reset_index(drop=True)
    if np.sum(np.abs(tab["alpha"] - 0.5) < 1e-12) != 1:
        raise ValueError("The path must contain alpha = 0.5.")
    tab = tab[tab["alpha"] >= 0.5 - 1e-12].reset_index(drop=True)
    i0 = int(np.where(np.abs(tab["alpha"] - 0.5) < 1e-12)[0][0])
    rho = tab["niveau2"][i0] / tab["contraste"][i0]
    dl = (tab["niveau2"][i0] - tab["niveau2"]) / tab["niveau2"][i0]
    coef = (1 - tab["alpha"]) / tab["alpha"]
    return pd.DataFrame({"alpha": tab["alpha"],
                         "delta_contraste": (tab["contraste"] - tab["contraste"][i0]) / tab["contraste"][i0],
                         "delta_niveau": dl, "borne_inf": coef * rho * dl,
                         "borne_sup": rho * dl, "rho": rho})
