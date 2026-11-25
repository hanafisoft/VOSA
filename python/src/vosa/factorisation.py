"""Canonical factorization X_k = U D_k V_k' + E_k of any fitted model."""

import numpy as np
import pandas as pd

from .resultats import VosaFactorisation, VosaFit


def factorisation_canonique(fit, blocs=None):
    """Write a fitted model in the canonical factorization X_k = U D_k V_k' + E_k.

    With orthonormal components U: d_kr = ||X_k' u_r|| (its square is the
    salience when the blocks have unit inertia), V_k = X_k' U D_k^{-1} (unit
    columns) and E_k = (I - UU') X_k (Bougeard et al., 2024).

    Parameters
    ----------
    fit : VosaFit
    blocs : list of the normalized data blocks; by default those kept in the
        fit (S built with :func:`blocs_vers_S` or :func:`bloc_W`).

    Returns
    -------
    VosaFactorisation
        ``U``, ``D``, ``V``, ``E``, ``part_restituee``, ``part_residuelle``,
        ``erreur_reconstruction``, ``methode``.
    """
    if not isinstance(fit, VosaFit):
        raise TypeError("fit must be a VosaFit object.")
    if blocs is None:
        blocs = fit["blocs"]
    if blocs is None:
        raise ValueError("No data blocks: build S with blocs_vers_S() or bloc_W(), or supply 'blocs'.")
    U = fit["U"]
    noms = list(fit["saliences"].index)
    dims = list(fit["saliences"].columns)
    if len(blocs) != len(noms):
        raise ValueError("The number of blocks does not match the fitted model.")
    R = U.shape[1]
    D = np.zeros((len(blocs), R))
    part = np.zeros((len(blocs), R))
    resid = np.zeros(len(blocs))
    V, E = {}, {}
    err = 0.0
    for k, Xk in enumerate(blocs):
        Xk = np.asarray(Xk, float)
        if Xk.shape[0] != U.shape[0]:
            raise ValueError("Blocks and components must have the same number of rows.")
        P = Xk.T @ U                                    # p_k x R, columns X_k' u_r
        d = np.sqrt(np.sum(P ** 2, axis=0))
        Vk = np.where(d > 0, P / np.where(d > 0, d, 1), 0.0)
        Ek = Xk - U @ P.T                               # = (I - UU') X_k
        tot = np.sum(Xk ** 2)
        D[k] = d
        part[k] = d ** 2 / tot
        resid[k] = np.sum(Ek ** 2) / tot
        V[noms[k]] = Vk
        E[noms[k]] = Ek
        err = max(err, np.max(np.abs(U @ np.diag(d) @ Vk.T + Ek - Xk)))
    return VosaFactorisation({
        "U": U,
        "D": pd.DataFrame(D, index=noms, columns=dims),
        "V": V, "E": E,
        "part_restituee": pd.DataFrame(part, index=noms, columns=dims),
        "part_residuelle": pd.Series(resid, index=noms),
        "erreur_reconstruction": err,
        "methode": fit["methode"]})
