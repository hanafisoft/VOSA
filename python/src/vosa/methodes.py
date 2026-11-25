"""Methods: CPCA, INDORT, the family g_alpha (SPA algorithm), VOSA, STATIS."""

import warnings

import numpy as np
import pandas as pd

from ._outils import (_depart_aleatoire, _dims, _eigh_decroissant, _fixer_signe,
                      _objet, _polar, _preparer_S, _resume, _rng, _saliences)
from .resultats import VosaChemin, VosaFit


def _verifier_R(R, n):
    if R < 1 or R > n:
        raise ValueError("R must be between 1 and n.")


def _sorties(S, U, R):
    A = _saliences(S, U)
    return pd.DataFrame(A, index=S.noms, columns=_dims(R))


def cpca(S, R=1, normaliser=True):
    """CPCA (alpha = 0): eigenvectors of the mean W_bar of the W_k.

    Parameters
    ----------
    S : list (or dict) of n x n scalar product matrices W_k, for instance from
        :func:`blocs_vers_S`.
    R : int
        Number of dimensions.
    normaliser : bool
        Scale each W_k to unit trace.

    Returns
    -------
    VosaFit
    """
    S = _preparer_S(S, normaliser)
    K, n = len(S), S[0].shape[0]
    _verifier_R(R, n)
    Wbar = sum(S) / K
    val, vec = _eigh_decroissant(Wbar)
    U = _fixer_signe(vec[:, :R])
    A = _sorties(S, U, R)
    return _objet(VosaFit, {"methode": "CPCA", "alpha": 0.0, "U": U, "saliences": A,
                            "resume": _resume(A.values), "valeurs_propres": val[:R],
                            "critere": K * val[:R].sum(), "Wbar": Wbar}, S, S.blocs)


def _indort_un_depart(S, U, maxit, tol):
    K = len(S)
    A = _saliences(S, U)
    Fv = np.sum(A ** 2)
    normeS = max(1.0, sum(np.max(np.abs(np.linalg.eigvalsh(Sk))) ** 2 for Sk in S))
    converge, it = False, 0
    for it in range(1, maxit + 1):
        M = sum(S[k] @ U * A[k] for k in range(K))
        sym = U.T @ M
        sym = (sym + sym.T) / 2
        if np.sqrt(np.sum((M - U @ sym) ** 2)) / normeS < tol:
            converge = True
            break
        U = _polar(M)
        A = _saliences(S, U)
        Fnew = np.sum(A ** 2)
        if Fnew < Fv - 1e-12 * max(1.0, abs(Fv)):
            warnings.warn("The criterion decreased: check that the W_k are positive semi-definite.")
        Fv = Fnew
    return {"U": U, "critere": Fv, "iterations": it, "converge": converge}


def indort(S, R=1, normaliser=True, nstart=10, maxit=1000, tol=1e-8, seed=None):
    """INDORT (orthogonal INDSCAL, simultaneous ComDim), alpha = 1/2.

    Maximizes the sum of the squared saliences with the classical polar
    iteration (monotone for positive semi-definite W_k), from the CPCA
    solution and ``nstart - 1`` random orthonormal starts.

    Returns
    -------
    VosaFit
    """
    S = _preparer_S(S, normaliser)
    K, n = len(S), S[0].shape[0]
    _verifier_R(R, n)
    rng = _rng(seed)
    departs = [cpca(S, R, normaliser=False)["U"]]
    departs += [_depart_aleatoire(rng, n, R) for _ in range(max(nstart - 1, 0))]
    res = [_indort_un_depart(S, U0, maxit, tol) for U0 in departs]
    crit = np.array([r["critere"] for r in res])
    best = res[int(np.argmax(crit))]
    U = best["U"]
    A = _saliences(S, U)
    o = np.argsort(-np.sum(A ** 2, axis=0), kind="stable")
    U = _fixer_signe(U[:, o])
    A = _sorties(S, U, R)
    return _objet(VosaFit, {"methode": "INDORT", "alpha": 0.5, "U": U, "saliences": A,
                            "resume": _resume(A.values), "critere": best["critere"],
                            "iterations": best["iterations"], "converge": best["converge"],
                            "criteres_departs": np.sort(crit)[::-1]}, S, S.blocs)


def _blocs_B(S, alpha):
    """Reduction to J = K + 1 blocks: B_k = sqrt(alpha/K)(W_k - W_bar), B_{K+1} = sqrt(1-alpha) W_bar."""
    K = len(S)
    Wbar = sum(S) / K
    B = [np.sqrt(alpha / K) * (Sk - Wbar) for Sk in S]
    B.append(np.sqrt(1 - alpha) * Wbar)
    return B


def _crit_B(B, U):
    return np.sum(_saliences(B, U) ** 2)


def _spa_un_depart(B, U, maxit, tol):
    J = len(B)
    ev = [np.linalg.eigvalsh(Bj)[[0, -1]] for Bj in B]
    c_dec = sum(max(r[1], 0) * max(-r[0], 0) for r in ev)                # SPA-0
    normeB = max(1.0, sum(max(abs(r[0]), abs(r[1])) ** 2 for r in ev))
    Fv = _crit_B(B, U)
    converge, it = False, 0
    for it in range(1, maxit + 1):
        b = _saliences(B, U)                                             # SPA-2
        M = sum(B[j] @ U * b[j] for j in range(J))
        sym = U.T @ M
        sym = (sym + sym.T) / 2
        if np.sqrt(np.sum((M - U @ sym) ** 2)) / normeB < tol:          # SPA-5
            converge = True
            break
        U0, Uc = _polar(M), _polar(c_dec * U + M)                        # SPA-3
        F0, Fc = _crit_B(B, U0), _crit_B(B, Uc)
        U, Fnew = (U0, F0) if F0 >= Fc else (Uc, Fc)                     # SPA-4
        if Fnew < Fv - 1e-10 * max(1.0, abs(Fv)):
            warnings.warn("The criterion decreased (should not happen).")
        Fv = Fnew
    return {"U": U, "F": Fv, "iterations": it, "converge": converge, "decalage": c_dec}


def spa(S, alpha, R=1, normaliser=True, nstart=10, U0=None, maxit=1000, tol=1e-6, seed=None):
    """Family g_alpha = (1 - alpha) level^2 + alpha contrast, Shifted Polar Ascent.

    Maximizes the mean over the R dimensions of g_alpha on the Stiefel
    manifold, with the SPA algorithm (monotone ascent for indefinite blocks).
    Starts: ``U0`` (an n x R matrix or a list of them), the CPCA solution and
    ``nstart`` random orthonormal matrices; the best solution is kept.

    Parameters
    ----------
    S : list (or dict) of the W_k.
    alpha : float in [0, 1].  0: CPCA, 0.5: INDORT, 1: VOSA.
    R : int, number of dimensions.
    normaliser : bool, scale each W_k to unit trace.
    nstart : int, number of random starts.
    U0 : array or list of arrays, optional starting values.
    maxit, tol : stopping rules (norm of the Riemannian gradient).
    seed : int or numpy Generator, optional.

    Returns
    -------
    VosaFit
    """
    if alpha < 0 or alpha > 1:
        raise ValueError("alpha must be in [0, 1].")
    S = _preparer_S(S, normaliser)
    n = S[0].shape[0]
    _verifier_R(R, n)
    rng = _rng(seed)
    B = _blocs_B(S, alpha)
    departs = [] if U0 is None else (list(U0) if isinstance(U0, (list, tuple)) else [U0])
    departs.append(cpca(S, R, normaliser=False)["U"])
    departs += [_depart_aleatoire(rng, n, R) for _ in range(max(nstart, 0))]
    res = [_spa_un_depart(B, np.asarray(U, float), maxit, tol) for U in departs]
    Fs = np.array([r["F"] for r in res])
    best = res[int(np.argmax(Fs))]
    U = best["U"]
    rs = _resume(_saliences(S, U))
    g = (1 - alpha) * rs["niveau2"].values + alpha * rs["contraste"].values
    o = np.argsort(-g, kind="stable")
    U = _fixer_signe(U[:, o])
    A = _sorties(S, U, R)
    rs = _resume(A.values)
    rs["g_alpha"] = (1 - alpha) * rs["niveau2"] + alpha * rs["contraste"]
    if alpha == 1:
        nom = "VOSA"
    elif alpha == 0.5:
        nom = "INDORT (SPA)"
    elif alpha == 0:
        nom = "CPCA (SPA)"
    else:
        nom = f"SPA (alpha = {alpha:g})"
    return _objet(VosaFit, {"methode": nom, "alpha": float(alpha), "U": U, "saliences": A,
                            "resume": rs, "G_alpha": rs["g_alpha"].mean(),
                            "iterations": best["iterations"], "converge": best["converge"],
                            "decalage": best["decalage"],
                            "criteres_departs": np.sort(Fs / R)[::-1]}, S, S.blocs)


def vosa(S, R=1, **kwargs):
    """VOSA (alpha = 1): maximizes the variance (contrast) of the saliences.

    Same arguments as :func:`spa`, with ``alpha = 1``.
    """
    return spa(S, alpha=1, R=R, **kwargs)


def chemin_spa(S, alphas=None, R=1, normaliser=True, nstart=10, maxit=1000, tol=1e-6, seed=None):
    """Path of the solutions of g_alpha over a grid of alpha.

    Forward sweep with warm starts (plus ``nstart`` random starts at each
    alpha), then a backward sweep; at each alpha the best solution is kept.

    Returns
    -------
    VosaChemin
        ``chemin`` (DataFrame: alpha, mean level^2, contrast, L2, G_alpha),
        ``solutions`` (list of VosaFit) and ``R``.
    """
    if alphas is None:
        alphas = np.round(np.arange(0, 1.0001, 0.05), 10)
    S = _preparer_S(S, normaliser)
    rng = _rng(seed)
    alphas = np.unique(np.asarray(alphas, float))
    m = len(alphas)
    sol = []
    for i in range(m):
        sol.append(spa(S, alphas[i], R, normaliser=False, nstart=nstart,
                       U0=sol[i - 1]["U"] if i > 0 else None, maxit=maxit, tol=tol, seed=rng))
    for i in range(m - 2, -1, -1):
        cand = spa(S, alphas[i], R, normaliser=False, nstart=0, U0=sol[i + 1]["U"],
                   maxit=maxit, tol=tol, seed=rng)
        if cand["G_alpha"] > sol[i]["G_alpha"] + 1e-12:
            sol[i] = cand
    tab = pd.DataFrame({"alpha": [s["alpha"] for s in sol],
                        "niveau2": [s["resume"]["niveau2"].mean() for s in sol],
                        "contraste": [s["resume"]["contraste"].mean() for s in sol],
                        "L2": [s["resume"]["L2"].mean() for s in sol],
                        "G_alpha": [s["G_alpha"] for s in sol]})
    return VosaChemin({"chemin": tab, "solutions": sol, "R": R})


def statis(S, R=1, normalisation="RV"):
    """STATIS, for comparison: eigenvectors of the compromise sum_k gamma_k W_k.

    Parameters
    ----------
    normalisation : {"RV", "trace"}
        ``"RV"``: the W_k are scaled to unit Frobenius norm before computing the
        interstructure (classical STATIS); ``"trace"``: unit trace.

    Returns
    -------
    VosaFit
        With the weights ``poids`` (and ``poids_somme1``), the interstructure,
        the RV matrix, ``part_RV`` and the compromise.
    """
    if normalisation not in ("RV", "trace"):
        raise ValueError("normalisation must be 'RV' or 'trace'.")
    S = _preparer_S(S, normaliser=True)
    K, n = len(S), S[0].shape[0]
    _verifier_R(R, n)
    Sf = [Sk / np.sqrt(np.sum(Sk ** 2)) for Sk in S] if normalisation == "RV" else list(S)
    P = np.array([[np.sum(Sf[k] * Sf[l]) for l in range(K)] for k in range(K)])
    nrm = np.sqrt(np.diag(P))
    RV = P / np.outer(nrm, nrm)
    _, vecP = _eigh_decroissant(P)
    gamma = np.abs(vecP[:, 0])
    Wc = sum(g * Sk for g, Sk in zip(gamma, Sf))
    val, vec = _eigh_decroissant(Wc)
    U = _fixer_signe(vec[:, :R])
    A = _sorties(S, U, R)
    eRV = np.linalg.eigvalsh(RV)[::-1]
    return _objet(VosaFit, {"methode": f"STATIS ({normalisation})", "alpha": np.nan, "U": U,
                            "saliences": A, "resume": _resume(A.values),
                            "poids": pd.Series(gamma, index=S.noms),
                            "poids_somme1": pd.Series(gamma / gamma.sum(), index=S.noms),
                            "interstructure": pd.DataFrame(P, index=S.noms, columns=S.noms),
                            "RV": pd.DataFrame(RV, index=S.noms, columns=S.noms),
                            "part_RV": eRV[0] / K, "compromis": Wc,
                            "valeurs_propres_compromis": val[:R]}, S, S.blocs)
