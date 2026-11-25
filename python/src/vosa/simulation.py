"""Simulated multiblock data with common and block-specific factors."""

import numpy as np

from ._outils import Resultat, _rng


def simuler_blocs(n=30, K=4, p=100, pi_c=(0.20, 0.15, 0.10), pi_d=(0.30, 0.30, 0, 0), seed=None):
    """Generator used in the simulation study.

    ``n`` individuals, ``K`` blocks of ``p`` variables; factors present in all
    blocks with shares of inertia ``pi_c``, and a factor t_d present in block
    k with share ``pi_d[k]``. Factors are centred and orthonormal; each block
    has its own random loadings; the remaining inertia is Gaussian noise. Each
    block is centred and scaled to unit inertia.

    Returns
    -------
    Resultat
        ``X`` (dict of the K blocks), ``t_d`` and ``T`` (all factors, the last
        column being t_d).
    """
    pi_c = np.asarray(pi_c, float)
    pi_d = np.asarray(pi_d, float)
    if len(pi_d) != K:
        raise ValueError("pi_d must have length K.")
    if np.any(pi_c.sum() + pi_d >= 1):
        raise ValueError("The shares of inertia must add up to less than 1.")
    rng = _rng(seed)
    q = len(pi_c) + 1
    Fac = rng.standard_normal((n, q))
    Fac = Fac - Fac.mean(axis=0)
    Fac, _ = np.linalg.qr(Fac)
    X = {}
    for k in range(K):
        parts = np.r_[pi_c, pi_d[k]]
        Xk = np.zeros((n, p))
        for j in range(q):
            v = rng.standard_normal(p)
            v = v / np.linalg.norm(v)
            Xk += np.sqrt(parts[j]) * np.outer(Fac[:, j], v)
        E = rng.standard_normal((n, p))
        E = E - E.mean(axis=0)
        E = E / np.sqrt(np.sum(E ** 2))
        Xk += np.sqrt(1 - parts.sum()) * E
        Xk = Xk - Xk.mean(axis=0)
        X[f"Block{k + 1}"] = Xk / np.sqrt(np.sum(Xk ** 2))
    return Resultat(X=X, t_d=Fac[:, q - 1], T=Fac)
