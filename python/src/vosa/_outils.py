"""Internal helpers (not part of the public interface)."""

import numpy as np
import pandas as pd


class Resultat(dict):
    """Dictionary whose items can also be read as attributes (``fit.saliences``)."""

    def __getattr__(self, nom):
        try:
            return self[nom]
        except KeyError:
            raise AttributeError(nom) from None

    def __setattr__(self, nom, valeur):
        self[nom] = valeur

    def __dir__(self):
        return list(self.keys()) + list(super().__dir__())


class ListeBlocs(list):
    """List of matrices W_k carrying block names and, optionally, the normalized blocks."""

    def __init__(self, matrices, noms=None, blocs=None):
        super().__init__(matrices)
        self.noms = list(noms) if noms is not None else [f"Block{k + 1}" for k in range(len(matrices))]
        self.blocs = blocs


def _rng(seed):
    """Random generator from a seed, a Generator, or None."""
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


def _noms_et_blocs(S):
    """Block names and normalized data blocks attached to S, if any."""
    if isinstance(S, dict):
        noms = list(S.keys())
        mats = list(S.values())
        blocs = None
    else:
        mats = list(S)
        noms = getattr(S, "noms", None) or [f"Block{k + 1}" for k in range(len(mats))]
        blocs = getattr(S, "blocs", None)
    if blocs is None:
        bk = [getattr(m, "bloc", None) for m in mats]
        if all(b is not None for b in bk):
            blocs = bk
    return noms, mats, blocs


def _preparer_S(S, normaliser=True):
    """Check the list, symmetrize, and scale each W_k to unit trace."""
    noms, mats, blocs = _noms_et_blocs(S)
    if len(mats) < 2:
        raise ValueError("S must contain at least two matrices.")
    n = np.asarray(mats[0]).shape[0]
    out = []
    for k, m in enumerate(mats):
        Sk = np.array(m, dtype=float)
        if Sk.shape != (n, n):
            raise ValueError("All matrices must be square and of the same size.")
        Sk = (Sk + Sk.T) / 2
        if normaliser:
            tr = np.trace(Sk)
            if tr <= 0:
                raise ValueError(f"Matrix {k + 1} has a non-positive trace.")
            Sk = Sk / tr
        out.append(Sk)
    return ListeBlocs(out, noms, blocs)


def _saliences(S, U):
    """K x R array of saliences a_k(u_r) = u_r' W_k u_r."""
    return np.array([np.sum(U * (Sk @ U), axis=0) for Sk in S])


def _resume(A):
    """Level / contrast summary of each dimension (variance with divisor K)."""
    niveau = A.mean(axis=0)
    contraste = ((A - niveau) ** 2).mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = niveau ** 2 / contraste
    return pd.DataFrame({"dimension": np.arange(1, A.shape[1] + 1),
                         "niveau": niveau, "niveau2": niveau ** 2,
                         "contraste": contraste, "L2": niveau ** 2 + contraste,
                         "rho": rho})


def _fixer_signe(U):
    """Sign convention: the largest (absolute) entry of each column is positive."""
    U = np.array(U, dtype=float)
    idx = np.argmax(np.abs(U), axis=0)
    s = np.sign(U[idx, np.arange(U.shape[1])])
    s[s == 0] = 1
    return U * s


def _polar(Y):
    """Orthonormal factor of the polar decomposition."""
    u, _, vt = np.linalg.svd(Y, full_matrices=False)
    return u @ vt


def _plafonds(S):
    """First eigenvalue of each W_k (ceiling of the saliences)."""
    return np.array([np.linalg.eigvalsh(Sk)[-1] for Sk in S])


def _eigh_decroissant(M):
    """Eigenvalues and eigenvectors in decreasing order."""
    val, vec = np.linalg.eigh(M)
    return val[::-1], vec[:, ::-1]


def _depart_aleatoire(rng, n, R):
    q, _ = np.linalg.qr(rng.standard_normal((n, R)))
    return q


def _dims(R):
    return [f"Dim{r + 1}" for r in range(R)]


def _objet(classe, contenu, S, blocs):
    """Build a fitted object."""
    contenu["n"] = S[0].shape[0]
    contenu["K"] = len(S)
    contenu["plafonds"] = pd.Series(_plafonds(S), index=S.noms)
    contenu["blocs"] = blocs
    return classe(contenu)


def _fusion_etiquettes(x, y, noms, tol=0.06):
    """Merge the labels of points that (almost) coincide."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rx = np.ptp(np.r_[0, x]) or 1.0
    ry = np.ptp(np.r_[0, y]) or 1.0
    groupe = list(range(len(x)))
    for i in range(len(x)):
        for j in range(i):
            if abs(x[i] - x[j]) / rx < tol and abs(y[i] - y[j]) / ry < tol:
                groupe[i] = groupe[j]
                break
    out = []
    for g in dict.fromkeys(groupe):
        idx = [i for i in range(len(x)) if groupe[i] == g]
        out.append((x[idx].mean(), y[idx].mean(), " / ".join(noms[i] for i in idx)))
    return out
