"""Preparation of the blocks: scalar product matrices W_k = X_k X_k'."""

import numpy as np

from ._outils import ListeBlocs


class MatriceW(np.ndarray):
    """Matrix W_k carrying the normalized block it comes from (attribute ``bloc``)."""

    def __new__(cls, W, bloc=None):
        obj = np.asarray(W, dtype=float).view(cls)
        obj.bloc = bloc
        return obj

    def __array_finalize__(self, obj):
        self.bloc = getattr(obj, "bloc", None)


def _normaliser_bloc(X, reduire=False):
    """Centre (and optionally scale) the columns, then scale the block to unit inertia."""
    Y = np.array(X, dtype=float)
    if Y.ndim == 1:
        Y = Y[:, None]
    Y = Y - Y.mean(axis=0)
    if reduire:
        sd = Y.std(axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            Y = np.where(sd > 0, Y / np.where(sd > 0, sd, 1), 0.0)
    tot = np.sqrt(np.sum(Y ** 2))
    if tot <= 0:
        raise ValueError("A block has zero inertia.")
    return Y / tot


def blocs_vers_S(X, reduire=False):
    """Scalar product matrices of a list (or dict) of data blocks.

    Each block is centred (and, if ``reduire``, its variables scaled to unit
    variance), then scaled to unit inertia, so that ``W_k = X_k X_k'`` has unit
    trace. The normalized blocks are kept (attribute ``blocs``) for
    :func:`factorisation_canonique`.

    Parameters
    ----------
    X : list or dict of 2-D arrays (or DataFrames), same number of rows.
    reduire : bool
        Scale the variables to unit variance before normalizing the blocks.

    Returns
    -------
    ListeBlocs
        List of the W_k, with attributes ``noms`` (block names) and ``blocs``.
    """
    if isinstance(X, dict):
        noms, mats = list(X.keys()), list(X.values())
    else:
        mats = list(X)
        noms = [f"Block{k + 1}" for k in range(len(mats))]
    if len(mats) < 2:
        raise ValueError("X must contain at least two blocks.")
    n = np.asarray(mats[0]).shape[0]
    if any(np.asarray(m).shape[0] != n for m in mats):
        raise ValueError("All blocks must have the same number of rows.")
    Y = [_normaliser_bloc(m, reduire) for m in mats]
    return ListeBlocs([y @ y.T for y in Y], noms, Y)


def bloc_W(X, sc=True):
    """Scalar product matrix of one block (centred, by default scaled, unit inertia).

    The normalized block is kept in the attribute ``bloc``, so that a list of
    such matrices can be passed to the methods and to
    :func:`factorisation_canonique`.
    """
    Y = _normaliser_bloc(X, reduire=sc)
    return MatriceW(Y @ Y.T, bloc=Y)
