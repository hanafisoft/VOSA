"""Tests of the Python package vosa (same properties as the R tests)."""

import numpy as np
import pandas as pd
import pytest

import vosa as vs

SIM = vs.simuler_blocs(n=25, K=4, p=40, pi_d=(0.3, 0, 0, 0), seed=11)
S = vs.blocs_vers_S(SIM.X)


def test_blocs_trace_unite_et_blocs_conserves():
    assert np.allclose([np.trace(W) for W in S], 1)
    assert S.blocs is not None and len(S.blocs) == 4
    assert S.noms == ["Block1", "Block2", "Block3", "Block4"]


def test_saliences_entre_0_et_plafond():
    fit = vs.vosa(S, R=2, nstart=3, seed=1)
    A = fit.saliences.values
    assert np.all(A >= -1e-10)
    assert np.all(A <= fit.plafonds.values[:, None] + 1e-10)


def test_proprietes_de_la_famille_R1():
    f0 = vs.cpca(S)
    f5 = vs.indort(S, nstart=10, seed=2)
    f1 = vs.vosa(S, nstart=10, seed=2)
    n0, n5, n1 = (f.resume.niveau2[0] for f in (f0, f5, f1))
    c0, c5, c1 = (f.resume.contraste[0] for f in (f0, f5, f1))
    l0, l5, l1 = (f.resume.L2[0] for f in (f0, f5, f1))
    assert n0 + 1e-10 >= n5 >= n1 - 1e-10
    assert c1 + 1e-10 >= c5 >= c0 - 1e-10
    assert l5 + 1e-10 >= max(l0, l1)


def test_spa_alpha_demi_atteint_indort():
    f5 = vs.indort(S, nstart=10, seed=3)
    s5 = vs.spa(S, alpha=0.5, nstart=10, seed=3)
    assert s5.resume.L2[0] == pytest.approx(f5.resume.L2[0], rel=1e-6)


def test_vosa_retrouve_le_facteur_propre():
    f1 = vs.vosa(S, nstart=10, seed=4)
    assert np.dot(f1.U[:, 0], SIM.t_d) ** 2 > 0.8


def test_spa_ne_baisse_jamais():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        for a in (0, 0.3, 0.5, 0.8, 1):
            vs.spa(S, alpha=a, R=2, nstart=5, seed=5)


def test_factorisation_canonique_reconstruit_les_blocs():
    fit = vs.indort(S, R=2, nstart=3, seed=5)
    fc = vs.factorisation_canonique(fit)
    assert fc.erreur_reconstruction < 1e-10
    assert np.allclose(fc.D.values ** 2, fit.saliences.values, atol=1e-10)
    assert np.allclose(fc.part_residuelle + fc.part_restituee.sum(axis=1), 1)


def test_part_heterogeneite_au_plus_1():
    W = sum(S)
    U = np.linalg.eigh(W)[1][:, ::-1]
    h = vs.heterogeneite(S, U)
    assert h.part_restituee.sum() <= 1 + 1e-10


def test_bornes_du_corollaire_2():
    ch = vs.chemin_spa(S, alphas=np.round(np.arange(0.5, 1.0001, 0.1), 10), nstart=5, seed=6)
    b = vs.borne_relative(ch)
    assert np.all(b.delta_contraste >= b.borne_inf - 1e-6)
    assert np.all(b.delta_contraste <= b.borne_sup + 1e-6)


def test_bloc_W_equivaut_a_blocs_vers_S():
    S2 = [vs.bloc_W(X, sc=False) for X in SIM.X.values()]
    assert np.allclose(vs.cpca(S2).saliences.values, vs.cpca(S).saliences.values)
    assert vs.factorisation_canonique(vs.cpca(S2)) is not None


def test_statis():
    f = vs.statis(S, R=2)
    assert f.poids_somme1.sum() == pytest.approx(1)
    assert np.allclose(np.diag(f.RV.values), 1)


def test_graphiques():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fit = vs.vosa(S, R=2, nstart=3, seed=7)
    fc = vs.factorisation_canonique(fit)
    fig, ax = plt.subplots(2, 4)
    fit.plot(ax=ax[0, 0])
    vs.plot_D(fc, ax=ax[0, 1])
    vs.plot_U(fc, groupes=SIM.t_d > 0, ax=ax[0, 2])
    vs.plot_U(fit, dims=1, ax=ax[0, 3])
    vs.plot_V(fc, bloc=0, dims=1, ax=ax[1, 0])
    vs.plot_V(fc, bloc="Block2", dims=(1, 2), ax=ax[1, 1])
    vs.plot_niveau_contraste({"CPCA": vs.cpca(S), "VOSA": fit}, ax=ax[1, 2])
    fc.plot(type="D", ax=ax[1, 3])
    plt.close(fig)
    assert np.allclose(vs.plot_D(fc).values, fit.saliences.values)
    plt.close("all")


def test_langue():
    ancienne = vs.vosa_langue()
    try:
        fit = vs.vosa(S, nstart=2, seed=3)
        vs.vosa_langue("fr")
        assert "Niveau et contraste" in repr(fit)
        assert list(vs.traduire(fit).keys()) == list(fit.keys())
        vs.vosa_langue("en")
        assert "Level and contrast" in repr(fit)
        tr = vs.traduire(fit)
        assert {"method", "summary", "converged"} <= set(tr.keys())
        assert "contrast" in vs.traduire(fit.resume).columns
        assert tr["saliences"].equals(fit.saliences)
        with pytest.raises(ValueError):
            vs.vosa_langue("de")
    finally:
        vs.vosa_langue(ancienne)
