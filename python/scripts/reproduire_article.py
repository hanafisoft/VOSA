"""Reproduction of the results of the article (Python version of
inst/scripts/reproduire_article.R).

    Ballah Houssein, Galharret and Hanafi (2026),
    "VOSA : analyse des saliences a variance optimisee dans le modele INDSCAL".

Data (not distributed with the package; put the files in DOSSIER_DONNEES):
  - Datafusionnetdata.mat : metabolomic data fusion, colorectal cancer / adenoma
      (Bro et al., 2013, Metabolomics 9, 3-8); requires scipy.
  - CRC_metadonnees.csv   : status and pairs of the 94 patients (rows i and i + 47
      are the adenoma and the cancer of the same pair).
  - potato_multibloc_full.csv : Potato data (Thygesen et al., 2001), as in the
      R package multiblock (all blocks, one column per variable).
  - breast_TCGA.csv : breast.TCGA training data of the R package mixOmics
      (columns 1-200 mRNA, 201-384 miRNA, 385-526 proteins).

Usage:  python reproduire_article.py [dossier_donnees] [dossier_resultats] [fr|en]

The random starts differ from those of R: the solutions of INDORT and VOSA are
the same up to the sign and to local optima; CPCA, STATIS and the canonical
factorization are identical.
"""

import os
import pickle
import platform
import sys

import numpy as np
import pandas as pd

import vosa as vs

DOSSIER_DONNEES = sys.argv[1] if len(sys.argv) > 1 else "donnees"
DOSSIER_RESULTATS = sys.argv[2] if len(sys.argv) > 2 else "resultats_python"
vs.vosa_langue(sys.argv[3] if len(sys.argv) > 3 else "fr")
GRAINE = 2026
os.makedirs(DOSSIER_RESULTATS, exist_ok=True)


def fichier(nom):
    return os.path.join(DOSSIER_DONNEES, nom)


def sortie(nom):
    return os.path.join(DOSSIER_RESULTATS, nom)


def ecrire_csv(tab, nom):
    vs.traduire(tab).to_csv(sortie(nom), index=False)


# ------------------------------------------------------------------ data ----
jeux = {}

# CRC (metabolomic data fusion): biomarkers and fluorescence scaled, NMR centred only
from scipy.io import loadmat  # noqa: E402

mat = loadmat(fichier("Datafusionnetdata.mat"), squeeze_me=True, struct_as_record=False)
Xcrc = np.asarray(mat["Data"].data, float)          # 94 x 476
jeux["CRC"] = {"Biomarqueurs": vs.bloc_W(Xcrc[:, 0:2]),
               "Fluorescence": vs.bloc_W(Xcrc[:, 2:21]),
               "RMN CPMG": vs.bloc_W(Xcrc[:, 21:222], sc=False),
               "RMN NOESY": vs.bloc_W(Xcrc[:, 222:476], sc=False)}
meta_crc = pd.read_csv(fichier("CRC_metadonnees.csv"))

# Potato: blocks centred only
pot = pd.read_csv(fichier("potato_multibloc_full.csv"), index_col=0)
blocs_pot = ["Chemical", "Compression", "NIRraw", "NIRcooked", "CPMGraw",
             "CPMGcooked", "FIDraw", "FIDcooked", "Sensory"]
jeux["Potato"] = {b: vs.bloc_W(pot.loc[:, [c for c in pot.columns if c.startswith(b + ".")]], sc=False)
                  for b in blocs_pot}

# breast.TCGA: variables scaled
tcga = pd.read_csv(fichier("breast_TCGA.csv")).values
jeux["breast.TCGA"] = {"ARNm": vs.bloc_W(tcga[:, 0:200]),
                       "miARN": vs.bloc_W(tcga[:, 200:384]),
                       "proteines": vs.bloc_W(tcga[:, 384:526])}


# ------------------------------------------------------------------ fits ----
def methodes(S, R):
    rng = np.random.default_rng(GRAINE)
    return {"CPCA": vs.cpca(S, R=R),
            "INDORT": vs.indort(S, R=R, nstart=20, seed=rng),
            "VOSA": vs.vosa(S, R=R, nstart=20, seed=rng),
            "STATIS": vs.statis(S, R=R, normalisation="trace")}


fits = {j: {"R1": methodes(S, 1), "R2": methodes(S, 2)} for j, S in jeux.items()}


# ------------------------------------------------------ tables of saliences ----
def tab_saliences(r):
    lignes = []
    for j, S in jeux.items():
        f = fits[j][f"R{r}"]
        tab = pd.DataFrame({"jeu": j, "bloc": list(S.keys()), "PI": f["CPCA"].plafonds.values})
        for m, fit in f.items():
            for d in range(r):
                tab[m if r == 1 else f"{m}{d + 1}"] = fit.saliences.iloc[:, d].values
        lignes.append(tab)
    return pd.concat(lignes, ignore_index=True)


ecrire_csv(tab_saliences(1), "saliences_R1.csv")
ecrire_csv(tab_saliences(2), "saliences_R2.csv")

nc = pd.DataFrame([{"jeu": j, "R": r, "methode": m,
                    "niveau2": fit.resume.niveau2.mean(), "contraste": fit.resume.contraste.mean(),
                    "rho_indort": fits[j]["R1"]["INDORT"].resume.rho[0]}
                   for j in jeux for r in (1, 2) for m, fit in fits[j][f"R{r}"].items()])
ecrire_csv(nc, "niveau_contraste.csv")

# ---------------------------------------------------------- heterogeneity ----
het = []
for j, S in jeux.items():
    for r in (1, 2):
        for m, fit in fits[j][f"R{r}"].items():
            h = vs.heterogeneite(S, fit)
            het.append({"jeu": j, "R": r, "methode": m, "H": h.H, "H_relative": h.H_relative,
                        "indice_RV": h.indice_RV, "part_restituee": h.part_restituee.sum()})
ecrire_csv(pd.DataFrame(het), "heterogeneite.csv")

diag_vosa = pd.concat([vs.diagnostic_blocs(fits[j]["R1"]["VOSA"], S).assign(jeu=j)
                       for j, S in jeux.items()], ignore_index=True)
ecrire_csv(diag_vosa, "diagnostic_blocs_VOSA.csv")

# --------------------------------------------- path and relative bounds ----
bornes = []
for j, S in jeux.items():
    ch = vs.chemin_spa(S, alphas=np.round(np.arange(0.5, 1.00001, 0.01), 10), nstart=10, seed=GRAINE)
    bornes.append(vs.borne_relative(ch).assign(jeu=j))
bornes = pd.concat(bornes, ignore_index=True)
ecrire_csv(bornes, "borne_relative.csv")

# --------------------------------------- CRC: paired t tests (47 pairs) ----
N_P = 47
assert (meta_crc["paire"].values[:N_P] == meta_crc["paire"].values[N_P:2 * N_P]).all()


def test_crc(f, libelle):
    from scipy import stats
    out = []
    for m, fit in f.items():
        s = np.sqrt(N_P) * fit.U[:, 0]
        tt = stats.ttest_rel(s[N_P:2 * N_P], s[:N_P])
        out.append({"blocs": libelle, "methode": m, "t": tt.statistic, "ddl": N_P - 1, "p": tt.pvalue})
    return out


crc4 = test_crc(fits["CRC"]["R1"], "4 blocs")
rng = np.random.default_rng(GRAINE)
S3 = {k: v for k, v in jeux["CRC"].items() if k != "Biomarqueurs"}
f3 = {"CPCA": vs.cpca(S3), "INDORT": vs.indort(S3, nstart=20, seed=rng),
      "VOSA": vs.vosa(S3, nstart=20, seed=rng), "STATIS": vs.statis(S3, normalisation="trace")}
crc3 = test_crc(f3, "sans biomarqueurs")
ecrire_csv(pd.DataFrame(crc4 + crc3), "crc_test_apparie.csv")

# ------------------------------------------ canonical factorization ----
fc = {j: {r: {m: vs.factorisation_canonique(fit) for m, fit in fr.items()}
          for r, fr in fj.items()} for j, fj in fits.items()}
with open(sortie("factorisation_canonique.pkl"), "wb") as fh:
    pickle.dump(fc, fh)
for j in fc:
    for r in (1, 2):
        D = pd.concat({m: f.D for m, f in fc[j][f"R{r}"].items()}, axis=1)
        D.columns = [f"{m}_{d}" for m, d in D.columns]
        ecrire_csv(D.reset_index(names="bloc"), f"factorisation_D_{j}_R{r}.csv")
    V = fc[j]["R1"]["VOSA"].V
    ecrire_csv(pd.concat([pd.DataFrame({"bloc": b, "variable": np.arange(1, len(v) + 1),
                                        "chargement": v[:, 0]}) for b, v in V.items()]),
               f"chargements_VOSA_{j}.csv")

# --------------------------------------------------------------- figures ----
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(2, 3, figsize=(13, 8))
for i, r in enumerate((1, 2)):
    for k, j in enumerate(jeux):
        vs.plot_niveau_contraste(fits[j][f"R{r}"], ax=ax[i, k],
                                 main=f"{j}, R = {r} (rho = {fits[j]['R1']['INDORT'].resume.rho[0]:.2f})")
fig.tight_layout()
fig.savefig(sortie("niveau_contraste.pdf"))

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for k, j in enumerate(jeux):
    vs.plot_borne_relative(bornes[bornes.jeu == j], ax=ax[k], main=j)
fig.tight_layout()
fig.savefig(sortie("borne_relative.pdf"))

fig, ax = plt.subplots(1, 3, figsize=(13, 4.5))
for k, j in enumerate(jeux):
    fits[j]["R1"]["VOSA"].plot(ax=ax[k], main=j)
fig.tight_layout()
fig.savefig(sortie("saliences_VOSA.pdf"))

fig, ax = plt.subplots(3, 3, figsize=(13, 13))
for k, j in enumerate(jeux):
    f = fc[j]["R2"]["VOSA"]
    vs.plot_D(f, ax=ax[k, 0], main=f"{j} - D_k (saliences)")
    vs.plot_U(f, groupes=meta_crc["statut"].values if j == "CRC" else None, ax=ax[k, 1], main=f"{j} - U")
    vs.plot_V(f, bloc=int(np.argmax(f.D.iloc[:, 0].values)), dims=1, ax=ax[k, 2])
fig.tight_layout()
fig.savefig(sortie("factorisation_VOSA.pdf"))
plt.close("all")

with open(sortie("versions.txt"), "w") as fh:
    fh.write(f"Python {platform.python_version()}, vosa {vs.__version__}, numpy {np.__version__}, "
             f"pandas {pd.__version__}, matplotlib {matplotlib.__version__}\n")
print("Results written in", os.path.abspath(DOSSIER_RESULTATS))
