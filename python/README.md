# vosa (Python) — Variance-Optimized Saliences Analysis for multiblock data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24531478.svg)](https://doi.org/10.5281/zenodo.24531478)

Python version of the R package `vosa` (same repository), accompanying the article *VOSA : analyse des saliences à variance optimisée dans le modèle INDSCAL* (Ballah Houssein, Galharret and Hanafi, 2026). Same functions, same arguments and same outputs as the R package.

Multiblock components are described by their **saliences**, the share of each block's inertia they reproduce. The family

  g_alpha = (1 − alpha) · level² + alpha · contrast

contains CPCA (alpha = 0), INDORT / orthogonal INDSCAL (alpha = 1/2, simultaneous ComDim) and **VOSA** (alpha = 1), which maximizes the variance of the saliences of a component. All criteria are computed by the Shifted Polar Ascent (SPA) algorithm, whose ascent is guaranteed for indefinite blocks.

## Installation

```bash
pip install "git+https://github.com/hanafisoft/VOSA#subdirectory=python"
```

Requirements: Python ≥ 3.9, numpy, pandas, matplotlib (scipy only to read the `.mat` file of the reproduction script).

## Example

```python
import vosa as vs

sim = vs.simuler_blocs(n=30, K=4, p=50, pi_d=(0.3, 0, 0, 0), seed=1)
S = vs.blocs_vers_S(sim.X)                 # W_k = X_k X_k', unit trace

fit = vs.vosa(S, R=2, nstart=10, seed=1)   # also cpca, indort, spa(alpha=...), statis
fit                                        # saliences, level and contrast by dimension
fit.summary()                              # with the ceilings of the blocks

fc = vs.factorisation_canonique(fit)       # X_k = U D_k V_k' + E_k
vs.plot_D(fc); vs.plot_U(fc); vs.plot_V(fc, bloc=0, dims=(1, 2))

ch = vs.chemin_spa(S, alphas=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0], nstart=5, seed=1)
vs.borne_relative(ch)                      # bounds of corollary 2
```

Results are attributes and keys at the same time (`fit.saliences` or `fit["saliences"]`). Saliences, summaries and the terms of the canonical factorization are pandas DataFrames; `U` and the loadings `V_k` are numpy arrays.

## Main functions

| Function | Role |
|---|---|
| `blocs_vers_S`, `bloc_W` | Scalar product matrices of the blocks (centred, optionally scaled, unit inertia) |
| `cpca`, `indort`, `spa`, `vosa`, `chemin_spa` | The family g_alpha and its path of solutions |
| `statis` | STATIS, for comparison |
| `factorisation_canonique` | Canonical factorization X_k = U D_k V_k' + E_k |
| `plot_saliences`, `plot_U`, `plot_D`, `plot_V` | Plots of the saliences, components, diagonal terms and loadings |
| `plot_niveau_contraste`, `plot_borne_relative` | Level/contrast plane, relative bounds |
| `plafonds`, `diagnostic_blocs`, `heterogeneite`, `borne_relative` | Diagnostics |
| `simuler_blocs` | Simulation generator |
| `vosa_langue`, `traduire` | Results in French or English |

## Language of the results

Printed results, plot labels and tables are in French by default:

```python
vs.vosa_langue("en")       # or the environment variable VOSA_LANGUE=en
fit                        # printed in English
vs.traduire(fit.resume).to_csv("summary.csv", index=False)
vs.vosa_langue("fr")
```

Attribute names (`fit.saliences`, `fit.resume`, ...) stay in French, so a script runs unchanged in both languages.

## Reproducing the article

```bash
python scripts/reproduire_article.py <data folder> <results folder> [fr|en]
```

The data files and their sources are listed at the top of the script. The tables of saliences (R = 1 and R = 2), the paired tests and the bounds of the article are reproduced to the third decimal. The random starts differ from those of R, so INDORT and VOSA may differ in sign or, in rare cases, reach another local optimum; CPCA, STATIS and the canonical factorization are identical.

## Tests

```bash
pip install -e ".[test]"
pytest tests
```

## Citation

Ballah Houssein, S., Galharret, J.-M. and Hanafi, M. (2026). vosa: Variance-Optimized Saliences Analysis for Multiblock Data. R and Python packages, version 1.0.0. https://github.com/hanafisoft/VOSA. doi:[10.5281/zenodo.24531478](https://doi.org/10.5281/zenodo.24531478)

## License

MIT. Ballah Houssein S., Galharret J.-M., Hanafi M. — StatSC (Oniris, INRAE), Oniris VetAgroBio Nantes, France.
