# vosa — Variance-Optimized Saliences Analysis for multiblock data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24531478.svg)](https://doi.org/10.5281/zenodo.24531478)

R package accompanying the article *VOSA : analyse des saliences à variance optimisée dans le modèle INDSCAL* (Ballah Houssein, Galharret and Hanafi, 2026).

Multiblock components are described by their **saliences**, the share of each block's inertia they reproduce. The family

  g_alpha = (1 − alpha) · level² + alpha · contrast

contains CPCA (alpha = 0), INDORT / orthogonal INDSCAL (alpha = 1/2, simultaneous ComDim) and **VOSA** (alpha = 1), which maximizes the variance of the saliences of a component. All criteria are computed by the Shifted Polar Ascent (SPA) algorithm, whose ascent is guaranteed for indefinite blocks.

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("hanafisoft/VOSA")

# or from the source archive
install.packages("vosa_1.0.0.tar.gz", repos = NULL, type = "source")
```

Source code: https://github.com/hanafisoft/VOSA (issues: https://github.com/hanafisoft/VOSA/issues).

### Python

A Python version with the same functions and outputs is in the folder [`python/`](python):

```bash
pip install "git+https://github.com/hanafisoft/VOSA#subdirectory=python"
```

## Example

```r
library(vosa)
sim <- simuler_blocs(pi_d = c(0.3, 0, 0, 0), seed = 1)   # factor specific to block 1
S   <- blocs_vers_S(sim$X)
fit <- vosa(S, R = 1, nstart = 10, seed = 1)
fit                                   # saliences, level and contrast
fc <- factorisation_canonique(fit)    # X_k = U D_k V_k' + E_k
plot_D(fc); plot_U(fc); plot_V(fc, bloc = 1)
diagnostic_blocs(fit, S)              # comparison with the block principal components
ch  <- chemin_spa(S, alphas = seq(0.5, 1, by = 0.05), seed = 1)
borne_relative(ch)                    # contrast gain and bounds (corollary 2)
```

## Main functions

| Function | Role |
|---|---|
| `blocs_vers_S()`, `bloc_W()` | scalar product matrices from data blocks (unit inertia) |
| `cpca()`, `indort()`, `vosa()`, `spa()`, `statis()` | methods |
| `chemin_spa()` | path of solutions over a grid of alpha |
| `factorisation_canonique()` | canonical factorization of any fit |
| `plot_U()`, `plot_D()`, `plot_V()` | plots of the components U, of the diagonal terms D_k and of the loadings V_k |
| `plafonds()`, `diagnostic_blocs()`, `heterogeneite()`, `borne_relative()` | diagnostics |
| `simuler_blocs()` | simulated data |

## Language of the results

Printed results, plot labels and tables are in French by default. To switch:

```r
vosa_langue("en")          # or options(vosa.langue = "en"), e.g. in .Rprofile
fit                        # printed in English
write.csv(traduire(fit$resume), "summary.csv", row.names = FALSE)
vosa_langue("fr")
```

Element names of the objects (`fit$saliences`, `fit$resume`, ...) stay in French, so a script runs unchanged in both languages; `traduire()` returns a copy with translated names.

## Reproducing the article

```r
file.edit(system.file("scripts", "reproduire_article.R", package = "vosa"))
```

The script lists the data files needed and writes all tables, figures and canonical factorizations of the article.

## Citation

Please cite:

Ballah Houssein, S., Galharret, J.-M. and Hanafi, M. (2026). vosa: Variance-Optimized Saliences Analysis for Multiblock Data. R and Python packages, version 1.0.0. https://github.com/hanafisoft/VOSA. doi:[10.5281/zenodo.24531478](https://doi.org/10.5281/zenodo.24531478)

In R: `citation("vosa")`.

## License

MIT.
