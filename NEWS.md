# vosa 1.0.1

* DOI of the archived version 1.0.0 (Zenodo): 10.5281/zenodo.24531478, added to the citation (`citation("vosa")`, `CITATION.cff`) and to the READMEs.

# vosa 1.0.0

* First release, accompanying the article "VOSA : analyse des saliences à variance optimisée dans le modèle INDSCAL".
* Methods: `cpca()`, `indort()`, `spa()` (family g_alpha, Shifted Polar Ascent), `vosa()`, `chemin_spa()`, `statis()`.
* Canonical factorization of any fit: `factorisation_canonique()`, with plots of the components `plot_U()`, of the diagonal terms `plot_D()` and of the loadings `plot_V()` (also `plot()` on the factorization).
* Diagnostics: `plafonds()`, `diagnostic_blocs()`, `heterogeneite()`, `borne_relative()`.
* Simulation generator `simuler_blocs()`; reproduction script in `inst/scripts/reproduire_article.R`.
* Results in French or English: `vosa_langue("fr")` / `vosa_langue("en")` (option `vosa.langue`, French by default) sets the language of printed output, plot labels and tables; `traduire()` returns a copy of an object or table with English or French names, for export. Element names of the objects stay in French, so scripts work in both languages.
* Plots: margins wide enough for axis titles and block names, labels of individuals (n <= 50) and of the variables with the largest loadings (`n_etiquettes`), labels of coinciding methods merged in `plot_niveau_contraste()` (e.g. "CPCA / STATIS").
* Python version of the package in `python/` (same functions and outputs; `pip install "git+https://github.com/hanafisoft/VOSA#subdirectory=python"`), with its tests and a script reproducing the results of the article.
