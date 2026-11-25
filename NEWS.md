# vosa 1.0.0

* First release, accompanying the article "VOSA : analyse des saliences à variance optimisée dans le modèle INDSCAL".
* Methods: `cpca()`, `indort()`, `spa()` (family g_alpha, Shifted Polar Ascent), `vosa()`, `chemin_spa()`, `statis()`.
* Canonical factorization of any fit: `factorisation_canonique()`, with plots of the components `plot_U()`, of the diagonal terms `plot_D()` and of the loadings `plot_V()` (also `plot()` on the factorization).
* Diagnostics: `plafonds()`, `diagnostic_blocs()`, `heterogeneite()`, `borne_relative()`.
* Simulation generator `simuler_blocs()`; reproduction script in `inst/scripts/reproduire_article.R`.
* Results in French or English: `vosa_langue("fr")` / `vosa_langue("en")` (option `vosa.langue`, French by default) sets the language of printed output, plot labels and tables; `traduire()` returns a copy of an object or table with English or French names, for export. Element names of the objects stay in French, so scripts work in both languages.
