#' vosa: Variance-Optimized Saliences Analysis for multiblock data
#'
#' The package computes, for a list of scalar product matrices
#' \eqn{W_k = X_k X_k'} (one per block, scaled to unit trace), components
#' described by their saliences \eqn{a_k(u) = u'W_ku}. The family
#' \eqn{g_\alpha = (1-\alpha)\,\mathrm{level}^2 + \alpha\,\mathrm{contrast}}
#' contains CPCA (\code{\link{cpca}}, \eqn{\alpha = 0}), INDORT
#' (\code{\link{indort}}, \eqn{\alpha = 1/2}) and VOSA (\code{\link{vosa}},
#' \eqn{\alpha = 1}); \code{\link{spa}} solves any \eqn{\alpha} and
#' \code{\link{chemin_spa}} the whole path. \code{\link{statis}} is provided for
#' comparison, \code{\link{factorisation_canonique}} writes any fit in the
#' canonical factorization \eqn{X_k = U D_k V_k' + E_k} (plots: \code{\link{plot_U}},
#' \code{\link{plot_D}}, \code{\link{plot_V}}), and
#' \code{\link{plafonds}}, \code{\link{diagnostic_blocs}},
#' \code{\link{heterogeneite}} and \code{\link{borne_relative}} give the
#' diagnostics of the article. The script reproducing the results of the
#' article is installed in \code{system.file("scripts",
#' "reproduire_article.R", package = "vosa")}.
#'
#' @keywords internal
"_PACKAGE"
