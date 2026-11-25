#' Plots of the canonical factorization: components U, diagonal terms D_k and loadings V_k
#'
#' Graphical representations of the three elements of the canonical
#' factorization \eqn{X_k = U D_k V_k' + E_k}:
#' \code{plot_U} draws the individuals on the components \eqn{U} (one or two
#' dimensions, optionally coloured by groups); \code{plot_D} draws the diagonal
#' terms \eqn{d_{kr}} of the \eqn{D_k}, or their squares (the saliences, shares
#' of the inertia of each block reproduced by each dimension), as a bar chart of
#' blocks by dimensions; \code{plot_V} draws the loadings \eqn{V_k} of one block
#' (as a profile against the variables for one dimension, for instance a
#' spectrum, or as a scatter plot of the variables for two dimensions).
#' Labels follow the language chosen with \code{\link{vosa_langue}}.
#' \code{plot()} on an object of class \code{"vosa_factorisation"} calls one of
#' these functions according to \code{type}.
#'
#' @param x an object of class \code{"vosa_factorisation"} (from
#'   \code{\link{factorisation_canonique}}) or, for \code{plot_U} and
#'   \code{plot_D}, an object of class \code{"vosa_fit"}.
#' @param dims one or two dimensions to plot.
#' @param groupes optional factor (length \eqn{n}) used to colour the individuals.
#' @param etiquettes logical; write the names (or numbers) of the individuals or
#'   of the variables.
#' @param carre logical; plot \eqn{d_{kr}^2} (saliences) rather than \eqn{d_{kr}}.
#' @param bloc block whose loadings are plotted (number or name).
#' @param axe optional numeric abscissa for the variables of the block (e.g.
#'   wavelengths or chemical shifts), used when one dimension is plotted.
#' @param main title of the plot (by default, the method, and the block for
#'   \code{plot_V}).
#' @param type for \code{plot()}: \code{"D"}, \code{"U"} or \code{"V"}.
#' @param ... further graphical arguments passed to \code{plot} or
#'   \code{barplot}, or, for \code{plot()}, to the selected function.
#' @return The plotted values, invisibly.
#' @examples
#' sim <- simuler_blocs(pi_d = c(0.3, 0, 0, 0), seed = 5)
#' fit <- vosa(blocs_vers_S(sim$X), R = 2, nstart = 3, seed = 5)
#' fc <- factorisation_canonique(fit)
#' plot_D(fc)
#' plot_U(fc, groupes = factor(sim$t_d > 0, labels = c("t_d < 0", "t_d > 0")))
#' plot_V(fc, bloc = 1, dims = 1)
#' plot_V(fc, bloc = 1, dims = 1:2)
#' plot(fc, type = "D")
#' @name graphiques-factorisation
NULL

#' @rdname graphiques-factorisation
#' @export
plot_U <- function(x, dims = c(1, 2), groupes = NULL, etiquettes = FALSE, main = NULL, ...) {
  U <- x$U
  if (max(dims) > ncol(U)) stop("dims exceeds the number of components.")
  n <- nrow(U)
  noms <- if (is.null(rownames(U))) as.character(seq_len(n)) else rownames(U)
  if (is.null(groupes)) {
    col <- rep("grey20", n)
  } else {
    groupes <- as.factor(groupes)
    if (length(groupes) != n) stop("groupes must have one value per individual.")
    pal <- grDevices::hcl.colors(max(2, nlevels(groupes)), "Dark 3")
    col <- pal[as.integer(groupes)]
  }
  titre <- if (!is.null(main)) main else if (!is.null(x$methode)) x$methode else ""
  if (length(dims) == 1) {
    s <- U[, dims]
    plot(seq_len(n), s, pch = 19, col = col, xlab = .tr("individu"),
         ylab = paste(.tr("composante"), dims), main = titre, ...)
    graphics::abline(h = 0, lty = 3)
    if (etiquettes) graphics::text(seq_len(n), s, noms, pos = 3, cex = 0.7)
    out <- s
  } else {
    s <- U[, dims[1:2], drop = FALSE]
    plot(s[, 1], s[, 2], pch = 19, col = col, asp = 1,
         xlab = paste(.tr("composante"), dims[1]), ylab = paste(.tr("composante"), dims[2]),
         main = titre, ...)
    graphics::abline(h = 0, v = 0, lty = 3)
    if (etiquettes) graphics::text(s[, 1], s[, 2], noms, pos = 3, cex = 0.7)
    out <- s
  }
  if (!is.null(groupes))
    graphics::legend("topright", levels(groupes), pch = 19,
                     col = grDevices::hcl.colors(max(2, nlevels(groupes)), "Dark 3")[seq_len(nlevels(groupes))],
                     bty = "n", cex = 0.8)
  invisible(out)
}

#' @rdname graphiques-factorisation
#' @export
plot_D <- function(x, carre = TRUE, main = NULL, ...) {
  if (is.null(main)) main <- if (!is.null(x$methode)) x$methode else ""
  D <- if (inherits(x, "vosa_fit")) sqrt(pmax(x$saliences, 0)) else x$D
  M <- if (carre) D^2 else D
  pal <- grDevices::hcl.colors(max(2, ncol(M)), "Blues 3")[seq_len(ncol(M))]
  graphics::barplot(t(M), beside = TRUE, col = pal, las = 2,
                    ylab = if (carre) .tr("d2") else .tr("d"),
                    main = main, ...)
  graphics::legend("topright", colnames(M), fill = pal, bty = "n", cex = 0.8)
  invisible(M)
}

#' @rdname graphiques-factorisation
#' @export
plot_V <- function(x, bloc = 1, dims = 1, axe = NULL, etiquettes = NULL, main = NULL, ...) {
  if (!inherits(x, "vosa_factorisation"))
    stop("x must be an object of class 'vosa_factorisation' (see factorisation_canonique).")
  Vk <- x$V[[bloc]]
  nom_bloc <- if (is.character(bloc)) bloc else names(x$V)[bloc]
  p <- nrow(Vk)
  noms <- if (is.null(rownames(Vk))) as.character(seq_len(p)) else rownames(Vk)
  if (is.null(etiquettes)) etiquettes <- p <= 40
  if (is.null(main)) main <- paste(nom_bloc, "-", x$methode)
  if (length(dims) == 1) {
    v <- Vk[, dims]
    ab <- if (is.null(axe)) seq_len(p) else axe
    if (length(ab) != p) stop("axe must have one value per variable of the block.")
    plot(ab, v, type = if (p > 40) "l" else "h", lwd = if (p > 40) 1.5 else 3,
         xlab = if (is.null(axe)) .tr("variable") else "", ylab = paste(.tr("chargement"), dims),
         main = main, ...)
    graphics::abline(h = 0, lty = 3)
    if (etiquettes && p <= 40) graphics::text(ab, v, noms, pos = ifelse(v >= 0, 3, 1), cex = 0.7)
    out <- v
  } else {
    v <- Vk[, dims[1:2], drop = FALSE]
    lim <- max(abs(v)) * 1.15
    plot(v[, 1], v[, 2], pch = 19, col = "grey30", asp = 1, xlim = c(-lim, lim), ylim = c(-lim, lim),
         xlab = paste(.tr("chargement"), dims[1]), ylab = paste(.tr("chargement"), dims[2]),
         main = main, ...)
    graphics::abline(h = 0, v = 0, lty = 3)
    if (etiquettes) graphics::text(v[, 1], v[, 2], noms, pos = 3, cex = 0.7)
    out <- v
  }
  invisible(out)
}

#' @rdname graphiques-factorisation
#' @export
plot.vosa_factorisation <- function(x, type = c("D", "U", "V"), ...) {
  type <- match.arg(type)
  switch(type, D = plot_D(x, ...), U = plot_U(x, ...), V = plot_V(x, ...))
}
