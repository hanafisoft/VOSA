#' Print, summary and plot methods
#'
#' The display follows the language chosen with \code{\link{vosa_langue}}.
#' \code{print} shows the method, the table of saliences and the level/contrast
#' summary; \code{summary} adds the ceilings of the blocks; \code{plot} draws the
#' saliences of one dimension as a bar chart, with the ceilings of the blocks
#' as points.
#'
#' @param x,object an object of class \code{"vosa_fit"}, \code{"vosa_chemin"}
#'   or \code{"vosa_factorisation"}.
#' @param digits number of digits printed.
#' @param dim dimension to plot.
#' @param main title of the plot (by default, the method and the dimension).
#' @param ... further arguments passed to \code{barplot}.
#' @return \code{print} and \code{plot} return \code{x} invisibly;
#'   \code{summary} returns a list.
#' @name methodes-vosa
NULL

#' @rdname methodes-vosa
#' @export
print.vosa_fit <- function(x, digits = 3, ...) {
  cat(x$methode, ": n =", x$n, ", K =", x$K, ", R =", ncol(x$U), "\n\n")
  cat(.tr("saliences"), "\n")
  print(round(x$saliences, digits))
  cat("\n", .tr("resume"), "\n", sep = "")
  print(traduire(round(x$resume, digits)), row.names = FALSE)
  invisible(x)
}

#' @rdname methodes-vosa
#' @export
summary.vosa_fit <- function(object, ...) {
  tab <- cbind(object$plafonds, object$saliences)
  colnames(tab)[1] <- .tr("plafond")
  out <- list(methode = object$methode, saliences = tab, resume = object$resume,
              plancher = 1 / (object$n - 1))
  print(round(tab, 3))
  cat(.tr("plancher"), round(out$plancher, 3), "\n")
  invisible(out)
}

#' @rdname methodes-vosa
#' @export
plot.vosa_fit <- function(x, dim = 1, main = NULL, ...) {
  if (is.null(main)) main <- paste(x$methode, "-", .tr("dimension"), dim)
  a <- x$saliences[, dim]
  pos <- graphics::barplot(a, ylim = c(0, max(c(a, x$plafonds)) * 1.1),
                           ylab = .tr("salience"), las = 2,
                           main = main, ...)
  graphics::points(pos, x$plafonds, pch = 25, bg = "grey40")
  graphics::abline(h = mean(a), lty = 2)
  graphics::legend("topright", c(.tr("plafond"), .tr("niveau")), pch = c(25, NA), lty = c(NA, 2),
                   pt.bg = "grey40", bty = "n")
  invisible(x)
}

#' @rdname methodes-vosa
#' @export
print.vosa_chemin <- function(x, digits = 4, ...) {
  cat(.tr("chemin"), x$R, "\n")
  print(traduire(round(x$chemin, digits)), row.names = FALSE)
  invisible(x)
}

#' @rdname methodes-vosa
#' @export
print.vosa_factorisation <- function(x, digits = 3, ...) {
  cat(.tr("fc_titre"), " (", x$methode, ")\n\n", sep = "")
  cat(.tr("fc_D"), "\n")
  print(round(x$D, digits))
  cat("\n", .tr("fc_part"), "\n", sep = "")
  print(round(x$part_restituee, digits))
  cat("\n", .tr("fc_resid"), "\n", sep = "")
  print(round(x$part_residuelle, digits))
  cat("\n", .tr("fc_err"), " ", signif(x$erreur_reconstruction, 3), "\n", sep = "")
  invisible(x)
}

#' Level/contrast plane
#'
#' Plots, for several fitted models, the mean squared level against the mean
#' contrast of their dimensions.
#'
#' @param fits a named list of objects of class \code{"vosa_fit"}.
#' @param ... further arguments passed to \code{plot}.
#' @return The data frame of the plotted points, invisibly.
#' @export
plot_niveau_contraste <- function(fits, ...) {
  if (is.null(names(fits))) names(fits) <- vapply(fits, `[[`, "", "methode")
  tab <- data.frame(methode = names(fits),
                    niveau2 = vapply(fits, function(f) mean(f$resume$niveau2), 0),
                    contraste = vapply(fits, function(f) mean(f$resume$contraste), 0))
  plot(tab$niveau2, tab$contraste, pch = 19, xlab = .tr("niveau2"),
       ylab = .tr("contraste"), xlim = c(0, max(tab$niveau2) * 1.15),
       ylim = c(0, max(tab$contraste) * 1.15), ...)
  graphics::text(tab$niveau2, tab$contraste, tab$methode, pos = 3, cex = 0.8)
  invisible(tab)
}

#' Plot of the relative bounds
#'
#' @param b data frame returned by \code{\link{borne_relative}}.
#' @param ... further arguments passed to \code{plot}.
#' @return \code{b}, invisibly.
#' @export
plot_borne_relative <- function(b, ...) {
  plot(b$alpha, b$delta_contraste, type = "n", xlab = expression(alpha),
       ylab = expression(delta[V]), ylim = range(c(0, b$borne_sup, b$delta_contraste)), ...)
  graphics::polygon(c(b$alpha, rev(b$alpha)), c(b$borne_inf, rev(b$borne_sup)),
                    col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  graphics::lines(b$alpha, b$delta_contraste, lwd = 2)
  invisible(b)
}
