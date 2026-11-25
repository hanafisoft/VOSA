#' Family g_alpha computed by Shifted Polar Ascent (SPA)
#'
#' Maximizes \eqn{G_\alpha(U) = R^{-1}\sum_r g_\alpha(u_r)}, with
#' \eqn{g_\alpha(u) = (1-\alpha)\,\bar a(u)^2 + \alpha\, V(u)}, where
#' \eqn{\bar a(u)} (level) and \eqn{V(u)} (contrast) are the mean and the
#' variance of the saliences \eqn{a_k(u) = u'W_ku}. The criterion is rewritten
#' as an INDORT criterion on the \eqn{K+1} blocks
#' \eqn{B_k = \sqrt{\alpha/K}(W_k - \bar W)} and
#' \eqn{B_{K+1} = \sqrt{1-\alpha}\,\bar W}, which are indefinite for
#' \eqn{\alpha > 0}. SPA computes at each iteration two candidates,
#' \eqn{\mathrm{polar}(M)} and \eqn{\mathrm{polar}(cU + M)}, and keeps the best;
#' the shift \eqn{c = \sum_j \lambda_{\max}^+(B_j)\lambda_{\min}^-(B_j)}
#' guarantees a monotone ascent.
#'
#' \code{vosa()} is \code{spa()} with \code{alpha = 1}: VOSA maximizes the
#' variance of the saliences (contrast).
#'
#' @inheritParams indort
#' @param alpha weight of the contrast, in [0, 1]: 0 = CPCA, 1/2 = INDORT,
#'   1 = VOSA.
#' @param U0 optional start (an \eqn{n \times R} matrix or a list of such
#'   matrices), e.g. the solution for a neighbouring value of \code{alpha}. The
#'   CPCA solution and \code{nstart} random starts are always added.
#' @param ... arguments passed from \code{vosa()} to \code{spa()}.
#' @return An object of class \code{"vosa_fit"}, a list with components
#'   \item{methode}{name of the method.}
#'   \item{alpha}{value of \eqn{\alpha}.}
#'   \item{U}{\eqn{n \times R} matrix of orthonormal components.}
#'   \item{saliences}{\eqn{K \times R} matrix of saliences \eqn{a_k(u_r)}.}
#'   \item{resume}{data frame with, for each dimension, the level, the squared
#'     level, the contrast, \eqn{L^2} = squared level + contrast, the ratio
#'     \eqn{\rho} = squared level / contrast, and \eqn{g_\alpha}.}
#'   \item{G_alpha}{value of the criterion.}
#'   \item{iterations, converge}{iterations and convergence of the best start.}
#'   \item{decalage}{shift \eqn{c}.}
#'   \item{criteres_departs}{criterion reached by each start (local maxima).}
#'   \item{n, K, plafonds}{sizes and ceilings \eqn{\pi_{k1}} of the blocks.}
#'   \item{blocs}{normalized data blocks if available, used by
#'     \code{\link{factorisation_canonique}}.}
#' @references Ballah Houssein, S., Galharret, J.-M. and Hanafi, M. (2026).
#'   VOSA: variance-optimized saliences analysis. Manuscript.
#' @examples
#' sim <- simuler_blocs(seed = 1)
#' S <- blocs_vers_S(sim$X)
#' fit <- vosa(S, R = 1, nstart = 5, seed = 1)
#' fit
#' factorisation_canonique(fit)
#' @export
spa <- function(S, alpha, R = 1, normaliser = TRUE, nstart = 10, U0 = NULL,
                maxit = 1000, tol = 1e-6, seed = NULL) {
  if (alpha < 0 || alpha > 1) stop("alpha must be in [0, 1].")
  blocs <- .extraire_blocs(S)
  S <- .preparer_S(S, normaliser)
  K <- length(S); n <- nrow(S[[1]])
  if (R < 1 || R > n) stop("R must be between 1 and n.")
  if (!is.null(seed)) set.seed(seed)
  B <- .blocs_B(S, alpha)
  departs <- if (is.null(U0)) list() else if (is.list(U0)) U0 else list(U0)
  departs <- c(departs, list(cpca(S, R, normaliser = FALSE)$U))
  if (nstart > 0) for (s in seq_len(nstart))
    departs[[length(departs) + 1]] <- qr.Q(qr(matrix(stats::rnorm(n * R), n, R)))
  res  <- lapply(departs, function(U) .spa_un_depart(B, U, maxit, tol))
  Fs   <- vapply(res, `[[`, 0, "F")
  best <- res[[which.max(Fs)]]
  U  <- best$U
  A  <- .saliences(S, U)
  rs <- .resume(A)
  g  <- (1 - alpha) * rs$niveau2 + alpha * rs$contraste
  o  <- order(g, decreasing = TRUE)
  U  <- .fixer_signe(U[, o, drop = FALSE])
  A  <- .saliences(S, U)
  rs <- .resume(A)
  rs$g_alpha <- (1 - alpha) * rs$niveau2 + alpha * rs$contraste
  dimnames(A) <- list(names(S), paste0("Dim", seq_len(R)))
  colnames(U) <- paste0("Dim", seq_len(R))
  nom <- if (alpha == 1) "VOSA" else if (alpha == 0.5) "INDORT (SPA)" else
    if (alpha == 0) "CPCA (SPA)" else sprintf("SPA (alpha = %g)", alpha)
  .objet(list(methode = nom,
              alpha = alpha,
              U = U,
              saliences = A,
              resume = rs,
              G_alpha = mean(rs$g_alpha),
              iterations = best$iterations,
              converge = best$converge,
              decalage = best$decalage,
              criteres_departs = sort(Fs / R, decreasing = TRUE)), S, blocs)
}

#' @rdname spa
#' @export
vosa <- function(S, R = 1, ...) spa(S, alpha = 1, R = R, ...)

.blocs_B <- function(S, alpha) {
  K <- length(S)
  Wbar <- Reduce(`+`, S) / K
  B <- lapply(S, function(Sk) sqrt(alpha / K) * (Sk - Wbar))
  B[[K + 1]] <- sqrt(1 - alpha) * Wbar
  B
}

.crit_B <- function(B, U) sum(.saliences(B, U)^2)

.spa_un_depart <- function(B, U, maxit, tol) {
  J  <- length(B)
  ev <- lapply(B, function(Bj) range(eigen(Bj, symmetric = TRUE, only.values = TRUE)$values))
  c_dec  <- sum(vapply(ev, function(r) max(r[2], 0) * max(-r[1], 0), 0))     # SPA-0
  normeB <- max(1, sum(vapply(ev, function(r) max(abs(r))^2, 0)))
  Fv <- .crit_B(B, U)
  converge <- FALSE
  it <- 0
  for (it in seq_len(maxit)) {
    b <- .saliences(B, U)                                                     # SPA-2
    M <- Reduce(`+`, lapply(seq_len(J), function(j) B[[j]] %*% U %*% diag(b[j, ], nrow = ncol(U))))
    sym <- crossprod(U, M); sym <- (sym + t(sym)) / 2
    if (sqrt(sum((M - U %*% sym)^2)) / normeB < tol) { converge <- TRUE; break } # SPA-5
    U0 <- .polar(M); Uc <- .polar(c_dec * U + M)                              # SPA-3
    F0 <- .crit_B(B, U0); Fc <- .crit_B(B, Uc)
    if (F0 >= Fc) { U <- U0; Fnew <- F0 } else { U <- Uc; Fnew <- Fc }        # SPA-4
    if (Fnew < Fv - 1e-10 * max(1, abs(Fv))) warning("The criterion decreased (should not happen).")
    Fv <- Fnew
  }
  list(U = U, F = Fv, iterations = it, converge = converge, decalage = c_dec)
}

#' Path of solutions over a grid of alpha
#'
#' Solves the family \eqn{g_\alpha} on a grid of values of \eqn{\alpha}, with an
#' increasing sweep (random starts plus the solution of the previous value) and
#' a decreasing sweep (the solution of the next value is kept if it improves the
#' criterion).
#'
#' @inheritParams spa
#' @param alphas grid of values of \eqn{\alpha}.
#' @return An object of class \code{"vosa_chemin"}: a list with
#'   \item{chemin}{data frame with, for each \eqn{\alpha}, the squared level, the
#'     contrast, \eqn{L^2} (means over the \code{R} dimensions) and
#'     \eqn{G_\alpha}.}
#'   \item{solutions}{list of the fitted objects.}
#' @examples
#' sim <- simuler_blocs(seed = 2)
#' ch <- chemin_spa(blocs_vers_S(sim$X), alphas = seq(0, 1, by = 0.25), nstart = 3, seed = 2)
#' ch
#' @export
chemin_spa <- function(S, alphas = seq(0, 1, by = 0.05), R = 1, normaliser = TRUE,
                       nstart = 10, maxit = 1000, tol = 1e-6, seed = NULL) {
  blocs <- .extraire_blocs(S)
  S <- .preparer_S(S, normaliser)
  if (!is.null(blocs)) attr(S, "blocs") <- blocs
  if (!is.null(seed)) set.seed(seed)
  alphas <- sort(unique(alphas))
  m <- length(alphas)
  sol <- vector("list", m)
  for (i in seq_len(m))
    sol[[i]] <- spa(S, alphas[i], R, normaliser = FALSE, nstart = nstart,
                    U0 = if (i > 1) sol[[i - 1]]$U else NULL, maxit = maxit, tol = tol)
  if (m > 1) for (i in (m - 1):1) {
    cand <- spa(S, alphas[i], R, normaliser = FALSE, nstart = 0,
                U0 = sol[[i + 1]]$U, maxit = maxit, tol = tol)
    if (cand$G_alpha > sol[[i]]$G_alpha + 1e-12) sol[[i]] <- cand
  }
  tab <- do.call(rbind, lapply(sol, function(s)
    data.frame(alpha     = s$alpha,
               niveau2   = mean(s$resume$niveau2),
               contraste = mean(s$resume$contraste),
               L2        = mean(s$resume$L2),
               G_alpha   = s$G_alpha)))
  structure(list(chemin = tab, solutions = sol, R = R), class = "vosa_chemin")
}
