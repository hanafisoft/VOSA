#' Ceilings and floor of the saliences
#'
#' The salience of block \eqn{k} on any unit-norm component cannot exceed the
#' share of inertia \eqn{\pi_{k1}} of its first principal direction (ceiling);
#' a random direction reproduces on average \eqn{1/(n-1)} of the inertia of a
#' centred block (floor).
#'
#' @inheritParams cpca
#' @return A data frame with the block names, the ceilings \eqn{\pi_{k1}} and
#'   the floor \eqn{1/(n-1)}.
#' @export
plafonds <- function(S) {
  S <- .preparer_S(S, TRUE)
  n <- nrow(S[[1]])
  data.frame(bloc = names(S), plafond = .plafonds(S), plancher = 1 / (n - 1),
             row.names = NULL)
}

#' Comparison of the components with the block principal components
#'
#' For each block and each dimension of a fitted model: salience, ceiling,
#' fraction of the ceiling reached, squared cosine with the first principal
#' component of the block, and whether the salience exceeds the level of the
#' dimension. A component close to the first principal component of one block
#' designates that block; close to those of several blocks, it is a direction
#' shared by this subset of blocks.
#'
#' @inheritParams cpca
#' @param fit an object of class \code{"vosa_fit"}.
#' @return A data frame (one row per block and dimension).
#' @export
diagnostic_blocs <- function(fit, S) {
  S <- .preparer_S(S, TRUE)
  U <- fit$U
  out <- list()
  for (k in seq_along(S)) {
    e  <- eigen(S[[k]], symmetric = TRUE)
    v1 <- e$vectors[, 1]
    for (r in seq_len(ncol(U))) {
      a <- fit$saliences[k, r]
      out[[length(out) + 1]] <- data.frame(
        bloc = names(S)[k], dimension = r, salience = a, plafond = e$values[1],
        fraction_plafond = a / e$values[1], cos2_cp1 = sum(U[, r] * v1)^2,
        au_dessus_niveau = a > mean(fit$saliences[, r]))
    }
  }
  do.call(rbind, out)
}

#' Heterogeneity of the block structures
#'
#' Global heterogeneity \eqn{H = K^{-1}\sum_k \lVert W_k - \bar W \rVert_F^2}
#' (interstructure inertia), its relative value, the RV-based index
#' \eqn{1 - \lambda_1(RV)/K}, and, if components are given, the share of
#' \eqn{H} reproduced by each of them, \eqn{V(u_r)/H}. For any orthonormal set
#' of components, these shares add up to at most 1.
#'
#' @inheritParams cpca
#' @param U optional \eqn{n \times R} matrix of orthonormal components, or a
#'   fitted object of class \code{"vosa_fit"}.
#' @return A list with \code{H}, \code{H_relative}, \code{indice_RV} and, if
#'   \code{U} is given, \code{part_restituee} (per dimension) and
#'   \code{part_cumulee}.
#' @export
heterogeneite <- function(S, U = NULL) {
  S <- .preparer_S(S, TRUE)
  K <- length(S)
  Wbar <- Reduce(`+`, S) / K
  H <- mean(vapply(S, function(Sk) sum((Sk - Wbar)^2), 0))
  Hrel <- H / mean(vapply(S, function(Sk) sum(Sk^2), 0))
  Sf <- lapply(S, function(Sk) Sk / sqrt(sum(Sk^2)))
  RV <- matrix(0, K, K)
  for (k in 1:K) for (l in 1:K) RV[k, l] <- sum(Sf[[k]] * Sf[[l]])
  out <- list(H = H, H_relative = Hrel,
              indice_RV = 1 - eigen(RV, symmetric = TRUE, only.values = TRUE)$values[1] / K)
  if (!is.null(U)) {
    if (inherits(U, "vosa_fit")) U <- U$U
    A <- .saliences(S, U)
    V <- colMeans(sweep(A, 2, colMeans(A))^2)
    out$part_restituee <- V / H
    out$part_cumulee <- cumsum(V) / H
  }
  out
}

#' Relative bounds on the contrast gain (corollary 2)
#'
#' From a path of solutions containing \eqn{\alpha = 1/2}, computes for each
#' \eqn{\alpha > 1/2} the relative gain of contrast
#' \eqn{\delta_V = \Delta V / V(u^*_{1/2})}, the relative loss of squared level
#' \eqn{\delta_{\bar a} = -\Delta \bar a^2 / \bar a(u^*_{1/2})^2}, and the
#' interval \eqn{[(1-\alpha)/\alpha\ \rho\,\delta_{\bar a},\ \rho\,\delta_{\bar a}]}
#' that contains \eqn{\delta_V}, with \eqn{\rho} computed at INDORT.
#'
#' @param chemin an object of class \code{"vosa_chemin"} (from
#'   \code{\link{chemin_spa}}) or its data frame \code{chemin}.
#' @return A data frame with \code{alpha}, \code{delta_contraste},
#'   \code{delta_niveau}, \code{borne_inf}, \code{borne_sup} and the value of
#'   \code{rho}.
#' @export
borne_relative <- function(chemin) {
  tab <- if (inherits(chemin, "vosa_chemin")) chemin$chemin else chemin
  i0 <- which(abs(tab$alpha - 0.5) < 1e-12)
  if (length(i0) != 1) stop("The path must contain alpha = 0.5.")
  tab <- tab[tab$alpha >= 0.5, ]
  i0 <- which(abs(tab$alpha - 0.5) < 1e-12)
  rho <- tab$niveau2[i0] / tab$contraste[i0]
  dl  <- (tab$niveau2[i0] - tab$niveau2) / tab$niveau2[i0]
  coef <- (1 - tab$alpha) / tab$alpha
  data.frame(alpha = tab$alpha,
             delta_contraste = (tab$contraste - tab$contraste[i0]) / tab$contraste[i0],
             delta_niveau = dl,
             borne_inf = coef * rho * dl,
             borne_sup = rho * dl,
             rho = rho, row.names = NULL)
}
