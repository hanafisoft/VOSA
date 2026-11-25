#' STATIS
#'
#' Interstructure (scalar products, or RV coefficients, between the
#' \eqn{W_k}), weights given by the first eigenvector, compromise
#' \eqn{W_c = \sum_k \gamma_k W_k}, and components given by the first
#' eigenvectors of the compromise. Saliences are always computed on the
#' unit-trace \eqn{W_k}, so that they are comparable with the other methods.
#'
#' @inheritParams cpca
#' @param normalisation \code{"RV"} (classical STATIS: each \eqn{W_k} has unit
#'   Frobenius norm for the interstructure) or \code{"trace"} (unit trace, as
#'   for CPCA and INDORT).
#' @return An object of class \code{"vosa_fit"} with, in addition,
#'   \code{poids}, \code{poids_somme1}, \code{interstructure}, \code{RV},
#'   \code{part_RV} (\eqn{\lambda_1(RV)/K}), \code{compromis} and
#'   \code{valeurs_propres_compromis}.
#' @export
statis <- function(S, R = 1, normalisation = c("RV", "trace")) {
  normalisation <- match.arg(normalisation)
  blocs <- .extraire_blocs(S)
  S <- .preparer_S(S, normaliser = TRUE)
  K <- length(S); n <- nrow(S[[1]])
  if (R < 1 || R > n) stop("R must be between 1 and n.")
  Sf <- if (normalisation == "RV") lapply(S, function(Sk) Sk / sqrt(sum(Sk^2))) else S
  P <- matrix(0, K, K, dimnames = list(names(S), names(S)))
  for (k in 1:K) for (l in k:K) P[k, l] <- P[l, k] <- sum(Sf[[k]] * Sf[[l]])
  nrm <- sqrt(diag(P))
  RV  <- P / outer(nrm, nrm)
  ei  <- eigen(P, symmetric = TRUE)
  gamma <- abs(ei$vectors[, 1])
  names(gamma) <- names(S)
  Wc <- Reduce(`+`, Map(`*`, Sf, gamma))
  e  <- eigen(Wc, symmetric = TRUE)
  U  <- .fixer_signe(e$vectors[, seq_len(R), drop = FALSE])
  A  <- .saliences(S, U)
  dimnames(A) <- list(names(S), paste0("Dim", seq_len(R)))
  colnames(U) <- paste0("Dim", seq_len(R))
  eRV <- eigen(RV, symmetric = TRUE, only.values = TRUE)$values
  .objet(list(methode = sprintf("STATIS (%s)", normalisation),
              alpha = NA,
              U = U,
              saliences = A,
              resume = .resume(A),
              poids = gamma,
              poids_somme1 = gamma / sum(gamma),
              interstructure = P,
              RV = RV,
              part_RV = eRV[1] / K,
              compromis = Wc,
              valeurs_propres_compromis = e$values[seq_len(R)]), S, blocs)
}
