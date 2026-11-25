#' Consensus PCA (CPCA)
#'
#' Maximizes the sum of the saliences of the components,
#' \eqn{\max_{U'U=I} \sum_r \sum_k a_k(u_r) = K\,\mathrm{tr}(U'\bar W U)}: the
#' components are the first eigenvectors of \eqn{\bar W = K^{-1}\sum_k W_k}
#' (point \eqn{\alpha = 0} of the family \eqn{g_\alpha}).
#'
#' @param S list of \eqn{n \times n} scalar product matrices \eqn{W_k}
#'   (see \code{\link{blocs_vers_S}}).
#' @param R number of components.
#' @param normaliser logical; scale each \eqn{W_k} to unit trace.
#' @return An object of class \code{"vosa_fit"}; see \code{\link{vosa}}.
#' @export
cpca <- function(S, R = 1, normaliser = TRUE) {
  blocs <- .extraire_blocs(S)
  S <- .preparer_S(S, normaliser)
  K <- length(S); n <- nrow(S[[1]])
  if (R < 1 || R > n) stop("R must be between 1 and n.")
  Wbar <- Reduce(`+`, S) / K
  e    <- eigen(Wbar, symmetric = TRUE)
  U    <- .fixer_signe(e$vectors[, seq_len(R), drop = FALSE])
  A    <- .saliences(S, U)
  dimnames(A) <- list(names(S), paste0("Dim", seq_len(R)))
  colnames(U) <- paste0("Dim", seq_len(R))
  .objet(list(methode = "CPCA",
              alpha = 0,
              U = U,
              saliences = A,
              resume = .resume(A),
              valeurs_propres = e$values[seq_len(R)],
              critere = K * sum(e$values[seq_len(R)]),
              Wbar = Wbar), S, blocs)
}
