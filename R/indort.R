#' Orthogonal INDSCAL (INDORT)
#'
#' Simultaneous extraction of \code{R} orthonormal components maximizing
#' \eqn{F(U) = \sum_r \sum_k a_k(u_r)^2} (orthogonal INDSCAL model
#' \eqn{W_k \approx U \Lambda_k U'}). For \code{R = 1}, this is the first
#' dimension of ComDim / CCSWA. Point \eqn{\alpha = 1/2} of the family
#' \eqn{g_\alpha}. The algorithm is the polar ascent \eqn{U \leftarrow
#' \mathrm{polar}(M)}, whose ascent is guaranteed since the \eqn{W_k} are
#' positive semi-definite.
#'
#' @inheritParams cpca
#' @param nstart number of starts: the CPCA solution plus \code{nstart - 1}
#'   random orthonormal matrices.
#' @param maxit maximum number of iterations per start.
#' @param tol tolerance on the relative norm of the Riemannian gradient.
#' @param seed optional random seed.
#' @return An object of class \code{"vosa_fit"}; see \code{\link{vosa}}.
#' @export
indort <- function(S, R = 1, normaliser = TRUE, nstart = 10,
                   maxit = 1000, tol = 1e-8, seed = NULL) {
  blocs <- .extraire_blocs(S)
  S <- .preparer_S(S, normaliser)
  K <- length(S); n <- nrow(S[[1]])
  if (R < 1 || R > n) stop("R must be between 1 and n.")
  if (!is.null(seed)) set.seed(seed)
  departs <- list(cpca(S, R, normaliser = FALSE)$U)
  if (nstart > 1) for (s in 2:nstart)
    departs[[s]] <- qr.Q(qr(matrix(stats::rnorm(n * R), n, R)))
  res  <- lapply(departs, function(U0) .indort_un_depart(S, U0, maxit, tol))
  crit <- vapply(res, `[[`, 0, "critere")
  best <- res[[which.max(crit)]]
  U <- best$U
  A <- .saliences(S, U)
  o <- order(colSums(A^2), decreasing = TRUE)
  U <- .fixer_signe(U[, o, drop = FALSE])
  A <- .saliences(S, U)
  dimnames(A) <- list(names(S), paste0("Dim", seq_len(R)))
  colnames(U) <- paste0("Dim", seq_len(R))
  .objet(list(methode = "INDORT",
              alpha = 0.5,
              U = U,
              saliences = A,
              resume = .resume(A),
              critere = best$critere,
              iterations = best$iterations,
              converge = best$converge,
              criteres_departs = sort(crit, decreasing = TRUE)), S, blocs)
}

.indort_un_depart <- function(S, U, maxit, tol) {
  K  <- length(S)
  A  <- .saliences(S, U)
  Fv <- sum(A^2)
  normeS <- max(1, sum(vapply(S, function(Sk)
    max(abs(eigen(Sk, symmetric = TRUE, only.values = TRUE)$values))^2, 0)))
  converge <- FALSE
  it <- 0
  for (it in seq_len(maxit)) {
    M <- Reduce(`+`, lapply(seq_len(K), function(k) S[[k]] %*% U %*% diag(A[k, ], nrow = ncol(U))))
    sym  <- crossprod(U, M); sym <- (sym + t(sym)) / 2
    grad <- sqrt(sum((M - U %*% sym)^2)) / normeS
    if (grad < tol) { converge <- TRUE; break }
    U    <- .polar(M)
    A    <- .saliences(S, U)
    Fnew <- sum(A^2)
    if (Fnew < Fv - 1e-12 * max(1, abs(Fv)))
      warning("The criterion decreased: check that the W_k are positive semi-definite.")
    Fv <- Fnew
  }
  list(U = U, critere = Fv, iterations = it, converge = converge)
}
