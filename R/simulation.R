#' Simulated multiblock data with common and block-specific factors
#'
#' Generator used in the simulation study: \code{n} individuals, \code{K} blocks
#' of \code{p} variables; factors \eqn{t_{c1}, t_{c2}, \dots} present in all
#' blocks with shares of inertia \code{pi_c}, and a factor \eqn{t_d} present in
#' block \eqn{k} with share \code{pi_d[k]}. Factors are centred and orthonormal;
#' each block has its own random loadings; the remaining inertia is Gaussian
#' noise. Each block is centred and scaled to unit inertia.
#'
#' @param n number of individuals.
#' @param K number of blocks.
#' @param p number of variables per block.
#' @param pi_c shares of inertia of the common factors.
#' @param pi_d vector of length \code{K}: share of inertia of \eqn{t_d} in each
#'   block.
#' @param seed optional random seed.
#' @return A list with \code{X} (list of the \code{K} blocks), \code{t_d} and
#'   \code{T} (matrix of all factors, the last column being \eqn{t_d}).
#' @examples
#' sim <- simuler_blocs(pi_d = c(0.3, 0, 0, 0), seed = 4)
#' fit <- vosa(blocs_vers_S(sim$X), nstart = 3, seed = 4)
#' sum(fit$U[, 1] * sim$t_d)^2      # squared cosine with the specific factor
#' @export
simuler_blocs <- function(n = 30, K = 4, p = 100, pi_c = c(0.20, 0.15, 0.10),
                          pi_d = c(0.30, 0.30, 0, 0), seed = NULL) {
  if (length(pi_d) != K) stop("pi_d must have length K.")
  if (any(sum(pi_c) + pi_d >= 1)) stop("The shares of inertia must add up to less than 1.")
  if (!is.null(seed)) set.seed(seed)
  q <- length(pi_c) + 1
  Fac <- matrix(stats::rnorm(n * q), n, q)
  Fac <- scale(Fac, center = TRUE, scale = FALSE)
  Fac <- qr.Q(qr(Fac))
  X <- vector("list", K)
  for (k in seq_len(K)) {
    parts <- c(pi_c, pi_d[k])
    Xk <- matrix(0, n, p)
    for (j in seq_len(q)) {
      v <- stats::rnorm(p); v <- v / sqrt(sum(v^2))
      Xk <- Xk + sqrt(parts[j]) * outer(Fac[, j], v)
    }
    E <- matrix(stats::rnorm(n * p), n, p)
    E <- scale(E, center = TRUE, scale = FALSE); E <- E / sqrt(sum(E^2))
    Xk <- Xk + sqrt(1 - sum(parts)) * E
    Xk <- scale(Xk, center = TRUE, scale = FALSE)
    attr(Xk, "scaled:center") <- NULL
    X[[k]] <- Xk / sqrt(sum(Xk^2))
  }
  names(X) <- paste0("Block", seq_len(K))
  list(X = X, t_d = Fac[, q], T = Fac)
}
