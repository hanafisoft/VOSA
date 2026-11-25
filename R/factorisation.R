#' Canonical factorization of a fitted model
#'
#' Writes a fitted multiblock model in the canonical factorization common to
#' multiblock methods (Bougeard et al., 2024):
#' \deqn{X_k = U D_k V_k' + E_k,}
#' where \eqn{U} (\eqn{n \times R}, \eqn{U'U = I}) holds the components,
#' \eqn{D_k} is diagonal with \eqn{d_{kr} = \lVert X_k' u_r \rVert},
#' \eqn{V_k = X_k' U D_k^{-1}} has unit-norm columns (the block loadings) and
#' \eqn{E_k = (I - UU')X_k} is the residual. When the blocks are scaled to unit
#' inertia (as by \code{\link{blocs_vers_S}} or \code{\link{bloc_W}}),
#' \eqn{d_{kr}^2} is the salience \eqn{a_k(u_r)}.
#'
#' @param fit an object of class \code{"vosa_fit"} (from \code{cpca},
#'   \code{indort}, \code{spa}, \code{vosa} or \code{statis}).
#' @param blocs list of data blocks; by default the normalized blocks stored in
#'   \code{fit} (available when \code{S} was built with \code{blocs_vers_S} or
#'   \code{bloc_W}).
#' @return An object of class \code{"vosa_factorisation"}, a list with
#'   \item{U}{components (\eqn{n \times R}).}
#'   \item{D}{\eqn{K \times R} matrix of the diagonal terms \eqn{d_{kr}}.}
#'   \item{V}{list of the \eqn{p_k \times R} loading matrices \eqn{V_k}.}
#'   \item{E}{list of the residual matrices \eqn{E_k}.}
#'   \item{part_restituee}{\eqn{K \times R} matrix of the share of the inertia
#'     of each block reproduced by each dimension, \eqn{d_{kr}^2 /
#'     \lVert X_k\rVert^2} (equal to the saliences for unit-inertia blocks).}
#'   \item{part_residuelle}{share of the inertia of each block left in
#'     \eqn{E_k}.}
#'   \item{erreur_reconstruction}{maximum absolute error of
#'     \eqn{U D_k V_k' + E_k - X_k} (numerical check).}
#' @references Bougeard, S., Peltier, C., Jaillais, B., Boulet, J.-C. and
#'   Hanafi, M. (2024). Benchmarking multiblock methods with canonical
#'   factorization. Chemometrics and Intelligent Laboratory Systems, 254, 105240.
#' @seealso \code{\link{plot_U}}, \code{\link{plot_D}}, \code{\link{plot_V}} for the plots of
#'   \eqn{U}, \eqn{D_k} and \eqn{V_k}.
#' @examples
#' sim <- simuler_blocs(seed = 3)
#' fit <- indort(blocs_vers_S(sim$X), R = 2, nstart = 3, seed = 3)
#' fc <- factorisation_canonique(fit)
#' fc$D^2      # equal to fit$saliences
#' plot(fc, type = "D")
#' @export
factorisation_canonique <- function(fit, blocs = NULL) {
  if (!inherits(fit, "vosa_fit")) stop("fit must be an object of class 'vosa_fit'.")
  if (is.null(blocs)) blocs <- fit$blocs
  if (is.null(blocs)) stop("No data blocks: build S with blocs_vers_S() or bloc_W(), or supply 'blocs'.")
  U <- fit$U
  if (length(blocs) != nrow(fit$saliences)) stop("The number of blocks does not match the fitted model.")
  R <- ncol(U)
  nm <- rownames(fit$saliences)
  D <- matrix(0, length(blocs), R, dimnames = list(nm, colnames(U)))
  V <- E <- vector("list", length(blocs))
  part <- D
  resid <- stats::setNames(numeric(length(blocs)), nm)
  err <- 0
  for (k in seq_along(blocs)) {
    Xk <- as.matrix(blocs[[k]])
    if (nrow(Xk) != nrow(U)) stop("Blocks and components must have the same number of rows.")
    P  <- crossprod(Xk, U)                       # p_k x R, columns X_k' u_r
    d  <- sqrt(colSums(P^2))
    Vk <- sweep(P, 2, ifelse(d > 0, d, 1), `/`)
    Vk[, d == 0] <- 0
    colnames(Vk) <- colnames(U)
    Ek <- Xk - U %*% t(P)                        # = (I - UU') X_k
    tot <- sum(Xk^2)
    D[k, ] <- d
    part[k, ] <- d^2 / tot
    resid[k] <- sum(Ek^2) / tot
    V[[k]] <- Vk
    E[[k]] <- Ek
    err <- max(err, max(abs(U %*% diag(d, nrow = R) %*% t(Vk) + Ek - Xk)))
  }
  names(V) <- names(E) <- nm
  structure(list(U = U, D = D, V = V, E = E, part_restituee = part,
                 part_residuelle = resid, erreur_reconstruction = err,
                 methode = fit$methode), class = "vosa_factorisation")
}
