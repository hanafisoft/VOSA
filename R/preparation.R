#' Scalar product matrices from data blocks
#'
#' Centres (and optionally scales) the columns of each block, scales each block
#' to unit total inertia, and returns the list of scalar product matrices
#' \eqn{W_k = X_k X_k'}. The normalized blocks are attached as attribute
#' \code{"blocs"}, so that fitted models can be returned in canonical
#' factorization form.
#'
#' @param X list of data blocks (matrices or data frames) with the same rows.
#' @param reduire logical; scale the variables to unit variance before
#'   computing the scalar products.
#' @return A named list of \eqn{n \times n} matrices of unit trace, with
#'   attribute \code{"blocs"} (the normalized blocks).
#' @export
blocs_vers_S <- function(X, reduire = FALSE) {
  if (!is.list(X) || length(X) < 2) stop("X must be a list of at least two blocks.")
  n <- nrow(as.matrix(X[[1]]))
  if (any(vapply(X, function(x) nrow(as.matrix(x)), 1L) != n))
    stop("All blocks must have the same number of rows.")
  Y <- lapply(X, .normaliser_bloc, reduire = reduire)
  S <- lapply(Y, tcrossprod)
  nm <- names(X)
  if (is.null(nm)) nm <- paste0("Block", seq_along(X))
  names(S) <- nm
  names(Y) <- nm
  attr(S, "blocs") <- Y
  S
}

#' Scalar product matrix of one block
#'
#' Centres (and by default scales) the columns of a block, scales it to unit
#' total inertia and returns \eqn{W = X X'}. The normalized block is attached as
#' attribute \code{"bloc"}.
#'
#' @param X a data block (matrix or data frame).
#' @param sc logical; scale the variables to unit variance.
#' @return An \eqn{n \times n} matrix of unit trace with attribute \code{"bloc"}.
#' @export
bloc_W <- function(X, sc = TRUE) {
  Y <- .normaliser_bloc(X, reduire = sc)
  W <- tcrossprod(Y)
  attr(W, "bloc") <- Y
  W
}

.normaliser_bloc <- function(X, reduire = FALSE) {
  Y <- scale(as.matrix(X), center = TRUE, scale = reduire)
  Y[is.na(Y)] <- 0                      # constant variables after scaling
  attr(Y, "scaled:center") <- NULL
  attr(Y, "scaled:scale") <- NULL
  tot <- sqrt(sum(Y^2))
  if (tot <= 0) stop("A block has zero inertia.")
  Y / tot
}
