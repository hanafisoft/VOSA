## Internal helpers (not exported)

## Collect the normalized data blocks attached to a list of matrices, if any:
## either attr(S, "blocs") (set by blocs_vers_S) or attr(S[[k]], "bloc") (set by bloc_W).
.extraire_blocs <- function(S) {
  b <- attr(S, "blocs")
  if (!is.null(b)) return(b)
  bk <- lapply(S, function(Sk) attr(Sk, "bloc"))
  if (all(!vapply(bk, is.null, logical(1)))) {
    if (!is.null(names(S))) names(bk) <- names(S)
    return(bk)
  }
  NULL
}

## Check the list, symmetrize, and scale each W_k to unit trace.
.preparer_S <- function(S, normaliser = TRUE) {
  if (!is.list(S) || length(S) < 2) stop("S must be a list of at least two matrices.")
  n <- nrow(as.matrix(S[[1]]))
  nm <- names(S)
  out <- vector("list", length(S))
  for (k in seq_along(S)) {
    Sk <- as.matrix(S[[k]])
    attr(Sk, "bloc") <- NULL
    if (nrow(Sk) != n || ncol(Sk) != n) stop("All matrices must be square and of the same size.")
    Sk <- (Sk + t(Sk)) / 2
    if (normaliser) {
      tr <- sum(diag(Sk))
      if (tr <= 0) stop(sprintf("Matrix %d has a non-positive trace.", k))
      Sk <- Sk / tr
    }
    out[[k]] <- Sk
  }
  names(out) <- if (is.null(nm)) paste0("Block", seq_along(S)) else nm
  out
}

## K x R matrix of saliences a_k(u_r) = u_r' W_k u_r
.saliences <- function(S, U) {
  matrix(vapply(S, function(Sk) colSums(U * (Sk %*% U)), numeric(ncol(U))),
         nrow = length(S), byrow = TRUE)
}

## Level / contrast summary for each dimension
.resume <- function(A) {
  niveau    <- colMeans(A)
  contraste <- colMeans(sweep(A, 2, niveau)^2)      # variance with divisor K
  data.frame(dimension = seq_len(ncol(A)),
             niveau    = niveau,
             niveau2   = niveau^2,
             contraste = contraste,
             L2        = niveau^2 + contraste,       # mean of squared saliences
             rho       = niveau^2 / contraste)
}

## Sign convention: the largest (absolute) entry of each column is positive
.fixer_signe <- function(U) {
  s <- apply(U, 2, function(u) sign(u[which.max(abs(u))]))
  s[s == 0] <- 1
  sweep(U, 2, s, `*`)
}

## Orthonormal factor of the polar decomposition
.polar <- function(Y) {
  sv <- svd(Y)
  sv$u %*% t(sv$v)
}

## First eigenvalue of each W_k (ceiling of the saliences)
.plafonds <- function(S) {
  vapply(S, function(Sk) eigen(Sk, symmetric = TRUE, only.values = TRUE)$values[1], numeric(1))
}

## Build a fitted object
.objet <- function(lst, S, blocs) {
  lst$n <- nrow(S[[1]])
  lst$K <- length(S)
  lst$plafonds <- .plafonds(S)
  lst$blocs <- blocs
  class(lst) <- "vosa_fit"
  lst
}
