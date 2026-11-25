sim <- simuler_blocs(n = 25, K = 4, p = 40, pi_d = c(0.3, 0, 0, 0), seed = 11)
S <- blocs_vers_S(sim$X)

test_that("blocks are scaled to unit trace and carry the normalized data", {
  expect_equal(unname(vapply(S, function(W) sum(diag(W)), 0)), rep(1, 4))
  expect_false(is.null(attr(S, "blocs")))
})

test_that("saliences lie between 0 and the ceiling of each block", {
  fit <- vosa(S, R = 2, nstart = 3, seed = 1)
  expect_true(all(fit$saliences >= -1e-10))
  expect_true(all(fit$saliences <= fit$plafonds + 1e-10))
})

test_that("level, contrast and norm follow the properties of the family (R = 1)", {
  f0 <- cpca(S)
  f5 <- indort(S, nstart = 10, seed = 2)
  f1 <- vosa(S, nstart = 10, seed = 2)
  expect_gte(f0$resume$niveau2 + 1e-10, f5$resume$niveau2)
  expect_gte(f5$resume$niveau2 + 1e-10, f1$resume$niveau2)
  expect_gte(f1$resume$contraste + 1e-10, f5$resume$contraste)
  expect_gte(f5$resume$contraste + 1e-10, f0$resume$contraste)
  expect_gte(f5$resume$L2 + 1e-10, f0$resume$L2)
  expect_gte(f5$resume$L2 + 1e-10, f1$resume$L2)
})

test_that("spa with alpha = 1/2 reaches the INDORT criterion", {
  f5 <- indort(S, nstart = 10, seed = 3)
  s5 <- spa(S, alpha = 0.5, nstart = 10, seed = 3)
  expect_equal(s5$resume$L2, f5$resume$L2, tolerance = 1e-6)
})

test_that("VOSA recovers the block-specific factor", {
  f1 <- vosa(S, nstart = 10, seed = 4)
  expect_gt(sum(f1$U[, 1] * sim$t_d)^2, 0.8)
})

test_that("the canonical factorization reconstructs the blocks", {
  fit <- indort(S, R = 2, nstart = 3, seed = 5)
  fc <- factorisation_canonique(fit)
  expect_lt(fc$erreur_reconstruction, 1e-10)
  expect_equal(unname(fc$D^2), unname(fit$saliences), tolerance = 1e-10)
  expect_equal(unname(fc$part_residuelle + rowSums(fc$part_restituee)), rep(1, 4), tolerance = 1e-10)
})

test_that("the shares of heterogeneity of an orthonormal basis add up to at most 1", {
  U <- eigen(Reduce(`+`, S), symmetric = TRUE)$vectors
  h <- heterogeneite(S, U)
  expect_lte(sum(h$part_restituee), 1 + 1e-10)
})

test_that("the contrast gain lies within the bounds of corollary 2", {
  ch <- chemin_spa(S, alphas = seq(0.5, 1, by = 0.1), nstart = 5, seed = 6)
  b <- borne_relative(ch)
  expect_true(all(b$delta_contraste >= b$borne_inf - 1e-6))
  expect_true(all(b$delta_contraste <= b$borne_sup + 1e-6))
})

test_that("bloc_W gives the same results as blocs_vers_S", {
  S2 <- lapply(sim$X, bloc_W, sc = FALSE)
  expect_equal(cpca(S2)$saliences, cpca(S)$saliences, tolerance = 1e-10)
  expect_false(is.null(factorisation_canonique(cpca(S2))))
})

test_that("the plots of U, D_k and V_k run", {
  fit <- vosa(S, R = 2, nstart = 3, seed = 7)
  fc <- factorisation_canonique(fit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot_D(fc))
  expect_silent(plot_D(fit, carre = FALSE))
  expect_silent(plot_U(fc, groupes = rep(c("a", "b"), length.out = nrow(fc$U))))
  expect_silent(plot_U(fit, dims = 1))
  expect_silent(plot_V(fc, bloc = 1, dims = 1))
  expect_silent(plot_V(fc, bloc = 2, dims = 1:2))
  expect_silent(plot(fc, type = "V", bloc = 1))
  expect_equal(unname(plot_D(fc)), unname(fit$saliences), tolerance = 1e-10)
})

test_that("language of the results", {
  ancienne <- vosa_langue()
  on.exit(vosa_langue(ancienne))
  fit <- vosa(S, R = 1, nstart = 2, seed = 3)
  vosa_langue("fr")
  expect_equal(vosa_langue(), "fr")
  expect_identical(names(traduire(fit)), names(fit))
  expect_output(print(fit), "Niveau et contraste")
  vosa_langue("en")
  expect_equal(vosa_langue(), "en")
  expect_output(print(fit), "Level and contrast")
  tr <- traduire(fit)
  expect_true(all(c("method", "summary", "converged") %in% names(tr)))
  expect_true("contrast" %in% names(traduire(fit$resume)))
  expect_equal(tr$saliences, fit$saliences)
  expect_error(vosa_langue("de"))
})
