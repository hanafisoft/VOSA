## =============================================================================
## Reproduction of the results of the article
##   Ballah Houssein, Galharret and Hanafi (2026),
##   "VOSA : analyse des saliences a variance optimisee dans le modele INDSCAL".
##
## Package vosa 1.0.0. Script location after installation:
##   system.file("scripts", "reproduire_article.R", package = "vosa")
##
## Data (not distributed with the package; put the files in 'dossier_donnees'):
##   - Datafusionnetdata.mat : metabolomic data fusion, colorectal cancer / adenoma
##       (Bro et al., 2013, Metabolomics 9, 3-8; "Metabolomic data fusion" data set,
##       University of Copenhagen); requires the package R.matlab.
##   - CRC_metadonnees.csv   : status and pairs of the 94 patients (rows i and i + 47
##       are the adenoma and the cancer of the same pair).
##   - potato_multibloc_full.csv : Potato data (Thygesen et al., 2001), as in the
##       R package multiblock (all blocks, one column per variable).
##   - breast_TCGA.csv : breast.TCGA training data of the R package mixOmics
##       (columns 1-200 mRNA, 201-384 miRNA, 385-526 proteins).
##
## Language: set 'langue' below to "fr" or "en"; printed results, plot labels and
## CSV column names follow it (see ?vosa_langue). Objects keep French element names.
##
## Outputs (in 'dossier_resultats'):
##   - saliences_R1.csv, saliences_R2.csv   : tables of saliences (ceilings PI included)
##   - niveau_contraste.csv                 : squared level and contrast of each method
##   - heterogeneite.csv                    : H and share of H reproduced by each method
##   - diagnostic_blocs_VOSA.csv            : VOSA component versus block principal components
##   - borne_relative.csv                   : relative contrast gain and bounds (corollary 2)
##   - crc_test_apparie.csv                 : paired t tests cancer / adenoma (CRC)
##   - factorisation_canonique.rds          : canonical factorizations X_k = U D_k V_k' + E_k
##                                            of every fit (all data sets, methods, R = 1, 2)
##   - factorisation_D_<data>_R<r>.csv      : diagonal terms d_kr of the D_k for each method
##   - chargements_VOSA_<data>.csv          : loadings V_k of the first VOSA dimension
##   - figures *.pdf (including factorisation_VOSA.pdf: D_k, U and V_k), and sessionInfo.txt
## =============================================================================

library(vosa)

dossier_donnees   <- "donnees"          # <- adapt
dossier_resultats <- "resultats"        # <- adapt
refaire_simulation <- FALSE             # TRUE: also rerun the simulation study (long)
graine <- 2026
langue <- "fr"                          # "fr" or "en": language of printed results,
                                        #   plot labels and CSV column names
vosa_langue(langue)
dir.create(dossier_resultats, showWarnings = FALSE, recursive = TRUE)
fichier <- function(...) file.path(dossier_donnees, ...)
sortie  <- function(...) file.path(dossier_resultats, ...)
ecrire_csv <- function(x, file, ...) utils::write.csv(traduire(x), file, ...)

## ---------------------------------------------------------------- data ----
jeux <- list()

## CRC (metabolomic data fusion): biomarkers and fluorescence scaled, NMR centred only
if (!requireNamespace("R.matlab", quietly = TRUE)) stop("Install the package R.matlab.")
mat <- R.matlab::readMat(fichier("Datafusionnetdata.mat"))
Xcrc <- mat$Data$data                    # 94 x 476 (as in the analysis script of the article)
jeux$CRC <- list(
  Biomarqueurs = bloc_W(Xcrc[, 1:2]),
  Fluorescence = bloc_W(Xcrc[, 3:21]),
  `RMN CPMG`   = bloc_W(Xcrc[, 22:222], sc = FALSE),
  `RMN NOESY`  = bloc_W(Xcrc[, 223:476], sc = FALSE))
meta_crc <- utils::read.csv(fichier("CRC_metadonnees.csv"))

## Potato: blocks centred only
pot <- utils::read.csv(fichier("potato_multibloc_full.csv"), row.names = 1)
blocs_pot <- c("Chemical", "Compression", "NIRraw", "NIRcooked", "CPMGraw",
               "CPMGcooked", "FIDraw", "FIDcooked", "Sensory")
jeux$Potato <- lapply(blocs_pot, function(b)
  bloc_W(pot[, startsWith(names(pot), paste0(b, ".")), drop = FALSE], sc = FALSE))
names(jeux$Potato) <- blocs_pot

## breast.TCGA: variables scaled
tcga <- utils::read.csv(fichier("breast_TCGA.csv"))
jeux$breast.TCGA <- list(ARNm      = bloc_W(tcga[, 1:200]),
                         miARN     = bloc_W(tcga[, 201:384]),
                         proteines = bloc_W(tcga[, 385:526]))

## ------------------------------------------------------------- fits ----
methodes <- function(S, R) {
  set.seed(graine)
  list(CPCA   = cpca(S, R = R),
       INDORT = indort(S, R = R, nstart = 20),
       VOSA   = vosa(S, R = R, nstart = 20),
       STATIS = statis(S, R = R, normalisation = "trace"))
}
fits <- lapply(jeux, function(S) list(R1 = methodes(S, 1), R2 = methodes(S, 2)))

## ------------------------------------------------- tables of saliences ----
tab_saliences <- function(r) {
  do.call(rbind, lapply(names(jeux), function(j) {
    f <- fits[[j]][[paste0("R", r)]]
    A <- do.call(cbind, lapply(names(f), function(m) {
      a <- f[[m]]$saliences
      colnames(a) <- if (r == 1) m else paste0(m, seq_len(r))
      a
    }))
    data.frame(jeu = j, bloc = names(jeux[[j]]), PI = f$CPCA$plafonds, A,
               row.names = NULL, check.names = FALSE)
  }))
}
ecrire_csv(tab_saliences(1), sortie("saliences_R1.csv"), row.names = FALSE)
ecrire_csv(tab_saliences(2), sortie("saliences_R2.csv"), row.names = FALSE)

nc <- do.call(rbind, lapply(names(jeux), function(j) do.call(rbind, lapply(1:2, function(r) {
  f <- fits[[j]][[paste0("R", r)]]
  do.call(rbind, lapply(names(f), function(m)
    data.frame(jeu = j, R = r, methode = m,
               niveau2 = mean(f[[m]]$resume$niveau2),
               contraste = mean(f[[m]]$resume$contraste),
               rho_indort = fits[[j]]$R1$INDORT$resume$rho)))
}))))
ecrire_csv(nc, sortie("niveau_contraste.csv"), row.names = FALSE)

## --------------------------------------------------------- heterogeneity ----
het <- do.call(rbind, lapply(names(jeux), function(j) do.call(rbind, lapply(1:2, function(r) {
  f <- fits[[j]][[paste0("R", r)]]
  do.call(rbind, lapply(names(f), function(m) {
    h <- heterogeneite(jeux[[j]], f[[m]])
    data.frame(jeu = j, R = r, methode = m, H = h$H, H_relative = h$H_relative,
               indice_RV = h$indice_RV, part_restituee = sum(h$part_restituee))
  }))
}))))
ecrire_csv(het, sortie("heterogeneite.csv"), row.names = FALSE)

diag_vosa <- do.call(rbind, lapply(names(jeux), function(j)
  cbind(jeu = j, diagnostic_blocs(fits[[j]]$R1$VOSA, jeux[[j]]))))
ecrire_csv(diag_vosa, sortie("diagnostic_blocs_VOSA.csv"), row.names = FALSE)

## ------------------------------------------- path and relative bounds ----
bornes <- do.call(rbind, lapply(names(jeux), function(j) {
  ch <- chemin_spa(jeux[[j]], alphas = seq(0.5, 1, by = 0.01), nstart = 10, seed = graine)
  cbind(jeu = j, borne_relative(ch))
}))
ecrire_csv(bornes, sortie("borne_relative.csv"), row.names = FALSE)

## ------------------------------------- CRC: paired t tests (47 pairs) ----
n_p <- 47
stopifnot(all(meta_crc$paire[1:n_p] == meta_crc$paire[n_p + 1:n_p]))
test_crc <- function(S, f, libelle) {
  do.call(rbind, lapply(names(f), function(m) {
    s <- sqrt(n_p) * f[[m]]$U[, 1]
    tt <- stats::t.test(s[n_p + 1:n_p], s[1:n_p], paired = TRUE)
    data.frame(blocs = libelle, methode = m, t = unname(tt$statistic),
               ddl = unname(tt$parameter), p = tt$p.value)
  }))
}
crc4 <- test_crc(jeux$CRC, fits$CRC$R1, "4 blocs")
set.seed(graine)
S3 <- jeux$CRC[-1]
f3 <- list(CPCA = cpca(S3), INDORT = indort(S3, nstart = 20),
           VOSA = vosa(S3, nstart = 20), STATIS = statis(S3, normalisation = "trace"))
crc3 <- test_crc(S3, f3, "sans biomarqueurs")
ecrire_csv(rbind(crc4, crc3), sortie("crc_test_apparie.csv"), row.names = FALSE)

## ----------------------------------------- canonical factorization ----
fc <- lapply(fits, function(fj) lapply(fj, function(fr) lapply(fr, factorisation_canonique)))
saveRDS(fc, sortie("factorisation_canonique.rds"))
for (j in names(fc)) for (r in 1:2) {
  fr <- fc[[j]][[paste0("R", r)]]
  D <- do.call(cbind, lapply(names(fr), function(m) {
    d <- fr[[m]]$D
    colnames(d) <- paste0(m, "_", colnames(d))
    d
  }))
  ecrire_csv(data.frame(bloc = rownames(D), D, check.names = FALSE),
                   sortie(sprintf("factorisation_D_%s_R%d.csv", j, r)), row.names = FALSE)
}
for (j in names(fc)) {
  V <- fc[[j]]$R1$VOSA$V
  ecrire_csv(do.call(rbind, lapply(names(V), function(b)
    data.frame(bloc = b, variable = seq_len(nrow(V[[b]])), chargement = V[[b]][, 1]))),
    sortie(sprintf("chargements_VOSA_%s.csv", j)), row.names = FALSE)
}

## -------------------------------------------------------------- figures ----
grDevices::pdf(sortie("niveau_contraste.pdf"), width = 10, height = 6.5)
graphics::par(mfrow = c(2, 3))
for (r in 1:2) for (j in names(jeux))
  plot_niveau_contraste(fits[[j]][[paste0("R", r)]],
                        main = sprintf("%s, R = %d (rho = %.2f)", j, r, fits[[j]]$R1$INDORT$resume$rho))
grDevices::dev.off()

grDevices::pdf(sortie("borne_relative.pdf"), width = 10, height = 3.5)
graphics::par(mfrow = c(1, 3))
for (j in names(jeux)) plot_borne_relative(bornes[bornes$jeu == j, ], main = j)
grDevices::dev.off()

grDevices::pdf(sortie("crc_indort_vosa.pdf"), width = 8, height = 4)
graphics::par(mfrow = c(1, 2))
statut <- factor(meta_crc$statut)
for (m in c("INDORT", "VOSA")) {
  s <- sqrt(n_p) * fits$CRC$R1[[m]]$U[, 1]
  graphics::boxplot(s ~ statut, main = m, ylab = "score x sqrt(n)")
  graphics::stripchart(s ~ statut, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.6)
}
grDevices::dev.off()

grDevices::pdf(sortie("saliences_VOSA.pdf"), width = 10, height = 4)
graphics::par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))
for (j in names(jeux)) plot(fits[[j]]$R1$VOSA, main = j)
grDevices::dev.off()


## figures of the canonical factorization (VOSA, R = 2): D_k, U and V_k
grDevices::pdf(sortie("factorisation_VOSA.pdf"), width = 10, height = 10)
graphics::par(mfrow = c(3, 3), mar = c(7, 4, 3, 1))
for (j in names(jeux)) {
  f <- fc[[j]]$R2$VOSA
  plot_D(f, main = paste(j, "- D_k (saliences)"))
  plot_U(f, groupes = if (j == "CRC") meta_crc$statut else NULL, main = paste(j, "- U"))
  b <- which.max(f$D[, 1])                       # block that carries the first VOSA dimension
  plot_V(f, bloc = b, dims = 1)
}
grDevices::dev.off()

## ------------------------------------------------ simulation (optional) ----
if (refaire_simulation) {
  ## Question 1: precision of the first dimension of g_alpha as a function of alpha
  scenarios <- list(S1 = c(0.15, 0.15, 0, 0), S2 = c(0.30, 0.30, 0, 0), S3 = c(0.30, 0, 0, 0))
  alphas <- seq(0, 1, by = 0.05)
  set.seed(graine)
  sim1 <- do.call(rbind, lapply(names(scenarios), function(sc) do.call(rbind, lapply(1:40, function(rep) {
    s <- simuler_blocs(pi_d = scenarios[[sc]])
    ch <- chemin_spa(blocs_vers_S(s$X), alphas = alphas, nstart = 8)
    data.frame(scenario = sc, rep = rep, alpha = alphas,
               cos2 = vapply(ch$solutions, function(f) sum(f$U[, 1] * s$t_d)^2, 0))
  }))))
  ecrire_csv(sim1, sortie("simulation_alpha.csv"), row.names = FALSE)
  ## Question 2: precision as a function of the heterogeneity H, R = 1 to 4
  sim2 <- do.call(rbind, lapply(c(1, 2), function(m) do.call(rbind, lapply(c(0, .05, .1, .15, .2, .3), function(pd)
    do.call(rbind, lapply(1:20, function(rep) {
      s <- simuler_blocs(pi_d = c(rep(pd, m), rep(0, 4 - m)))
      S <- blocs_vers_S(s$X)
      do.call(rbind, lapply(1:4, function(R) {
        U <- list(CPCA = cpca(S, R)$U, STATIS = statis(S, R)$U,
                  INDORT = indort(S, R, nstart = 8)$U, VOSA = vosa(S, R, nstart = 8)$U)
        data.frame(m = m, pi_d = pd, rep = rep, H = heterogeneite(S)$H, R = R,
                   methode = names(U), precision = vapply(U, function(u) sum(crossprod(u, s$t_d)^2), 0))
      }))
    }))))))
  ecrire_csv(sim2, sortie("simulation_heterogeneite.csv"), row.names = FALSE)
}

writeLines(utils::capture.output(utils::sessionInfo()), sortie("sessionInfo.txt"))
message("Results written in ", normalizePath(dossier_resultats))
