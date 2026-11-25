#' Language of the results (French or English)
#'
#' The package can display its results in French or in English: headings of
#' the printed results, column names of the printed tables, and labels of the
#' plots. \code{vosa_langue()} returns the current language and, with an
#' argument, sets it for the session (option \code{vosa.langue}; French by
#' default). \code{traduire()} returns a copy of any result of the package
#' (fitted model, factorization, path, table) whose element and column names
#' are translated, for instance before exporting it to a file.
#'
#' The names of the elements of the objects returned by the functions of the
#' package are not changed by the language (they remain those used in the code
#' of the article, e.g. \code{fit$saliences}, \code{fit$resume$contraste}), so
#' that scripts work in both languages; only the display, the plots and the
#' copies made by \code{traduire()} are translated.
#'
#' @param langue \code{"fr"} (French) or \code{"en"} (English); for
#'   \code{vosa_langue()}, \code{NULL} returns the current language.
#' @param x a result of the package: object of class \code{"vosa_fit"},
#'   \code{"vosa_factorisation"} or \code{"vosa_chemin"}, a data frame or a list.
#' @return \code{vosa_langue()} returns the current language (invisibly when it
#'   is set); \code{traduire()} returns a plain list or a data frame with
#'   translated names.
#' @examples
#' sim <- simuler_blocs(seed = 1)
#' fit <- vosa(blocs_vers_S(sim$X), nstart = 3, seed = 1)
#' ancienne <- vosa_langue("en")
#' fit                                   # printed in English
#' names(traduire(fit$resume))           # level, level2, contrast, ...
#' vosa_langue("fr")
#' fit                                   # printed in French
#' vosa_langue(ancienne)
#' @name langue
NULL

#' @rdname langue
#' @export
vosa_langue <- function(langue = NULL) {
  if (is.null(langue)) return(.langue())
  langue <- match.arg(langue, c("fr", "en"))
  options(vosa.langue = langue)
  invisible(langue)
}

#' @rdname langue
#' @export
traduire <- function(x, langue = vosa_langue()) {
  langue <- match.arg(langue, c("fr", "en"))
  if (is.data.frame(x)) {
    names(x) <- .tr_noms(names(x), langue)
    return(x)
  }
  if (is.list(x)) {
    out <- lapply(x, function(el) if (is.list(el)) traduire(el, langue) else el)
    if (!is.null(names(x))) names(out) <- .tr_noms(names(x), langue)
    return(out)
  }
  x
}

## ---- internal ---------------------------------------------------------------

.langue <- function() {
  l <- getOption("vosa.langue", "fr")
  if (!l %in% c("fr", "en")) "fr" else l
}

## Names of elements and columns: French (as in the code) -> English
.dico_noms <- c(
  methode = "method", alpha = "alpha", U = "U", saliences = "saliences",
  resume = "summary", dimension = "dimension", niveau = "level",
  niveau2 = "level2", contraste = "contrast", L2 = "L2", rho = "rho",
  g_alpha = "g_alpha", G_alpha = "G_alpha", valeurs_propres = "eigenvalues",
  critere = "criterion", Wbar = "Wbar", iterations = "iterations",
  converge = "converged", decalage = "shift", criteres_departs = "start_criteria",
  n = "n", K = "K", plafonds = "ceilings", blocs = "blocks",
  poids = "weights", poids_somme1 = "weights_sum1",
  interstructure = "interstructure", RV = "RV", part_RV = "RV_share",
  compromis = "compromise", valeurs_propres_compromis = "compromise_eigenvalues",
  chemin = "path", solutions = "solutions", R = "R",
  D = "D", V = "V", E = "E", part_restituee = "reproduced_share",
  part_residuelle = "residual_share", erreur_reconstruction = "reconstruction_error",
  bloc = "block", plafond = "ceiling", plancher = "floor", salience = "salience",
  fraction_plafond = "ceiling_fraction", cos2_cp1 = "cos2_pc1",
  au_dessus_niveau = "above_level", H = "H", H_relative = "H_relative",
  indice_RV = "RV_index", part_cumulee = "cumulative_share",
  delta_contraste = "contrast_gain", delta_niveau = "level_loss",
  borne_inf = "lower_bound", borne_sup = "upper_bound",
  X = "X", t_d = "t_d", T = "T", jeu = "dataset", variable = "variable",
  chargement = "loading")

.tr_noms <- function(noms, langue = .langue()) {
  if (is.null(noms) || langue == "fr") return(noms)
  out <- unname(.dico_noms[noms])
  ifelse(is.na(out), noms, out)
}

## Messages and labels
.dico_messages <- list(
  saliences      = c(fr = "Saliences :", en = "Saliences:"),
  resume         = c(fr = "Niveau et contraste par dimension :", en = "Level and contrast by dimension:"),
  plancher       = c(fr = "Plancher 1/(n-1) =", en = "Floor 1/(n-1) ="),
  chemin         = c(fr = "Chemin des solutions de g_alpha, R =", en = "Path of solutions of g_alpha, R ="),
  fc_titre       = c(fr = "Factorisation canonique X_k = U D_k V_k' + E_k", en = "Canonical factorization X_k = U D_k V_k' + E_k"),
  fc_D           = c(fr = "Termes diagonaux d_kr des D_k :", en = "Diagonal terms d_kr of D_k:"),
  fc_part        = c(fr = "Part d'inertie restitu\u00e9e (d_kr^2 / ||X_k||^2) :", en = "Share of inertia reproduced (d_kr^2 / ||X_k||^2):"),
  fc_resid       = c(fr = "Part r\u00e9siduelle ||E_k||^2 / ||X_k||^2 :", en = "Residual share ||E_k||^2 / ||X_k||^2:"),
  fc_err         = c(fr = "Erreur de reconstruction :", en = "Reconstruction error:"),
  salience       = c(fr = "salience", en = "salience"),
  plafond        = c(fr = "plafond", en = "ceiling"),
  niveau         = c(fr = "niveau", en = "level"),
  niveau2        = c(fr = "niveau au carr\u00e9", en = "squared level"),
  contraste      = c(fr = "contraste", en = "contrast"),
  dimension      = c(fr = "dimension", en = "dimension"),
  individu       = c(fr = "individu", en = "individual"),
  composante     = c(fr = "composante", en = "component"),
  variable       = c(fr = "variable", en = "variable"),
  chargement     = c(fr = "chargement, dim.", en = "loading, dim."),
  d2             = c(fr = "d_kr^2 (salience)", en = "d_kr^2 (salience)"),
  d              = c(fr = "d_kr", en = "d_kr"),
  score          = c(fr = "score x racine(n)", en = "score x sqrt(n)")
)

.tr <- function(cle, langue = .langue()) {
  m <- .dico_messages[[cle]]
  if (is.null(m)) cle else unname(m[langue])
}
