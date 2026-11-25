% Please edit documentation in R/langue.R
\name{langue}
\alias{langue}
\alias{vosa_langue}
\alias{traduire}
\title{Language of the results (French or English)}
\usage{
vosa_langue(langue = NULL)

traduire(x, langue = vosa_langue())
}
\arguments{
\item{langue}{\code{"fr"} (French, default) or \code{"en"} (English). In
\code{vosa_langue()}, \code{NULL} returns the current language.}

\item{x}{a fit (class \code{"vosa_fit"}), a path (\code{"vosa_chemin"}), a
canonical factorization (\code{"vosa_factorisation"}), a list or a data frame.}
}
\value{
\code{vosa_langue()} returns the current language, or invisibly the language
set. \code{traduire()} returns \code{x} with element names and column names
translated; the values and the class are unchanged.
}
\description{
\code{vosa_langue()} reads or sets the language used by the package. The
language is stored in the option \code{vosa.langue}; it can also be set with
\code{options(vosa.langue = "en")}, for example in a \file{.Rprofile} file.

The language governs the printed output (\code{print}, \code{summary}), the
axis labels and legends of the plots (\code{plot}, \code{plot_U},
\code{plot_D}, \code{plot_V}, \code{plot_niveau_contraste},
\code{plot_borne_relative}) and the column names of the tables shown.

The element names of the objects returned by the functions stay in French
(\code{fit$saliences}, \code{fit$resume}, \code{fc$part_restituee}...), so that
a script works whatever the language. \code{traduire()} returns a copy with
translated names, for export (for example with \code{write.csv}).
}
\examples{
S <- blocs_vers_S(simuler_blocs(n = 30, p = 40, seed = 1)$X)
fit <- vosa(S, R = 1, nstart = 3, seed = 1)
vosa_langue()            # "fr"
fit
vosa_langue("en")
fit
names(traduire(fit))
head(traduire(fit$resume))
vosa_langue("fr")
}
