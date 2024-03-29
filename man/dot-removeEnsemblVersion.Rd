% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stripEnsemblVersions.R
\name{.removeEnsemblVersion}
\alias{.removeEnsemblVersion}
\title{Regex off Ensembl identifier version numbers}
\usage{
.removeEnsemblVersion(x)
}
\arguments{
\item{x}{\code{character()} vector of Ensembl gene identifiers
to remove the version numbers from.}
}
\value{
\code{character()} vector of Ensembl gene identifiers with the version
numbers removed.
}
\description{
Converts Ensembl identifiers from 'ENS*.XX' to 'ENS*', where 'XX' is any
version number and '*' is the rest of the identifier text.
}
\details{
Specific regex pattern removed is '\\..*$'.
}
