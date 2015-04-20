\name{siESR1}
\docType{data}
\alias{siESR1}
\alias{siESR1limma}

\title{Knockdown dataset: siESR1 in MCF-7 cells.}

\description{
  The data consists of 6 microarray samples after knocking down ESR1 in MCF-7 cells,
  retrived from GEO dataset (Series GSE18431). The data have been normalized (quantile 
  normalization and log2 transformatoin) as described in the orifinal dataset, and are 
  presented here in the form of an \code{exprSet} object.
}

\usage{
data(siESR1)
}
\format{
  The format is:
  An \code{ExpressionSet} object with covariates representing experimental conditions:
\itemize{
\item \code{Sample}: Sample IDs
\item \code{Target}: Group names for differential expression analysis (e.g. targets for limma). The 
corresponding constrasts can be retrived by \code{notes(siESR1)}
}
}

\source{
Park YY, Kim K, Kim SB, Hennessy BT et al. Reconstruction of nuclear receptor network reveals that NR2E3 is a novel upstream regulator of ESR1 in breast cancer. EMBO Mol Med 2012 Jan;4(1):52-67.
}

\note{
The differential expression analysis documented in the vignette is available at 'siESR1limma'.
}

\examples{
  data(siESR1)
  #notes(siESR1)
  
  #limma pre-processed dataset
  #data(siESR1limma)
}

\keyword{datasets}

