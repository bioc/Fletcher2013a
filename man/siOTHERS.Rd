\name{siOTHERS}
\docType{data}
\alias{siOTHERS}

\title{Knockdown dataset: siPTTG1, siSPDEF, siE2F2 and siELF3 in MCF-7 cells.}

\description{
  The data consists of 30 microarray samples after knocking down PTTG1, SPDEF, 
  E2F2 or ELF3 in MCF-7 cells. The data have been normalized (RMA algorithm) 
  and are presented in the form of an \code{exprSet} object. The experiment 
  was carried out on 3 Humanv4 BeadChips using 30 arrays (1 samples per BeadChip). 
  These arrays interrogate 48107 randomly-distributed bead-types, and in this experiment
  there was a mean of 22 beads per bead-type (Standard Deviation of 5). In total, there
  are 21642 genes being interrogated, with 7872 genes being interrogate by more than 
  one bead-type and 13167 bead-types not being assigned to a gene symbol.
 }

\usage{
data(siOTHERS)
}

\format{
  The format is:
  An \code{ExpressionSet} object with covariates representing experimental conditions:
\itemize{
\item \code{Sample}: Sample IDs.
\item \code{Target}: Group names for differential expression analysis (e.g. targets for limma). The 
corresponding constrasts can be retrived by \code{notes(siOTHERS)}.
}
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Journal paper (submitted). 2012.
}

\examples{
  data(siOTHERS)
  #notes(siOTHERS)
}

\keyword{datasets}

