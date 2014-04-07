\name{Exp3}
\alias{Exp3}
\alias{Exp3limma}
\docType{data}

\title{Exp3 dataset: FGFR2b over-expression experiments.}

\description{
 The data consists of 125 microarray samples from MCF-7 cells treated under 
 different conditions, at 5 time points (0, 3, 6, 12 and 24 h). The data have 
 been normalized (RMA algorithm) and are presented in the form of an \code{exprSet} 
 object. The experiment was carried out on 11 Humanv4 BeadChips using 12 samples
 per BeadChip. The original arrays contains 48324 features, with a mean of 22 beads 
 per feature (Standard Deviation of 5). Additional description about sample 
 groups can be retrieved from \code{phenoData} slot (see examples).
}

\usage{
data(Exp3)
}

\format{
  The format is:
  An \code{ExpressionSet} object with covariates representing experimental conditions:
\itemize{
\item \code{Sample}: Sample IDs.
\item \code{Group}: Major experimental groups (MinusTet and PlusTet).
\item \code{Time}: Treatment time.
\item \code{Replicates}: Biological replicates (numerical sequence).
\item \code{MinusTet.UT}: Control group (vehicle for MinusTet group).
\item \code{MinusTet.E2}: E2 treatment group.
\item \code{MinusTet.E2.FGF10}: E2+FGF10 treatment group.
\item \code{PlusTet.UT}: Control group (vehicle for PlusTet group).
\item \code{PlusTet.E2}: E2 treatment group.
\item \code{PlusTet.E2.FGF10}: E2+FGF10 treatment group.
\item \code{Treatment}: Group names (treatments).
\item \code{Target}: Group names for differential expression analysis (e.g. targets for limma). The 
corresponding constrasts can be retrived by \code{notes(Exp3)}.
}
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.
}

\note{
The differential expression analysis documented in the vignette is available at 'Exp3limma'.
}

\examples{
  data(Exp3)
  #gexp<-exprs(Exp3)
  #geneids<-fData(Exp3)
  #targets<-pData(Exp3)
  #mycontrasts<-notes(Exp3)$contrasts
  
  #limma pre-processed dataset
  #data(Exp3limma)
}

\keyword{datasets}

