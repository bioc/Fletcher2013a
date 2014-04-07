\name{Exp2}
\alias{Exp2}
\alias{Exp2limma}
\docType{data}

\title{Exp2 dataset: iF2 construct experiments.}

\description{
 The data consists of 71 microarray samples from MCF-7 cells treated under 
 different conditions, at 3 time points (0, 6 and 24 h). The data have been 
 normalized (RMA algorithm) and are presented in the form of an \code{exprSet} 
 object. The experiment was carried out on 6 Humanv4 BeadChips using 12 samples 
 per BeadChip. The original arrays contain 48324 features, with a mean of 22 beads 
 per feature (Standard Deviation of 5). Additional description about sample 
 groups can be retrieved from \code{phenoData} slot (see examples).
}

\usage{
data(Exp2)
}

\format{
  The format is:
  An \code{ExpressionSet} object with covariates representing experimental conditions:
\itemize{
\item \code{Sample}: Sample IDs.
\item \code{Time}: Treatment time.
\item \code{Replicates}: Biological replicates (numerical sequence).
\item \code{UT}: Control group (vehicle).
\item \code{E2}: E2 treatment group.
\item \code{E2.AP20187}: E2+AP20187 treatment group.
\item \code{E2.AP20187.PD}: E2+AP20187+PD173074 treatment group.
\item \code{E2.FGF10}: E2+FGF10 treatment group.
\item \code{Treatment}: Group names.
\item \code{TreatmentGroups}: Group names with common starting point.
\item \code{Target}: Group names for differential expression analysis (e.g. targets for limma). The 
corresponding constrasts can be retrived by \code{notes(Exp2)}.
\item \code{TecRep}: Technical replicates.
\item \code{isOriginal}: Simple vector mapping all samples (excluding technical replicates).
}
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.
}

\note{
The differential expression analysis documented in the vignette is available at 'Exp2limma'.
}

\examples{
  data(Exp2)
  #gexp<-exprs(Exp2)
  #geneids<-fData(Exp2)
  #targets<-pData(Exp2)
  #mycontrasts<-notes(Exp2)$contrasts
  
  #limma pre-processed dataset
  #data(Exp2limma)
}
\keyword{datasets}

