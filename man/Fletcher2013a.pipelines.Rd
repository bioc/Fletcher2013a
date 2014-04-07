\name{Fletcher2013a.pipelines}
\alias{Fletcher2013pipeline.limma}
\alias{Fletcher2013pipeline.pca}
\alias{Fletcher2013pipeline.deg}
\alias{Fletcher2013pipeline.supp}

\title{
A pipeline to reproduce results for Fletcher et al. 2013.
}

\description{
Pipeline functions to reproduce results in Fletcher et al. 2013.
}

\usage{
 Fletcher2013pipeline.limma(exprSet, par=list(method="global", adjust.method="BH", p.value=1e-2, lfc=0))
 Fletcher2013pipeline.pca(exprSet)
 Fletcher2013pipeline.deg(what="Exp1", idtype="probeid", response="all", mode="all")
 Fletcher2013pipeline.supp()
}

\arguments{
\item{exprSet}{
Any \code{ExpressionSet} object available in the package \code{Fletcher2013a}: 
\code{\link[Fletcher2013a:Exp1]{Exp1}}, \code{\link[Fletcher2013a:Exp2]{Exp2}} and 
\code{\link[Fletcher2013a:Exp3]{Exp3}}.
}

\item{par}{
a list with parameters for limma analysis (see \code{\link[limma:decideTests]{decideTests}}) for details.
}

\item{what}{
a single character value specifying an experimental system: "Exp1", "Exp2" or "Exp3".  
}

\item{idtype}{
a single character value specifying a signature report option: "probeid" or "entrez".
}

\item{response}{
a single character value specifying how signatures should be consolidated across time points: "early", "late" or "all". (see \code{Details} for a detailed time-course dataset).
}

\item{mode}{
a single character value specifying how signatures should be consolidated regarding the mode of action: "positive", "negative" or "all". (see \code{Details} for a detailed time-course dataset).
}

}

\details{

These functions reproduce results of the differential expression analysis in Fletcher et al., 2013:

(1) \code{Fletcher2013pipeline.limma}: main function to run the limma analysis on the gene expression datasets.

(2) \code{Fletcher2013pipeline.pca}: complementary function to run a principal components analysis on the gene expression datasets.

(3) \code{Fletcher2013pipeline.deg}: function to extract retults computed in \code{Fletcher2013pipeline.limma}. 
Usefull to retrieve consolidated gene lists from the differential expression analyses. For detailed 
time-course information, please use the full pre-processed datasets 
\code{\link[Fletcher2013a:Exp1limma]{Exp1limma}}, 
\code{\link[Fletcher2013a:Exp2limma]{Exp2limma}} and 
\code{\link[Fletcher2013a:Exp3limma]{Exp3limma}}.

(4) \code{Fletcher2013pipeline.supp}: this function generates additional figures for the vignette.

}

\value{ 
All results will be saved in the current work directory.
}

\source{
Michael NC Fletcher, Mauro AA Castro, Suet-Feung Chin, Oscar Rueda, Xin Wang, Carlos Caldas, Bruce AJ Ponder, Florian Markowetz, Kerstin B Meyer. Master regulators of FGFR2 signalling and breast cancer risk. Nature Communications, 4:2464, 2013.
}

\author{
Mauro Castro \email{mauro.a.castro@gmail.com}
}

\examples{
\dontrun{
data(Exp1)
#Fletcher2013pipeline.limma(Exp1)
#Fletcher2013pipeline.pca(Exp1)
}
}












