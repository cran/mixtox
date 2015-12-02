\name{antibiotox}
\alias{antibiotox}
\docType{data}
\title{Toxicity of Seven Antibiotics on Photobacteria}
\description{
Seven antibiotics include: Neomycin sulfate(NEO), streptomycin sulfate(STR), kanamycin sulfate(KAN),
spectinomycin dihydrochloridehydrate(SPE), paromomycin sulfate(PAR), 
dihydrostreptomycin sesquisulfate hydrate(DIH), and gentamycin sulfate(GEN). Their toxicity
on photobacteria were tested using microtoxicity analysis. \cr
The concentration-response data include the toxicity of seven antibiotics and ten mixtures 
designed by udcr. The fitting information of seven antibiotics and 10 udcr mixtures. 
}
\usage{data(antibiotox)}
\format{
  A list with concentration-response data of 7 antibiotics, 10 udcr mixtures, and associated
  fitting information.
  \describe{
    \item{\code{PAR$x}}{a numeric vector of concentrations}
    \item{\code{PAR$y}}{a numeric matrix of responses}
	\item{\code{PAR$name}}{name of test substance}
	\item{\code{PAR$concNum}}{the number of concentrations}
	\item{\code{PAR$tierNum}}{the number of repetitions}
	\item{\code{PAR$type}}{type of test substance}
	
	\item{\code{sgl$model}}{model used to fit the concentration-response data
		}
	\item{\code{sgl$param}}{fitted coefficients of concentration-response curves
		}
	
	\item{\code{udcr.mix$model}}{model used to fit the concentration-response data of udcr mixtures on 
		photobacteria
		}
	\item{\code{udcr.mix$param}}{fitted coefficients of concentration-response curves
		}
	\item{\code{udcr.pct}}{the percentage of individual chemicals in the udcr mixtures}
	
	\item{\code{eecr.mix$model}}{model used to fit the concentration-response data of eecr mixtures on 
		photobacteria
		}
	\item{\code{eecr.mix$param}}{fitted coefficients of concentration-response curves
		}
	\item{\code{eecr.pct}}{the percentage of individual chemicals in the eecr mixtures}
	
  }
}
\details{
The monotonic concentration-response data need to be scaled into [0, 1] to suit some equaitons(
e.g., Weibull and Logit) using the following equation:
\deqn{E = \frac{{{I_0} - {I_i}}}{{{I_0}}}}
where \eqn{I_0} is the average of controls for inhibition test or the average 
of the maximum effect for the stimulation test and \eqn{I_i} the average 
effect of the \eqn{i^{th}} treatment.
}
\examples{
# example 1
## Retrive the toxicity data of PAR on photobacteria.
PAR <- antibiotox$PAR
}
\keyword{concentration-response data}
\keyword{curve fitting}