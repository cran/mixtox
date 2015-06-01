\name{cytotox}
\alias{cytotox}
\docType{data}
\title{Cytotoxicity of Heavy Metals and Ionic Liquids on MCF-7}
\description{
Chemicals include four heavy metals: 
\eqn{{\rm{NiN}}{{\rm{O}}_{\rm{3}}}\cdot 6{{\rm{H}}_{\rm{2}}}{\rm{O}}}(Ni),
\eqn{{\rm{ZnS}}{{\rm{O}}_{\rm{4}}} \cdot {\rm{7}}{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Zn), 
\eqn{{\rm{CuS}}{{\rm{O}}_{\rm{4}}} \cdot 5{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Cu), and 
\eqn{{\rm{MnC}}{{\rm{l}}_{\rm{2}}} \cdot 4{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Mn); 
four ionic liquids are 1-Octyl-3-methylimidazolium chloride (Omim), 
1-Dodecyl-3-methylimidazolium chloride (Dmim), 1-Ethyl-3-methylimidazolium tetrafluoroborate 
(Emim), and 1-Hexyl-3-Methylimidazolium tetrafluoroborate (Hmim). \cr
The concentration-response data include the cytotoxicity of these eight compounds 
and ten mixtures designed by udcr.\cr
The fitting information of eight chemicals and 10 uncr mixtures. 
}
\usage{data(cytotox)}
\format{
  A list with concentration-response data of 8 chemicals and 10 uncr mixtures and their associated
  fitting information.
  \describe{
    \item{\code{Ni$x}}{a numeric vector of concentrations.}
    \item{\code{Ni$y}}{a numeric matrix of response values, that is relative change in volume.}
	\item{\code{Ni$name}}{name of test substance}
	\item{\code{Ni$concNum}}{the number of concentrations}
	\item{\code{Ni$tierNum}}{the number of repetitions}
	\item{\code{Ni$type}}{type of test substance}
	
	\item{\code{sgl$model}}{model used to fit the concentration-response data of individual 
		chemicals on MCF-7
		}
	\item{\code{sgl$param}}{fitted coefficients of the concentration-response curves corresponding 
		to sgl$model
		}
	
	\item{\code{mix$model}}{model used to fit the concentration-response data of udcr mixtures on 
		MCF-7
		}
	\item{\code{mix$param}}{fitted coefficients of the concentration-response curves corresponding
		to mix$model
		}
	\item{\code{$pct}}{the percent of individual chemicals in the uncr mixtures }	
  }
}
\details{
The monotonic concentration-response data need to be scaled into [0, 1] using the 
following equation:
\deqn{E = \frac{{{I_0} - {I_i}}}{{{I_0}}}}
where \eqn{I_0} is the average of controls for inhibition test or the average 
of the maximum effect for the stimulation test and \eqn{I_i} the average 
effect of the \eqn{i^{th}} treatment.
}
\source{ The cytotoxicity experiments were conducted in our lab}
\examples{
data(cytotox)
Ni <- cytotox$Ni
}
\keyword{concentration-response data}
\keyword{curve fitting}