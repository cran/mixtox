\name{nmECx}
\alias{nmECx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect Concentration Calculation for Non-monotonic Curves}
\description{
Effect concentrations are calculated at particular effects based on the fitting coefficients 
of 13 monotonic equations.
}
\usage{nmECx(model, param, effv, minx)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{character vector of equations:("Brain_Consens", "BCV", "Beckon", "Biphasic",
               "Hill_six")}
  \item{param}{numeric matrix of fitting coefficients with rownames (equation selected) and 
		colnames (ALpha, Beta, Gamma, Delta, and Epsilon).
		}
  \item{effv}{numeric vector with single or multiple effect values (miny ~ 1)}
  \item{minx}{numeric vector with single or multiple concentrations that 
               induce maximun stimulation}
}

\details{
effect concentrations will be calculated with provided equations(model), associated fitting 
parameters (param), and effects (effv). Effect (effv) should be a value(s) between miny ~ 1. For 
example, \eqn{effv} should be 0.5 if we want to calculate a concentration causes 50\% effect.
\eqn{minx} should be calculated using curveFit
}

\value{ a numeric vector of effect concentration(s)}
\references{
Zhu X-W, Liu S-S, Qin L-T, Chen F, Liu H-L. 2013. Modeling non-monotonic dose-response 
relationships: Model evaluation and hormetic quantities exploration. Ecotoxicology and 
Environmental Safety 89:130-136.\cr
Di Veroli GY, Fornari C, Goldlust I, Mills G, Koh SB, Bramhall JL, et al. 2015. An automated 
fitting procedure and software for dose-response curves with multiphasic features. 
Scitific Report 5: 14701.
Reference to {curveFit}
}
\author{xiangwei zhu}
\note{For equations with less parameters, the default can be set as zero or any other 
 numeric values. 
	}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{CEx}}
\code{\link{curveFit}}
}
\examples{
## example 1
# calculate ECL-10, ECR-10, EC5, and EC50 of the four hormesis curves 
model <- hormesis$sgl$model
param <- hormesis$sgl$param
minx <- hormesis$sgl$minx
nmECx(model, param, effv = c(-0.10, 0.05, 0.50), minx)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{effect concentration}
\keyword{ non-monotonic concentration-response curve}