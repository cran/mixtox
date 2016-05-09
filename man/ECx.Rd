\name{ECx}
\alias{ECx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect Concentration Calculation for Sigmoidal Models}
\description{
Effect concentrations are calculated at particular effects based on the fitting coefficients 
of 13 Sigmoidal models.
}
\usage{ECx(model, param, effv)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{character vector of equations:("Hill", "Hill_two", "Hill_three", "Hill_four", 
 "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
 "BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)")}
  \item{param}{numeric matrix of fitting coefficients with rownames (equations) and 
		colnames (Alpha, Beta, Gamma, Delta, and Epsilon).
		}
  \item{effv}{numeric vector with single or multiple effect values}
}

\details{effect concentrations will be calculated with provided equations (model), associated 
fitting parameters (param), and effects (effv). 
\eqn{effv} should be 0.5 if we want to calculate a concentration 
causes 50\% effect.\cr 
As an example, the inverse function of the six sigmoidal equations are listed as follows:\cr	
inverse Hill_two: \deqn{{c = \beta E/\left( {\alpha  - E} \right)}}
inverse Weibull: \deqn{c = {10^{\left( {\ln ( - \ln (1 - E)) - \alpha } \right)/\beta }}}
inverse Logit: \deqn{c = {10^{\left( {\ln (E/(1 - E)) - \alpha } \right)/\beta }}}
inverse BCW: \deqn{c = {\left( {(\gamma /\beta )(\ln ( - \ln (1 - E)) - \alpha )
				+ 1} \right)^{1/\gamma }}}
inverse BCL: \deqn{c = {((\gamma /\beta )( - \ln ((1 - E)/E) - \alpha ) + 1)^{1/\gamma }}}
inverse GL: \deqn{c = {10^{(( - \ln ({{(1/E)}^{1/\gamma }} - 1) - \alpha )/\beta )}}}
where \eqn{E} is effect and \eqn{c} is the concentration.
}

\value{ 
\item{ecx}{a numeric vector of effect concentration(s)}
}
\references{
Hill equation (biochemistry) http://en.wikipedia.org/wiki/Hill_equation_(biochemistry)\cr
Reference to {curveFit}
}
\author{xiangwei zhu}
\note{For equations with less parameters, the default can be set as zero or any other 
 numeric values. }

\examples{
## example 1
# calculate EC5 and EC50 of seven antibiotics on the photobacteria
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
effv <- c(0.05, 0.5)
ECx(model, param, effv = c(0.05, 0.50))

## example 2
# calculate EC5 and EC50 of four heavy metals and four ionic liquids on the MCF-7 cells
model <- cytotox$sgl$model
param <- cytotox$sgl$param
ECx(model, param, effv = c(0.05, 0.50))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{effect concentration}
\keyword{ concentration-response curve}
