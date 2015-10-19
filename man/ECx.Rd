\name{ECx}
\alias{ECx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect Concentration Calculation}
\description{
Effect concentrations are calculated at particular effects based on the fitting coefficients 
of six equations (Hill, Weibull, Logit, BCW, BCL, GL)
}
\usage{ECx(model, param, effv)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{character vector of equations: Hill, Weibull, Logit, BCW, BCL, GL}
  \item{param}{numeric matrix of fitting coefficients with rownames (equation selected) and 
		colnames (ALpha, Beta, and Gamma).For equations with two parameters, Gamma can be set 
		as zero or any other numeric value.
		}
  \item{effv}{numeric vector with single or multiple effect values (0 ~ 1)}
}

\details{
effect concentrations will be calculated with provided equations (model), associated fitting 
parameters (param), and effects (effv). Effect (effv) should be a value(s) between 0~1. For 
example, \eqn{effv} should be 0.5 if we want to calculate a concentration causes 50\% effect.
}

\value{ a numeric vector of effect concentration(s)}
\references{
Hill equation (biochemistry) http://en.wikipedia.org/wiki/Hill_equation_(biochemistry)\cr
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the 
Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry 20(2):448-457.
}
\author{xiangwei zhu}
\note{More equations are needed if we want to calculate effect concentrations based on 
	concentration-response curves fitted by other equations
	}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{CEx}}
}
\examples{
# data(cytotox)
## example 1
# calculate EC50 of heavy metal Ni(2+) on the MCF-7 cells
model.ni <- cytotox$sgl$model
param.ni <- cytotox$sgl$param
ECx(model.ni, param.ni, effv = 0.5) 

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
