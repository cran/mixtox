\name{CEx}
\alias{CEx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect Calculation}
\description{
Effects can be calculated at particular concentrations based on the inverse function of Hill, 
Weibull, Logit, BCW, BCL, GL.
}
\usage{
CEx(model, param, conc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{character vector of equation names: Hill, Weibull, Logit, BCW, BCL, GL}
  \item{param}{numeric matrix of fitting coefficients with rownames (equation selected) and 
  colnames (ALpha, Beta, and Gamma). For equations with two parameters, Gamma can be set as 
  zero or any other numeric value.
  }
  \item{conc}{numeric vector with single or multiple concentrations.}
}
\details{effects will be calculated with provided equations (model), associated fitting parameters
		(param), and concentrations(conc). The inverse function of those six sigmoidal equations
		are as follows:\cr		
		inverse Hill: \deqn{{c = \beta E/\left( {\alpha  - E} \right)}}
		inverse Weibull: \deqn{c = {10^{\left( {\ln ( - \ln (1 - E)) - \alpha } \right)/\beta }}}
		inverse Logit: \deqn{c = {10^{\left( {\ln (E/(1 - E)) - \alpha } \right)/\beta }}}
		inverse BCW: \deqn{c = {\left( {(\gamma /\beta )(\ln ( - \ln (1 - E)) - \alpha )
						+ 1} \right)^{1/\gamma }}}
		inverse BCL: \deqn{c = {((\gamma /\beta )( - \ln ((1 - E)/E) - \alpha ) + 1)^{1/\gamma }}}
		inverse GL: \deqn{c = {10^{(( - \ln ({{(1/E)}^{1/\gamma }} - 1) - \alpha )/\beta )}}}
		where \eqn{E} is effect and \eqn{c} is the concentration.	
		}
\value{ a numeric vector of effect(s)}
\references{
Hill equation (biochemistry) http://en.wikipedia.org/wiki/Hill_equation_(biochemistry)\cr
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the
 Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry 20(2):448-457.
}
\author{xiangwei zhu}
\note{More equations are needed if we want to calculate effect concentrations based on 
concentration-response curves fitted by other equations
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ECx}}
}
\examples{
# data(cytotox)

## example 1
# calculate the effect of heavy metal Ni(2+) at the concentration of 0.00005 mol/L on the MCF-7
model.ni <- cytotox$sgl$model
param.ni <- cytotox$sgl$param
CEx(model.ni, param.ni, conc = 0.00005)

## example 2
# calculate the effect caused by four heavy metals and four ionic liquids at the concentration of
# 0.00001 and 0.00002 mol/L on the MCF-7 cells
model <- cytotox$sgl$model
param <- cytotox$sgl$param
CEx(model, param, conc = c(0.00001, 0.00002)) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ effect concentration }
\keyword{inverse function}
