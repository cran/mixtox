\name{gcaPred}
\alias{gcaPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture Toxicity Prediction Using GCA (General)}
\description{Predict the mixture toxicity based on individual 
			concentration-response information. Thirteen monotonic(sigmoidal) models
			("Hill", "Hill_two", "Hill_three", "Hill_four",
		"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
		"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") are 
		incorporated to construct the GCA equation. The dichotomy 
			technique is used to solve the constructed equation. 
			Three optional mixture design methods are provided. One is the arbitrary 
			concentration ratio (acr) for mixture components. Users can arbitrarily 
			deign a random ratio for each component in the 
			mixture. Other two options are equal effect concentration 
			ratio (eecr) and uniform design concentration ratio (udcr).
			}
\usage{gcaPred(model, param, mixType = c("acr", "eecr", "udcr"), effv, refEffv = 
				c(0.1, 0.50, 0.80), lb = 1E-7, ub = 0.9)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{vector of equation names}
  \item{param}{numeric matrix of fitting coefficients with rownames (equation selected) and 
  colnames (ALpha, Beta, and Gamma). For equations with two parameters, Gamma can be set as 
  zero or any other numeric value.
  }
  \item{mixType}{experimental design of the mixture. acr: arbitrary concentration ratio; 
		eecr: equal effect concentration ratio; udcr: uniform design concentration ratio.
		}
  \item{effv}{numeric vector with single or multiple effect values (0 ~ 1).}
  \item{refEffv}{to determine the concentration ranges for predicting effect.}
  \item{lb}{lower bound for solving constructed IA using diochotomy method.}
  \item{ub}{upper bound bound for solving constructed IA using diochotomy method.}
}
\details{The generalized concentration addition (GCA) model was proposed by Howard and Webster to 
 predict mixtures containing partial agonists (Howard and Webster, 2009). 
 \deqn{\sum\limits_{i = 1}^n {\frac{{{c_i}}}{{f_i^{ - 1}\left( E \right)}}}  = 1}
 Empirical data were used to fit concentration-response function, and then predict the mixture 
 response using the inverse function. Previous studies used Hill_two function to fit individual 
 concentration response curves in the GCA prediction (Hadrup et al., 2013; Howard et al., 2010). 
 Here, we incorporated thirteen functions to construct
 the GCA equation and the dichotomy technique is used to solve the constructed equation. 
 }
\value{
\item{x}{ a series of concentrations}
\item{e}{a series of effects caused by the concentrations (x) as predicted by gca}
\item{pct}{the concentration ratio (percent) of every component in the mixture}
\item{uniTab}{the uniform design table used to construct the mixture when mixType is udcr}
}
\references{
Howard, G.J., Schlezinger, J.J., Hahn, M.E., Webster, T.F., 2010. Generalized Concentration 
Addition Predicts Joint Effects of Aryl Hydrocarbon Receptor Agonists with Partial Agonists 
and Competitive Antagonists. Environ. Health Perspect. 118, 666-672.\cr
Howard, G.J., Webster, T.F., 2009. Generalized concentration addition: A method for examining 
mixtures containing partial agonists. J. Theor. Biol. 259, 469-477.
}
\author{xiangwei zhu}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{gcaHill}}
}
\examples{
## example 1
# using GCA to predict the mixtures designed by equal effect concentration ratio (eecr) at the 
# effect concentration of EC05 and EC50
# the eecr mixture design is based on seven antibiotics(seven factors).
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
gcaPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))

## example 2
# using GCA to predict the mixtures designed by uniform design concentration ratio (udcr)
# the udcr mixture design is based on seven antibiotics(seven factors).
# five levels (EC05, EC10, EC20, EC30, and EC50 ) are allocated in the uniform table using the 
# pseudo-level technique (Liang et al., 2001)
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
effv <- c(0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.30, 0.30, 0.50, 0.50)
gcaPred(model, param, mixType = "udcr", effv)

## example 3
# using GCA to predict the mixtures designed by arbitrary concentration ratio (acr)
# the udcr mixture design is based on four heavy metals (four factors).
# the every component in the mixture shares exactly the same ratio (0.25) 
model <- cytotox$sgl$model[1 : 4]
param <- cytotox$sgl$param[1 : 4, ]
effv <- c(0.25, 0.25, 0.25, 0.25)
gcaPred(model, param, mixType = "acr", effv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized concentration addition }
\keyword{ equal effect concentration ratio }
\keyword{ uniform design concentration ratio }
\keyword{ arbitrary concentration ratio }
\keyword{ uniform design table}
\keyword{ mixture effect}
