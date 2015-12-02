\name{caPred}
\alias{caPred}
\title{Mixture Toxicity Prediction Using Concentration Addition}
\description{Predict the mixture toxicity based on individual concentration-response 
			information. Three optional mixture design methods are provided. One 
			is the arbitrary concentration ratio (acr) for mixture components. Users can 
			deign a random ratio for each
			component in the mixture. Other two options are equal effect 
			concentration ratio (eecr) and uniform design concentration ratio (udcr).
}
\usage{caPred(model, param, mixType = c("acr", "eecr", "udcr"), effv)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{character vector of equation names: Hill, Weibull, Logit, BCW, BCL, GL}
  \item{param}{numeric matrix of fitting coefficients with row names (equation selected) and 
		column names (Alpha, Beta, and Gamma). For equations with only two parameters, Gamma can 
		be set as zero or any other numeric value.}
  \item{mixType}{experimental design of the mixture. 
			acr: arbitrary concentration ratio; 
			eecr: equal effect concentration ratio; 
			udcr: uniform design concentration ratio.}
  \item{effv}{numeric vector with single or multiple effect values (0 ~ 1).}
}
\details{Concentration addition (CA) is designed for mixtures of chemicals that have similar 
		mechanisms of action. For a well-defined mixture (e.g., a mixture of n components), 
		CA is expressed mathematically as: 
		\deqn{\sum\limits_{i = 1}^n {\frac{{{c_i}}}{{EC{x_i}}}}  = 1}
		where \eqn{EC{x_i}} is the effect concentration of the \eqn{i^{th}} compound that 
		causes x\% effect when applied individually at \eqn{{c_i}}. The \eqn{{c_i}} can be 
		computed from the following equation:
		\deqn{{c_i} = {p_i} \cdot {c_{mix}} = {p_i} \cdot E{C_{x,mix}}}
		where \eqn{p_i} is the proportion of \eqn{i^{th}} component in the mixture, 
		\eqn{c_{mix}} the mixture concentration and \eqn{E{C_{x,mix}}} the 
		concentration of the mixture at a specific effect x\%. The prediction 
		of combined effects of mixture-components based on 
		CA can then be expressed as:
		\deqn{E{C_{x,mix}} = {\left( {\sum\limits_{i = 1}^n {\frac{{{p_i}}}{{E{C_{x,i}}}}}} 
				\right)^{ - 1}}}
}
\value{
\item{ca}{ a series of effect concentrations predicted by CA}
\item{e}{a series of effects associated with the effect concentrations in ca}
\item{pct}{the concentration ratio (percent) of every component in the mixture}
\item{uniTab}{the uniform design table used to construct the mixture when mixType is 'udcr'}
}
\references{
Liang, Yi-zeng, Kai-tai Fang, and Qing-song Xu. 2001. Uniform Design and Its Applications in 
Chemistry and Chemical Engineering. Chemometrics and Intelligent Laboratory Systems 
58(1):43-57.\cr
Backhaus, T., Faust, M., 2012. Predictive environmental risk assessment of chemical 
mixtures: A conceptual framework. Environmental Science and Technology. 46, 2564-2573.
}
\author{xiangwei zhu}
\note{
Note that effv is dependent on the mixType. 
if the mixType is acr, the length of effv should be the same as that of model. 
if the mixType is eecr, effv should contain at least one value 
if the mixTpe is udcr, elements in effv are levels, the length of effv is the number of runs. 
the number of runs should be in accordance with the length of model based on the uniform design 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{iaPred}}
}
\examples{
## example 1
# using CA to predict the mixtures designed by equal effect concentration ratio (eecr) at the
# effect concentration of EC05 and EC50
# the eecr mixture design is based on seven antibiotics(seven factors).
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
caPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))

## example 2
# using CA to predict the mixtures designed by uniform design concentration ratio (udcr)
# the udcr mixture design is based on four heavy metals and four ionic liquids (eight factors).
# five levels (EC05, EC10, EC20, EC30, and EC50 ) are allocated in the uniform table using the
# pseudo-level technique (Liang et al., 2001)
model <- cytotox$sgl$model
param <- cytotox$sgl$param
effv <- c(0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.30, 0.30, 0.50, 0.50)
caPred(model, param, mixType = "udcr", effv)

## example 3
# using CA to predict the mixtures designed by arbitrary concentration ratio (acr)
# the udcr mixture design is based on five antibiotics (five factors).
# the every component in the mixture shares exactly the same ratio (0.20) 
model <- antibiotox$sgl$model[1 : 5]
param <- antibiotox$sgl$param[1 : 5, ]
effv <- c(0.2, 0.2, 0.2, 0.2, 0.2)
caPred(model, param, mixType = "acr", effv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ concentration addition}
\keyword{ equal effect concentration ratio }
\keyword{ uniform design concentration ratio }
\keyword{ arbitrary concentration ratio }
\keyword{ uniform design table}
\keyword{ mixture effects}
