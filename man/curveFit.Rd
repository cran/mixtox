\name{curveFit}
\alias{curveFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Curve Fitting}
\description{
Six sigmoidal equations (Hill, Logit, Weibull, Box-Cox-Weibull (BCW), Box-Cox-Logit (BCL), 
Generalised Logit (GL)) are provided to fit the concentration-response data. The goodness of fit 
is evaluated by the following statistics: coefficient of 
determination (\eqn{R^2}), adjusted coefficient of determination (\eqn{R_{adj}^2}), 
root mean squared error (RMSE), mean absolute error (MAE), 
and Akaike information criterion (AIC).
}
\usage{
curveFit(x, expr, eq = c("Hill", "Weibull", "Logit", "BCW", "BCL", "GL"), 
param, 	effv, sigLev = 0.05, noec = TRUE, fig = TRUE, algo = default)
}

\arguments{
  \item{x}{a numeric vector of experimental concentrations}
  \item{expr}{a numeric matrix with one or more columns. each column represents one experimental
		repetition.}
  \item{eq}{equations to be used: Hill, Weibull, Logit, BCW, BCL, GL.}
  \item{param}{starting values for the curve fitting. Two values for Hill, Weibull, Logit, and 
		three values for BCW, BCL, GL.}
  \item{effv}{numeric vector with single or multiple effect values (0 ~ 1).}
  \item{sigLev}{The significant level for confidence intervals and Dunnett\'s test. The default
		is 0.05.}
  \item{noec}{a logical value (TRUE of FALSE). Whether to calculate NOEC and LOEC for the 
		concentration-response data.}
  \item{fig}{a logical value (TRUE of FALSE). Whether to show the concentration-response curve.}
  \item{algo}{algorithm used in the non-linear least squares fitting. 
		if package 'nls2' is installed on your R platform. The following choices 
		are available : brute-force"(alternately called 
		"grid-search"), "random-search", "plinear-brute" and "plinear-random". Reference to the 
		package 'nls2'.
		}
}
\details{The curve fitting is conditionally dependent on the package 'nls2'. 
		It will decide whether to use the built-in nonlinear least squares fitting 
		function nls or its robust counterpart 
		'nls2' (http://cran.r-project.org/web/packages/nls2/index.html).\cr
		Six sigmoidal equations are as follows:\cr
		Hill: \deqn{E = \alpha c/\left( {\beta  + c} \right)}
		Weibull: \deqn{E = 1 - \exp ( - \exp (\alpha  + \beta \lg (c)))}
		Logit: \deqn{E = {(1 + \exp ( - \alpha  - \beta \lg (c)))^{ - 1}}}
		BCW: \deqn{E = 1 - \exp \left( { - \exp \left( {\alpha  + \beta \left( 
				{\frac{{{c^\gamma } - 1}}{\gamma }} \right)} \right)} \right)}
		BCL: \deqn{E = {(1 + \exp ( - \alpha  - \beta (({c^\gamma } - 1)/\gamma )))^{ - 1}}}
		GL: \deqn{E = 1/{(1 + \exp ( - \alpha  - \beta \lg (c)))^\gamma }}
		where \eqn{E} is effect and \eqn{c} is the concentration.
		}
\value{
\item{fitInfo }{curve fitting information including the formula used to fit the concentration-
		response data. The fitted coefficients with standard errors, t test value, and p value. 
		Residual standard error and the degree of freedom are also provided.}
\item{p }{fitted coefficients of the formula}
\item{sta }{Statistics about the goodness of fit (\eqn{R^2}, \eqn{R_{adj}^2}, MAE, RMSE, AIC)}
\item{crcInfo }{a numeric matrix with the experimental concentration (x), fitted response (yhat), 
		experimental responses, lower and upper bounds of observation-based confidence intervals
		(OCI.low and OCI.up), and lower and upper bounds of function-based confidence intervals
		(FCI.low and FCI.up)
		}
\item{eci }{The confidence intervals of effect concentration at the response of effv}
\item{effvci }{The confidence intervals of the response effv}
\item{noecInfo$mat }{a matrix of experimental concentrations, Student's t-statistic, 
		F distribution at the sigLev,and sign (-1 or 1)}
\item{noecInfo$no }{non-observed effect concentration (NOEC)}
\item{noecInfo$lo }{least-observed effect concentration (LOEC)}

}
\references{
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the
Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry 
20(2):448-457.\cr
Spiess, A.-N., Neumeyer, N., 2010. An evaluation of R2 as an inadequate measure for nonlinear
models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol. 
10, 11.\cr
Gryze, S. De, Langhans, I., Vandebroek, M., 2007. Using the correct intervals for prediction: A
tutorial on tolerance intervals for ordinary least-squares regression. Chemom. Intell. Lab. 
Syst. 87, 147-154.
}
\author{xiangwei zhu}
\note{
	More sigmoidal equations (Scholze et al., 2001) can still be added into this function. Note
	that the starting values for parameters are not easy to pick in curve fitting even though at
	the	presence of nls2 package. Generally, two-parameter equations are more robust to mediocre
	starting values that three-parameter equations. One rule of thumb is to choose equations with
	less parameters at the same level of accuracy.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(cytotox)
## example 1
# Fit the concentration-response data of heavy metal Ni(2+) on MCF-7 cells
# Logit equation is selected as the fitting function with a starting value of 12 and 3.
x <- cytotox$Ni$x
expr <- cytotox$Ni$y
curveFit(x, expr, eq = 'Logit', param = c(12, 3), effv = 0.5)

## example 2
# Fit the concentration-response data of heavy metal Ni(2+) on MCF-7 cells. 
# Calculate the concentrations that cause 5\% and 50\% inhibition of the growth of MCF-7 and
# corresponding confidence intervals.
x <- cytotox$Ni$x
expr <- cytotox$Ni$y
curveFit(x, expr, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5))

## example 3
# Fit the concentration-response data of ion liquid Omim on MCF-7 cells
# Weibull equation is selected as the fitting function with a starting value of 6 and 2.
# Calculate the concentrations that cause 5\% and 50\% inhibition of the growth of MCF-7 and
# corresponding confidence intervals.
x <- cytotox$Omim$x
expr <- cytotox$Omim$y
curveFit(x, expr, eq = 'Weibull', param = c(6, 2), effv = c(0.05, 0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ curve fitting }
\keyword{concentration-response curve}
\keyword{goodness of fit}
\keyword{confidence intervals}
\keyword{function-based confidence intervals}
\keyword{observation-based confidence intervals}
\keyword{effect concentration}
