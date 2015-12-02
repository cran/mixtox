\name{curveFit}
\alias{curveFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Curve Fitting}
\description{
Thirteen monotonic(sigmoidal) equations ("Hill", "Hill_two", "Hill_three", "Hill_four",
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and six
non-monotonic(J-shaped) equations ("Brain_Consens", "BCV", "Beckon", "Biphasic", "Hill_six")
are provided to fit concentration-response data. The statistics for goodness of fit
is evaluated by the following statistics: coefficient of
determination (\eqn{R^2}), adjusted coefficient of determination (\eqn{R_{adj}^2}),
root mean squared error (RMSE), mean absolute error (MAE), Akaike information criterion (AIC),
bias-corrected Akaike information criterion(AICc), and Bayesian information criterion (BIC).
}
\usage{
curveFit(x, expr, eq , param, effv, fig = TRUE, ylimit, 
         xlabel = 'lg[concentration, mol/L]', ylabel = 'Inhibition [\%]', 
		 sigLev = 0.05, noec = TRUE, algo = "default")

}

\arguments{
  \item{x}{a numeric vector of experimental concentrations}
  \item{expr}{a numeric matrix with one or more columns. each column represents one experimental 
		repetition.}
  \item{eq}{equations to be called: "Hill", "Hill_two", "Hill_three", "Hill_four",
            "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
            "BCW", "BCL", "GL", "Brain_Consens", "BCV", "Beckon", "Biphasic", "Hill_six".}
  \item{param}{starting values for curve fitting.}
  \item{effv}{numeric vector with single or multiple effect values.}
  \item{fig}{a logical value (TRUE of FALSE). Whether to show the concentration-response curve.}
  \item{ylimit}{a two value vector defines the lower and upper limits of the figure.}
  \item{xlabel}{define x-axis label.}
  \item{ylabel}{define y-axis label.}
  \item{sigLev}{The significant level for confidence intervals and Dunnett\'s test. The default
        is 0.05.}
  \item{noec}{a logical value (TRUE of FALSE) to determine the calculation of NOEC and LOEC.}
  \item{algo}{algorithm used in the non-linear least squares fitting.
		if package 'nls2' is installed on your R platform. The following choices
		are available : brute-force"(alternately called
		"grid-search"), "random-search", "plinear-brute" and "plinear-random". Reference to the
		package 'nls2'.
		}
}
\details{Curve fitting is dependent on the package 'nls2'
        (http://cran.r-project.org/web/packages/nls2/index.html).\cr
		Monotonic(sigmoidal) equations are listed as follows:\cr
		Hill: \deqn{E = 1/\left( {1 + {{\left( {\alpha /c} \right)}^\beta }} \right)}
		Hill_two: \deqn{E = \beta c/\left( {\alpha  + c} \right)}
		Hill_three: \deqn{E = \gamma /\left( {1 + {{\left( {\alpha /c} \right)}^\beta }} \right)}
		Hill_four: \deqn{E = \delta  + \left( {\gamma  - \delta } \right)/\left( {1 +
							{{\left( {\alpha /c} \right)}^\beta }} \right)}
		where \eqn{alpha} = EC50, \eqn{beta} = m(Hill coefficient), \eqn{gamma} = Top,
		and \eqn{delta} = Bottom \cr

		Weibull: \deqn{E = 1 - \exp ( - \exp (\alpha  + \beta \lg (c)))}
		Weibull_three: 	\deqn{E = \gamma \left( {1 - \exp \left( { - \exp \left( {\alpha  +
		                       \beta \lg \left( c \right)} \right)} \right)} \right)}
		Weibull_four: \deqn{E = \gamma  + \left( {\delta  - \gamma } \right)\exp \left(
		                   { - \exp \left({\alpha  + \beta \lg \left( c \right)} \right)} \right)}
		Logit: \deqn{E = {(1 + \exp ( - \alpha  - \beta \lg (c)))^{ - 1}}}
		Logit_three: \deqn{E = \gamma /\left( {1 + \exp \left( {\left( { - \alpha } \right)
		                - \beta \lg \left( c \right)} \right)} \right)}
		Logit_four: \deqn{E = \delta  + \left( {\gamma  - \delta } \right)/\left( {1 + \exp \left
		            ( {\left( { - \alpha } \right) - \beta \lg \left( c \right)} \right)} \right)}

		where \eqn{alpha} is the location parameter and \eqn{beta} slope parameter.
		\eqn{gamma} = Top, and \eqn{delta} = Bottom \cr

		BCW: \deqn{E = 1 - \exp \left( { - \exp \left( {\alpha  + \beta \left(
				{\frac{{{c^\gamma } - 1}}{\gamma }} \right)} \right)} \right)}
		BCL: \deqn{E = {(1 + \exp ( - \alpha  - \beta (({c^\gamma } - 1)/\gamma )))^{ - 1}}}
		GL: \deqn{E = 1/{(1 + \exp ( - \alpha  - \beta \lg (c)))^\gamma }}

		Non-monotonic(J-shaped) equations:\cr
		Hill_six: \deqn{E = \left( {\gamma /{{\left( {1 + \alpha /c} \right)}^\beta }} \right)/
			                   (gamma\_one/{\left( {1 + alpha\_one/c} \right)^{beta\_one}}}

		Brain_Consens: \deqn{E = 1 - \left( {1 + \alpha c} \right)/\left( {1 + \exp \left(
		                      {\beta \gamma } \right){c^\beta }} \right)}
		where alpha is the initial rate of increase at low concentration, beta the way in which
		response decreases with concentration, and gamma no simple interpretation.\cr

		BCV: \deqn{E = 1 - \alpha \left( {1 + \beta c} \right)/\left( {1 + \left( {1 + 2\beta
					\gamma } \right){{\left( {c/\gamma } \right)}^\delta }} \right)}
		where alpha is untreated control, beta the initial rate of increase at low concentration,
		gamma the concentration cause 50\% inhibition, and delta no simple interpretation.\cr

		Cedergreen: \deqn{E = 1 - \left( {1 + \alpha \exp \left( { - 1/\left( {{c^\beta }}
						\right)} \right)} \right)/\left( {1 + \exp \left( {\gamma \left({\ln
						\left( c \right) - \ln \left( \delta  \right)} \right)} \right)} \right)}
		where alpha the initial rate of increase at low concentration, beta the rat of the
		hormetic effect manifests itself, gamma the steepness of the curve after the
		maximum hormetic effect, and delta the lower bound on the EC50 level.\cr

		Beckon: \deqn{E = \left( {\alpha  + \left( {1 - \alpha } \right)/\left( {1 + {{\left(
						{\beta /c} \right)}^\gamma }} \right)} \right)/\left( {1 + {{\left(
						{c/\delta } \right)}^\varepsilon }} \right)}
		where alpha is the minimum effect that would be approached by the downslope in the absence
		of the upslope, beta the concentration at the midpoint of the falling slope, gamma the
		steepness of the rising(positive) slope, delta the concentration at the midpoint of the
		rising slope, and epsilon the steepness of the falling(negative) slope.\cr

		Biphasic: \deqn{E = \alpha  - \alpha /\left( {1 + {{10}^{\left( {\left( {c - \beta }
						\right)\gamma } \right)}}} \right) + \left( {1 - \alpha } \right)/\left
						( {1 + {{10}^{\left( {\left( {\delta  - c} \right)\varepsilon } \right)}}}
						\right)}
		where alpha is the minim effect that would be approached by the downslope in the absence
		of the upslope, beta the concentration at the midpoint of the falling slope, gamma the
		steepness of the rising(positive) slope, delta the concentration at the midpoint of the
		rising slope, and epsilon the steepness of the falling(negative) slope.\cr
		In all, \eqn{E} is effect and \eqn{c} is the concentration.
	}
\value{
\item{fitInfo }{curve fitting information including the formula used to fit the concentration-
		response data. The fitted coefficients with standard errors, t test value, and p value.
		Residual standard error and the degree of freedom are also provided.}
\item{p }{fitted coefficients of the formula}
\item{sta }{Statistics about the goodness of fit (\eqn{R^2}, \eqn{R_{adj}^2}, MAE, RMSE, AIC,
			AICc, and BIC)}
\item{crcInfo }{a numeric matrix with the experimental concentration (x), fitted response (yhat),
		experimental responses, lower and upper bounds of observation-based confidence intervals
		(OCI.low and OCI.up), and lower and upper bounds of function-based confidence intervals
		(FCI.low and FCI.up)
		}
\item{eci }{confidence intervals of effect concentration at the response of effv}
\item{effvci }{confidence intervals of the response effv for monotonic equations}
\item{minx }{concentration to induce the largest stimulation for non-monotonic equations}
\item{miny }{the largest stimulation, only for non-monotonic equations}
\item{noecInfo$mat }{a matrix of experimental concentrations, Student's t-statistic,
		F distribution at the sigLev,and sign (-1 or 1)}
\item{noecInfo$no }{non-observed effect concentration (NOEC)}
\item{noecInfo$lo }{least-observed effect concentration (LOEC)}

}
\references{
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the
Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry
20(2):448-457.\cr
Zhu X-W, et.al. 2013. Modeling non-monotonic dose-response relationships: Model evaluation
 and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136.\cr
Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures 
containing partial agonists. J. Theor. Biol. 259:469-477.\cr	
Spiess, A.-N., Neumeyer, N., 2010. An evaluation of R2 as an inadequate measure for nonlinear
models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol.
10, 11.\cr
Di Veroli GY, Fornari C, Goldlust I, Mills G, Koh SB, Bramhall JL, et al. 2015. An automated 
fitting procedure and software for dose-response curves with multiphasic features. 
Scitific Report 5: 14701.\cr
Gryze, S. De, Langhans, I., Vandebroek, M., 2007. Using the correct intervals for prediction: A
tutorial on tolerance intervals for ordinary least-squares regression. Chemom. Intell. Lab.
Syst. 87, 147-154.
}
\author{xiangwei zhu}
\note{
	Note that the starting values for parameters are not easy to pick in curve fitting 
	Generally, two-parameter 
	equations are more robust to bad 
	starting values than three-parameter equations. One rule of thumb is to choose equations with
	less parameters at the same level of accuracy.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# Fit the non-monotonic concentration-response data

x <- hormesis$OmimCl$x
expr <- hormesis$OmimCl$y
curveFit(x, expr, eq = 'Biphasic', param = c(-0.34, 0.001, 884, 0.01, 128), effv = 0.5)

x <- hormesis$HmimCl$x
expr <- hormesis$HmimCl$y
curveFit(x, expr, eq = 'Biphasic', param = c(-0.59, 0.001, 160,0.05, 19),  effv = c(0.05, 0.5))

x <- hormesis$ACN$x
expr <- hormesis$ACN$y
curveFit(x, expr, eq = 'Brain_Consens', param = c(2.5, 2.8, 0.6, 2.44),  effv = c(0.05, 0.5))

x <- hormesis$Acetone$x
expr <- hormesis$Acetone$y
curveFit(x, expr, eq = 'BCV', param = c(1.0, 3.8, 0.6, 2.44),  effv = c(0.05, 0.5))

## example 2
# Fit the concentration-response data of heavy metal Ni(2+) on MCF-7 cells.
# Calculate the concentrations that cause 5\% and 50\% inhibition of the growth of MCF-7 and
# corresponding confidence intervals.

x <- cytotox$Ni$x
expr <- cytotox$Ni$y
curveFit(x, expr, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5))

## example 3
# Fit the concentration-response data of Paromomycin Sulfate (PAR) on photobacteria.
# Calculate the concentrations that cause 5\% and 50\% inhibition of the growth of photobacteria 
# and corresponding confidence intervals.

x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
curveFit(x, expr, eq = 'Logit', param = c(26, 4), effv = c(0.05, 0.5))
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
