\name{tuneFit}
\alias{tuneFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Optimal Starting values for Curve Fitting}
\description{ This generic function first seaches optimal starting values based on 
trial and error. The the concentration-responses will be fitted using the optimal 
starting values. The statistics for goodness of fit
is evaluated by the following statistics: coefficient of
determination (\eqn{R^2}), adjusted coefficient of determination (\eqn{R_{adj}^2}),
root mean squared error (RMSE), mean absolute error (MAE), Akaike information criterion (AIC),
bias-corrected Akaike information criterion(AICc), and Bayesian information criterion (BIC).
Thirteen sigmoidal models ("Hill", "Hill_two", "Hill_three", "Hill_four",
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and four
J-shaped models ("Brain_Consens", "BCV", "Biphasic", "Hill_six")
are provided to fit concentration-response data. 
}
\usage{
tuneFit(conc, rspn, eq = 'Weibull', effv, highBar = 5000, bar = 1000, sav = FALSE)
}

\arguments{
  \item{conc}{a numeric vector (matrix) of experimental concentrations}
  \item{rspn}{a numeric vector (matrix) of responses corresponding to conc, it should have the
  same length (or rows and columns for matrix) as conc.}
  \item{eq}{models for curve fitting: "Hill", "Hill_two", "Hill_three", "Hill_four",
            "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
            "BCW", "BCL", "GL", "Brain_Consens", "BCV", "Biphasic", "Hill_six".}
  \item{effv}{one numeric value [0, 1].}
  \item{highBar}{if the number of starting values exceeds highBar, A random sample of 
					starting values will be taken.}
  \item{bar}{the number of random samples}
  \item{sav}{save fitting results in a text file (TRUE/FALSE).}
}
\details{The starting values are indispensable for nonlinear optimization algorithms. 
Guessing appropriate starting values requires significant experience and statistical 
 insight. However, users in the area of toxicology always get confused and frustrated 
 in trying so many times to fit a concentration-response curve successfully. 
	tuneFit provides high frequency trial and error approach to find appropriate starting 
	values for users. It will deploy those starting values one by one until 
	finding the right one.
	}
\value{
\item{sta }{Statistics about the goodness of fit (\eqn{R^2}, \eqn{R_{adj}^2}, MAE, RMSE, AIC,
			AICc, and BIC)}
}

\author{xiangwei zhu}
\note{
	tuneFit will load the file staval.rda which contains hundreds of starting values. However, 
	those starting values are also limited. We encourage users to send fitted coefficients to us 
	to enlarge the coverage of staval.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# Fit the non-monotonic concentration-response data

x <- hormesis$OmimCl$x
expr <- hormesis$OmimCl$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Biphasic')

x <- hormesis$HmimCl$x
expr <- hormesis$HmimCl$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Biphasic', effv = 0.05)


## example 2
# Fit the concentration-response data of heavy metal Ni(2+) on MCF-7 cells.
# Calculate the concentrations that cause 5\% inhibition of the growth of MCF-7

x <- cytotox$Ni$x
expr <- cytotox$Ni$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Logit', effv = 0.05)

## example 3
# Fit the concentration-response data of Paromomycin Sulfate (PAR) on photobacteria.
# Calculate the concentrations that cause 50\% inhibition of the growth of photobacteria 

x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Logit', effv = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ curve fitting }
\keyword{goodness of fit}
\keyword{trial and error}
\keyword{effect concentration}
