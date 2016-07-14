\name{figPlot}
\alias{figPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot concentration response curve}
\description{Producing concentration response curves with/without confidence intervals.
}
\usage{figPlot(crcInfo, ylimit, xlabel = 'log[concentration, mol/L]', 
			ylabel = 'Inhibition [\%]', lgd = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{crcInfo}{a numeric matrix with the experimental concentration (x), fitted response (yhat),
		experimental responses, lower and upper bounds of (non-simultaneous) prediction intervals
		(PI.low and PI.up), and lower and upper bounds of (non-simultaneous) confidence intervals
		(CI.low and CI.up)}
  \item{ylimit}{range of the y axis.}
  \item{xlabel, ylabel}{plot labels.}
  \item{lgd}{legend of the plot}
}
\details{Plot the concentration response curves, experimental data, CI, and PI in one figure.
}

\author{xiangwei zhu}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## example 1
# 
x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
f <- curveFit(x, expr, eq = 'Logit', param = c(26, 4), effv = c(0.05, 0.5))
figPlot(f$crcInfo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{concentration response curve}
