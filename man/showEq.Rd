\name{showEq}
\alias{showEq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{List the Expression of Six Sigmoidal Equations}
\description{
Six equations are employed to fit the sigmoidal concentration-response data.
}
\usage{
showEq()
}
%- maybe also 'usage' for other objects documented here.
\details{
In Hill equation with slope parameter \eqn{n = 1}, parameter Beta is equivalent to the effect 
concentration causing 50% of maximal response (EC50). In Weibull and Logit equations, 
Alpha and Beta are the parameters (location and slope, respectively) to be estimated. 
}
\value{The formula of six equations (with abbr.) will be listed on the screen.}
\references{
Scholze, M., Boedeker, W., Faust, M., Backhaus, T., Altenburger, R., Grimme, L.H., 2001. A general
best-fit method for concentration-response curves and the estimation of low-effect concentrations.
Environ. Toxicol. Chem. 20, 448-457
}
\author{xiangwei zhu}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{curveFit}}
}
\examples{
## show those equations
showEq()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sigmoidal equation }
