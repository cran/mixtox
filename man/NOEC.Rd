\name{NOEC}
\alias{NOEC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{NOEC and LOEC Calculation}
\description{Dunnett's test is used to calculate the NOEC and LOEC}
\usage{NOEC(x, expr, sigLev = 0.05)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector of experimental concentrations}
  \item{expr}{ a numeric matrix. each column is an experimental repetition and at least three
		repetitions are required.}
  \item{sigLev}{The significant level for Dunnett's test. The default is 0.05.}
}
\details{
Dunnett's test (Dunnett, 1964) is performed by computing a Student's t-statistic for 
each experimental, or treatment, group where the statistic compares the treatment 
group to a single control group. Since each comparison has the same control in common, 
the procedure incorporates the dependencies between these comparisons. In particular, 
the t-statistics are all derived from the same estimate of the error variance 
which is obtained by pooling the sums of squares for error across all 
(treatment and control) groups. The formal test statistic for Dunnett's test 
is either the largest in absolute value of these t-statistics (if a 
two-tailed test is required), or the most negative or most positive of the 
t-statistics (if a one-tailed test is required).
}
\value{
\item{mat }{a matrix of experimental concentrations, Student's t-statistic, \eqn{F} distribution 
			at the sigLev, and sign (-1 or 1)}
\item{no }{NOEC}
\item{lo }{LOEC}
}
\references{
Dunnett, C.W., 1964. New tables for multiple comparisons with a control. Biometrics 30, 482-491.
}
\author{xiangwei zhu}
\note{
At least three repetitions are required to calculate the NOEC and LOEC
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{curveFit}}
}
\examples{
data(cytotox)

## example 1
# calcualte the NOEC and LOEC of heavy metal Ni(2+) on the MCF-7 cells at the default significance 
# level of 0.05
x <- cytotox$Ni$x
expr <- cytotox$Ni$y
NOEC(x, expr)

## example 2
# calcualte the NOEC and LOEC of heavy metal Cu(2+) on the MCF-7 cells at the default significance
# level of 0.01
x <- cytotox$Cu$x
expr <- cytotox$Cu$y
NOEC(x, expr, sigLev = 0.01)
}
\keyword{NOEC}
\keyword{LOEC}
