\name{rsq}
\alias{rsq}
\title{R squared and pseudo R squared}
\description{
  This function computes the R squared for multiple hurdle models. The
  measure is an extension of the R squared of linear model or may be based on
  the likelihood and may measure the goodness of fit for the zero part,
  the positive part of the model or the whole model.
}
\usage{
rsq(object, type = c("coefdet", "lratio"),
            adj = FALSE, r2pos = c("rss", "ess", "cor"))
}
\arguments{
  \item{object}{an object of class \code{"mhurdle"},}
  \item{type}{one of \code{"coefdet"} or \code{"lratio"} to select a
    linear measure of the R squared or a Mc Fadden like measure based on
    the likelihood function,}
  \item{adj}{if \code{TRUE} a correction for the degrees of freedom is
    performed,}
  \item{r2pos}{only for linear R squared, should the positive part of
    the R squared be computed using the residual sum of squares
    (\code{"rss"}), the explained sum of squares (\code{"ess"}) or the
    coefficient of correlation between the fitted values and the
    response (\code{cor}).
    }
}
\value{
  a numerical value
}
%\details{
%}
% \seealso{
% }

\references{
  
 McFadden D (1974). The Measurement of Urban Travel Demand. Journal of
 Public Economics, 3, 303-328.
 
}
\examples{

# A dependent double hurdle example
data8d <- simhurdle(n = 1000, h3 = FALSE, h1 = TRUE, corr = "h1", dist = "n")
z8d <- mhurdle(y ~ x1 | x2 | 0, data = data8d, dist = "n", corr = "h1", method = "bfgs")

rsq(z8d, type = "lratio")
rsq(z8d, type = "coefdet", adj = TRUE, r2pos = "rss")

}
\keyword{htest}
