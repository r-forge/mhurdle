\name{r.squared}
\alias{r.squared}
\title{R squared and pseudo R squared}
\description{
  This function computes the R squared for multiple hurdle models. The
  measure is linear or may be based on the likelihood and may measure
  the goodness of fit for the zero part, the positive part of the model
  or the whole model.
}
\usage{
r.squared(object, which = c("all", "zero", "positive"))
}
\arguments{
  \item{object}{an object of class \code{"mhurdle"},}
  \item{which}{which part of the model the goodness of fit has to be
    computed, one of \code{"all"} (the whole model), \code{"zero"} (null
    versus positive values), \code{"positive"} (values for the truncated
    sample of positive values),}
}
\value{
  a numerical value
}
\details{
  When \code{type="linear"}, a linear $R^2$ is computed, when
  \code{type="mcfadden"}, a pseudo-R$^2$ based on the likelihood of the
  ``naive'' model is computed.
}
\seealso{
}
\examples{
data("tobin", package = "survival")
dhi <- mhurdle(durable ~ age + quant | age + quant | 0, tobin, 
            dist = "n", corr = FALSE)
r.squared(dhi)
r.squared(dhi, "zero")
r.squared(dhi, "positive")

}
\keyword{}
