\name{sim.mhurdle}
\alias{sim.mhurdle}
\title{Function to simulate limited dependent variables}
\description{
  
  This function enables the simulation of limited dependent variables
  
}
\usage{
sim.mhurdle(reg = c(0.5, 0.5, 1), sel = NULL, ifr = NULL, rho = 0,
           n = 200, dist = c("l", "t", "n"))
}
\arguments{

  \item{reg}{the parameters of the regression equation (of length three
    : the intercept, the slope and the standard deviation of the error,}
  \item{sel}{the parameters of the selection equation (of lenght two :
    the intercept and the slope). If null, there is no selection
    process,}
  \item{ifr}{the parameters of the purchase equation (of lenght two :
    the intercept and the slope). If null, there is no purchase
    process,}
  \item{rho}{the coefficient of correlation between the regression and
    the selection equations : only relevant if there is a selection process,}
  \item{n}{the size of the sample,}
  \item{dist}{the distribution of the error of the consumption equation:
    one of \code{"n"} (normal), \code{"l"} (log-normal) or \code{"t"}
    (truncated normal),}
}
\value{
  a \code{data.frame} containing \code{y} the dependent variable,
  \code{x} the explanatory variable in the consumption, \code{s} the
  explanatory variable in the selection equation and \code{p} the
  explanatory variable in the infrequency equation.
}

\details{
  The explanatory variables are standard normal variables.
}
\examples{

# simulation of a data set of 1000 observations for Cragg's dependent
# double hurdle model

d <- sim.mhurdle(c(0.5, 0.5, 2), sel = c(0.5, 0.5), n = 1000, dist = "n")

}
\author{}
\keyword{}

