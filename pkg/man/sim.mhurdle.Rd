\name{sim.mhurdle}
\alias{sim.mhurdle}
\title{Function to simulate limited dependent variables}
\description{
  
  This function enables the simulation of limited dependent variables
  
}
\usage{
sim.mhurdle(param = c(0.5, 0.5, 0.5, 0.5, 0.5), rho = 0, n = 200, 
         sel = TRUE, ifr = FALSE,
         dist = c("l", "t", "n"), res = FALSE) 
}
\arguments{

  \item{param}{the parameters of the simulation, see details,}
  \item{rho}{the coefficient of correlation,}
  \item{n}{the size of the sample,}
  \item{sel}{a boolean, \code{TRUE} if 0 may occurs in a selection
    process,}
  \item{ifr}{a boolean, \code{TRUE} if 0 may occurs because the good is
    bought infrequently,}
  \item{dist}{the distribution of the error of the consumption
    equation: one of \code{"n"} (normal), \code{"l"} (log-normal) or \code{"t"} (truncated
    normal),}
  \item{res}{a boolean, \code{TRUE} if 0 may occurs as a corner
    solution,}
}
\value{
  a \code{data.frame} containing \code{y} the dependent variable, \code{x} the
  explanatory variable in the consumption equation and \code{s} the
  explanatory variable in the selection/infrequency equation.
}

\details{
  
  The \code{param} argument is a vector of length 5 containing the
  intercept and the coefficient for the sel/ifruency equation,
  the intercept and the coefficient for the consumption equation, and
  the standard deviation of the error term in the consumption equation.
  The two explanatory variables are standard normal variables.
}
\references{
}
\seealso{
}
\examples{

# simulation of a data set of 1000 observations for Cragg's dependent
# double hurdle model

d <- sim.mhurdle(c(0,1,1,2), rho = 0.6, n = 1000, sel = TRUE,
              dist = "n", res = TRUE)

}
\author{}
\keyword{}

