\name{vuongtest}
\alias{vuongtest}
\title{Vuoung test for non-nested models}
\description{
  The Vuong test is suitable to discriminate between two non-nested models.
}
\usage{
vuongtest(x, y)
}
\arguments{
  \item{x}{a first fitted model of class \code{"mhurdle"},}
  \item{y}{a second fitted model of class \code{"mhurdle"},}
}
\value{
  an object of class \code{"htest"}
}
\details{
}
\references{
  Vuong Q.H. (1989) Likelihood ratio tests for model selection and
  non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
}
\seealso{
  \code{vuong} in package \code{pscl}.
}
\examples{
## dependent simple hurdle ptobit vs an
## independent double hurdle ptobit

data("tobin",package="survival")
dhi <- mhurdle(durable ~ age + quant | age + quant, tobin, 
            sel = FALSE, ifr = TRUE, dist = "n", res = TRUE, corr = FALSE)
dhs <- mhurdle(durable ~ age + quant | age + quant, tobin, 
            sel = FALSE, ifr = TRUE,  dist = "l", res = FALSE, corr = TRUE)
vuongtest(dhi, dhs)

}
\author{}
\keyword{}
