\name{vuongtest}
\alias{vuongtest}
\title{Vuoung test for non-nested models}
\description{
  The Vuong test is suitable to discriminate between two non-nested models.
}
\usage{
vuongtest(x, y,
          type = c("non.nested", "nested", "overlapping"),
          hyp = FALSE)
}
\arguments{
  \item{x}{a first fitted model of class \code{"mhurdle"},}
  \item{y}{a second fitted model of class \code{"mhurdle"},}
  \item{type}{the kind of test to be computed,}
  \item{hyp}{a boolean, \code{TRUE} if one of the models is asumed to be
    the true model.}
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
dhi <- mhurdle(durable ~ 0 | age + quant | age + quant, tobin, 
             dist = "n", corr = FALSE)
dhs <- mhurdle(durable ~ 0 | age + quant | age + quant, tobin, 
              dist = "l", corr = FALSE)
vuongtest(dhi, dhs)

}
\author{}
\keyword{}
