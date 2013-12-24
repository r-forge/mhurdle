\name{vuongtest}
\alias{vuongtest}
\title{Vuoung test for non-nested models}
\description{
  The Vuong test is suitable to discriminate between two non-nested models.
}
\usage{
vuongtest(x, y,
          type = c("non-nested", "nested", "overlapping"),
          hyp = FALSE,
          variance = c("centered", "uncentered"),
          matrix = c("large", "reduced")
         )
}
\arguments{
  \item{x}{a first fitted model of class \code{"mhurdle"},}
  \item{y}{a second fitted model of class \code{"mhurdle"},}
  \item{type}{the kind of test to be computed,}
  \item{hyp}{a boolean, \code{TRUE} if one of the models is asumed to be
    the true model,}
  \item{variance}{the variance is estimated using the \code{centered} or
    \code{uncentered} expression,}
  \item{matrix}{the W matrix can be computed using the general
    expression \code{large} or the reduced matrix \code{reduced} (only
    relevant for the nested case),}
}
\value{
  an object of class \code{"htest"}
}
% \details{
% }
\references{
  Vuong Q.H. (1989) Likelihood ratio tests for model selection and
  non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
}
\seealso{
  \code{vuong} in package \code{pscl}.
}
\examples{
data("tobin", package = "survival")
# selection double hurdle model
model110i <- mhurdle(durable ~ age |  quant | 0, tobin,  dist = "n")
# double-hurdle p-tobit
model011i <- mhurdle(durable ~ 0 |  quant | age, tobin,  dist = "n")
# Vuong test for strictly non-nested models
vuongtest(model011i, model110i)
}
\keyword{htest}
