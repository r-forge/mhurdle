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
data("Interview", package = "mhurdle")
# dependent double hurdle model
dhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
              dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)

# a double hurdle p-tobit model
ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
              dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
vuongtest(dhm, ptm)
}
\keyword{htest}
