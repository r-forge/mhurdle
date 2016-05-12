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
# a small subset of the Intview data
set.seed(1)
data("Intview", package = "mhurdle")
SIntview <- Intview[sample(1:nrow(Intview), 1E+03), ]

# dependent double hurdle model
dhm <- mhurdle(vacations ~ car + size | lincu + lincu2 + lincu3 | 0, SIntview,
              dist = "ln2", method = "bfgs", corr = TRUE)

# a double hurdle p-tobit model
ptm <- mhurdle(vacations ~ 0 | lincu + lincu2 + lincu3 | car + size, SIntview,
              dist = "ln2", method = "bfgs", corr = TRUE)
vuongtest(dhm, ptm)
}
\keyword{htest}
