\name{mhurdle}
\alias{mhurdle}
\alias{coef.mhurdle}
\alias{vcov.mhurdle}
\alias{logLik.mhurdle}
\alias{print.mhurdle}
\alias{summary.mhurdle}
\alias{print.summary.mhurdle}
\alias{predict.mhurdle}
\alias{update.mhurdle}
\alias{fitted.mhurdle}

\title{Estimation of limited dependent variable models}
\description{mhurdle fits a large set of models relevant when the dependent
  variable is 0 for a part of the sample. 
}
\usage{
 
mhurdle(formula, data, subset, weights, na.action,
     start = NULL,
     dist = c("l","t","n"),
     corr = FALSE, ...)

\method{coef}{mhurdle}(object,
                    which = c("all", "sel", "ifr", "reg",
                              "other", "sigma", "rho"),
                     ...)
\method{vcov}{mhurdle}(object,
                    which = c("all", "sel", "ifr", "reg",
                              "other", "sigma", "rho"),
                     ...)
\method{logLik}{mhurdle}(object,
                    which = c("all", "zero", "positive"),
                    naive = FALSE,
                    ...)
\method{print}{mhurdle}(x, digits = max(3, getOption("digits") - 2),
                     width = getOption("width"), ...)
\method{summary}{mhurdle}(object, ...)
\method{print}{summary.mhurdle}(x, digits = max(3, getOption("digits") - 2),
                             width = getOption("width"), ...)

\method{fitted}{mhurdle}(object,
                         which = c("all", "zero", "positive"),
                         ...)
\method{predict}{mhurdle}(object,
                       newdata = NULL, ...)
\method{update}{mhurdle}(object, new, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fitted,}
  \item{data}{a \code{data.frame},}
  \item{newdata}{a \code{data.frame} for which the predictions should be computed,}
  \item{subset}{see \code{\link{lm}},}
  \item{weights}{see \code{\link{lm}},}
  \item{na.action}{see \code{\link{lm}},}
  \item{start}{starting values,}
  \item{dist}{the distribution of the error of the consumption
    equation: one of \code{"n"} (normal), \code{"l"} (log-normal) or \code{"t"} (truncated
    normal),}
  \item{corr}{a boolean, \code{TRUE} if the two error terms are
    correlated,}
  \item{naive}{a boolean, it \code{TRUE}, the likelihood of the naive
    model is returned,}
  \item{object,x}{an object of class \code{"mhurdle"},}
  \item{new}{an updated formula for the \code{update} method,}
  \item{digits}{see \code{\link{print}},}
  \item{width}{see \code{\link{print}},}
  \item{which}{which coefficients or covariances should be extracted ? Those of the
    selection (\code{"sel"}), infrequency (\code{"ifr"}) or regression
    (\code{"reg"}) equation, the other coefficients \code{"other"} (the
    standard error and the coefficient of corr), the standard error
    (\code{"sigma"}) or the coefficient of correlation (\code{rho}),}
  \item{\dots}{further arguments.}
  
}
\value{
  

  an object of class \code{c("mhurdle","maxLik")}.
    
  A \code{"mhurdle"} object has the following elements :
  
  %  \describe{
  \item{coefficients}{the vector of coefficients,}
  \item{vcov}{the covariance matrix of the coefficients,}
  \item{fitted.values}{a matrix of fitted.values, the first column being
  the probability of 0 and the second one the mean values for the
  positive observations,}
  \item{logLik}{the log-likelihood,}
  \item{gradient}{the gradient at convergence,}
  \item{model}{a data.frame containing the variables used for the
    estimation,}
  \item{coef.names}{a list containing the names of the coefficients in
    the sel equation, in the dist equation and the other
    coefficients (the standard deviation of the error term and the
    corr coefficient if \code{corr = TRUE}),}
  \item{formula}{the model formula, an object of class \code{Formula},}
  \item{call}{the call,}
  \item{rho}{the lagrange multiplier test of no correlation.}
    
}
\details{
}
\references{
}
\seealso{
}
\examples{
# A dependent double hurdle example
data("tobin", package = "survival")
x <- mhurdle(durable ~ age+quant | age+quant | 0, data = tobin, dist = "n",
          corr = TRUE, method = "bfgs")
# A p-tobit with a log-normal distribution
x <- mhurdle(durable ~ 0 | age+quant | age+quant, data = tobin, dist = "l",
          method = "bfgs")

}
\author{}
\keyword{}
