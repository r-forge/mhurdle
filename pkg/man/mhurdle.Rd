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
     dist = c("ln","tn","n", "bc", "ihs", "ln2", "bc2"),
     corr = NULL, ...)
\method{coef}{mhurdle}(object,
   which = c("all", "h1", "h2", "h3", "sd", "corr", "tr", "mu"), ...)
\method{vcov}{mhurdle}(object,
   which = c("all", "h1", "h2", "h3", "sd", "corr", "tr", "mu"), ...)
\method{logLik}{mhurdle}(object, naive = FALSE, ...)
\method{print}{mhurdle}(x, digits = max(3, getOption("digits") - 2),
                     width = getOption("width"), ...)
\method{summary}{mhurdle}(object, ...)
\method{print}{summary.mhurdle}(x, digits = max(3, getOption("digits") - 2),
   width = getOption("width"), ...)

\method{fitted}{mhurdle}(object,
   which = c("all", "zero", "positive"), ...)
\method{predict}{mhurdle}(object, newdata = NULL, ...)
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
  \item{corr}{indicates whether the errors of the different equations
   are correlated. For models with two equations, this can be either
   \code{"d"} for dependent and \code{"i"} for independent. For models
   with three equations, this should be a character of length three
   containing values of "i" and "d",}
  \item{naive}{a boolean, it \code{TRUE}, the likelihood of the naive
    model is returned,}
  \item{object,x}{an object of class \code{"mhurdle"},}
  \item{new}{an updated formula for the \code{update} method,}
  \item{digits}{see \code{\link{print}},}
  \item{width}{see \code{\link{print}},}
  \item{which}{which coefficients or covariances should be extracted ? Those of the
    selection (\code{"h1"}), consumption (\code{"h2"}) or purchase
    (\code{"h3"}) equation, the other coefficients \code{"other"} (the
    standard error and the coefficient of corr), the standard error
    (\code{"sigma"}) or the coefficient of correlation (\code{"rho"}),}
  \item{\dots}{further arguments.}
  
}
\value{
  

  an object of class \code{c("mhurdle", "maxLik")}.
    
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
    the selection equation, the regression equation, the infrequency of
    purchase equation and the other coefficients (the standard deviation
    of the error term and the coefficient of correlation if \code{corr = TRUE}),}
  \item{formula}{the model formula, an object of class \code{Formula},}
  \item{call}{the call,}
  \item{rho}{the lagrange multiplier test of no correlation.}
    
}

\details{
  
  \code{mhurdle} fits models for which the dependent variable is zero for
a part of the sample. Null values of the dependent variable may occurs
because of one or several mechanisms : good rejection, lack of
ressources and purchase infrequency. The model is described using a
three-parts formula : the first part describes the selection process if
any, the second part the regression equation and the third part the
purchase infrequency process. \code{y ~ 0 | x1 + x2 | z1 + z2} means
that there is no selection process. \code{y ~ w1 + w2 | x1 + x2 | 0} and
\code{y ~ w1 + w2 | x1 + x2} describe the same model with no purchase
infrequency process. The second part is mandatory, it explains the
positive values of the dependant variable. The \code{dist} argument
indicates the distribution of the error term. If \code{dist = "n"}, the
error term is normal and (at least part of) the zero observations are
also explained by the second part as the result of a corner
solution. Several models described in the litterature are obtained as
special cases :

A model with a formula like \code{y~0|x1+x2} and \code{dist="n"} is the
Tobit model proposed by Tobin (1958).

\code{y~w1+w2|x1+x2} and \code{dist="l"} or \code{dist="t"} is the
single hurdle model proposed by Cragg (1971). With \code{dist="n"}, the
double hurdle model also proposed by Cragg (1971) is obtained. With
\code{corr="h1"} we get the correlated version of this model described
by Blundell (1987). 

\code{y~0|x1+x2|z1+z2} is the P-Tobit model of Deaton and Irish (1984),
which can be a single hurdle model if \code{dist="t"} or \code{dist="l"}
or a double hurdle model if \code{dist="n"}.

}
 
\references{
  
  Blundell R, Meghir C (1987). Bivariate Alternatives to the Tobit
  Model. Journal of Econometrics, 34, 179-200.

  Cragg JG (1971). Some Statistical Models for Limited Dependent
  Variables with Applications for the Demand for Durable
  Goods. Econometrica, 39(5), 829-44.

  Deaton A, Irish M (1984). A Statistical Model for Zero Expenditures in
  Household Budgets.  Journal of Public Economics, 23, 59-80.

  Tobin J (1958). Estimation of Relationships for Limited Dependent
  Variables. Econometrica, 26(1), 24-36.

}
% \seealso{
% }
\examples{
data("tobin", package = "survival")
#  tobit model
model010 <- mhurdle(durable ~ 0 | age + quant | 0, tobin, dist = "n")
#  independent double hurdle model
model110i <- mhurdle(durable ~ age |  quant | 0, tobin,  dist = "n")
# Cragg log-normal single hurdle model
model100il <- mhurdle(durable ~ age |  quant | 0, tobin, dist = "ln")
# Cragg truncated-normal single hurdle model
model100it <- update(model100il, dist = "tn")
# a double-hurdle p-tobit
model011i <- mhurdle(durable ~ 0 |  quant | age, tobin,  dist = "n")
}
% \author{}
\keyword{regression}
