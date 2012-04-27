\name{Comics}
\docType{data}
\alias{Comics}
\title{Comics}
\description{
  a cross section from 2003

  \emph{number of observations} : 5159

  \emph{observation} : individuals

  \emph{country} : France
}
\usage{
data(Comics)
}
\format{A dataframe containing :
  \describe{
    
  \item{area}{one of \code{rural}, \code{small}, \code{medium},
    \code{large} and \code{paris},}
  \item{income}{the income of the household (in thousands of euros per
    year),}
  \item{uc}{the number of consumption unit (one for the first two
    adults, one half for other members of the household),}
  \item{size}{the number of persons in the household,}
  \item{age}{the age of the person,}
  \item{empl}{the kind of job, a factor with 9 levels,}
  \item{gender}{one of \code{male} and \code{female},}
  \item{couple}{does the person live in couple ?  A factor with levels
    \code{yes} and \code{no},}
  \item{educ}{the number of years of education.}
}
}
\source{
  Web site of the ``Institut National de la Statistique et des Etudes
    Economiques'', \url{http://www.insee.fr/}
}

\keyword{datasets}
