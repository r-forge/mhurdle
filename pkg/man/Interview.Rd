\name{Interview}
\docType{data}
\alias{Interview}
\title{Interview}
\description{
  a cross section from 2014

  \emph{number of observations} : 1000

  \emph{observation} : households

  \emph{country} : United-States
}
\usage{
data(Interview)
}
\format{A dataframe containing :
  \describe{

    \item{month}{the month of the interview,}
    
    \item{size}{the number of person in the household,}

    \item{cu}{the number of consumption units in the household,}

    \item{income}{the income of the household for the 12 month before
      the interview,}

    \item{linc}{the logarithme of the net income per consumption unit
      divided by its mean,}

    \item{linc2}{the square of \code{link},}

    \item{smsa}{does the household live in a SMSA (\code{yes} or \code{no}),}

    \item{sex}{the sex of the reference person of the household (\code{male} and \code{female}),}

    \item{race}{the race of the head of the household, one of \code{white}, \code{black}, \code{indian}, \code{asian}, \code{pacific} and \code{multirace},}

    \item{hispanic}{is the reference person of the household is hispanic (\code{no} or \code{yes}),}

    \item{educ}{the number of year of education of the reference person of the household,}

    \item{age}{the age of the reference person of the household - 50,}

    \item{age2}{the square of \code{age}}

    \item{car}{cars in the household,}

    \item{food}{food,}
    \item{alcool}{,}
    \item{housing}{,}
    \item{apparel}{,}
    \item{transport}{,}
    \item{health}{,}
    \item{entertainment}{,}
    \item{perscare}{,}
    \item{reading}{,}
    \item{education}{,}
    \item{tobacco}{,}
    \item{miscexp}{,}
    \item{cashcont}{,}
    \item{insurance}{,}
    \item{shows}{,}
    \item{foodaway}{,}
    \item{vacations}{.}
  }    
}
\source{ Consumer Expenditure Survey (CE), program of the US Bureau of
    Labor Statistics \url{http://www.bls.gov/cex/}, diary survey.
}

\keyword{datasets}
