% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqmbgx.r
\name{seqmbgx}
\alias{seqmbgx}
\title{Maurer-Bretz sequential graphical approach}
\usage{
seqmbgx(
  xm = qnorm(matrix(rep(c(0.03, 0.04, 0.01), times = 4), ncol = 3, nrow = 4)),
  informationm = matrix(rep(c(0.4, 0.8, 1), each = 4), ncol = 3, nrow = 4),
  spending = rep("OBF", 4),
  param.spending = rep(1, 4),
  alpha = 0.025,
  sided = -1,
  W = c(0.5, 0.5, 0, 0),
  G = rbind(c(0, 0, 1, 0), c(0, 0, 0, 1), c(0, 1, 0, 0), c(1, 0, 0, 0)),
  tol = 1e-10,
  retrospective = 0
)
}
\arguments{
\item{xm}{Numeric matrix of test statistics for each endpoint (in row) and
each time point (in column).}

\item{informationm}{Numeric matrix of information fractions for the
statistics \code{xm}.}

\item{spending}{Character vector for the type(s) of the spending function for
each endpoint.}

\item{param.spending}{parameter in the spending function}

\item{alpha}{overall family-wise error rate}

\item{sided}{Integer scalar indicating the side of the test:
\itemize{
\item \code{-1}: Reject if test statistic is smaller than or equal to the critical value (one-sided)
\item \code{1}: Reject if test statistic is greater or equal to the critical value (one-sided)
\item \code{0}: Reject if the absolute value of the test statistic is greater than the critical value (two-sided)
}}

\item{W}{Numeric vector of the weights of the graph.}

\item{G}{Numeric transition matrix of the graph.}

\item{tol}{Numeric scalar of tolerance level for computing the critical
values.}

\item{retrospective}{Integer scalar with the following potential values
\itemize{
\item \code{0}: (default) only compares the current test statistic with the updated
critical value;
\item \code{1}: compares all the test statistics up to the current one with the updated
critical values.
}

Even though retrospectively looking at the values is statistically valid
in terms of control the type-1 error, not retrospectively looking at the
past comparisons avoids the dilemma of retrospectively increasing the alpha
level for the un-rejected hypothesis in the past.}
}
\value{
List with elements
\itemize{
\item \code{Hrej}: rejected hypotheses
\item \code{rejected}: the index set of rejected hypotheses
\item \code{decisionsm}: rejection decision for each endpoint (row) at each timepoint (column)
\item \code{cumdecisionsm}: cumulative rejection decision for each endpoint (row) at
each timepoint (column)
}
}
\description{
\code{seqmbgx()} conducts group-sequential testing for multiple hypotheses based
on Maurer-Bretz approach.
}
\examples{
seqmbgx(
  xm = qnorm(matrix(rep(c(0.03, 0.04, 0.01), times = 4), ncol = 3, nrow = 4)),
  informationm = matrix(rep(c(0.4, 0.8, 1), each = 4), ncol = 3, nrow = 4),
  spending = rep("OBF", 4),
  param.spending = rep(1, 4),
  alpha = 0.025,
  W = c(0.5, 0.5, 0, 0),
  G = rbind(c(0, 0, 1, 0), c(0, 0, 0, 1), c(0, 1, 0, 0), c(1, 0, 0, 0)),
  retrospective = 0
)
}
\author{
Xiaodong Luo
}
\concept{Maurer-Bretz}
\concept{graphical procedure}
\concept{group-sequential}
