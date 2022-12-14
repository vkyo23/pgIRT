% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_init.R
\name{make_init}
\alias{make_init}
\title{Generating initial values for pgIRT model}
\usage{
make_init(data, model = c("default", "dynamic"), T = NULL, constraint = NULL)
}
\arguments{
\item{data}{a matrix object, roll-call matrix}

\item{model}{string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)}

\item{T}{an integer, indicating the number of sessions. For dynamic model, you must specify this argument.}

\item{constraint}{an integer or integer vector (for dynamic model, T-length), indicating the voter whose ideal point is always set positive.}
}
\description{
\code{make_init} is an auxiliary function for pgIRT model and generates optimal initial values based on roll-call matrix.
}
\examples{
\dontrun{
data(m_data_dyn)
mat <- make_rollcall(m_data_dyn,
                     unit = "unit",
                     bill = "bill",
                     vote_col = "vote")
mat <- clean_rollcall(mat)
init <- make_init(mat,
                  model = "default",
                  constraint = 1)
                  
}
}
