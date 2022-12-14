% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_rollcall.R
\name{clean_rollcall}
\alias{clean_rollcall}
\title{Cleaning roll-call matrix for pgIRT model}
\usage{
clean_rollcall(data)
}
\arguments{
\item{data}{a matrix, data.frame or tbl_df of roll-call (i.e. individual-item) matrix (I x J).}
}
\description{
\code{clean_rolllcall} is an auxiliary function for pgIRT model and cleans a roll-call matrix
by removing noisy items and individuals who have no response records.
}
\examples{
\dontrun{
data(m_data_dyn)
mat <- make_rollcall(m_data_dyn,
                     unit = "unit",
                     bill = "bill",
                     vote_col = "vote")
mat <- clean_rollcall(mat)
}
}
