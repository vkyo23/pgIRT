% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_dyn_options.R
\name{make_dyn_options}
\alias{make_dyn_options}
\title{Generating options for dynamic pgIRT model}
\usage{
make_dyn_options(
  dataframe,
  unit_id = NULL,
  bill_id = NULL,
  time_id = NULL,
  vote_col = NULL,
  add_matched_bill = NULL,
  clean = FALSE
)
}
\arguments{
\item{dataframe}{data.frame, tbl_df or matrix object, dataframe (long-format, i.e. voter-bill unit) of voting data.}

\item{unit_id}{a string, column name which indicates voters. Integer or numeric column is allowed.}

\item{bill_id}{a string, column name which indicates bills. Only integer or numeric column is allowed.}

\item{time_id}{a string, column name which indicates sessions. Only integer or numeric column is allowed.}

\item{vote_col}{a string, column name which indicates the votes. Only integer or numeric column is allowed.}

\item{add_matched_bill}{(Optional). An integer vector which indicates matched bills (the location of matched bills, starts from 1)
and is the same length as the number of bills. If a bill does not have matches, please fill NA. Using \link{make_bill_match} is strongly recommended.}

\item{clean}{a bool, whether removing noisy bills or not. Default is FALSE.}
}
\description{
\code{make_dyn_options} is an auxiliary function for dynamic pgIRT model and generates options for such models.
}
\examples{
\dontrun{
data(m_data_dyn)
ops <- make_dyn_options(m_data_dyn,
                        unit_id = "unit",
                        bill_id = "bill",
                        time_id = "time",
                        vote_col = "vote")
}         
}
