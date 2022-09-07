#' @title Simulated voting data (multinomial)
#' @description \code{m_data} is simulated voting dataset
#'
#' @format a matrix, rows are voters (I = 100) and cols are bills (J = 15). 
"m_data"

#' @title Simulated voting long-format data (dynamic multinomial)
#' @description \code{m_data_dyn} is simulated voting long-format dataset
#'
#' @format a tbl_df with 100 voters, 120 bills and 10 sessions. Observations are 11,597 (NAs are dropped) with 4 columns:
#'  \describe{
#'   \item{unit}{IDs of voters}
#'   \item{time}{IDs of sessions}
#'   \item{bill}{IDs of bills}
#'   \item{vote}{Voting choice of a voter}
#' }
"m_data_dyn"

#' @title Simulated matched bills (dynamic multinomial)
#' @description \code{sim_match} is simulated matched bills for dynamic estimation
#'
#' @format a tbl_df with 20 matched votes (40 x 2):
#'  \describe{
#'   \item{rc1}{IDs of bills which have a matched bill}
#'   \item{rc2}{IDs of bills which are matched with \code{rc1}}
#' }
"sim_match"