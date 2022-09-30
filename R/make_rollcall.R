#' @title Generating roll-call matrix for pgIRT model
#' @description \code{make_rolllcall} is an auxiliary function for pgIRT model and generates roll-call matrix from long-dataframe.
#' 
#' @param dataframe data.frame, tbl_df or matrix object, dataframe (long-format, i.e. voter-bill unit) of voting data.
#' @param unit_id a string, column name which indicates voters. Integer or numeric column is allowed.
#' @param bill_id a string, column name which indicates bills. Only integer or numeric column is allowed.
#' @param vote_col a string, column name which indicates the votes. Only integer or numeric column is allowed.
#' 
#' @importFrom dplyr select arrange %>%
#' @importFrom rlang sym enquo !! :=
#' @importFrom tidyr pivot_wider
#' @export
#' 
#' @examples
#' \dontrun{
#' data(m_data_dyn)
#' mat <- make_rollcall(m_data_dyn,
#'                      unit = "unit",
#'                      bill = "bill",
#'                      vote_col = "vote")
#' }

make_rollcall <- function(dataframe, 
                          unit_id = NULL, 
                          bill_id = NULL, 
                          vote_col = NULL) {
  
  if (!class(dataframe)[1] %in% c("data.frame", "tbl_df")) {
    stop("`dataframe` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    dataframe <- as_tibble(dataframe)
  }
  if (is.null(unit_id)) stop("Please specify `unit_id`.")
  if (is.null(bill_id)) stop("Please specify `bill_id`.")
  if (is.null(vote_col)) stop("Please specify `vote_col`.")
  unit_id <- rlang::sym(unit_id)
  bill_id <- rlang::enquo(bill_id)
  vote_col <- rlang::enquo(vote_col)

  temp <- dataframe %>%
    dplyr::select(!!unit_id, !!bill_id, !!vote_col) %>%
    tidyr::pivot_wider(names_from = !!bill_id, values_from = !!vote_col) %>%
    dplyr::arrange(!!unit_id) %>%
    as.matrix()

  rname <- temp[, 1]

  temp <- apply(temp[, -1], 2, as.numeric)

  rownames(temp) <- rname
  cat('* Created', nrow(temp), 'x', ncol(temp), 'matrix.\n')
  return(temp)
}


