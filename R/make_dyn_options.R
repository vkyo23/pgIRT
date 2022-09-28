#' @title Generating options for dynamic pgIRT model
#' @description \code{make_dyn_options} is an auxiliary function for dynamic pgIRT model and generates options for such models.
#'
#' @param dataframe data.frame, tbl_df or matrix object, dataframe (long-format, i.e. voter-bill unit) of voting data.
#' @param unit_id a string, column name which indicates voters. Integer or numeric column is allowed.
#' @param bill_id a string, column name which indicates bills. Only integer or numeric column is allowed.
#' @param time_id a string, column name which indicates sessions. Only integer or numeric column is allowed.
#' @param vote_col a string, column name which indicates the votes. Only integer or numeric column is allowed.
#' @param add_matched_bill (Optional). An integer vector which indicates matched bills (the location of matched bills, starts from 1) 
#' and is the same length as the number of bills. If a bill does not have matches, please fill NA. Using \link{make_bill_match} is strongly recommended.
#' @param clean a bool, whether removing noisy bills or not. Default is FALSE.
#' 
#' @importFrom rlang sym enquo !! :=
#' @importFrom dplyr select arrange mutate filter distinct %>%
#' @importFrom tidyr pivot_wider
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' data(m_data_dyn)
#' ops <- make_dyn_options(m_data_dyn,
#'                         unit_id = "unit",
#'                         bill_id = "bill",
#'                         time_id = "time",
#'                         vote_col = "vote")
#' }         

make_dyn_options <- function(dataframe, 
                             unit_id = NULL, 
                             bill_id = NULL, 
                             time_id = NULL, 
                             vote_col = NULL,
                             add_matched_bill = NULL, 
                             clean = FALSE) {
  
  if (!class(dataframe)[1] %in% c("data.frame", "tbl_df")) {
    stop("`dataframe` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    dataframe <- as_tibble(dataframe)
  }
  if (is.null(unit_id)) stop("Please specify `unit_id`.")
  if (is.null(bill_id)) stop("Please specify `bill_id`.")
  if (is.null(time_id)) stop("Please specify `time_id`.")
  if (is.null(vote_col)) stop("Please specify `vote_col`.")
  unit_id <- rlang::sym(unit_id)
  bill_id <- rlang::sym(bill_id)
  time_id <- rlang::sym(time_id)
  vote_col <- rlang::enquo(vote_col)
  
  `%!in%` <- function(x, y) !x %in% y
  
  temp <- dataframe %>%
    dplyr::select(!!unit_id, !!bill_id, !!vote_col) %>%
    tidyr::pivot_wider(names_from = !!bill_id, values_from = !!vote_col) %>%
    dplyr::arrange(!!unit_id) %>% 
    as.matrix()
  
  rname <- temp[, 1]
  
  temp <- apply(temp[, -1], 2, as.numeric)
  
  rownames(temp) <- rname
  
  
  
  if (clean) {
    category <- apply(temp, 2, function(x) unique(x[!is.na(x)]))
    num_cat <- lapply(category, length) %>% unlist()
    drop_res <- names(num_cat)[which(num_cat == 1)]
    
    no_attend <- which(apply(temp, 1, function(x) all(is.na(x))))
    
    if (length(drop_res) != 0) {
      cat("Remove some bills because they are unanimous votings:", drop_res, "\n")
      dataframe <- dataframe %>%
        filter(!!bill_id %!in% drop_res)
    }
    dr <- which(!apply(temp, 2, function(x) min(x, na.rm = TRUE) == 1))
    if (length(dr) != 0) {
      cat("Remove some bills which do not start from 1:", names(dr), "\n")
      dataframe <- dataframe %>%
        dplyr::filter(!!bill_id %!in% names(dr))
    }
    if (length(no_attend) != 0) {
      cat("Remove some units who have no voting record:", no_attend, "\n")
      dataframe <- dataframe %>%
        filter(!!unit_id %!in% no_attend)
    }
  }
  
  ## Time-individual map (matrix) -----
  tmp1 <- dataframe %>%
    dplyr::arrange(!!unit_id) %>%
    dplyr::select(!!unit_id, !!bill_id, !!time_id) %>%
    dplyr::mutate(id := paste0(!!unit_id, "-", !!time_id),
                  time := !!time_id - min(!!time_id)) %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::arrange(!!unit_id) %>%
    dplyr::select(!!time_id, !!unit_id) %>%
    as.matrix()
  tmp1[, 1] <- as.numeric(tmp1[, 1]) - as.numeric(min(tmp1[, 1])) + 1
  
  tmp2 <- dataframe %>%
    dplyr::mutate(time := !!time_id - min(!!time_id)) %>%
    dplyr::distinct(!!bill_id, .keep_all = TRUE) %>%
    dplyr::select(!!bill_id, !!time_id) %>%
    as.matrix()
  tmp2 <- as.numeric(tmp2[, 2]) - as.numeric(min(tmp2[, 2])) + 1
  
  L <- list(session_individual = tmp1,
            bill_session = tmp2)
  
  if (!is.null(add_matched_bill)) {
    L$matched_bill <- add_matched_bill
  } else {
    L$matched_bill <- NULL
  }
  return(L)
}


