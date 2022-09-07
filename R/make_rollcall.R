#' @title Generating roll-call matrix for pgIRT model
#' @description \code{make_rolllcall} is an auxiliary function for pgIRT model and generates roll-call matrix from long-dataframe.
#' 
#' @param dataframe data.frame, tbl_df or matrix object, dataframe (long-format, i.e. voter-bill unit) of voting data.
#' @param unit_id a string, column name which indicates voters. Integer or numeric column is allowed.
#' @param bill_id a string, column name which indicates bills. Only integer or numeric column is allowed.
#' @param vote_col a string, column name which indicates the votes. Only integer or numeric column is allowed.
#' @param drop_unanimous a bool, whether removing unanimous bills or not. Default is FALSE.
#' 
#' @importFrom dplyr select arrange %>%
#' @importFrom rlang sym enquo !! :=
#' @importFrom tidyr pivot_wider
#' @export

make_rollcall <- function(dataframe, 
                          unit_id = NULL, 
                          bill_id = NULL, 
                          vote_col = NULL, 
                          drop_unanimous = FALSE) {
  
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

  if (drop_unanimous) {
    category <- apply(temp, 2, function(x) unique(x[!is.na(x)]))
    num_cat <- lapply(category, length) %>% unlist()
    max_cat <- lapply(category, max) %>% unlist()
    drop_res <- names(max_cat)[which(num_cat == 1)]

    if (length(drop_res) != 0) {
      cat("Remove some bills because they are unanimous votings:", drop_res, "\n")
      temp <- temp[, which(!colnames(temp) %in% drop_res)]
      num_cat <- num_cat[which(!names(num_cat) %in% drop_res)]
      max_cat <- max_cat[which(!names(max_cat) %in% drop_res)]
    }
  }
  no_attend <- which(apply(temp, 1, function(x) all(is.na(x))))
  if (length(no_attend) != 0) {
    cat("Remove some units who have no voting record:", no_attend, "\n")
    temp <- temp[-no_attend, ]
    rname <- rname[-no_attend]
    rownames(temp) <- rname
  }
  temp <- temp[, (apply(temp, 2, function(x) min(x, na.rm = TRUE) == 1) | apply(temp, 2, function(x) min(x, na.rm = TRUE) == 0))]
  return(temp)
}


